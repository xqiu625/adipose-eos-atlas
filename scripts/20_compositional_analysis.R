#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))
if (nzchar(Sys.getenv("RERUN_R_LIB"))) {
  .rerun_r_lib_resolved <- normalizePath(.libPaths()[1], mustWork = FALSE)
  .rerun_r_lib_expected <- normalizePath(Sys.getenv("RERUN_R_LIB"), mustWork = FALSE)
  if (!identical(.rerun_r_lib_resolved, .rerun_r_lib_expected)) {
    stop(sprintf(
      "RERUN_R_LIB audit guard: .libPaths()[1] resolved to '%s' but RERUN_R_LIB='%s' -- the shared/legacy library would take precedence over the pinned isolated library. Refusing to proceed.",
      .rerun_r_lib_resolved, .rerun_r_lib_expected))
  }
  if (!nzchar(Sys.getenv("CMDSTAN"))) {
    stop("RERUN_R_LIB audit guard: this script loads sccomp, which requires CMDSTAN to be set when running under the audited rerun.")
  }
  cat("=== Pinned R environment audit guard ===\n")
  cat("R.version.string:", R.version.string, "\n")
  cat("RERUN_R_LIB:", Sys.getenv("RERUN_R_LIB"), "\n")
  cat("CMDSTAN:", Sys.getenv("CMDSTAN"), "\n")
  cat(".libPaths():\n")
  print(.libPaths())
  cat("=========================================\n")
}

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --job-name=comp_analysis
#SBATCH --output=comp_analysis_%j.out

# ==============================================================================
# 20_compositional_analysis.R
# Bayesian compositional analysis of cluster proportions across conditions
# using sccomp (Bayesian compositional analysis)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
comp_dir <- file.path(BASE, "composition")
fig_dir  <- file.path(BASE, "figures")
dir.create(comp_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# eos_colors loaded from 00_theme.R

cat("========================================\n")
cat("20_compositional_analysis.R\n")
cat("Compositional analysis of cluster proportions\n")
cat("========================================\n\n")

# --- Load data ----------------------------------------------------------------
cat("Loading eos_subclustered.rds...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
# Remove C6 if still present (too few cells, n=14)
if ("C6" %in% levels(Idents(eos))) {
  eos <- subset(eos, idents = "C6", invert = TRUE)
  Idents(eos) <- droplevels(Idents(eos))
  eos$eos_cluster <- droplevels(eos$eos_cluster)
  cat(sprintf("  Removed C6. Remaining: %d eos (C0-C5)\n", ncol(eos)))
}
cat(sprintf("  %d cells, %d clusters\n", ncol(eos), length(unique(eos$eos_cluster))))

# --- Build composition table --------------------------------------------------
cat("Building composition table...\n")

comp_data <- eos@meta.data %>%
  group_by(orig.ident, eos_cluster) %>%
  summarise(count = n(), .groups = "drop")

# Add sample-level metadata (one row per sample)
sample_meta <- eos@meta.data %>%
  select(orig.ident, condition, sex, genotype, infection) %>%
  distinct()

# Ensure all cluster-sample combinations exist. group_by/summarise above only
# emits rows for observed combinations, so a sample with zero cells in a given
# cluster is silently absent rather than modeled as a true biological zero.
all_combos <- expand.grid(
  orig.ident  = unique(comp_data$orig.ident),
  eos_cluster = levels(eos$eos_cluster),
  stringsAsFactors = FALSE
)
comp_data <- all_combos %>%
  left_join(comp_data, by = c("orig.ident", "eos_cluster")) %>%
  mutate(count = replace_na(count, 0)) %>%
  left_join(sample_meta, by = "orig.ident")

# Compute proportions per sample
comp_data <- comp_data %>%
  group_by(orig.ident) %>%
  mutate(
    total = sum(count),
    proportion = count / total
  ) %>%
  ungroup()

cat(sprintf("  %d samples, %d cluster-sample combinations\n",
            length(unique(comp_data$orig.ident)), nrow(comp_data)))
cat(sprintf("  Conditions: %s\n", paste(unique(comp_data$condition), collapse = ", ")))

# --- Try sccomp; fall back to per-cluster glm if unavailable ----------------------------
use_sccomp <- requireNamespace("sccomp", quietly = TRUE)

# Detect which sccomp API is available (v1: sccomp_glm, v2: sccomp_estimate)
sccomp_func <- NULL
if (use_sccomp) {
  library(sccomp)
  if (exists("sccomp_estimate", mode = "function")) {
    sccomp_func <- "sccomp_estimate"
    cat("\nUsing sccomp (v2 API: sccomp_estimate)...\n")
  } else if (exists("sccomp_glm", mode = "function")) {
    sccomp_func <- "sccomp_glm"
    cat("\nUsing sccomp (v1 API: sccomp_glm)...\n")
  } else {
    cat("\nsccomp loaded but no known function found. Falling back to per-cluster glm.\n")
    use_sccomp <- FALSE
  }
}

# Helper to call the right sccomp function
run_sccomp <- function(data, formula_comp, formula_var = ~ 1) {
  if (sccomp_func == "sccomp_estimate") {
    data %>% sccomp_estimate(
      formula_composition = formula_comp,
      formula_variability = formula_var,
      .sample = orig.ident,
      .cell_group = eos_cluster,
      .count = count
    )
  } else {
    data %>% sccomp_glm(
      formula_composition = formula_comp,
      formula_variability = formula_var,
      .sample = orig.ident,
      .cell_group = eos_cluster,
      .count = count
    )
  }
}

if (use_sccomp) {

  # --- Test 1: Sex effect (WT, HFD only) ---------------------------------------
  cat("\n--- Test 1: Sex effect (WT, HFD only) ---\n")
  data_sex <- comp_data %>% filter(infection == "HFD", genotype == "WT")
  cat(sprintf("  Samples: %d, Cells in total: %d\n",
              length(unique(data_sex$orig.ident)), sum(data_sex$count)))

  res_sex <- tryCatch({
    run_sccomp(data_sex, ~ sex)
  }, error = function(e) {
    cat(sprintf("  sccomp error: %s\n", e$message))
    NULL
  })

  # Helper to run sccomp_test and save results as flat CSV
  save_sccomp <- function(res, test_name) {
    if (is.null(res)) return()
    # Run hypothesis testing
    res_tested <- tryCatch({
      sccomp_test(res)
    }, error = function(e) {
      cat(sprintf("  sccomp_test error: %s\n", e$message))
      res  # fall back to untested result
    })
    # Convert to data frame, dropping list columns
    res_df <- as.data.frame(res_tested)
    list_cols <- sapply(res_df, is.list)
    if (any(list_cols)) {
      cat(sprintf("  Dropping list columns: %s\n", paste(names(res_df)[list_cols], collapse=", ")))
      res_df <- res_df[, !list_cols, drop = FALSE]
    }
    out_path <- file.path(comp_dir, paste0("sccomp_", test_name, ".csv"))
    write.csv(res_df, out_path, row.names = FALSE)
    cat(sprintf("  Saved: sccomp_%s.csv\n", test_name))
    # Print results
    cat(sprintf("  Results for %s:\n", test_name))
    print(res_df)
    invisible(res_df)
  }

  save_sccomp(res_sex, "sex")

  # --- Test 2: Genotype effect (within females, HFD) --------------------------
  cat("\n--- Test 2: Genotype effect (Females, HFD) ---\n")
  data_geno <- comp_data %>% filter(sex == "Female", infection == "HFD")
  cat(sprintf("  Samples: %d, Cells in total: %d\n",
              length(unique(data_geno$orig.ident)), sum(data_geno$count)))

  res_geno <- tryCatch({
    run_sccomp(data_geno, ~ genotype)
  }, error = function(e) {
    cat(sprintf("  sccomp error: %s\n", e$message))
    NULL
  })

  save_sccomp(res_geno, "genotype")

  # --- Test 3: Infection effect (within females) ------------------------------
  cat("\n--- Test 3: Infection effect (Females) ---\n")
  data_inf <- comp_data %>% filter(sex == "Female")
  cat(sprintf("  Samples: %d, Cells in total: %d\n",
              length(unique(data_inf$orig.ident)), sum(data_inf$count)))

  res_inf <- tryCatch({
    run_sccomp(data_inf, ~ infection)
  }, error = function(e) {
    cat(sprintf("  sccomp error: %s\n", e$message))
    NULL
  })

  save_sccomp(res_inf, "infection")

} else {
  cat("\nsccomp not available. Falling back to per-cluster binomial GLM (glm)...\n")
  cat("WARNING: this is an approximate frequentist substitute, not a repeat of\n")
  cat("the sccomp Bayesian compositional model. Outputs are flagged with a\n")
  cat("'method' column and should not be treated as equivalent to sccomp results.\n")

  run_glm_composition <- function(data, formula_str, test_name) {
    clusters <- sort(unique(data$eos_cluster))
    covariate_cols <- setdiff(colnames(data),
                               c("eos_cluster", "count", "total", "proportion"))
    results <- list()

    for (cl in clusters) {
      cat(sprintf("  Cluster %s...\n", cl))

      # One row per sample: successes = cells in cluster cl, failures = cells
      # in every other cluster, so successes + failures = total cells in that
      # sample (sum-constrained). This is the correct binomial representation
      # of "cluster cl's share of the sample"; the previous construction gave
      # each cluster-row its own denominator, which is not sum-constrained and
      # does not represent cluster count vs. total-minus-cluster count.
      cl_data <- data %>%
        group_by(orig.ident) %>%
        summarise(
          successes = sum(count[eos_cluster == cl]),
          failures  = sum(count[eos_cluster != cl]),
          across(all_of(covariate_cols), dplyr::first),
          .groups = "drop"
        )

      # One row per sample means there is no repeated-measures structure left
      # to justify a (1|orig.ident) random effect, so this is a plain glm.
      fit <- tryCatch({
        glm(
          as.formula(paste("cbind(successes, failures) ~", formula_str)),
          data = cl_data,
          family = binomial
        )
      }, error = function(e) {
        cat(sprintf("    glm error: %s\n", e$message))
        NULL
      })

      if (!is.null(fit)) {
        coefs <- summary(fit)$coefficients
        if (nrow(coefs) >= 2) {
          results[[cl]] <- data.frame(
            cluster = cl,
            term = rownames(coefs)[2],
            estimate = coefs[2, "Estimate"],
            std_error = coefs[2, "Std. Error"],
            z_value = coefs[2, "z value"],
            p_value = coefs[2, "Pr(>|z|)"],
            method = "glm_fallback",
            stringsAsFactors = FALSE
          )
        }
      }
    }

    if (length(results) > 0) {
      res_df <- do.call(rbind, results)
      res_df$fdr <- p.adjust(res_df$p_value, method = "BH")
      rownames(res_df) <- NULL
      return(res_df)
    }
    return(NULL)
  }

  # Test 1: Sex effect
  cat("\n--- Test 1: Sex effect (WT, HFD only) ---\n")
  data_sex <- comp_data %>% filter(infection == "HFD", genotype == "WT")
  res_sex_df <- run_glm_composition(data_sex, "sex", "sex")
  if (!is.null(res_sex_df)) {
    write.csv(res_sex_df, file.path(comp_dir, "sccomp_sex.csv"), row.names = FALSE)
    cat("  Saved: sccomp_sex.csv (method=glm_fallback)\n")
  }

  # Test 2: Genotype effect
  cat("\n--- Test 2: Genotype effect (Females, HFD) ---\n")
  data_geno <- comp_data %>% filter(sex == "Female", infection == "HFD")
  res_geno_df <- run_glm_composition(data_geno, "genotype", "genotype")
  if (!is.null(res_geno_df)) {
    write.csv(res_geno_df, file.path(comp_dir, "sccomp_genotype.csv"), row.names = FALSE)
    cat("  Saved: sccomp_genotype.csv (method=glm_fallback)\n")
  }

  # Test 3: Infection effect
  cat("\n--- Test 3: Infection effect (Females) ---\n")
  data_inf <- comp_data %>% filter(sex == "Female")
  res_inf_df <- run_glm_composition(data_inf, "infection", "infection")
  if (!is.null(res_inf_df)) {
    write.csv(res_inf_df, file.path(comp_dir, "sccomp_infection.csv"), row.names = FALSE)
    cat("  Saved: sccomp_infection.csv (method=glm_fallback)\n")
  }
}

# ==============================================================================
# Figures
# ==============================================================================
cat("\n========================================\n")
cat("Generating figures...\n")
cat("========================================\n")

# --- Figure 1: Stacked bar plot (cluster proportions per sample, grouped) -----
cat("Figure 1: Stacked bar plot (proportions per sample)...\n")

# Order samples by condition for cleaner grouping
sample_order <- sample_meta %>%
  arrange(condition) %>%
  pull(orig.ident)

prop_df <- comp_data %>%
  mutate(orig.ident = factor(orig.ident, levels = sample_order))

p1 <- ggplot(prop_df, aes(x = orig.ident, y = proportion, fill = eos_cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = eos_colors) +
  labs(x = "Sample", y = "Proportion", fill = "Cluster",
       title = "Eosinophil Cluster Proportions per Sample") +
  theme_pub(base_size = 12) +
  theme(
    axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
    legend.position = "right"
  ) +
  facet_grid(. ~ condition, scales = "free_x", space = "free_x")

pdf(file.path(fig_dir, "20_stacked_bar_composition.pdf"), width = 14, height = 6)
print(p1)
dev.off()
cat("  Saved: 20_stacked_bar_composition.pdf\n")

# --- Figure 2: Forest plot (effect sizes with credible/confidence intervals) --
cat("Figure 2: Forest plot (effect sizes)...\n")

# Combine results from all tests for forest plot
all_results <- list()
for (test_name in c("sex", "genotype", "infection")) {
  csv_path <- file.path(comp_dir, paste0("sccomp_", test_name, ".csv"))
  if (file.exists(csv_path)) {
    df <- read.csv(csv_path)
    # Identify the cluster column (sccomp uses cell_group, glm fallback uses cluster)
    cluster_col <- intersect(c("eos_cluster", "cell_group", "cluster"), colnames(df))
    estimate_col <- intersect(c("c_effect", "estimate"), colnames(df))
    lower_col <- intersect(c("c_lower", "lower"), colnames(df))
    upper_col <- intersect(c("c_upper", "upper"), colnames(df))

    if (length(cluster_col) > 0 && length(estimate_col) > 0) {
      plot_df <- data.frame(
        cluster = df[[cluster_col[1]]],
        estimate = df[[estimate_col[1]]],
        test = test_name,
        stringsAsFactors = FALSE
      )
      # Add CI if available (sccomp), otherwise compute from SE (glm fallback)
      if (length(lower_col) > 0 && length(upper_col) > 0) {
        plot_df$lower <- df[[lower_col[1]]]
        plot_df$upper <- df[[upper_col[1]]]
      } else if ("std_error" %in% colnames(df)) {
        plot_df$lower <- df$estimate - 1.96 * df$std_error
        plot_df$upper <- df$estimate + 1.96 * df$std_error
      } else {
        plot_df$lower <- plot_df$estimate
        plot_df$upper <- plot_df$estimate
      }
      # Add significance
      if ("fdr" %in% colnames(df)) {
        plot_df$significant <- df$fdr < 0.1
      } else if ("c_FDR" %in% colnames(df)) {
        plot_df$significant <- df$c_FDR < 0.1
      } else {
        plot_df$significant <- FALSE
      }
      all_results[[test_name]] <- plot_df
    }
  }
}

if (length(all_results) > 0) {
  forest_df <- do.call(rbind, all_results)
  rownames(forest_df) <- NULL

  forest_df$test <- factor(forest_df$test,
                           levels = c("sex", "genotype", "infection"),
                           labels = c("Sex (M vs F)", "Genotype (KO vs WT)",
                                      "Infection (Infected vs HFD)"))

  p2 <- ggplot(forest_df, aes(x = estimate, y = cluster, color = significant)) +
    geom_vline(xintercept = 0, linetype = "dashed", color = "grey50") +
    geom_errorbarh(aes(xmin = lower, xmax = upper), height = 0.2, linewidth = 0.6) +
    geom_point(size = 2.5) +
    scale_color_manual(values = c("TRUE" = "#E84A5F", "FALSE" = "grey60"),
                       labels = c("TRUE" = "FDR < 0.1", "FALSE" = "NS")) +
    facet_wrap(~ test, scales = "free_x") +
    labs(x = "Log fold change in composition", y = "Cluster",
         color = "Significance",
         title = "Compositional Changes Across Conditions") +
    theme_pub(base_size = 12) +
    theme(
      legend.position = "bottom"
    )

  pdf(file.path(fig_dir, "20_forest_composition.pdf"), width = 12, height = 5)
  print(p2)
  dev.off()
  cat("  Saved: 20_forest_composition.pdf\n")
} else {
  cat("  Skipping forest plot (no results available)\n")
}

# --- Figure 3: Dot plot of significant composition changes --------------------
cat("Figure 3: Dot plot (significant changes, FDR < 0.1)...\n")

if (length(all_results) > 0) {
  sig_df <- forest_df %>% filter(significant == TRUE)

  if (nrow(sig_df) > 0) {
    p3 <- ggplot(sig_df, aes(x = cluster, y = test)) +
      geom_point(aes(size = abs(estimate), color = estimate)) +
      scale_color_gradient2(low = "#2166AC", mid = "white", high = "#B2182B",
                            midpoint = 0, name = "Effect size") +
      scale_size_continuous(range = c(2, 8), name = "|Effect size|") +
      labs(x = "Cluster", y = "Comparison",
           title = "Significant Compositional Changes (FDR < 0.1)") +
      theme_pub(base_size = 12) +
      theme(
        axis.text.x = element_text(angle = 45, hjust = 1)
      )

    pdf(file.path(fig_dir, "20_dotplot_significant.pdf"), width = 8, height = 4)
    print(p3)
    dev.off()
    cat("  Saved: 20_dotplot_significant.pdf\n")
  } else {
    cat("  No significant changes at FDR < 0.1; skipping dot plot\n")
  }
} else {
  cat("  Skipping dot plot (no results available)\n")
}

# --- Figure 4: Boxplot of cluster proportions across conditions ---------------
cat("Figure 4: Boxplot (per-cluster proportions by condition)...\n")

# Create condition groups for faceting
box_df <- comp_data %>%
  mutate(
    comparison_sex = ifelse(infection == "HFD" & genotype == "WT", sex, NA),
    comparison_geno = ifelse(infection == "HFD" & sex == "Female", genotype, NA),
    comparison_inf = ifelse(sex == "Female", infection, NA)
  )

# Panel A: Sex effect
box_sex <- box_df %>%
  filter(!is.na(comparison_sex)) %>%
  ggplot(aes(x = eos_cluster, y = proportion, fill = sex)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, position = position_dodge(0.75)) +
  geom_point(aes(color = sex), position = position_jitterdodge(jitter.width = 0.1,
                                                                dodge.width = 0.75),
             size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = c("Female" = "#E84A5F", "Male" = "#6CA6CD")) +
  scale_color_manual(values = c("Female" = "#E84A5F", "Male" = "#6CA6CD")) +
  labs(x = "Cluster", y = "Proportion", title = "Sex Effect (WT, HFD only)") +
  theme_pub(base_size = 11) +
  theme(legend.position = "bottom")

# Panel B: Genotype effect
box_geno <- box_df %>%
  filter(!is.na(comparison_geno)) %>%
  ggplot(aes(x = eos_cluster, y = proportion, fill = genotype)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, position = position_dodge(0.75)) +
  geom_point(aes(color = genotype), position = position_jitterdodge(jitter.width = 0.1,
                                                                     dodge.width = 0.75),
             size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = c("WT" = "#5BA4A4", "KO" = "#FF6F1F")) +
  scale_color_manual(values = c("WT" = "#5BA4A4", "KO" = "#FF6F1F")) +
  labs(x = "Cluster", y = "Proportion", title = "Genotype Effect (Females, HFD)") +
  theme_pub(base_size = 11) +
  theme(legend.position = "bottom")

# Panel C: Infection effect
box_inf <- box_df %>%
  filter(!is.na(comparison_inf)) %>%
  ggplot(aes(x = eos_cluster, y = proportion, fill = infection)) +
  geom_boxplot(outlier.size = 0.5, alpha = 0.7, position = position_dodge(0.75)) +
  geom_point(aes(color = infection), position = position_jitterdodge(jitter.width = 0.1,
                                                                      dodge.width = 0.75),
             size = 1.5, alpha = 0.8) +
  scale_fill_manual(values = c("HFD" = "#CBAACB", "HFDInfected" = "#B03A48")) +
  scale_color_manual(values = c("HFD" = "#CBAACB", "HFDInfected" = "#B03A48")) +
  labs(x = "Cluster", y = "Proportion", title = "Infection Effect (Females)") +
  theme_pub(base_size = 11) +
  theme(legend.position = "bottom")

p4 <- box_sex / box_geno / box_inf +
  plot_annotation(
    title = "Cluster Proportions Across Conditions",
    subtitle = "Each dot = one sample",
    theme = theme(plot.title = element_text(face = "bold", size = 14))
  )

pdf(file.path(fig_dir, "20_boxplot_proportions.pdf"), width = 10, height = 14)
print(p4)
dev.off()
cat("  Saved: 20_boxplot_proportions.pdf\n")

# --- Summary ------------------------------------------------------------------
cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
cat(sprintf("  Cells:    %d\n", ncol(eos)))
cat(sprintf("  Samples:  %d\n", length(unique(eos$orig.ident))))
cat(sprintf("  Clusters: %d\n", length(unique(eos$eos_cluster))))
cat(sprintf("  Method:   %s\n", ifelse(use_sccomp, "sccomp (Bayesian)", "glm (frequentist, fallback)")))
cat("\nOutput files:\n")
cat(sprintf("  %s/sccomp_sex.csv\n", comp_dir))
cat(sprintf("  %s/sccomp_genotype.csv\n", comp_dir))
cat(sprintf("  %s/sccomp_infection.csv\n", comp_dir))
cat(sprintf("  %s/20_stacked_bar_composition.pdf\n", fig_dir))
cat(sprintf("  %s/20_forest_composition.pdf\n", fig_dir))
cat(sprintf("  %s/20_dotplot_significant.pdf\n", fig_dir))
cat(sprintf("  %s/20_boxplot_proportions.pdf\n", fig_dir))
cat("\nDone! 20_compositional_analysis.R completed successfully.\n")
