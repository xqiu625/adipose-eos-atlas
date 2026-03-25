#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=8:00:00
#SBATCH --job-name=deg_eos
#SBATCH --output=deg_eos_%j.out

# ==============================================================================
# 07_deg_eos.R
# MAST mixed-effects DEG analysis for eosinophils (4 comparisons)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(MAST)
  library(SingleCellExperiment)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
  library(data.table)
})

SCRIPTS <- dirname(sys.frame(1)$ofile)
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
deg_dir <- file.path(BASE, "deg")
fig_dir <- file.path(BASE, "figures")
dir.create(deg_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper: run MAST with mixed effects --------------------------------------
run_mast <- function(seu, group_col, test_level, ref_level, sample_col = "orig.ident") {
  # Subset to the two groups
  cells_use <- colnames(seu)[seu@meta.data[[group_col]] %in% c(test_level, ref_level)]
  seu_sub <- subset(seu, cells = cells_use)

  cat(sprintf("    Cells: %d test (%s), %d ref (%s)\n",
              sum(seu_sub@meta.data[[group_col]] == test_level), test_level,
              sum(seu_sub@meta.data[[group_col]] == ref_level), ref_level))

  # Get normalized data and create SCE
  norm_mat <- GetAssayData(seu_sub, slot = "data")
  sce <- SingleCellExperiment(
    assays = list(logcounts = norm_mat),
    colData = seu_sub@meta.data
  )

  # Convert to SingleCellAssay
  sca <- SceToSingleCellAssay(sce, check_sanity = FALSE)

  # Add covariates
  colData(sca)$cngeneson <- scale(colSums(assay(sca) > 0))[, 1]
  colData(sca)$group <- factor(colData(sca)[[group_col]],
                               levels = c(ref_level, test_level))
  colData(sca)$sample <- factor(colData(sca)[[sample_col]])

  # Filter lowly-expressed genes (expressed in >= 10% of cells)
  freq <- freq(sca)
  sca <- sca[freq > 0.1, ]
  cat(sprintf("    Genes tested: %d\n", nrow(sca)))

  # Run MAST with random intercept per sample
  cat("    Fitting zlm model (glmer, may take a while)...\n")
  zlm_fit <- zlm(~ group + cngeneson + (1 | sample), sca,
                  method = "glmer", ebayes = FALSE, strictConvergence = FALSE)

  # Extract results for the group contrast
  coef_name <- paste0("group", test_level)
  cat(sprintf("    Extracting LRT for coefficient: %s\n", coef_name))
  summary_dt <- summary(zlm_fit, doLRT = coef_name)
  summary_table <- summary_dt$datatable

  # Get hurdle p-values and log2FC
  hurdle_pval <- summary_table[component == "H" & contrast == coef_name,
                               .(primerid, `Pr(>Chisq)`)]
  logfc_vals  <- summary_table[component == "logFC" & contrast == coef_name,
                               .(primerid, coef, ci.hi, ci.lo)]

  res <- merge(hurdle_pval, logfc_vals, by = "primerid")
  setnames(res, c("gene", "pval", "log2FC", "ci_hi", "ci_lo"))
  res$padj <- p.adjust(res$pval, method = "BH")

  # Sort by adjusted p-value
  res <- res[order(res$padj), ]
  return(as.data.frame(res))
}

# --- Helper: volcano plot -----------------------------------------------------
plot_volcano <- function(deg_df, title, filename, top_n = 20) {
  deg_df <- deg_df %>%
    mutate(
      sig = case_when(
        padj < 0.05 & log2FC >  0.5 ~ "Up",
        padj < 0.05 & log2FC < -0.5 ~ "Down",
        TRUE                         ~ "NS"
      )
    )

  # Label top genes by significance
  top_genes <- deg_df %>%
    filter(sig != "NS") %>%
    arrange(padj) %>%
    head(top_n)

  p <- ggplot(deg_df, aes(x = log2FC, y = -log10(pval), color = sig)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text_repel(data = top_genes,
                    aes(label = gene), size = 2.5, max.overlaps = 25,
                    fontface = "italic", color = "black") +
    scale_color_manual(values = vol_colors) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    labs(x = "log2 Fold Change",
         y = expression(-log[10](italic(p)-value)),
         title = title, color = "Significance") +
    theme_pub()

  pdf(file.path(fig_dir, filename), width = 8, height = 7)
  print(p)
  dev.off()
}

# --- Load eosinophil object ---------------------------------------------------
cat("Loading eos_subclustered.rds...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))

# --- Define comparisons -------------------------------------------------------
comparisons <- list(
  sex_effect = list(
    name      = "sex_effect",
    desc      = "WT_Female_HFD vs WT_Male_HFD (sex effect)",
    subset_conditions = c("WT_Female_HFD", "WT_Male_HFD"),
    group_col = "condition",
    test      = "WT_Female_HFD",
    ref       = "WT_Male_HFD"
  ),
  genotype_female = list(
    name      = "genotype_female",
    desc      = "KO_Female_HFD vs WT_Female_HFD (genotype effect, female)",
    subset_conditions = c("KO_Female_HFD", "WT_Female_HFD"),
    group_col = "condition",
    test      = "KO_Female_HFD",
    ref       = "WT_Female_HFD"
  ),
  genotype_male = list(
    name      = "genotype_male",
    desc      = "KO_Male_HFD vs WT_Male_HFD (genotype effect, male)",
    subset_conditions = c("KO_Male_HFD", "WT_Male_HFD"),
    group_col = "condition",
    test      = "KO_Male_HFD",
    ref       = "WT_Male_HFD"
  ),
  infection_female = list(
    name      = "infection_female",
    desc      = "Infected vs non-infected females",
    subset_conditions = c("WT_Female_HFDInfected", "KO_Female_HFDInfected",
                          "WT_Female_HFD", "KO_Female_HFD"),
    group_col = "infection_status",
    test      = "Infected",
    ref       = "NonInfected"
  )
)

# --- Run MAST for each comparison --------------------------------------------
for (comp in comparisons) {
  cat(sprintf("\n=== Comparison: %s ===\n", comp$desc))

  # Subset to relevant conditions
  eos_sub <- subset(eos, condition %in% comp$subset_conditions)

  # For infection comparison, create a binary group variable
  if (comp$name == "infection_female") {
    eos_sub$infection_status <- ifelse(
      grepl("Infected", eos_sub$condition), "Infected", "NonInfected"
    )
  }

  cat(sprintf("  Total cells in comparison: %d\n", ncol(eos_sub)))

  # Run MAST
  tryCatch({
    deg_res <- run_mast(eos_sub,
                        group_col = comp$group_col,
                        test_level = comp$test,
                        ref_level  = comp$ref)

    # Filter significant DEGs
    sig_res <- deg_res %>%
      filter(padj < 0.05, abs(log2FC) > 0.5)
    cat(sprintf("  Significant DEGs (padj < 0.05, |log2FC| > 0.5): %d\n",
                nrow(sig_res)))
    cat(sprintf("    Up: %d, Down: %d\n",
                sum(sig_res$log2FC > 0), sum(sig_res$log2FC < 0)))

    # Save full results
    out_file <- file.path(deg_dir, sprintf("deg_%s.csv", comp$name))
    write.csv(deg_res, out_file, row.names = FALSE)
    cat(sprintf("  Saved: %s\n", out_file))

    # Volcano plot
    plot_volcano(deg_res,
                 title    = comp$desc,
                 filename = sprintf("07_volcano_%s.pdf", comp$name))
    cat(sprintf("  Volcano plot saved.\n"))

  }, error = function(e) {
    cat(sprintf("  ERROR in comparison %s: %s\n", comp$name, e$message))
  })
}

cat("\nDone! 07_deg_eos.R completed successfully.\n")
