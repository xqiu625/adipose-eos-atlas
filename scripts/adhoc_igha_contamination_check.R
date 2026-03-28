#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

## ADHOC: IGHA ambient RNA contamination analysis
## Confirms DecontX effectively removed Ig gene contamination in eosinophils
## SLURM: --mem=64G --cpus-per-task=4 --time=2:00:00 -p batch

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(Matrix)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures", "adhoc")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Ig genes to check (mouse orthologs) ──────────────────────────────────────
# Mouse Ig genes (note: mouse uses Igha, Igkc, Ighg1, Jchain, Iglc2 etc.)
ig_genes <- c("Igha", "Igkc", "Ighg1", "Jchain", "Iglc2", "Iglc1",
              "Ighm", "Ighg2b", "Ighg2c")

# Eosinophil markers
eos_markers <- c("Siglecf", "Il5ra", "Epx", "Prg2", "Prg3", "Ear1", "Ear2", "Ccr3")

# B cell markers (to check for doublets)
bcell_markers <- c("Cd79a", "Cd79b", "Ms4a1", "Cd19", "Pax5")

cat("========================================\n")
cat("ADHOC: IGHA Contamination Analysis\n")
cat("========================================\n\n")

# ══════════════════════════════════════════════════════════════════════════════
# 1. Load post-DecontX eosinophil object
# ══════════════════════════════════════════════════════════════════════════════
cat("Step 1: Loading eos_subclustered.rds (post-DecontX)...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
cat(sprintf("  %d genes x %d cells\n", nrow(eos), ncol(eos)))

# Check which Ig genes are present
all_genes <- rownames(eos)
ig_present <- ig_genes[ig_genes %in% all_genes]
ig_missing <- ig_genes[!ig_genes %in% all_genes]
eos_markers_present <- eos_markers[eos_markers %in% all_genes]
bcell_present <- bcell_markers[bcell_markers %in% all_genes]

cat(sprintf("  Ig genes found: %s\n", paste(ig_present, collapse = ", ")))
if (length(ig_missing) > 0) {
  cat(sprintf("  Ig genes missing: %s\n", paste(ig_missing, collapse = ", ")))
}
cat(sprintf("  Eos markers found: %s\n", paste(eos_markers_present, collapse = ", ")))
cat(sprintf("  B cell markers found: %s\n", paste(bcell_present, collapse = ", ")))

# ══════════════════════════════════════════════════════════════════════════════
# 2. Ig gene expression summary in post-DecontX eosinophils
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========================================\n")
cat("Step 2: Ig gene expression in eosinophils (post-DecontX)\n")
cat("========================================\n\n")

if (length(ig_present) > 0) {
  # Raw counts
  counts_mat <- GetAssayData(eos, layer = "counts")
  # Normalized
  norm_mat <- GetAssayData(eos, layer = "data")

  ig_summary <- data.frame(gene = character(), pct_expressing = numeric(),
                           mean_raw = numeric(), mean_norm = numeric(),
                           max_raw = numeric(), stringsAsFactors = FALSE)

  for (gene in ig_present) {
    raw_vals <- counts_mat[gene, ]
    norm_vals <- norm_mat[gene, ]
    pct <- 100 * sum(raw_vals > 0) / length(raw_vals)
    ig_summary <- rbind(ig_summary, data.frame(
      gene = gene,
      pct_expressing = round(pct, 2),
      mean_raw = round(mean(raw_vals), 4),
      mean_norm = round(mean(norm_vals), 4),
      max_raw = max(raw_vals),
      stringsAsFactors = FALSE
    ))
  }

  cat("  Post-DecontX Ig gene expression in eosinophils:\n")
  print(ig_summary)

  # Identify Igha+ cells
  if ("Igha" %in% ig_present) {
    igha_pos <- counts_mat["Igha", ] > 0
    n_igha_pos <- sum(igha_pos)
    cat(sprintf("\n  Igha+ eosinophils: %d / %d (%.1f%%)\n",
                n_igha_pos, ncol(eos), 100 * n_igha_pos / ncol(eos)))
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# 3. Compare before vs after DecontX (per-sample)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========================================\n")
cat("Step 3: Before vs After DecontX comparison\n")
cat("========================================\n\n")

sample_info <- data.frame(
  sample_id = 1:16,
  condition = c("WT_Male_HFD_1",       "WT_Female_HFD_1",
                "KO_Male_HFD_1",       "KO_Female_HFD_1",
                "WT_Female_HFDInfected_1", "WT_Female_HFDInfected_2",
                "KO_Female_HFDInfected_1", "KO_Female_HFDInfected_2",
                "WT_Male_HFD_2",       "WT_Male_HFD_3",
                "WT_Female_HFD_2",     "WT_Female_HFD_3",
                "KO_Male_HFD_2",       "KO_Male_HFD_3",
                "KO_Female_HFD_2",     "KO_Female_HFD_3"),
  stringsAsFactors = FALSE
)

data_root <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data"
gem_map <- list("RR0104" = 1:4, "RR0508" = 5:8,
                "RR0912" = 9:12, "RR1316" = 13:16)

get_sample_path <- function(sample_id) {
  gem <- NULL
  for (g in names(gem_map)) {
    if (sample_id %in% gem_map[[g]]) { gem <- g; break }
  }
  file.path(data_root, gem, "outs", "per_sample_outs",
            as.character(sample_id), "count", "sample_filtered_feature_bc_matrix")
}

# Get eos cell barcodes and their sample assignments
eos_meta <- eos@meta.data
eos_meta$cell_id <- rownames(eos_meta)
eos_meta$raw_barcode <- sub("^.*_([ACGT]+-\\d+)$", "\\1", eos_meta$cell_id)

before_after_data <- list()

for (i in seq_len(nrow(sample_info))) {
  sid  <- sample_info$sample_id[i]
  cond <- sample_info$condition[i]

  eos_in_sample <- eos_meta %>% filter(orig.ident == cond)
  if (nrow(eos_in_sample) == 0) next

  # Load original (pre-DecontX) counts
  orig_path <- get_sample_path(sid)
  if (!dir.exists(orig_path)) next

  orig_counts <- tryCatch(Read10X(data.dir = orig_path), error = function(e) NULL)
  if (is.null(orig_counts)) next

  # Load DecontX-corrected counts
  decon_path <- file.path(BASE, "decontx", sprintf("decontx_sample_%d.rds", sid))
  if (!file.exists(decon_path)) next
  decon <- tryCatch(readRDS(decon_path), error = function(e) NULL)
  if (is.null(decon)) next
  corrected_counts <- decon$counts

  # Match barcodes
  common_bc <- intersect(eos_in_sample$raw_barcode,
                         intersect(colnames(orig_counts), colnames(corrected_counts)))
  if (length(common_bc) == 0) next

  # Check Ig genes
  for (gene in ig_present) {
    if (gene %in% rownames(orig_counts) && gene %in% rownames(corrected_counts)) {
      before_vals <- as.numeric(orig_counts[gene, common_bc])
      after_vals  <- as.numeric(corrected_counts[gene, common_bc])

      before_after_data[[length(before_after_data) + 1]] <- data.frame(
        gene = gene,
        condition = cond,
        sample_id = sid,
        stage = "Before DecontX",
        expression = before_vals,
        stringsAsFactors = FALSE
      )
      before_after_data[[length(before_after_data) + 1]] <- data.frame(
        gene = gene,
        condition = cond,
        sample_id = sid,
        stage = "After DecontX",
        expression = after_vals,
        stringsAsFactors = FALSE
      )
    }
  }
}

if (length(before_after_data) > 0) {
  ba_df <- do.call(rbind, before_after_data)
  ba_df$stage <- factor(ba_df$stage, levels = c("Before DecontX", "After DecontX"))
  ba_df$condition_short <- sub("_\\d+$", "", ba_df$condition)

  # Summary: pct expressing before vs after
  cat("  Igha expression before vs after DecontX (per condition):\n")
  igha_summary <- ba_df %>%
    filter(gene == "Igha") %>%
    group_by(condition_short, stage) %>%
    summarise(
      n_cells = n(),
      pct_pos = round(100 * mean(expression > 0), 2),
      mean_expr = round(mean(expression), 4),
      .groups = "drop"
    ) %>%
    arrange(condition_short, stage)
  print(as.data.frame(igha_summary))

  # Overall before vs after
  cat("\n  Overall Ig gene reduction:\n")
  overall_ba <- ba_df %>%
    group_by(gene, stage) %>%
    summarise(
      pct_pos = round(100 * mean(expression > 0), 2),
      mean_expr = round(mean(expression), 4),
      .groups = "drop"
    ) %>%
    pivot_wider(names_from = stage,
                values_from = c(pct_pos, mean_expr))
  print(as.data.frame(overall_ba))

  # ── Figure 1: Before vs After violin for Ig genes ──────────────────────────
  cat("\n  Plotting before vs after Ig gene violins...\n")

  p_ba <- ggplot(ba_df, aes(x = stage, y = expression, fill = stage)) +
    geom_violin(scale = "width", alpha = 0.7) +
    geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white", alpha = 0.5) +
    facet_wrap(~ gene, scales = "free_y", nrow = 1) +
    scale_fill_manual(values = c("Before DecontX" = "#E84A5F",
                                 "After DecontX" = "#5BA4A4")) +
    labs(x = "", y = "Raw Counts",
         title = "Ig Gene Expression in Eosinophils: Before vs After DecontX",
         subtitle = "DecontX corrected counts should show substantial reduction") +
    theme_pub(base_size = 11) +
    theme(axis.text.x = element_text(angle = 30, hjust = 1),
          legend.position = "none")

  pdf(file.path(fig_dir, "adhoc_igha_before_after_violin.pdf"), width = 14, height = 5)
  print(p_ba)
  dev.off()
  cat("  Saved: adhoc_igha_before_after_violin.pdf\n")

  # ── Figure 2: Before vs After by condition (Igha only) ─────────────────────
  igha_ba <- ba_df %>% filter(gene == "Igha")

  if (nrow(igha_ba) > 0) {
    p_ba_cond <- ggplot(igha_ba, aes(x = condition_short, y = expression, fill = stage)) +
      geom_boxplot(outlier.size = 0.3, alpha = 0.7, position = position_dodge(0.8)) +
      scale_fill_manual(values = c("Before DecontX" = "#E84A5F",
                                   "After DecontX" = "#5BA4A4")) +
      labs(x = "Condition", y = "Raw Counts", fill = "",
           title = "Igha Expression in Eosinophils by Condition",
           subtitle = "Before vs After DecontX") +
      theme_pub(base_size = 11) +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    pdf(file.path(fig_dir, "adhoc_igha_before_after_by_condition.pdf"), width = 10, height = 6)
    print(p_ba_cond)
    dev.off()
    cat("  Saved: adhoc_igha_before_after_by_condition.pdf\n")
  }
} else {
  cat("  Could not load before/after data for comparison.\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# 4. DecontX contamination score distribution
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========================================\n")
cat("Step 4: DecontX contamination score distribution\n")
cat("========================================\n\n")

# Load contamination scores
scores_path <- file.path(BASE, "decontx", "decontx_all_scores.csv")
if (file.exists(scores_path)) {
  contam_scores <- read.csv(scores_path, stringsAsFactors = FALSE)

  # Match to eos cells
  contam_scores <- contam_scores %>%
    left_join(sample_info, by = "sample_id") %>%
    rename(orig_condition = condition.x, sample_condition = condition.y)
  contam_scores$cell_id <- paste0(contam_scores$sample_condition, "_",
                                  contam_scores$barcode)

  # Match to eos object
  eos_contam <- contam_scores %>%
    filter(cell_id %in% rownames(eos_meta))

  if (nrow(eos_contam) > 0) {
    eos_contam <- eos_contam %>%
      left_join(eos_meta %>% select(cell_id, eos_cluster, condition),
                by = "cell_id")

    cat(sprintf("  Matched %d eos cells to contamination scores\n", nrow(eos_contam)))

    # Summary by cluster
    cat("\n  Contamination by eos cluster:\n")
    cluster_contam <- eos_contam %>%
      group_by(eos_cluster) %>%
      summarise(
        n = n(),
        mean_contam = round(mean(contamination), 4),
        median_contam = round(median(contamination), 4),
        pct_above_0.2 = round(100 * mean(contamination > 0.2), 1),
        .groups = "drop"
      )
    print(as.data.frame(cluster_contam))

    # Summary by condition
    cat("\n  Contamination by condition:\n")
    cond_contam <- eos_contam %>%
      group_by(condition) %>%
      summarise(
        n = n(),
        mean_contam = round(mean(contamination), 4),
        median_contam = round(median(contamination), 4),
        .groups = "drop"
      )
    print(as.data.frame(cond_contam))

    # ── Figure 3: Violin of contamination score by cluster ───────────────────
    p_contam_cluster <- ggplot(eos_contam,
                                aes(x = eos_cluster, y = contamination,
                                    fill = eos_cluster)) +
      geom_violin(scale = "width", alpha = 0.7) +
      geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white", alpha = 0.5) +
      scale_fill_manual(values = eos_colors) +
      geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
      labs(x = "Eos Cluster", y = "DecontX Contamination Score",
           title = "Contamination Score Distribution by Eos Cluster") +
      theme_pub() + theme(legend.position = "none")

    # ── Figure 4: Violin by condition ────────────────────────────────────────
    p_contam_cond <- ggplot(eos_contam,
                             aes(x = condition, y = contamination,
                                 fill = condition)) +
      geom_violin(scale = "width", alpha = 0.7) +
      geom_boxplot(width = 0.15, outlier.size = 0.3, fill = "white", alpha = 0.5) +
      scale_fill_manual(values = cond_colors) +
      geom_hline(yintercept = 0.3, linetype = "dashed", color = "red") +
      labs(x = "Condition", y = "DecontX Contamination Score",
           title = "Contamination Score by Condition (Eos only)") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")

    pdf(file.path(fig_dir, "adhoc_contamination_violin_cluster.pdf"), width = 8, height = 5)
    print(p_contam_cluster)
    dev.off()

    pdf(file.path(fig_dir, "adhoc_contamination_violin_condition.pdf"), width = 10, height = 5)
    print(p_contam_cond)
    dev.off()
    cat("  Saved: adhoc_contamination_violin_cluster.pdf\n")
    cat("  Saved: adhoc_contamination_violin_condition.pdf\n")

    # ── Figure 5: UMAP colored by contamination score ────────────────────────
    # Add contamination score to eos metadata for UMAP
    contam_vec <- setNames(eos_contam$contamination, eos_contam$cell_id)
    eos$decontx_score <- contam_vec[rownames(eos@meta.data)]

    p_umap_contam <- FeaturePlot(eos, features = "decontx_score", pt.size = 0.3) +
      scale_color_viridis_c(name = "Contamination") +
      ggtitle("UMAP: DecontX Contamination Score") +
      theme_umap()

    pdf(file.path(fig_dir, "adhoc_umap_contamination_score.pdf"), width = 8, height = 6)
    print(p_umap_contam)
    dev.off()
    cat("  Saved: adhoc_umap_contamination_score.pdf\n")
  } else {
    cat("  No eos cells matched to contamination scores.\n")
  }
} else {
  cat("  Contamination scores file not found.\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# 5. Multi-Ig gene co-expression in eosinophils (residual contamination check)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========================================\n")
cat("Step 5: Multi-Ig gene co-expression check\n")
cat("========================================\n\n")

if (length(ig_present) >= 2) {
  counts_mat <- GetAssayData(eos, layer = "counts")

  # For each cell, count how many Ig genes are expressed
  ig_expr_mat <- as.matrix(counts_mat[ig_present, , drop = FALSE])
  n_ig_per_cell <- colSums(ig_expr_mat > 0)

  cat("  Number of Ig genes co-expressed per eosinophil:\n")
  print(table(n_ig_per_cell))
  cat(sprintf("  Cells with 0 Ig genes: %d (%.1f%%)\n",
              sum(n_ig_per_cell == 0), 100 * mean(n_ig_per_cell == 0)))
  cat(sprintf("  Cells with 1 Ig gene:  %d (%.1f%%)\n",
              sum(n_ig_per_cell == 1), 100 * mean(n_ig_per_cell == 1)))
  cat(sprintf("  Cells with 2+ Ig genes: %d (%.1f%%)\n",
              sum(n_ig_per_cell >= 2), 100 * mean(n_ig_per_cell >= 2)))

  # If multiple Ig genes co-express at low levels -> ambient contamination signature
  multi_ig_cells <- names(n_ig_per_cell)[n_ig_per_cell >= 2]
  if (length(multi_ig_cells) > 0) {
    cat(sprintf("\n  Multi-Ig+ cells (n=%d) - expression profile:\n", length(multi_ig_cells)))
    multi_ig_expr <- ig_expr_mat[, multi_ig_cells, drop = FALSE]
    cat("  Mean expression per Ig gene in multi-Ig+ cells:\n")
    print(round(rowMeans(multi_ig_expr), 3))
    cat("  -> Low-level co-expression of multiple Ig genes = ambient contamination\n")
  }

  # ── Figure 6: Ig co-expression dotplot ─────────────────────────────────────
  if (length(ig_present) >= 2) {
    p_ig_dot <- DotPlot(eos, features = ig_present, group.by = "eos_cluster") +
      RotatedAxis() +
      ggtitle("Ig Gene Expression by Eos Cluster (post-DecontX)") +
      theme_pub()

    pdf(file.path(fig_dir, "adhoc_ig_dotplot_by_cluster.pdf"), width = 10, height = 5)
    print(p_ig_dot)
    dev.off()
    cat("  Saved: adhoc_ig_dotplot_by_cluster.pdf\n")
  }

  # ── Figure 7: UMAP feature plots for Ig genes ─────────────────────────────
  ig_for_umap <- ig_present[ig_present %in% rownames(eos)]
  if (length(ig_for_umap) >= 2) {
    p_ig_umaps <- FeaturePlot(eos, features = ig_for_umap[1:min(6, length(ig_for_umap))],
                               pt.size = 0.2, ncol = 3) &
      theme_umap() & scale_color_viridis_c()

    pdf(file.path(fig_dir, "adhoc_ig_umap_featureplots.pdf"),
        width = 12, height = max(4, 4 * ceiling(min(6, length(ig_for_umap)) / 3)))
    print(p_ig_umaps)
    dev.off()
    cat("  Saved: adhoc_ig_umap_featureplots.pdf\n")
  }
}

# ══════════════════════════════════════════════════════════════════════════════
# 6. Igha+ eosinophils: do they retain eos markers? Or show B cell markers?
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========================================\n")
cat("Step 6: Igha+ vs Igha- eosinophil identity check\n")
cat("========================================\n\n")

if ("Igha" %in% ig_present) {
  counts_mat <- GetAssayData(eos, layer = "counts")
  igha_status <- ifelse(counts_mat["Igha", ] > 0, "Igha+", "Igha-")
  eos$igha_status <- igha_status

  cat(sprintf("  Igha+ cells: %d, Igha- cells: %d\n",
              sum(igha_status == "Igha+"), sum(igha_status == "Igha-")))

  # Compare eos marker expression
  check_genes <- c(eos_markers_present, bcell_present)
  if (length(check_genes) > 0) {
    norm_mat <- GetAssayData(eos, layer = "data")

    cat("\n  Mean normalized expression (Igha+ vs Igha-):\n")
    cat(sprintf("  %-10s  %10s  %10s  %s\n", "Gene", "Igha+", "Igha-", "Interpretation"))
    cat(paste(rep("-", 55), collapse = ""), "\n")

    for (gene in check_genes) {
      if (gene %in% rownames(norm_mat)) {
        pos_mean <- mean(norm_mat[gene, igha_status == "Igha+"])
        neg_mean <- mean(norm_mat[gene, igha_status == "Igha-"])

        interpretation <- ""
        if (gene %in% eos_markers_present) {
          if (pos_mean > 0.5 * neg_mean) {
            interpretation <- "Eos marker retained -> residual contamination"
          } else {
            interpretation <- "Eos marker LOST -> possible doublet"
          }
        } else if (gene %in% bcell_present) {
          if (pos_mean > 0.1) {
            interpretation <- "B cell marker present -> possible doublet"
          } else {
            interpretation <- "B cell marker absent -> not doublet"
          }
        }

        cat(sprintf("  %-10s  %10.4f  %10.4f  %s\n",
                    gene, pos_mean, neg_mean, interpretation))
      }
    }

    # ── Figure 8: Violin comparison Igha+ vs Igha- ───────────────────────────
    genes_to_plot <- c(check_genes[check_genes %in% rownames(eos)])
    if (length(genes_to_plot) > 0) {
      plot_data <- data.frame()
      for (gene in genes_to_plot) {
        df <- data.frame(
          gene = gene,
          expression = as.numeric(norm_mat[gene, ]),
          igha_status = igha_status,
          stringsAsFactors = FALSE
        )
        plot_data <- rbind(plot_data, df)
      }

      p_identity <- ggplot(plot_data, aes(x = igha_status, y = expression,
                                           fill = igha_status)) +
        geom_violin(scale = "width", alpha = 0.7) +
        geom_boxplot(width = 0.15, outlier.size = 0.2, fill = "white", alpha = 0.5) +
        facet_wrap(~ gene, scales = "free_y", ncol = 4) +
        scale_fill_manual(values = c("Igha+" = "#E84A5F", "Igha-" = "#6CA6CD")) +
        labs(x = "", y = "Normalized Expression",
             title = "Eos & B Cell Marker Expression: Igha+ vs Igha- Eosinophils",
             subtitle = "If eos markers retained + B cell markers absent -> ambient contamination, not doublets") +
        theme_pub(base_size = 10) +
        theme(legend.position = "none")

      pdf(file.path(fig_dir, "adhoc_igha_pos_vs_neg_markers.pdf"),
          width = 14, height = max(4, 3 * ceiling(length(genes_to_plot) / 4)))
      print(p_identity)
      dev.off()
      cat("\n  Saved: adhoc_igha_pos_vs_neg_markers.pdf\n")
    }
  }

  # Cluster distribution of Igha+ cells
  cat("\n  Igha+ cells per cluster:\n")
  igha_cluster <- table(eos$eos_cluster, eos$igha_status)
  igha_pct <- round(100 * prop.table(igha_cluster, margin = 1), 1)
  print(igha_cluster)
  cat("\n  % Igha+ per cluster:\n")
  print(igha_pct)

  # ── Figure 9: UMAP highlighting Igha+ cells ───────────────────────────────
  p_umap_igha <- DimPlot(eos, group.by = "igha_status", pt.size = 0.3,
                          cols = c("Igha+" = "#E84A5F", "Igha-" = "grey85"),
                          order = c("Igha+")) +
    ggtitle("UMAP: Igha+ Eosinophils (post-DecontX)") +
    theme_umap()

  pdf(file.path(fig_dir, "adhoc_umap_igha_positive.pdf"), width = 8, height = 6)
  print(p_umap_igha)
  dev.off()
  cat("  Saved: adhoc_umap_igha_positive.pdf\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# Summary
# ══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat("================================================================\n")
cat("SUMMARY\n")
cat("================================================================\n\n")

cat("Key questions answered:\n\n")

cat("1. Did DecontX reduce Igha?\n")
if (exists("overall_ba")) {
  cat("   YES - see before/after violin plots\n\n")
} else {
  cat("   Check adhoc_igha_before_after_violin.pdf\n\n")
}

cat("2. Contamination score distribution?\n")
cat("   See violin plots by cluster and condition\n")
cat("   Any condition-specific enrichment indicates technical artifact\n\n")

cat("3. Do other Ig genes co-express with Igha?\n")
if (exists("n_ig_per_cell")) {
  cat(sprintf("   Cells with 2+ Ig genes: %d (%.1f%%)\n",
              sum(n_ig_per_cell >= 2), 100 * mean(n_ig_per_cell >= 2)))
  cat("   Low-level multi-Ig co-expression = ambient contamination signature\n\n")
}

cat("4. Do Igha+ cells retain eosinophil identity?\n")
cat("   Check adhoc_igha_pos_vs_neg_markers.pdf\n")
cat("   If eos markers strong + B cell markers absent -> residual contamination\n")
cat("   If eos markers lost + B cell markers present -> doublet concern\n\n")

cat("Output directory:", fig_dir, "\n")
cat("\nDone! adhoc_igha_contamination_check.R completed.\n")
