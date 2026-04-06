#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 30_cross_dataset_similarity.R
# Cross-dataset cluster similarity: Qiu C0-C5 vs Hu clusters vs Gurtner organs
# Method: Spearman correlation of average gene expression profiles
# (Similar to approach in Hu et al. 2026)
# Input:  qiu_hu_integrated.rds, qiu_gurtner_integrated.rds
# Output: correlation heatmaps + CSVs
# SLURM:  batch, 64G, 4 CPUs, 2h
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(patchwork)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
outdir <- file.path(BASE, "integration", "similarity")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# ==============================================================================
# Helper: compute cluster-to-cluster Spearman correlation
# ==============================================================================
compute_cluster_correlation <- function(obj, group_col, groups_a, groups_b,
                                        n_hvg = 2000) {
  # Use variable features or find them
  if (length(VariableFeatures(obj)) > 0) {
    hvg <- VariableFeatures(obj)
  } else {
    obj <- FindVariableFeatures(obj, nfeatures = n_hvg, verbose = FALSE)
    hvg <- VariableFeatures(obj)
  }
  hvg <- hvg[hvg %in% rownames(obj)]

  # Average expression per group
  Idents(obj) <- group_col
  avg <- AverageExpression(obj, features = hvg, assays = "RNA", layer = "data")$RNA
  avg <- as.matrix(avg)

  # Keep only requested groups
  all_groups <- c(groups_a, groups_b)
  keep <- colnames(avg) %in% all_groups
  avg <- avg[, keep, drop = FALSE]

  # Spearman correlation
  cor_mat <- cor(avg, method = "spearman")
  return(list(cor = cor_mat, avg = avg))
}

# ==============================================================================
# 1. Qiu C0-C5 vs Hu clusters
# ==============================================================================
cat("=== Qiu vs Hu cluster similarity ===\n")

hu_int_file <- file.path(BASE, "integration", "qiu_hu_integrated.rds")
if (file.exists(hu_int_file)) {
  cat("Loading Qiu-Hu integrated object...\n")
  hu_int <- readRDS(hu_int_file)
  if (inherits(hu_int[["RNA"]], "Assay5")) hu_int[["RNA"]] <- JoinLayers(hu_int[["RNA"]])

  # Create dataset_cluster label
  hu_int$dataset_cluster <- ifelse(
    hu_int$dataset == "Qiu_2025",
    paste0("Qiu_", hu_int$original_cluster),
    paste0("Hu_", hu_int$original_cluster)
  )

  qiu_clusters <- paste0("Qiu_C", 0:5)
  hu_clusters <- sort(unique(hu_int$dataset_cluster[hu_int$dataset == "Hu_Reference"]))
  cat("  Qiu clusters:", paste(qiu_clusters, collapse = ", "), "\n")
  cat("  Hu clusters:", paste(hu_clusters, collapse = ", "), "\n")

  res_hu <- compute_cluster_correlation(hu_int, "dataset_cluster",
                                         qiu_clusters, hu_clusters)

  # Save correlation matrix
  write.csv(res_hu$cor, file.path(outdir, "qiu_hu_spearman_correlation.csv"))

  # Heatmap: Qiu rows vs Hu columns
  qiu_in <- intersect(qiu_clusters, rownames(res_hu$cor))
  hu_in <- intersect(hu_clusters, colnames(res_hu$cor))

  if (length(qiu_in) > 0 && length(hu_in) > 0) {
    sub_cor <- res_hu$cor[qiu_in, hu_in, drop = FALSE]

    pdf(file.path(outdir, "qiu_hu_similarity_heatmap.pdf"), width = 10, height = 6)
    pheatmap(sub_cor,
             color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
             breaks = seq(0.5, 1, length.out = 101),
             cluster_rows = TRUE, cluster_cols = TRUE,
             display_numbers = TRUE, number_format = "%.2f",
             fontsize_number = 8,
             fontsize_row = 10, fontsize_col = 10,
             main = "Spearman Correlation: Qiu C0-C5 vs Hu Clusters",
             border_color = NA)
    dev.off()
    cat("  Saved qiu_hu_similarity_heatmap.pdf\n")
  }

  rm(hu_int); gc()
} else {
  cat("  WARNING: qiu_hu_integrated.rds not found, skipping\n")
}

# ==============================================================================
# 2. Qiu C0-C5 vs Gurtner BM / Blood
# ==============================================================================
cat("\n=== Qiu vs Gurtner cluster similarity ===\n")

gurtner_int_file <- file.path(BASE, "integration", "qiu_gurtner_integrated.rds")
if (file.exists(gurtner_int_file)) {
  cat("Loading Qiu-Gurtner integrated object...\n")
  gurtner_int <- readRDS(gurtner_int_file)
  if (inherits(gurtner_int[["RNA"]], "Assay5")) gurtner_int[["RNA"]] <- JoinLayers(gurtner_int[["RNA"]])

  # Create dataset_organ label for Gurtner, dataset_cluster for Qiu
  gurtner_int$similarity_group <- case_when(
    gurtner_int$dataset == "Qiu_2025_Adipose" ~ paste0("Qiu_", gurtner_int$original_cluster),
    gurtner_int$organ == "Bone Marrow" ~ "Gurtner_BM",
    gurtner_int$organ == "Blood" ~ "Gurtner_Blood",
    TRUE ~ paste0("Gurtner_", gurtner_int$organ)
  )

  qiu_clusters <- paste0("Qiu_C", 0:5)
  gurtner_groups <- sort(unique(gurtner_int$similarity_group[
    gurtner_int$dataset != "Qiu_2025_Adipose"]))
  cat("  Qiu clusters:", paste(qiu_clusters, collapse = ", "), "\n")
  cat("  Gurtner groups:", paste(gurtner_groups, collapse = ", "), "\n")

  res_gurtner <- compute_cluster_correlation(gurtner_int, "similarity_group",
                                              qiu_clusters, gurtner_groups)

  write.csv(res_gurtner$cor, file.path(outdir, "qiu_gurtner_spearman_correlation.csv"))

  qiu_in <- intersect(qiu_clusters, rownames(res_gurtner$cor))
  gurtner_in <- intersect(gurtner_groups, colnames(res_gurtner$cor))

  if (length(qiu_in) > 0 && length(gurtner_in) > 0) {
    sub_cor <- res_gurtner$cor[qiu_in, gurtner_in, drop = FALSE]

    pdf(file.path(outdir, "qiu_gurtner_similarity_heatmap.pdf"), width = 10, height = 6)
    pheatmap(sub_cor,
             color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
             breaks = seq(0.5, 1, length.out = 101),
             cluster_rows = TRUE, cluster_cols = TRUE,
             display_numbers = TRUE, number_format = "%.2f",
             fontsize_number = 8,
             fontsize_row = 10, fontsize_col = 10,
             main = "Spearman Correlation: Qiu C0-C5 vs Gurtner Tissues",
             border_color = NA)
    dev.off()
    cat("  Saved qiu_gurtner_similarity_heatmap.pdf\n")
  }

  rm(gurtner_int); gc()
} else {
  cat("  WARNING: qiu_gurtner_integrated.rds not found, skipping\n")
}

# ==============================================================================
# 3. Combined: all datasets
# ==============================================================================
cat("\n=== Three-way similarity (if three-way integration exists) ===\n")

three_way_file <- file.path(BASE, "integration", "three_way", "three_way_integrated.rds")
if (file.exists(three_way_file)) {
  cat("Loading three-way integrated object...\n")
  tw <- readRDS(three_way_file)
  if (inherits(tw[["RNA"]], "Assay5")) tw[["RNA"]] <- JoinLayers(tw[["RNA"]])

  # Create unified group labels
  tw$similarity_group <- case_when(
    tw$dataset == "Qiu_2025_Adipose" ~ paste0("Qiu_", tw$original_cluster),
    tw$dataset == "Hu_2026_Healthy" ~ paste0("Hu_", tw$original_cluster),
    tw$organ == "Bone Marrow" ~ "Gurtner_BM",
    tw$organ == "Blood" ~ "Gurtner_Blood",
    TRUE ~ paste0("Gurtner_", tw$organ)
  )

  all_groups <- sort(unique(tw$similarity_group))
  cat("  All groups:", paste(all_groups, collapse = ", "), "\n")

  # Filter to groups with >= 30 cells
  group_sizes <- table(tw$similarity_group)
  keep_groups <- names(group_sizes[group_sizes >= 30])
  cat("  Groups with >= 30 cells:", length(keep_groups), "\n")

  tw_sub <- subset(tw, cells = colnames(tw)[tw$similarity_group %in% keep_groups])
  tw_sub$similarity_group <- factor(tw_sub$similarity_group)

  res_tw <- compute_cluster_correlation(tw_sub, "similarity_group",
                                         keep_groups, keep_groups)

  write.csv(res_tw$cor, file.path(outdir, "three_way_spearman_correlation.csv"))

  pdf(file.path(outdir, "three_way_similarity_heatmap.pdf"), width = 12, height = 10)
  # Annotation for dataset source
  ann_row <- data.frame(
    Dataset = case_when(
      grepl("^Qiu_", rownames(res_tw$cor)) ~ "Qiu",
      grepl("^Hu_", rownames(res_tw$cor)) ~ "Hu",
      TRUE ~ "Gurtner"
    ),
    row.names = rownames(res_tw$cor)
  )
  ann_colors <- list(Dataset = c("Qiu" = "#E63946", "Hu" = "#4DAF4A", "Gurtner" = "#1D3557"))

  pheatmap(res_tw$cor,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           breaks = seq(0.5, 1, length.out = 101),
           cluster_rows = TRUE, cluster_cols = TRUE,
           annotation_row = ann_row, annotation_col = ann_row,
           annotation_colors = ann_colors,
           display_numbers = TRUE, number_format = "%.2f",
           fontsize_number = 6,
           fontsize_row = 9, fontsize_col = 9,
           main = "Cross-Dataset Cluster Similarity (Spearman)",
           border_color = NA)
  dev.off()
  cat("  Saved three_way_similarity_heatmap.pdf\n")

  rm(tw, tw_sub); gc()
} else {
  cat("  Three-way integrated RDS not found, skipping\n")
}

cat("\n========== Cross-dataset similarity analysis complete ==========\n")
cat("Output directory:", outdir, "\n")
