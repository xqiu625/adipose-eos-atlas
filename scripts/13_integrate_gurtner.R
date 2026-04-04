#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --job-name=integrate_gurtner
#SBATCH --output=integrate_gurtner_%j.out

# =============================================================================
# 13_integrate_gurtner.R
# Harmony integration of our adipose eosinophils with Gurtner reference
# (GSE182001, blood & bone marrow eosinophils)
# Input:  $BASE/eos_subclustered.rds, gurtner_eos.rds (candidate paths)
# Output: qiu_gurtner_integrated.rds, PDFs
# =============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(harmony)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
outdir <- file.path(BASE, "integration")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Load our eosinophil data ------------------------------------------------
cat("Loading our eosinophil data ...\n")
eos <- tryCatch(
  readRDS(file.path(BASE, "eos_subclustered.rds")),
  error = function(e) { cat("ERROR: Cannot load eos_subclustered.rds:", e$message, "\n"); quit(status = 1) }
)
# Remove C6 if still present (too few cells)
if ("C6" %in% levels(Idents(eos))) {
  eos <- subset(eos, idents = "C6", invert = TRUE)
  Idents(eos) <- droplevels(Idents(eos))
  cat(sprintf("  Removed C6. Remaining: %d eos (C0-C5)\n", ncol(eos)))
}
cat("  Our cells:", ncol(eos), "\n")

# --- Load Gurtner reference data ---------------------------------------------
cat("Searching for Gurtner reference data (GSE182001) ...\n")
gurtner_candidates <- c(
  "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/Gurtner_eosinophil_subset.rds",
  "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/GSE182001_eosinophil_subset.rds",
  "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/Gurtner_merged_tissues_seurat_object.rds"
)

gurtner <- NULL
for (path in gurtner_candidates) {
  cat("  Trying:", path, "\n")
  if (file.exists(path)) {
    cat("  Found! Loading ...\n")
    gurtner <- tryCatch(
      readRDS(path),
      error = function(e) { cat("  ERROR loading:", e$message, "\n"); NULL }
    )
    if (!is.null(gurtner)) break
  }
}

if (is.null(gurtner)) {
  cat("ERROR: Gurtner reference data not found at any candidate path.\n")
  cat("  Searched:\n")
  for (p in gurtner_candidates) cat("    ", p, "\n")
  quit(status = 1)
}
cat("  Gurtner cells:", ncol(gurtner), "\n")

# --- Add dataset labels -------------------------------------------------------
cat("Adding dataset labels ...\n")
eos$dataset     <- "Qiu_2025_Adipose"
gurtner$dataset <- "Gurtner_2022"

# Preserve original cluster identities
eos$original_cluster <- as.character(Idents(eos))
if ("seurat_clusters" %in% colnames(gurtner@meta.data)) {
  gurtner$original_cluster <- as.character(gurtner$seurat_clusters)
} else {
  gurtner$original_cluster <- as.character(Idents(gurtner))
}

# Preserve organ labels from Gurtner (Bone Marrow, Blood, Colon, etc.)
eos$organ <- "Adipose"
if (!"organ" %in% colnames(gurtner@meta.data)) {
  gurtner$organ <- "Unknown"
}
cat("  Gurtner organ distribution:\n")
print(table(gurtner$organ))

# --- Merge objects ------------------------------------------------------------
cat("Merging datasets ...\n")
combined <- merge(eos, y = gurtner, add.cell.ids = c("Qiu", "Gurtner"),
                  project = "Qiu_Gurtner_Integration")
cat("  Combined cells:", ncol(combined), "\n")

# --- Standard Seurat workflow -------------------------------------------------
cat("Normalizing ...\n")
combined <- NormalizeData(combined, verbose = FALSE)

cat("Finding variable features ...\n")
combined <- FindVariableFeatures(combined, selection.method = "vst",
                                  nfeatures = 2000, verbose = FALSE)

cat("Scaling data ...\n")
combined <- ScaleData(combined, verbose = FALSE)

cat("Running PCA ...\n")
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

# --- Harmony integration ------------------------------------------------------
cat("Running Harmony integration (group.by = dataset) ...\n")
combined <- RunHarmony(combined, group.by.vars = "dataset",
                        reduction.use = "pca", verbose = TRUE)

# --- UMAP on Harmony embeddings ----------------------------------------------
cat("Running UMAP on Harmony embeddings ...\n")
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30, verbose = FALSE)

# --- Clustering on integrated space ------------------------------------------
cat("Finding neighbors and clusters on Harmony space ...\n")
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, verbose = FALSE)
cat("  Integrated clusters:", length(unique(Idents(combined))), "\n")

# --- Figures ------------------------------------------------------------------
cat("Generating figures ...\n")

# Figure 1: UMAP colored by dataset (combined + split)
pdf(file.path(outdir, "gurtner_integration_umap_dataset.pdf"), width = 14, height = 6)
tryCatch({
  p1 <- DimPlot(combined, reduction = "umap", group.by = "dataset",
                 cols = c("Qiu_2025_Adipose" = "#E63946",
                          "Gurtner_2022" = "#1D3557"),
                 pt.size = 0.3) +
    ggtitle("Combined: Adipose vs Blood/BM eosinophils") +
    theme_umap()

  p2 <- DimPlot(combined, reduction = "umap", group.by = "dataset",
                 split.by = "dataset",
                 cols = c("Qiu_2025_Adipose" = "#E63946",
                          "Gurtner_2022" = "#1D3557"),
                 pt.size = 0.3) +
    ggtitle("Split by dataset") +
    theme_umap()

  print(p1 + p2)
}, error = function(e) cat("  Warning: dataset UMAP failed:", e$message, "\n"))
dev.off()

# Figure 2: UMAP colored by integrated cluster
pdf(file.path(outdir, "gurtner_integration_umap_cluster.pdf"), width = 8, height = 6)
tryCatch({
  p <- DimPlot(combined, reduction = "umap", label = TRUE, label.size = 4,
               pt.size = 0.3) +
    ggtitle("Integrated clusters (Qiu + Gurtner)") +
    theme_umap()
  print(p)
}, error = function(e) cat("  Warning: cluster UMAP failed:", e$message, "\n"))
dev.off()

# Figure 3: UMAP colored by original cluster identity per dataset
pdf(file.path(outdir, "gurtner_integration_umap_original_clusters.pdf"), width = 14, height = 6)
tryCatch({
  p1 <- DimPlot(subset(combined, dataset == "Qiu_2025_Adipose"),
                 reduction = "umap", group.by = "original_cluster",
                 cols = eos_colors,
                 label = TRUE, pt.size = 0.3) +
    ggtitle("Qiu_2025 adipose eos clusters") +
    theme_umap()

  p2 <- DimPlot(subset(combined, dataset == "Gurtner_2022"),
                 reduction = "umap", group.by = "original_cluster",
                 label = TRUE, pt.size = 0.3) +
    ggtitle("Gurtner blood/BM eos clusters") +
    theme_umap()

  print(p1 + p2)
}, error = function(e) cat("  Warning: original cluster UMAP failed:", e$message, "\n"))
dev.off()

# Figure 4: Cluster comparison - confusion matrix
pdf(file.path(outdir, "gurtner_integration_confusion_matrix.pdf"), width = 10, height = 8)
tryCatch({
  conf_df <- combined@meta.data %>%
    select(dataset, integrated_cluster = seurat_clusters, original_cluster) %>%
    group_by(dataset, integrated_cluster, original_cluster) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(dataset, integrated_cluster) %>%
    mutate(frac = n / sum(n))

  p <- ggplot(conf_df, aes(x = integrated_cluster, y = original_cluster, fill = frac)) +
    geom_tile() +
    geom_text(aes(label = n), size = 2.5) +
    scale_fill_gradient(low = "white", high = "#264653", name = "Fraction") +
    facet_wrap(~dataset, scales = "free_y") +
    labs(x = "Integrated cluster", y = "Original cluster",
         title = "Cluster correspondence: Adipose vs Blood/BM eos") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
}, error = function(e) cat("  Warning: confusion matrix failed:", e$message, "\n"))
dev.off()

# Figure 5: Composition comparison
pdf(file.path(outdir, "gurtner_integration_composition.pdf"), width = 10, height = 6)
tryCatch({
  comp_df <- combined@meta.data %>%
    group_by(seurat_clusters, dataset) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(frac = n / sum(n))

  p <- ggplot(comp_df, aes(x = seurat_clusters, y = frac, fill = dataset)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Qiu_2025_Adipose" = "#E63946",
                                  "Gurtner_2022" = "#1D3557")) +
    labs(x = "Integrated cluster", y = "Fraction",
         title = "Dataset composition per integrated cluster") +
    theme_pub()
  print(p)
}, error = function(e) cat("  Warning: composition plot failed:", e$message, "\n"))
dev.off()

# Figure 6: Tissue-specific marker expression across integrated clusters
pdf(file.path(outdir, "gurtner_integration_markers.pdf"), width = 12, height = 8)
tryCatch({
  # Key eosinophil markers to compare across tissue sources
  eos_markers <- c("Siglecf", "Epx", "Prg2", "Ccr3", "Il5ra",
                    "Cd63", "Cd9", "Il4", "Il13", "Ccl6")
  # Keep only markers present in the data
  eos_markers <- eos_markers[eos_markers %in% rownames(combined)]

  if (length(eos_markers) >= 3) {
    p <- DotPlot(combined, features = eos_markers, group.by = "seurat_clusters",
                  split.by = "dataset") +
      RotatedAxis() +
      ggtitle("Marker expression: Adipose vs Blood/BM eos") +
      theme_pub()
    print(p)
  } else {
    cat("  Fewer than 3 markers found in data - skipping marker dotplot\n")
  }
}, error = function(e) cat("  Warning: marker dotplot failed:", e$message, "\n"))
dev.off()

# --- Save integrated object ---------------------------------------------------
cat("Saving integrated object ...\n")
saveRDS(combined, file.path(outdir, "qiu_gurtner_integrated.rds"))
cat("  Saved qiu_gurtner_integrated.rds\n")

# --- Export for CellRank (Python) ---------------------------------------------
cat("Exporting data for CellRank ...\n")
export_dir <- file.path(outdir, "for_cellrank")
dir.create(export_dir, showWarnings = FALSE, recursive = TRUE)

# Join layers for Seurat v5 (merged object has split layers)
combined <- JoinLayers(combined)

# Counts matrix (sparse)
counts_mat <- GetAssayData(combined, layer = "counts")
Matrix::writeMM(counts_mat, file.path(export_dir, "counts.mtx"))

# Gene names and cell barcodes
writeLines(rownames(combined), file.path(export_dir, "genes.txt"))
writeLines(colnames(combined), file.path(export_dir, "barcodes.txt"))

# Metadata
write.csv(combined@meta.data, file.path(export_dir, "metadata.csv"))

# Harmony embeddings
write.csv(Embeddings(combined, "harmony"), file.path(export_dir, "harmony_embeddings.csv"))

# UMAP embeddings
write.csv(Embeddings(combined, "umap"), file.path(export_dir, "umap_embeddings.csv"))

cat("  Exported to:", export_dir, "\n")

cat("\n========== Gurtner integration complete ==========\n")
cat("Output directory:", outdir, "\n")
sessionInfo()
