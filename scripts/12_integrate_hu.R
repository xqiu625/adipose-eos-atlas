#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --job-name=integrate_hu
#SBATCH --output=integrate_hu_%j.out

# =============================================================================
# 12_integrate_hu.R
# Harmony integration of our adipose eosinophils with Hu reference dataset
# Input:  $BASE/eos_subclustered.rds, hu_adipose_eos.rds (candidate paths)
# Output: qiu_hu_integrated.rds, PDFs
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
cat("  Our cells:", ncol(eos), "\n")

# --- Load Hu reference data --------------------------------------------------
cat("Searching for Hu reference data ...\n")
hu_candidates <- c(
  "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/hu_2026/hu_adipose_eos.rds",
  "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/integration/hu_adipose_eos.rds",
  "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/hu_eos/hu_adipose_eos.rds"
)

hu <- NULL
for (path in hu_candidates) {
  cat("  Trying:", path, "\n")
  if (file.exists(path)) {
    cat("  Found! Loading ...\n")
    hu <- tryCatch(
      readRDS(path),
      error = function(e) { cat("  ERROR loading:", e$message, "\n"); NULL }
    )
    if (!is.null(hu)) break
  }
}

if (is.null(hu)) {
  cat("ERROR: Hu reference data not found at any candidate path.\n")
  cat("  Searched:\n")
  for (p in hu_candidates) cat("    ", p, "\n")
  quit(status = 1)
}
cat("  Hu cells:", ncol(hu), "\n")

# --- Add dataset labels -------------------------------------------------------
cat("Adding dataset labels ...\n")
eos$dataset <- "Qiu_2025"
hu$dataset  <- "Hu_Reference"

# Preserve original cluster identities
eos$original_cluster <- as.character(Idents(eos))
if ("seurat_clusters" %in% colnames(hu@meta.data)) {
  hu$original_cluster <- as.character(hu$seurat_clusters)
} else {
  hu$original_cluster <- as.character(Idents(hu))
}

# --- Merge objects ------------------------------------------------------------
cat("Merging datasets ...\n")
combined <- merge(eos, y = hu, add.cell.ids = c("Qiu", "Hu"), project = "Qiu_Hu_Integration")
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
combined <- RunHarmony(combined, group.by.vars = "dataset", reduction = "pca",
                        dims.use = 1:30, verbose = TRUE)

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

# Figure 1: UMAP colored by dataset (split panel)
pdf(file.path(outdir, "hu_integration_umap_dataset.pdf"), width = 14, height = 6)
tryCatch({
  p1 <- DimPlot(combined, reduction = "umap", group.by = "dataset",
                 cols = c("Qiu_2025" = "#E63946", "Hu_Reference" = "#457B9D"),
                 pt.size = 0.3) +
    ggtitle("Combined") +
    theme_umap()

  p2 <- DimPlot(combined, reduction = "umap", group.by = "dataset",
                 split.by = "dataset",
                 cols = c("Qiu_2025" = "#E63946", "Hu_Reference" = "#457B9D"),
                 pt.size = 0.3) +
    ggtitle("Split by dataset") +
    theme_umap()

  print(p1 + p2)
}, error = function(e) cat("  Warning: dataset UMAP failed:", e$message, "\n"))
dev.off()

# Figure 2: UMAP colored by integrated cluster
pdf(file.path(outdir, "hu_integration_umap_cluster.pdf"), width = 8, height = 6)
tryCatch({
  p <- DimPlot(combined, reduction = "umap", label = TRUE, label.size = 4,
               pt.size = 0.3) +
    ggtitle("Integrated clusters (Qiu + Hu)") +
    theme_umap()
  print(p)
}, error = function(e) cat("  Warning: cluster UMAP failed:", e$message, "\n"))
dev.off()

# Figure 3: UMAP colored by original cluster identity
pdf(file.path(outdir, "hu_integration_umap_original_clusters.pdf"), width = 14, height = 6)
tryCatch({
  p1 <- DimPlot(subset(combined, dataset == "Qiu_2025"),
                 reduction = "umap", group.by = "original_cluster",
                 cols = eos_colors,
                 label = TRUE, pt.size = 0.3) +
    ggtitle("Qiu_2025 original clusters") +
    theme_umap()

  p2 <- DimPlot(subset(combined, dataset == "Hu_Reference"),
                 reduction = "umap", group.by = "original_cluster",
                 label = TRUE, pt.size = 0.3) +
    ggtitle("Hu_Reference original clusters") +
    theme_umap()

  print(p1 + p2)
}, error = function(e) cat("  Warning: original cluster UMAP failed:", e$message, "\n"))
dev.off()

# Figure 4: Confusion matrix - our clusters vs Hu labels
pdf(file.path(outdir, "hu_integration_confusion_matrix.pdf"), width = 10, height = 8)
tryCatch({
  # Build confusion matrix: integrated cluster vs original cluster, faceted by dataset
  conf_df <- combined@meta.data %>%
    select(dataset, integrated_cluster = seurat_clusters, original_cluster) %>%
    group_by(dataset, integrated_cluster, original_cluster) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(dataset, integrated_cluster) %>%
    mutate(frac = n / sum(n))

  p <- ggplot(conf_df, aes(x = integrated_cluster, y = original_cluster, fill = frac)) +
    geom_tile() +
    geom_text(aes(label = n), size = 2.5) +
    scale_fill_gradient(low = "white", high = "#2D6A4F", name = "Fraction") +
    facet_wrap(~dataset, scales = "free_y") +
    labs(x = "Integrated cluster", y = "Original cluster",
         title = "Cluster correspondence: Qiu vs Hu") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
  print(p)
}, error = function(e) cat("  Warning: confusion matrix failed:", e$message, "\n"))
dev.off()

# Figure 5: Composition comparison (proportion of each dataset per integrated cluster)
pdf(file.path(outdir, "hu_integration_composition.pdf"), width = 10, height = 6)
tryCatch({
  comp_df <- combined@meta.data %>%
    group_by(seurat_clusters, dataset) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(frac = n / sum(n))

  p <- ggplot(comp_df, aes(x = seurat_clusters, y = frac, fill = dataset)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = c("Qiu_2025" = "#E63946", "Hu_Reference" = "#457B9D")) +
    labs(x = "Integrated cluster", y = "Fraction",
         title = "Dataset composition per integrated cluster") +
    theme_pub()
  print(p)
}, error = function(e) cat("  Warning: composition plot failed:", e$message, "\n"))
dev.off()

# --- Save integrated object ---------------------------------------------------
cat("Saving integrated object ...\n")
saveRDS(combined, file.path(outdir, "qiu_hu_integrated.rds"))
cat("  Saved qiu_hu_integrated.rds\n")

cat("\n========== Hu integration complete ==========\n")
cat("Output directory:", outdir, "\n")
sessionInfo()
