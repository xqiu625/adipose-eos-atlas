#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# =============================================================================
# 24_integrate_qiu_hu_gurtner.R
# Three-way Harmony integration: Qiu (HFD adipose) + Hu (healthy adipose) +
# Gurtner (blood/bone marrow) eosinophils.
# Exports integrated object + CellRank-ready files.
#
# Input:  eos_subclustered.rds, hu_adipose_eos.rds, Gurtner eos .rds
# Output: three_way_integrated.rds, for_cellrank/ exports, PDFs
# SLURM:  epyc, 64G, 8 CPUs, 4h
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
outdir <- file.path(BASE, "integration", "three_way")
fig_dir <- file.path(outdir, "figures")
export_dir <- file.path(outdir, "for_cellrank")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(export_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load all three datasets
# =============================================================================

# --- Our eosinophils (HFD adipose, C0-C6) ---
cat("Loading Qiu adipose eosinophils ...\n")
qiu <- readRDS(file.path(BASE, "eos_subclustered.rds"))
cat(sprintf("  Qiu: %d cells\n", ncol(qiu)))

# --- Hu healthy adipose eosinophils ---
cat("Loading Hu healthy adipose eosinophils ...\n")
hu_candidates <- c(
  "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/hu_2026/hu_adipose_eos.rds",
  "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/integration/hu_adipose_eos.rds",
  "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/hu_eos/hu_adipose_eos.rds"
)
hu <- NULL
for (path in hu_candidates) {
  if (file.exists(path)) {
    cat(sprintf("  Found: %s\n", path))
    hu <- readRDS(path)
    break
  }
}
if (is.null(hu)) { cat("ERROR: Hu data not found\n"); quit(status = 1) }
cat(sprintf("  Hu: %d cells\n", ncol(hu)))

# --- Gurtner blood/BM eosinophils ---
cat("Loading Gurtner blood/BM eosinophils ...\n")
gurtner_candidates <- c(
  "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/Gurtner_eosinophil_subset.rds",
  "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/GSE182001_eosinophil_subset.rds",
  "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature/Gurtner_merged_tissues_seurat_object.rds"
)
gurtner <- NULL
for (path in gurtner_candidates) {
  if (file.exists(path)) {
    cat(sprintf("  Found: %s\n", path))
    gurtner <- readRDS(path)
    break
  }
}
if (is.null(gurtner)) { cat("ERROR: Gurtner data not found\n"); quit(status = 1) }
cat(sprintf("  Gurtner: %d cells\n", ncol(gurtner)))

# =============================================================================
# 2. Harmonize metadata across datasets
# =============================================================================
cat("\nPreparing metadata ...\n")

# Dataset labels
qiu$dataset <- "Qiu_2025_Adipose"
hu$dataset  <- "Hu_2026_Healthy"
gurtner$dataset <- "Gurtner_2022_Blood_BM"

# Original cluster identities
qiu$original_cluster <- as.character(Idents(qiu))
if ("seurat_clusters" %in% colnames(hu@meta.data)) {
  hu$original_cluster <- as.character(hu$seurat_clusters)
} else {
  hu$original_cluster <- as.character(Idents(hu))
}
if ("seurat_clusters" %in% colnames(gurtner@meta.data)) {
  gurtner$original_cluster <- as.character(gurtner$seurat_clusters)
} else {
  gurtner$original_cluster <- as.character(Idents(gurtner))
}

# Organ / tissue labels
qiu$organ <- "Adipose_HFD"
hu$organ  <- "Adipose_Healthy"
if (!"organ" %in% colnames(gurtner@meta.data)) {
  gurtner$organ <- "Unknown"
}

# Condition
if (!"condition" %in% colnames(hu@meta.data)) {
  hu$condition <- "Healthy"
}
if (!"condition" %in% colnames(gurtner@meta.data)) {
  gurtner$condition <- gurtner$organ
}

# Study label (for display)
qiu$study <- "Qiu (HFD adipose)"
hu$study  <- "Hu (healthy adipose)"
gurtner$study <- "Gurtner (blood/BM)"

cat("  Qiu organs:    ", paste(unique(qiu$organ), collapse = ", "), "\n")
cat("  Hu organs:     ", paste(unique(hu$organ), collapse = ", "), "\n")
cat("  Gurtner organs:", paste(unique(gurtner$organ), collapse = ", "), "\n")
cat("  Gurtner organ distribution:\n")
print(table(gurtner$organ))

# =============================================================================
# 3. Merge and integrate
# =============================================================================
cat("\nMerging all three datasets ...\n")
combined <- merge(qiu, y = list(hu, gurtner),
                  add.cell.ids = c("Qiu", "Hu", "Gurtner"),
                  project = "Three_Way_Integration")
cat(sprintf("  Combined: %d cells\n", ncol(combined)))

cat("Normalizing ...\n")
combined <- NormalizeData(combined, verbose = FALSE)

cat("Finding variable features ...\n")
combined <- FindVariableFeatures(combined, selection.method = "vst",
                                  nfeatures = 2000, verbose = FALSE)

cat("Scaling ...\n")
combined <- ScaleData(combined, verbose = FALSE)

cat("Running PCA ...\n")
combined <- RunPCA(combined, npcs = 30, verbose = FALSE)

cat("Running Harmony (group.by = dataset) ...\n")
combined <- RunHarmony(combined, group.by.vars = "dataset",
                        reduction.use = "pca", verbose = TRUE)

cat("Running UMAP on Harmony ...\n")
combined <- RunUMAP(combined, reduction = "harmony", dims = 1:30, verbose = FALSE)

cat("Finding neighbors and clusters ...\n")
combined <- FindNeighbors(combined, reduction = "harmony", dims = 1:30, verbose = FALSE)
combined <- FindClusters(combined, resolution = 0.5, verbose = FALSE)
cat(sprintf("  Integrated clusters: %d\n", length(unique(Idents(combined)))))

# =============================================================================
# 4. Figures
# =============================================================================
cat("\nGenerating figures ...\n")

dataset_cols <- c("Qiu_2025_Adipose" = "#E63946",
                   "Hu_2026_Healthy" = "#4DAF4A",
                   "Gurtner_2022_Blood_BM" = "#1D3557")

# Fig 1: UMAP by dataset (combined + split)
pdf(file.path(fig_dir, "three_way_umap_dataset.pdf"), width = 18, height = 6)
tryCatch({
  p1 <- DimPlot(combined, reduction = "umap", group.by = "dataset",
                 cols = dataset_cols, pt.size = 0.2) +
    ggtitle("Three-way Integration") + theme_umap()
  p2 <- DimPlot(combined, reduction = "umap", group.by = "dataset",
                 split.by = "dataset", cols = dataset_cols, pt.size = 0.2) +
    ggtitle("Split by Dataset") + theme_umap()
  print(p1 + p2 + plot_layout(widths = c(1, 2)))
}, error = function(e) cat("  Fig 1 failed:", e$message, "\n"))
dev.off()

# Fig 2: UMAP by integrated cluster
pdf(file.path(fig_dir, "three_way_umap_cluster.pdf"), width = 8, height = 6)
tryCatch({
  print(DimPlot(combined, reduction = "umap", label = TRUE, label.size = 4,
                 pt.size = 0.2) +
    ggtitle("Integrated Clusters") + theme_umap())
}, error = function(e) cat("  Fig 2 failed:", e$message, "\n"))
dev.off()

# Fig 3: UMAP by original cluster (split by dataset)
pdf(file.path(fig_dir, "three_way_umap_original_clusters.pdf"), width = 20, height = 6)
tryCatch({
  p1 <- DimPlot(subset(combined, dataset == "Qiu_2025_Adipose"),
                 reduction = "umap", group.by = "original_cluster",
                 cols = eos_colors, label = TRUE, pt.size = 0.3) +
    ggtitle("Qiu: Eos Clusters (C0-C6)") + theme_umap()
  p2 <- DimPlot(subset(combined, dataset == "Hu_2026_Healthy"),
                 reduction = "umap", group.by = "original_cluster",
                 label = TRUE, pt.size = 0.3) +
    ggtitle("Hu: Healthy Eos") + theme_umap()
  p3 <- DimPlot(subset(combined, dataset == "Gurtner_2022_Blood_BM"),
                 reduction = "umap", group.by = "original_cluster",
                 label = TRUE, pt.size = 0.3) +
    ggtitle("Gurtner: Blood/BM Eos") + theme_umap()
  print(p1 + p2 + p3)
}, error = function(e) cat("  Fig 3 failed:", e$message, "\n"))
dev.off()

# Fig 4: UMAP by organ
pdf(file.path(fig_dir, "three_way_umap_organ.pdf"), width = 10, height = 6)
tryCatch({
  print(DimPlot(combined, reduction = "umap", group.by = "organ",
                 pt.size = 0.2) +
    ggtitle("Tissue / Organ") + theme_umap())
}, error = function(e) cat("  Fig 4 failed:", e$message, "\n"))
dev.off()

# Fig 5: Dataset composition per integrated cluster
pdf(file.path(fig_dir, "three_way_composition.pdf"), width = 10, height = 6)
tryCatch({
  comp_df <- combined@meta.data %>%
    group_by(seurat_clusters, dataset) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(seurat_clusters) %>%
    mutate(frac = n / sum(n))
  print(ggplot(comp_df, aes(x = seurat_clusters, y = frac, fill = dataset)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = dataset_cols) +
    labs(x = "Integrated Cluster", y = "Fraction",
         title = "Dataset Composition per Cluster") +
    theme_pub())
}, error = function(e) cat("  Fig 5 failed:", e$message, "\n"))
dev.off()

# Fig 6: Confusion matrix
pdf(file.path(fig_dir, "three_way_confusion_matrix.pdf"), width = 14, height = 8)
tryCatch({
  conf_df <- combined@meta.data %>%
    select(dataset, integrated_cluster = seurat_clusters, original_cluster) %>%
    group_by(dataset, integrated_cluster, original_cluster) %>%
    summarise(n = n(), .groups = "drop") %>%
    group_by(dataset, integrated_cluster) %>%
    mutate(frac = n / sum(n))
  print(ggplot(conf_df, aes(x = integrated_cluster, y = original_cluster, fill = frac)) +
    geom_tile() +
    geom_text(aes(label = n), size = 2) +
    scale_fill_gradient(low = "white", high = "#264653", name = "Fraction") +
    facet_wrap(~dataset, scales = "free_y") +
    labs(x = "Integrated Cluster", y = "Original Cluster",
         title = "Cluster Correspondence Across Datasets") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)))
}, error = function(e) cat("  Fig 6 failed:", e$message, "\n"))
dev.off()

# =============================================================================
# 5. Save integrated object
# =============================================================================
cat("\nSaving integrated object ...\n")
saveRDS(combined, file.path(outdir, "three_way_integrated.rds"))

# =============================================================================
# 6. Export for CellRank (Python)
# =============================================================================
cat("Exporting for CellRank ...\n")
combined <- JoinLayers(combined)

# Counts
counts_mat <- GetAssayData(combined, layer = "counts")
Matrix::writeMM(counts_mat, file.path(export_dir, "counts.mtx"))

# Gene names and barcodes
writeLines(rownames(combined), file.path(export_dir, "genes.txt"))
writeLines(colnames(combined), file.path(export_dir, "barcodes.txt"))

# Metadata
write.csv(combined@meta.data, file.path(export_dir, "metadata.csv"))

# Harmony embeddings
write.csv(Embeddings(combined, "harmony"), file.path(export_dir, "harmony_embeddings.csv"))

# UMAP embeddings
write.csv(Embeddings(combined, "umap"), file.path(export_dir, "umap_embeddings.csv"))

cat(sprintf("  Exported to: %s\n", export_dir))

cat("\n========== Three-way integration complete ==========\n")
cat("Output directory:", outdir, "\n")
cat("Cells per dataset:\n")
print(table(combined$dataset))
cat("\nCells per organ:\n")
print(table(combined$organ))
sessionInfo()
