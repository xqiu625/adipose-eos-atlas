#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

## All cell types UMAP + clustering
## SLURM: --mem=64G --cpus-per-task=4 --time=2:00:00 -p batch

# ── Libraries ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

# ── Shared theme ─────────────────────────────────────────────────────────────
SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

# ── Configuration ──────────────────────────────────────────────────────────────
BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load data ──────────────────────────────────────────────────────────────────
cat("========================================\n")
cat("All cell types UMAP + clustering\n")
cat("========================================\n\n")

input_file <- file.path(BASE, "seurat_qc_filtered.rds")
cat("Loading:", input_file, "\n")
sobj <- readRDS(input_file)
cat(sprintf("Loaded: %d genes x %d cells\n\n", nrow(sobj), ncol(sobj)))

# ── Standard Seurat workflow ───────────────────────────────────────────────────
cat("NormalizeData...\n")
sobj <- NormalizeData(sobj, verbose = FALSE)

cat("FindVariableFeatures (nfeatures = 2000)...\n")
sobj <- FindVariableFeatures(sobj, nfeatures = 2000, verbose = FALSE)

cat("ScaleData...\n")
sobj <- ScaleData(sobj, verbose = FALSE)

cat("RunPCA (npcs = 30)...\n")
sobj <- RunPCA(sobj, npcs = 30, verbose = FALSE)

cat("FindNeighbors (dims = 1:20, k.param = 30)...\n")
sobj <- FindNeighbors(sobj, dims = 1:20, k.param = 30, verbose = FALSE)

cat("RunUMAP (dims = 1:20, n.neighbors = 50, min.dist = 0.1)...\n")
sobj <- RunUMAP(sobj, dims = 1:20, n.neighbors = 50L, min.dist = 0.1,
                spread = 1, verbose = FALSE)

cat("FindClusters (res = 0.5)...\n")
sobj <- FindClusters(sobj, resolution = 0.5, verbose = FALSE)

cat("FindClusters (res = 0.8)...\n")
sobj <- FindClusters(sobj, resolution = 0.8, verbose = FALSE)

cat("FindClusters (res = 1.0)...\n")
sobj <- FindClusters(sobj, resolution = 1.0, verbose = FALSE)

# Set default identity to resolution 0.8
Idents(sobj) <- "RNA_snn_res.0.8"
sobj$seurat_clusters <- sobj$RNA_snn_res.0.8

cat(sprintf("\nClusters at res 0.5: %d\n", length(unique(sobj$RNA_snn_res.0.5))))
cat(sprintf("Clusters at res 0.8: %d\n", length(unique(sobj$RNA_snn_res.0.8))))
cat(sprintf("Clusters at res 1.0: %d\n", length(unique(sobj$RNA_snn_res.1))))
cat("\nCluster sizes (res 0.8):\n")
print(sort(table(sobj$seurat_clusters), decreasing = TRUE))

# ── Figures ────────────────────────────────────────────────────────────────────
cat("\n========================================\n")
cat("Generating figures\n")
cat("========================================\n\n")

# 1. UMAP by clusters (res 0.8)
cat("UMAP colored by seurat_clusters (res 0.8)...\n")
p1 <- DimPlot(sobj, group.by = "seurat_clusters", label = TRUE, repel = TRUE,
              pt.size = 0.3) +
  ggtitle("Clusters (res 0.8)") + theme_umap() + NoLegend()

save_pdf(p1, file.path(fig_dir, "umap_clusters_res0.8.pdf"),
         width = 8, height = 7)

# 2. UMAP by condition
cat("UMAP colored by condition...\n")
p2 <- DimPlot(sobj, group.by = "condition", pt.size = 0.3,
              cols = cond_colors) +
  ggtitle("Condition") + theme_umap()

save_pdf(p2, file.path(fig_dir, "umap_condition.pdf"),
         width = 10, height = 7)

# 3. UMAP by orig.ident (sample)
cat("UMAP colored by orig.ident (sample)...\n")
p3 <- DimPlot(sobj, group.by = "orig.ident", pt.size = 0.3) +
  ggtitle("Sample") + theme_umap()

save_pdf(p3, file.path(fig_dir, "umap_sample.pdf"),
         width = 12, height = 7)

# 4. UMAP by gem_reaction (batch check)
cat("UMAP colored by gem_reaction (batch check)...\n")
p4 <- DimPlot(sobj, group.by = "gem_reaction", pt.size = 0.3) +
  ggtitle("GEM Reaction (Batch)") + theme_umap()

save_pdf(p4, file.path(fig_dir, "umap_gem_reaction.pdf"),
         width = 9, height = 7)

# Combined panel
cat("Saving combined panel...\n")
p_combined <- (p1 | p2) / (p3 | p4) +
  plot_annotation(title = "Adipose Atlas - All Cells",
                  theme = theme(plot.title = element_text(size = 16, face = "bold",
                                                         family = "Arial")))

save_pdf(p_combined, file.path(fig_dir, "umap_combined_panel.pdf"),
         width = 18, height = 14)

# ── Save ───────────────────────────────────────────────────────────────────────
cat("\n========================================\n")
cat("Saving results\n")
cat("========================================\n\n")

out_file <- file.path(BASE, "all_cells_clustered.rds")
cat("Saving:", out_file, "\n")
saveRDS(sobj, out_file)

cat(sprintf("Figures saved to: %s\n", fig_dir))
cat("Done.\n")
