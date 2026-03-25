#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --job-name=eos_subcluster
#SBATCH --output=eos_subcluster_%j.out

# ==============================================================================
# 04_eos_subcluster.R
# Subset eosinophils, re-cluster, find markers, and generate composition plots
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
})

SCRIPTS <- dirname(sys.frame(1)$ofile)
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"

# --- Load annotated object and subset eosinophils ----------------------------
cat("Loading annotated object...\n")
obj <- readRDS(file.path(BASE, "all_cells_annotated.rds"))

cat("Subsetting eosinophils...\n")
eos <- subset(obj, broad_celltype == "Eosinophil")
cat(sprintf("  %d eosinophils extracted\n", ncol(eos)))
rm(obj); gc()

# --- Re-process eosinophils ---------------------------------------------------
cat("Re-normalizing and finding variable features...\n")
eos <- NormalizeData(eos, verbose = FALSE)
eos <- FindVariableFeatures(eos, nfeatures = 2000, verbose = FALSE)
eos <- ScaleData(eos, verbose = FALSE)

cat("Running PCA...\n")
eos <- RunPCA(eos, npcs = 30, verbose = FALSE)

cat("Running UMAP...\n")
eos <- RunUMAP(eos, dims = 1:20, verbose = FALSE)

# --- Clustering ---------------------------------------------------------------
cat("Finding neighbors and clustering (resolution 0.8)...\n")
eos <- FindNeighbors(eos, dims = 1:20, verbose = FALSE)
eos <- FindClusters(eos, resolution = 0.8, verbose = FALSE)

# Label clusters as C0, C1, ...
eos$eos_cluster <- paste0("C", as.integer(as.character(Idents(eos))))
Idents(eos) <- "eos_cluster"
cat(sprintf("  Found %d clusters: %s\n", length(unique(eos$eos_cluster)),
            paste(sort(unique(eos$eos_cluster)), collapse = ", ")))

# --- Find cluster markers -----------------------------------------------------
cat("Finding cluster markers (FindAllMarkers)...\n")
markers <- FindAllMarkers(eos, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)
cat(sprintf("  %d marker genes identified\n", nrow(markers)))

# --- Figures ------------------------------------------------------------------
fig_dir <- file.path(BASE, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Determine which colors to use based on clusters found
clusters_present <- sort(unique(eos$eos_cluster))
eos_colors_use <- eos_colors[clusters_present]

# UMAP by cluster
cat("Plotting UMAP by cluster...\n")
pdf(file.path(fig_dir, "04_eos_umap_cluster.pdf"), width = 7, height = 6)
print(DimPlot(eos, group.by = "eos_cluster", cols = eos_colors_use, pt.size = 0.5) +
        ggtitle("Eosinophil Subclusters") + theme_umap())
dev.off()

# UMAP by condition
cat("Plotting UMAP by condition...\n")
pdf(file.path(fig_dir, "04_eos_umap_condition.pdf"), width = 7, height = 6)
print(DimPlot(eos, group.by = "condition", pt.size = 0.5) +
        ggtitle("Eosinophils by Condition") + theme_umap())
dev.off()

# UMAP by sample
cat("Plotting UMAP by sample...\n")
pdf(file.path(fig_dir, "04_eos_umap_sample.pdf"), width = 8, height = 6)
print(DimPlot(eos, group.by = "orig.ident", pt.size = 0.5) +
        ggtitle("Eosinophils by Sample") + theme_umap())
dev.off()

# Stacked bar: composition by condition
cat("Plotting composition by condition...\n")
comp_df <- eos@meta.data %>%
  group_by(condition, eos_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(frac = n / sum(n))

pdf(file.path(fig_dir, "04_eos_composition_condition.pdf"), width = 8, height = 6)
print(ggplot(comp_df, aes(x = condition, y = frac, fill = eos_cluster)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = eos_colors_use) +
        labs(x = "Condition", y = "Fraction", fill = "Cluster",
             title = "Eosinophil Cluster Composition by Condition") +
        theme_pub() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

# Top 5 marker heatmap
cat("Plotting top 5 marker heatmap...\n")
top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

pdf(file.path(fig_dir, "04_eos_marker_heatmap.pdf"), width = 12, height = 10)
print(DoHeatmap(eos, features = top5$gene, group.colors = eos_colors_use) +
        ggtitle("Top 5 Markers per Eosinophil Cluster") +
        theme(axis.text.y = element_text(size = 8)))
dev.off()

# --- Save outputs -------------------------------------------------------------
cat("Saving eos_subclustered.rds...\n")
saveRDS(eos, file.path(BASE, "eos_subclustered.rds"))

cat("Saving eos_cluster_markers.csv...\n")
write.csv(markers, file.path(BASE, "eos_cluster_markers.csv"), row.names = FALSE)

cat("Done! 04_eos_subcluster.R completed successfully.\n")
