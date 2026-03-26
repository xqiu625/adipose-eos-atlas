#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

## Eosinophil subclustering with resolution justification (clustree + silhouette)
## Addresses Reviewer 1 point 3d
## SLURM: --mem=128G --cpus-per-task=4 --time=2:00:00 -p batch

# ── Libraries ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(RColorBrewer)
  library(clustree)
  library(cluster)
  library(patchwork)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures")
supp_dir <- file.path(fig_dir, "supplement")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load annotated object and subset eosinophils ─────────────────────────────
cat("Loading annotated object...\n")
obj <- readRDS(file.path(BASE, "all_cells_annotated.rds"))

cat("Subsetting eosinophils...\n")
eos <- subset(obj, celltype == "Eosinophils")
cat(sprintf("  %d eosinophils extracted\n", ncol(eos)))
rm(obj); gc()

# ── Re-process eosinophils ───────────────────────────────────────────────────
cat("Re-normalizing and finding variable features...\n")
eos <- NormalizeData(eos, verbose = FALSE)
eos <- FindVariableFeatures(eos, nfeatures = 2000, verbose = FALSE)
eos <- ScaleData(eos, verbose = FALSE)

cat("Running PCA...\n")
eos <- RunPCA(eos, npcs = 30, verbose = FALSE)

cat("Running UMAP...\n")
eos <- RunUMAP(eos, dims = 1:20, verbose = FALSE)

cat("Finding neighbors...\n")
eos <- FindNeighbors(eos, dims = 1:20, verbose = FALSE)

# ── Step 1: Multi-resolution clustering + clustree ───────────────────────────
cat("\n========================================\n")
cat("Step 1: Multi-resolution clustering\n")
cat("========================================\n\n")

resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5)
cat(sprintf("Clustering at %d resolutions: %s\n",
            length(resolutions), paste(resolutions, collapse = ", ")))

for (res in resolutions) {
  eos <- FindClusters(eos, resolution = res, verbose = FALSE)
  n_clust <- length(unique(eos@meta.data[[paste0("RNA_snn_res.", res)]]))
  cat(sprintf("  res %.1f -> %d clusters\n", res, n_clust))
}

# Clustree diagram
cat("Generating clustree diagram...\n")
pdf(file.path(supp_dir, "fig_clustree.pdf"), width = 12, height = 10)
print(
  clustree(eos, prefix = "RNA_snn_res.") +
    ggtitle("Clustree: Cluster Stability Across Resolutions") +
    theme(plot.title = element_text(size = 14, face = "bold"))
)
dev.off()

# ── Step 2: Silhouette analysis ──────────────────────────────────────────────
cat("\n========================================\n")
cat("Step 2: Silhouette analysis\n")
cat("========================================\n\n")

pca_emb <- Embeddings(eos, "pca")[, 1:20]

# Subsample for silhouette if needed
n_cells <- nrow(pca_emb)
if (n_cells > 5000) {
  cat(sprintf("Subsampling %d -> 5000 cells for silhouette...\n", n_cells))
  set.seed(42)
  sil_idx <- sample(n_cells, 5000)
} else {
  sil_idx <- seq_len(n_cells)
}
dist_mat <- dist(pca_emb[sil_idx, ])

sil_results <- data.frame(
  resolution = numeric(), n_clusters = integer(),
  mean_sil = numeric(), median_sil = numeric()
)

for (res in resolutions) {
  col_name <- paste0("RNA_snn_res.", res)
  clusters <- as.integer(as.character(eos@meta.data[sil_idx, col_name]))
  n_clust  <- length(unique(clusters))

  if (n_clust < 2) {
    sil_results <- rbind(sil_results, data.frame(
      resolution = res, n_clusters = n_clust, mean_sil = NA, median_sil = NA))
    next
  }

  sil <- silhouette(clusters, dist_mat)
  mean_s   <- mean(sil[, 3])
  median_s <- median(sil[, 3])
  cat(sprintf("  res %.1f: %d clusters, mean sil = %.4f\n", res, n_clust, mean_s))

  sil_results <- rbind(sil_results, data.frame(
    resolution = res, n_clusters = n_clust,
    mean_sil = mean_s, median_sil = median_s))
}

# Silhouette + cluster count plot
p_sil <- ggplot(sil_results, aes(x = resolution, y = mean_sil)) +
  geom_line(color = "#2166AC", linewidth = 1) +
  geom_point(color = "#2166AC", size = 3) +
  geom_point(data = sil_results %>% filter(resolution == 0.8),
             color = "red", size = 5, shape = 18) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "red", alpha = 0.5) +
  labs(x = "Resolution", y = "Mean Silhouette Score",
       title = "Silhouette Score vs Resolution") +
  theme_pub()

p_nclust <- ggplot(sil_results, aes(x = resolution, y = n_clusters)) +
  geom_line(color = "#D6604D", linewidth = 1) +
  geom_point(color = "#D6604D", size = 3) +
  geom_point(data = sil_results %>% filter(resolution == 0.8),
             color = "red", size = 5, shape = 18) +
  geom_vline(xintercept = 0.8, linetype = "dashed", color = "red", alpha = 0.5) +
  scale_y_continuous(breaks = seq(0, max(sil_results$n_clusters) + 2, by = 1)) +
  labs(x = "Resolution", y = "Number of Clusters",
       title = "Number of Clusters vs Resolution") +
  theme_pub()

pdf(file.path(supp_dir, "fig_silhouette_resolution.pdf"), width = 12, height = 5)
print(p_sil + p_nclust + plot_layout(ncol = 2))
dev.off()

# UMAP panels at selected resolutions
cat("Plotting UMAP at selected resolutions...\n")
umap_resolutions <- c(0.5, 0.8, 1.0, 1.5)
umap_plots <- list()
for (res in umap_resolutions) {
  col_name <- paste0("RNA_snn_res.", res)
  eos@meta.data[[paste0("label_", res)]] <- paste0("C",
    as.integer(as.character(eos@meta.data[[col_name]])))
  n_clust <- length(unique(eos@meta.data[[col_name]]))
  p <- DimPlot(eos, group.by = paste0("label_", res), pt.size = 0.3,
               label = TRUE, label.size = 3) +
    ggtitle(sprintf("Resolution %.1f (%d clusters)", res, n_clust)) +
    theme_umap() + theme(legend.position = "none")
  if (res == 0.8) p <- p + theme(plot.background = element_rect(
    color = "red", linewidth = 2, fill = NA))
  umap_plots[[as.character(res)]] <- p
}

pdf(file.path(supp_dir, "fig_umap_resolutions.pdf"), width = 12, height = 10)
print(wrap_plots(umap_plots, ncol = 2) +
  plot_annotation(title = "UMAP at Multiple Resolutions (red = chosen)", tag_levels = "A"))
dev.off()

# Save metrics
write.csv(sil_results, file.path(supp_dir, "resolution_metrics.csv"), row.names = FALSE)

# Print summary
cat("\nResolution Selection Summary:\n")
print(sil_results)
res08 <- sil_results %>% filter(resolution == 0.8)
cat(sprintf("\nChosen resolution 0.8: %d clusters, mean sil = %.4f\n",
            res08$n_clusters, res08$mean_sil))

# ── Step 3: Set final clusters at res 0.8 ────────────────────────────────────
cat("\n========================================\n")
cat("Step 3: Final clustering (res 0.8)\n")
cat("========================================\n\n")

eos$eos_cluster <- paste0("C", as.integer(as.character(
  eos@meta.data[["RNA_snn_res.0.8"]])))
Idents(eos) <- "eos_cluster"
cat(sprintf("  Found %d clusters: %s\n", length(unique(eos$eos_cluster)),
            paste(sort(unique(eos$eos_cluster)), collapse = ", ")))

# ── Step 4: Find cluster markers ─────────────────────────────────────────────
cat("Finding cluster markers (FindAllMarkers)...\n")
markers <- FindAllMarkers(eos, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)
cat(sprintf("  %d marker genes identified\n", nrow(markers)))

# ── Step 5: Figures ──────────────────────────────────────────────────────────
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

# Composition by condition
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

# ── Save outputs ─────────────────────────────────────────────────────────────
cat("\nSaving eos_subclustered.rds...\n")
saveRDS(eos, file.path(BASE, "eos_subclustered.rds"))

cat("Saving eos_cluster_markers.csv...\n")
write.csv(markers, file.path(BASE, "eos_cluster_markers.csv"), row.names = FALSE)

cat("Done! 04_eos_subcluster.R completed successfully.\n")
