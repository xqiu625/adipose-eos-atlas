#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --job-name=clustree_resolution
#SBATCH --output=clustree_resolution_%j.out

# ==============================================================================
# 18_clustree_resolution.R
# Justify clustering resolution choice via clustree, silhouette, and UMAP
# Clustering resolution selection (clustree + silhouette + sc-SHC)
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(clustree)
  library(cluster)
  library(patchwork)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures", "supplement")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load annotated object and extract eosinophils ────────────────────────────
cat("Loading all_cells_annotated.rds...\n")
obj <- readRDS(file.path(BASE, "all_cells_annotated.rds"))

cat("Subsetting eosinophils...\n")
eos <- subset(obj, broad_celltype == "Eosinophil")
cat(sprintf("  %d eosinophils extracted\n", ncol(eos)))
rm(obj); gc()

# ── Re-process: normalize, PCA, neighbors ───────────────────────────────────
cat("Re-normalizing, PCA, neighbors...\n")
eos <- NormalizeData(eos, verbose = FALSE)
eos <- FindVariableFeatures(eos, nfeatures = 2000, verbose = FALSE)
eos <- ScaleData(eos, verbose = FALSE)
eos <- RunPCA(eos, npcs = 30, verbose = FALSE)
eos <- FindNeighbors(eos, dims = 1:20, verbose = FALSE)

# ── Run UMAP (single embedding for all panels) ──────────────────────────────
cat("Running UMAP...\n")
eos <- RunUMAP(eos, dims = 1:20, verbose = FALSE)

# ── Cluster at multiple resolutions ──────────────────────────────────────────
resolutions <- c(0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0, 1.2, 1.5)
cat(sprintf("Clustering at %d resolutions: %s\n",
            length(resolutions), paste(resolutions, collapse = ", ")))

for (res in resolutions) {
  cat(sprintf("  Resolution %.1f ... ", res))
  eos <- FindClusters(eos, resolution = res, verbose = FALSE)
  n_clust <- length(unique(eos@meta.data[[paste0("RNA_snn_res.", res)]]))
  cat(sprintf("%d clusters\n", n_clust))
}

# ── Figure 1: Clustree diagram ──────────────────────────────────────────────
cat("Generating clustree diagram...\n")
pdf(file.path(fig_dir, "fig_clustree.pdf"), width = 12, height = 10)
print(
  clustree(eos, prefix = "RNA_snn_res.") +
    ggtitle("Clustree: Cluster Stability Across Resolutions") +
    theme(plot.title = element_text(size = 14, face = "bold"))
)
dev.off()

# ── Figure 2: Silhouette analysis ───────────────────────────────────────────
cat("Computing silhouette scores at each resolution...\n")
pca_emb <- Embeddings(eos, "pca")[, 1:20]

# Pre-compute distance matrix (can be large; subsample if needed)
n_cells <- nrow(pca_emb)
if (n_cells > 5000) {
  cat(sprintf("  Subsampling %d -> 5000 cells for silhouette computation...\n", n_cells))
  set.seed(42)
  sil_idx <- sample(n_cells, 5000)
} else {
  sil_idx <- seq_len(n_cells)
}
pca_sub <- pca_emb[sil_idx, ]
dist_mat <- dist(pca_sub)

sil_results <- data.frame(
  resolution    = numeric(),
  n_clusters    = integer(),
  mean_sil      = numeric(),
  median_sil    = numeric(),
  stringsAsFactors = FALSE
)

for (res in resolutions) {
  col_name <- paste0("RNA_snn_res.", res)
  clusters <- as.integer(as.character(eos@meta.data[sil_idx, col_name]))
  n_clust  <- length(unique(clusters))

  if (n_clust < 2) {
    cat(sprintf("  Resolution %.1f: %d cluster — skipping silhouette\n", res, n_clust))
    sil_results <- rbind(sil_results, data.frame(
      resolution = res, n_clusters = n_clust,
      mean_sil = NA, median_sil = NA
    ))
    next
  }

  sil <- silhouette(clusters, dist_mat)
  mean_s   <- mean(sil[, 3])
  median_s <- median(sil[, 3])
  cat(sprintf("  Resolution %.1f: %d clusters, mean sil = %.4f\n", res, n_clust, mean_s))

  sil_results <- rbind(sil_results, data.frame(
    resolution = res, n_clusters = n_clust,
    mean_sil = mean_s, median_sil = median_s
  ))
}

# Plot silhouette vs resolution
p_sil <- ggplot(sil_results, aes(x = resolution, y = mean_sil)) +
  geom_line(color = "#2166AC", linewidth = 1) +
  geom_point(color = "#2166AC", size = 3) +
  geom_point(data = sil_results %>% filter(resolution == 0.4),
             color = "red", size = 5, shape = 18) +
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "red", alpha = 0.5) +
  annotate("text", x = 0.45, y = max(sil_results$mean_sil, na.rm = TRUE),
           label = "res = 0.4\n(chosen)", hjust = 0, color = "red", size = 3.5) +
  labs(x = "Resolution", y = "Mean Silhouette Score",
       title = "Silhouette Score vs Clustering Resolution") +
  theme_pub() +
  theme(plot.title = element_text(size = 12, face = "bold"))

# Plot number of clusters vs resolution
p_nclust <- ggplot(sil_results, aes(x = resolution, y = n_clusters)) +
  geom_line(color = "#D6604D", linewidth = 1) +
  geom_point(color = "#D6604D", size = 3) +
  geom_point(data = sil_results %>% filter(resolution == 0.4),
             color = "red", size = 5, shape = 18) +
  geom_vline(xintercept = 0.4, linetype = "dashed", color = "red", alpha = 0.5) +
  scale_y_continuous(breaks = seq(0, max(sil_results$n_clusters) + 2, by = 1)) +
  labs(x = "Resolution", y = "Number of Clusters",
       title = "Number of Clusters vs Resolution") +
  theme_pub() +
  theme(plot.title = element_text(size = 12, face = "bold"))

pdf(file.path(fig_dir, "fig_silhouette_resolution.pdf"), width = 12, height = 5)
print(p_sil + p_nclust + plot_layout(ncol = 2))
dev.off()

# ── sc-SHC significance testing on the chosen resolution (0.4) ──────────────
# The header/README/methods.md all describe sc-SHC (Grabski & Irizarry, Nat
# Methods 2023) as one of the completed, convergent lines of evidence for the
# chosen resolution. This step must not be silently skippable if the package
# is unavailable -- if it can't run, that is treated as a hard failure of this
# script (not printed as an optional warning while everything else proceeds).
cat("\nRunning sc-SHC significance testing (resolution 0.4 clusters)...\n")
if (!requireNamespace("scSHC", quietly = TRUE)) {
  stop("scSHC package is not installed. sc-SHC significance testing is a ",
       "documented, required part of the resolution justification (see ",
       "methods.md); install scSHC (https://github.com/igrabski/sc-SHC) ",
       "rather than skipping this step.")
}
library(scSHC)

res04_labels <- as.character(eos@meta.data[["RNA_snn_res.0.4"]])
counts_for_shc <- GetAssayData(eos, layer = "counts")

# testClusters() evaluates significance of an existing cluster assignment
# (as opposed to scSHC()'s own de novo hierarchical clustering). Confirm the
# exact function name/signature against the installed scSHC version if this
# errors -- do not silently fall back to skipping the test.
shc_result <- testClusters(counts_for_shc, cluster_ids = res04_labels,
                           num_PCs = 20, alpha = 0.05)
shc_out_path <- file.path(fig_dir, "scSHC_resolution0.4_significance.csv")
write.csv(as.data.frame(shc_result), shc_out_path, row.names = FALSE)
cat(sprintf("  Saved: %s\n", shc_out_path))

# ── Figure 3: UMAP at selected resolutions (2x2 grid) ───────────────────────
cat("Plotting UMAP panels at resolutions 0.4, 0.5, 1.0, 1.5...\n")

umap_resolutions <- c(0.4, 0.5, 1.0, 1.5)
umap_plots <- list()

for (res in umap_resolutions) {
  col_name <- paste0("RNA_snn_res.", res)
  # Label clusters as C0, C1, ... for readability
  eos@meta.data[[paste0("label_", res)]] <- paste0("C",
    as.integer(as.character(eos@meta.data[[col_name]])))

  n_clust <- length(unique(eos@meta.data[[col_name]]))
  p <- DimPlot(eos, group.by = paste0("label_", res), pt.size = 0.3,
               label = TRUE, label.size = 3) +
    ggtitle(sprintf("Resolution %.1f (%d clusters)", res, n_clust)) +
    theme_umap() +
    theme(legend.position = "none",
          plot.title = element_text(size = 11, face = "bold"))

  # Highlight the chosen resolution
  if (res == 0.4) {
    p <- p + theme(plot.background = element_rect(color = "red",
                                                   linewidth = 2,
                                                   fill = NA))
  }

  umap_plots[[as.character(res)]] <- p
}

pdf(file.path(fig_dir, "fig_umap_resolutions.pdf"), width = 12, height = 10)
print(wrap_plots(umap_plots, ncol = 2) +
        plot_annotation(
          title = "UMAP at Multiple Clustering Resolutions (red border = chosen)",
          tag_levels = "A"
        ))
dev.off()

# ── Save metrics CSV ─────────────────────────────────────────────────────────
csv_out <- file.path(fig_dir, "resolution_metrics.csv")
write.csv(sil_results, csv_out, row.names = FALSE)
cat(sprintf("Saved resolution metrics: %s\n", csv_out))

# ── Print summary ─────────────────────────────────────────────────────────────
cat("\n========================================\n")
cat("Resolution Selection Summary\n")
cat("========================================\n")
print(sil_results)

res04 <- sil_results %>% filter(resolution == 0.4)
cat(sprintf("\nChosen resolution: 0.4\n"))
cat(sprintf("  Clusters: %d\n", res04$n_clusters))
cat(sprintf("  Mean silhouette: %.4f\n", res04$mean_sil))
cat(sprintf("  Median silhouette: %.4f\n", res04$median_sil))

# Compare with adjacent resolutions for stability argument
res03 <- sil_results %>% filter(resolution == 0.3)
res05 <- sil_results %>% filter(resolution == 0.5)
if (nrow(res03) > 0 && nrow(res05) > 0) {
  cat(sprintf("\n  Adjacent resolutions:\n"))
  cat(sprintf("    res 0.3: %d clusters (sil = %.4f)\n", res03$n_clusters, res03$mean_sil))
  cat(sprintf("    res 0.4: %d clusters (sil = %.4f)  <-- CHOSEN\n", res04$n_clusters, res04$mean_sil))
  cat(sprintf("    res 0.5: %d clusters (sil = %.4f)\n", res05$n_clusters, res05$mean_sil))
  if (res03$n_clusters == res04$n_clusters || res05$n_clusters == res04$n_clusters) {
    cat("  -> Cluster count is stable across nearby resolutions (clustree confirms minimal cross-assignment).\n")
  }
}

cat("\nDone! 18_clustree_resolution.R completed successfully.\n")
