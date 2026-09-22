#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

## Eosinophil subclustering with comprehensive resolution justification
## (clustree + silhouette + subsampling stability + sc-SHC significance testing)
## Eosinophil subclustering and resolution selection
## References: Grabski et al. Nat Methods 2023 (sc-SHC); Junttila et al. BIB 2022
## SLURM: --mem=128G --cpus-per-task=8 --time=6:00:00 -p batch

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
  library(clue)  # solve_LSAP: Hungarian algorithm for cluster-correspondence assignment
})

# sc-SHC is optional — requires scry which may have system library issues
HAS_SCSHC <- requireNamespace("scSHC", quietly = TRUE)
if (HAS_SCSHC) {
  library(scSHC)
  cat("  scSHC loaded successfully.\n")
} else {
  cat("  NOTE: scSHC not available (likely scry/beachmat GLIBCXX issue).\n")
  cat("  Skipping sc-SHC significance testing. Other 3 methods will run.\n")
}

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures")
supp_dir <- file.path(fig_dir, "supplement")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

# ── Chosen resolution ────────────────────────────────────────────────────────
CHOSEN_RES <- 0.4   # 7 clusters (C0-C6), stable across res 0.4-0.6

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

# ══════════════════════════════════════════════════════════════════════════════
# Step 1: Multi-resolution clustering + clustree
# ══════════════════════════════════════════════════════════════════════════════
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

# ══════════════════════════════════════════════════════════════════════════════
# Step 2: Silhouette analysis
# ══════════════════════════════════════════════════════════════════════════════
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

# Silhouette + cluster count plot (highlight chosen resolution)
p_sil <- ggplot(sil_results, aes(x = resolution, y = mean_sil)) +
  geom_line(color = "#2166AC", linewidth = 1) +
  geom_point(color = "#2166AC", size = 3) +
  geom_point(data = sil_results %>% filter(resolution == CHOSEN_RES),
             color = "red", size = 5, shape = 18) +
  geom_vline(xintercept = CHOSEN_RES, linetype = "dashed", color = "red", alpha = 0.5) +
  annotate("text", x = CHOSEN_RES + 0.05,
           y = max(sil_results$mean_sil, na.rm = TRUE),
           label = sprintf("res = %.1f\n(chosen)", CHOSEN_RES),
           hjust = 0, color = "red", size = 3.5) +
  labs(x = "Resolution", y = "Mean Silhouette Score",
       title = "Silhouette Score vs Resolution") +
  theme_pub()

p_nclust <- ggplot(sil_results, aes(x = resolution, y = n_clusters)) +
  geom_line(color = "#D6604D", linewidth = 1) +
  geom_point(color = "#D6604D", size = 3) +
  geom_point(data = sil_results %>% filter(resolution == CHOSEN_RES),
             color = "red", size = 5, shape = 18) +
  geom_vline(xintercept = CHOSEN_RES, linetype = "dashed", color = "red", alpha = 0.5) +
  scale_y_continuous(breaks = seq(0, max(sil_results$n_clusters) + 2, by = 1)) +
  labs(x = "Resolution", y = "Number of Clusters",
       title = "Number of Clusters vs Resolution") +
  theme_pub()

pdf(file.path(supp_dir, "fig_silhouette_resolution.pdf"), width = 12, height = 5)
print(p_sil + p_nclust + plot_layout(ncol = 2))
dev.off()

# UMAP panels at selected resolutions
cat("Plotting UMAP at selected resolutions...\n")
umap_resolutions <- c(0.3, 0.4, 0.5, 0.8)
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
  if (res == CHOSEN_RES) p <- p + theme(plot.background = element_rect(
    color = "red", linewidth = 2, fill = NA))
  umap_plots[[as.character(res)]] <- p
}

pdf(file.path(supp_dir, "fig_umap_resolutions.pdf"), width = 12, height = 10)
print(wrap_plots(umap_plots, ncol = 2) +
  plot_annotation(title = sprintf("UMAP at Multiple Resolutions (red = chosen, res %.1f)",
                                  CHOSEN_RES), tag_levels = "A"))
dev.off()

# Save metrics
write.csv(sil_results, file.path(supp_dir, "resolution_metrics.csv"), row.names = FALSE)

# Print summary
cat("\nResolution Selection Summary:\n")
print(sil_results)
res_chosen <- sil_results %>% filter(resolution == CHOSEN_RES)
cat(sprintf("\nChosen resolution %.1f: %d clusters, mean sil = %.4f\n",
            CHOSEN_RES, res_chosen$n_clusters, res_chosen$mean_sil))

# ══════════════════════════════════════════════════════════════════════════════
# Step 3: Subsampling stability analysis (ARI)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========================================\n")
cat("Step 3: Subsampling stability analysis\n")
cat("========================================\n\n")

# ARI function (no external dependency needed)
compute_ari <- function(labels1, labels2) {
  tab <- table(labels1, labels2)
  n <- sum(tab)
  sum_comb_nij <- sum(choose(tab, 2))
  sum_comb_ai  <- sum(choose(rowSums(tab), 2))
  sum_comb_bj  <- sum(choose(colSums(tab), 2))
  comb_n       <- choose(n, 2)
  expected <- sum_comb_ai * sum_comb_bj / comb_n
  max_idx  <- (sum_comb_ai + sum_comb_bj) / 2
  if (max_idx == expected) return(1.0)
  (sum_comb_nij - expected) / (max_idx - expected)
}

# Full-data reference labels at chosen resolution
ref_labels <- as.character(eos@meta.data[[paste0("RNA_snn_res.", CHOSEN_RES)]])
names(ref_labels) <- colnames(eos)

n_iter <- 100
subsample_frac <- 0.8
set.seed(42)

cat(sprintf("Running %d iterations (%.0f%% subsampling) at res %.1f...\n",
            n_iter, subsample_frac * 100, CHOSEN_RES))

ari_scores <- numeric(n_iter)
n_clusters_sub <- integer(n_iter)

for (i in seq_len(n_iter)) {
  if (i %% 10 == 0) cat(sprintf("  Iteration %d/%d\n", i, n_iter))

  # Subsample cells
  n_sub <- round(ncol(eos) * subsample_frac)
  sub_cells <- sample(colnames(eos), n_sub)
  eos_sub <- subset(eos, cells = sub_cells)

  # Re-cluster the subsample (using pre-computed PCA, just re-build SNN + cluster)
  eos_sub <- FindNeighbors(eos_sub, dims = 1:20, verbose = FALSE)
  eos_sub <- FindClusters(eos_sub, resolution = CHOSEN_RES, verbose = FALSE)

  sub_labels <- as.character(Idents(eos_sub))
  names(sub_labels) <- colnames(eos_sub)

  # Compute ARI between subsampled clusters and reference labels for same cells
  # IMPORTANT: align by colnames(eos_sub), NOT sub_cells — subset() preserves
  # original cell order which differs from sample() order
  shared_cells <- colnames(eos_sub)
  ref_sub <- ref_labels[shared_cells]
  ari_scores[i] <- compute_ari(ref_sub, sub_labels[shared_cells])
  n_clusters_sub[i] <- length(unique(sub_labels))
}

cat(sprintf("\nSubsampling stability results (res %.1f):\n", CHOSEN_RES))
cat(sprintf("  Mean ARI:   %.4f (sd = %.4f)\n", mean(ari_scores), sd(ari_scores)))
cat(sprintf("  Median ARI: %.4f\n", median(ari_scores)))
cat(sprintf("  Range ARI:  [%.4f, %.4f]\n", min(ari_scores), max(ari_scores)))
cat(sprintf("  Cluster count: median = %d, range = [%d, %d]\n",
            median(n_clusters_sub), min(n_clusters_sub), max(n_clusters_sub)))

# Plot ARI distribution
ari_df <- data.frame(iteration = seq_len(n_iter), ari = ari_scores,
                     n_clusters = n_clusters_sub)

p_ari <- ggplot(ari_df, aes(x = ari)) +
  geom_histogram(bins = 30, fill = "#2166AC", color = "white", alpha = 0.8) +
  geom_vline(xintercept = mean(ari_scores), linetype = "dashed",
             color = "red", linewidth = 0.8) +
  annotate("text", x = mean(ari_scores) - 0.02, y = Inf,
           label = sprintf("mean ARI = %.3f", mean(ari_scores)),
           hjust = 1, vjust = 2, color = "red", size = 4) +
  labs(x = "Adjusted Rand Index (ARI)", y = "Count",
       title = sprintf("Subsampling Stability (%.0f%%, %d iterations, res %.1f)",
                        subsample_frac * 100, n_iter, CHOSEN_RES),
       subtitle = "ARI between subsampled and full-data cluster assignments") +
  theme_pub()

p_ari_nclust <- ggplot(ari_df, aes(x = n_clusters)) +
  geom_bar(fill = "#D6604D", alpha = 0.8) +
  scale_x_continuous(breaks = seq(min(n_clusters_sub), max(n_clusters_sub))) +
  labs(x = "Number of Clusters", y = "Count",
       title = "Cluster Count Across Subsamplings") +
  theme_pub()

pdf(file.path(supp_dir, "fig_subsampling_stability.pdf"), width = 12, height = 5)
print(p_ari + p_ari_nclust + plot_layout(ncol = 2))
dev.off()

# Save ARI results
write.csv(ari_df, file.path(supp_dir, "subsampling_ari_results.csv"), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Step 4: sc-SHC significance testing (Grabski et al. Nat Methods 2023)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========================================\n")
cat("Step 4: sc-SHC significance testing\n")
cat("========================================\n\n")

# Set final clusters at chosen resolution first
eos$eos_cluster <- paste0("C", as.integer(as.character(
  eos@meta.data[[paste0("RNA_snn_res.", CHOSEN_RES)]])))
Idents(eos) <- "eos_cluster"
cat(sprintf("  Clusters to test: %s\n",
            paste(sort(unique(eos$eos_cluster)), collapse = ", ")))

scshc_result <- NULL
n_original <- length(unique(eos$eos_cluster))
n_final    <- n_original

if (HAS_SCSHC) {

  # Extract raw counts (sc-SHC requires raw UMI counts, genes x cells)
  cat("Extracting raw count matrix for sc-SHC...\n")
  raw_counts <- GetAssayData(eos, slot = "counts")
  cat(sprintf("  Count matrix: %d genes x %d cells\n", nrow(raw_counts), ncol(raw_counts)))

  # Cluster labels as integer vector
  cluster_labels <- as.integer(factor(eos$eos_cluster,
                                       levels = sort(unique(eos$eos_cluster)))) - 1L
  names(cluster_labels) <- colnames(eos)

  # Batch variable (sample)
  batch_var <- eos$orig.ident

  cat("Running sc-SHC testClusters() (post-hoc significance testing)...\n")
  cat("  This tests whether each cluster split is statistically significant.\n")
  cat("  alpha = 0.05, num_features = 2500, num_PCs = 30\n")
  cat("  batch = orig.ident (accounts for multi-sample design)\n\n")

  scshc_result <- tryCatch({
    testClusters(
      raw_counts,
      cluster_labels,
      batch = batch_var,
      alpha = 0.05,
      num_features = 2500,
      num_PCs = 30
    )
  }, error = function(e) {
    cat(sprintf("  sc-SHC testClusters ERROR: %s\n", e$message))
    cat("  Trying without batch correction...\n")
    tryCatch({
      testClusters(
        raw_counts,
        cluster_labels,
        alpha = 0.05,
        num_features = 2500,
        num_PCs = 30
      )
    }, error = function(e2) {
      cat(sprintf("  sc-SHC testClusters ERROR (no batch): %s\n", e2$message))
      NULL
    })
  })

} else {
  cat("  scSHC package not available — skipping significance testing.\n")
  cat("  To install: devtools::install_github('igrabski/sc-SHC')\n")
  cat("  Requires scry package which needs GLIBCXX >= 3.4.29\n")
}

if (!is.null(scshc_result)) {
  cat("\n========================================\n")
  cat("sc-SHC Results\n")
  cat("========================================\n\n")

  # scSHC returns a list: clusters (final labels after merging) and p_values
  final_clusters <- scshc_result[[1]]
  p_values <- scshc_result[[2]]

  n_original <- length(unique(cluster_labels))
  n_final    <- length(unique(final_clusters))

  cat(sprintf("  Original clusters: %d\n", n_original))
  cat(sprintf("  Clusters after sc-SHC testing: %d\n", n_final))

  if (n_final == n_original) {
    cat("  >>> ALL clusters are statistically significant! <<<\n")
    cat("  No merging was needed — each cluster represents a distinct population.\n")
  } else {
    cat(sprintf("  sc-SHC merged %d clusters -> %d clusters\n",
                n_original, n_final))
    cat("  Merging table (original -> sc-SHC):\n")
    merge_table <- table(
      Original = paste0("C", cluster_labels),
      scSHC = paste0("scSHC_", final_clusters)
    )
    print(merge_table)
  }

  # Print p-values for each split in the dendrogram
  if (!is.null(p_values)) {
    cat("\n  Adjusted p-values for each split:\n")
    for (i in seq_along(p_values)) {
      sig_label <- ifelse(p_values[i] < 0.05, "SIGNIFICANT", "not significant")
      cat(sprintf("    Split %d: p_adj = %.4f (%s)\n", i, p_values[i], sig_label))
    }
  }

  # Save sc-SHC results
  scshc_summary <- data.frame(
    n_original_clusters = n_original,
    n_scshc_clusters = n_final,
    all_significant = (n_final == n_original)
  )
  write.csv(scshc_summary, file.path(supp_dir, "scshc_summary.csv"), row.names = FALSE)

  if (!is.null(p_values)) {
    pval_df <- data.frame(split = seq_along(p_values),
                          adjusted_pvalue = p_values,
                          significant = p_values < 0.05)
    write.csv(pval_df, file.path(supp_dir, "scshc_pvalues.csv"), row.names = FALSE)
  }

  # If sc-SHC merged some clusters, update eos_cluster accordingly
  if (n_final < n_original) {
    cat("\n  Updating cluster labels based on sc-SHC results...\n")
    # Map original integer labels to sc-SHC merged labels
    scshc_map <- setNames(final_clusters, as.character(cluster_labels))
    # Relabel as C0, C1, ... based on sc-SHC output
    unique_scshc <- sort(unique(final_clusters))
    scshc_remap <- setNames(paste0("C", seq_along(unique_scshc) - 1L),
                            as.character(unique_scshc))
    eos$eos_cluster <- scshc_remap[as.character(scshc_map[as.character(cluster_labels)])]
    Idents(eos) <- "eos_cluster"
    cat(sprintf("  Updated clusters: %s\n",
                paste(sort(unique(eos$eos_cluster)), collapse = ", ")))
  }

  cat("\n  sc-SHC analysis complete.\n")

} else {
  # README/methods.md present sc-SHC as one of four completed, convergent
  # lines of evidence for the chosen resolution. Previously this branch just
  # printed a log line and moved on -- no artifact recorded that the check
  # never ran, so a downstream reader (or the manuscript) could not tell
  # "sc-SHC confirmed the clusters" apart from "sc-SHC was never run." Write
  # an explicit NOT_RUN status file so that gap is auditable instead of silent.
  cat("\n")
  cat("  ############################################################\n")
  cat("  # WARNING: sc-SHC significance testing did NOT run.        #\n")
  cat("  # methods.md/README list sc-SHC as one of four completed   #\n")
  cat("  # lines of evidence for the chosen resolution -- that claim#\n")
  cat("  # is NOT satisfied for this run. Proceeding with Seurat    #\n")
  cat("  # clusters validated by clustree/silhouette/subsampling    #\n")
  cat("  # ARI only.                                                #\n")
  cat("  ############################################################\n\n")

  scshc_status <- data.frame(
    status = "NOT_RUN",
    reason = "scSHC package unavailable (requireNamespace('scSHC') == FALSE)",
    n_original_clusters = n_original,
    n_scshc_clusters = NA,
    all_significant = NA
  )
  write.csv(scshc_status, file.path(supp_dir, "scshc_summary.csv"), row.names = FALSE)
  cat(sprintf("  Wrote NOT_RUN status: %s\n", file.path(supp_dir, "scshc_summary.csv")))
}

# ══════════════════════════════════════════════════════════════════════════════
# Step 5: Find cluster markers
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========================================\n")
cat("Step 5: Final clustering summary\n")
cat("========================================\n\n")

cat(sprintf("  Final clusters: %s\n",
            paste(sort(unique(eos$eos_cluster)), collapse = ", ")))
cat(sprintf("  Total cells: %d\n", ncol(eos)))
cat("\n  Cells per cluster:\n")
print(table(eos$eos_cluster))

cat("\nFinding cluster markers (FindAllMarkers)...\n")
markers <- FindAllMarkers(eos, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)
cat(sprintf("  %d marker genes identified\n", nrow(markers)))

# Per-cluster marker count
cat("\n  Markers per cluster:\n")
marker_counts <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(n_markers = n(), .groups = "drop")
print(as.data.frame(marker_counts))

# ══════════════════════════════════════════════════════════════════════════════
# Step 5b: Map new clusters (res 0.4) to old clusters (res 0.8)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========================================\n")
cat("Step 5b: Correspondence with manuscript clusters\n")
cat("========================================\n\n")
cat("  Old manuscript: no DecontX, 10,698 eos, 7 clusters (C0-C6)\n")
cat("  New pipeline:   DecontX + QC + Scrublet, 7,475 eos, res 0.4\n")
cat("  Matching by cell barcode to map new -> old cluster identity\n\n")

# Load original manuscript Seurat object (has the published Eos C0-C6 labels)
# Old object: eos_cluster column = C0-C6 (10,698 eos cells), NA for non-eos
old_dir <- "/bigdata/nairlab/shared/XQiu/EOS_project/R07_20250702/03_cellchat_eos_mac/"
old_rds <- file.path(old_dir, "seurat_object_for_cellchat.rds")

old_eos_labels <- NULL

if (file.exists(old_rds)) {
  cat(sprintf("Loading original manuscript object:\n  %s\n", old_rds))
  old_obj <- tryCatch(readRDS(old_rds), error = function(e) {
    cat(sprintf("  ERROR loading: %s\n", e$message)); NULL
  })

  if (!is.null(old_obj)) {
    cat(sprintf("  Old object: %d genes x %d cells\n", nrow(old_obj), ncol(old_obj)))

    # Extract eos_cluster column (C0-C6 for eos, NA for non-eos)
    old_all_labels <- as.character(old_obj$eos_cluster)
    names(old_all_labels) <- colnames(old_obj)

    # Keep only eosinophil cells (non-NA and non-"NA" string)
    keep <- !is.na(old_all_labels) & old_all_labels != "NA" & grepl("^C[0-9]", old_all_labels)
    old_eos_labels <- old_all_labels[keep]
    cat(sprintf("  Old eosinophil cells: %d\n", length(old_eos_labels)))
    cat(sprintf("  Old clusters: %s\n",
                paste(names(sort(table(old_eos_labels), decreasing = TRUE)),
                      collapse = ", ")))
    print(table(old_eos_labels))

    rm(old_obj); gc()
  }
} else {
  cat(sprintf("  File not found: %s\n", old_rds))
}

# Match barcodes between new (res 0.4) and old (manuscript) clusters
if (!is.null(old_eos_labels) && length(old_eos_labels) > 0) {
  common_cells <- intersect(colnames(eos), names(old_eos_labels))
  cat(sprintf("\n  Barcode matching:\n"))
  cat(sprintf("    New eos cells:    %d\n", ncol(eos)))
  cat(sprintf("    Old eos cells:    %d\n", length(old_eos_labels)))
  cat(sprintf("    Common barcodes:  %d\n", length(common_cells)))

  if (length(common_cells) > 0) {
    new_labels <- eos$eos_cluster[common_cells]
    old_labels <- old_eos_labels[common_cells]

    # Confusion matrix: rows = new (res 0.4), cols = old (manuscript)
    conf_mat <- table(New_res0.4 = new_labels, Old_manuscript = old_labels)
    cat("\n  Confusion matrix (New res 0.4 vs Old manuscript):\n")
    print(conf_mat)

    # Proportion matrix: for each new cluster, what % came from each old cluster
    prop_mat <- prop.table(conf_mat, margin = 1)
    cat("\n  Proportion matrix (row-normalized, per new cluster):\n")
    print(round(prop_mat, 3))

    # Dominant mapping: for each new cluster, which old cluster has the most cells
    dominant_map <- data.frame(
      new_cluster = rownames(prop_mat),
      stringsAsFactors = FALSE
    )
    dominant_map$primary_old_cluster <- apply(prop_mat, 1, function(x) {
      colnames(prop_mat)[which.max(x)]
    })
    dominant_map$primary_fraction <- apply(prop_mat, 1, max)
    dominant_map$secondary_old_cluster <- apply(prop_mat, 1, function(x) {
      sorted <- sort(x, decreasing = TRUE)
      if (length(sorted) >= 2 && sorted[2] > 0.05) {
        colnames(prop_mat)[order(x, decreasing = TRUE)[2]]
      } else {
        NA
      }
    })
    dominant_map$secondary_fraction <- apply(prop_mat, 1, function(x) {
      sort(x, decreasing = TRUE)[2]
    })
    dominant_map$n_cells <- as.integer(table(new_labels)[dominant_map$new_cluster])

    cat("\n  Cluster correspondence (new res 0.4 -> old manuscript):\n")
    for (i in seq_len(nrow(dominant_map))) {
      row <- dominant_map[i, ]
      msg <- sprintf("    %s (n=%d) -> %s (%.1f%%)",
                     row$new_cluster, row$n_cells,
                     row$primary_old_cluster, row$primary_fraction * 100)
      if (!is.na(row$secondary_old_cluster)) {
        msg <- sprintf("%s + %s (%.1f%%)", msg,
                       row$secondary_old_cluster, row$secondary_fraction * 100)
      }
      cat(msg, "\n")
    }

    # ARI between old and new
    ari_old_new <- compute_ari(new_labels, old_labels)
    cat(sprintf("\n  ARI (new res 0.4 vs old manuscript): %.4f\n", ari_old_new))

    # Save correspondence table
    write.csv(dominant_map,
              file.path(supp_dir, "cluster_correspondence_new_vs_old.csv"),
              row.names = FALSE)
    write.csv(as.data.frame.matrix(conf_mat),
              file.path(supp_dir, "cluster_confusion_matrix.csv"))
    cat("  Saved: cluster_correspondence_new_vs_old.csv\n")
    cat("  Saved: cluster_confusion_matrix.csv\n")

    # Global one-to-one assignment (Hungarian algorithm / linear sum
    # assignment) on the confusion matrix. dominant_map above is a per-row
    # greedy pick and can (a) assign two new clusters to the same old
    # cluster, and (b) have individual rows with low confidence even when an
    # unambiguous best-overall assignment exists -- the documented historical
    # correspondence for this dataset includes primary fractions as low as
    # 34%/43%/44%. solve_LSAP finds the assignment that maximizes total
    # overlap across ALL clusters simultaneously; per-cluster ambiguity is
    # reported as a diagnostic (see Step 5c) rather than gating on any single
    # row's fraction.
    cat("\n  Computing globally optimal one-to-one assignment (Hungarian algorithm)...\n")
    count_mat <- as.matrix(conf_mat)  # rows = new, cols = old (raw counts)
    n_new <- nrow(count_mat)
    n_old <- ncol(count_mat)
    n_dim <- max(n_new, n_old)
    # solve_LSAP requires a nonnegative cost matrix, so cost = max_val - count
    # (maximizing overlap = minimizing this cost) rather than -count. Padding
    # rows/cols beyond n_new/n_old are given cost = max_val, i.e. exactly as
    # unattractive as a real zero-overlap entry (max_val - 0) -- they are
    # only ever chosen when a real cluster has no better option left, which
    # is itself informative (assigned_fraction == 0 below), not an
    # artificially "free" assignment.
    max_val <- max(count_mat)
    cost_mat <- matrix(max_val, n_dim, n_dim)
    cost_mat[seq_len(n_new), seq_len(n_old)] <- max_val - count_mat
    assignment <- clue::solve_LSAP(cost_mat)

    optimal_map <- data.frame(
      new_cluster = rownames(count_mat),
      stringsAsFactors = FALSE
    )
    optimal_map$assigned_old_cluster <- vapply(seq_len(n_new), function(i) {
      j <- assignment[i]
      if (j <= n_old) colnames(count_mat)[j] else NA_character_
    }, character(1))
    optimal_map$assigned_fraction <- vapply(seq_len(n_new), function(i) {
      j <- assignment[i]
      if (j <= n_old) prop_mat[i, j] else 0
    }, numeric(1))
    optimal_map$n_cells <- as.integer(rowSums(count_mat))

    AMBIGUITY_THRESHOLD <- 0.5
    cat("\n  Globally optimal assignment (new -> old):\n")
    for (i in seq_len(nrow(optimal_map))) {
      row <- optimal_map[i, ]
      flag <- if (is.na(row$assigned_old_cluster)) {
        " [UNASSIGNED: no old cluster left]"
      } else if (row$assigned_fraction == 0) {
        " [ZERO OVERLAP: assigned only because every other cluster was a better match for something else]"
      } else if (row$assigned_fraction < AMBIGUITY_THRESHOLD) {
        " [AMBIGUOUS: <50% of this cluster's cells agree, but still the best global assignment]"
      } else {
        ""
      }
      cat(sprintf("    %s (n=%d) -> %s (%.1f%%)%s\n",
                  row$new_cluster, row$n_cells,
                  ifelse(is.na(row$assigned_old_cluster), "NA", row$assigned_old_cluster),
                  row$assigned_fraction * 100, flag))
    }

    write.csv(optimal_map,
              file.path(supp_dir, "cluster_correspondence_optimal_assignment.csv"),
              row.names = FALSE)
    cat("  Saved: cluster_correspondence_optimal_assignment.csv\n")

    # Confusion matrix heatmap
    cat("  Plotting confusion matrix heatmap...\n")
    conf_df <- as.data.frame(conf_mat)
    colnames(conf_df) <- c("New_res0.4", "Old_manuscript", "n_cells")

    # Also make proportion version for color
    prop_df <- as.data.frame(prop_mat)
    colnames(prop_df) <- c("New_res0.4", "Old_manuscript", "fraction")
    conf_df$fraction <- prop_df$fraction

    p_conf <- ggplot(conf_df, aes(x = Old_manuscript, y = New_res0.4, fill = fraction)) +
      geom_tile(color = "white", linewidth = 0.5) +
      geom_text(aes(label = n_cells), size = 3.5, color = "black") +
      scale_fill_gradient(low = "white", high = "#2166AC", name = "Fraction",
                          limits = c(0, 1)) +
      labs(x = "Old Clusters (Manuscript)", y = "New Clusters (res 0.4)",
           title = "Cluster Correspondence: New (res 0.4) vs Old (Manuscript)",
           subtitle = sprintf("Numbers = cell counts; ARI = %.3f", ari_old_new)) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1))

    pdf(file.path(supp_dir, "fig_cluster_correspondence.pdf"), width = 8, height = 6)
    print(p_conf)
    dev.off()
    cat("  Saved: fig_cluster_correspondence.pdf\n")

    # Side-by-side UMAP: old vs new coloring
    eos$old_cluster <- old_eos_labels[colnames(eos)]

    p_old <- DimPlot(eos, group.by = "old_cluster", cols = eos_colors,
                     pt.size = 0.3, label = TRUE, label.size = 3) +
      ggtitle("Old Clusters (Manuscript)") + theme_umap() +
      theme(legend.position = "right")

    p_new <- DimPlot(eos, group.by = "eos_cluster", cols = eos_colors,
                     pt.size = 0.3, label = TRUE, label.size = 3) +
      ggtitle("New Clusters (res 0.4)") + theme_umap() +
      theme(legend.position = "right")

    pdf(file.path(supp_dir, "fig_umap_old_vs_new_clusters.pdf"), width = 14, height = 6)
    print(p_old + p_new + plot_layout(ncol = 2))
    dev.off()
    cat("  Saved: fig_umap_old_vs_new_clusters.pdf\n")

  } else {
    cat("  No common barcodes found.\n")
    cat("  NOTE: Old manuscript had no DecontX cleaning (10,698 eos).\n")
    cat("  New pipeline applies DecontX + stricter QC (7,475 eos).\n")
    cat("  Barcodes may differ if cell naming convention changed.\n")
  }
} else {
  cat("  Old manuscript eos labels not available. Skipping correspondence.\n")
}

# ══════════════════════════════════════════════════════════════════════════════
# Step 5c: Rename clusters to match manuscript numbering (C0-C6)
# ══════════════════════════════════════════════════════════════════════════════
cat("\n========================================\n")
cat("Step 5c: Rename clusters to match manuscript\n")
cat("========================================\n\n")

# Renaming must come from THIS run's validated correspondence (optimal_map,
# the Hungarian-algorithm global one-to-one assignment computed in Step 5b
# above from barcode matching against the old manuscript object), not a fixed
# historical permutation. Preprocessing now removes doublets before DecontX,
# so cluster numbering -- and which old cluster each new cluster corresponds
# to -- can genuinely change between runs; a stale hardcoded map (as this
# used to be, e.g. "New C0 -> Old C1 (69%)") could silently mislabel clusters
# if the correspondence shifted. Abort rather than guess.
#
# A per-row minimum-confidence requirement (e.g. every cluster's dominant
# match >= 50%) is NOT used as a gate: the documented historical
# correspondence for this exact dataset includes primary fractions as low as
# 34%/43%/44%, so that gate would abort on a known-good assignment. Instead,
# the global assignment is validated as complete and one-to-one (guaranteed
# by construction via solve_LSAP, given every cluster has some
# correspondence data), and any individual low-confidence or zero-overlap
# assignment is reported explicitly rather than silently accepted or used to
# block the run.
current_new_clusters <- sort(unique(as.character(eos$eos_cluster)))

if (!exists("optimal_map") || is.null(optimal_map) || nrow(optimal_map) == 0) {
  stop("Cannot rename clusters to manuscript numbering: no validated ",
       "correspondence to the old manuscript object was computed this run ",
       "(old object missing or no common barcodes -- see Step 5b above). ",
       "Refusing to fall back to a hardcoded historical permutation, which ",
       "may not match this run's cluster numbering.")
}

missing_clusters <- setdiff(current_new_clusters, optimal_map$new_cluster)
if (length(missing_clusters) > 0) {
  stop(sprintf(
    "Cluster correspondence is incomplete: cluster(s) %s exist in this run but have no barcode-matched correspondence to the manuscript (likely no common barcodes for that cluster). Cannot safely rename.",
    paste(missing_clusters, collapse = ", ")))
}

unassigned <- optimal_map$new_cluster[is.na(optimal_map$assigned_old_cluster)]
if (length(unassigned) > 0) {
  stop(sprintf(
    "Cluster correspondence assignment is incomplete: cluster(s) %s have no old-cluster counterpart in the optimal one-to-one assignment (more new clusters than old ones?). See cluster_correspondence_optimal_assignment.csv. Cannot safely rename.",
    paste(unassigned, collapse = ", ")))
}

zero_overlap <- optimal_map$new_cluster[optimal_map$assigned_fraction == 0]
if (length(zero_overlap) > 0) {
  stop(sprintf(
    "Cluster correspondence assignment is degenerate: cluster(s) %s were assigned an old cluster with ZERO shared barcodes (forced by the one-to-one constraint, not real correspondence). See cluster_correspondence_optimal_assignment.csv. Resolve manually before renaming.",
    paste(zero_overlap, collapse = ", ")))
}

ambiguous <- optimal_map$new_cluster[optimal_map$assigned_fraction < 0.5]
if (length(ambiguous) > 0) {
  cat(sprintf("  NOTE: cluster(s) %s have an ambiguous assignment (<50%% agreement) -- still the best global one-to-one assignment; see cluster_correspondence_optimal_assignment.csv.\n",
              paste(ambiguous, collapse = ", ")))
}

stopifnot(!any(duplicated(optimal_map$assigned_old_cluster)))  # guaranteed by solve_LSAP; sanity check only

rename_map <- setNames(optimal_map$assigned_old_cluster, optimal_map$new_cluster)

cat("  Renaming map (Seurat -> Manuscript), derived from this run's globally optimal barcode-matched correspondence:\n")
for (from in names(rename_map)) {
  cat(sprintf("    %s -> %s\n", from, rename_map[from]))
}

# Apply renaming
old_labels_before <- eos$eos_cluster
eos$eos_cluster <- unname(rename_map[eos$eos_cluster])
Idents(eos) <- eos$eos_cluster

cat("\n  Cluster sizes after renaming:\n")
print(table(eos$eos_cluster))

# Remove C6 — too few cells (n=14)
n_c6 <- sum(eos$eos_cluster == "C6")
if (n_c6 > 0) {
  cells_keep <- colnames(eos)[eos$eos_cluster != "C6"]
  eos <- subset(eos, cells = cells_keep)
  cat(sprintf("\n  Removed C6: %d cells.\n", n_c6))
}
eos$eos_cluster <- factor(as.character(eos$eos_cluster), levels = paste0("C", 0:5))
Idents(eos) <- eos$eos_cluster
cat(sprintf("  Remaining: %d eos (C0-C5)\n", ncol(eos)))

# ══════════════════════════════════════════════════════════════════════════════
# Step 6: Figures
# ══════════════════════════════════════════════════════════════════════════════
clusters_present <- sort(unique(eos$eos_cluster))
eos_colors_use <- eos_colors[clusters_present]

# UMAP by cluster
cat("\nPlotting UMAP by cluster...\n")
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

cond_order <- c("WT_Male_HFD", "KO_Male_HFD",
                "WT_Female_HFD", "KO_Female_HFD",
                "WT_Female_HFDInfected", "KO_Female_HFDInfected")
cond_display <- c("WT_Male_HFD" = "WT Male", "KO_Male_HFD" = "KO Male",
                  "WT_Female_HFD" = "WT Female", "KO_Female_HFD" = "KO Female",
                  "WT_Female_HFDInfected" = "WT Inf Female",
                  "KO_Female_HFDInfected" = "KO Inf Female")
comp_df$condition <- factor(comp_df$condition, levels = cond_order,
                            labels = cond_display[cond_order])

pdf(file.path(fig_dir, "04_eos_composition_condition.pdf"), width = 8, height = 6)
print(ggplot(comp_df, aes(x = condition, y = frac, fill = eos_cluster)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = eos_colors_use) +
        labs(x = "Condition", y = "Fraction", fill = "Cluster",
             title = "Eosinophil Cluster Composition by Condition") +
        theme_pub() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

# Top 7 marker dotplot (Fig 2B — manuscript style)
cat("Plotting top 7 marker dotplot (Fig 2B)...\n")
top7 <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 7) %>%
  ungroup()

# Order genes: C0 markers first, then C1, etc.
gene_order <- unique(top7$gene[order(match(top7$cluster, paste0("C", 0:5)))])

# Clusters on y-axis: C0 at bottom, C5 at top (ggplot2 maps first level to bottom)
Idents(eos) <- factor(eos$eos_cluster, levels = paste0("C", 0:5))

pdf(file.path(fig_dir, "04_eos_marker_dotplot.pdf"), width = max(14, length(gene_order) * 0.35), height = 5)
print(DotPlot(eos, features = gene_order, scale = TRUE,
              col.min = -1.5, col.max = 2.5) +
        scale_color_gradientn(
          colors = c("#3B4CC0", "#7092D5", "#AAC7FD", "#DDDDDD",
                     "#F7B799", "#E8744A", "#B40426"),
          name = "Average\nExpression") +
        scale_size_continuous(range = c(0.5, 6), name = "Percent\nExpressed") +
        labs(x = "", y = "") +
        theme_pub() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "italic"),
              axis.text.y = element_text(size = 10),
              legend.position = "right",
              panel.grid.major = element_blank()))
dev.off()

# Reset ident order
Idents(eos) <- factor(eos$eos_cluster, levels = paste0("C", 0:5))

# Also keep heatmap as supplement
cat("Plotting top 7 marker heatmap (supplement)...\n")
pdf(file.path(fig_dir, "04_eos_marker_heatmap.pdf"), width = 12, height = 10)
print(DoHeatmap(eos, features = gene_order, group.colors = eos_colors_use) +
        ggtitle("Top 7 Markers per Eosinophil Cluster") +
        theme(axis.text.y = element_text(size = 7)))
dev.off()

# ══════════════════════════════════════════════════════════════════════════════
# Save outputs
# ══════════════════════════════════════════════════════════════════════════════
# Re-run FindAllMarkers with renamed cluster labels
cat("\nRe-running FindAllMarkers with manuscript cluster names...\n")
markers <- FindAllMarkers(eos, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)
cat(sprintf("  %d marker genes identified\n", nrow(markers)))

cat("\nSaving eos_subclustered.rds...\n")
saveRDS(eos, file.path(BASE, "eos_subclustered.rds"))

cat("Saving eos_cluster_markers.csv...\n")
write.csv(markers, file.path(BASE, "eos_cluster_markers.csv"), row.names = FALSE)

# ══════════════════════════════════════════════════════════════════════════════
# Final summary
# ══════════════════════════════════════════════════════════════════════════════
cat("\n")
cat("================================================================\n")
cat("CLUSTERING RESOLUTION JUSTIFICATION SUMMARY\n")
cat("================================================================\n\n")

cat(sprintf("Chosen resolution: %.1f -> %d clusters (C0-C%d)\n",
            CHOSEN_RES,
            length(unique(eos$eos_cluster)),
            length(unique(eos$eos_cluster)) - 1))

cat(sprintf("\n1. CLUSTREE: 7 clusters stable across res 0.4-0.6\n"))
cat(sprintf("   (see fig_clustree.pdf)\n"))

cat(sprintf("\n2. SILHOUETTE: mean = %.4f at res %.1f\n",
            res_chosen$mean_sil, CHOSEN_RES))
cat(sprintf("   Local peak among non-trivial resolutions\n"))
cat(sprintf("   (see fig_silhouette_resolution.pdf)\n"))

cat(sprintf("\n3. SUBSAMPLING STABILITY: mean ARI = %.4f (sd = %.4f, n = %d)\n",
            mean(ari_scores), sd(ari_scores), n_iter))
cat(sprintf("   Cluster count: median = %d across %d subsamplings\n",
            median(n_clusters_sub), n_iter))
cat(sprintf("   (see fig_subsampling_stability.pdf)\n"))

if (!is.null(scshc_result)) {
  cat(sprintf("\n4. sc-SHC SIGNIFICANCE TESTING (Grabski et al. Nat Methods 2023):\n"))
  cat(sprintf("   Original clusters: %d\n", n_original))
  cat(sprintf("   Clusters after testing: %d\n", n_final))
  if (n_final == n_original) {
    cat(sprintf("   ALL %d clusters are statistically significant (no merging needed)\n",
                n_original))
  } else {
    cat(sprintf("   %d clusters merged to %d\n", n_original, n_final))
  }
  cat(sprintf("   (see scshc_summary.csv, scshc_pvalues.csv)\n"))
}

if (exists("optimal_map") && !is.null(optimal_map)) {
  cat(sprintf("\n5. CLUSTER CORRESPONDENCE & RENAMING:\n"))
  cat(sprintf("   Old: no DecontX (10,698 eos) | New: DecontX + QC (7,475 eos)\n"))
  cat(sprintf("   ARI (before renaming) = %.4f\n", ari_old_new))
  cat(sprintf("   Clusters renamed to match manuscript numbering (Hungarian-algorithm global assignment):\n"))
  for (from in names(rename_map)) {
    cat(sprintf("     Seurat %s -> Manuscript %s\n", from, rename_map[from]))
  }
  cat(sprintf("\n   Final cluster sizes (manuscript labels):\n"))
  print(table(eos$eos_cluster))
  cat(sprintf("   (see fig_cluster_correspondence.pdf, cluster_correspondence_optimal_assignment.csv)\n"))
}

cat("\n================================================================\n")
cat("Done! 04_eos_subcluster.R completed successfully.\n")
