#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

## 04b: Rename Seurat res 0.4 clusters to match manuscript C0-C6 numbering
## Based on barcode correspondence with original manuscript object
## SLURM: --mem=64G --cpus-per-task=4 --time=1:00:00 -p batch

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures")
supp_dir <- file.path(fig_dir, "supplement")
dir.create(supp_dir, showWarnings = FALSE, recursive = TRUE)

# ── Load current eos object (res 0.4, Seurat labels) ────────────────────────
cat("Loading eos_subclustered.rds...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
cat(sprintf("  %d cells, clusters: %s\n", ncol(eos),
            paste(sort(unique(eos$eos_cluster)), collapse = ", ")))
cat("  Before renaming:\n")
print(table(eos$eos_cluster))

# ── Rename map (from 04 barcode correspondence analysis) ─────────────────────
# New C0 (n=1949) -> Old C1 (69%)  => rename to C1
# New C1 (n=1739) -> Old C0 (34%)  => rename to C0
# New C2 (n=1658) -> Old C4 (44%)  => rename to C4
# New C3 (n=1156) -> Old C2 (60%)  => rename to C2
# New C4 (n=756)  -> Old C5 (82%)  => rename to C5
# New C5 (n=203)  -> Old C3 (43%)  => rename to C3
# New C6 (n=14)   -> Old C6 (100%) => keep C6

rename_map <- c(
  "C0" = "C1",
  "C1" = "C0",
  "C2" = "C4",
  "C3" = "C2",
  "C4" = "C5",
  "C5" = "C3",
  "C6" = "C6"
)

cat("\n  Renaming map (Seurat res 0.4 -> Manuscript):\n")
for (from in names(rename_map)) {
  cat(sprintf("    %s -> %s\n", from, rename_map[from]))
}

# ── Apply renaming ───────────────────────────────────────────────────────────
new_labels <- rename_map[eos$eos_cluster]
names(new_labels) <- colnames(eos)
eos$eos_cluster <- new_labels
Idents(eos) <- factor(new_labels, levels = paste0("C", 0:6))  # C6 removed below

cat("\n  After renaming:\n")
print(table(eos$eos_cluster))

# ── Remove C6 — too few cells (Cldn5/Rasip1, no mature eos markers)
n_c6 <- sum(Idents(eos) == "C6")
eos <- subset(eos, idents = "C6", invert = TRUE)
eos$eos_cluster <- droplevels(factor(eos$eos_cluster, levels = paste0("C", 0:5)))
Idents(eos) <- eos$eos_cluster
cat(sprintf("\n  Removed C6: %d cells. Remaining: %d eos (C0-C5)\n", n_c6, ncol(eos)))

# ── Re-run FindAllMarkers with renamed labels ────────────────────────────────
cat("\nFinding cluster markers with manuscript labels...\n")
markers <- FindAllMarkers(eos, only.pos = TRUE, min.pct = 0.25,
                          logfc.threshold = 0.25, verbose = FALSE)
cat(sprintf("  %d marker genes identified\n", nrow(markers)))

cat("\n  Markers per cluster:\n")
marker_counts <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  summarise(n_markers = n(), .groups = "drop")
print(as.data.frame(marker_counts))

# ── Figures ──────────────────────────────────────────────────────────────────
clusters_present <- sort(unique(eos$eos_cluster))
eos_colors_use <- eos_colors[clusters_present]

cat("\nPlotting UMAP by cluster (manuscript labels)...\n")
pdf(file.path(fig_dir, "04b_eos_umap_cluster.pdf"), width = 7, height = 6)
print(DimPlot(eos, group.by = "eos_cluster", cols = eos_colors_use, pt.size = 0.5) +
        ggtitle("Eosinophil Subclusters (Manuscript Labels)") + theme_umap())
dev.off()

cat("Plotting composition by condition...\n")
comp_df <- eos@meta.data %>%
  group_by(condition, eos_cluster) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(frac = n / sum(n))

pdf(file.path(fig_dir, "04b_eos_composition_condition.pdf"), width = 8, height = 6)
print(ggplot(comp_df, aes(x = condition, y = frac, fill = eos_cluster)) +
        geom_bar(stat = "identity", position = "stack") +
        scale_fill_manual(values = eos_colors_use) +
        labs(x = "Condition", y = "Fraction", fill = "Cluster",
             title = "Eosinophil Cluster Composition by Condition") +
        theme_pub() +
        theme(axis.text.x = element_text(angle = 45, hjust = 1)))
dev.off()

cat("Plotting top 5 marker heatmap...\n")
top5 <- markers %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 5)

pdf(file.path(fig_dir, "04b_eos_marker_heatmap.pdf"), width = 12, height = 10)
print(DoHeatmap(eos, features = top5$gene, group.colors = eos_colors_use) +
        ggtitle("Top 5 Markers per Eosinophil Cluster") +
        theme(axis.text.y = element_text(size = 8)))
dev.off()

# ── Save ─────────────────────────────────────────────────────────────────────
cat("\nSaving eos_subclustered.rds (overwrite with renamed labels)...\n")
saveRDS(eos, file.path(BASE, "eos_subclustered.rds"))

cat("Saving eos_cluster_markers.csv...\n")
write.csv(markers, file.path(BASE, "eos_cluster_markers.csv"), row.names = FALSE)

cat("\nDone! 04b_rename_clusters.R completed successfully.\n")
