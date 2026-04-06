#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 31_fig2B_dotplot.R
# Regenerate Fig 2B marker dotplot matching previous submission format:
#   - C0 at bottom, C5 at top (NO rev())
#   - Top 7 markers per cluster, ordered C0→C5
#   - RdBu diverging color scale
# Input:  eos_subclustered.rds, eos_cluster_markers.csv
# Output: Fig2B_marker_dotplot.pdf
# SLURM:  batch, 64G, 4 CPUs, 1h
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures")

# --- Load data ----------------------------------------------------------------
cat("Loading eos_subclustered.rds...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
if (inherits(eos[["RNA"]], "Assay5")) eos[["RNA"]] <- JoinLayers(eos[["RNA"]])

cat(sprintf("  %d cells, clusters: %s\n", ncol(eos),
            paste(sort(unique(eos$eos_cluster)), collapse = ", ")))

# --- Load markers -------------------------------------------------------------
markers <- read.csv(file.path(BASE, "eos_cluster_markers.csv"))
cat(sprintf("  %d marker genes loaded\n", nrow(markers)))

# Top 7 markers per cluster (by avg_log2FC)
top7 <- markers %>%
  filter(p_val_adj < 0.05) %>%
  group_by(cluster) %>%
  slice_max(order_by = avg_log2FC, n = 7) %>%
  ungroup()

# Order genes: C0 markers first, then C1, ..., C5
gene_order <- unique(top7$gene[order(match(top7$cluster, paste0("C", 0:5)))])
cat(sprintf("  %d genes in dotplot\n", length(gene_order)))

# --- Set cluster order: C0 at BOTTOM, C5 at TOP (no rev!) --------------------
Idents(eos) <- factor(eos$eos_cluster, levels = paste0("C", 0:5))

# --- Plot (matching previous submission format) -------------------------------
cat("Plotting Fig 2B dotplot...\n")
pdf(file.path(fig_dir, "Fig2B_marker_dotplot.pdf"),
    width = max(14, length(gene_order) * 0.35), height = 5)
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

cat("Saved: Fig2B_marker_dotplot.pdf\n")
cat("  Gene order:", paste(gene_order, collapse = ", "), "\n")
