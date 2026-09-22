#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))
if (nzchar(Sys.getenv("RERUN_R_LIB"))) {
  .rerun_r_lib_resolved <- normalizePath(.libPaths()[1], mustWork = FALSE)
  .rerun_r_lib_expected <- normalizePath(Sys.getenv("RERUN_R_LIB"), mustWork = FALSE)
  if (!identical(.rerun_r_lib_resolved, .rerun_r_lib_expected)) {
    stop(sprintf(
      "RERUN_R_LIB audit guard: .libPaths()[1] resolved to '%s' but RERUN_R_LIB='%s' -- the shared/legacy library would take precedence over the pinned isolated library. Refusing to proceed.",
      .rerun_r_lib_resolved, .rerun_r_lib_expected))
  }
  cat("=== Pinned R environment audit guard ===\n")
  cat("R.version.string:", R.version.string, "\n")
  cat("RERUN_R_LIB:", Sys.getenv("RERUN_R_LIB"), "\n")
  cat("CMDSTAN:", Sys.getenv("CMDSTAN"), "\n")
  cat(".libPaths():\n")
  print(.libPaths())
  cat("=========================================\n")
}

# ==============================================================================
# 28_gurtner_bm_blood_heatmap.R
# Integrated heatmap + dotplot: Qiu C0-C5 vs Gurtner BM vs Gurtner Blood
# Input:  $BASE/integration/qiu_gurtner_integrated.rds
# Output: PDFs in $BASE/integration/figures/
# SLURM:  batch, 64G, 4 CPUs, 2h
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(patchwork)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
outdir <- file.path(BASE, "integration", "figures")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Load integrated object ---------------------------------------------------
cat("Loading Qiu-Gurtner integrated object...\n")
combined <- readRDS(file.path(BASE, "integration", "qiu_gurtner_integrated.rds"))
if (inherits(combined[["RNA"]], "Assay5")) combined[["RNA"]] <- JoinLayers(combined[["RNA"]])
cat(sprintf("  %d cells\n", ncol(combined)))

# --- Create tissue_group: Qiu C0-C5, Gurtner BM, Gurtner Blood ---------------
cat("Creating tissue groups...\n")
combined$tissue_group <- case_when(
  combined$dataset == "Qiu_2025_Adipose" ~ paste0("Adipose_", combined$original_cluster),
  combined$organ == "Bone Marrow"         ~ "Gurtner_BM",
  combined$organ == "Blood"               ~ "Gurtner_Blood",
  TRUE                                    ~ paste0("Gurtner_", combined$organ)
)
cat("  Group sizes:\n")
print(sort(table(combined$tissue_group), decreasing = TRUE))

# Focus on main groups: Qiu C0-C5 + BM + Blood
main_groups <- c(paste0("Adipose_C", 0:5), "Gurtner_BM", "Gurtner_Blood")
main_cells <- combined$tissue_group %in% main_groups
cat(sprintf("  Keeping %d cells in main groups\n", sum(main_cells)))
sub <- subset(combined, cells = colnames(combined)[main_cells])
sub$tissue_group <- factor(sub$tissue_group, levels = main_groups)
Idents(sub) <- "tissue_group"

# --- Colors -------------------------------------------------------------------
group_colors <- c(
  eos_colors[paste0("C", 0:5)],
  "Gurtner_BM"    = "#1D3557",
  "Gurtner_Blood" = "#457B9D"
)
names(group_colors)[1:6] <- paste0("Adipose_C", 0:5)

# --- Key eosinophil markers ---------------------------------------------------
marker_genes <- c(
  # Pan-eos
  "Siglecf", "Epx", "Prg2", "Ccr3", "Il5ra",
  # Maturation / activation
  "Cd63", "Cd9", "Itgam", "Itgb2",
  # Cytokines / effectors
  "Il4", "Il13", "Ccl6", "Retnla",
  # MHC-II (C0)
  "H2-Ab1", "H2-Aa", "Cd74",
  # Migration (C1)
  "Sell", "Cxcr4",
  # Stress (C2)
  "Atf6", "Ddit3",
  # Tissue resident (C3)
  "Cd44",
  # Progenitor
  "Kit", "Cd34"
)
marker_genes <- marker_genes[marker_genes %in% rownames(sub)]
cat(sprintf("  %d markers found in data\n", length(marker_genes)))

# ==============================================================================
# Figure 1: Dotplot — markers across Qiu C0-C5 + BM + Blood
# ==============================================================================
cat("Plotting dotplot...\n")
pdf(file.path(outdir, "gurtner_bm_blood_dotplot.pdf"),
    width = max(14, length(marker_genes) * 0.4), height = 6)
tryCatch({
  p <- DotPlot(sub, features = marker_genes, scale = TRUE,
               col.min = -1.5, col.max = 2.5) +
    scale_color_gradientn(
      colors = c("#3B4CC0", "#7092D5", "#AAC7FD", "#DDDDDD",
                 "#F7B799", "#E8744A", "#B40426"),
      name = "Scaled\nExpression") +
    scale_size_continuous(range = c(0.5, 6), name = "Percent\nExpressed") +
    labs(x = "", y = "",
         title = "Eosinophil Markers: Adipose Subtypes vs BM & Blood") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "italic"),
          axis.text.y = element_text(size = 10))
  print(p)
}, error = function(e) cat("  Dotplot failed:", e$message, "\n"))
dev.off()

# ==============================================================================
# Figure 2: Average expression heatmap
# ==============================================================================
cat("Plotting average expression heatmap...\n")

avg_expr <- AverageExpression(sub, features = marker_genes, group.by = "tissue_group",
                              assays = "RNA", layer = "data")$RNA
avg_mat <- as.matrix(avg_expr)

# Scale by gene (row)
avg_scaled <- t(scale(t(avg_mat)))

# Annotation colors
ann_col <- data.frame(
  Source = ifelse(grepl("^Adipose", colnames(avg_scaled)), "Qiu Adipose",
           ifelse(grepl("BM", colnames(avg_scaled)), "Gurtner BM", "Gurtner Blood")),
  row.names = colnames(avg_scaled)
)
ann_colors <- list(
  Source = c("Qiu Adipose" = "#E63946", "Gurtner BM" = "#1D3557", "Gurtner Blood" = "#457B9D")
)

pdf(file.path(outdir, "gurtner_bm_blood_heatmap.pdf"), width = 10, height = 10)
tryCatch({
  pheatmap(avg_scaled,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           cluster_cols = FALSE,
           cluster_rows = TRUE,
           annotation_col = ann_col,
           annotation_colors = ann_colors,
           fontsize_row = 9,
           fontsize_col = 10,
           main = "Average Expression: Adipose Subtypes vs BM & Blood",
           border_color = NA,
           breaks = seq(-2, 2, length.out = 101))
}, error = function(e) cat("  Heatmap failed:", e$message, "\n"))
dev.off()

# ==============================================================================
# Figure 3: UMAP highlighting BM vs Blood vs Adipose
# ==============================================================================
cat("Plotting UMAP by tissue source...\n")

combined$source <- case_when(
  combined$dataset == "Qiu_2025_Adipose" ~ "Adipose (Qiu)",
  combined$organ == "Bone Marrow"         ~ "Bone Marrow (Gurtner)",
  combined$organ == "Blood"               ~ "Blood (Gurtner)",
  TRUE                                    ~ paste0("Other (", combined$organ, ")")
)

pdf(file.path(outdir, "gurtner_bm_blood_umap.pdf"), width = 16, height = 6)
tryCatch({
  p1 <- DimPlot(combined, group.by = "source", pt.size = 0.3,
                cols = c("Adipose (Qiu)" = "#E63946",
                         "Bone Marrow (Gurtner)" = "#1D3557",
                         "Blood (Gurtner)" = "#457B9D")) +
    ggtitle("Tissue Source") + theme_umap()

  # Split into 3 panels
  bm_cells <- colnames(combined)[combined$organ == "Bone Marrow"]
  blood_cells <- colnames(combined)[combined$organ == "Blood"]
  adipose_cells <- colnames(combined)[combined$dataset == "Qiu_2025_Adipose"]

  p2 <- DimPlot(combined, cells.highlight = list(
    "Bone Marrow" = bm_cells, "Blood" = blood_cells, "Adipose" = adipose_cells),
    cols.highlight = c("#1D3557", "#457B9D", "#E63946"),
    sizes.highlight = 0.3, pt.size = 0.1) +
    ggtitle("Highlighted by Source") + theme_umap() + NoLegend()

  # Qiu clusters on integrated UMAP
  qiu_sub <- subset(combined, dataset == "Qiu_2025_Adipose")
  p3 <- DimPlot(qiu_sub, group.by = "original_cluster", cols = eos_colors,
                label = TRUE, pt.size = 0.3) +
    ggtitle("Qiu Eos Subtypes") + theme_umap()

  print(p1 + p3)
}, error = function(e) cat("  UMAP failed:", e$message, "\n"))
dev.off()

cat("\n========== Gurtner BM vs Blood heatmap complete ==========\n")
cat("Output directory:", outdir, "\n")
