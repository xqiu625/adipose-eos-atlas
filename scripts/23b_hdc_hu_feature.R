#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 23b_hdc_hu_feature.R
# Show Hdc is genuinely expressed in Hu FACS-sorted eosinophils
# Panel A: Feature plot on Hu UMAP — Hdc expression
# Panel B: Dotplot — Hdc alongside canonical eos markers by depot
# SLURM: batch, 64G, 4 CPUs, 1h
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

BASE    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures", "supplement")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# --- Load Hu data ---
cat("Loading Hu eos...\n")
hu <- readRDS("/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/hu_2026/hu_adipose_eos.rds")
if (inherits(hu[["RNA"]], "Assay5")) hu[["RNA"]] <- JoinLayers(hu[["RNA"]])
cat(sprintf("  %d cells\n", ncol(hu)))

# Find tissue/depot column
depot_col <- NULL
for (col in colnames(hu@meta.data)) {
  vals <- as.character(unique(hu@meta.data[[col]]))
  if (any(c("BAT", "eWAT", "scWAT") %in% vals)) {
    depot_col <- col
    break
  }
}
if (!is.null(depot_col)) {
  cat(sprintf("  Depot column: %s\n", depot_col))
  Idents(hu) <- depot_col
} else {
  cat("  No depot column found, using default idents\n")
}

# Ensure normalized
hu <- NormalizeData(hu, verbose = FALSE)

# Run UMAP if not present
if (!"umap" %in% Reductions(hu)) {
  cat("  Computing UMAP...\n")
  hu <- FindVariableFeatures(hu, nfeatures = 2000, verbose = FALSE)
  hu <- ScaleData(hu, verbose = FALSE)
  hu <- RunPCA(hu, npcs = 20, verbose = FALSE)
  hu <- RunUMAP(hu, dims = 1:20, verbose = FALSE)
}

# --- Panel A: Feature plot — Hdc on Hu UMAP ---
cat("Panel A: Hdc feature plot...\n")
p_feat <- FeaturePlot(hu, features = "Hdc", pt.size = 0.5,
                       order = TRUE, min.cutoff = 0) +
  scale_color_gradientn(
    colors = c("grey90", "#fee0d2", "#fc9272", "#de2d26", "#a50f15"),
    name = "Hdc") +
  ggtitle(expression(bold(italic("Hdc") ~ "Expression in Hu et al. FACS-sorted Eosinophils")),
          subtitle = "CD11b+SiglecF+ sorted, healthy C57BL/6 adipose tissue") +
  theme_umap()

# Add % expressing annotation
hdc_pct <- round(100 * mean(GetAssayData(hu, layer = "data")["Hdc", ] > 0), 1)
p_feat <- p_feat + annotate("text", x = Inf, y = -Inf,
                              label = sprintf("%.1f%% cells express Hdc", hdc_pct),
                              hjust = 1.1, vjust = -0.5, size = 4, fontface = "bold")

save_pdf(p_feat, file.path(fig_dir, "fig_23b_hdc_featureplot_hu.pdf"), 8, 7)

# --- Panel B: Dotplot — Hdc + eos markers by depot ---
cat("Panel B: Dotplot Hdc + eos markers by depot...\n")
markers <- c("Hdc", "Siglecf", "Ccr3", "Epx", "Prg2", "Il5ra")
markers <- markers[markers %in% rownames(hu)]

if (!is.null(depot_col)) {
  p_dot <- DotPlot(hu, features = markers, group.by = depot_col) +
    scale_color_gradientn(
      colors = c("#f7f7f7", "#fee0d2", "#fc9272", "#de2d26", "#a50f15"),
      name = "Avg Expr") +
    scale_size_continuous(range = c(1, 8), name = "% Expressed") +
    labs(x = "", y = "Adipose Depot",
         title = expression(bold(italic("Hdc") ~ "Among Canonical Eosinophil Markers")),
         subtitle = "Hu et al. 2026 FACS-sorted (CD11b+SiglecF+) eosinophils") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 11),
          axis.text.y = element_text(size = 11))

  save_pdf(p_dot, file.path(fig_dir, "fig_23b_hdc_dotplot_hu.pdf"), 8, 4)
}

# --- Panel C: Split feature plot by depot ---
cat("Panel C: Split feature by depot...\n")
if (!is.null(depot_col)) {
  p_split <- FeaturePlot(hu, features = "Hdc", pt.size = 0.5,
                          split.by = depot_col, order = TRUE, min.cutoff = 0) &
    scale_color_gradientn(
      colors = c("grey90", "#fee0d2", "#fc9272", "#de2d26", "#a50f15"),
      name = "Hdc") &
    theme_umap()

  save_pdf(p_split, file.path(fig_dir, "fig_23b_hdc_featureplot_hu_split.pdf"), 16, 5)
}

# --- Panel D: Bar chart — % expressing per depot ---
cat("Panel D: % expressing bar chart...\n")
if (!is.null(depot_col)) {
  depot_stats <- data.frame()
  for (dep in unique(hu@meta.data[[depot_col]])) {
    mask <- hu@meta.data[[depot_col]] == dep
    expr <- GetAssayData(hu, layer = "data")["Hdc", mask]
    depot_stats <- rbind(depot_stats, data.frame(
      depot = dep,
      pct = 100 * mean(expr > 0),
      n_cells = sum(mask),
      stringsAsFactors = FALSE
    ))
  }

  p_bar <- ggplot(depot_stats, aes(x = depot, y = pct, fill = depot)) +
    geom_bar(stat = "identity", width = 0.6, show.legend = FALSE) +
    geom_text(aes(label = sprintf("%.1f%%\n(n=%d)", pct, n_cells)),
              vjust = -0.3, size = 3.5, fontface = "bold") +
    scale_fill_manual(values = c("eWAT" = "#E63946", "scWAT" = "#457B9D", "BAT" = "#2A9D8F")) +
    scale_y_continuous(limits = c(0, max(depot_stats$pct) * 1.3)) +
    labs(x = "Adipose Depot", y = "% Cells Expressing Hdc",
         title = expression(bold(italic("Hdc") ~ "Detection in FACS-sorted Eosinophils")),
         subtitle = "Hu et al. 2026 (CD11b+SiglecF+ sorted, healthy C57BL/6)") +
    theme_pub()

  save_pdf(p_bar, file.path(fig_dir, "fig_23b_hdc_pct_bar_hu.pdf"), 5, 5)
}

cat("\n23b complete.\n")
