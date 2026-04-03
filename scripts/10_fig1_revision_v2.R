#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 10_fig1_revision_v2.R
# Revised Figure 1 panels addressing Rebecca's requests & reviewer comments
#
# Changes from v1 (10_fig1_global.R):
#   Fig 1B: All filled circles, different color per group (Rev 2 Comment 3)
#           Both immune cell AND eos-only PCA/MDS
#   Fig 1C: Updated UMAPs in PDF (immune only + immune with gray background)
#   Fig 1D: Dot plot of canonical immune cell markers (replicating manuscript)
#   Fig 1E-H: Volcano + DEG heatmap (Hu-style, Rev 2 Comment 5) + GO enrichment
#
# Input:  all_cells_annotated.rds, harmony_tier1_clustered.rds, DEG CSVs
# Output: PDF figures in revision/figures/fig1/
# SLURM:  batch, 128G, 4 CPUs, 4h
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(tidyr)
  library(patchwork)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

# --- Paths ---
PNAS     <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision"
BASE     <- file.path(PNAS, "atlas")
deg_dir  <- file.path(BASE, "deg")  # pseudobulk DESeq2 results
out_dir  <- file.path(PNAS, "revision/figures/fig1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# Fig 1B: Pseudobulk PCA/MDS — all filled circles, distinct colors
# Both immune cells and eos-only (Rev 2 Comment 3)
# =============================================================================
cat("=== Fig 1B: Pseudobulk PCA/MDS (all circles, distinct colors) ===\n")

obj <- readRDS(file.path(BASE, "all_cells_annotated.rds"))
cat(sprintf("  Loaded: %d genes x %d cells\n", nrow(obj), ncol(obj)))

# Handle Seurat v5 layers
if (inherits(obj[["RNA"]], "Assay5")) {
  obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
}

# PCA/MDS helper: DESeq2 rlog -> top 500 variable genes -> PCA + MDS
run_fig1b <- function(sobj, label, prefix) {
  cat(sprintf("\n--- Fig 1B: %s (%d cells) ---\n", label, ncol(sobj)))

  # Aggregate raw counts per sample
  counts_mat <- GetAssayData(sobj, layer = "counts")
  samples <- unique(sobj$orig.ident)

  pb_list <- lapply(samples, function(s) {
    cells <- which(sobj$orig.ident == s)
    Matrix::rowSums(counts_mat[, cells, drop = FALSE])
  })
  pb_counts <- do.call(cbind, pb_list)
  colnames(pb_counts) <- samples
  pb_counts <- pb_counts[rowSums(pb_counts) > 0, ]

  # Sample metadata
  col_data <- data.frame(
    sample    = samples,
    condition = sapply(samples, function(s) {
      unique(sobj$condition[sobj$orig.ident == s])[1]
    }),
    row.names = samples
  )

  # DESeq2 rlog
  dds <- DESeqDataSetFromMatrix(countData = round(pb_counts),
                                colData   = col_data,
                                design    = ~ 1)
  rlog_mat <- assay(rlog(dds, blind = TRUE))

  # Top 500 variable genes
  gene_vars <- apply(rlog_mat, 1, var)
  ntop <- min(500, sum(gene_vars > 0))
  top_genes <- names(sort(gene_vars, decreasing = TRUE))[1:ntop]

  # PCA
  pca_res <- prcomp(t(rlog_mat[top_genes, ]), scale. = FALSE)
  pct_var <- round(100 * pca_res$sdev^2 / sum(pca_res$sdev^2), 1)
  pca_df <- data.frame(
    Dim1 = pca_res$x[, 1], Dim2 = pca_res$x[, 2],
    sample = rownames(pca_res$x),
    condition = col_data[rownames(pca_res$x), "condition"]
  )

  # MDS (1 - Spearman correlation)
  cor_mat <- cor(rlog_mat[top_genes, ], method = "spearman")
  mds_res <- cmdscale(as.dist(1 - cor_mat), k = 2, eig = TRUE)
  mds_df <- data.frame(
    Dim1 = mds_res$points[, 1], Dim2 = mds_res$points[, 2],
    sample = rownames(mds_res$points),
    condition = col_data[rownames(mds_res$points), "condition"]
  )

  # Plot helper: ALL FILLED CIRCLES (shape = 16), different color per group
  make_fig1b <- function(df, xlab, ylab, title_text) {
    present <- intersect(names(cond_colors), unique(df$condition))
    df$condition <- factor(df$condition, levels = names(cond_colors))

    ggplot(df, aes(x = Dim1, y = Dim2, color = condition)) +
      geom_point(shape = 16, size = 5) +
      scale_color_manual(values = cond_colors[present],
                         labels = cond_labels[present],
                         name = "Group") +
      geom_text_repel(aes(label = sample), size = 2.5, color = "grey30",
                      max.overlaps = 20, segment.color = "grey70") +
      geom_hline(yintercept = 0, color = "grey80", linewidth = 0.3) +
      geom_vline(xintercept = 0, color = "grey80", linewidth = 0.3) +
      labs(title = title_text, x = xlab, y = ylab) +
      theme_pub() +
      theme(legend.position = "right",
            plot.title = element_text(face = "bold", size = 14))
  }

  # PCA
  p_pca <- make_fig1b(pca_df,
                       sprintf("PC1 (%.1f%%)", pct_var[1]),
                       sprintf("PC2 (%.1f%%)", pct_var[2]),
                       sprintf("Pseudobulk PCA — %s", label))
  save_pdf(p_pca, file.path(out_dir, sprintf("fig1B_%s_pca.pdf", prefix)), 9, 7)

  # MDS
  p_mds <- make_fig1b(mds_df,
                       "MDS Dimension 1", "MDS Dimension 2",
                       sprintf("Pseudobulk MDS — %s", label))
  save_pdf(p_mds, file.path(out_dir, sprintf("fig1B_%s_mds.pdf", prefix)), 9, 7)
}

# Immune cells PCA/MDS
immune_types <- c("Eosinophils", "Macrophages", "Dendritic cells",
                  "B cells", "NK cells", "Neutrophils")
immune <- subset(obj, celltype %in% immune_types)
run_fig1b(immune, "Immune Cells", "immune")
rm(immune); gc()

# Eosinophils only PCA/MDS
eos <- subset(obj, celltype == "Eosinophils")
run_fig1b(eos, "Eosinophils", "eos")
rm(eos); gc()

# =============================================================================
# Fig 1C: Updated UMAPs in PDF
# (1) Immune cells only  (2) Immune cells + others in gray
# =============================================================================
cat("\n=== Fig 1C: UMAPs ===\n")

ct_colors <- broad_celltype_colors[immune_types]

# (1) Immune cells only
immune_cells <- colnames(obj)[obj$celltype %in% immune_types]
sobj_immune <- subset(obj, cells = immune_cells)
sobj_immune$celltype <- factor(sobj_immune$celltype, levels = immune_types)

p_immune <- DimPlot(sobj_immune, group.by = "celltype", label = TRUE,
                    repel = TRUE, pt.size = 0.3, cols = ct_colors) +
  ggtitle("Immune Cells") + theme_umap()
save_pdf(p_immune, file.path(out_dir, "fig1C_umap_immune_only.pdf"), 10, 8)

# (2) All cells: immune colored, others in gray
all_ct <- sort(unique(obj$celltype))
all_colors <- setNames(rep("grey85", length(all_ct)), all_ct)
all_colors[names(ct_colors)] <- ct_colors

# Reorder cells so immune are plotted on top
obj$plot_order <- ifelse(obj$celltype %in% immune_types, 1, 0)
cell_order <- order(obj$plot_order)

p_all <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE,
                 pt.size = 0.2, cols = all_colors,
                 order = rev(immune_types)) +
  ggtitle("All Cells (non-immune in gray)") + theme_umap()
save_pdf(p_all, file.path(out_dir, "fig1C_umap_all_cells_immune_highlight.pdf"), 10, 8)

rm(sobj_immune); gc()

# =============================================================================
# Fig 1D: Dot plot of canonical immune cell marker genes
# Replicating original manuscript Fig 1D style
# =============================================================================
cat("\n=== Fig 1D: Immune cell marker dot plot ===\n")

immune_sub <- subset(obj, celltype %in% immune_types)

# Set cell type order to match manuscript (bottom to top)
immune_sub$celltype <- factor(immune_sub$celltype,
  levels = rev(c("Dendritic cells", "B cells",
                 "NK cells", "Neutrophils", "Macrophages", "Eosinophils")))

# Canonical marker genes — manuscript set + Nkg7 for NK cells
marker_genes <- c("Siglecf", "Ccr3", "Il5ra",             # Eosinophils
                   "Mertk", "Mrc1", "Fcgr1", "Retnla",    # Macrophages
                   "S100a8",                                # Neutrophils
                   "Nkg7",                                  # NK cells
                   "Cd19",                                  # B cells
                   "Itgax")                                 # Dendritic cells

# Filter to genes present in the data
marker_genes <- marker_genes[marker_genes %in% rownames(immune_sub)]
cat(sprintf("  Using %d / %d marker genes\n", length(marker_genes),
            length(c("Siglecf","Ccr3","Il5ra","Mertk","Mrc1","Fcgr1","Retnla",
                     "S100a8","Nkg7","Cd19","Itgax"))))

p_dot <- DotPlot(immune_sub, features = marker_genes, group.by = "celltype",
                 scale = TRUE) +
  scale_color_gradientn(
    colors = c("#3B4CC0", "#7092D5", "#AAC7FD", "#DDDDDD",
               "#F7B799", "#E8744A", "#B40426"),
    name = "Average\nExpression") +
  scale_size_continuous(range = c(1, 8), name = "Percent\nExpressed") +
  labs(x = "", y = "") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 10, face = "italic"),
        axis.text.y = element_text(size = 11),
        legend.position = "right",
        panel.grid.major = element_line(color = "grey90", linewidth = 0.3))

save_pdf(p_dot, file.path(out_dir, "fig1D_dotplot_markers.pdf"), 12, 6)
cat("  Saved: fig1D_dotplot_markers.pdf\n")

rm(immune_sub); gc()

# =============================================================================
# Fig 1E-H: Volcano plots + DEG heatmaps (Hu-style) + GO enrichment
# =============================================================================
cat("\n=== Fig 1E-H: Volcanos + Heatmaps + GO ===\n")

# Load Tier 1 object for heatmaps
tier1 <- readRDS(file.path(PNAS, "revision/harmony_tier1_clustered.rds"))
if (inherits(tier1[["RNA"]], "Assay5")) {
  tier1[["RNA"]] <- JoinLayers(tier1[["RNA"]])
}

# Helper: make volcano (same as v1)
ko_exclude <- c("Retnla", "Igha")

make_volcano <- function(deg, title, pdf_path, label_left, label_right,
                         exclude_genes = NULL) {
  deg$gene <- rownames(deg)
  if (!is.null(exclude_genes)) deg <- deg[!deg$gene %in% exclude_genes, ]

  deg$sig <- "NS"
  deg$sig[deg$p_val_adj < 0.05 & deg$avg_log2FC > 0.25]  <- "Up"
  deg$sig[deg$p_val_adj < 0.05 & deg$avg_log2FC < -0.25] <- "Down"
  deg$sig <- factor(deg$sig, levels = c("Down", "NS", "Up"))

  top_pos <- head(deg[deg$avg_log2FC > 0, ][order(deg$p_val_adj[deg$avg_log2FC > 0]), ], 10)
  top_neg <- head(deg[deg$avg_log2FC < 0, ][order(deg$p_val_adj[deg$avg_log2FC < 0]), ], 10)
  label_df <- rbind(top_pos, top_neg)

  deg$neglog10p <- -log10(pmax(deg$p_val_adj, 1e-300))
  sig_colors <- c("Down" = "#2166AC", "NS" = "grey70", "Up" = "#B2182B")
  x_range <- range(deg$avg_log2FC, na.rm = TRUE)
  y_max <- max(deg$neglog10p, na.rm = TRUE) * 1.05

  p <- ggplot(deg, aes(x = avg_log2FC, y = neglog10p, color = sig)) +
    geom_point(size = 1.2, alpha = 0.7) +
    scale_color_manual(values = sig_colors, guide = "none") +
    geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    geom_text_repel(data = label_df,
      aes(x = avg_log2FC, y = -log10(pmax(p_val_adj, 1e-300)), label = gene),
      size = 3.2, max.overlaps = 25, segment.color = "grey50",
      fontface = "italic", show.legend = FALSE, min.segment.length = 0.2) +
    annotate("text", x = x_range[1] * 0.7, y = -y_max * 0.03,
             label = paste("UP in", label_left),
             hjust = 0.5, size = 4, fontface = "bold", color = "#2166AC") +
    annotate("text", x = x_range[2] * 0.7, y = -y_max * 0.03,
             label = paste("UP in", label_right),
             hjust = 0.5, size = 4, fontface = "bold", color = "#B2182B") +
    labs(title = title, x = expression(Log[2]*"(FC)"),
         y = expression(-Log[10]*"(adj. p value)")) +
    coord_cartesian(clip = "off") +
    theme_classic(base_size = 13) +
    theme(plot.title = element_text(face = "bold", size = 14),
          plot.margin = margin(t = 10, r = 10, b = 25, l = 10))

  ggsave(pdf_path, p, width = 8, height = 6)
  cat(sprintf("  Saved volcano: %s\n", pdf_path))
}

# Helper: DEG heatmap (Hu-style — shows top significant genes across conditions)
# condition_order: optional vector of condition names to control column ordering
make_deg_heatmap <- function(deg, seu_obj, group_col, groups, title, pdf_path,
                             exclude_genes = NULL, max_genes = 50,
                             condition_order = NULL) {
  deg$gene <- rownames(deg)
  if (!is.null(exclude_genes)) deg <- deg[!deg$gene %in% exclude_genes, ]

  sig_genes <- deg %>%
    filter(p_val_adj < 0.05, abs(avg_log2FC) > 0.25) %>%
    arrange(p_val_adj)

  if (nrow(sig_genes) == 0) {
    # Relax threshold
    sig_genes <- deg %>%
      filter(p_val_adj < 0.1) %>%
      arrange(p_val_adj)
    cat(sprintf("  Warning: relaxed threshold to padj<0.1 (%d genes)\n", nrow(sig_genes)))
  }

  if (nrow(sig_genes) == 0) {
    cat(sprintf("  SKIP heatmap: no significant genes for %s\n", title))
    return(invisible(NULL))
  }

  genes_use <- head(sig_genes$gene, max_genes)
  cat(sprintf("  Heatmap: %d genes for %s\n", length(genes_use), title))

  # Get scaled expression per sample
  mat <- GetAssayData(seu_obj, layer = "data")[genes_use, , drop = FALSE]
  samples <- unique(seu_obj$orig.ident)

  # Pseudobulk mean expression per sample
  pb_mat <- sapply(samples, function(s) {
    cells <- which(seu_obj$orig.ident == s)
    rowMeans(as.matrix(mat[, cells, drop = FALSE]))
  })
  if (!is.matrix(pb_mat)) pb_mat <- matrix(pb_mat, nrow = length(genes_use),
                                            dimnames = list(genes_use, samples))

  # Scale genes (z-score across samples) — skip for single gene
  if (nrow(pb_mat) == 1) {
    pb_scaled <- pb_mat
  } else {
    pb_scaled <- t(scale(t(pb_mat)))
    if (!is.matrix(pb_scaled)) pb_scaled <- matrix(pb_scaled, nrow = nrow(pb_mat),
                                                    dimnames = dimnames(pb_mat))
  }

  # Annotation
  sample_cond <- sapply(samples, function(s) {
    unique(seu_obj@meta.data[seu_obj$orig.ident == s, "condition"])[1]
  })
  anno_col <- data.frame(Condition = sample_cond, row.names = samples)
  present <- intersect(names(cond_colors), unique(sample_cond))
  anno_colors <- list(Condition = cond_colors[present])

  # Order columns: custom condition_order or alphabetical
  if (!is.null(condition_order)) {
    cond_rank <- setNames(seq_along(condition_order), condition_order)
    col_order <- order(cond_rank[sample_cond])
  } else {
    col_order <- order(sample_cond)
  }
  pb_scaled <- pb_scaled[, col_order, drop = FALSE]
  anno_col <- anno_col[col_order, , drop = FALSE]

  # Separate up/down genes
  up_genes <- sig_genes$gene[sig_genes$avg_log2FC > 0]
  down_genes <- sig_genes$gene[sig_genes$avg_log2FC < 0]
  row_order <- c(intersect(up_genes, genes_use), intersect(down_genes, genes_use))
  if (length(row_order) > 0) pb_scaled <- pb_scaled[row_order, , drop = FALSE]

  # Row annotation (direction)
  row_anno <- data.frame(
    Direction = ifelse(rownames(pb_scaled) %in% up_genes, "Up", "Down"),
    row.names = rownames(pb_scaled)
  )
  anno_colors$Direction <- c("Up" = "#B2182B", "Down" = "#2166AC")

  # Truncate extreme z-scores for visualization
  pb_scaled[pb_scaled > 2.5] <- 2.5
  pb_scaled[pb_scaled < -2.5] <- -2.5

  cairo_pdf(pdf_path, width = max(8, ncol(pb_scaled) * 0.6 + 4),
            height = max(6, nrow(pb_scaled) * 0.25 + 3))
  pheatmap(pb_scaled,
           cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_col = anno_col,
           annotation_row = row_anno,
           annotation_colors = anno_colors,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           fontsize_row = 7, fontsize_col = 8,
           main = title, border_color = NA,
           show_colnames = TRUE)
  dev.off()
  cat(sprintf("  Saved heatmap: %s\n", pdf_path))
}

# Helper: GO enrichment — butterfly bar chart combining UP + DOWN (manuscript style)
run_go_enrichment <- function(deg, title, prefix, out_dir,
                               label_left, label_right,
                               exclude_genes = NULL) {
  deg$gene <- rownames(deg)
  if (!is.null(exclude_genes)) deg <- deg[!deg$gene %in% exclude_genes, ]

  sig_up <- deg$gene[deg$p_val_adj < 0.05 & deg$avg_log2FC > 0]
  sig_down <- deg$gene[deg$p_val_adj < 0.05 & deg$avg_log2FC < 0]

  go_results <- list()

  for (direction in c("up", "down")) {
    genes <- if (direction == "up") sig_up else sig_down
    if (length(genes) < 3) {
      cat(sprintf("  GO %s: only %d genes, skipping\n", direction, length(genes)))
      next
    }

    ego <- tryCatch({
      enrichGO(gene     = genes,
               OrgDb    = org.Mm.eg.db,
               keyType  = "SYMBOL",
               ont      = "BP",
               pAdjustMethod = "BH",
               pvalueCutoff  = 0.05,
               qvalueCutoff  = 0.2,
               readable = TRUE)
    }, error = function(e) {
      cat(sprintf("  GO %s error: %s\n", direction, e$message))
      NULL
    })

    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      go_df <- as.data.frame(ego)
      go_df$direction <- direction
      go_results[[direction]] <- go_df
      write.csv(go_df, file.path(out_dir,
                sprintf("%s_GO_BP_%s.csv", prefix, direction)),
                row.names = FALSE)
      cat(sprintf("  GO %s: %d enriched terms\n", direction, nrow(go_df)))
    } else {
      cat(sprintf("  GO %s: no enriched terms\n", direction))
    }
  }

  # --- Combined butterfly bar chart (manuscript Fig 1E-H style) ---
  if (length(go_results) == 0) return(go_results)

  n_top <- 5
  frames <- list()

  # "down" = neg log2FC = UP in label_left → bars extend LEFT
  if ("down" %in% names(go_results)) {
    d <- head(go_results[["down"]], n_top)
    frames[["down"]] <- data.frame(
      Description = as.character(d$Description),
      gene_count  = as.numeric(sub("/.*", "", d$GeneRatio)),
      neg_log10p  = -log10(d$pvalue),
      side        = "left",
      stringsAsFactors = FALSE
    )
    frames[["down"]]$plot_x <- -frames[["down"]]$neg_log10p
  }

  # "up" = pos log2FC = UP in label_right → bars extend RIGHT
  if ("up" %in% names(go_results)) {
    u <- head(go_results[["up"]], n_top)
    frames[["up"]] <- data.frame(
      Description = as.character(u$Description),
      gene_count  = as.numeric(sub("/.*", "", u$GeneRatio)),
      neg_log10p  = -log10(u$pvalue),
      side        = "right",
      stringsAsFactors = FALSE
    )
    frames[["up"]]$plot_x <- frames[["up"]]$neg_log10p
  }

  plot_df <- do.call(rbind, frames)
  cat(sprintf("  Butterfly: %d rows (%s)\n", nrow(plot_df),
              paste(names(frames), sapply(frames, nrow), sep = "=", collapse = ", ")))
  plot_df$label <- stringr::str_wrap(plot_df$Description, width = 25)

  # Y ordering: left (down) terms bottom, right (up) terms top
  left_desc  <- if ("down" %in% names(frames))
    rev(frames[["down"]]$Description) else character(0)
  right_desc <- if ("up" %in% names(frames))
    frames[["up"]]$Description else character(0)
  plot_df$Description <- factor(plot_df$Description, levels = c(left_desc, right_desc))

  # Text near center axis
  plot_df$text_x    <- ifelse(plot_df$side == "left", -0.15, 0.15)
  plot_df$text_hjust <- ifelse(plot_df$side == "left", 1, 0)

  x_max <- max(abs(plot_df$plot_x)) * 1.15

  # Use distinct fill color per side so both are always visible
  side_colors <- c("left" = "#2166AC", "right" = "#B2182B")

  p <- ggplot(plot_df, aes(x = plot_x, y = Description, fill = side)) +
    geom_col(width = 0.85, alpha = 0.85) +
    geom_text(aes(x = text_x, label = label, hjust = text_hjust),
              size = 3, lineheight = 0.85, color = "black", show.legend = FALSE) +
    scale_fill_manual(values = side_colors, guide = "none") +
    scale_x_continuous(labels = abs, limits = c(-x_max, x_max)) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    labs(x = expression(-Log[10] ~ "(P-value)"), y = "", title = "GO BP") +
    annotate("text", x = -x_max * 0.5, y = nrow(plot_df) + 0.7,
             label = paste("UP in", label_left), hjust = 0.5,
             size = 3.5, fontface = "bold", color = "#2166AC") +
    annotate("text", x = x_max * 0.5, y = nrow(plot_df) + 0.7,
             label = paste("UP in", label_right), hjust = 0.5,
             size = 3.5, fontface = "bold", color = "#B2182B") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(t = 25, r = 10, b = 10, l = 10),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_cartesian(clip = "off")

  n_rows <- nrow(plot_df)
  save_pdf(p, file.path(out_dir, sprintf("%s_GO_BP_butterfly.pdf", prefix)),
           8, max(5, n_rows * 0.55 + 2))
  cat(sprintf("  Saved: %s_GO_BP_butterfly.pdf\n", prefix))

  return(go_results)
}

# Helper: GSEA using the full ranked gene list (all genes, not just significant)
run_gsea <- function(deg_full, title, prefix, out_dir,
                     label_left, label_right, exclude_genes = NULL) {
  deg_full$gene <- rownames(deg_full)
  if (!is.null(exclude_genes)) deg_full <- deg_full[!deg_full$gene %in% exclude_genes, ]

  # Rank by DESeq2 stat (= log2FC / lfcSE); positive = up in label_right
  gene_list <- setNames(deg_full$stat, deg_full$gene)
  gene_list <- sort(gene_list, decreasing = TRUE)
  gene_list <- gene_list[!is.na(gene_list)]

  gsea_res <- tryCatch({
    gseGO(geneList  = gene_list,
          OrgDb     = org.Mm.eg.db,
          keyType   = "SYMBOL",
          ont       = "BP",
          minGSSize = 10,
          maxGSSize = 500,
          pvalueCutoff = 0.05,
          pAdjustMethod = "BH",
          seed      = TRUE)
  }, error = function(e) {
    cat(sprintf("  GSEA error: %s\n", e$message))
    NULL
  })

  if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
    cat(sprintf("  GSEA: no enriched terms for %s\n", prefix))
    return(invisible(NULL))
  }

  gsea_df <- as.data.frame(gsea_res)
  write.csv(gsea_df, file.path(out_dir, sprintf("%s_GSEA_BP.csv", prefix)),
            row.names = FALSE)
  cat(sprintf("  GSEA: %d enriched BP terms for %s\n", nrow(gsea_df), prefix))

  # --- Butterfly bar chart: top activated (right) vs suppressed (left) ---
  n_top <- 5
  activated  <- gsea_df[gsea_df$NES > 0, ]
  suppressed <- gsea_df[gsea_df$NES < 0, ]

  frames <- list()
  if (nrow(suppressed) > 0) {
    s <- head(suppressed[order(suppressed$pvalue), ], n_top)
    frames[["left"]] <- data.frame(
      Description = as.character(s$Description),
      NES         = s$NES,
      neg_log10p  = -log10(s$pvalue),
      plot_x      = -(-log10(s$pvalue)),
      side        = "left",
      stringsAsFactors = FALSE
    )
  }
  if (nrow(activated) > 0) {
    a <- head(activated[order(activated$pvalue), ], n_top)
    frames[["right"]] <- data.frame(
      Description = as.character(a$Description),
      NES         = a$NES,
      neg_log10p  = -log10(a$pvalue),
      plot_x      = -log10(a$pvalue),
      side        = "right",
      stringsAsFactors = FALSE
    )
  }

  if (length(frames) == 0) return(invisible(gsea_df))

  plot_df <- do.call(rbind, frames)
  plot_df$label <- stringr::str_wrap(plot_df$Description, width = 25)

  left_desc  <- if ("left" %in% names(frames))
    rev(frames[["left"]]$Description) else character(0)
  right_desc <- if ("right" %in% names(frames))
    frames[["right"]]$Description else character(0)
  plot_df$Description <- factor(plot_df$Description, levels = c(left_desc, right_desc))

  plot_df$text_x    <- ifelse(plot_df$side == "left", -0.15, 0.15)
  plot_df$text_hjust <- ifelse(plot_df$side == "left", 1, 0)

  x_max <- max(abs(plot_df$plot_x)) * 1.15
  side_colors <- c("left" = "#2166AC", "right" = "#B2182B")

  p <- ggplot(plot_df, aes(x = plot_x, y = Description, fill = side)) +
    geom_col(width = 0.85, alpha = 0.85) +
    geom_text(aes(x = text_x, label = label, hjust = text_hjust),
              size = 3, lineheight = 0.85, color = "black", show.legend = FALSE) +
    scale_fill_manual(values = side_colors, guide = "none") +
    scale_x_continuous(labels = abs, limits = c(-x_max, x_max)) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.5) +
    labs(x = expression(-Log[10] ~ "(P-value)"), y = "", title = "GSEA GO BP") +
    annotate("text", x = -x_max * 0.5, y = nrow(plot_df) + 0.7,
             label = paste("UP in", label_left), hjust = 0.5,
             size = 3.5, fontface = "bold", color = "#2166AC") +
    annotate("text", x = x_max * 0.5, y = nrow(plot_df) + 0.7,
             label = paste("UP in", label_right), hjust = 0.5,
             size = 3.5, fontface = "bold", color = "#B2182B") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 12),
      panel.grid.major.y = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.8),
      plot.margin = margin(t = 25, r = 10, b = 10, l = 10),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_cartesian(clip = "off")

  n_rows <- nrow(plot_df)
  save_pdf(p, file.path(out_dir, sprintf("%s_GSEA_BP_butterfly.pdf", prefix)),
           8, max(5, n_rows * 0.55 + 2))
  cat(sprintf("  Saved: %s_GSEA_BP_butterfly.pdf\n", prefix))

  return(invisible(gsea_df))
}

# --- Define comparison panels (using pseudobulk DESeq2 results) ---
volcanos <- list(
  list(csv = "deg_pseudobulk_sex_effect.csv",
       prefix = "fig1E", title = "Eos DEG: WT Female vs WT Male",
       label_right = "WT Female", label_left = "WT Male",
       excl = c("Ddx3y"),
       groups = c("WT_Female_HFD", "WT_Male_HFD"),
       cond_order = NULL),
  list(csv = "deg_pseudobulk_genotype_female.csv",
       prefix = "fig1F", title = "Eos DEG: KO Female vs WT Female",
       label_right = "KO Female", label_left = "WT Female",
       excl = ko_exclude,
       groups = c("KO_Female_HFD", "WT_Female_HFD"),
       cond_order = NULL),
  list(csv = "deg_pseudobulk_genotype_male.csv",
       prefix = "fig1G", title = "Eos DEG: KO Male vs WT Male",
       label_right = "KO Male", label_left = "WT Male",
       excl = ko_exclude,
       groups = c("KO_Male_HFD", "WT_Male_HFD"),
       cond_order = NULL),
  list(csv = "deg_pseudobulk_infection_female.csv",
       prefix = "fig1H", title = "Eos DEG: Infected vs Naive Female",
       label_right = "Inf. Female", label_left = "Naive Female",
       excl = NULL,
       groups = c("WT_Female_HFD", "WT_Female_HFDInfected",
                  "KO_Female_HFD", "KO_Female_HFDInfected"),
       cond_order = c("WT_Female_HFD", "WT_Female_HFDInfected",
                      "KO_Female_HFD", "KO_Female_HFDInfected"))
)

for (v in volcanos) {
  csv_path <- file.path(deg_dir, v$csv)
  if (!file.exists(csv_path)) {
    cat(sprintf("  SKIP: %s not found\n", v$csv))
    next
  }
  deg <- read.csv(csv_path)
  # Standardize DESeq2 column names to match code expectations
  deg <- deg %>%
    filter(!is.na(padj)) %>%
    rename(avg_log2FC = log2FC, p_val = pval, p_val_adj = padj)
  rownames(deg) <- deg$gene
  if (!is.null(v$excl)) {
    deg <- deg[!deg$gene %in% v$excl, ]
    rownames(deg) <- deg$gene
    cat(sprintf("  Excluded genes: %s\n", paste(v$excl, collapse = ", ")))
  }
  n_sig <- sum(deg$p_val_adj < 0.05 & abs(deg$avg_log2FC) > 0.25, na.rm = TRUE)
  cat(sprintf("\n%s: %d genes, %d significant (padj<0.05, |FC|>0.25)\n",
      v$csv, nrow(deg), n_sig))

  # Volcano
  make_volcano(deg, v$title,
               file.path(out_dir, paste0(v$prefix, "_volcano.pdf")),
               label_left = v$label_left, label_right = v$label_right,
               exclude_genes = v$excl)

  # DEG heatmap (Hu-style)
  tier1_sub <- subset(tier1, condition %in% v$groups)
  make_deg_heatmap(deg, tier1_sub, "condition", v$groups,
                   v$title,
                   file.path(out_dir, paste0(v$prefix, "_heatmap.pdf")),
                   exclude_genes = v$excl,
                   condition_order = v$cond_order)

  # GO enrichment (ORA butterfly bar chart)
  run_go_enrichment(deg, v$title, v$prefix, out_dir,
                    label_left = v$label_left, label_right = v$label_right,
                    exclude_genes = v$excl)

  # GSEA (full ranked list)
  run_gsea(deg, v$title, v$prefix, out_dir,
           label_left = v$label_left, label_right = v$label_right,
           exclude_genes = v$excl)
}

rm(obj, tier1); gc()
cat(sprintf("\nDone! All Fig 1 panels saved to %s\n", out_dir))
