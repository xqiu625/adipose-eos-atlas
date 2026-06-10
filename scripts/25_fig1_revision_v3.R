#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 25_fig1_revision_v3.R
# Revised Figure 1 panels addressing Rebecca's 2026-04-03 email:
#   1. Fig 1C: Updated UMAP with C6 (CEos) removed
#   2. Fig 1D: Updated dotplot with C6 removed
#   3. GSEA for all DEG comparisons (fig1E-H + multifactor + infection splits)
#
# Based on 10_fig1_revision_v2.R — only regenerates panels that need updating.
# Previously completed panels (fig1B PCA/MDS, volcanos, heatmaps, ORA GO)
# are NOT re-run here; see 10_fig1_revision_v2.R for those.
#
# Input:  all_cells_annotated.rds, eos_subclustered.rds, DESeq2 CSVs
# Output: PDF figures + GSEA CSVs in revision/figures/fig1/
# SLURM:  batch, 128G, 4 CPUs, 4h
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(stringr)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

# --- Paths ---
PNAS    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision"
BASE    <- file.path(PNAS, "atlas")
deg_dir <- file.path(BASE, "deg")
out_dir <- file.path(PNAS, "revision/figures/fig1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)


# =============================================================================
# Fig 1C: Updated UMAPs — C6 (CEos) removed
# Rebecca: "could you provide the following new figures (the umap ... since
#           now we omitted CEos)"
# =============================================================================
cat("=== Fig 1C: UMAPs (C6 removed) ===\n")

# Identify C6 barcodes: eosinophils in all_cells_annotated that are NOT in
# eos_subclustered (which had C6 removed by 04b_rename_clusters.R)
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
eos_barcodes <- colnames(eos)  # C0-C5 only (C6 already removed by 04b)
rm(eos); gc()

# Load all-cells object first to find the C6 cells
obj <- readRDS(file.path(BASE, "all_cells_annotated.rds"))
if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

all_eos_barcodes <- colnames(obj)[obj$celltype == "Eosinophils"]
c6_barcodes <- setdiff(all_eos_barcodes, eos_barcodes)
cat(sprintf("  C6 barcodes to remove: %d (eos in all_cells but not in eos_subclustered)\n",
            length(c6_barcodes)))

cat(sprintf("  All cells before: %d\n", ncol(obj)))

cells_keep <- setdiff(colnames(obj), c6_barcodes)
n_removed <- ncol(obj) - length(cells_keep)
obj <- subset(obj, cells = cells_keep)
cat(sprintf("  After C6 removal: %d cells (%d removed)\n", ncol(obj), n_removed))

immune_types <- c("Eosinophils", "Macrophages", "Dendritic cells",
                  "B cells", "NK cells", "Neutrophils", "Mast cells")
ct_colors <- broad_celltype_colors[immune_types]

# (1) Immune cells only
sobj_immune <- subset(obj, celltype %in% immune_types)
sobj_immune$celltype <- factor(sobj_immune$celltype, levels = immune_types)
cat(sprintf("  Immune cells: %d\n", ncol(sobj_immune)))

p_immune <- DimPlot(sobj_immune, group.by = "celltype", label = TRUE,
                    repel = TRUE, pt.size = 0.3, cols = ct_colors) +
  ggtitle("Immune Cells") + theme_umap()
save_pdf(p_immune, file.path(out_dir, "fig1C_umap_immune_only.pdf"), 10, 8)

# (2) All cells: immune colored, others in gray
all_ct <- sort(unique(obj$celltype))
all_colors <- setNames(rep("grey85", length(all_ct)), all_ct)
all_colors[names(ct_colors)] <- ct_colors

p_all <- DimPlot(obj, group.by = "celltype", label = TRUE, repel = TRUE,
                 pt.size = 0.2, cols = all_colors,
                 order = rev(immune_types)) +
  ggtitle("All Cells (non-immune in gray)") + theme_umap()
save_pdf(p_all, file.path(out_dir, "fig1C_umap_all_cells_immune_highlight.pdf"), 10, 8)


# =============================================================================
# Fig 1D: Dot plot of canonical immune cell markers — C6 removed
# Rebecca: "and also the dot plot of immune cell marker genes, since now
#           we omitted CEos"
# =============================================================================
cat("\n=== Fig 1D: Immune cell marker dot plot (C6 removed) ===\n")

sobj_immune$celltype <- factor(sobj_immune$celltype,
  levels = rev(c("Mast cells", "Dendritic cells", "B cells",
                 "NK cells", "Neutrophils", "Macrophages", "Eosinophils")))

marker_genes <- c("Siglecf", "Ccr3", "Il5ra",             # Eosinophils
                   "Mertk", "Mrc1", "Fcgr1", "Retnla",    # Macrophages
                   "S100a8",                                # Neutrophils
                   "Nkg7",                                  # NK cells
                   "Cd19",                                  # B cells
                   "Itgax",                                 # Dendritic cells
                   "Tpsb2")                                 # Mast cells
marker_genes <- marker_genes[marker_genes %in% rownames(sobj_immune)]
cat(sprintf("  Using %d marker genes\n", length(marker_genes)))

p_dot <- DotPlot(sobj_immune, features = marker_genes, group.by = "celltype",
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

rm(obj, sobj_immune); gc()


# =============================================================================
# GSEA helper — butterfly bar chart (same style as 10_fig1_revision_v2.R)
# =============================================================================
run_gsea <- function(deg_csv, prefix, title, label_left, label_right,
                     out_dir, exclude_genes = NULL) {
  deg <- read.csv(deg_csv)
  if (!is.null(exclude_genes)) deg <- deg[!deg$gene %in% exclude_genes, ]

  # Use DESeq2 stat column for ranking (= log2FC / lfcSE)
  if ("stat" %in% colnames(deg)) {
    gene_list <- setNames(deg$stat, deg$gene)
  } else {
    cat(sprintf("  WARNING: no 'stat' column in %s, using log2FC\n", basename(deg_csv)))
    gene_list <- setNames(deg$log2FC, deg$gene)
  }
  gene_list <- sort(gene_list[!is.na(gene_list)], decreasing = TRUE)
  cat(sprintf("  GSEA ranked list: %d genes from %s\n", length(gene_list), basename(deg_csv)))

  # Single fixed GSEA: rank by DESeq2 stat, BH-adjusted, report terms at FDR < 0.05
  gsea_res <- tryCatch({
    gseGO(geneList      = gene_list,
          OrgDb         = org.Mm.eg.db,
          keyType       = "SYMBOL",
          ont           = "BP",
          minGSSize     = 10,
          maxGSSize     = 500,
          pvalueCutoff  = 0.05,
          pAdjustMethod = "BH",
          seed          = TRUE)
  }, error = function(e) {
    cat(sprintf("  GSEA error: %s\n", e$message))
    NULL
  })

  if (is.null(gsea_res) || nrow(as.data.frame(gsea_res)) == 0) {
    cat(sprintf("  GSEA: no terms at FDR < 0.05 for %s\n", prefix))
    write.csv(data.frame(note = "No significant GSEA terms (FDR < 0.05)"),
              file.path(out_dir, sprintf("%s_GSEA_BP.csv", prefix)),
              row.names = FALSE)
    return(invisible(NULL))
  }

  gsea_df <- as.data.frame(gsea_res)
  write.csv(gsea_df, file.path(out_dir, sprintf("%s_GSEA_BP.csv", prefix)),
            row.names = FALSE)
  cat(sprintf("  GSEA: %d enriched BP terms for %s\n", nrow(gsea_df), prefix))

  # --- Butterfly bar chart ---
  n_top <- 5
  activated  <- gsea_df[gsea_df$NES > 0, ]
  suppressed <- gsea_df[gsea_df$NES < 0, ]

  frames <- list()
  if (nrow(suppressed) > 0) {
    s <- head(suppressed[order(suppressed$p.adjust), ], n_top)
    frames[["left"]] <- data.frame(
      Description = as.character(s$Description),
      NES = s$NES, neg_log10p = -log10(s$p.adjust),
      plot_x = -(-log10(s$p.adjust)), side = "left",
      stringsAsFactors = FALSE)
  }
  if (nrow(activated) > 0) {
    a <- head(activated[order(activated$p.adjust), ], n_top)
    frames[["right"]] <- data.frame(
      Description = as.character(a$Description),
      NES = a$NES, neg_log10p = -log10(a$p.adjust),
      plot_x = -log10(a$p.adjust), side = "right",
      stringsAsFactors = FALSE)
  }
  if (length(frames) == 0) return(invisible(gsea_df))

  plot_df <- do.call(rbind, frames)
  plot_df$label <- str_wrap(plot_df$Description, width = 25)

  left_desc  <- if ("left" %in% names(frames)) rev(frames[["left"]]$Description) else character(0)
  right_desc <- if ("right" %in% names(frames)) frames[["right"]]$Description else character(0)
  plot_df$Description <- factor(plot_df$Description, levels = c(left_desc, right_desc))
  plot_df$text_x    <- ifelse(plot_df$side == "left", -0.15, 0.15)
  plot_df$text_hjust <- ifelse(plot_df$side == "left", 1, 0)
  x_max <- max(abs(plot_df$plot_x)) * 1.15
  side_colors <- c("left" = "#2166AC", "right" = "#B2182B")

  # 5.5cm x 7.5cm = 2.17in x 2.95in
  w_in <- 5.5 / 2.54
  h_in <- 7.5 / 2.54

  p <- ggplot(plot_df, aes(x = plot_x, y = Description, fill = side)) +
    geom_col(width = 0.85, alpha = 0.85) +
    geom_text(aes(x = text_x, label = label, hjust = text_hjust),
              size = 1.8, lineheight = 0.8, color = "black", show.legend = FALSE) +
    scale_fill_manual(values = side_colors, guide = "none") +
    scale_x_continuous(labels = abs, limits = c(-x_max, x_max)) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
    labs(x = expression(-Log[10] ~ "(adjusted P)"), y = "", title = "GSEA GO BP") +
    annotate("text", x = -x_max * 0.5, y = nrow(plot_df) + 0.7,
             label = paste("UP in", label_left), hjust = 0.5,
             size = 2, fontface = "bold", color = "#2166AC") +
    annotate("text", x = x_max * 0.5, y = nrow(plot_df) + 0.7,
             label = paste("UP in", label_right), hjust = 0.5,
             size = 2, fontface = "bold", color = "#B2182B") +
    theme_minimal(base_size = 7) +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 6),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 7),
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(t = 15, r = 5, b = 5, l = 5),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_cartesian(clip = "off")

  save_pdf(p, file.path(out_dir, sprintf("%s_GSEA_BP_butterfly.pdf", prefix)),
           w_in, h_in)

  return(invisible(gsea_df))
}


# =============================================================================
# Re-generate ORA GO butterfly plots at 5.5cm x 7.5cm (Rebecca's size request)
# Reads existing GO CSV results and re-plots at new size
# =============================================================================
cat("\n=== Re-generating GO butterfly plots at 5.5cm x 7.5cm ===\n")

go_w <- 5.5 / 2.54  # inches
go_h <- 7.5 / 2.54

replot_go_butterfly <- function(prefix, label_left, label_right, out_dir) {
  up_csv <- file.path(out_dir, sprintf("%s_GO_BP_up.csv", prefix))
  dn_csv <- file.path(out_dir, sprintf("%s_GO_BP_down.csv", prefix))

  frames <- list()
  n_top <- 5

  if (file.exists(dn_csv)) {
    d <- read.csv(dn_csv)
    if (nrow(d) > 0) {
      d <- head(d[order(d$pvalue), ], n_top)
      frames[["left"]] <- data.frame(
        Description = as.character(d$Description),
        plot_x = -(-log10(d$pvalue)), side = "left",
        stringsAsFactors = FALSE)
    }
  }
  if (file.exists(up_csv)) {
    u <- read.csv(up_csv)
    if (nrow(u) > 0) {
      u <- head(u[order(u$pvalue), ], n_top)
      frames[["right"]] <- data.frame(
        Description = as.character(u$Description),
        plot_x = -log10(u$pvalue), side = "right",
        stringsAsFactors = FALSE)
    }
  }

  if (length(frames) == 0) { cat(sprintf("  SKIP %s: no GO results\n", prefix)); return() }

  plot_df <- do.call(rbind, frames)
  plot_df$label <- str_wrap(plot_df$Description, width = 25)

  left_desc  <- if ("left" %in% names(frames)) rev(frames[["left"]]$Description) else character(0)
  right_desc <- if ("right" %in% names(frames)) frames[["right"]]$Description else character(0)
  plot_df$Description <- factor(plot_df$Description, levels = c(left_desc, right_desc))
  plot_df$text_x    <- ifelse(plot_df$side == "left", -0.15, 0.15)
  plot_df$text_hjust <- ifelse(plot_df$side == "left", 1, 0)
  x_max <- max(abs(plot_df$plot_x)) * 1.15
  side_colors <- c("left" = "#2166AC", "right" = "#B2182B")

  p <- ggplot(plot_df, aes(x = plot_x, y = Description, fill = side)) +
    geom_col(width = 0.85, alpha = 0.85) +
    geom_text(aes(x = text_x, label = label, hjust = text_hjust),
              size = 1.8, lineheight = 0.8, color = "black", show.legend = FALSE) +
    scale_fill_manual(values = side_colors, guide = "none") +
    scale_x_continuous(labels = abs, limits = c(-x_max, x_max)) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
    labs(x = expression(-Log[10] ~ "(P-value)"), y = "", title = "GO BP") +
    annotate("text", x = -x_max * 0.5, y = nrow(plot_df) + 0.7,
             label = paste("UP in", label_left), hjust = 0.5,
             size = 2, fontface = "bold", color = "#2166AC") +
    annotate("text", x = x_max * 0.5, y = nrow(plot_df) + 0.7,
             label = paste("UP in", label_right), hjust = 0.5,
             size = 2, fontface = "bold", color = "#B2182B") +
    theme_minimal(base_size = 7) +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 6),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 7),
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(t = 15, r = 5, b = 5, l = 5),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_cartesian(clip = "off")

  save_pdf(p, file.path(out_dir, sprintf("%s_GO_BP_butterfly.pdf", prefix)), go_w, go_h)
}

go_panels <- list(
  list(prefix = "fig1E", left = "WT Male", right = "WT Female"),
  list(prefix = "fig1F", left = "WT Female", right = "KO Female"),
  list(prefix = "fig1G", left = "WT Male", right = "KO Male"),
  list(prefix = "fig1H", left = "Naive Female", right = "Inf. Female")
)

for (g in go_panels) {
  cat(sprintf("  Replotting %s GO butterfly...\n", g$prefix))
  replot_go_butterfly(g$prefix, g$left, g$right, out_dir)
}


# =============================================================================
# Fig 1E-H GSEA: re-run on all 4 pseudobulk DESeq2 comparisons
# Rebecca: "And also maybe GSEA analysis as well?"
# =============================================================================
cat("\n=== GSEA: Fig 1E-H (pseudobulk DESeq2 comparisons) ===\n")

ko_exclude <- c("Retnla", "Igha")

gsea_jobs <- list(
  list(csv = file.path(deg_dir, "deg_pseudobulk_sex_effect.csv"),
       prefix = "fig1E", title = "Sex: WT Female vs WT Male",
       left = "WT Male", right = "WT Female",
       excl = c("Ddx3y")),
  list(csv = file.path(deg_dir, "deg_pseudobulk_genotype_female.csv"),
       prefix = "fig1F", title = "Genotype: KO Female vs WT Female",
       left = "WT Female", right = "KO Female",
       excl = ko_exclude),
  list(csv = file.path(deg_dir, "deg_pseudobulk_genotype_male.csv"),
       prefix = "fig1G", title = "Genotype: KO Male vs WT Male",
       left = "WT Male", right = "KO Male",
       excl = ko_exclude),
  list(csv = file.path(deg_dir, "deg_pseudobulk_infection_female.csv"),
       prefix = "fig1H", title = "Infection: Infected vs Naive (all females)",
       left = "Naive Female", right = "Inf. Female",
       excl = NULL)
)

for (job in gsea_jobs) {
  cat(sprintf("\n--- GSEA: %s ---\n", job$title))
  if (file.exists(job$csv)) {
    run_gsea(job$csv, job$prefix, job$title,
             job$left, job$right, out_dir, job$excl)
  } else {
    cat(sprintf("  SKIP: %s not found\n", basename(job$csv)))
  }
}


# =============================================================================
# Multifactor GSEA: genotype + infection as covariates (all 10 female samples)
# Rebecca agreed with this approach in her 2026-04-03 email
# =============================================================================
cat("\n=== GSEA: Multifactor (all 10 females, ~ genotype + infection) ===\n")

mf_jobs <- list(
  list(csv = file.path(deg_dir, "deg_pseudobulk_multifactor_genotype.csv"),
       prefix = "multifactor_genotype",
       title = "Genotype (KO vs WT, ctrl for infection)",
       left = "WT", right = "KO"),
  list(csv = file.path(deg_dir, "deg_pseudobulk_multifactor_infection.csv"),
       prefix = "multifactor_infection",
       title = "Infection (Inf vs Naive, ctrl for genotype)",
       left = "Naive", right = "Infected")
)

for (job in mf_jobs) {
  cat(sprintf("\n--- GSEA: %s ---\n", job$title))
  if (file.exists(job$csv)) {
    run_gsea(job$csv, job$prefix, job$title,
             job$left, job$right, out_dir)
  } else {
    cat(sprintf("  SKIP: %s not found\n", basename(job$csv)))
  }
}


cat(sprintf("\nDone! All Fig 1 updates saved to %s\n", out_dir))
