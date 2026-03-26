#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

## Cell type annotation via label transfer from previous round (R02)
## SLURM: --mem=128G --cpus-per-task=4 --time=2:00:00 -p batch

# ── Libraries ──────────────────────────────────────────────────────────────────
suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(patchwork)
  library(dplyr)
})

# ── Shared theme ─────────────────────────────────────────────────────────────
SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

# ── Configuration ──────────────────────────────────────────────────────────────
BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# Previous labeled object (13 cell types, 97,061 cells)
PREV_LABELED <- "/bigdata/nairlab/shared/XQiu/EOS_project/R02_draft_v1/data/processed/seurat_object_labeled.rds"

# Broad category mapping
broad_map <- c(
  Eosinophils              = "Granulocyte",
  Macrophages              = "Myeloid",
  "Endothelial cells"      = "Stromal",
  "Dendritic cells"        = "Myeloid",
  "Mesenchymal stem cells" = "Stromal",
  Pericytes                = "Stromal",
  "B cells"                = "Lymphoid",
  "NK cells"               = "Lymphoid",
  "Mesothelial cells"      = "Stromal",
  "Mast cells"             = "Granulocyte",
  Adipocytes               = "Stromal",
  Neutrophils              = "Granulocyte",
  "Epithelial cells"       = "Stromal"
)

# ── Load data ────────────────────────────────────────────────────────────────
cat("========================================\n")
cat("Cell type annotation via label transfer\n")
cat("========================================\n\n")

input_file <- file.path(BASE, "all_cells_clustered.rds")
cat("Loading:", input_file, "\n")
sobj <- readRDS(input_file)
cat(sprintf("Loaded: %d genes x %d cells\n", nrow(sobj), ncol(sobj)))
cat(sprintf("Clusters: %d\n\n", length(unique(Idents(sobj)))))

# ── Step 1: FindAllMarkers (for reference, saved to CSV) ─────────────────────
cat("========================================\n")
cat("Step 1: FindAllMarkers\n")
cat("========================================\n\n")

cat("Running FindAllMarkers (only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25)...\n")
markers <- tryCatch({
  FindAllMarkers(sobj, only.pos = TRUE, min.pct = 0.25,
                 logfc.threshold = 0.25, verbose = FALSE)
}, error = function(e) {
  cat(sprintf("ERROR in FindAllMarkers: %s\n", e$message))
  NULL
})

if (!is.null(markers)) {
  cat(sprintf("Found %d marker genes across %d clusters\n",
              nrow(markers), length(unique(markers$cluster))))

  top5 <- markers %>%
    group_by(cluster) %>%
    top_n(5, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))
  cat("\nTop 5 markers per cluster:\n")
  for (cl in sort(unique(top5$cluster))) {
    genes <- top5$gene[top5$cluster == cl]
    cat(sprintf("  Cluster %s: %s\n", cl, paste(genes, collapse = ", ")))
  }
}

# ── Step 2: Label transfer from previous round ───────────────────────────────
cat("\n========================================\n")
cat("Step 2: Label transfer from R02 labeled object\n")
cat("========================================\n\n")

if (!file.exists(PREV_LABELED)) {
  stop("Previous labeled object not found: ", PREV_LABELED)
}

cat("Loading:", PREV_LABELED, "\n")
prev_obj <- readRDS(PREV_LABELED)
cat(sprintf("  Previous object: %d cells, %d cell types\n",
            ncol(prev_obj), length(unique(na.omit(as.character(prev_obj$cell_type))))))
cat("  Cell types: ", paste(sort(unique(na.omit(as.character(prev_obj$cell_type)))), collapse = ", "), "\n")

# Build barcode -> cell_type map (as.character to avoid factor levels)
prev_labels <- setNames(as.character(prev_obj$cell_type), colnames(prev_obj))
rm(prev_obj); gc(verbose = FALSE)

# Match barcodes
shared_barcodes <- intersect(colnames(sobj), names(prev_labels))
unmatched_barcodes <- setdiff(colnames(sobj), names(prev_labels))
cat(sprintf("\n  Shared barcodes: %d / %d (%.1f%%)\n",
            length(shared_barcodes), ncol(sobj),
            100 * length(shared_barcodes) / ncol(sobj)))
cat(sprintf("  Unmatched barcodes: %d\n", length(unmatched_barcodes)))

# Transfer labels
sobj@meta.data$celltype <- NA_character_
sobj@meta.data[shared_barcodes, "celltype"] <- prev_labels[shared_barcodes]

cat("\nTransferred cell type labels:\n")
print(table(sobj$celltype, useNA = "ifany"))

# For unmatched cells, label as "Unassigned"
n_na <- sum(is.na(sobj$celltype))
if (n_na > 0) {
  cat(sprintf("\n%d cells without previous labels -> labeled 'Unassigned'\n", n_na))
  sobj@meta.data$celltype[is.na(sobj@meta.data$celltype)] <- "Unassigned"
}

# Assign broad cell type
sobj@meta.data$broad_celltype <- broad_map[sobj$celltype]
sobj@meta.data$broad_celltype[is.na(sobj$broad_celltype)] <- "Other"

# ── Step 3: Cross-tabulate clusters vs cell types ────────────────────────────
cat("\n========================================\n")
cat("Step 3: Cluster vs cell type mapping\n")
cat("========================================\n\n")

# Majority cell type per cluster
cluster_ct <- table(Cluster = as.character(Idents(sobj)), CellType = sobj$celltype)
cat("Majority cell type per cluster:\n")
for (cl in rownames(cluster_ct)) {
  majority <- names(which.max(cluster_ct[cl, ]))
  pct <- 100 * max(cluster_ct[cl, ]) / sum(cluster_ct[cl, ])
  cat(sprintf("  Cluster %s: %s (%.1f%%, n=%d)\n", cl, majority, pct, sum(cluster_ct[cl, ])))
}

# ── Summary ──────────────────────────────────────────────────────────────────
cat("\n========================================\n")
cat("Cell type summary\n")
cat("========================================\n\n")

cat("Cells per cell type:\n")
print(sort(table(sobj$celltype), decreasing = TRUE))

cat("\nCells per broad cell type:\n")
print(sort(table(sobj$broad_celltype), decreasing = TRUE))

cat("\nCell type proportions per condition:\n")
prop_table <- prop.table(table(sobj$condition, sobj$celltype), margin = 1)
print(round(prop_table, 3))

# ── Step 4: Figures (Fig 1C style — immune cells only) ────────────────────────
cat("\n========================================\n")
cat("Step 4: Generating figures\n")
cat("========================================\n\n")

# 7 immune cell types shown in manuscript Fig 1C
immune_types <- c("Eosinophils", "Macrophages", "Dendritic cells",
                  "B cells", "NK cells", "Neutrophils", "Mast cells")

# Subset to immune cells only for UMAP (matching manuscript Fig 1C)
immune_cells <- colnames(sobj)[sobj$celltype %in% immune_types]
sobj_immune <- subset(sobj, cells = immune_cells)
cat(sprintf("Immune cells for Fig 1C UMAP: %d / %d\n",
            ncol(sobj_immune), ncol(sobj)))

# Color mapping (from broad_celltype_colors in 00_theme.R)
ct_colors <- broad_celltype_colors[immune_types]

# UMAP colored by cell type (immune only, Fig 1C style)
cat("UMAP colored by celltype (immune only, Fig 1C)...\n")
sobj_immune$celltype <- factor(sobj_immune$celltype, levels = immune_types)
p1 <- DimPlot(sobj_immune, group.by = "celltype", label = TRUE, repel = TRUE,
              pt.size = 0.3, cols = ct_colors) +
  ggtitle("Immune Cells") + theme_umap()
save_pdf(p1, file.path(fig_dir, "umap_celltype.pdf"), width = 10, height = 8)

# UMAP with ALL cell types (for reference, non-immune in grey)
cat("UMAP colored by celltype (all cells, non-immune grey)...\n")
all_ct <- sort(unique(sobj$celltype))
all_colors <- setNames(rep("grey85", length(all_ct)), all_ct)
all_colors[names(ct_colors)] <- ct_colors
p1_all <- DimPlot(sobj, group.by = "celltype", label = TRUE, repel = TRUE,
                  pt.size = 0.2, cols = all_colors) +
  ggtitle("All Cells (non-immune in grey)") + theme_umap()
save_pdf(p1_all, file.path(fig_dir, "umap_celltype_all.pdf"), width = 10, height = 8)

# Stacked bar plot by condition (immune cells only)
cat("Stacked bar by condition (immune only)...\n")
ct_counts <- as.data.frame(table(sobj_immune$condition, sobj_immune$celltype))
colnames(ct_counts) <- c("condition", "celltype", "count")

p2 <- ggplot(ct_counts, aes(x = condition, y = count, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = ct_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Cell Type Composition by Condition",
       x = "Condition", y = "Proportion", fill = "Cell Type") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))
save_pdf(p2, file.path(fig_dir, "barplot_celltype_by_condition.pdf"), width = 10, height = 7)

# Combined panel
p_combined <- p1 | p2
save_pdf(p_combined, file.path(fig_dir, "celltype_annotation_panel.pdf"), width = 18, height = 8)

# ── Save ─────────────────────────────────────────────────────────────────────
cat("\n========================================\n")
cat("Saving results\n")
cat("========================================\n\n")

out_rds <- file.path(BASE, "all_cells_annotated.rds")
cat("Saving:", out_rds, "\n")
saveRDS(sobj, out_rds)

if (!is.null(markers)) {
  markers_file <- file.path(BASE, "celltype_markers.csv")
  cat("Saving:", markers_file, "\n")
  write.csv(markers, markers_file, row.names = FALSE)
}

cat(sprintf("\nFigures saved to: %s\n", fig_dir))
cat("Done.\n")
