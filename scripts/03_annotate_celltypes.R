#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

## Cell type annotation using marker genes
## SLURM: --mem=64G --cpus-per-task=4 --time=2:00:00 -p batch

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

# Cell type marker gene sets
marker_sets <- list(
  Eosinophil    = c("Epx", "Prg2", "Prg3", "Ear1", "Ear2", "Ear6",
                    "Ccr3", "Siglecf", "Il5ra"),
  Macrophage    = c("Adgre1", "Csf1r", "Cd68", "Lyz2", "C1qa", "C1qb"),
  T_cell        = c("Cd3d", "Cd3e", "Cd3g", "Trac"),
  B_cell        = c("Cd79a", "Cd79b", "Ms4a1"),
  NK_cell       = c("Nkg7", "Klrb1c", "Ncr1"),
  Dendritic     = c("Itgax", "H2-Ab1", "Flt3"),
  Neutrophil    = c("S100a8", "S100a9", "Ly6g"),
  Other_immune  = c("Ptprc")
)

# Broad category mapping
broad_map <- c(
  Eosinophil   = "Granulocyte",
  Macrophage   = "Myeloid",
  T_cell       = "Lymphoid",
  B_cell       = "Lymphoid",
  NK_cell      = "Lymphoid",
  Dendritic    = "Myeloid",
  Neutrophil   = "Granulocyte",
  Other_immune = "Other"
)

# ── Load data ──────────────────────────────────────────────────────────────────
cat("========================================\n")
cat("Cell type annotation using marker genes\n")
cat("========================================\n\n")

input_file <- file.path(BASE, "all_cells_clustered.rds")
cat("Loading:", input_file, "\n")
sobj <- readRDS(input_file)
cat(sprintf("Loaded: %d genes x %d cells\n", nrow(sobj), ncol(sobj)))
cat(sprintf("Clusters: %d\n\n", length(unique(Idents(sobj)))))

# ── Step 1: FindAllMarkers ─────────────────────────────────────────────────────
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

  # Top 5 markers per cluster
  cat("\nTop 5 markers per cluster:\n")
  top5 <- markers %>%
    group_by(cluster) %>%
    top_n(5, wt = avg_log2FC) %>%
    arrange(cluster, desc(avg_log2FC))
  for (cl in sort(unique(top5$cluster))) {
    genes <- top5$gene[top5$cluster == cl]
    cat(sprintf("  Cluster %s: %s\n", cl, paste(genes, collapse = ", ")))
  }
}

# ── Step 2: AddModuleScore for each cell type ─────────────────────────────────
cat("\n========================================\n")
cat("Step 2: Cell type module scoring\n")
cat("========================================\n\n")

# Check which markers are present in the dataset
all_genes <- rownames(sobj)
for (ct in names(marker_sets)) {
  present <- marker_sets[[ct]][marker_sets[[ct]] %in% all_genes]
  missing <- marker_sets[[ct]][!marker_sets[[ct]] %in% all_genes]
  cat(sprintf("  %s: %d/%d markers present", ct,
              length(present), length(marker_sets[[ct]])))
  if (length(missing) > 0) {
    cat(sprintf(" (missing: %s)", paste(missing, collapse = ", ")))
  }
  cat("\n")
}
cat("\n")

# Compute module scores
for (ct in names(marker_sets)) {
  genes_present <- marker_sets[[ct]][marker_sets[[ct]] %in% all_genes]
  if (length(genes_present) < 1) {
    cat(sprintf("Skipping %s: no markers found in dataset\n", ct))
    next
  }
  cat(sprintf("Scoring %s (%d genes)...\n", ct, length(genes_present)))
  score_name <- paste0(ct, "_score")
  tryCatch({
    sobj <- AddModuleScore(sobj, features = list(genes_present),
                           name = score_name, verbose = FALSE)
    # AddModuleScore appends "1" to the name
    colnames(sobj@meta.data)[colnames(sobj@meta.data) == paste0(score_name, "1")] <- score_name
  }, error = function(e) {
    cat(sprintf("  ERROR scoring %s: %s\n", ct, e$message))
  })
}

# ── Step 3: Assign clusters to cell types ──────────────────────────────────────
cat("\n========================================\n")
cat("Step 3: Assign clusters to cell types\n")
cat("========================================\n\n")

score_cols <- paste0(names(marker_sets), "_score")
score_cols <- score_cols[score_cols %in% colnames(sobj@meta.data)]

# Compute mean score per cluster for each cell type
clusters <- sort(unique(Idents(sobj)))
score_matrix <- matrix(NA, nrow = length(clusters), ncol = length(score_cols),
                       dimnames = list(as.character(clusters), score_cols))

for (cl in clusters) {
  cells <- WhichCells(sobj, idents = cl)
  for (sc in score_cols) {
    score_matrix[as.character(cl), sc] <- mean(sobj@meta.data[cells, sc],
                                                na.rm = TRUE)
  }
}

cat("Mean module scores per cluster:\n")
print(round(score_matrix, 4))
cat("\n")

# Assign cell type to each cluster based on highest mean score
cluster_assignments <- data.frame(
  cluster  = as.character(clusters),
  celltype = NA_character_,
  max_score = NA_real_,
  stringsAsFactors = FALSE
)

for (i in seq_along(clusters)) {
  cl <- as.character(clusters[i])
  scores <- score_matrix[cl, ]
  best_idx <- which.max(scores)
  best_type <- sub("_score$", "", names(scores)[best_idx])
  cluster_assignments$celltype[i]  <- best_type
  cluster_assignments$max_score[i] <- scores[best_idx]
}

cat("Cluster -> Cell type assignments:\n")
for (i in seq_len(nrow(cluster_assignments))) {
  cat(sprintf("  Cluster %s -> %s (score: %.4f)\n",
              cluster_assignments$cluster[i],
              cluster_assignments$celltype[i],
              cluster_assignments$max_score[i]))
}

# Add celltype and broad_celltype to Seurat object
cluster_to_type <- setNames(cluster_assignments$celltype,
                            cluster_assignments$cluster)
sobj$celltype <- cluster_to_type[as.character(Idents(sobj))]
sobj$broad_celltype <- broad_map[sobj$celltype]

# ── Step 3b: Transfer previous round annotations by barcode matching ──────────
cat("\n========================================\n")
cat("Step 3b: Transfer previous round annotations\n")
cat("========================================\n\n")

# Use the fully labeled object from R02 draft (13 cell types, 97,061 cells)
prev_labeled_file <- file.path(dirname(BASE), "R02_draft_v1/data/processed/seurat_object_labeled.rds")

if (file.exists(prev_labeled_file)) {
  cat("Loading previous labeled object:", prev_labeled_file, "\n")
  prev_obj <- readRDS(prev_labeled_file)
  cat(sprintf("  Previous object: %d cells, %d cell types\n",
              ncol(prev_obj), length(unique(na.omit(prev_obj$cell_type)))))

  # Build barcode -> cell_type map from previous object
  prev_labels <- setNames(prev_obj$cell_type, colnames(prev_obj))
  rm(prev_obj); gc(verbose = FALSE)

  # Match barcodes
  shared_barcodes <- intersect(colnames(sobj), names(prev_labels))
  cat(sprintf("  Shared barcodes: %d / %d (%.1f%%)\n",
              length(shared_barcodes), ncol(sobj),
              100 * length(shared_barcodes) / ncol(sobj)))

  sobj@meta.data$prev_celltype <- NA_character_
  if (length(shared_barcodes) > 0) {
    sobj@meta.data[shared_barcodes, "prev_celltype"] <- prev_labels[shared_barcodes]
  }
  cat("\nPrevious cell type labels transferred:\n")
  print(table(sobj$prev_celltype, useNA = "ifany"))

  # Cross-tabulate: module-score annotation vs previous labels
  cat("\nCross-tabulation: module-score annotation vs previous labels:\n")
  print(table(ModuleScore = sobj$celltype, Previous = sobj$prev_celltype, useNA = "no"))

  # Use previous labels as primary where available, fall back to module score
  sobj@meta.data$celltype_final <- ifelse(
    !is.na(sobj@meta.data$prev_celltype),
    sobj@meta.data$prev_celltype,
    sobj@meta.data$celltype
  )
  # Keep module-score annotation as secondary
  sobj@meta.data$celltype_module <- sobj@meta.data$celltype

  # Replace celltype with final consensus
  sobj@meta.data$celltype <- sobj@meta.data$celltype_final

  cat("\nFinal cell type labels (previous + module-score fallback):\n")
  print(sort(table(sobj$celltype), decreasing = TRUE))
} else {
  cat("Previous labeled object not found:", prev_labeled_file, "\n")
  cat("Using module-score annotations only.\n")
}

# Update broad_celltype based on final celltype
broad_map_full <- c(
  Eosinophils            = "Granulocyte",
  Macrophages            = "Myeloid",
  "Endothelial cells"    = "Stromal",
  "Dendritic cells"      = "Myeloid",
  "Mesenchymal stem cells" = "Stromal",
  Pericytes              = "Stromal",
  "B cells"              = "Lymphoid",
  "NK cells"             = "Lymphoid",
  "Mesothelial cells"    = "Stromal",
  "Mast cells"           = "Granulocyte",
  Adipocytes             = "Stromal",
  Neutrophils            = "Granulocyte",
  "Epithelial cells"     = "Stromal",
  # Module-score fallbacks
  Eosinophil   = "Granulocyte",
  Macrophage   = "Myeloid",
  T_cell       = "Lymphoid",
  B_cell       = "Lymphoid",
  NK_cell      = "Lymphoid",
  Dendritic    = "Myeloid",
  Neutrophil   = "Granulocyte",
  Other_immune = "Other"
)
sobj@meta.data$broad_celltype <- broad_map_full[sobj$celltype]
sobj@meta.data$broad_celltype[is.na(sobj$broad_celltype)] <- "Other"

# ── Summary ────────────────────────────────────────────────────────────────────
cat("\n========================================\n")
cat("Cell type summary\n")
cat("========================================\n\n")

cat("Cells per cell type:\n")
ct_table <- sort(table(sobj$celltype), decreasing = TRUE)
print(ct_table)

cat("\nCells per broad cell type:\n")
print(sort(table(sobj$broad_celltype), decreasing = TRUE))

cat("\nCell type proportions per condition:\n")
prop_table <- prop.table(table(sobj$condition, sobj$celltype), margin = 1)
print(round(prop_table, 3))

# ── Step 4: Figures ────────────────────────────────────────────────────────────
cat("\n========================================\n")
cat("Step 4: Generating figures\n")
cat("========================================\n\n")

# UMAP colored by celltype (final labels)
cat("UMAP colored by celltype (final labels)...\n")

# Build color mapping for cell types present
ct_levels <- sort(unique(sobj$celltype))
ct_colors <- broad_celltype_colors[ct_levels]
# Fall back to grey for any unmatched types
ct_colors[is.na(ct_colors)] <- "grey50"

p1 <- DimPlot(sobj, group.by = "celltype", label = TRUE, repel = TRUE,
              pt.size = 0.3, cols = ct_colors) +
  ggtitle("Cell Type Annotation (Final)") + theme_umap()

save_pdf(p1, file.path(fig_dir, "umap_celltype.pdf"),
         width = 10, height = 8)

# UMAP colored by broad celltype
cat("UMAP colored by broad celltype...\n")
p1b <- DimPlot(sobj, group.by = "broad_celltype", label = TRUE, repel = TRUE,
               pt.size = 0.3) +
  ggtitle("Broad Cell Type") + theme_umap()

save_pdf(p1b, file.path(fig_dir, "umap_broad_celltype.pdf"),
         width = 10, height = 8)

# Stacked bar plot by condition
cat("Stacked bar by condition...\n")
ct_counts <- as.data.frame(table(sobj$condition, sobj$celltype))
colnames(ct_counts) <- c("condition", "celltype", "count")

p2 <- ggplot(ct_counts, aes(x = condition, y = count, fill = celltype)) +
  geom_bar(stat = "identity", position = "fill") +
  scale_fill_manual(values = ct_colors) +
  scale_y_continuous(labels = scales::percent_format()) +
  labs(title = "Cell Type Composition by Condition",
       x = "Condition", y = "Proportion", fill = "Cell Type") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 9))

save_pdf(p2, file.path(fig_dir, "barplot_celltype_by_condition.pdf"),
         width = 10, height = 7)

# Combined panel
p_combined <- p1 | p2
save_pdf(p_combined, file.path(fig_dir, "celltype_annotation_panel.pdf"),
         width = 18, height = 8)

# ── Save ───────────────────────────────────────────────────────────────────────
cat("\n========================================\n")
cat("Saving results\n")
cat("========================================\n\n")

# Save annotated object
out_rds <- file.path(BASE, "all_cells_annotated.rds")
cat("Saving:", out_rds, "\n")
saveRDS(sobj, out_rds)

# Save markers CSV
if (!is.null(markers)) {
  markers_file <- file.path(BASE, "celltype_markers.csv")
  cat("Saving:", markers_file, "\n")
  write.csv(markers, markers_file, row.names = FALSE)
}

cat(sprintf("\nFigures saved to: %s\n", fig_dir))
cat("Done.\n")
