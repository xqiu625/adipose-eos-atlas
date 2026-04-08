#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 25_monocle3_threeway.R
# Monocle3 trajectory on three-way integrated data:
# Qiu HFD adipose + Hu healthy adipose (by tissue) + Gurtner blood/BM
# Root = Gurtner Bone Marrow cells (most immature)
# ==============================================================================

suppressPackageStartupMessages({
  library(Biobase)
  library(SingleCellExperiment)
  library(monocle3)
  library(Seurat)
  library(ggplot2)
  library(dplyr)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
outdir <- file.path(BASE, "monocle_threeway")
fig_dir <- file.path(outdir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load data (use saved CDS if available, otherwise build from scratch) --
# Force recompute trajectory (root cell was wrong in previous run)
SKIP_TO_FIGURES <- FALSE

if (!SKIP_TO_FIGURES) {
cat("Loading three-way integrated object...\n")
integrated_path <- file.path(BASE, "integration", "three_way", "three_way_integrated.rds")
combined <- readRDS(integrated_path)

if (inherits(combined[["RNA"]], "Assay5")) {
  combined[["RNA"]] <- JoinLayers(combined[["RNA"]])
}
cat(sprintf("  Loaded: %d cells\n", ncol(combined)))

# --- 2. Convert to Monocle3 CDS (manual — avoids SeuratWrappers version issues)
cat("Converting to Monocle3 cell_data_set...\n")
combined <- JoinLayers(combined)
expr_mat <- GetAssayData(combined, layer = "counts")
cell_meta <- combined@meta.data
gene_meta <- data.frame(gene_short_name = rownames(combined), row.names = rownames(combined))
cds <- new_cell_data_set(expr_mat,
                         cell_metadata = cell_meta,
                         gene_metadata = gene_meta)

# Transfer UMAP coordinates
umap_coords <- Embeddings(combined, "umap")
umap_coords <- umap_coords[colnames(cds), ]
reducedDims(cds)[["UMAP"]] <- umap_coords

# --- 3. Preprocess and learn trajectory ---------------------------------------
cat("Preprocessing and learning trajectory...\n")
cds <- preprocess_cds(cds, num_dim = 30)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

# --- 4. Set root = Gurtner Bone Marrow cells (BM ONLY, not Blood) -------------
cat("Setting root to Gurtner Bone Marrow cells...\n")
bm_mask <- colData(cds)$organ == "Bone Marrow"
cat(sprintf("  Bone Marrow cells: %d\n", sum(bm_mask)))

if (sum(bm_mask) == 0) {
  cat("  WARNING: No BM cells found. Trying any Gurtner cell.\n")
  bm_mask <- colData(cds)$dataset == "Gurtner_2022_Blood_BM"
}

# Pick BM cell with highest CytoTRACE (most progenitor-like)
if ("ct_score" %in% colnames(colData(cds))) {
  bm_ct <- colData(cds)$ct_score[bm_mask]
  root_idx <- which.max(bm_ct)
  root_cell <- colnames(cds)[bm_mask][root_idx]
} else {
  root_cell <- colnames(cds)[bm_mask][1]
}
cat(sprintf("  Selected root: %s (organ: %s)\n", root_cell,
            colData(cds)[root_cell, "organ"]))

# Map to principal graph node
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
root_pr_node <- paste0("Y_", closest_vertex[root_cell, ])
cat(sprintf("  Root cell: %s -> graph node: %s\n", root_cell, root_pr_node))

# Verify node exists in graph; fallback to interactive-style root selection
graph_nodes <- igraph::V(principal_graph(cds)[["UMAP"]])$name
cat(sprintf("  Graph has %d nodes. Example: %s\n", length(graph_nodes),
            paste(head(graph_nodes, 5), collapse = ", ")))

if (root_pr_node %in% graph_nodes) {
  cds <- order_cells(cds, root_pr_nodes = root_pr_node)
  cat(sprintf("  Root set successfully: %s\n", root_pr_node))
} else {
  # Fallback: find graph node closest to BM cells centroid
  cat("  Direct node mapping failed. Finding closest graph node to BM centroid...\n")
  bm_umap <- reducedDims(cds)[["UMAP"]][bm_mask, , drop = FALSE]
  bm_centroid <- colMeans(bm_umap)

  # Get graph node positions
  graph_coords <- cds@principal_graph_aux[["UMAP"]]$dp_mst
  if (!is.null(graph_coords)) {
    dists <- apply(t(graph_coords), 1, function(x) sqrt(sum((x - bm_centroid)^2)))
    root_pr_node <- names(which.min(dists))
    cds <- order_cells(cds, root_pr_nodes = root_pr_node)
    cat(sprintf("  Root set via centroid: %s\n", root_pr_node))
  } else {
    cat("  ERROR: Cannot find graph node positions. Skipping order_cells.\n")
  }
}

} # end if (!SKIP_TO_FIGURES)

# --- 5. Create unified group labels -------------------------------------------
cat("Creating group labels...\n")

ds <- as.character(colData(cds)$dataset)
org <- as.character(colData(cds)$organ)
cl <- as.character(colData(cds)$original_cluster)

group_label <- rep(NA_character_, ncol(cds))

# Qiu: C0-C5
qiu_idx <- which(ds == "Qiu_2025_Adipose")
group_label[qiu_idx] <- cl[qiu_idx]

# Gurtner: BM and Blood only (main figure)
bm_idx <- which(org == "Bone Marrow")
group_label[bm_idx] <- "Gurtner BM"
blood_idx <- which(org == "Blood")
group_label[blood_idx] <- "Gurtner Blood"

# Hu: depot from barcode prefix (supplement only)
hu_idx <- which(ds == "Hu_2026_Healthy")
hu_barcodes <- colnames(cds)[hu_idx]
hu_depot <- sub("^Hu_([^_]+)_.*", "\\1", hu_barcodes)
group_label[hu_idx] <- paste0("Hu ", hu_depot)

# Other Gurtner organs → NA (excluded from main figure)
colData(cds)$group_label <- group_label
cat("  All group sizes:\n")
print(table(group_label, useNA = "ifany"))

# Main figure groups: BM + Blood + C0-C5 (no Hu, no other Gurtner)
main_groups <- c("Gurtner BM", "Gurtner Blood",
                 "C0", "C1", "C2", "C3", "C4", "C5")
colData(cds)$main_group <- ifelse(group_label %in% main_groups, group_label, NA)

# Supplement groups: adds Hu depots
supp_groups <- c(main_groups, "Hu eWAT", "Hu scWAT", "Hu BAT")
colData(cds)$supp_group <- ifelse(group_label %in% supp_groups, group_label, NA)

cat("  Main figure cells:", sum(!is.na(colData(cds)$main_group)), "\n")
cat("  Supplement cells:", sum(!is.na(colData(cds)$supp_group)), "\n")

# --- Group colors ---------------------------------------------------------------
group_colors <- c(
  "C0" = "#6CA6CD", "C1" = "#E84A5F", "C2" = "#FF6F1F",
  "C3" = "#CBAACB", "C4" = "#B03A48", "C5" = "#5BA4A4",
  "Gurtner BM"    = "#2E8B57",
  "Gurtner Blood" = "#DAA520",
  "Hu eWAT"  = "#FF8C00",
  "Hu scWAT" = "#8B4513",
  "Hu BAT"   = "#DC143C"
)

# --- 6. Generate figures ------------------------------------------------------
cat("Generating figures...\n")

# Helper: shared monocle3 theme
monocle_theme <- theme(
  plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
  legend.text = element_text(size = 11),
  legend.position = "right",
  axis.title = element_text(size = 12, face = "bold"),
  axis.text = element_blank(),
  axis.ticks = element_blank(),
  panel.background = element_blank()
)

# === MAIN FIGURE: BM + Blood + C0-C5 only ===
cat("  Main figure (BM + Blood + C0-C5)...\n")
main_cells <- colnames(cds)[!is.na(colData(cds)$main_group)]
cds_main <- cds[, main_cells]

pdf(file.path(fig_dir, "monocle3_main_groups.pdf"), width = 10, height = 8)
p <- plot_cells(cds_main, color_cells_by = "main_group",
                label_cell_groups = TRUE, label_leaves = FALSE,
                label_branch_points = FALSE, graph_label_size = 4,
                cell_size = 0.5, group_label_size = 5) +
  scale_color_manual(values = group_colors, name = "") +
  ggtitle("Monocle3 Trajectory") + monocle_theme
print(p)
dev.off()

pdf(file.path(fig_dir, "monocle3_main_pseudotime.pdf"), width = 10, height = 8)
p <- plot_cells(cds_main, color_cells_by = "pseudotime",
                label_cell_groups = FALSE, label_leaves = FALSE,
                label_branch_points = FALSE, cell_size = 0.5) +
  scale_color_viridis_c(name = "Pseudotime", option = "inferno") +
  ggtitle("Pseudotime (Root = Gurtner BM)") + monocle_theme
print(p)
dev.off()

# Main combined (stacked)
pdf(file.path(fig_dir, "monocle3_main_combined.pdf"), width = 10, height = 16)
p1 <- plot_cells(cds_main, color_cells_by = "main_group",
                 label_cell_groups = TRUE, label_leaves = FALSE,
                 label_branch_points = FALSE, graph_label_size = 4,
                 cell_size = 0.5, group_label_size = 5) +
  scale_color_manual(values = group_colors, name = "") +
  ggtitle("Monocle3 Trajectory") + monocle_theme
p2 <- plot_cells(cds_main, color_cells_by = "pseudotime",
                 label_cell_groups = FALSE, label_leaves = FALSE,
                 label_branch_points = FALSE, cell_size = 0.5) +
  scale_color_viridis_c(name = "Pseudotime", option = "inferno") +
  ggtitle("Pseudotime") + monocle_theme
print(p1 / p2)
dev.off()

# Main pseudotime boxplot
pdf(file.path(fig_dir, "monocle3_main_pseudotime_boxplot.pdf"), width = 10, height = 6)
pt_df <- data.frame(
  pseudotime = pseudotime(cds_main),
  group = colData(cds_main)$main_group,
  stringsAsFactors = FALSE
)
pt_df <- pt_df[is.finite(pt_df$pseudotime) & !is.na(pt_df$group), ]
med_order <- pt_df %>% group_by(group) %>% summarise(med = median(pseudotime)) %>% arrange(med)
pt_df$group <- factor(pt_df$group, levels = med_order$group)
p <- ggplot(pt_df, aes(x = group, y = pseudotime, fill = group)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Pseudotime by Group (BM + Blood + C0-C5)", x = "", y = "Pseudotime") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
print(p)
dev.off()

# === SUPPLEMENT: adds Hu depots ===
cat("  Supplement (+ Hu eWAT/scWAT/BAT)...\n")
supp_cells <- colnames(cds)[!is.na(colData(cds)$supp_group)]
cds_supp <- cds[, supp_cells]

pdf(file.path(fig_dir, "monocle3_supp_groups.pdf"), width = 10, height = 8)
p <- plot_cells(cds_supp, color_cells_by = "supp_group",
                label_cell_groups = TRUE, label_leaves = FALSE,
                label_branch_points = FALSE, graph_label_size = 4,
                cell_size = 0.5, group_label_size = 5) +
  scale_color_manual(values = group_colors, name = "") +
  ggtitle("Monocle3 Trajectory (with Hu Adipose)") + monocle_theme
print(p)
dev.off()

pdf(file.path(fig_dir, "monocle3_supp_pseudotime_boxplot.pdf"), width = 12, height = 6)
pt_df <- data.frame(
  pseudotime = pseudotime(cds_supp),
  group = colData(cds_supp)$supp_group,
  stringsAsFactors = FALSE
)
pt_df <- pt_df[is.finite(pt_df$pseudotime) & !is.na(pt_df$group), ]
med_order <- pt_df %>% group_by(group) %>% summarise(med = median(pseudotime)) %>% arrange(med)
pt_df$group <- factor(pt_df$group, levels = med_order$group)
p <- ggplot(pt_df, aes(x = group, y = pseudotime, fill = group)) +
  geom_boxplot(outlier.size = 0.5) +
  scale_fill_manual(values = group_colors) +
  labs(title = "Pseudotime by Group (with Hu Adipose)", x = "", y = "Pseudotime") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1), legend.position = "none")
print(p)
dev.off()

# Save the CDS
saveRDS(cds, file.path(outdir, "monocle3_threeway_cds.rds"))
cat("Saved monocle3_threeway_cds.rds\n")

cat("\n=== Three-way Monocle3 analysis complete ===\n")
cat("Root: Gurtner Bone Marrow\n")
cat("Output directory:", outdir, "\n")