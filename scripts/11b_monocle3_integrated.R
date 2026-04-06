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
  library(monocle3)
  library(Seurat)
  library(SeuratWrappers)
  library(ggplot2)
  library(dplyr)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
outdir <- file.path(BASE, "monocle_threeway")
fig_dir <- file.path(outdir, "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- 1. Load three-way integrated object -------------------------------------
cat("Loading three-way integrated object...\n")
integrated_path <- file.path(BASE, "integration", "three_way", "three_way_integrated.rds")
combined <- readRDS(integrated_path)

if (inherits(combined[["RNA"]], "Assay5")) {
  combined[["RNA"]] <- JoinLayers(combined[["RNA"]])
}
cat(sprintf("  Loaded: %d cells\n", ncol(combined)))

# --- 2. Convert to Monocle3 CDS ----------------------------------------------
cat("Converting to Monocle3 cell_data_set...\n")
cds <- as.cell_data_set(combined)

# Transfer UMAP coordinates
umap_coords <- Embeddings(combined, "umap")
umap_coords <- umap_coords[colnames(cds), ]
reducedDims(cds)[["UMAP"]] <- umap_coords

# --- 3. Preprocess and learn trajectory ---------------------------------------
cat("Preprocessing and learning trajectory...\n")
cds <- preprocess_cds(cds, num_dim = 30)
cds <- cluster_cells(cds, reduction_method = "UMAP")
cds <- learn_graph(cds)

# --- 4. Set root = Gurtner Bone Marrow cells ----------------------------------
cat("Setting root to Gurtner Bone Marrow cells...\n")
bm_mask <- colData(cds)$organ == "Bone Marrow" | 
           (colData(cds)$dataset == "Gurtner_2022_Blood_BM" & 
            colData(cds)$organ %in% c("Bone Marrow", "Blood"))

if (sum(bm_mask) == 0) {
  cat("WARNING: No Gurtner BM cells found. Using first Gurtner cell as fallback.\n")
  bm_mask <- colData(cds)$dataset == "Gurtner_2022_Blood_BM"
}

# Prefer highest CytoTRACE score if available (most immature)
if ("ct_score" %in% colnames(colData(cds))) {
  bm_ct <- colData(cds)$ct_score[bm_mask]
  root_idx <- which.max(bm_ct)
  root_cell <- colnames(cds)[bm_mask][root_idx]
} else {
  root_cell <- colnames(cds)[bm_mask][1]
}

# Map to principal graph node
closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
root_pr_node <- closest_vertex[root_cell, ]

cds <- order_cells(cds, root_pr_nodes = root_pr_node)
cat(sprintf("  Root cell set: %s (Gurtner BM)\n", root_cell))

# --- 5. Generate figures ------------------------------------------------------
cat("Generating figures...\n")

# Pseudotime on UMAP
pdf(file.path(fig_dir, "monocle3_threeway_pseudotime.pdf"), width = 9, height = 7)
p <- plot_cells(cds, color_cells_by = "pseudotime",
                label_cell_groups = FALSE,
                label_leaves = FALSE,
                label_branch_points = FALSE,
                cell_size = 0.6) +
  ggtitle("Pseudotime Trajectory (Root = Gurtner Bone Marrow)") +
  theme_umap()
print(p)
dev.off()

# By dataset
pdf(file.path(fig_dir, "monocle3_threeway_dataset.pdf"), width = 9, height = 7)
p <- plot_cells(cds, color_cells_by = "dataset", label_cell_groups = TRUE, cell_size = 0.6) +
  ggtitle("Dataset") + theme_umap()
print(p)
dev.off()

# By organ/tissue (shows Hu eWAT, scWAT, BAT separately)
pdf(file.path(fig_dir, "monocle3_threeway_organ.pdf"), width = 9, height = 7)
p <- plot_cells(cds, color_cells_by = "organ", label_cell_groups = TRUE, cell_size = 0.6) +
  ggtitle("Tissue / Organ") + theme_umap()
print(p)
dev.off()

# By original cluster
pdf(file.path(fig_dir, "monocle3_threeway_cluster.pdf"), width = 9, height = 7)
p <- plot_cells(cds, color_cells_by = "original_cluster", label_cell_groups = TRUE, cell_size = 0.6) +
  ggtitle("Original Cluster") + theme_umap()
print(p)
dev.off()

# Pseudotime boxplot by dataset
pdf(file.path(fig_dir, "monocle3_threeway_pseudotime_by_dataset.pdf"), width = 10, height = 6)
pt_df <- data.frame(
  pseudotime = pseudotime(cds),
  dataset = colData(cds)$dataset,
  stringsAsFactors = FALSE
)
pt_df <- pt_df[is.finite(pt_df$pseudotime), ]
ggplot(pt_df, aes(x = dataset, y = pseudotime, fill = dataset)) +
  geom_boxplot(outlier.size = 0.5) +
  labs(title = "Pseudotime by Dataset") +
  theme_pub()
dev.off()

# Save the CDS
saveRDS(cds, file.path(outdir, "monocle3_threeway_cds.rds"))
cat("Saved monocle3_threeway_cds.rds\n")

cat("\n=== Three-way Monocle3 analysis complete ===\n")
cat("Root: Gurtner Bone Marrow\n")
cat("Output directory:", outdir, "\n")