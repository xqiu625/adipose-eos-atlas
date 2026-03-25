#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --job-name=monocle3
#SBATCH --output=monocle3_%j.out

# =============================================================================
# 11_monocle3.R
# Monocle3 trajectory / pseudotime analysis on eosinophil subclusters
# Input:  $BASE/eos_subclustered.rds
# Output: monocle3_cds.rds, trajectory_genes.csv, PDFs
# =============================================================================

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
outdir <- file.path(BASE, "monocle")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

# --- Load data ---------------------------------------------------------------
cat("Loading eos_subclustered.rds ...\n")
eos <- tryCatch(
  readRDS(file.path(BASE, "eos_subclustered.rds")),
  error = function(e) { cat("ERROR: Cannot load input:", e$message, "\n"); quit(status = 1) }
)
cat("  Cells:", ncol(eos), " Clusters:", paste(sort(unique(Idents(eos))), collapse = ", "), "\n")

# --- Convert Seurat to Monocle3 cell_data_set --------------------------------
cat("Converting Seurat object to cell_data_set ...\n")
cds <- tryCatch(
  as.cell_data_set(eos),
  error = function(e) { cat("ERROR converting to CDS:", e$message, "\n"); quit(status = 1) }
)
cat("  CDS created with", ncol(cds), "cells and", nrow(cds), "genes\n")

# --- Preprocess ---------------------------------------------------------------
cat("Preprocessing (num_dim = 20) ...\n")
cds <- preprocess_cds(cds, num_dim = 20)

# --- Transfer Seurat UMAP coordinates ----------------------------------------
cat("Transferring Seurat UMAP embeddings ...\n")
umap_coords <- Embeddings(eos, "umap")
# Ensure cell order matches
umap_coords <- umap_coords[colnames(cds), ]
reducedDims(cds)[["UMAP"]] <- umap_coords
cat("  UMAP coordinates set (", nrow(umap_coords), "cells )\n")

# --- Cluster cells ------------------------------------------------------------
cat("Clustering cells in UMAP space ...\n")
cds <- cluster_cells(cds, reduction_method = "UMAP")

# --- Learn trajectory graph ---------------------------------------------------
cat("Learning trajectory graph ...\n")
cds <- learn_graph(cds)
cat("  Trajectory graph learned\n")

# --- Identify root: cluster with lowest Siglecf expression (most immature) ---
cat("Identifying root based on Siglecf expression (lowest = most immature) ...\n")

# Helper: get_earliest_principal_node using marker gene
get_earliest_principal_node <- function(cds, marker_gene = "Siglecf") {
  # Get expression of marker gene per cell
  cell_ids <- colnames(cds)

  if (marker_gene %in% rownames(cds)) {
    marker_expr <- as.numeric(assay(cds)[marker_gene, ])
    names(marker_expr) <- cell_ids
  } else {
    cat("  WARNING: Siglecf not found in gene list. Using first cluster as root.\n")
    # Fall back: pick cluster with most cells as root
    cluster_tab <- table(clusters(cds))
    root_cluster <- names(which.max(cluster_tab))
    cell_ids_root <- names(clusters(cds))[clusters(cds) == root_cluster]

    closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
    closest_vertex <- as.numeric(closest_vertex[cell_ids_root, ])
    root_pr_node <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
      which.max(table(closest_vertex))
    ]
    return(root_pr_node)
  }

  # Compute mean expression per cluster — root = lowest Siglecf (most immature/progenitor-like)
  cell_clusters <- clusters(cds)
  cluster_means <- tapply(marker_expr, cell_clusters, mean)
  root_cluster <- names(which.min(cluster_means))
  cat("  Root cluster (lowest Siglecf = most immature):", root_cluster, "\n")

  # Map root cluster cells to closest principal graph node
  cell_ids_root <- names(cell_clusters)[cell_clusters == root_cluster]
  closest_vertex <- cds@principal_graph_aux[["UMAP"]]$pr_graph_cell_proj_closest_vertex
  closest_vertex <- as.numeric(closest_vertex[cell_ids_root, ])
  root_pr_node <- igraph::V(principal_graph(cds)[["UMAP"]])$name[
    which.max(table(closest_vertex))
  ]
  cat("  Root principal graph node:", root_pr_node, "\n")
  return(root_pr_node)
}

root_node <- get_earliest_principal_node(cds, marker_gene = "Siglecf")

cat("Ordering cells by pseudotime ...\n")
cds <- order_cells(cds, root_pr_nodes = root_node)
cat("  Pseudotime range:", round(min(pseudotime(cds)), 2), "-",
    round(max(pseudotime(cds)[is.finite(pseudotime(cds))]), 2), "\n")

# --- Graph test: genes varying along trajectory ------------------------------
cat("Running graph_test for trajectory-variable genes ...\n")
traj_genes <- tryCatch({
  gt_res <- graph_test(cds, neighbor_graph = "principal_graph", cores = 4)
  gt_res <- gt_res %>%
    filter(q_value < 0.05) %>%
    arrange(q_value)
  cat("  Significant trajectory genes:", nrow(gt_res), "\n")
  gt_res
}, error = function(e) {
  cat("  WARNING: graph_test failed:", e$message, "\n")
  NULL
})

# Save trajectory genes
if (!is.null(traj_genes) && nrow(traj_genes) > 0) {
  write.csv(traj_genes, file.path(outdir, "trajectory_genes.csv"), row.names = TRUE)
  cat("  Saved trajectory_genes.csv\n")
} else {
  cat("  No trajectory genes to save\n")
}

# --- Figures ------------------------------------------------------------------
cat("Generating figures ...\n")

# Store Seurat cluster labels in CDS colData for plotting
colData(cds)$seurat_cluster <- as.character(Idents(eos)[colnames(cds)])
if ("condition" %in% colnames(eos@meta.data)) {
  colData(cds)$condition <- eos$condition[colnames(cds)]
}

# Figure 1: UMAP by pseudotime
pdf(file.path(outdir, "monocle3_pseudotime_umap.pdf"), width = 8, height = 6)
tryCatch({
  p <- plot_cells(cds,
                  color_cells_by = "pseudotime",
                  label_cell_groups = FALSE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  graph_label_size = 3,
                  cell_size = 0.5) +
    ggtitle("Eosinophil pseudotime") +
    theme_umap()
  print(p)
}, error = function(e) cat("  Warning: pseudotime UMAP failed:", e$message, "\n"))
dev.off()

# Figure 2: UMAP by cluster
pdf(file.path(outdir, "monocle3_cluster_umap.pdf"), width = 8, height = 6)
tryCatch({
  p <- plot_cells(cds,
                  color_cells_by = "seurat_cluster",
                  label_cell_groups = TRUE,
                  label_leaves = FALSE,
                  label_branch_points = FALSE,
                  cell_size = 0.5) +
    scale_color_manual(values = eos_colors) +
    ggtitle("Eosinophil clusters with trajectory") +
    theme_umap()
  print(p)
}, error = function(e) cat("  Warning: cluster UMAP failed:", e$message, "\n"))
dev.off()

# Figure 3: Pseudotime boxplot by cluster
pdf(file.path(outdir, "monocle3_pseudotime_boxplot.pdf"), width = 8, height = 5)
tryCatch({
  pt_df <- data.frame(
    pseudotime = pseudotime(cds),
    cluster = colData(cds)$seurat_cluster,
    stringsAsFactors = FALSE
  )
  pt_df <- pt_df[is.finite(pt_df$pseudotime), ]

  p <- ggplot(pt_df, aes(x = reorder(cluster, pseudotime, FUN = median),
                          y = pseudotime, fill = cluster)) +
    geom_boxplot(outlier.size = 0.3) +
    scale_fill_manual(values = eos_colors) +
    labs(x = "Cluster", y = "Pseudotime",
         title = "Pseudotime distribution by cluster") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
  print(p)
}, error = function(e) cat("  Warning: pseudotime boxplot failed:", e$message, "\n"))
dev.off()

# Figure 4: Pseudotime boxplot by condition (if condition exists)
if ("condition" %in% colnames(colData(cds))) {
  pdf(file.path(outdir, "monocle3_pseudotime_by_condition.pdf"), width = 10, height = 5)
  tryCatch({
    pt_df <- data.frame(
      pseudotime = pseudotime(cds),
      condition = colData(cds)$condition,
      stringsAsFactors = FALSE
    )
    pt_df <- pt_df[is.finite(pt_df$pseudotime), ]

    p <- ggplot(pt_df, aes(x = condition, y = pseudotime, fill = condition)) +
      geom_boxplot(outlier.size = 0.3) +
      scale_fill_manual(values = cond_colors) +
      labs(x = "Condition", y = "Pseudotime",
           title = "Pseudotime distribution by condition") +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            legend.position = "none")
    print(p)
  }, error = function(e) cat("  Warning: condition boxplot failed:", e$message, "\n"))
  dev.off()
}

# Figure 5: Top trajectory-variable genes heatmap
if (!is.null(traj_genes) && nrow(traj_genes) >= 10) {
  pdf(file.path(outdir, "monocle3_trajectory_gene_heatmap.pdf"), width = 10, height = 12)
  tryCatch({
    top_genes <- head(rownames(traj_genes), 50)
    # Ensure genes are in the CDS
    top_genes <- top_genes[top_genes %in% rownames(cds)]

    if (length(top_genes) >= 5) {
      pt_matrix <- as.matrix(assay(cds)[top_genes, , drop = FALSE])
      # Order cells by pseudotime
      pt_order <- order(pseudotime(cds))
      pt_order <- pt_order[is.finite(pseudotime(cds)[pt_order])]
      pt_matrix <- pt_matrix[, pt_order, drop = FALSE]

      # Smooth and scale for visualization
      pt_matrix <- t(scale(t(pt_matrix)))
      pt_matrix[pt_matrix > 3] <- 3
      pt_matrix[pt_matrix < -3] <- -3

      # Use base heatmap if pheatmap not available
      if (requireNamespace("pheatmap", quietly = TRUE)) {
        pheatmap::pheatmap(pt_matrix,
                           cluster_cols = FALSE,
                           cluster_rows = TRUE,
                           show_colnames = FALSE,
                           color = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                           main = "Top trajectory-variable genes ordered by pseudotime")
      } else {
        heatmap(pt_matrix, Colv = NA, scale = "none",
                labCol = NA,
                col = colorRampPalette(c("navy", "white", "firebrick3"))(100),
                main = "Top trajectory-variable genes")
      }
    }
  }, error = function(e) cat("  Warning: gene heatmap failed:", e$message, "\n"))
  dev.off()
}

# --- Save CDS ----------------------------------------------------------------
cat("Saving monocle3 CDS ...\n")
saveRDS(cds, file.path(outdir, "monocle3_cds.rds"))
cat("  Saved monocle3_cds.rds\n")

cat("\n========== Monocle3 trajectory analysis complete ==========\n")
cat("Output directory:", outdir, "\n")
sessionInfo()
