#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 72_monocle_integrated.R   (INTEGRATION-based trajectory — Fig 2G)
# Monocle3 pseudotime on the ComBat-integrated eosinophils (i0-i5), using the
# ComBat UMAP as the reduced embedding. Reports pseudotime ordering across the
# integrated clusters for comparison with the published trajectory.
# Input : $REDO/eos_combat010.rds ; Output: $REDO/monocle/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(monocle3); library(Matrix) })

REDO <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/revision/redo"
out  <- file.path(REDO, "monocle"); dir.create(out, recursive=TRUE, showWarnings=FALSE)
eos <- readRDS(file.path(REDO, "eos_combat010.rds")); DefaultAssay(eos) <- "decontX"

expr <- GetAssayData(eos, layer="counts")
cd <- data.frame(row.names=colnames(eos), cluster=as.character(eos$new_cluster),
                 condition=eos$condition, sex=eos$sex, genotype=eos$genotype)
gd <- data.frame(row.names=rownames(eos), gene_short_name=rownames(eos))
cds <- new_cell_data_set(expr, cell_metadata=cd, gene_metadata=gd)
cds <- preprocess_cds(cds, num_dim=30)
# inject the ComBat UMAP as the reduced-dim used for graph learning
um <- Embeddings(eos, "umap_combat"); colnames(um) <- c("UMAP_1","UMAP_2")
reducedDims(cds)[["UMAP"]] <- um[colnames(cds), ]
cds <- cluster_cells(cds, reduction_method="UMAP")
cds <- learn_graph(cds, use_partition=FALSE)

# root: cluster with highest Klf2/Ccr3 (most circulating-like / C3-like = i4)
root_cl <- "i4"
root_cells <- colnames(cds)[cds@colData$cluster == root_cl]
cds <- order_cells(cds, root_cells = root_cells)
pt <- pseudotime(cds)
eos$pseudotime_combat <- pt[colnames(eos)]

df <- data.frame(cluster=as.character(eos$new_cluster), pseudotime=eos$pseudotime_combat)
sm <- aggregate(pseudotime ~ cluster, df[is.finite(df$pseudotime), ], function(x) round(median(x),3))
write.csv(sm, file.path(out, "pseudotime_by_integrated_cluster.csv"), row.names=FALSE)
cat("=== median pseudotime by integrated cluster (root", root_cl, ") ===\n"); print(sm)
saveRDS(cds, file.path(out, "monocle_integrated_cds.rds"))

pdf(file.path(out, "monocle_integrated_umap.pdf"), width=7, height=6)
print(plot_cells(cds, color_cells_by="cluster", label_groups_by_cluster=FALSE,
                 label_leaves=FALSE, label_branch_points=FALSE))
print(plot_cells(cds, color_cells_by="pseudotime", label_leaves=FALSE, label_branch_points=FALSE))
dev.off()
cat("\nDone 72.\n")
