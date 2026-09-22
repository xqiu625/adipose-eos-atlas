#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"; .libPaths(c(lib_path, .libPaths()))
# Per-cell data for the integration-based Fig 2: UMAP, cluster, MHC-II program
# score, pseudotime. Output one CSV for figure-making.
suppressPackageStartupMessages({ library(Seurat); library(monocle3) })
REDO <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/revision/redo"
out  <- file.path(REDO, "fig2"); dir.create(out, showWarnings=FALSE)
eos <- readRDS(file.path(REDO, "eos_combat010.rds")); DefaultAssay(eos) <- "decontX"
mhc <- intersect(c("H2-Aa","H2-Ab1","H2-Eb1","Cd74"), rownames(eos))
eos <- AddModuleScore(eos, features=list(mhc), name="MHCII", seed=1)
um <- Embeddings(eos, "umap_combat")
cds <- readRDS(file.path(REDO, "monocle/monocle_integrated_cds.rds"))
pt <- pseudotime(cds)
df <- data.frame(barcode=colnames(eos),
                 UMAP_1=um[,1], UMAP_2=um[,2],
                 cluster=as.character(eos$new_cluster),
                 sex=eos$sex, condition=eos$condition,
                 MHCII=round(eos$MHCII1,3),
                 pseudotime=round(pt[colnames(eos)],3))
write.csv(df, file.path(out, "fig2_percell.csv"), row.names=FALSE)
# top markers per cluster for the dotplot (decont)
mk <- read.csv(file.path(REDO,"combat_core/combat_cluster_markers.csv"))
cat("wrote fig2_percell.csv:", nrow(df), "cells; clusters:", paste(sort(unique(df$cluster)),collapse=","), "\n")
