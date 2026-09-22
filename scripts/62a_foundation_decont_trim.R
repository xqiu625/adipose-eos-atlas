#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 62a_foundation_decont_trim.R   (Integration+decontamination redo — STAGE 1)
# Build the decontaminated, C0-trimmed eosinophil matrix that all downstream
# re-analyses will use.
#   1. Cluster-aware DecontX on all cells (z = celltype; macrophages etc. are the
#      ambient source) -> per-cell decontaminated counts + per-cell contamination.
#   2. Subset to eosinophils (C0-C5).
#   3. Drop C0 cells with per-cell DecontX contamination > 0.30.
#   4. Export to h5ad for Stage 2 (ComBat integration + re-clustering, 62b):
#        X            = raw counts (for DESeq2 pseudobulk / SCENIC)
#        layers[decontX] = decontaminated counts (for ComBat / expression)
#        obs: orig.ident, condition, sex, genotype, infection, gem_reaction,
#             old_eos_cluster (original C0-C5), decontx_contam
# Input : all_cells_annotated.rds, eos_subclustered.rds
# Output: $BASE/revision/redo/eos_decont_trimmed.h5ad  (+ summary CSV)
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(celda); library(Matrix); library(reticulate) })

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out  <- file.path(BASE, "revision", "redo"); dir.create(out, recursive=TRUE, showWarnings=FALSE)
use_condaenv("/bigdata/godziklab/shared/Xinru/cellrank_env")
ad <- import("anndata"); np <- import("numpy"); sp <- import("scipy.sparse"); pd <- import("pandas")

# ── 1. cluster-aware DecontX on all cells ────────────────────────────────────
cat("Loading all_cells_annotated.rds...\n")
allc <- readRDS(file.path(BASE, "all_cells_annotated.rds"))
if (inherits(allc[["RNA"]],"Assay5")) allc[["RNA"]] <- JoinLayers(allc[["RNA"]])
ctcol <- intersect(c("celltype","broad_celltype"), colnames(allc@meta.data))[1]
z <- as.character(allc[[ctcol]][,1])
counts_all <- GetAssayData(allc, layer="counts")
cat("Running DecontX (z = celltype) on", ncol(allc), "cells...\n")
res <- decontX(x = counts_all, z = z)
decont_all <- res$decontXcounts
contam_all <- setNames(res$contamination, colnames(counts_all))
cat("  mean contamination (all):", round(mean(contam_all),3), "\n")

# ── 2. eosinophils + original C0-C5 labels ───────────────────────────────────
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
if ("C6" %in% levels(Idents(eos))) { eos <- subset(eos, idents="C6", invert=TRUE); Idents(eos)<-droplevels(Idents(eos)) }
eos_cl <- setNames(as.character(eos$eos_cluster), colnames(eos))
eos_bc <- intersect(names(eos_cl), colnames(counts_all))
cat("Eosinophils mapped:", length(eos_bc), "\n")

# ── 3. drop C0 cells with per-cell contamination > 0.30 ──────────────────────
contam_eos <- contam_all[eos_bc]; cl_eos <- eos_cl[eos_bc]
drop <- (cl_eos == "C0") & (contam_eos > 0.30)
cat(sprintf("C0 cells total: %d ; C0 with contam>0.30 (dropped): %d\n",
            sum(cl_eos=="C0"), sum(drop)))
cat("Per-cluster contam>0.30 counts (for reference):\n")
print(tapply(contam_eos > 0.30, cl_eos, sum))
keep_bc <- eos_bc[!drop]
cat(sprintf("Eosinophils retained: %d (from %d)\n", length(keep_bc), length(eos_bc)))

# ── 4. assemble + export h5ad ────────────────────────────────────────────────
raw_t    <- t(counts_all[, keep_bc, drop=FALSE])         # cells x genes
decont_t <- t(round(decont_all[, keep_bc, drop=FALSE]))  # round decont to integers
to_csc <- function(m) sp$csc_matrix(tuple(np$array(m@x), np$array(as.integer(m@i)), np$array(as.integer(m@p))),
                                    shape=tuple(as.integer(nrow(m)), as.integer(ncol(m))))
md <- eos@meta.data[keep_bc, ]
obs <- data.frame(
  orig.ident = md$orig.ident, condition = md$condition,
  sex = md$sex, genotype = md$genotype, infection = md$infection,
  gem_reaction = md$gem_reaction, old_eos_cluster = cl_eos[keep_bc],
  decontx_contam = round(contam_eos[keep_bc], 4), row.names = keep_bc)
obs[] <- lapply(obs, function(x) if (is.factor(x)) as.character(x) else x)
obs_pd <- pd$DataFrame(obs); obs_pd$index <- keep_bc
var_pd <- pd$DataFrame(index = rownames(counts_all))
adata <- ad$AnnData(X = to_csc(raw_t), obs = obs_pd, var = var_pd)
adata$layers["decontX"] <- to_csc(decont_t)
h5 <- file.path(out, "eos_decont_trimmed.h5ad")
adata$write_h5ad(h5)
cat("Wrote", h5, "shape", paste(adata$shape, collapse=" x "), "\n")

write.csv(data.frame(metric=c("eos_before","C0_before","C0_dropped_contam_gt_0.3","eos_after"),
                     value=c(length(eos_bc), sum(cl_eos=="C0"), sum(drop), length(keep_bc))),
          file.path(out, "stage1_trim_summary.csv"), row.names=FALSE)
cat("\nDone 62a.\n")
