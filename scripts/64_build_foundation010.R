#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 64_build_foundation010.R   (decont + >0.10 trim — downstream foundation)
# Build the eosinophil object for the re-run: decontaminated counts, per-cell
# DecontX contamination > 0.10 removed, PUBLISHED C0-C5 labels carried over
# (un-integrated primary; ComBat kept separately as sensitivity). Reuses the
# decont output from 62a (eos_decont_trimmed.h5ad: X=raw, layers[decontX]=decont,
# obs[decontx_contam], obs[old_eos_cluster]).
# Output: $REDO/eos_decont010.rds  (RNA assay = raw counts; decontX assay = decont;
#         Idents = C0-C5; standard Seurat processing on decont lognorm)
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix); library(reticulate) })

REDO <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/revision/redo"
use_condaenv("/bigdata/godziklab/shared/Xinru/cellrank_env")
adm <- import("anndata")

cat("Reading eos_decont_trimmed.h5ad...\n")
ad <- adm$read_h5ad(file.path(REDO, "eos_decont_trimmed.h5ad"))
obs <- ad$obs
bc  <- rownames(obs); genes <- rownames(ad$var)
raw    <- t(as(ad$X, "CsparseMatrix"));          dimnames(raw)    <- list(genes, bc)
decont <- t(as(ad$layers["decontX"], "CsparseMatrix")); dimnames(decont) <- list(genes, bc)

# >0.10 contamination filter
keep <- which(as.numeric(obs$decontx_contam) <= 0.10)
cat(sprintf("Cells: %d -> %d after DecontX contam <= 0.10 (removed %d)\n",
            length(bc), length(keep), length(bc)-length(keep)))
raw <- raw[, keep]; decont <- decont[, keep]; obs <- obs[keep, ]

cat("Per-cluster retained:\n"); print(table(obs$old_eos_cluster))

# build Seurat: RNA = raw counts (for DESeq2/SCENIC); decontX assay = decont
eos <- CreateSeuratObject(counts = raw, meta.data = obs)
eos[["decontX"]] <- CreateAssayObject(counts = decont)
eos$eos_cluster <- obs$old_eos_cluster
Idents(eos) <- factor(eos$eos_cluster, levels=paste0("C",0:5))

# standard processing on decontaminated expression for markers/embedding
DefaultAssay(eos) <- "decontX"
eos <- NormalizeData(eos, verbose=FALSE)
eos <- FindVariableFeatures(eos, nfeatures=2000, verbose=FALSE)
eos <- ScaleData(eos, verbose=FALSE)
eos <- RunPCA(eos, npcs=30, verbose=FALSE)
eos <- RunUMAP(eos, dims=1:20, verbose=FALSE)
# normalize RNA (raw) too for any raw-based steps
DefaultAssay(eos) <- "RNA"; eos <- NormalizeData(eos, verbose=FALSE)

saveRDS(eos, file.path(REDO, "eos_decont010.rds"))
write.csv(data.frame(cluster=names(table(eos$eos_cluster)),
                     n=as.integer(table(eos$eos_cluster))),
          file.path(REDO, "foundation010_cluster_counts.csv"), row.names=FALSE)
cat(sprintf("\nSaved eos_decont010.rds: %d cells x %d genes (assays: RNA raw, decontX)\n",
            ncol(eos), nrow(eos)))
cat("Done 64.\n")
