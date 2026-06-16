#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 52_ambient_correction.R   (Reviewer #1 — is C0 just ambient macrophage RNA?)
# Cluster-aware DecontX on the all-cells object (uses final cell types as the
# contamination sources), then for the eosinophil C0-C5 clusters compare, raw vs
# decontaminated: macrophage ambient genes (C1q/Apoe/Lyz2/Folr2) vs MHC-II
# (H2-Aa/H2-Ab1/Cd74). Expectation if C0 is a real MHC-II state with ambient mac
# RNA: mac genes fall after decontamination, MHC-II persists.
# Input : $BASE/all_cells_annotated.rds, $BASE/eos_subclustered.rds
# Output: $BASE/revision/figures/cluster_robustness/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(celda); library(Matrix); library(dplyr) })

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out  <- file.path(BASE, "revision", "figures", "cluster_robustness")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

all <- readRDS(file.path(BASE, "all_cells_annotated.rds"))
if (inherits(all[["RNA"]], "Assay5")) all[["RNA"]] <- JoinLayers(all[["RNA"]])
ctcol <- intersect(c("celltype","broad_celltype","cell_type","annotation"), colnames(all@meta.data))[1]
z <- as.character(all[[ctcol]][,1])
counts <- GetAssayData(all, layer = "counts")
cat("All cells:", ncol(all), " types:", length(unique(z)), "\n")

# ── cluster-aware DecontX ─────────────────────────────────────────────────────
cat("Running DecontX (z = cell type)...\n")
res <- decontX(x = counts, z = z)
decont <- res$decontXcounts
contam <- res$contamination
names(contam) <- colnames(counts)
cat("Mean residual contamination:", round(mean(contam),3), "\n")

# ── eosinophil C0-C5 labels (by barcode) ──────────────────────────────────────
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
if ("C6" %in% levels(Idents(eos))) { eos <- subset(eos, idents="C6", invert=TRUE); Idents(eos) <- droplevels(Idents(eos)) }
eos_cl <- setNames(as.character(Idents(eos)), colnames(eos))
common <- intersect(names(eos_cl), colnames(counts))
cat("Eos cells mapped:", length(common), "\n")

lognorm <- function(m) { cs <- Matrix::colSums(m); cs[cs==0] <- 1; log1p(sweep(m, 2, cs, "/") * 1e4) }
raw_e   <- lognorm(counts[, common, drop=FALSE])
dec_e   <- lognorm(decont[, common, drop=FALSE])
cl      <- eos_cl[common]
mac_g <- intersect(c("C1qa","C1qb","Apoe","Lyz2","Folr2","F13a1"), rownames(counts))
mhc_g <- intersect(c("H2-Aa","H2-Ab1","Cd74","H2-Eb1"), rownames(counts))

meanby <- function(mat, genes) sapply(levels(factor(cl)), function(k)
  round(mean(Matrix::colMeans(mat[genes, cl==k, drop=FALSE])), 3))

cmp <- data.frame(
  cluster      = levels(factor(cl)),
  contam       = round(tapply(contam[common], cl, mean), 3),
  mac_raw      = meanby(raw_e, mac_g),
  mac_decont   = meanby(dec_e, mac_g),
  mhc2_raw     = meanby(raw_e, mhc_g),
  mhc2_decont  = meanby(dec_e, mhc_g)
) %>% mutate(mac_drop_pct = round(100*(mac_raw - mac_decont)/pmax(mac_raw,1e-6),0),
             mhc2_drop_pct = round(100*(mhc2_raw - mhc2_decont)/pmax(mhc2_raw,1e-6),0))
write.csv(cmp, file.path(out, "ambient_correction_by_cluster.csv"), row.names = FALSE)
cat("\n=== raw vs DecontX, eos clusters (mac ambient genes vs MHC-II) ===\n"); print(cmp)

# per-gene C0 detail
c0 <- cl == "C0"
detail <- data.frame(gene = c(mac_g, mhc_g),
  type = c(rep("mac_ambient", length(mac_g)), rep("MHC-II", length(mhc_g))),
  C0_raw = round(c(Matrix::rowMeans(raw_e[mac_g, c0, drop=FALSE]), Matrix::rowMeans(raw_e[mhc_g, c0, drop=FALSE])), 3),
  C0_decont = round(c(Matrix::rowMeans(dec_e[mac_g, c0, drop=FALSE]), Matrix::rowMeans(dec_e[mhc_g, c0, drop=FALSE])), 3))
write.csv(detail, file.path(out, "ambient_correction_C0_genes.csv"), row.names = FALSE)
cat("\n=== C0 per-gene raw vs decont ===\n"); print(detail)

cat("\nDone 52_ambient_correction.R\n")
