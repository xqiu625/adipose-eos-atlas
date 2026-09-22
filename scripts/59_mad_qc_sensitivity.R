#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))
if (nzchar(Sys.getenv("RERUN_R_LIB"))) {
  .rerun_r_lib_resolved <- normalizePath(.libPaths()[1], mustWork = FALSE)
  .rerun_r_lib_expected <- normalizePath(Sys.getenv("RERUN_R_LIB"), mustWork = FALSE)
  if (!identical(.rerun_r_lib_resolved, .rerun_r_lib_expected)) {
    stop(sprintf(
      "RERUN_R_LIB audit guard: .libPaths()[1] resolved to '%s' but RERUN_R_LIB='%s' -- the shared/legacy library would take precedence over the pinned isolated library. Refusing to proceed.",
      .rerun_r_lib_resolved, .rerun_r_lib_expected))
  }
  cat("=== Pinned R environment audit guard ===\n")
  cat("R.version.string:", R.version.string, "\n")
  cat("RERUN_R_LIB:", Sys.getenv("RERUN_R_LIB"), "\n")
  cat("CMDSTAN:", Sys.getenv("CMDSTAN"), "\n")
  cat(".libPaths():\n")
  print(.libPaths())
  cat("=========================================\n")
}

# ==============================================================================
# 59_mad_qc_sensitivity.R   (Reviewer / Handler-Gurtner 2026 — QC strategy check)
# We used fixed-threshold QC (nFeature 200-5000, %mt<=10, DecontX<=0.2, Scrublet).
# This checks robustness to an adaptive median-absolute-deviation (MAD) QC:
#   (1) per-sample MAD-derived bounds vs our fixed thresholds
#   (2) how many of the retained eosinophils MAD would flag as outliers
#   (3) re-cluster after MAD filtering and compare to C0-C5 (ARI)
# Input : $BASE/eos_subclustered.rds ; Output: $BASE/revision/figures/mad_qc/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out  <- file.path(BASE, "revision", "figures", "mad_qc"); dir.create(out, recursive=TRUE, showWarnings=FALSE)
set.seed(1)

eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
if ("C6" %in% levels(Idents(eos))) { eos <- subset(eos, idents="C6", invert=TRUE); Idents(eos)<-droplevels(Idents(eos)) }
if (inherits(eos[["RNA"]],"Assay5")) eos[["RNA"]] <- JoinLayers(eos[["RNA"]])
md <- eos@meta.data
md$cluster <- as.character(Idents(eos)); md$sample <- md$orig.ident

# ── MAD outlier flag (per sample) ─────────────────────────────────────────────
is_out <- function(x, nmads, type=c("both","higher","lower")) {
  type <- match.arg(type); m <- median(x); d <- mad(x)
  if (d == 0) return(rep(FALSE, length(x)))
  switch(type, both = x < m - nmads*d | x > m + nmads*d,
               higher = x > m + nmads*d, lower = x < m - nmads*d)
}
flag_mad <- function(nmads_lib=5, nmads_mt=3) {
  o <- logical(nrow(md))
  for (s in unique(md$sample)) {
    ix <- which(md$sample == s)
    o[ix] <- is_out(log1p(md$nCount_RNA[ix]), nmads_lib) |
             is_out(log1p(md$nFeature_RNA[ix]), nmads_lib) |
             is_out(md$percent.mt[ix], nmads_mt, "higher")
  }
  o
}
mad5 <- flag_mad(5,3)
mad3 <- flag_mad(3,3)
md$mad5_outlier <- mad5; md$mad3_outlier <- mad3

cat(sprintf("Retained eosinophils: %d\n", nrow(md)))
cat(sprintf("Flagged as outliers by MAD (5 MAD lib / 3 MAD mt): %d (%.1f%%)\n", sum(mad5), 100*mean(mad5)))
cat(sprintf("Flagged by stricter MAD (3 MAD / 3 MAD):           %d (%.1f%%)\n", sum(mad3), 100*mean(mad3)))
cat("\nMAD-flagged by cluster (5 MAD):\n"); print(round(100*tapply(mad5, md$cluster, mean), 1))

# ── per-sample MAD-derived bounds vs fixed thresholds ─────────────────────────
bounds <- md %>% group_by(sample) %>% summarise(
  nFeat_median = median(nFeature_RNA),
  nFeat_MAD_lower = round(expm1(median(log1p(nFeature_RNA)) - 5*mad(log1p(nFeature_RNA)))),
  nFeat_MAD_upper = round(expm1(median(log1p(nFeature_RNA)) + 5*mad(log1p(nFeature_RNA)))),
  mt_MAD_upper = round(median(percent.mt) + 3*mad(percent.mt), 2), .groups="drop")
bounds$fixed_lower <- 200; bounds$fixed_upper <- 5000; bounds$fixed_mt <- 10
write.csv(bounds, file.path(out, "mad_bounds_vs_fixed.csv"), row.names=FALSE)
cat("\n=== per-sample MAD bounds vs fixed (200-5000, mt<=10) ===\n"); print(as.data.frame(bounds))

# ── re-cluster after MAD filtering; ARI vs C0-C5 ──────────────────────────────
keep <- rownames(md)[!mad5]
sub <- DietSeurat(subset(eos, cells = keep), assays = "RNA")  # drop stale graphs/reductions
if (inherits(sub[["RNA"]],"Assay5")) sub[["RNA"]] <- JoinLayers(sub[["RNA"]])
DefaultAssay(sub) <- "RNA"
sub <- NormalizeData(sub, verbose=FALSE)
sub <- FindVariableFeatures(sub, nfeatures=2000, verbose=FALSE)
sub <- ScaleData(sub, verbose=FALSE)
sub <- RunPCA(sub, npcs=30, verbose=FALSE)
sub <- FindNeighbors(sub, reduction="pca", dims=1:20, verbose=FALSE)
sub <- FindClusters(sub, resolution=0.4, verbose=FALSE)
compute_ari <- function(a,b){ t<-table(a,b); n<-sum(t)
  sij<-sum(choose(t,2)); ai<-sum(choose(rowSums(t),2)); bj<-sum(choose(colSums(t),2)); cn<-choose(n,2)
  e<-ai*bj/cn; mx<-(ai+bj)/2; if(mx==e) 1 else (sij-e)/(mx-e) }
ari <- compute_ari(md[keep,"cluster"], sub$seurat_clusters)
nmi_tab <- table(md[keep,"cluster"], sub$seurat_clusters)
res <- data.frame(metric=c("eos_retained","mad5_flagged_pct","ari_MADrecluster_vs_C0C5",
                           "n_clusters_orig","n_clusters_MAD"),
                  value=c(nrow(md), round(100*mean(mad5),1), round(ari,3),
                          length(unique(md$cluster)), length(unique(sub$seurat_clusters))))
write.csv(res, file.path(out, "mad_qc_summary.csv"), row.names=FALSE)
write.csv(as.data.frame.matrix(nmi_tab), file.path(out, "mad_recluster_vs_C0C5_crosstab.csv"))
cat(sprintf("\n=== ARI(re-cluster after MAD filtering vs original C0-C5) = %.3f ===\n", ari))
cat("(high ARI -> clusters are robust to fixed-vs-MAD QC strategy)\n")
cat("\nDone 59_mad_qc_sensitivity.R\n")
