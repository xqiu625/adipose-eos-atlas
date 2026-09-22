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
# Resolution-matched MAD re-clustering: re-cluster the MAD-filtered eosinophils at
# the resolution that yields 6 clusters and compare to C0-C5 (removes the
# resolution confound from 59's fixed res=0.4, which gave 4 clusters).
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })
BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out  <- file.path(BASE, "revision", "figures", "mad_qc"); set.seed(1)

eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
if ("C6" %in% levels(Idents(eos))) { eos <- subset(eos, idents="C6", invert=TRUE); Idents(eos)<-droplevels(Idents(eos)) }
if (inherits(eos[["RNA"]],"Assay5")) eos[["RNA"]] <- JoinLayers(eos[["RNA"]])
md <- eos@meta.data; md$cluster <- as.character(Idents(eos)); md$sample <- md$orig.ident
is_out <- function(x,n,t="both"){ m<-median(x); d<-mad(x); if(d==0) return(rep(FALSE,length(x)))
  if(t=="higher") x>m+n*d else x<m-n*d | x>m+n*d }
o <- logical(nrow(md))
for (s in unique(md$sample)) { ix<-which(md$sample==s)
  o[ix] <- is_out(log1p(md$nCount_RNA[ix]),5) | is_out(log1p(md$nFeature_RNA[ix]),5) | is_out(md$percent.mt[ix],3,"higher") }
keep <- rownames(md)[!o]

sub <- DietSeurat(subset(eos, cells=keep), assays="RNA")
if (inherits(sub[["RNA"]],"Assay5")) sub[["RNA"]] <- JoinLayers(sub[["RNA"]])
sub <- NormalizeData(sub,verbose=FALSE); sub <- FindVariableFeatures(sub,nfeatures=2000,verbose=FALSE)
sub <- ScaleData(sub,verbose=FALSE); sub <- RunPCA(sub,npcs=30,verbose=FALSE)
sub <- FindNeighbors(sub,reduction="pca",dims=1:20,verbose=FALSE)
compute_ari <- function(a,b){ t<-table(a,b); n<-sum(t); sij<-sum(choose(t,2)); ai<-sum(choose(rowSums(t),2))
  bj<-sum(choose(colSums(t),2)); e<-ai*bj/choose(n,2); mx<-(ai+bj)/2; if(mx==e)1 else (sij-e)/(mx-e) }
rows <- list()
for (r in c(0.4,0.6,0.8,1.0,1.2)) {
  sub <- FindClusters(sub, resolution=r, verbose=FALSE)
  nk <- length(unique(sub$seurat_clusters)); ari <- compute_ari(md[keep,"cluster"], sub$seurat_clusters)
  rows[[as.character(r)]] <- c(res=r, n_clusters=nk, ARI_vs_C0C5=round(ari,3))
  cat(sprintf("res %.1f -> %d clusters, ARI vs C0-C5 = %.3f\n", r, nk, ari))
}
res <- as.data.frame(do.call(rbind, rows))
write.csv(res, file.path(out, "mad_recluster_resolution_sweep.csv"), row.names=FALSE)
best <- res[which.min(abs(res$n_clusters-6)), ]
cat(sprintf("\n=== matched: res %.1f gives %d clusters, ARI vs C0-C5 = %.3f ===\n", best$res, best$n_clusters, best$ARI_vs_C0C5))
cat("Done 59b.\n")
