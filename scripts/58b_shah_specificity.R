#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))
# Specificity control for the Shah concordance: our EOSINOPHIL pseudobulk should
# correlate with Shah bulk eosinophils MORE than our macrophage/neutrophil
# pseudobulk does — otherwise the ~0.88 is just shared housekeeping. Also report
# on the most-variable genes (housekeeping removed), which is the sensitive test.
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

EXT <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/external/shah2023"
BASE<- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out <- file.path(BASE, "revision", "figures", "shah_concordance")

t2g <- read.table(file.path(EXT, "tx2gene.tsv"), sep="\t", stringsAsFactors=FALSE, quote="")
tx2g <- setNames(t2g$V2, t2g$V1)
sf <- list.files(file.path(EXT, "raw"), pattern="\\.sf\\.gz$", full.names=TRUE)
read_sf <- function(f){ d<-read.table(gzfile(f),header=TRUE,sep="\t"); d$g<-tx2g[sub("\\.\\d+$","",d$Name)]; d<-d[!is.na(d$g),]; tapply(d$NumReads,d$g,sum) }
shl <- lapply(sf, read_sf); gs <- Reduce(union, lapply(shl, names))
shah <- sapply(shl, function(v){x<-setNames(rep(0,length(gs)),gs); x[names(v)]<-v; x})
shah_eos <- rowMeans(shah)                     # mean across all 12 Shah eosinophil samples
shah_var <- apply(shah, 1, var)                # for HVG selection

# our pseudobulks: eos (C0-C5 object) + macrophage/neutrophil (all-cells)
eos <- readRDS(file.path(BASE,"eos_subclustered.rds")); if (inherits(eos[["RNA"]],"Assay5")) eos[["RNA"]]<-JoinLayers(eos[["RNA"]])
our_eos <- Matrix::rowSums(GetAssayData(eos,layer="counts")); rm(eos); gc()
ac <- readRDS(file.path(BASE,"all_cells_annotated.rds")); if (inherits(ac[["RNA"]],"Assay5")) ac[["RNA"]]<-JoinLayers(ac[["RNA"]])
acc <- GetAssayData(ac, layer="counts"); ct <- ac$celltype
our_mac    <- Matrix::rowSums(acc[, ct=="Macrophages", drop=FALSE])
our_neutro <- Matrix::rowSums(acc[, ct=="Neutrophils", drop=FALSE])
rm(ac, acc); gc()

cpm <- function(v) log1p(v/sum(v)*1e6)
common <- Reduce(intersect, list(names(shah_eos), names(our_eos), names(our_mac), names(our_neutro)))
common <- common[our_eos[common] > 0]
# HVG: top 2000 most-variable genes in Shah (removes invariant housekeeping)
hvg <- intersect(names(sort(shah_var, decreasing=TRUE))[1:3000], common)

S <- function(a,b,g) round(cor(cpm(a[g]), cpm(b[g]), method="spearman"), 3)
res <- data.frame(
  our_pseudobulk = c("Eosinophils (C0-C5)","Macrophages","Neutrophils"),
  vs_Shah_eos_allgenes = c(S(our_eos,shah_eos,common), S(our_mac,shah_eos,common), S(our_neutro,shah_eos,common)),
  vs_Shah_eos_top3000var = c(S(our_eos,shah_eos,hvg), S(our_mac,shah_eos,hvg), S(our_neutro,shah_eos,hvg)),
  n_genes_all = length(common), n_genes_hvg = length(hvg))
write.csv(res, file.path(out, "shah_specificity_control.csv"), row.names=FALSE)
cat("=== Specificity: our pseudobulk vs Shah bulk EOSINOPHILS ===\n"); print(res)
cat("\n(Eosinophils should correlate higher than macrophages/neutrophils, especially on variable genes.)\n")
# gene-level for eos vs mac scatter on HVG
gl <- data.frame(gene=hvg, shah_eos=cpm(shah_eos[hvg]), our_eos=cpm(our_eos[hvg]), our_mac=cpm(our_mac[hvg]))
write.csv(gl, file.path(out, "shah_specificity_genelevel_hvg.csv"), row.names=FALSE)
cat("Done 58b.\n")
