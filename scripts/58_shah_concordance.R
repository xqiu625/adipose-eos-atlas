#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 58_shah_concordance.R   (Reviewer #1 — fixation/representativeness, bulk substitute)
# Compare the pseudobulk of our scRNA eosinophils against the PUBLISHED bulk
# RNA-seq of sorted eosinophils from Shah, Knights et al. 2023 (GSE217551:
# gonadal/visceral adipose, subcutaneous adipose, blood). High concordance with
# the ADIPOSE bulk (and higher than with blood) shows our FLEX/fixed transcriptome
# is representative of bona fide adipose eosinophils — addressing the reviewer's
# concern that fixation/processing distorted the population, without new bulk.
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

EXT <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/external/shah2023"
BASE<- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out <- file.path(BASE, "revision", "figures", "shah_concordance"); dir.create(out, recursive=TRUE, showWarnings=FALSE)

# ── tx2gene ───────────────────────────────────────────────────────────────────
t2g <- read.table(file.path(EXT, "tx2gene.tsv"), sep="\t", stringsAsFactors=FALSE, quote="")
tx2g <- setNames(t2g$V2, t2g$V1)

# ── Shah salmon .sf -> gene-level counts ──────────────────────────────────────
sf <- list.files(file.path(EXT, "raw"), pattern="\\.sf\\.gz$", full.names=TRUE)
grp <- ifelse(grepl("scEos", sf), "subcutaneous", ifelse(grepl("gEos", sf), "gonadal", "blood"))
read_sf <- function(f) {
  d <- read.table(gzfile(f), header=TRUE, sep="\t", stringsAsFactors=FALSE)
  d$gene <- tx2g[sub("\\.\\d+$", "", d$Name)]
  d <- d[!is.na(d$gene), ]
  tapply(d$NumReads, d$gene, sum)
}
cat("Reading", length(sf), "Shah samples...\n")
shah_list <- lapply(sf, read_sf)
genes_sh <- Reduce(union, lapply(shah_list, names))
shah <- sapply(shah_list, function(v) { x <- setNames(rep(0, length(genes_sh)), genes_sh); x[names(v)] <- v; x })
colnames(shah) <- paste0(grp, "_", seq_along(sf))
cat("Shah gene x sample:", dim(shah)[1], "x", dim(shah)[2], " groups:", paste(table(grp), names(table(grp))), "\n")

# ── Our scRNA eosinophil pseudobulk ───────────────────────────────────────────
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
if ("C6" %in% levels(Idents(eos))) { eos <- subset(eos, idents="C6", invert=TRUE); Idents(eos)<-droplevels(Idents(eos)) }
if (inherits(eos[["RNA"]], "Assay5")) eos[["RNA"]] <- JoinLayers(eos[["RNA"]])
cnt <- GetAssayData(eos, layer="counts")
our_all <- Matrix::rowSums(cnt)
our_cl  <- sapply(paste0("C",0:5), function(c) Matrix::rowSums(cnt[, WhichCells(eos, idents=c), drop=FALSE]))

# ── CPM + log; common genes ───────────────────────────────────────────────────
cpm <- function(v) log1p(v / sum(v) * 1e6)
shah_g <- list(gonadal = rowMeans(shah[, grp=="gonadal"]),
               subcutaneous = rowMeans(shah[, grp=="subcutaneous"]),
               blood = rowMeans(shah[, grp=="blood"]))
common <- Reduce(intersect, c(list(names(our_all)), lapply(shah_g, names)))
# keep genes expressed in our eos
common <- common[our_all[common] > 0]
cat("Common expressed genes:", length(common), "\n")

S <- function(a, b) round(cor(cpm(a[common]), cpm(b[common]), method="spearman"), 3)
rows <- list()
rows[["Ours_all_eos"]] <- sapply(shah_g, function(g) S(our_all, g))
for (c in paste0("C",0:5)) rows[[c]] <- sapply(shah_g, function(g) S(our_cl[,c], g))
res <- do.call(rbind, rows)
write.csv(res, file.path(out, "shah_concordance_spearman.csv"))
cat("\n=== Spearman: our eos pseudobulk vs Shah bulk (rows=ours, cols=Shah tissue) ===\n"); print(res)

# ── gene-level table for figure (our all-eos vs Shah gonadal & blood) ─────────
gl <- data.frame(gene = common,
                 ours = cpm(our_all[common]),
                 shah_gonadal = cpm(shah_g$gonadal[common]),
                 shah_blood = cpm(shah_g$blood[common]))
write.csv(gl, file.path(out, "shah_concordance_genelevel.csv"), row.names=FALSE)
cat("\nDone 58_shah_concordance.R. Outputs in", out, "\n")
