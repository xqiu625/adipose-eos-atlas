#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 70_sccomp_both.R   (composition: un-integrated C0-C5 vs integrated i0-i5)
# Robust sccomp (single-core; drop list-columns before writing). Runs the three
# manuscript contrasts on BOTH the decont+>0.10 un-integrated object (C0-C5) and
# the ComBat-integrated object (i0-i5), so we can compare the sex/RELMa/infection
# compositional shifts side by side.
# Output: $REDO/sccomp/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(sccomp) })
options(lifecycle_verbosity = "quiet")
REDO <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/revision/redo"
out  <- file.path(REDO, "sccomp"); dir.create(out, recursive=TRUE, showWarnings=FALSE)

run_sccomp <- function(md, formula, tag) {
  sub <- md[, c("sample","cluster", all.vars(formula))]
  est <- tryCatch(sccomp_estimate(sub, formula_composition=formula, .sample=sample,
                                  .cell_group=cluster, cores=1, verbose=FALSE),
                  error=function(e){cat("  [",tag,"] est error:",conditionMessage(e),"\n");NULL})
  if (is.null(est)) return(NULL)
  res <- tryCatch(sccomp_test(est), error=function(e){cat("  test err:",conditionMessage(e),"\n"); est})
  res <- as.data.frame(res)
  res <- res[, vapply(res, function(c) !is.list(c), logical(1)), drop=FALSE]   # drop list-cols
  write.csv(res, file.path(out, paste0(tag, "_full.csv")), row.names=FALSE)
  pa <- intersect(c("parameter","factor"), colnames(res))[1]
  cg <- intersect(c("cell_group","cluster"), colnames(res))[1]
  ef <- intersect(c("c_effect","effect"), colnames(res))[1]
  fd <- intersect(c("c_FDR","FDR"), colnames(res))[1]
  k <- res[!grepl("ntercept", res[[pa]]), c(cg,pa,ef,fd)]
  colnames(k) <- c("cluster","parameter","effect","FDR"); k$sig <- ifelse(k$FDR<0.1,"*","")
  write.csv(k, file.path(out, paste0(tag, ".csv")), row.names=FALSE)
  cat(sprintf("=== %s ===\n", tag)); print(k[order(k$cluster), ]); k
}

for (obj in c("unintegrated","combat")) {
  rds <- if (obj=="unintegrated") "eos_decont010.rds" else "eos_combat010.rds"
  cat("\n############", obj, "(", rds, ") ############\n")
  eos <- readRDS(file.path(REDO, rds))
  cl  <- if (obj=="unintegrated") as.character(Idents(eos)) else as.character(eos$new_cluster)
  md  <- data.frame(sample=as.character(eos$orig.ident), cluster=cl,
                    sex=as.character(eos$sex), genotype=as.character(eos$genotype),
                    condition=as.character(eos$condition), stringsAsFactors=FALSE)
  md$infected <- grepl("Inf", md$condition)
  run_sccomp(md[md$genotype=="WT" & !md$infected, ], ~ sex,      paste0(obj,"_sex_WT"))
  run_sccomp(md[md$sex=="Female" & !md$infected, ],  ~ genotype, paste0(obj,"_genotype_F"))
  run_sccomp(md[md$sex=="Female", ],                 ~ infected, paste0(obj,"_infection_F"))
  rm(eos); gc()
}
cat("\nDone 70.\n")
