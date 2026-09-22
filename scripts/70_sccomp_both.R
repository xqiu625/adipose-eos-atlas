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
  if (!nzchar(Sys.getenv("CMDSTAN"))) {
    stop("RERUN_R_LIB audit guard: this script loads sccomp, which requires CMDSTAN to be set when running under the audited rerun.")
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
