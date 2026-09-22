#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 68_sccomp_decont.R   (re-run composition on decont + >0.10 object)
# Bayesian compositional analysis (sccomp) on the decontaminated, contamination
# ->0.10-filtered eosinophils with the published C0-C5 labels. Confirms the
# sex / RELMa / infection cluster-proportion shifts (Fig 2D) hold after the
# contamination removal. Robust output extraction (full tibble saved).
# Input : $REDO/eos_decont010.rds ; Output: $REDO/sccomp/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(sccomp); library(tidyr) })

REDO <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/revision/redo"
out  <- file.path(REDO, "sccomp"); dir.create(out, recursive=TRUE, showWarnings=FALSE)
eos <- readRDS(file.path(REDO, "eos_decont010.rds"))
md <- eos@meta.data
md$cluster <- as.character(Idents(eos))
md$sample  <- as.character(md$orig.ident)
md$infected <- grepl("Inf", md$condition)

run_sccomp <- function(sub, formula, label) {
  sub <- sub[, c("sample","cluster", all.vars(formula))]
  est <- tryCatch(
    sccomp_estimate(sub, formula_composition = formula,
                    .sample = sample, .cell_group = cluster,
                    cores = 4, verbose = FALSE),
    error = function(e){ cat("  [",label,"] estimate error:", conditionMessage(e), "\n"); NULL })
  if (is.null(est)) return(invisible())
  res <- tryCatch(sccomp_test(est),
                  error = function(e){ cat("  [",label,"] test error:", conditionMessage(e), "\n"); est })
  df <- as.data.frame(res)
  write.csv(df, file.path(out, paste0("sccomp_", label, "_full.csv")), row.names = FALSE)
  # extract the effect of interest (non-intercept parameter)
  eff_col <- intersect(c("c_effect","effect"), colnames(df))[1]
  fdr_col <- intersect(c("c_FDR","FDR"), colnames(df))[1]
  cg_col  <- intersect(c("cell_group","cluster"), colnames(df))[1]
  par_col <- intersect(c("parameter","factor"), colnames(df))[1]
  if (!is.na(eff_col)) {
    keep <- df[!grepl("ntercept", df[[par_col]]), c(cg_col, par_col, eff_col, fdr_col)]
    colnames(keep) <- c("cluster","parameter","effect","FDR")
    keep$contrast <- label
    write.csv(keep, file.path(out, paste0("sccomp_", label, ".csv")), row.names = FALSE)
    cat(sprintf("\n=== [%s] composition shifts (FDR<0.1 flagged) ===\n", label))
    keep$sig <- ifelse(keep$FDR < 0.1, "*", "")
    print(keep[order(keep$cluster), ])
  }
}

cat("Contrast 1: sex (WT female vs WT male, HFD only)\n")
run_sccomp(md[md$genotype=="WT" & !md$infected, ], ~ sex, "sex_WT")
cat("\nContrast 2: genotype (KO vs WT, females HFD only)\n")
run_sccomp(md[md$sex=="Female" & !md$infected, ], ~ genotype, "genotype_F")
cat("\nContrast 3: infection (infected vs naive, all females)\n")
run_sccomp(md[md$sex=="Female", ], ~ infected, "infection_F")

cat("\nDone 68.\n")
