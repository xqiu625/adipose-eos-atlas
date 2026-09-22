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
# 53_c0_coordination.R   (Reviewer #1 — MHC-II is an intrinsic program in C0,
# not ambient macrophage RNA). Within C0 cells, MHC-II genes should co-vary with
# each other (coordinated program) but NOT with the ambient macrophage module
# (C1q/Apoe/Lyz2). Low MHC-II<->mac cross-correlation = MHC-II is intrinsic.
# Input: $BASE/eos_subclustered.rds ; Output: $BASE/revision/figures/cluster_robustness/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out  <- file.path(BASE, "revision", "figures", "cluster_robustness")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
if ("C6" %in% levels(Idents(eos))) { eos <- subset(eos, idents="C6", invert=TRUE); Idents(eos) <- droplevels(Idents(eos)) }
if (inherits(eos[["RNA"]], "Assay5")) eos[["RNA"]] <- JoinLayers(eos[["RNA"]])

mhc <- intersect(c("H2-Aa","H2-Ab1","H2-Eb1","Cd74"), rownames(eos))
mac <- intersect(c("C1qa","C1qb","Apoe","Lyz2","Folr2"), rownames(eos))
genes <- c(mhc, mac)

c0 <- WhichCells(eos, idents = "C0")
expr <- as.matrix(GetAssayData(eos, layer = "data")[genes, c0])
cm <- cor(t(expr), method = "spearman")
write.csv(round(cm, 3), file.path(out, "c0_coordination_cormat.csv"))
cat("=== C0 within-cluster Spearman correlation (MHC-II + mac ambient genes) ===\n"); print(round(cm, 2))

offblock <- function(M, a, b) { v <- M[a, b]; mean(v[is.finite(v)]) }
withinblock <- function(M, a) { sub <- M[a, a]; mean(sub[upper.tri(sub)]) }
summ <- data.frame(
  block = c("within MHC-II", "within mac-ambient", "MHC-II x mac (cross)"),
  mean_spearman = round(c(withinblock(cm, mhc), withinblock(cm, mac), offblock(cm, mhc, mac)), 3)
)
write.csv(summ, file.path(out, "c0_coordination_summary.csv"), row.names = FALSE)
cat("\n=== block-correlation summary ===\n"); print(summ)
cat("\nIntrinsic MHC-II program if: within-MHC-II HIGH, MHC-II x mac (cross) LOW.\n")
cat("Done 53_c0_coordination.R\n")
