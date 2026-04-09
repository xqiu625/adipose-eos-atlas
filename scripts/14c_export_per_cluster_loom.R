#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(reticulate)
})

use_condaenv("/bigdata/godziklab/shared/Xinru/cellrank_env")
loompy <- import("loompy")
np     <- import("numpy")
scipy  <- import("scipy.sparse")

BASE    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out_dir <- file.path(BASE, "scenic_per_cluster")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("14c_export_per_cluster_loom.R\n")
cat("Export each Eos cluster as loom for per-cluster pySCENIC\n\n")

cat("Loading eos_subclustered.rds...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
if ("C6" %in% levels(Idents(eos))) {
  eos <- subset(eos, idents = "C6", invert = TRUE)
  Idents(eos) <- droplevels(Idents(eos))
}
cat(sprintf("  %d cells, clusters: %s\n", ncol(eos),
            paste(levels(Idents(eos)), collapse = ", ")))

clusters <- levels(Idents(eos))

for (cl in clusters) {
  cat(sprintf("\nExporting %s...\n", cl))
  sub <- subset(eos, idents = cl)
  cat(sprintf("  %d cells\n", ncol(sub)))

  cl_dir <- file.path(out_dir, cl)
  dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)
  loom_path <- file.path(cl_dir, "input.loom")

  # Raw counts (genes x cells)
  raw <- GetAssayData(sub, layer = "counts")

  # Convert to scipy sparse (CSC, cells x genes for loom = genes x cells)
  raw_csc <- scipy$csc_matrix(
    tuple(np$array(raw@x),
          np$array(as.integer(raw@i)),
          np$array(as.integer(raw@p))),
    shape = tuple(as.integer(nrow(raw)), as.integer(ncol(raw)))
  )

  # Build loom via loompy.create
  row_attrs <- list(Gene = rownames(raw))
  col_attrs <- list(
    CellID    = colnames(raw),
    cluster   = rep(cl, ncol(raw)),
    condition = as.character(sub$condition),
    sex       = as.character(sub$sex)
  )

  loompy$create(
    loom_path,
    raw_csc,
    row_attrs = reticulate::r_to_py(row_attrs),
    col_attrs = reticulate::r_to_py(col_attrs)
  )
  cat(sprintf("  Saved: %s\n", loom_path))
}

cat("\n14c complete.\n")
cat(sprintf("Output: %s\n", out_dir))
