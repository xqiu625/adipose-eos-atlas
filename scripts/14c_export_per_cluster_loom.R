#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(reticulate)
})

use_condaenv("/bigdata/godziklab/shared/Xinru/cellrank_env")
ad    <- import("anndata")
np    <- import("numpy")
scipy <- import("scipy.sparse")

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

  # Raw counts (genes x cells) -> transpose to cells x genes for AnnData
  raw <- GetAssayData(sub, layer = "counts")

  # Convert to scipy CSR (cells x genes)
  raw_csr <- scipy$csr_matrix(
    tuple(np$array(raw@x),
          np$array(as.integer(raw@i)),
          np$array(as.integer(raw@p))),
    shape = tuple(as.integer(nrow(raw)), as.integer(ncol(raw)))
  )$T   # transpose: cells x genes

  # Build AnnData and write as loom (avoids loompy np.string_ bug)
  obs_df <- reticulate::r_to_py(
    data.frame(
      row.names = colnames(raw),
      cluster   = rep(cl, ncol(raw)),
      condition = as.character(sub$condition),
      sex       = as.character(sub$sex)
    )
  )
  var_df <- reticulate::r_to_py(
    data.frame(row.names = rownames(raw))
  )

  adata <- ad$AnnData(X = raw_csr, obs = obs_df, var = var_df)
  adata$obs_names <- reticulate::r_to_py(as.list(colnames(raw)))
  adata$var_names <- reticulate::r_to_py(as.list(rownames(raw)))
  adata$write_loom(loom_path, write_obsm_varm = FALSE)
  cat(sprintf("  Saved: %s\n", loom_path))
}

cat("\n14c complete.\n")
cat(sprintf("Output: %s\n", out_dir))
