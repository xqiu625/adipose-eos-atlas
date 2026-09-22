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

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

BASE    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out_dir <- file.path(BASE, "scenic_per_cluster")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("14c_export_per_cluster_loom.R\n")
cat("Export each Eos cluster as mtx+metadata for per-cluster pySCENIC\n\n")

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
  n_cells <- ncol(sub)
  cat(sprintf("  %d cells\n", n_cells))

  cl_dir <- file.path(out_dir, cl)
  dir.create(cl_dir, recursive = TRUE, showWarnings = FALSE)

  # Raw counts matrix (genes x cells)
  raw <- GetAssayData(sub, layer = "counts")

  # Write sparse matrix (genes x cells)
  writeMM(raw, file.path(cl_dir, "matrix.mtx"))

  # Write gene names and cell barcodes
  write.table(rownames(raw), file.path(cl_dir, "features.tsv"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)
  write.table(colnames(raw), file.path(cl_dir, "barcodes.tsv"),
              quote = FALSE, row.names = FALSE, col.names = FALSE)

  # Write cell metadata
  meta <- data.frame(
    barcode   = colnames(raw),
    cluster   = rep(cl, n_cells),
    condition = as.character(sub$condition),
    sex       = as.character(sub$sex),
    stringsAsFactors = FALSE
  )
  write.csv(meta, file.path(cl_dir, "metadata.csv"), row.names = FALSE)

  cat(sprintf("  Written: matrix.mtx (%d genes x %d cells), features.tsv, barcodes.tsv, metadata.csv\n",
              nrow(raw), n_cells))
}

cat("\n14c complete.\n")
cat(sprintf("Output: %s\n", out_dir))
cat("Next: run 14d_create_looms.py to convert to loom format\n")
