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
# Probe reference + our objects for the broad-reference (Tier 1) build.
suppressPackageStartupMessages({ library(Seurat) })

GD  <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature"
B   <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"

probe <- function(path, label_cols = NULL, tag) {
  if (!file.exists(path)) { cat("\n###", tag, "MISSING:", path, "\n"); return(invisible()) }
  o <- readRDS(path)
  cat("\n############", tag, "############\n")
  cat("cells:", ncol(o), " genes:", nrow(o), "\n")
  cat("assays:", paste(Assays(o), collapse=","), "\n")
  cat("meta cols:", paste(colnames(o@meta.data), collapse=", "), "\n")
  cat("idents:", paste(head(levels(Idents(o)), 30), collapse=", "), "\n")
  # auto-detect likely label columns
  cand <- grep("type|anno|cluster|subset|ident|tissue|organ|celltype|stage|label|state",
               colnames(o@meta.data), value = TRUE, ignore.case = TRUE)
  for (c in unique(c(cand, label_cols))) {
    v <- o@meta.data[[c]]
    if (!is.null(v) && length(unique(v)) <= 40) {
      cat(sprintf("  [%s]: ", c)); print(table(v))
    }
  }
  cat("example barcodes:", paste(head(colnames(o), 3), collapse=" | "), "\n")
  rm(o); gc()
}

probe(file.path(GD, "Gurtner_eosinophil_subset.rds"), tag = "GURTNER_EOS_SUBSET")
probe(file.path(GD, "Gurtner_merged_tissues_seurat_object.rds"), tag = "GURTNER_MERGED_TISSUES")
probe(file.path(B, "all_cells_annotated.rds"), label_cols = c("celltype","broad_celltype"), tag = "OUR_ALL_CELLS")
probe(file.path(B, "eos_subclustered.rds"), label_cols = c("eos_cluster"), tag = "OUR_EOS")
cat("\nDone probe.\n")
