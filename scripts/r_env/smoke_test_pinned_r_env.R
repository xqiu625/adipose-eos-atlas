#!/usr/bin/env Rscript
# Package-level smoke tests for the pinned R environment (see
# setup_pinned_r_env.R). Loads each version-sensitive package and runs one
# trivial synthetic-data call to catch API breakage before the R branch of
# the audited rerun is released to run actual pipeline scripts.

rerun_root <- Sys.getenv("RERUN_ROOT")
if (rerun_root == "") stop("RERUN_ROOT must be set")
lib_dir <- file.path(rerun_root, "r_env", "renv_library")
.libPaths(c(lib_dir, .libPaths()))

results <- list()
smoke <- function(name, expr) {
  ok <- tryCatch({ force(expr); TRUE }, error = function(e) {
    message("SMOKE TEST FAILED: ", name, " -- ", conditionMessage(e))
    FALSE
  })
  results[[name]] <<- ok
  cat(name, ":", if (ok) "PASS" else "FAIL", "\n")
}

set.seed(1)

# Seurat: build a tiny object, normalize, cluster
smoke("Seurat", {
  library(Seurat)
  m <- matrix(rpois(200 * 20, 5), nrow = 200, ncol = 20)
  rownames(m) <- paste0("gene", 1:200); colnames(m) <- paste0("cell", 1:20)
  so <- CreateSeuratObject(counts = m)
  so <- NormalizeData(so, verbose = FALSE)
  so <- FindVariableFeatures(so, verbose = FALSE)
  stopifnot(ncol(so) == 20)
})

# harmony: run on tiny synthetic PCA-like matrix with 2 batches
smoke("harmony", {
  library(harmony)
  pcs <- matrix(rnorm(40 * 5), nrow = 40, ncol = 5)
  meta <- data.frame(batch = rep(c("a", "b"), each = 20))
  rownames(pcs) <- rownames(meta) <- paste0("cell", 1:40)
  colnames(pcs) <- paste0("PC", 1:5)
  h <- RunHarmony(pcs, meta, vars_use = "batch", verbose = FALSE)
  stopifnot(nrow(h) == 40)
})

# CellChat: just confirm the package loads and createCellChat exists
smoke("CellChat", {
  library(CellChat)
  stopifnot(exists("createCellChat"))
  stopifnot(as.character(packageVersion("CellChat")) == "1.6.1")
})

# monocle3: build a tiny cell_data_set
smoke("monocle3", {
  library(monocle3)
  m <- matrix(rpois(50 * 10, 5), nrow = 50, ncol = 10)
  rownames(m) <- paste0("gene", 1:50); colnames(m) <- paste0("cell", 1:10)
  gene_meta <- data.frame(gene_short_name = rownames(m), row.names = rownames(m))
  cell_meta <- data.frame(row.names = colnames(m))
  cds <- new_cell_data_set(m, cell_metadata = cell_meta, gene_metadata = gene_meta)
  stopifnot(ncol(cds) == 10)
  stopifnot(as.character(packageVersion("monocle3")) == "1.2.9")
})

# sccomp: confirm it loads and its main entrypoint exists (full model fit is
# too slow for a smoke test)
smoke("sccomp", {
  library(sccomp)
  stopifnot(exists("sccomp_estimate"))
  cat("  sccomp version actually installed:", as.character(packageVersion("sccomp")), "\n")
})

# scSHC: confirm it loads and its main entrypoint exists
smoke("scSHC", {
  library(scSHC)
  stopifnot(exists("scSHC"))
})

# DESeq2: tiny synthetic count matrix through DESeqDataSetFromMatrix
smoke("DESeq2", {
  suppressMessages(library(DESeq2))
  cts <- matrix(rpois(6 * 20, 50), nrow = 20, ncol = 6)
  colnames(cts) <- paste0("s", 1:6)
  coldata <- data.frame(condition = rep(c("a", "b"), each = 3), row.names = colnames(cts))
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~condition)
  stopifnot(ncol(dds) == 6)
})

cat("\n=== Smoke test summary ===\n")
n_pass <- sum(unlist(results))
n_total <- length(results)
cat(n_pass, "/", n_total, "passed\n")
for (nm in names(results)) cat("-", nm, ":", if (results[[nm]]) "PASS" else "FAIL", "\n")

if (n_pass < n_total) quit(status = 1) else quit(status = 0)
