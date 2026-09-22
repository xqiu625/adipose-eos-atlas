#!/usr/bin/env Rscript
# Package-level smoke tests for the pinned R environment (see
# setup_pinned_r_env.R). Loads each version-sensitive package and runs a
# small synthetic-data call to catch API breakage before the R branch of the
# audited rerun is released to run actual pipeline scripts. Per investigator
# requirement: sccomp and scSHC get real FUNCTIONAL calls through their core
# entrypoints (not just library()/exists() checks), since successful install
# does not verify API compatibility for the two packages with the most
# consequential version-sensitive behavior in this rerun.

rerun_root <- Sys.getenv("RERUN_ROOT")
if (rerun_root == "") stop("RERUN_ROOT must be set")
lib_dir <- file.path(rerun_root, "r_env", "renv_library")
.libPaths(c(lib_dir, .libPaths()))

results <- list()
smoke <- function(name, expr) {
  warned <- character(0)
  ok <- withCallingHandlers(
    tryCatch({ force(expr); TRUE }, error = function(e) {
      message("SMOKE TEST FAILED: ", name, " -- ", conditionMessage(e))
      FALSE
    }),
    warning = function(w) {
      warned <<- c(warned, conditionMessage(w))
      invokeRestart("muffleWarning")
    }
  )
  results[[name]] <<- list(ok = ok, warnings = warned)
  cat(name, ":", if (ok) "PASS" else "FAIL",
      if (length(warned) > 0) paste0(" (", length(warned), " warning(s))") else "", "\n")
  if (length(warned) > 0) for (w in warned) cat("    warning:", w, "\n")
}

set.seed(1)

# ---- Load-only tests for the remaining required packages ------------------
smoke("clusterProfiler (load)", { library(clusterProfiler) })
smoke("celda (load)", { library(celda) })
smoke("ComplexHeatmap (load)", { library(ComplexHeatmap) })
smoke("SeuratWrappers (load)", { library(SeuratWrappers) })
smoke("clustree (load)", { library(clustree) })

# ---- Seurat: build a tiny object, normalize, cluster -----------------------
smoke("Seurat", {
  library(Seurat)
  m <- matrix(rpois(200 * 20, 5), nrow = 200, ncol = 20)
  rownames(m) <- paste0("gene", 1:200); colnames(m) <- paste0("cell", 1:20)
  so <- CreateSeuratObject(counts = m)
  so <- NormalizeData(so, verbose = FALSE)
  so <- FindVariableFeatures(so, verbose = FALSE)
  stopifnot(ncol(so) == 20)
  stopifnot(as.character(packageVersion("Seurat")) == "5.1.0")
})

# ---- harmony: run on tiny synthetic PCA-like matrix with 2 batches --------
smoke("harmony", {
  library(harmony)
  pcs <- matrix(rnorm(40 * 5), nrow = 40, ncol = 5)
  meta <- data.frame(batch = rep(c("a", "b"), each = 20))
  rownames(pcs) <- rownames(meta) <- paste0("cell", 1:40)
  colnames(pcs) <- paste0("PC", 1:5)
  h <- RunHarmony(pcs, meta, vars_use = "batch", verbose = FALSE)
  stopifnot(nrow(h) == 40)
  stopifnot(as.character(packageVersion("harmony")) == "1.2.1")
})

# ---- CellChat: load + version check (full analysis too heavy for a smoke
# test; createCellChat existing + version match is the meaningful check here) --
smoke("CellChat", {
  library(CellChat)
  stopifnot(exists("createCellChat"))
  stopifnot(as.character(packageVersion("CellChat")) == "1.6.1")
})

# ---- monocle3: build a tiny cell_data_set ----------------------------------
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

# ---- DESeq2: tiny synthetic count matrix through DESeqDataSetFromMatrix ---
smoke("DESeq2", {
  suppressMessages(library(DESeq2))
  cts <- matrix(rpois(6 * 20, 50), nrow = 20, ncol = 6)
  colnames(cts) <- paste0("s", 1:6)
  coldata <- data.frame(condition = rep(c("a", "b"), each = 3), row.names = colnames(cts))
  dds <- DESeqDataSetFromMatrix(countData = cts, colData = coldata, design = ~condition)
  stopifnot(ncol(dds) == 6)
})

# ---- sccomp: REAL functional call through sccomp_estimate() on tiny
# synthetic compositional data (sccomp 1.x long-format API) --------------------
smoke("sccomp FUNCTIONAL (sccomp_estimate)", {
  suppressPackageStartupMessages(library(sccomp))
  stopifnot(as.character(packageVersion("sccomp")) == "1.7.16")
  set.seed(42)
  n_samples <- 6
  cell_groups <- paste0("cg", 1:3)
  synth <- do.call(rbind, lapply(1:n_samples, function(s) {
    counts <- rpois(3, lambda = c(50, 30, 20) + ifelse(s <= 3, c(10, 0, -5), 0))
    counts <- as.integer(pmax(counts, 1))
    data.frame(
      sample = paste0("sample", s),
      cell_group = cell_groups,
      count = counts,
      condition = if (s <= 3) "A" else "B",
      stringsAsFactors = FALSE
    )
  }))
  synth$sample <- factor(synth$sample)
  synth$cell_group <- factor(synth$cell_group)
  synth$condition <- factor(synth$condition)

  fit <- sccomp_estimate(
    synth,
    formula_composition = ~condition,
    .sample = sample,
    .cell_group = cell_group,
    .count = count,
    cores = 1,
    variational_inference = TRUE
  )
  stopifnot(!is.null(fit))
  cat("  sccomp_estimate() returned an object of class:", paste(class(fit), collapse = ", "), "\n")
})

# ---- scSHC: REAL functional call through scSHC() on tiny synthetic counts --
smoke("scSHC FUNCTIONAL (scSHC())", {
  library(scSHC)
  stopifnot(exists("scSHC"))
  set.seed(7)
  n_genes <- 40; n_cells <- 30
  m <- matrix(rpois(n_genes * n_cells, lambda = 3), nrow = n_genes, ncol = n_cells)
  # inject a clear 2-cluster structure so the test isn't degenerate
  m[1:10, 1:15] <- m[1:10, 1:15] + 15
  rownames(m) <- paste0("gene", 1:n_genes)
  colnames(m) <- paste0("cell", 1:n_cells)
  res <- scSHC(m, num_features = n_genes, num_PCs = 5, parallel = FALSE)
  stopifnot(!is.null(res))
  cat("  scSHC() returned an object of class:", paste(class(res), collapse = ", "), "\n")
})

cat("\n=== Smoke test summary ===\n")
n_pass <- sum(vapply(results, function(r) r$ok, logical(1)))
n_total <- length(results)
cat(n_pass, "/", n_total, "passed\n")
for (nm in names(results)) cat("-", nm, ":", if (results[[nm]]$ok) "PASS" else "FAIL", "\n")

if (n_pass < n_total) quit(status = 1) else quit(status = 0)
