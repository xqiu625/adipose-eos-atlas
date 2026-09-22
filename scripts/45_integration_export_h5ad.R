#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --job-name=int_export
#SBATCH --output=int_export_%j.out

# ==============================================================================
# 45_integration_export_h5ad.R
# Export the eosinophil object to h5ad as the input for the integration-method
# benchmark (46_integration_benchmark.py).
#
# Writes raw counts in X plus the cell metadata needed to score integration:
#   - cluster       : manuscript C0-C5 (bio-conservation label_key)
#   - sample        : orig.ident, the per-library batch
#   - gem_reaction  : GEM1-GEM4, the 10x Flex multiplex reaction (technical batch)
#   - condition / sex / genotype : biological covariates (held out, for plots)
#
# The benchmark uses `gem_reaction` (or `sample`) as the batch_key and `cluster`
# as the label_key. Counts are exported raw so scVI/scANVI can model them and so
# every method re-derives HVGs/PCA identically.
#
# Input:  $BASE/eos_subclustered.rds (C0-C6; C6 dropped here)
# Output: $BASE/integration/benchmark/eos_for_integration.h5ad
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(dplyr)
  library(reticulate)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

use_condaenv("/bigdata/godziklab/shared/Xinru/cellrank_env")
anndata <- import("anndata")
np      <- import("numpy")
scipy   <- import("scipy.sparse")
pd      <- import("pandas")

BASE    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out_dir <- file.path(BASE, "integration", "benchmark")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# ── Load eosinophil object ───────────────────────────────────────────────────
cat("Loading eos_subclustered.rds...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
if ("C6" %in% levels(Idents(eos))) {
  eos <- subset(eos, idents = "C6", invert = TRUE)
  Idents(eos)      <- droplevels(Idents(eos))
  eos$eos_cluster  <- droplevels(eos$eos_cluster)
  cat(sprintf("  Removed C6. Remaining: %d eos (C0-C5)\n", ncol(eos)))
}
if (inherits(eos[["RNA"]], "Assay5")) eos[["RNA"]] <- JoinLayers(eos[["RNA"]])
cat(sprintf("  %d genes x %d cells\n", nrow(eos), ncol(eos)))

# ── Derive GEM reaction (technical batch), mirrors 19_batch_effect_evidence.R ──
cat("Assigning GEM reactions...\n")
sm <- eos@meta.data %>%
  select(orig.ident, condition) %>%
  distinct() %>%
  mutate(
    replicate_num = as.integer(sub(".*_(\\d+)$", "\\1", orig.ident)),
    is_infected   = grepl("Infected", condition),
    genotype_g    = sub("^(WT|KO)_.*", "\\1", condition),
    gem_reaction  = case_when(
      is_infected                                                ~ "GEM2",
      !is_infected & replicate_num == 1                          ~ "GEM1",
      !is_infected & genotype_g == "WT" & replicate_num %in% 2:3 ~ "GEM3",
      !is_infected & genotype_g == "KO" & replicate_num %in% 2:3 ~ "GEM4",
      TRUE                                                       ~ "Unknown"
    )
  )
gem_map          <- setNames(sm$gem_reaction, sm$orig.ident)
# unname(): gem_map[...] is named by sample ID; Seurat's `$<-` would try to match
# those names to cell barcodes and error ("No cell overlap"). Strip names so the
# vector is assigned positionally, in cell order.
eos$gem_reaction <- unname(gem_map[eos$orig.ident])
print(table(eos$orig.ident, eos$gem_reaction))

# ── Raw counts (cells x genes) -> scipy CSC ───────────────────────────────────
cat("Extracting raw counts...\n")
raw_t <- t(GetAssayData(eos, layer = "counts"))
to_scipy_csc <- function(mat_t) {
  scipy$csc_matrix(
    tuple(np$array(mat_t@x),
          np$array(as.integer(mat_t@i)),
          np$array(as.integer(mat_t@p))),
    shape = tuple(as.integer(nrow(mat_t)), as.integer(ncol(mat_t)))
  )
}
raw_scipy <- to_scipy_csc(raw_t)
cat(sprintf("  Raw counts (cells x genes): %d x %d\n", nrow(raw_t), ncol(raw_t)))

# ── Cell metadata ─────────────────────────────────────────────────────────────
wanted    <- c("eos_cluster", "orig.ident", "gem_reaction", "condition", "sex", "genotype")
available <- wanted[wanted %in% colnames(eos@meta.data)]
obs_df    <- eos@meta.data[, available, drop = FALSE]
names(obs_df)[names(obs_df) == "eos_cluster"] <- "cluster"
names(obs_df)[names(obs_df) == "orig.ident"]  <- "sample"
# ensure character (categoricals carried as strings into AnnData)
obs_df[] <- lapply(obs_df, as.character)
cat(sprintf("  obs columns: %s\n", paste(colnames(obs_df), collapse = ", ")))

obs_pd       <- pd$DataFrame(obs_df)
obs_pd$index <- colnames(eos)
var_pd       <- pd$DataFrame(index = rownames(eos))

adata <- anndata$AnnData(X = raw_scipy, obs = obs_pd, var = var_pd)

h5ad_path <- file.path(out_dir, "eos_for_integration.h5ad")
cat(sprintf("Saving h5ad: %s\n", h5ad_path))
adata$write_h5ad(h5ad_path)

cat("\nCluster distribution:\n");      print(table(obs_df$cluster))
cat("\nGEM reaction distribution:\n"); print(table(obs_df$gem_reaction))
cat("\nDone! 45_integration_export_h5ad.R completed successfully.\n")
