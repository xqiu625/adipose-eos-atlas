#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=1:00:00
#SBATCH --job-name=scenic_export
#SBATCH --output=scenic_export_%j.out

# ==============================================================================
# 14_scenic_export.R
# Export eosinophil Seurat object to h5ad format for pySCENIC
# Transfers raw counts, normalized counts, metadata, and UMAP embeddings
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(reticulate)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

use_condaenv("/bigdata/godziklab/shared/Xinru/cellrank_env")
anndata <- import("anndata")
np      <- import("numpy")
scipy   <- import("scipy.sparse")

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"

# -- Output directory ----------------------------------------------------------
out_dir <- file.path(BASE, "scenic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# -- Load eosinophil object ----------------------------------------------------
cat("========================================\n")
cat("14_scenic_export.R\n")
cat("Export eosinophil data for pySCENIC\n")
cat("========================================\n\n")

cat("Loading eos_subclustered.rds...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
# Remove C6 if still present (too few cells, n=14)
if ("C6" %in% levels(Idents(eos))) {
  eos <- subset(eos, idents = "C6", invert = TRUE)
  Idents(eos) <- droplevels(Idents(eos))
  cat(sprintf("  Removed C6. Remaining: %d eos (C0-C5)\n", ncol(eos)))
}
cat(sprintf("  %d genes x %d cells\n", nrow(eos), ncol(eos)))

# -- Extract raw counts --------------------------------------------------------
cat("Extracting raw counts...\n")
raw_counts <- GetAssayData(eos, layer = "counts")

# Transpose: AnnData expects cells x genes
raw_t <- t(raw_counts)

cat(sprintf("  Raw counts matrix (cells x genes): %d x %d\n",
            nrow(raw_t), ncol(raw_t)))

# -- Extract normalized counts -------------------------------------------------
cat("Extracting normalized counts...\n")
norm_counts <- GetAssayData(eos, layer = "data")
norm_t <- t(norm_counts)

# -- Convert to scipy sparse (CSC) --------------------------------------------
cat("Converting matrices to scipy sparse format...\n")

to_scipy_csc <- function(mat_t) {
  # mat_t is a dgCMatrix (cells x genes after transpose)
  # scipy is already imported as scipy.sparse
  scipy$csc_matrix(
    tuple(
      np$array(mat_t@x),
      np$array(as.integer(mat_t@i)),
      np$array(as.integer(mat_t@p))
    ),
    shape = tuple(as.integer(nrow(mat_t)), as.integer(ncol(mat_t)))
  )
}

raw_scipy  <- to_scipy_csc(raw_t)
norm_scipy <- to_scipy_csc(norm_t)

# -- Build cell metadata -------------------------------------------------------
cat("Building cell metadata...\n")

meta_cols <- c("eos_cluster", "condition", "sex", "genotype", "orig.ident")
available <- meta_cols[meta_cols %in% colnames(eos@meta.data)]
missing   <- meta_cols[!meta_cols %in% colnames(eos@meta.data)]

if (length(missing) > 0) {
  cat(sprintf("  Warning: metadata columns not found: %s\n",
              paste(missing, collapse = ", ")))
}

obs_df <- eos@meta.data[, available, drop = FALSE]

# Rename orig.ident -> sample for clarity
if ("orig.ident" %in% colnames(obs_df)) {
  colnames(obs_df)[colnames(obs_df) == "orig.ident"] <- "sample"
}

# Rename eos_cluster -> cluster for clarity
if ("eos_cluster" %in% colnames(obs_df)) {
  colnames(obs_df)[colnames(obs_df) == "eos_cluster"] <- "cluster"
}

cat(sprintf("  Metadata columns: %s\n", paste(colnames(obs_df), collapse = ", ")))
cat(sprintf("  %d cells\n", nrow(obs_df)))

# -- Extract UMAP coordinates --------------------------------------------------
cat("Extracting UMAP coordinates...\n")
umap_coords <- Embeddings(eos, reduction = "umap")
cat(sprintf("  UMAP dimensions: %d x %d\n", nrow(umap_coords), ncol(umap_coords)))

# -- Build AnnData object ------------------------------------------------------
cat("Creating AnnData object...\n")

# Convert obs to pandas DataFrame
pd <- import("pandas")
obs_pd <- pd$DataFrame(obs_df)
obs_pd$index <- colnames(eos)

# Gene names as var
var_pd <- pd$DataFrame(index = rownames(eos))

adata <- anndata$AnnData(
  X    = raw_scipy,
  obs  = obs_pd,
  var  = var_pd
)

# Add normalized counts as a layer
adata$layers[["normalized"]] <- norm_scipy

# Add UMAP to obsm
adata$obsm[["X_umap"]] <- np$array(umap_coords)

cat(sprintf("  AnnData: %s\n", adata$shape))

# -- Save h5ad -----------------------------------------------------------------
h5ad_path <- file.path(out_dir, "eos_for_scenic.h5ad")
cat(sprintf("Saving h5ad: %s\n", h5ad_path))
adata$write_h5ad(h5ad_path)

# -- Export cell metadata as CSV -----------------------------------------------
csv_path <- file.path(out_dir, "eos_cell_metadata.csv")
cat(sprintf("Saving cell metadata CSV: %s\n", csv_path))
obs_out <- obs_df
obs_out$barcode <- colnames(eos)
obs_out$UMAP_1  <- umap_coords[, 1]
obs_out$UMAP_2  <- umap_coords[, 2]
write.csv(obs_out, csv_path, row.names = FALSE)

# -- Summary -------------------------------------------------------------------
cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
cat(sprintf("Cells exported:   %d\n", ncol(eos)))
cat(sprintf("Genes exported:   %d\n", nrow(eos)))
cat(sprintf("h5ad file:        %s\n", h5ad_path))
cat(sprintf("Metadata CSV:     %s\n", csv_path))
cat(sprintf("Metadata columns: %s\n", paste(colnames(obs_df), collapse = ", ")))

cat("\nCluster distribution:\n")
print(table(obs_df$cluster))

# -- Export mac data for SCENIC heatmap (script 16) ---------------------------
mac_rds <- file.path(BASE, "mac_subtyped.rds")
if (file.exists(mac_rds)) {
  cat("\nExporting mac data for SCENIC heatmap...\n")
  mac <- readRDS(mac_rds)
  if (inherits(mac[["RNA"]], "Assay5")) mac[["RNA"]] <- JoinLayers(mac[["RNA"]])

  # Metadata CSV
  mac_meta <- data.frame(
    barcode = colnames(mac),
    mac_subtype = mac$mac_subtype,
    row.names = colnames(mac)
  )
  write.csv(mac_meta, file.path(BASE, "mac_subtyped_metadata.csv"))
  cat(sprintf("  Saved mac_subtyped_metadata.csv (%d cells)\n", nrow(mac_meta)))

  # Export mac as h5ad for AUCell scoring in script 16
  mac <- NormalizeData(mac, verbose = FALSE)
  mac_counts <- GetAssayData(mac, layer = "data")
  mac_t <- t(mac_counts)
  mac_scipy <- scipy$csc_matrix(
    tuple(np$array(mac_t@x), np$array(as.integer(mac_t@i)), np$array(as.integer(mac_t@p))),
    shape = tuple(as.integer(nrow(mac_t)), as.integer(ncol(mac_t)))
  )

  pd <- import("pandas")
  mac_obs <- pd$DataFrame(mac_meta)
  mac_obs$index <- colnames(mac)
  mac_var <- pd$DataFrame(index = rownames(mac))

  mac_adata <- anndata$AnnData(X = mac_scipy, obs = mac_obs, var = mac_var)
  mac_h5ad <- file.path(out_dir, "mac_for_scenic_scoring.h5ad")
  mac_adata$write_h5ad(mac_h5ad)
  cat(sprintf("  Saved mac h5ad: %s (%d cells)\n", mac_h5ad, ncol(mac)))
  rm(mac); gc()
}

cat("\nDone! 14_scenic_export.R completed successfully.\n")
