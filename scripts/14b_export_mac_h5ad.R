#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

suppressPackageStartupMessages({
  library(Seurat)
  library(reticulate)
})

use_condaenv("/bigdata/godziklab/shared/Xinru/cellrank_env")
anndata <- import("anndata")
np      <- import("numpy")
scipy   <- import("scipy.sparse")

BASE    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out_dir <- file.path(BASE, "scenic")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

cat("14b_export_mac_h5ad.R\n")
cat("Export macrophage cells for regulon scoring\n\n")

# Try loading Mac-specific object first, fall back to full atlas
mac_rds <- file.path(BASE, "mac_subtyped.rds")
atlas_rds <- file.path(BASE, "eos_mac_atlas.rds")

if (file.exists(mac_rds)) {
  cat(sprintf("Loading %s...\n", mac_rds))
  mac <- readRDS(mac_rds)
} else if (file.exists(atlas_rds)) {
  cat(sprintf("Mac RDS not found; loading atlas: %s\n", atlas_rds))
  mac <- readRDS(atlas_rds)
  # Keep only Mac subtypes
  mac_idents <- c("Mac_M1", "Mac_M2_Func", "Mac_M2_Surf")
  mac <- subset(mac, idents = mac_idents)
} else {
  # Try harmony_clustered if atlas not found
  alt_rds <- file.path(BASE, "harmony_clustered.rds")
  if (!file.exists(alt_rds)) {
    # Look in parent directory
    alt_rds <- file.path(dirname(BASE), "harmony_clustered.rds")
  }
  cat(sprintf("Loading %s...\n", alt_rds))
  obj <- readRDS(alt_rds)
  mac_idents <- c("Mac_M1", "Mac_M2_Func", "Mac_M2_Surf")
  Idents(obj) <- obj@meta.data[["cell_type"]]
  mac <- subset(obj, idents = mac_idents)
  rm(obj); gc()
}

cat(sprintf("  Mac cells: %d, genes: %d\n", ncol(mac), nrow(mac)))
cat(sprintf("  Mac subtypes: %s\n", paste(unique(Idents(mac)), collapse = ", ")))

# Extract raw counts
cat("Extracting counts...\n")
raw_counts <- GetAssayData(mac, layer = "counts")
raw_t <- t(raw_counts)

to_scipy_csc <- function(mat_t) {
  scipy$csc_matrix(
    tuple(np$array(mat_t@x),
          np$array(as.integer(mat_t@i)),
          np$array(as.integer(mat_t@p))),
    shape = tuple(as.integer(nrow(mat_t)), as.integer(ncol(mat_t)))
  )
}
raw_scipy <- to_scipy_csc(raw_t)

# Metadata: include mac_subtype column
obs_df <- mac@meta.data
# Standardize subtype column name
if (!"mac_subtype" %in% colnames(obs_df)) {
  obs_df$mac_subtype <- as.character(Idents(mac))
}
obs_df <- obs_df[, intersect(c("mac_subtype", "condition", "sex", "genotype", "orig.ident"),
                              colnames(obs_df)), drop = FALSE]

# UMAP
umap_key <- grep("^umap", colnames(mac@reductions), ignore.case = TRUE, value = TRUE)[1]
if (!is.na(umap_key)) {
  umap_coords <- mac@reductions[[umap_key]]@cell.embeddings[, 1:2]
  obsm_list <- list("X_umap" = umap_coords)
} else {
  obsm_list <- NULL
}

# Build AnnData
cat("Building AnnData...\n")
var_df <- data.frame(row.names = rownames(mac))
r_obs_df <- reticulate::r_to_py(obs_df)

adata <- anndata$AnnData(
  X    = raw_scipy,
  obs  = r_obs_df,
  var  = reticulate::r_to_py(var_df)
)
adata$obs_names <- reticulate::r_to_py(colnames(mac))
adata$var_names <- reticulate::r_to_py(rownames(mac))

if (!is.null(obsm_list)) {
  adata$obsm[["X_umap"]] <- np$array(umap_coords)
}

out_path <- file.path(out_dir, "mac_for_scenic_scoring.h5ad")
adata$write_h5ad(out_path)
cat(sprintf("Saved: %s\n", out_path))

# Also save metadata CSV for easy loading
meta_path <- file.path(BASE, "mac_subtyped_metadata.csv")
write.csv(obs_df, meta_path)
cat(sprintf("Saved metadata: %s\n", meta_path))
cat("\n14b complete.\n")
