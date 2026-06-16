#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))
# Follow-up on the Tier-1 broad reference: (a) cleaner Q2 = pseudobulk Spearman
# correlation of our C0-C5 vs Hu Eos1-8; (b) export per-cell UMAP + meta for figures.
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })

OUT <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/revision/figures/broad_reference"
obj <- readRDS(file.path(OUT, "broad_reference_tier1.rds"))
if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

# ── (a) pseudobulk correlation: our C0-C5 vs Hu Eos1-8 ───────────────────────
hvg  <- VariableFeatures(obj)
expr <- as.matrix(GetAssayData(obj, layer = "data")[hvg, ])
grp  <- obj$grp
ours_eos <- paste0("C", 0:5)
hu_sub <- ifelse(obj$dataset == "Hu", as.character(obj$ref_sub), NA)

pb <- function(mask) rowMeans(expr[, mask, drop = FALSE])
ours_pb <- sapply(ours_eos, function(c) pb(grp == c & obj$dataset == "Ours"))
hu_levels <- sort(unique(na.omit(hu_sub)))
hu_pb <- sapply(hu_levels, function(h) pb(which(hu_sub == h)))

cormat <- cor(ours_pb, hu_pb, method = "spearman")
write.csv(round(cormat, 3), file.path(OUT, "q2_pseudobulk_correlation_C0C5_vs_Hu.csv"))
cat("=== Q2 pseudobulk Spearman correlation: C0-C5 (rows) x Hu Eos (cols) ===\n")
print(round(cormat, 3))
cat("\nBest Hu match per cluster:\n")
for (c in ours_eos) cat(sprintf("  %s -> %s (r=%.3f)\n", c, colnames(cormat)[which.max(cormat[c,])], max(cormat[c,])))

# ── (b) export per-cell UMAP + meta (downsample non-eos for plotting) ─────────
um <- as.data.frame(Embeddings(obj, "umap")); colnames(um) <- c("UMAP_1","UMAP_2")
um$dataset <- obj$dataset; um$grp <- obj$grp; um$frac_ref_eos <- obj$frac_ref_eos
# lineage category for plotting
um$category <- dplyr::case_when(
  um$dataset == "Hu" ~ "Ref eos (Hu)",
  um$dataset == "Gurtner" ~ "Ref eos (Gurtner)",
  grepl("^C[0-5]$", um$grp) ~ "Our eosinophils",
  um$grp == "Macrophages" ~ "Our macrophages",
  TRUE ~ "Our other lineages")
set.seed(1)
# keep all eos + refs; downsample big non-eos groups to <=3000 for the figure
keep <- unlist(lapply(split(seq_len(nrow(um)), um$grp), function(ix)
  if (length(ix) > 3000 && !grepl("^C[0-5]$", um$grp[ix[1]])) sample(ix, 3000) else ix))
write.csv(um[sort(keep), ], file.path(OUT, "umap_cells_for_fig.csv"), row.names = FALSE)
cat("\nExported", length(keep), "cells for figure.\nDone 57.\n")
