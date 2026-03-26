#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))
library(Seurat); library(Matrix)

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"

cat("=== Post-DecontX annotated object ===\n")
obj <- readRDS(file.path(BASE, "all_cells_annotated.rds"))

for (g in c("Igha", "Retnla")) {
  cat(sprintf("\n=== %s ===\n", g))
  if (g %in% rownames(obj)) {
    vals <- GetAssayData(obj, layer = "counts")[g, ]
    cat(sprintf("Total counts: %.1f, Cells expressing: %.1f%%\n",
                sum(vals), 100 * mean(vals > 0)))
    ct <- obj$celltype
    for (t in sort(unique(ct))) {
      cells <- vals[ct == t]
      cat(sprintf("  %-20s n=%5d mean=%.3f %%expr=%.1f%%\n",
                  t, length(cells), mean(cells), 100 * mean(cells > 0)))
    }
  } else {
    cat(g, "NOT in gene list\n")
  }
}

cat("\n=== Before vs After DecontX (sample 2, WT_Female_HFD_1) ===\n")
decon <- readRDS(file.path(BASE, "decontx/decontx_sample_2.rds"))
orig <- Read10X("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR0104/outs/per_sample_outs/2/count/sample_filtered_feature_bc_matrix")
corr <- decon$counts
common <- intersect(colnames(orig), colnames(corr))

for (g in c("Igha", "Retnla")) {
  if (g %in% rownames(orig) && g %in% rownames(corr)) {
    b <- mean(orig[g, common])
    a <- mean(corr[g, common])
    cat(sprintf("%s: before=%.4f after=%.4f reduction=%.1f%%\n",
                g, b, a, 100 * (b - a) / max(b, 1e-9)))
  } else {
    cat(g, "not found in one of the matrices\n")
  }
}

cat("\nDone.\n")
