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

# ==============================================================================
# 66_c2_robustness_sexgene_gem.R   (Reviewer #1 / Meera — defend C2)
# C2 is >90% HFD male. Two tests:
#   (A) GEM check — are the male HFD samples (and the C2 cells) confined to one
#       GEM reaction? If C2 spans multiple GEMs, its male enrichment is biology,
#       not a batch artifact.
#   (B) Sex-gene-removal — remove Y-linked / X-inactivation genes, re-cluster, and
#       test whether a C2-equivalent cluster still emerges. If C2 survives, it is
#       not merely "maleness" driven by a handful of sex genes.
#   (C) How many of C2's markers are sex-linked genes?
# Input : $REDO/eos_decont010.rds (decont + >0.10; C0-C5)
# Output: $REDO/c2/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr) })

REDO <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/revision/redo"
out  <- file.path(REDO, "c2"); dir.create(out, recursive=TRUE, showWarnings=FALSE)
eos <- readRDS(file.path(REDO, "eos_decont010.rds"))
DefaultAssay(eos) <- "decontX"
eos$cluster <- as.character(Idents(eos))

sex_genes <- intersect(c("Ddx3y","Uty","Eif2s3y","Kdm5d","Uba1y","Zfy1","Zfy2","Usp9y",
                         "Xist","Tsix"), rownames(eos))
cat("Sex-linked genes present:", paste(sex_genes, collapse=", "), "\n")

# ── (A) GEM check ────────────────────────────────────────────────────────────
cat("\n=== (A) GEM reaction distribution ===\n")
cat("Male HFD samples -> GEM:\n")
male_hfd <- eos@meta.data %>% filter(sex=="Male", !grepl("Inf",condition)) %>%
  select(orig.ident, gem_reaction) %>% distinct() %>% arrange(orig.ident)
print(male_hfd)
cat("\nC2 cells by GEM reaction:\n")
print(table(eos$gem_reaction[eos$cluster=="C2"]))
cat("\nC2 cells by sample:\n")
print(table(eos$orig.ident[eos$cluster=="C2"]))
cat("\nAll clusters x GEM (cells):\n")
print(table(eos$cluster, eos$gem_reaction))
write.csv(as.data.frame.matrix(table(eos$cluster, eos$gem_reaction)),
          file.path(out, "cluster_by_gem.csv"))
write.csv(male_hfd, file.path(out, "male_hfd_gem_map.csv"), row.names=FALSE)

# ── (B) sex-gene-removal re-clustering ───────────────────────────────────────
cat("\n=== (B) Re-cluster after removing sex-linked genes ===\n")
sub <- DietSeurat(eos, assays="decontX")
sub <- NormalizeData(sub, verbose=FALSE)
sub <- FindVariableFeatures(sub, nfeatures=2000, verbose=FALSE)
VariableFeatures(sub) <- setdiff(VariableFeatures(sub), sex_genes)   # drop sex genes
sub <- ScaleData(sub, verbose=FALSE)
sub <- RunPCA(sub, npcs=30, verbose=FALSE)
sub <- FindNeighbors(sub, dims=1:20, verbose=FALSE)
sub <- FindClusters(sub, resolution=0.4, verbose=FALSE)
ct <- table(orig=eos$cluster, new=sub$seurat_clusters)
cat("original cluster (rows) x new cluster without sex genes (cols):\n"); print(ct)
write.csv(as.data.frame.matrix(ct), file.path(out, "c2_nosexgene_crosstab.csv"))
# does a new cluster correspond to C2?  (C2 cells' dominant new cluster + its purity)
c2_new <- sub$seurat_clusters[eos$cluster=="C2"]
dom <- names(sort(table(c2_new), decreasing=TRUE))[1]
purity_c2_stays <- mean(c2_new == dom)
new_cluster_c2frac <- mean(eos$cluster[sub$seurat_clusters==dom] == "C2")
cat(sprintf("\nC2 cells: %.0f%% land in new cluster %s; that new cluster is %.0f%% original-C2\n",
            100*purity_c2_stays, dom, 100*new_cluster_c2frac))

# ── (C) C2 markers: how many are sex-linked? ─────────────────────────────────
cat("\n=== (C) C2 markers vs sex genes ===\n")
DefaultAssay(eos) <- "decontX"
c2mk <- FindMarkers(eos, ident.1="C2", only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, verbose=FALSE)
c2mk$gene <- rownames(c2mk); c2sig <- c2mk[c2mk$p_val_adj<0.05, ]
cat(sprintf("C2 significant markers: %d ; of which sex-linked: %d (%s)\n",
            nrow(c2sig), sum(c2sig$gene %in% sex_genes),
            paste(intersect(c2sig$gene, sex_genes), collapse=", ")))
write.csv(c2sig, file.path(out, "c2_markers.csv"), row.names=FALSE)
cat("Top 15 C2 markers:", paste(head(c2sig$gene[order(-c2sig$avg_log2FC)],15), collapse=", "), "\n")

cat("\nDone 66.\n")
