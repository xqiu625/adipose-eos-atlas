#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=highmem
#SBATCH --mem=180G
#SBATCH --cpus-per-task=8
#SBATCH --time=6:00:00

# ==============================================================================
# 56_broad_reference_tier1.R   (Reviewer #1 — broad-reference co-clustering, Tier 1)
# Co-embed OUR cells (eosinophils C0-C5 + macrophages / mast / neutrophils / DC /
# NK / B as negative-control lineages) with EXTERNAL eosinophils (Hu, Gurtner),
# then test:
#   Q1 IDENTITY  — do our eosinophils sit among reference eosinophils while our
#                  other lineages do not? (fraction of each cell's neighbours that
#                  are reference eosinophils, summarised per our lineage)
#   Q2 SUBCLUSTER— do C0-C5 map onto Hu Eos1-8 subgroups? (nearest-reference label
#                  transfer in the integrated space)
# Integration: Harmony on `dataset` (Ours / Hu / Gurtner).
# Output: $BASE/revision/figures/broad_reference/
# ==============================================================================
suppressPackageStartupMessages({
  library(Seurat); library(harmony); library(dplyr); library(ggplot2); library(RANN)
})

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
GD   <- "/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/Nature"
HU   <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/hu_2026/hu_adipose_eos.rds"
out  <- file.path(BASE, "revision", "figures", "broad_reference")
dir.create(out, recursive = TRUE, showWarnings = FALSE)
set.seed(1)

mk_lognorm <- function(o) { DefaultAssay(o) <- "RNA"; if (inherits(o[["RNA"]], "Assay5")) o[["RNA"]] <- JoinLayers(o[["RNA"]]); o }

# ── OURS: all-cells, immune lineages; eosinophils relabelled C0-C5 ────────────
cat("Loading our all-cells...\n")
ours <- mk_lognorm(readRDS(file.path(BASE, "all_cells_annotated.rds")))
keep_lin <- c("Eosinophils","Macrophages","Mast cells","Neutrophils","Dendritic cells","NK cells","B cells")
ours <- subset(ours, subset = celltype %in% keep_lin)
# map C0-C5 onto eosinophils by barcode
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
ec  <- setNames(as.character(eos$eos_cluster), colnames(eos)); rm(eos); gc()
lin <- as.character(ours$celltype)
isE <- lin == "Eosinophils"
lin[isE] <- ifelse(colnames(ours)[isE] %in% names(ec), ec[colnames(ours)[isE]], NA)
ours$grp <- ifelse(is.na(lin), as.character(ours$celltype), lin)
ours <- subset(ours, cells = colnames(ours)[!is.na(lin) | !isE])   # drop eos w/o C0-C5 (the C6/14)
ours$dataset <- "Ours"
ours <- DietSeurat(ours, assays = "RNA")
cat("  ours:", ncol(ours), " grp:", paste(names(table(ours$grp)), collapse=","), "\n")

# ── Hu eosinophils (Eos1-8) ───────────────────────────────────────────────────
cat("Loading Hu...\n")
hu <- mk_lognorm(readRDS(HU))
hu$grp <- "Hu_Eos"; hu$ref_sub <- as.character(hu$New_clusters); hu$dataset <- "Hu"
hu <- DietSeurat(hu, assays = "RNA")

# ── Gurtner eosinophils (organ subtypes) ──────────────────────────────────────
cat("Loading Gurtner...\n")
gu <- mk_lognorm(readRDS(file.path(GD, "Gurtner_eosinophil_subset.rds")))
gu$grp <- "Gurtner_Eos"; gu$ref_sub <- paste0("Gurtner_", as.character(gu$organ)); gu$dataset <- "Gurtner"
gu <- DietSeurat(gu, assays = "RNA")

# ── Common genes + merge ──────────────────────────────────────────────────────
common <- Reduce(intersect, list(rownames(ours), rownames(hu), rownames(gu)))
cat("Common genes:", length(common), "\n")
ours <- ours[common, ]; hu <- hu[common, ]; gu <- gu[common, ]
obj <- merge(ours, list(hu, gu)); rm(ours, hu, gu); gc()
if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])
obj$is_ref_eos <- obj$dataset %in% c("Hu", "Gurtner")
cat("Merged:", ncol(obj), "cells\n"); print(table(obj$dataset))

# ── Standard pipeline + Harmony on dataset ────────────────────────────────────
obj <- NormalizeData(obj, verbose = FALSE) |> FindVariableFeatures(nfeatures = 2000, verbose = FALSE) |>
  ScaleData(verbose = FALSE) |> RunPCA(npcs = 30, verbose = FALSE)
obj <- RunHarmony(obj, group.by.vars = "dataset", dims.use = 1:30, verbose = FALSE)
obj <- RunUMAP(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)

# ── Q1 IDENTITY: fraction of reference-eosinophil neighbours, per our lineage ──
cat("Q1: reference-eos neighbourhood fraction...\n")
emb <- Embeddings(obj, "harmony")[, 1:30]
k <- 30
nn <- RANN::nn2(emb, k = k + 1)$nn.idx[, -1]
refmat <- matrix(obj$is_ref_eos[nn], ncol = k)
obj$frac_ref_eos <- rowMeans(refmat)
our <- obj$dataset == "Ours"
q1 <- data.frame(lineage = obj$grp[our], frac = obj$frac_ref_eos[our]) %>%
  group_by(lineage) %>%
  summarise(n = n(), mean_frac_ref_eos = round(mean(frac), 3), median = round(median(frac), 3)) %>%
  arrange(desc(mean_frac_ref_eos))
write.csv(q1, file.path(out, "q1_identity_refeos_neighbor_fraction.csv"), row.names = FALSE)
cat("\n=== Q1: fraction of neighbours that are reference eosinophils (per our lineage) ===\n"); print(as.data.frame(q1))

# ── Q2 SUBCLUSTER: C0-C5 -> nearest Hu Eos1-8 ─────────────────────────────────
cat("Q2: C0-C5 -> Hu Eos subgroup (nearest reference in Harmony space)...\n")
hu_idx <- which(obj$dataset == "Hu")
our_eos <- which(our & grepl("^C[0-5]$", obj$grp))
nn_hu <- RANN::nn2(emb[hu_idx, ], query = emb[our_eos, ], k = 15)$nn.idx
hu_sub <- obj$ref_sub[hu_idx]
pred_hu <- apply(nn_hu, 1, function(ix) names(sort(table(hu_sub[ix]), decreasing = TRUE))[1])
q2 <- as.data.frame.matrix(table(our_C = obj$grp[our_eos], Hu = pred_hu))
write.csv(q2, file.path(out, "q2_C0C5_to_Hu_Eos.csv"))
cat("\n=== Q2: C0-C5 (rows) vs nearest Hu Eos subgroup (cols) ===\n"); print(q2)

# ── Joint clustering composition (positive control: lineages separate) ────────
obj <- FindNeighbors(obj, reduction = "harmony", dims = 1:30, verbose = FALSE)
obj <- FindClusters(obj, resolution = 0.5, verbose = FALSE)
jc <- as.data.frame.matrix(table(joint = obj$seurat_clusters, grp = obj$grp))
write.csv(jc, file.path(out, "joint_cluster_composition.csv"))

# ── Figures ───────────────────────────────────────────────────────────────────
save_pdf <- function(p, f, w=8, h=6.5) ggsave(file.path(out, f), p, width=w, height=h)
save_pdf(DimPlot(obj, group.by="dataset") + ggtitle("Dataset"), "umap_dataset.pdf")
save_pdf(DimPlot(obj, group.by="grp", label=TRUE, repel=TRUE) + ggtitle("Lineage / cluster") + NoLegend(), "umap_lineage.pdf", 9, 8)
obj$plot_refeos <- obj$frac_ref_eos
save_pdf(FeaturePlot(obj, "plot_refeos") + ggtitle("Fraction of neighbours = reference eosinophils"), "umap_refeos_fraction.pdf")
q1$lineage <- factor(q1$lineage, levels = q1$lineage)
save_pdf(ggplot(q1, aes(lineage, mean_frac_ref_eos, fill = grepl("^C[0-5]$", lineage))) +
  geom_col() + scale_fill_manual(values = c("FALSE"="#999999","TRUE"="#C0392B"), guide="none") +
  coord_flip() + theme_classic(base_size=12) +
  labs(y="mean fraction of neighbours that are reference eosinophils", x=NULL,
       title="Identity: our eosinophil clusters (C0-C5) sit among reference eosinophils",
       subtitle="other lineages (grey) do not"), "q1_identity_bar.pdf", 8, 5)

saveRDS(obj, file.path(out, "broad_reference_tier1.rds"))
cat("\nDone 56_broad_reference_tier1.R. Outputs in", out, "\n")
