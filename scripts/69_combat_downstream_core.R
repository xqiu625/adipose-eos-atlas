#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 69_combat_downstream_core.R   (INTEGRATION-based downstream — core)
# Build the ComBat-integrated eosinophil object and re-run the core analyses on
# the INTEGRATED clustering (not the published C0-C5):
#   - attach ComBat clusters (new_cluster) + ComBat embedding + ComBat-UMAP
#   - FindAllMarkers on the integrated clusters (to characterise/name them)
#   - sccomp composition (sex / RELMa / infection) on the integrated clusters
#   - save eos_combat010.rds for the heavier integration-downstream (CellChat,
#     SCENIC, trajectory)
# Inputs : $REDO/eos_decont010.rds, $REDO/combat010_labels.csv,
#          $REDO/eos_combat010.h5ad (obsm X_combat, X_umap)
# Output : $REDO/combat_core/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(sccomp) })
options(lifecycle_verbosity = "quiet")

REDO <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/revision/redo"
out  <- file.path(REDO, "combat_core"); dir.create(out, recursive=TRUE, showWarnings=FALSE)

eos <- readRDS(file.path(REDO, "eos_decont010.rds"))
lab <- read.csv(file.path(REDO, "combat010_labels.csv"), stringsAsFactors=FALSE)
rownames(lab) <- lab$barcode
common <- intersect(colnames(eos), lab$barcode)
eos <- subset(eos, cells = common)
eos$new_cluster <- factor(paste0("i", lab[colnames(eos), "new_cluster"]))   # integrated cluster ids i0..
eos$new_C       <- lab[colnames(eos), "new_C"]
cat("Integration-based clusters:\n"); print(table(eos$new_cluster))

# attach ComBat embedding + ComBat-UMAP from CSV (exported by 67b)
ce <- read.csv(file.path(REDO, "combat010_embedding.csv"), row.names=1, check.names=FALSE)
ce <- ce[colnames(eos), ]
emb <- as.matrix(ce[, grep("^combat_", colnames(ce))]);  colnames(emb) <- paste0("combat_", seq_len(ncol(emb)))
um  <- as.matrix(ce[, c("UMAP_1","UMAP_2")]);             colnames(um)  <- paste0("UMAPcombat_", 1:2)
eos[["combat"]]      <- CreateDimReducObject(embeddings=emb, key="combat_", assay="decontX")
eos[["umap_combat"]] <- CreateDimReducObject(embeddings=um,  key="UMAPcombat_", assay="decontX")

Idents(eos) <- eos$new_cluster
DefaultAssay(eos) <- "decontX"

# ── markers on integrated clusters ───────────────────────────────────────────
cat("\n[1] FindAllMarkers on integrated clusters...\n")
mk <- FindAllMarkers(eos, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, verbose=FALSE)
write.csv(mk, file.path(out, "combat_cluster_markers.csv"), row.names=FALSE)
top <- mk %>% filter(p_val_adj<0.05) %>% group_by(cluster) %>% slice_max(avg_log2FC, n=10)
for (cl in levels(Idents(eos))) cat(sprintf("  %s: %s\n", cl, paste(top$gene[top$cluster==cl], collapse=", ")))
print(as.data.frame(table(cluster = mk$cluster[mk$p_val_adj<0.05])))

# save the integration-based object NOW (before sccomp, which is fragile)
saveRDS(eos, file.path(REDO, "eos_combat010.rds"))
cat("Saved eos_combat010.rds (integration-based object).\n")

# ── sccomp on integrated clusters ────────────────────────────────────────────
cat("\n[2] sccomp composition on integrated clusters...\n")
md <- eos@meta.data; md$cluster <- as.character(eos$new_cluster); md$sample <- as.character(md$orig.ident)
md$infected <- grepl("Inf", md$condition)
run_sccomp <- function(sub, formula, label) {
  sub <- sub[, c("sample","cluster", all.vars(formula))]
  est <- tryCatch(sccomp_estimate(sub, formula_composition=formula, .sample=sample,
                                  .cell_group=cluster, cores=1, verbose=FALSE),
                  error=function(e){cat("  [",label,"] est error:",conditionMessage(e),"\n");NULL})
  if (is.null(est)) return()
  res <- tryCatch(sccomp_test(est), error=function(e) est)
  df <- as.data.frame(res); write.csv(df, file.path(out, paste0("sccomp_", label, "_full.csv")), row.names=FALSE)
  cg <- intersect(c("cell_group","cluster"),colnames(df))[1]; pa<-intersect(c("parameter","factor"),colnames(df))[1]
  ef <- intersect(c("c_effect","effect"),colnames(df))[1]; fd<-intersect(c("c_FDR","FDR"),colnames(df))[1]
  if (!is.na(ef)) { k <- df[!grepl("ntercept",df[[pa]]), c(cg,pa,ef,fd)]; colnames(k)<-c("cluster","parameter","effect","FDR")
    k$sig <- ifelse(k$FDR<0.1,"*",""); cat(sprintf("=== [%s] ===\n",label)); print(k[order(k$cluster),]) }
}
tryCatch({
  run_sccomp(md[md$genotype=="WT" & !md$infected,], ~ sex, "sex_WT")
  run_sccomp(md[md$sex=="Female" & !md$infected,], ~ genotype, "genotype_F")
  run_sccomp(md[md$sex=="Female",], ~ infected, "infection_F")
}, error=function(e) cat("sccomp section failed (object already saved):", conditionMessage(e), "\n"))

saveRDS(eos, file.path(REDO, "eos_combat010.rds"))
cat("\nSaved eos_combat010.rds (integration-based object). Done 69.\n")
