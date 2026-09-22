#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 65_core_downstream.R   (re-run core analyses on decont + >0.10 object)
# Verifies the paragraph's claim that the six-cluster structure and the
# sex/RELMa/infection composition are unchanged after decontamination + filtering.
#   1. Cluster markers (FindAllMarkers) on the decontaminated data
#   2. sccomp compositional shifts: sex (WT F vs M), genotype (KO vs WT, F),
#      infection (HFDInf vs HFD, F)
#   3. Pseudobulk DESeq2 DEG counts for the same three contrasts (whole-eos)
# Input : $REDO/eos_decont010.rds ; Output: $REDO/core/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(dplyr); library(DESeq2); library(sccomp) })

REDO <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/revision/redo"
out  <- file.path(REDO, "core"); dir.create(out, recursive=TRUE, showWarnings=FALSE)

eos <- readRDS(file.path(REDO, "eos_decont010.rds"))
DefaultAssay(eos) <- "decontX"
cat("Loaded eos_decont010.rds:", ncol(eos), "cells\n"); print(table(Idents(eos)))

# ── 1. cluster markers ───────────────────────────────────────────────────────
cat("\n[1] FindAllMarkers...\n")
mk <- FindAllMarkers(eos, only.pos=TRUE, min.pct=0.25, logfc.threshold=0.25, verbose=FALSE)
write.csv(mk, file.path(out, "markers_decont010.csv"), row.names=FALSE)
top <- mk %>% filter(p_val_adj<0.05) %>% group_by(cluster) %>% slice_max(avg_log2FC, n=8)
cat("Top markers per cluster:\n")
for (cl in levels(Idents(eos))) cat(sprintf("  %s: %s\n", cl,
    paste(top$gene[top$cluster==cl], collapse=", ")))
nsig <- as.data.frame(table(cluster = mk$cluster[mk$p_val_adj < 0.05]))
cat("Significant markers per cluster:\n"); print(nsig)

# ── 2. sccomp compositional analysis ─────────────────────────────────────────
cat("\n[2] sccomp compositional shifts...\n")
md <- eos@meta.data
md$cluster <- as.character(Idents(eos))
md$sample  <- md$orig.ident
run_sccomp <- function(sub, formula, label) {
  res <- tryCatch({
    sccomp_estimate(sub, formula_composition = formula, .sample = sample,
                    .cell_group = cluster, cores = 4, verbose = FALSE) |>
      sccomp_test()
  }, error=function(e){ cat("  sccomp", label, "error:", conditionMessage(e),"\n"); NULL })
  if (!is.null(res)) {
    df <- as.data.frame(res)[, intersect(c("cluster","parameter","c_effect","c_FDR"), colnames(as.data.frame(res)))]
    df$contrast <- label
    write.csv(df, file.path(out, paste0("sccomp_", label, ".csv")), row.names=FALSE)
    cat("  ", label, "written (", nrow(df), "rows)\n")
    print(df)
  }
}
sex_sub <- md[md$genotype=="WT" & !grepl("Inf",md$condition), ]
run_sccomp(sex_sub, ~ sex, "sex_WT")
gen_sub <- md[md$sex=="Female" & !grepl("Inf",md$condition), ]
run_sccomp(gen_sub, ~ genotype, "genotype_F")
inf_sub <- md[md$sex=="Female", ]; inf_sub$infected <- grepl("Inf", inf_sub$condition)
run_sccomp(inf_sub, ~ infected, "infection_F")

# ── 3. pseudobulk DESeq2 DEG counts (whole-eos, three contrasts) ─────────────
cat("\n[3] pseudobulk DESeq2...\n")
DefaultAssay(eos) <- "RNA"
cnt <- GetAssayData(eos, layer="counts")
samples <- unique(eos$orig.ident)
pb <- sapply(samples, function(s) Matrix::rowSums(cnt[, eos$orig.ident==s, drop=FALSE]))
colnames(pb) <- samples
meta <- eos@meta.data[!duplicated(eos$orig.ident), ]; rownames(meta) <- meta$orig.ident
meta <- meta[samples, ]
meta$infected <- grepl("Inf", meta$condition)
deg_n <- function(idx, var, label) {
  cd <- meta[idx, ]; m <- pb[, rownames(cd)]
  m <- m[rowSums(m)>10, ]
  cd[[var]] <- factor(cd[[var]])
  if (nlevels(cd[[var]]) < 2) { cat("  ",label,"skip\n"); return() }
  dds <- DESeqDataSetFromMatrix(round(m), cd, as.formula(paste("~", var)))
  dds <- DESeq(dds, quiet=TRUE)
  r <- results(dds)
  n <- sum(r$padj < 0.05 & abs(r$log2FoldChange) > 1, na.rm=TRUE)
  cat(sprintf("  %s: %d sig DEGs (padj<0.05, |LFC|>1)\n", label, n))
  write.csv(as.data.frame(r), file.path(out, paste0("deg_", label, ".csv")))
}
deg_n(meta$genotype=="WT" & !meta$infected, "sex", "sex_WT")
deg_n(meta$sex=="Female" & !meta$infected, "genotype", "genotype_F")
deg_n(meta$sex=="Female", "infected", "infection_F")

cat("\nDone 65.\n")
