#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 51_contamination_check.R   (Reviewer #1 — is C0/C2 macrophage contamination?)
# On the FULL-gene eos object, per cluster:
#   - eosinophil score (Siglecf/Il5ra/Ccr3/Il4)
#   - macrophage genes split into ABUNDANT/ambient-prone (C1q/Apoe/Lyz2/Folr2/F13a1)
#     vs CELL-INTRINSIC surface (Mertk/Fcgr1/Cd68/Csf1r)
#   - MHC-II score (H2-Aa/H2-Ab1/Cd74) — does C0's stated identity hold?
#   - Scrublet doublet score & %doublet, library size, %mito
# Verdict logic:
#   real contamination/doublets -> intrinsic-mac high + scrublet high + nCount high
#   ambient only                -> abundant-mac high but intrinsic/scrublet normal
# Input : $BASE/eos_subclustered.rds ; Output: $BASE/revision/figures/cluster_robustness/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr); library(tidyr) })

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out  <- file.path(BASE, "revision", "figures", "cluster_robustness")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
if ("C6" %in% levels(Idents(eos))) { eos <- subset(eos, idents="C6", invert=TRUE); Idents(eos) <- droplevels(Idents(eos)) }
if (inherits(eos[["RNA"]], "Assay5")) eos[["RNA"]] <- JoinLayers(eos[["RNA"]])
eos$cl <- Idents(eos)

g <- function(v) intersect(v, rownames(eos))
sets <- list(
  eos_score        = g(c("Siglecf","Il5ra","Ccr3","Il4")),
  mac_abundant     = g(c("C1qa","C1qb","C1qc","Apoe","Lyz2","Folr2","F13a1","Selenop","Pltp")),
  mac_intrinsic    = g(c("Mertk","Fcgr1","Cd68","Csf1r")),
  mhc2_score       = g(c("H2-Aa","H2-Ab1","H2-Eb1","Cd74"))
)
for (nm in names(sets)) { cat(nm, ":", paste(sets[[nm]], collapse=","), "\n") }
for (nm in names(sets)) eos <- AddModuleScore(eos, features = list(sets[[nm]]), name = nm, seed = 1)

# scrublet doublet flag -> logical
sd <- eos$scrublet_doublet
is_dbl <- if (is.logical(sd)) sd else grepl("doublet|true|1", as.character(sd), ignore.case = TRUE)

expr <- GetAssayData(eos, layer = "data")
pct_pos <- function(gene) if (gene %in% rownames(eos)) tapply(expr[gene,] > 0, eos$cl, mean)*100 else NA

summ <- data.frame(cluster = levels(eos$cl)) %>%
  rowwise() %>% mutate(
    n            = sum(eos$cl == cluster),
    eos_score    = round(mean(eos$eos_score1[eos$cl==cluster]), 3),
    mac_abundant = round(mean(eos$mac_abundant1[eos$cl==cluster]), 3),
    mac_intrinsic= round(mean(eos$mac_intrinsic1[eos$cl==cluster]), 3),
    mhc2_score   = round(mean(eos$mhc2_score1[eos$cl==cluster]), 3),
    scrublet     = round(mean(eos$scrublet_score[eos$cl==cluster], na.rm=TRUE), 3),
    pct_doublet  = round(100*mean(is_dbl[eos$cl==cluster]), 1),
    med_nCount   = round(median(eos$nCount_RNA[eos$cl==cluster])),
    med_nFeature = round(median(eos$nFeature_RNA[eos$cl==cluster])),
    pct_mt       = round(mean(eos$percent.mt[eos$cl==cluster]), 2),
    pct_Mertk    = round(pct_pos("Mertk")[cluster], 1),
    pct_Fcgr1    = round(pct_pos("Fcgr1")[cluster], 1)
  ) %>% ungroup()
write.csv(summ, file.path(out, "contamination_by_cluster.csv"), row.names = FALSE)
cat("\n=== per-cluster contamination summary ===\n"); print(as.data.frame(summ))

# ── Figures: module scores + QC by cluster ───────────────────────────────────
md <- eos@meta.data
md$cl <- eos$cl
long <- md %>% select(cl, eos_score1, mac_abundant1, mac_intrinsic1, mhc2_score1) %>%
  pivot_longer(-cl, names_to = "score", values_to = "value")
p1 <- ggplot(long, aes(cl, value, fill = cl)) + geom_violin(scale="width") +
  geom_boxplot(width=0.1, outlier.size=0.2, fill="white") +
  facet_wrap(~score, scales = "free_y") + theme_classic(base_size = 11) +
  theme(legend.position = "none") +
  labs(x = NULL, y = "module score", title = "Identity / contamination scores by eosinophil cluster")
ggsave(file.path(out, "contamination_scores_violin.pdf"), p1, width = 9, height = 7)

p2 <- DotPlot(eos, features = c(sets$eos_score, sets$mhc2_score, sets$mac_abundant, sets$mac_intrinsic)) +
  RotatedAxis() + labs(title = "Eos / MHC-II / macrophage (abundant vs intrinsic) genes by cluster")
ggsave(file.path(out, "contamination_dotplot.pdf"), p2, width = 12, height = 5)

cat("\nDone 51_contamination_check.R\n")
