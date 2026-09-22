#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 47_c2_robustness.R   (Reviewer #1, Concern #1)
# Is cluster C2 (>90% WT_Male) a genuine eosinophil state, or a male-specific
# transcriptional / stress response (immediate-early + UPR genes)?
#   - C2 marker genes (C2 vs rest)
#   - fraction of each cluster's markers that are IEG / stress genes
#   - per-cluster IEG/stress module score
# Input : $BASE/eos_subclustered.rds (C0-C5)
# Output: $BASE/revision/figures/cluster_robustness/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr) })

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out  <- file.path(BASE, "revision", "figures", "cluster_robustness")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
if ("C6" %in% levels(Idents(eos))) { eos <- subset(eos, idents = "C6", invert = TRUE); Idents(eos) <- droplevels(Idents(eos)) }
if (inherits(eos[["RNA"]], "Assay5")) eos[["RNA"]] <- JoinLayers(eos[["RNA"]])
clusters <- levels(Idents(eos))
cat("Clusters:", paste(clusters, collapse=","), " n=", ncol(eos), "\n")

# ── IEG / stress gene set (AP-1 / immediate-early + heat-shock / UPR) ──────────
ieg <- c("Fos","Fosb","Fosl1","Fosl2","Jun","Junb","Jund","Egr1","Egr2","Atf3",
         "Ier2","Ier3","Dusp1","Zfp36","Nr4a1","Nr4a2","Nr4a3","Btg2","Gadd45b",
         "Socs3","Klf2","Klf4","Klf6","Cebpb",
         "Hspa1a","Hspa1b","Hspb1","Dnajb1","Dnaja1","Hsph1","Hspa8",
         "Atf4","Ddit3","Atf6","Xbp1","Herpud1")
ieg <- intersect(ieg, rownames(eos))
cat("IEG/stress genes present:", length(ieg), "\n")

# ── Per-cluster markers + IEG fraction ────────────────────────────────────────
allmark <- FindAllMarkers(eos, only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
write.csv(allmark, file.path(out, "eos_cluster_markers.csv"), row.names = FALSE)

frac <- lapply(clusters, function(cl) {
  m  <- allmark %>% filter(cluster == cl, p_val_adj < 0.05) %>% arrange(desc(avg_log2FC))
  topm <- head(m$gene, 50)
  data.frame(cluster = cl, n_sig_markers = nrow(m),
             n_top = length(topm),
             n_ieg_in_top = sum(topm %in% ieg),
             ieg_frac_top = round(mean(topm %in% ieg), 3),
             top_ieg = paste(intersect(topm, ieg), collapse=";"))
}) %>% bind_rows()
write.csv(frac, file.path(out, "c2_ieg_fraction_by_cluster.csv"), row.names = FALSE)
cat("\n=== IEG fraction among top-50 markers per cluster ===\n"); print(frac[,1:5])

# ── C2 markers explicitly ─────────────────────────────────────────────────────
c2 <- FindMarkers(eos, ident.1 = "C2", only.pos = TRUE, min.pct = 0.25, logfc.threshold = 0.25, verbose = FALSE)
c2$gene <- rownames(c2); c2$is_ieg <- c2$gene %in% ieg
write.csv(c2, file.path(out, "c2_markers.csv"), row.names = FALSE)
cat(sprintf("\nC2: %d sig markers; %d (%.0f%%) are IEG/stress; top-20 IEG: %s\n",
    sum(c2$p_val_adj<0.05), sum(c2$is_ieg & c2$p_val_adj<0.05),
    100*mean(head(c2$gene[c2$p_val_adj<0.05],20) %in% ieg),
    paste(head(c2$gene[c2$is_ieg],20), collapse=", ")))

# ── Per-cluster IEG/stress module score ──────────────────────────────────────
eos <- AddModuleScore(eos, features = list(ieg), name = "IEG_stress", seed = 1)
df <- data.frame(cluster = Idents(eos), score = eos$IEG_stress1)
sm <- df %>% group_by(cluster) %>% summarise(mean = mean(score), median = median(score)) %>% arrange(desc(mean))
write.csv(sm, file.path(out, "ieg_module_score_by_cluster.csv"), row.names = FALSE)
cat("\n=== IEG/stress module score by cluster ===\n"); print(as.data.frame(sm))

p <- ggplot(df, aes(cluster, score, fill = cluster)) +
  geom_violin(scale="width", trim=TRUE) + geom_boxplot(width=0.12, outlier.size=0.2, fill="white") +
  labs(y = "IEG / stress module score", x = NULL,
       title = "Immediate-early / stress signature by eosinophil cluster",
       subtitle = "High C2 score would support reviewer's 'transcriptional response, not a state'") +
  theme_classic(base_size = 12) + theme(legend.position = "none")
ggsave(file.path(out, "ieg_module_violin.pdf"), p, width = 7, height = 5)

cat("\nDone 47_c2_robustness.R. Outputs in", out, "\n")
