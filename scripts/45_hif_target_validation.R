#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 45_hif_target_validation.R
# Validate HIF sex switch: compare known HIF-1α and HIF-2α target genes
# across WT♀, WT♂, KO♀, KO♂ eosinophils (naive/HFD-only samples)
#
# 1. Differential expression of individual HIF target genes (Wilcoxon)
# 2. Module scoring for HIF-1α vs HIF-2α target gene sets
# 3. GSEA using HALLMARK_HYPOXIA + curated HIF-1/2 gene sets
#
# SLURM: --partition=batch --mem=32G --cpus-per-task=4 --time=01:00:00
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

BASE    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("========================================\n")
cat("45_hif_target_validation.R\n")
cat("========================================\n\n")

# --- Known HIF target genes ---------------------------------------------------
# HIF-1α targets: glycolytic/pro-inflammatory
hif1a_targets <- c("Ldha", "Pdk1", "Pgk1", "Eno1", "Slc2a1", "Aldoa",
                   "Gapdh", "Pkm", "Pfkfb3", "Bnip3", "Vegfa",
                   "Il1b", "Nos2", "Tnf", "Hmox1",
                   "Hk2", "Pfkl", "Tpi1", "Gpi1")

# HIF-2α (Epas1) targets: tissue-protective/anti-inflammatory
hif2a_targets <- c("Arg1", "Cited2", "Epo", "Sox2",
                   "Oct4", "Tgfa", "Ccnd1", "Mmp9",
                   "Vegfa", "Angpt2", "Tie1",
                   "Il10", "Mrc1", "Cd163", "Retnla",
                   "Chi3l3", "Clec10a", "Fn1")

# Shared hypoxia response (from MSigDB HALLMARK_HYPOXIA top genes)
hallmark_hypoxia <- c("Vegfa", "Slc2a1", "Ldha", "Pgk1", "Aldoa", "Eno1",
                      "Gapdh", "Pkm", "Hk2", "Pdk1", "Bnip3", "Pfkfb3",
                      "Ndrg1", "P4ha1", "P4ha2", "Plod2", "Lox",
                      "Egln1", "Egln3", "Adm", "Adora2b",
                      "Hmox1", "Serpine1", "Tgfbi", "Bhlhe40",
                      "Ddit4", "Ero1l", "Plin2", "Ak4",
                      "Gys1", "Pfkp", "Tpi1", "Gpi1")

cat(sprintf("HIF-1a targets: %d genes\n", length(hif1a_targets)))
cat(sprintf("HIF-2a targets: %d genes\n", length(hif2a_targets)))
cat(sprintf("Hallmark hypoxia: %d genes\n", length(hallmark_hypoxia)))

# --- Load data ----------------------------------------------------------------
cat("\nLoading eos_subclustered.rds...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))

# Remove C6
if ("C6" %in% levels(Idents(eos))) {
  eos <- subset(eos, idents = "C6", invert = TRUE)
  Idents(eos) <- droplevels(Idents(eos))
}

# Filter to naive (HFD-only, exclude infected)
naive_conditions <- c("WT_Male_HFD", "WT_Female_HFD", "KO_Male_HFD", "KO_Female_HFD")
eos <- subset(eos, subset = condition %in% naive_conditions)
cat(sprintf("  Naive eosinophils: %d cells\n", ncol(eos)))

# Create sex_geno group
eos$sex_geno <- paste0(eos$genotype, " ", ifelse(eos$sex == "Female", "\u2640", "\u2642"))
eos$sex_geno <- factor(eos$sex_geno, levels = c("WT \u2640", "WT \u2642", "KO \u2640", "KO \u2642"))

cat("Cells per group:\n")
print(table(eos$sex_geno))

# Check which target genes are present
all_genes <- rownames(eos)
hif1a_present <- intersect(hif1a_targets, all_genes)
hif2a_present <- intersect(hif2a_targets, all_genes)
hypoxia_present <- intersect(hallmark_hypoxia, all_genes)

cat(sprintf("\nHIF-1a targets in data: %d / %d\n", length(hif1a_present), length(hif1a_targets)))
cat(sprintf("HIF-2a targets in data: %d / %d\n", length(hif2a_present), length(hif2a_targets)))
cat(sprintf("Hallmark hypoxia in data: %d / %d\n", length(hypoxia_present), length(hallmark_hypoxia)))

# =============================================================================
# 1. Module scores
# =============================================================================
cat("\n--- Module scoring ---\n")

eos <- AddModuleScore(eos, features = list(hif1a_present), name = "HIF1a_score", seed = 42)
eos <- AddModuleScore(eos, features = list(hif2a_present), name = "HIF2a_score", seed = 42)
eos <- AddModuleScore(eos, features = list(hypoxia_present), name = "Hypoxia_score", seed = 42)

# Compute HIF balance: HIF2a - HIF1a (positive = HIF2a dominant)
eos$HIF_balance <- eos$HIF2a_score1 - eos$HIF1a_score1

# Summary stats
cat("\nModule score means per group:\n")
score_summary <- eos@meta.data %>%
  group_by(sex_geno) %>%
  summarise(
    HIF1a = mean(HIF1a_score1),
    HIF2a = mean(HIF2a_score1),
    Hypoxia = mean(Hypoxia_score1),
    HIF_balance = mean(HIF_balance),
    n = n()
  )
print(as.data.frame(score_summary))

# --- Violin plots for module scores ---
group_colors <- c("WT \u2640" = "#B2182B", "WT \u2642" = "#2166AC",
                  "KO \u2640" = "#E08080", "KO \u2642" = "#6BAED6")

p_hif1a <- VlnPlot(eos, features = "HIF1a_score1", group.by = "sex_geno",
                   cols = group_colors, pt.size = 0) +
  labs(title = "HIF-1\u03b1 Target Score", y = "Module Score", x = "") +
  theme_pub() + theme(legend.position = "none")

p_hif2a <- VlnPlot(eos, features = "HIF2a_score1", group.by = "sex_geno",
                   cols = group_colors, pt.size = 0) +
  labs(title = "HIF-2\u03b1 Target Score", y = "Module Score", x = "") +
  theme_pub() + theme(legend.position = "none")

p_balance <- VlnPlot(eos, features = "HIF_balance", group.by = "sex_geno",
                     cols = group_colors, pt.size = 0) +
  labs(title = "HIF Balance (HIF-2\u03b1 \u2212 HIF-1\u03b1)",
       y = "Balance Score", x = "") +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey40") +
  theme_pub() + theme(legend.position = "none")

p_hypoxia <- VlnPlot(eos, features = "Hypoxia_score1", group.by = "sex_geno",
                     cols = group_colors, pt.size = 0) +
  labs(title = "Hallmark Hypoxia Score", y = "Module Score", x = "") +
  theme_pub() + theme(legend.position = "none")

p_modules <- (p_hif1a | p_hif2a) / (p_balance | p_hypoxia) +
  plot_annotation(title = "HIF Target Gene Module Scores by Sex and Genotype",
                  theme = theme(plot.title = element_text(face = "bold", size = 14)))

save_pdf(p_modules, file.path(fig_dir, "45_hif_module_scores.pdf"),
         width = 10, height = 8)

# =============================================================================
# 2. Wilcoxon tests for individual HIF targets
# =============================================================================
cat("\n--- Individual gene testing ---\n")

# Test key genes across groups
key_genes <- unique(c("Hif1a", "Epas1",
                      "Ldha", "Pdk1", "Pgk1", "Eno1", "Slc2a1", "Hk2",
                      "Arg1", "Cited2", "Mrc1", "Retnla", "Il10",
                      "Vegfa", "Bnip3", "Hmox1"))
key_genes <- intersect(key_genes, all_genes)
cat(sprintf("Testing %d key genes\n", length(key_genes)))

# WT Female vs WT Male
Idents(eos) <- "sex_geno"
deg_sex <- FindMarkers(eos, ident.1 = "WT \u2640", ident.2 = "WT \u2642",
                       features = key_genes, min.pct = 0, logfc.threshold = 0,
                       test.use = "wilcox")
deg_sex$gene <- rownames(deg_sex)
deg_sex$comparison <- "WT_F_vs_M"

cat("\nWT Female vs WT Male (key HIF genes):\n")
deg_sex_sig <- deg_sex %>% arrange(p_val_adj)
for (i in seq_len(nrow(deg_sex_sig))) {
  r <- deg_sex_sig[i, ]
  dir <- ifelse(r$avg_log2FC > 0, "Female", "Male")
  sig <- ifelse(r$p_val_adj < 0.05, "*", "")
  cat(sprintf("  %s: log2FC=%+.3f, padj=%.2e %s (%s-enriched)\n",
              r$gene, r$avg_log2FC, r$p_val_adj, sig, dir))
}

# KO Female vs KO Male
deg_ko <- FindMarkers(eos, ident.1 = "KO \u2640", ident.2 = "KO \u2642",
                      features = key_genes, min.pct = 0, logfc.threshold = 0,
                      test.use = "wilcox")
deg_ko$gene <- rownames(deg_ko)
deg_ko$comparison <- "KO_F_vs_M"

cat("\nKO Female vs KO Male (key HIF genes):\n")
deg_ko_sig <- deg_ko %>% arrange(p_val_adj)
for (i in seq_len(nrow(deg_ko_sig))) {
  r <- deg_ko_sig[i, ]
  dir <- ifelse(r$avg_log2FC > 0, "Female", "Male")
  sig <- ifelse(r$p_val_adj < 0.05, "*", "")
  cat(sprintf("  %s: log2FC=%+.3f, padj=%.2e %s (%s-enriched)\n",
              r$gene, r$avg_log2FC, r$p_val_adj, sig, dir))
}

# =============================================================================
# 3. Dot plot of key HIF targets
# =============================================================================
cat("\n--- Dot plot ---\n")

# Split into HIF-1a and HIF-2a target groups
hif1a_key <- intersect(c("Hif1a", "Ldha", "Pdk1", "Pgk1", "Eno1", "Slc2a1", "Hk2", "Bnip3", "Hmox1"), all_genes)
hif2a_key <- intersect(c("Epas1", "Arg1", "Cited2", "Mrc1", "Retnla", "Il10", "Vegfa"), all_genes)
dot_genes <- c(hif1a_key, hif2a_key)

p_dot <- DotPlot(eos, features = dot_genes, group.by = "sex_geno",
                 cols = c("#2166AC", "#B2182B"), dot.scale = 6) +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 9),
        axis.text.y = element_text(face = "bold")) +
  labs(title = "HIF Target Gene Expression by Sex and Genotype") +
  geom_vline(xintercept = length(hif1a_key) + 0.5, linetype = "dashed", color = "grey40")

# Add category labels
p_dot <- p_dot +
  annotate("text", x = length(hif1a_key) / 2, y = 5.2,
           label = "HIF-1\u03b1 Targets", fontface = "bold", size = 3.5, color = "#2166AC") +
  annotate("text", x = length(hif1a_key) + length(hif2a_key) / 2, y = 5.2,
           label = "HIF-2\u03b1 Targets", fontface = "bold", size = 3.5, color = "#B2182B")

save_pdf(p_dot, file.path(fig_dir, "45_hif_target_dotplot.pdf"),
         width = 10, height = 4.5)

# =============================================================================
# 4. Statistical summary: module score comparisons
# =============================================================================
cat("\n--- Statistical tests on module scores ---\n")

# Wilcoxon tests for module scores
for (score_name in c("HIF1a_score1", "HIF2a_score1", "HIF_balance", "Hypoxia_score1")) {
  cat(sprintf("\n%s:\n", score_name))

  # WT F vs WT M
  wt_f <- eos@meta.data[eos$sex_geno == "WT \u2640", score_name]
  wt_m <- eos@meta.data[eos$sex_geno == "WT \u2642", score_name]
  p_wt <- wilcox.test(wt_f, wt_m)$p.value
  cat(sprintf("  WT F vs M: p=%.2e (F mean=%.4f, M mean=%.4f)\n",
              p_wt, mean(wt_f), mean(wt_m)))

  # KO F vs KO M
  ko_f <- eos@meta.data[eos$sex_geno == "KO \u2640", score_name]
  ko_m <- eos@meta.data[eos$sex_geno == "KO \u2642", score_name]
  p_ko <- wilcox.test(ko_f, ko_m)$p.value
  cat(sprintf("  KO F vs M: p=%.2e (F mean=%.4f, M mean=%.4f)\n",
              p_ko, mean(ko_f), mean(ko_m)))

  # WT F vs KO F
  p_f <- wilcox.test(wt_f, ko_f)$p.value
  cat(sprintf("  WT F vs KO F: p=%.2e\n", p_f))

  # WT M vs KO M
  p_m <- wilcox.test(wt_m, ko_m)$p.value
  cat(sprintf("  WT M vs KO M: p=%.2e\n", p_m))
}

# =============================================================================
# 5. Save results
# =============================================================================
cat("\n--- Saving results ---\n")

# Combine DEG results
all_deg <- rbind(deg_sex, deg_ko)
write.csv(all_deg, file.path(fig_dir, "45_hif_target_deg_results.csv"), row.names = FALSE)

# Save module scores per cell
module_df <- eos@meta.data %>%
  select(sex_geno, condition, sex, genotype, eos_cluster,
         HIF1a_score1, HIF2a_score1, Hypoxia_score1, HIF_balance)
write.csv(module_df, file.path(fig_dir, "45_hif_module_scores_per_cell.csv"), row.names = FALSE)

cat("\nOutputs:\n")
cat("  45_hif_module_scores.pdf         — Violin plots of HIF module scores\n")
cat("  45_hif_target_dotplot.pdf        — Dot plot of key HIF target genes\n")
cat("  45_hif_target_deg_results.csv    — Wilcoxon DEG results for key genes\n")
cat("  45_hif_module_scores_per_cell.csv — Module scores per cell\n")
cat("\n45 complete.\n")
