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

#SBATCH --partition=batch
#SBATCH --mem=128G
#SBATCH --cpus-per-task=4
#SBATCH --time=4:00:00
#SBATCH --job-name=mac_subtype
#SBATCH --output=mac_subtype_%j.out

# ==============================================================================
# 06_mac_subtype.R  (v2 — macrophage subtype analysis)
#
# Original: M1/M2 module scores + subtype assignment
# Added:
#   - Retnla (RELMa) expression across mac subtypes + conditions
#   - Mac DEG: WT vs KO (RELMa KO impact on macrophages)
#   - Mac subtype composition stats by condition
#   - Cross-cell-type marker comparison (MHC-II, LFA-1, CD124, Cd14)
#   - Mac subtypes mapped onto all-cells UMAP
#
# SLURM: batch, 128G, 4 CPUs, 4h
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(ggrepel)
  library(pheatmap)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# Module score gene sets
m1_genes      <- c("Itgax", "Cd80", "Tnf", "Il6", "Ccl2")
m2_surf_genes <- c("Mrc1", "Lyve1", "Cd301b", "Cd163")
m2_func_genes <- c("Chil3", "Il10", "Arg1", "Tgfb1", "Clec10a")

# --- Load annotated object and subset macrophages ----------------------------
cat("Loading annotated object...\n")
obj <- readRDS(file.path(BASE, "all_cells_annotated.rds"))
if (inherits(obj[["RNA"]], "Assay5")) obj[["RNA"]] <- JoinLayers(obj[["RNA"]])

cat("Subsetting macrophages...\n")
mac <- subset(obj, celltype == "Macrophages")
cat(sprintf("  %d macrophages extracted\n", ncol(mac)))

# --- Normalize and compute module scores -------------------------------------
cat("Normalizing data...\n")
mac <- NormalizeData(mac, verbose = FALSE)

cat("Computing module scores...\n")
mac <- AddModuleScore(mac, features = list(m1_genes),      name = "M1_score",
                      ctrl = 50, verbose = FALSE)
mac <- AddModuleScore(mac, features = list(m2_surf_genes), name = "M2_Surf_score",
                      ctrl = 50, verbose = FALSE)
mac <- AddModuleScore(mac, features = list(m2_func_genes), name = "M2_Func_score",
                      ctrl = 50, verbose = FALSE)

# Seurat appends "1" to module score names
mac$M1_score      <- mac$M1_score1
mac$M2_Surf_score <- mac$M2_Surf_score1
mac$M2_Func_score <- mac$M2_Func_score1

# --- Assign each cell to the highest-scoring subtype -------------------------
cat("Assigning macrophage subtypes...\n")
score_mat <- data.frame(
  M1      = mac$M1_score,
  M2_Surf = mac$M2_Surf_score,
  M2_Func = mac$M2_Func_score
)
max_idx <- apply(score_mat, 1, which.max)
subtype_labels <- c("Mac_M1", "Mac_M2_Surf", "Mac_M2_Func")
mac$mac_subtype <- subtype_labels[max_idx]

cat("  Subtype distribution:\n")
print(table(mac$mac_subtype))

# --- Figures ------------------------------------------------------------------

# UMAP colored by subtype
cat("Plotting UMAP by subtype...\n")
pdf(file.path(fig_dir, "06_mac_umap_subtype.pdf"), width = 7, height = 6)
print(DimPlot(mac, group.by = "mac_subtype", cols = mac_colors, pt.size = 0.5) +
        ggtitle("Macrophage Subtypes") + theme_umap())
dev.off()

# Violin plot of module scores
cat("Plotting module score violins...\n")
violin_df <- mac@meta.data %>%
  select(mac_subtype, M1_score, M2_Surf_score, M2_Func_score) %>%
  pivot_longer(cols = c(M1_score, M2_Surf_score, M2_Func_score),
               names_to = "module", values_to = "score") %>%
  mutate(module = gsub("_score", "", module))

pdf(file.path(fig_dir, "06_mac_module_violins.pdf"), width = 10, height = 5)
print(
  ggplot(violin_df, aes(x = mac_subtype, y = score, fill = mac_subtype)) +
    geom_violin(scale = "width", trim = TRUE) +
    geom_boxplot(width = 0.1, outlier.size = 0.3, fill = "white") +
    facet_wrap(~ module, scales = "free_y") +
    scale_fill_manual(values = mac_colors) +
    labs(x = "Subtype", y = "Module Score",
         title = "Macrophage Module Scores by Subtype") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1),
          legend.position = "none")
)
dev.off()

# Composition by condition
cat("Plotting composition by condition...\n")
comp_df <- mac@meta.data %>%
  group_by(condition, mac_subtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(condition) %>%
  mutate(frac = n / sum(n))

pdf(file.path(fig_dir, "06_mac_composition_condition.pdf"), width = 8, height = 6)
print(
  ggplot(comp_df, aes(x = condition, y = frac, fill = mac_subtype)) +
    geom_bar(stat = "identity", position = "stack") +
    scale_fill_manual(values = mac_colors) +
    labs(x = "Condition", y = "Fraction", fill = "Subtype",
         title = "Macrophage Subtype Composition by Condition") +
    theme_pub() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1))
)
dev.off()

# --- Save mac object ----------------------------------------------------------
cat("Saving mac_subtyped.rds...\n")
saveRDS(mac, file.path(BASE, "mac_subtyped.rds"))


# =============================================================================
# SUPPLEMENTARY ANALYSES
# =============================================================================

# Ensure metadata fields exist
if (!"genotype" %in% colnames(mac@meta.data)) {
  mac$genotype <- ifelse(grepl("^WT_|^WT", mac$orig.ident), "WT", "KO")
}
if (!"Sex" %in% colnames(mac@meta.data)) {
  mac$Sex <- ifelse(grepl("_Male_|Male", mac$orig.ident), "Male", "Female")
}
if (!"infection_status" %in% colnames(mac@meta.data)) {
  mac$infection_status <- ifelse(grepl("Infected", mac$orig.ident), "Infected", "NonInfected")
}

# =============================================================================
# R1: Retnla (RELMa) expression in macrophage subtypes + conditions
# RELMa expression across macrophage subtypes
# =============================================================================
cat("\n=== R1: Retnla (RELMa) expression in macrophages ===\n")

# Violin: Retnla expression by mac subtype
p_retnla_sub <- VlnPlot(mac, features = "Retnla", group.by = "mac_subtype",
                         cols = mac_colors, pt.size = 0) +
  ggtitle("Retnla (RELMa) expression by macrophage subtype") +
  theme_pub() + theme(legend.position = "none")
save_pdf(p_retnla_sub, file.path(fig_dir, "06_mac_retnla_by_subtype.pdf"), 7, 5)

# Violin: Retnla by genotype within each subtype
mac$subtype_genotype <- paste(mac$mac_subtype, mac$genotype, sep = "_")
p_retnla_geno <- VlnPlot(mac, features = "Retnla", group.by = "subtype_genotype",
                          pt.size = 0) +
  ggtitle("Retnla by subtype x genotype") +
  theme_pub() + theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(p_retnla_geno, file.path(fig_dir, "06_mac_retnla_subtype_genotype.pdf"), 10, 5)

# Retnla expression summary stats
retnla_stats <- mac@meta.data %>%
  mutate(retnla_expr = GetAssayData(mac, layer = "data")["Retnla", ]) %>%
  group_by(mac_subtype, genotype) %>%
  summarise(
    n_cells = n(),
    pct_expressing = mean(retnla_expr > 0) * 100,
    mean_expr = mean(retnla_expr),
    .groups = "drop"
  )
cat("  Retnla expression summary:\n")
print(as.data.frame(retnla_stats))
write.csv(retnla_stats, file.path(fig_dir, "06_mac_retnla_stats.csv"), row.names = FALSE)


# =============================================================================
# R2: Mac subtype composition changes: WT vs KO (RELMa KO impact)
# Effect of RELMa absence on macrophage populations
# =============================================================================
cat("\n=== R2: Mac subtype composition — WT vs KO ===\n")

# Per-sample composition (for statistical testing)
sample_comp <- mac@meta.data %>%
  group_by(orig.ident, condition, genotype, Sex, infection_status, mac_subtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(orig.ident) %>%
  mutate(frac = n / sum(n))

write.csv(sample_comp, file.path(fig_dir, "06_mac_composition_per_sample.csv"),
          row.names = FALSE)

# Chi-squared test: subtype distribution WT vs KO (HFD females)
for (sex_val in c("Female", "Male")) {
  wt_cells <- mac@meta.data %>%
    filter(genotype == "WT", Sex == sex_val, infection_status == "NonInfected")
  ko_cells <- mac@meta.data %>%
    filter(genotype == "KO", Sex == sex_val, infection_status == "NonInfected")

  if (nrow(wt_cells) > 0 && nrow(ko_cells) > 0) {
    wt_tab <- table(wt_cells$mac_subtype)
    ko_tab <- table(ko_cells$mac_subtype)
    all_types <- union(names(wt_tab), names(ko_tab))
    cont_table <- rbind(
      WT = wt_tab[all_types], KO = ko_tab[all_types]
    )
    cont_table[is.na(cont_table)] <- 0

    chi_res <- chisq.test(cont_table)
    cat(sprintf("  Chi-sq %s WT vs KO: X2=%.2f, df=%d, p=%.4g\n",
                sex_val, chi_res$statistic, chi_res$parameter, chi_res$p.value))

    # Per-subtype proportion comparison
    for (st in all_types) {
      wt_pct <- round(100 * wt_tab[st] / sum(wt_tab), 1)
      ko_pct <- round(100 * ko_tab[st] / sum(ko_tab), 1)
      cat(sprintf("    %s: WT=%.1f%%, KO=%.1f%%\n", st, wt_pct, ko_pct))
    }
  }
}

# Stacked bar: Mac subtype by genotype (HFD only, split by sex)
hfd_mac <- mac@meta.data %>% filter(infection_status == "NonInfected")
comp_geno <- hfd_mac %>%
  group_by(genotype, Sex, mac_subtype) %>%
  summarise(n = n(), .groups = "drop") %>%
  group_by(genotype, Sex) %>%
  mutate(frac = n / sum(n))

p_comp_geno <- ggplot(comp_geno, aes(x = interaction(genotype, Sex), y = frac,
                                      fill = mac_subtype)) +
  geom_bar(stat = "identity", position = "stack") +
  scale_fill_manual(values = mac_colors) +
  labs(x = "", y = "Fraction", fill = "Subtype",
       title = "Mac subtype: WT vs KO by sex (HFD only)") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
save_pdf(p_comp_geno, file.path(fig_dir, "06_mac_composition_genotype.pdf"), 8, 6)


# =============================================================================
# R3: Pseudobulk DEG — WT vs KO macrophages (RELMa KO impact)
# RELMa KO effect on macrophages
# =============================================================================
cat("\n=== R3: Pseudobulk DESeq2 — WT vs KO macrophages ===\n")

run_mac_deseq2 <- function(mac_sub, test_level, ref_level, group_col,
                            sample_col = "orig.ident") {
  mat <- GetAssayData(mac_sub, layer = "counts")
  meta <- mac_sub@meta.data
  meta$group <- meta[[group_col]]
  meta$sample <- meta[[sample_col]]
  samples <- unique(meta$sample)

  pb_list <- lapply(samples, function(s) {
    cells <- which(meta$sample == s)
    Matrix::rowSums(mat[, cells, drop = FALSE])
  })
  pb_counts <- do.call(cbind, pb_list)
  colnames(pb_counts) <- samples
  pb_counts <- pb_counts[rowSums(pb_counts) >= 10, ]

  col_data <- data.frame(
    sample = samples,
    group  = sapply(samples, function(s) unique(meta$group[meta$sample == s])[1]),
    row.names = samples
  )
  col_data$group <- factor(col_data$group, levels = c(ref_level, test_level))

  cat(sprintf("    Pseudobulk: %d genes, %d samples (%d test, %d ref)\n",
              nrow(pb_counts), length(samples),
              sum(col_data$group == test_level),
              sum(col_data$group == ref_level)))

  dds <- DESeqDataSetFromMatrix(countData = round(pb_counts),
                                colData = col_data, design = ~ group)
  dds <- DESeq(dds)
  res <- as.data.frame(results(dds, contrast = c("group", test_level, ref_level)))
  res$gene <- rownames(res)
  res <- res %>% rename(log2FC = log2FoldChange, pval = pvalue) %>% arrange(padj)
  return(res)
}

# Female: KO vs WT macrophages
mac_fem_hfd <- subset(mac, Sex == "Female" & infection_status == "NonInfected")
cat(sprintf("  Female HFD macs: %d cells\n", ncol(mac_fem_hfd)))

tryCatch({
  deg_mac_fem <- run_mac_deseq2(mac_fem_hfd, "KO", "WT", "genotype")
  n_sig <- sum(!is.na(deg_mac_fem$padj) & deg_mac_fem$padj < 0.05 &
               abs(deg_mac_fem$log2FC) > 0.5)
  cat(sprintf("  Mac DEGs Female (KO vs WT): %d (padj<0.05, |FC|>0.5)\n", n_sig))
  write.csv(deg_mac_fem, file.path(fig_dir, "06_mac_deg_genotype_female.csv"),
            row.names = FALSE)

  # Volcano
  deg_mac_fem$sig <- "NS"
  deg_mac_fem$sig[!is.na(deg_mac_fem$padj) & deg_mac_fem$padj < 0.05 &
                   deg_mac_fem$log2FC > 0.5] <- "Up"
  deg_mac_fem$sig[!is.na(deg_mac_fem$padj) & deg_mac_fem$padj < 0.05 &
                   deg_mac_fem$log2FC < -0.5] <- "Down"

  top_genes <- bind_rows(
    head(deg_mac_fem %>% filter(log2FC > 0) %>% arrange(padj), 10),
    head(deg_mac_fem %>% filter(log2FC < 0) %>% arrange(padj), 10)
  )

  p_vol <- ggplot(deg_mac_fem %>% filter(!is.na(padj)),
                  aes(x = log2FC, y = -log10(pval), color = sig)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 2.5, max.overlaps = 30, color = "black",
                    fontface = "italic") +
    scale_color_manual(values = vol_colors) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    labs(title = "Mac DEG: KO vs WT Female (pseudobulk)",
         x = "log2 Fold Change", y = expression(-log[10](p))) +
    theme_pub()
  save_pdf(p_vol, file.path(fig_dir, "06_mac_volcano_genotype_female.pdf"), 8, 7)
}, error = function(e) {
  cat(sprintf("  ERROR mac DEG female: %s\n", e$message))
})

# Male: KO vs WT macrophages
mac_male_hfd <- subset(mac, Sex == "Male")
cat(sprintf("  Male HFD macs: %d cells\n", ncol(mac_male_hfd)))

tryCatch({
  deg_mac_male <- run_mac_deseq2(mac_male_hfd, "KO", "WT", "genotype")
  n_sig <- sum(!is.na(deg_mac_male$padj) & deg_mac_male$padj < 0.05 &
               abs(deg_mac_male$log2FC) > 0.5)
  cat(sprintf("  Mac DEGs Male (KO vs WT): %d (padj<0.05, |FC|>0.5)\n", n_sig))
  write.csv(deg_mac_male, file.path(fig_dir, "06_mac_deg_genotype_male.csv"),
            row.names = FALSE)

  # Volcano
  deg_mac_male$sig <- "NS"
  deg_mac_male$sig[!is.na(deg_mac_male$padj) & deg_mac_male$padj < 0.05 &
                    deg_mac_male$log2FC > 0.5] <- "Up"
  deg_mac_male$sig[!is.na(deg_mac_male$padj) & deg_mac_male$padj < 0.05 &
                    deg_mac_male$log2FC < -0.5] <- "Down"

  top_genes <- bind_rows(
    head(deg_mac_male %>% filter(log2FC > 0) %>% arrange(padj), 10),
    head(deg_mac_male %>% filter(log2FC < 0) %>% arrange(padj), 10)
  )

  p_vol <- ggplot(deg_mac_male %>% filter(!is.na(padj)),
                  aes(x = log2FC, y = -log10(pval), color = sig)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 2.5, max.overlaps = 30, color = "black",
                    fontface = "italic") +
    scale_color_manual(values = vol_colors) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    labs(title = "Mac DEG: KO vs WT Male (pseudobulk)",
         x = "log2 Fold Change", y = expression(-log[10](p))) +
    theme_pub()
  save_pdf(p_vol, file.path(fig_dir, "06_mac_volcano_genotype_male.pdf"), 8, 7)
}, error = function(e) {
  cat(sprintf("  ERROR mac DEG male: %s\n", e$message))
})

rm(mac_fem_hfd, mac_male_hfd); gc()


# =============================================================================
# R4: Cross-cell-type marker comparison — Eos vs Mac
# MHC2, LFA-1, CD124 as classical macrophage markers
# Cd14 as an ambient-RNA / macrophage signal (not eosinophil)
# =============================================================================
cat("\n=== R4: Eos vs Mac marker comparison ===\n")

# Subset eos from all-cells object (obj still in memory)
eos_obj <- subset(obj, celltype == "Eosinophils")
cat(sprintf("  Eosinophils: %d,  Macrophages: %d\n", ncol(eos_obj), ncol(mac)))

# Markers of interest
shared_markers <- c(
  "H2-Ab1", "H2-Eb1", "H2-Aa", "Cd74",    # MHC-II
  "Itgb2",                                   # LFA-1 (CD18)
  "Itgal",                                   # LFA-1 (CD11a)
  "Il4ra",                                   # CD124
  "Cd14",                                    # monocyte/mac marker (ambient RNA?)
  "Siglecf", "Ccr3",                         # eos-specific
  "Mertk", "Mrc1",                           # mac-specific
  "Retnla"                                   # RELMa
)
shared_markers <- shared_markers[shared_markers %in% rownames(obj)]

# Merge eos + mac for comparison
eos_obj$broad <- "Eosinophil"
mac$broad <- "Macrophage"
eos_mac <- merge(eos_obj, mac)
eos_mac <- NormalizeData(eos_mac, verbose = FALSE)
eos_mac$broad <- factor(eos_mac$broad, levels = c("Eosinophil", "Macrophage"))

# Dot plot: shared markers in eos vs mac
p_shared <- DotPlot(eos_mac, features = shared_markers, group.by = "broad",
                    scale = TRUE) +
  scale_color_gradientn(
    colors = c("#3B4CC0", "#7092D5", "#AAC7FD", "#DDDDDD",
               "#F7B799", "#E8744A", "#B40426"),
    name = "Avg Expr") +
  scale_size_continuous(range = c(2, 10), name = "% Expr") +
  labs(x = "", y = "", title = "Marker expression: Eosinophils vs Macrophages") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic"))
save_pdf(p_shared, file.path(fig_dir, "06_eos_vs_mac_marker_dotplot.pdf"), 12, 5)

# Violin: Cd14 specifically (ambient RNA concern)
p_cd14 <- VlnPlot(eos_mac, features = "Cd14", group.by = "broad",
                   pt.size = 0) +
  ggtitle("Cd14 expression: Eosinophils vs Macrophages") +
  theme_pub() + theme(legend.position = "none")
save_pdf(p_cd14, file.path(fig_dir, "06_cd14_eos_vs_mac.pdf"), 6, 5)

# Quantify: % cells expressing each marker
expr_stats <- data.frame()
for (g in shared_markers) {
  eos_expr <- GetAssayData(eos_obj, layer = "data")[g, ]
  mac_expr <- GetAssayData(mac, layer = "data")[g, ]
  expr_stats <- rbind(expr_stats, data.frame(
    gene = g,
    eos_pct = round(mean(eos_expr > 0) * 100, 1),
    eos_mean = round(mean(eos_expr), 3),
    mac_pct = round(mean(mac_expr > 0) * 100, 1),
    mac_mean = round(mean(mac_expr), 3)
  ))
}
cat("  Marker expression summary (eos vs mac):\n")
print(expr_stats)
write.csv(expr_stats, file.path(fig_dir, "06_eos_vs_mac_marker_stats.csv"),
          row.names = FALSE)

# Eos subcluster-level Cd14 check
eos_sub <- readRDS(file.path(BASE, "eos_subclustered.rds"))
if (inherits(eos_sub[["RNA"]], "Assay5")) eos_sub[["RNA"]] <- JoinLayers(eos_sub[["RNA"]])
eos_sub <- NormalizeData(eos_sub, verbose = FALSE)

p_cd14_eos <- VlnPlot(eos_sub, features = "Cd14", group.by = "eos_cluster",
                       cols = eos_colors, pt.size = 0) +
  ggtitle("Cd14 by eosinophil subcluster") +
  theme_pub() + theme(legend.position = "none")
save_pdf(p_cd14_eos, file.path(fig_dir, "06_cd14_by_eos_cluster.pdf"), 8, 5)

# MHC-II in eos subclusters (relevant for C1 characterization)
for (g in c("H2-Ab1", "H2-Eb1", "Cd74", "Itgb2", "Il4ra")) {
  if (g %in% rownames(eos_sub)) {
    p <- VlnPlot(eos_sub, features = g, group.by = "eos_cluster",
                 cols = eos_colors, pt.size = 0) +
      ggtitle(paste(g, "by eosinophil subcluster")) +
      theme_pub() + theme(legend.position = "none")
    g_safe <- gsub("-", "_", g)
    save_pdf(p, file.path(fig_dir, sprintf("06_%s_by_eos_cluster.pdf", g_safe)), 8, 5)
  }
}
rm(eos_sub); gc()


# =============================================================================
# R5: Mac subtypes on the all-cells UMAP
# Macrophage subtypes used in Figs 3 and 4
# =============================================================================
cat("\n=== R5: Mac subtypes on all-cells UMAP ===\n")

# Transfer mac subtype labels onto the all-cells object
obj$mac_subtype <- NA
obj$mac_subtype[colnames(obj) %in% colnames(mac)] <-
  mac$mac_subtype[match(
    colnames(obj)[colnames(obj) %in% colnames(mac)],
    colnames(mac)
  )]

# Color only macrophages, everything else gray
obj$mac_highlight <- ifelse(is.na(obj$mac_subtype), "Other", obj$mac_subtype)
highlight_colors <- c(mac_colors, "Other" = "grey90")

# Reorder so mac cells plot on top
cell_order <- order(is.na(obj$mac_subtype))

p_mac_on_all <- DimPlot(obj, group.by = "mac_highlight",
                         cols = highlight_colors, pt.size = 0.3,
                         order = c("Mac_M1", "Mac_M2_Func", "Mac_M2_Surf")) +
  ggtitle("Macrophage subtypes on all-cells UMAP") + theme_umap()
save_pdf(p_mac_on_all, file.path(fig_dir, "06_mac_subtypes_on_allcells_umap.pdf"), 10, 8)

# Split by condition
p_mac_split <- DimPlot(obj, group.by = "mac_highlight",
                        cols = highlight_colors, pt.size = 0.2,
                        split.by = "condition", ncol = 3,
                        order = c("Mac_M1", "Mac_M2_Func", "Mac_M2_Surf")) +
  ggtitle("Mac subtypes by condition") + theme_umap() +
  theme(strip.text = element_text(size = 9))
save_pdf(p_mac_split, file.path(fig_dir, "06_mac_subtypes_on_allcells_split.pdf"), 18, 10)

rm(obj, eos_obj, eos_mac); gc()

cat("\nDone! 06_mac_subtype.R (v2) completed successfully.\n")
