#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --job-name=mac_subtype
#SBATCH --output=mac_subtype_%j.out

# ==============================================================================
# 06_mac_subtype.R
# Macrophage subtyping via M1/M2 module scores
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
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

cat("Subsetting macrophages...\n")
mac <- subset(obj, broad_celltype == "Macrophage")
cat(sprintf("  %d macrophages extracted\n", ncol(mac)))
rm(obj); gc()

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

# --- Save output --------------------------------------------------------------
cat("Saving mac_subtyped.rds...\n")
saveRDS(mac, file.path(BASE, "mac_subtyped.rds"))

cat("Done! 06_mac_subtype.R completed successfully.\n")
