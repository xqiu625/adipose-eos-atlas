#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 43_composition_biorep.R
# Cluster composition by condition showing individual biological replicates
# Supplement to Figure 2D
#
# SLURM: --partition=batch --mem=32G --cpus-per-task=4 --time=00:30:00
# Conda env: /bigdata/godziklab/shared/Xinru/cellrank_env
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(tidyr)
  library(patchwork)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
fig_dir <- file.path(BASE, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("========================================\n")
cat("43_composition_biorep.R\n")
cat("Cluster composition with biological replicates\n")
cat("========================================\n\n")

# --- Load data ----------------------------------------------------------------
cat("Loading eos_subclustered.rds...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))

if ("C6" %in% levels(Idents(eos))) {
  eos <- subset(eos, idents = "C6", invert = TRUE)
  Idents(eos) <- droplevels(Idents(eos))
  eos$eos_cluster <- droplevels(eos$eos_cluster)
}
cat(sprintf("  %d cells, %d clusters\n", ncol(eos), length(unique(eos$eos_cluster))))

# --- Per-sample composition ---------------------------------------------------
cat("Computing per-sample composition...\n")

comp <- eos@meta.data %>%
  group_by(orig.ident, eos_cluster) %>%
  summarise(count = n(), .groups = "drop")

# Add metadata
sample_meta <- eos@meta.data %>%
  select(orig.ident, condition, sex, genotype, infection) %>%
  distinct()

comp <- comp %>%
  left_join(sample_meta, by = "orig.ident") %>%
  group_by(orig.ident) %>%
  mutate(total = sum(count), proportion = count / total) %>%
  ungroup()

# Ensure all cluster-sample combinations exist
all_combos <- expand.grid(
  orig.ident  = unique(comp$orig.ident),
  eos_cluster = paste0("C", 0:5),
  stringsAsFactors = FALSE
)
comp <- all_combos %>%
  left_join(comp, by = c("orig.ident", "eos_cluster")) %>%
  left_join(sample_meta, by = c("orig.ident"), suffix = c("", ".y")) %>%
  mutate(
    condition  = coalesce(condition, condition.y),
    sex        = coalesce(sex, sex.y),
    genotype   = coalesce(genotype, genotype.y),
    infection  = coalesce(infection, infection.y),
    count      = replace_na(count, 0),
    total      = replace_na(total, 0),
    proportion = replace_na(proportion, 0)
  ) %>%
  select(-ends_with(".y"))

cond_order <- c("WT_Male_HFD", "KO_Male_HFD",
                "WT_Female_HFD", "KO_Female_HFD",
                "WT_Female_HFDInfected", "KO_Female_HFDInfected")
comp$condition <- factor(comp$condition, levels = cond_order)
comp$eos_cluster <- factor(comp$eos_cluster, levels = paste0("C", 0:5))

# Short condition labels
comp$cond_short <- factor(comp$condition, levels = cond_order,
                          labels = c("WT \u2642", "KO \u2642",
                                     "WT \u2640", "KO \u2640",
                                     "WT \u2640 Inf", "KO \u2640 Inf"))

cat(sprintf("  %d samples across %d conditions\n",
            length(unique(comp$orig.ident)), length(cond_order)))

# --- Panel A: Stacked bar per replicate ---------------------------------------
cat("Generating Panel A: stacked bar per replicate...\n")

# Create replicate label (short)
comp <- comp %>%
  group_by(condition) %>%
  mutate(rep_label = paste0(cond_short, "\nRep", dense_rank(orig.ident))) %>%
  ungroup()

# Order replicates by condition
rep_order <- comp %>%
  select(condition, rep_label) %>%
  distinct() %>%
  arrange(condition) %>%
  pull(rep_label)
comp$rep_label <- factor(comp$rep_label, levels = rep_order)

p_stacked <- ggplot(comp, aes(x = rep_label, y = proportion, fill = eos_cluster)) +
  geom_bar(stat = "identity", position = "stack", width = 0.8) +
  scale_fill_manual(values = eos_colors, name = "Cluster") +
  scale_y_continuous(labels = scales::percent, expand = c(0, 0)) +
  labs(x = "", y = "Proportion",
       title = "Cluster Composition per Biological Replicate") +
  theme_pub(base_size = 11) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7))

# Add condition separators
n_per_cond <- comp %>% select(condition, orig.ident) %>% distinct() %>%
  group_by(condition) %>% summarise(n = n())
sep_positions <- cumsum(n_per_cond$n) + 0.5
sep_positions <- sep_positions[-length(sep_positions)]
for (sp in sep_positions) {
  p_stacked <- p_stacked +
    geom_vline(xintercept = sp, linetype = "dashed", color = "grey50", linewidth = 0.4)
}

# --- Panel B: Boxplot per cluster faceted by condition ------------------------
cat("Generating Panel B: boxplot with individual points...\n")

p_boxplot <- ggplot(comp, aes(x = cond_short, y = proportion, fill = cond_short)) +
  geom_boxplot(outlier.shape = NA, alpha = 0.6, linewidth = 0.4) +
  geom_jitter(aes(color = cond_short), width = 0.15, size = 1.5, alpha = 0.8) +
  facet_wrap(~ eos_cluster, nrow = 1, scales = "free_y") +
  scale_fill_manual(values = cond_colors, guide = "none") +
  scale_color_manual(values = cond_colors, guide = "none") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "", y = "Proportion per sample",
       title = "Cluster Proportion by Condition (Biological Replicates)") +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 7),
        strip.text = element_text(face = "bold", size = 11))

# --- Panel C: Dot plot — mean + individual reps -------------------------------
cat("Generating Panel C: dot plot with mean and replicates...\n")

# Compute condition means
cond_means <- comp %>%
  group_by(condition, cond_short, eos_cluster) %>%
  summarise(mean_prop = mean(proportion),
            sd_prop = sd(proportion),
            n = n(), .groups = "drop") %>%
  mutate(se_prop = sd_prop / sqrt(n))

p_dot <- ggplot() +
  # Individual replicates (small, transparent)
  geom_point(data = comp, aes(x = eos_cluster, y = proportion, color = cond_short),
             position = position_dodge(width = 0.7), size = 1.2, alpha = 0.4) +
  # Condition means (larger, solid)
  geom_point(data = cond_means, aes(x = eos_cluster, y = mean_prop, color = cond_short),
             position = position_dodge(width = 0.7), size = 3, shape = 18) +
  # Error bars (SE)
  geom_errorbar(data = cond_means,
                aes(x = eos_cluster, ymin = mean_prop - se_prop,
                    ymax = mean_prop + se_prop, color = cond_short),
                position = position_dodge(width = 0.7), width = 0.2, linewidth = 0.5) +
  scale_color_manual(values = cond_colors, name = "Condition") +
  scale_y_continuous(labels = scales::percent) +
  labs(x = "Cluster", y = "Proportion",
       title = "Cluster Proportion: Mean \u00b1 SE with Biological Replicates") +
  theme_pub(base_size = 11) +
  theme(legend.position = "right")

# --- Save all panels ----------------------------------------------------------
cat("Saving figures...\n")

# Combined
p_combined <- p_stacked / p_boxplot / p_dot +
  plot_annotation(tag_levels = "A",
                  title = "Supplementary: Cluster Composition with Biological Replicates") &
  theme(plot.tag = element_text(face = "bold", size = 14))

save_pdf(p_combined, file.path(fig_dir, "43_composition_biorep_combined.pdf"),
         width = 14, height = 16)

# Individual panels
save_pdf(p_stacked, file.path(fig_dir, "43_composition_biorep_stacked.pdf"),
         width = 12, height = 5)
save_pdf(p_boxplot, file.path(fig_dir, "43_composition_biorep_boxplot.pdf"),
         width = 14, height = 4)
save_pdf(p_dot, file.path(fig_dir, "43_composition_biorep_dotplot.pdf"),
         width = 10, height = 5)

# --- Print summary ------------------------------------------------------------
cat("\n========================================\n")
cat("Summary\n")
cat("========================================\n")
cat("Samples per condition:\n")
print(table(sample_meta$condition))
cat("\nCells per sample:\n")
print(sort(table(eos$orig.ident)))
cat("\n43_composition_biorep.R complete.\n")
