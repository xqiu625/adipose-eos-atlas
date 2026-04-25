#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 45b_hif_scenic_targets.R
# Validate Fig 5C using ACTUAL SCENIC-identified regulon target genes
# (not generic literature targets)
#
# Extracts the real target genes for each TF in Fig 5C from the per-cluster
# pySCENIC regulons.csv files, then:
#   1. Module scoring with SCENIC-identified targets
#   2. Violin plots across WT♀/WT♂/KO♀/KO♂ (naive only)
#   3. Statistical tests
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

BASE       <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
SCENIC_DIR <- file.path(BASE, "scenic_per_cluster")
fig_dir    <- file.path(BASE, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

EOS_CLUSTERS <- c("C0", "C1", "C2", "C3", "C4", "C5")

# The 10 TFs from Fig 5C
FIG5C_TFS <- c("Tbx21", "Epas1", "Runx1", "Cebpa", "Stat5a",
               "Cebpe", "Hif1a", "Gata2", "Xbp1", "Stat3")
FEMALE_TFS <- c("Tbx21", "Epas1", "Runx1", "Cebpa", "Stat5a")
MALE_TFS   <- c("Cebpe", "Hif1a", "Gata2", "Xbp1", "Stat3")

cat("========================================\n")
cat("45b_hif_scenic_targets.R\n")
cat("Validate Fig 5C with SCENIC regulon targets\n")
cat("========================================\n\n")

# =============================================================================
# 1. Extract SCENIC regulon target genes from regulons.csv
# =============================================================================
cat("Step 1: Extracting SCENIC regulon targets...\n")

# Load pySCENIC regulons using raw deserialization to get gene lists
# The regulons.csv from pySCENIC ctx can be loaded with pyscenic utils
# But since we're in R, we'll parse the AUCell column names + use the
# adjacency/regulon files

# Alternative approach: load from Python-generated regulon data
# For each cluster, load regulons.csv and extract target genes
# pySCENIC regulons.csv is a complex format — instead, let's use
# the adj.csv (GRNBoost2 adjacencies) filtered by the regulons

# Best approach: use Python to extract, save as simple CSV, then load in R
# Let's write a helper Python script inline

py_script <- tempfile(fileext = ".py")
writeLines(c(
  "import os, sys, json",
  "import pandas as pd",
  "",
  sprintf("SCENIC_DIR = '%s'", SCENIC_DIR),
  sprintf("FIG5C_TFS = %s", paste0("[", paste(sprintf("'%s'", FIG5C_TFS), collapse=","), "]")),
  "EOS_CLUSTERS = ['C0','C1','C2','C3','C4','C5']",
  "",
  "# Try loading regulons from each cluster",
  "all_regulons = {}  # tf -> set of target genes",
  "",
  "for cl in EOS_CLUSTERS:",
  "    reg_path = os.path.join(SCENIC_DIR, cl, 'regulons.csv')",
  "    if not os.path.exists(reg_path):",
  "        continue",
  "    try:",
  "        from pyscenic.utils import load_motifs",
  "        from pyscenic.prune import df2regulons",
  "        motifs = load_motifs(reg_path)",
  "        regs = df2regulons(motifs)",
  "        for r in regs:",
  "            tf = r.name.split('(')[0].strip()",
  "            if tf in FIG5C_TFS:",
  "                genes = list(r.genes)",
  "                if tf not in all_regulons or len(genes) > len(all_regulons[tf]):",
  "                    all_regulons[tf] = genes",
  "                    print(f'  {cl}: {tf} -> {len(genes)} targets', file=sys.stderr)",
  "    except Exception as e:",
  "        print(f'  {cl}: error: {e}', file=sys.stderr)",
  "",
  "# Save as simple CSV: tf, targets (comma-separated)",
  sprintf("out_path = '%s'", file.path(fig_dir, "45b_scenic_regulon_targets.csv")),
  "rows = []",
  "for tf in FIG5C_TFS:",
  "    if tf in all_regulons:",
  "        rows.append({'tf': tf, 'n_targets': len(all_regulons[tf]),",
  "                     'targets': ';'.join(sorted(all_regulons[tf]))})",
  "        print(f'{tf}: {len(all_regulons[tf])} targets')",
  "    else:",
  "        print(f'{tf}: NOT FOUND in any cluster regulon')",
  "pd.DataFrame(rows).to_csv(out_path, index=False)",
  "print(f'Saved: {out_path}')"
), py_script)

# Run Python to extract regulon targets
python_bin <- "/bigdata/godziklab/shared/Xinru/cellrank_env/bin/python"
cat("  Running Python to extract regulon targets...\n")
system2(python_bin, py_script, env = c(
  paste0("LD_LIBRARY_PATH=/bigdata/godziklab/shared/Xinru/cellrank_env/lib:",
         "/opt/linux/rocky/8.x/x86_64/pkgs/gcc/12.2.0/lib64:",
         Sys.getenv("LD_LIBRARY_PATH"))
))

# Load the extracted targets
targets_file <- file.path(fig_dir, "45b_scenic_regulon_targets.csv")
if (!file.exists(targets_file)) {
  stop("Failed to extract regulon targets. Check Python output above.")
}

targets_df <- read.csv(targets_file)
cat(sprintf("\n  Extracted targets for %d TFs\n", nrow(targets_df)))

# Parse into named list
regulon_genes <- list()
for (i in seq_len(nrow(targets_df))) {
  tf <- targets_df$tf[i]
  genes <- strsplit(targets_df$targets[i], ";")[[1]]
  regulon_genes[[tf]] <- genes
  cat(sprintf("  %s: %d target genes\n", tf, length(genes)))
}

# =============================================================================
# 2. Load Seurat object, filter to naive
# =============================================================================
cat("\nStep 2: Loading eos data...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))

if ("C6" %in% levels(Idents(eos))) {
  eos <- subset(eos, idents = "C6", invert = TRUE)
  Idents(eos) <- droplevels(Idents(eos))
}

naive_conditions <- c("WT_Male_HFD", "WT_Female_HFD", "KO_Male_HFD", "KO_Female_HFD")
eos <- subset(eos, subset = condition %in% naive_conditions)

eos$sex_geno <- paste0(eos$genotype, " ", ifelse(eos$sex == "Female", "\u2640", "\u2642"))
eos$sex_geno <- factor(eos$sex_geno, levels = c("WT \u2640", "WT \u2642", "KO \u2640", "KO \u2642"))

cat(sprintf("  %d naive eosinophils\n", ncol(eos)))
print(table(eos$sex_geno))

# Check gene overlap
all_genes <- rownames(eos)
for (tf in names(regulon_genes)) {
  present <- intersect(regulon_genes[[tf]], all_genes)
  cat(sprintf("  %s: %d / %d targets in expression data\n",
              tf, length(present), length(regulon_genes[[tf]])))
  regulon_genes[[tf]] <- present  # keep only present genes
}

# =============================================================================
# 3. Module scoring with SCENIC targets
# =============================================================================
cat("\nStep 3: Module scoring with SCENIC-identified targets...\n")

# Score each TF's regulon
for (tf in names(regulon_genes)) {
  if (length(regulon_genes[[tf]]) < 3) {
    cat(sprintf("  Skipping %s (< 3 targets)\n", tf))
    next
  }
  score_name <- paste0(tf, "_regulon")
  eos <- AddModuleScore(eos, features = list(regulon_genes[[tf]]),
                        name = score_name, seed = 42)
  # Rename the auto-numbered column
  col_auto <- paste0(score_name, "1")
  if (col_auto %in% colnames(eos@meta.data)) {
    eos@meta.data[[score_name]] <- eos@meta.data[[col_auto]]
  }
  cat(sprintf("  %s: scored with %d targets\n", tf, length(regulon_genes[[tf]])))
}

# =============================================================================
# 4. Violin plots — one per TF
# =============================================================================
cat("\nStep 4: Generating violin plots...\n")

group_colors <- c("WT \u2640" = "#B2182B", "WT \u2642" = "#2166AC",
                  "KO \u2640" = "#E08080", "KO \u2642" = "#6BAED6")

plot_list <- list()
for (tf in FIG5C_TFS) {
  score_col <- paste0(tf, "_regulon")
  if (!score_col %in% colnames(eos@meta.data)) next

  p <- VlnPlot(eos, features = score_col, group.by = "sex_geno",
               cols = group_colors, pt.size = 0) +
    labs(title = paste0(tf, " Regulon\n(", length(regulon_genes[[tf]]), " targets)"),
         y = "Module Score", x = "") +
    theme_pub(base_size = 9) +
    theme(legend.position = "none",
          plot.title = element_text(size = 10, face = "bold.italic"))

  plot_list[[tf]] <- p
}

# Arrange: female-enriched on top row, male-enriched on bottom
n_plots <- length(plot_list)
if (n_plots > 0) {
  p_combined <- wrap_plots(plot_list, ncol = 5, nrow = 2) +
    plot_annotation(
      title = "Fig 5C Validation: SCENIC Regulon Target Expression",
      subtitle = "Module scores computed from SCENIC-identified target genes (naive samples only)",
      theme = theme(
        plot.title = element_text(face = "bold", size = 14),
        plot.subtitle = element_text(size = 10, color = "grey40")
      )
    )

  save_pdf(p_combined, file.path(fig_dir, "45b_scenic_target_violins.pdf"),
           width = 16, height = 7)
}

# =============================================================================
# 5. Statistical tests
# =============================================================================
cat("\nStep 5: Statistical tests...\n")

stats_results <- list()
for (tf in FIG5C_TFS) {
  score_col <- paste0(tf, "_regulon")
  if (!score_col %in% colnames(eos@meta.data)) next

  wt_f <- eos@meta.data[eos$sex_geno == "WT \u2640", score_col]
  wt_m <- eos@meta.data[eos$sex_geno == "WT \u2642", score_col]
  ko_f <- eos@meta.data[eos$sex_geno == "KO \u2640", score_col]
  ko_m <- eos@meta.data[eos$sex_geno == "KO \u2642", score_col]

  p_sex_wt <- wilcox.test(wt_f, wt_m)$p.value
  p_sex_ko <- wilcox.test(ko_f, ko_m)$p.value
  p_geno_f <- wilcox.test(wt_f, ko_f)$p.value
  p_geno_m <- wilcox.test(wt_m, ko_m)$p.value

  expected_dir <- ifelse(tf %in% FEMALE_TFS, "Female", "Male")
  observed_dir <- ifelse(mean(wt_f) > mean(wt_m), "Female", "Male")
  consistent <- expected_dir == observed_dir

  cat(sprintf("\n  %s (expected: %s-enriched):\n", tf, expected_dir))
  cat(sprintf("    WT♀=%.4f  WT♂=%.4f  -> %s-enriched %s\n",
              mean(wt_f), mean(wt_m), observed_dir,
              ifelse(consistent, "CONSISTENT", "*** INCONSISTENT ***")))
  cat(sprintf("    WT F vs M: p=%.2e\n", p_sex_wt))
  cat(sprintf("    KO F vs M: p=%.2e\n", p_sex_ko))
  cat(sprintf("    WT F vs KO F: p=%.2e\n", p_geno_f))
  cat(sprintf("    WT M vs KO M: p=%.2e\n", p_geno_m))

  stats_results[[tf]] <- data.frame(
    tf = tf, expected = expected_dir, observed = observed_dir,
    consistent = consistent,
    mean_WT_F = mean(wt_f), mean_WT_M = mean(wt_m),
    mean_KO_F = mean(ko_f), mean_KO_M = mean(ko_m),
    p_sex_WT = p_sex_wt, p_sex_KO = p_sex_ko,
    p_geno_F = p_geno_f, p_geno_M = p_geno_m
  )
}

stats_df <- do.call(rbind, stats_results)
write.csv(stats_df, file.path(fig_dir, "45b_scenic_validation_stats.csv"), row.names = FALSE)

# =============================================================================
# 6. Summary dot plot — all TFs, SCENIC targets
# =============================================================================
cat("\nStep 6: Summary dot plot...\n")

score_cols <- paste0(FIG5C_TFS, "_regulon")
score_cols <- score_cols[score_cols %in% colnames(eos@meta.data)]

# Rename columns for cleaner plotting
for (tf in FIG5C_TFS) {
  sc <- paste0(tf, "_regulon")
  if (sc %in% colnames(eos@meta.data)) {
    eos@meta.data[[tf]] <- eos@meta.data[[sc]]
  }
}

p_dot <- DotPlot(eos, features = FIG5C_TFS[FIG5C_TFS %in% colnames(eos@meta.data)],
                 group.by = "sex_geno",
                 cols = c("#2166AC", "#B2182B"), dot.scale = 6) +
  theme_pub(base_size = 10) +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, face = "italic", size = 9),
        axis.text.y = element_text(face = "bold")) +
  labs(title = "SCENIC Regulon Target Scores by Sex and Genotype") +
  geom_vline(xintercept = length(FEMALE_TFS) + 0.5, linetype = "dashed", color = "grey40")

save_pdf(p_dot, file.path(fig_dir, "45b_scenic_target_dotplot.pdf"),
         width = 10, height = 4.5)

# =============================================================================
# 7. Print consistency summary
# =============================================================================
cat("\n========================================\n")
cat("CONSISTENCY SUMMARY\n")
cat("========================================\n")
n_consistent <- sum(stats_df$consistent)
n_total <- nrow(stats_df)
cat(sprintf("\n  %d / %d TFs show CONSISTENT direction between\n", n_consistent, n_total))
cat("  SCENIC regulon activity (AUCell) and target gene expression (module score)\n\n")

for (i in seq_len(nrow(stats_df))) {
  r <- stats_df[i, ]
  icon <- ifelse(r$consistent, "OK", "XX")
  cat(sprintf("  [%s] %s: expected=%s, observed=%s (p=%.2e)\n",
              icon, r$tf, r$expected, r$observed, r$p_sex_WT))
}

cat(sprintf("\nTarget gene lists saved: %s\n", targets_file))
cat("\nOutputs:\n")
cat("  45b_scenic_target_violins.pdf     — Violin plots per TF\n")
cat("  45b_scenic_target_dotplot.pdf     — Summary dot plot\n")
cat("  45b_scenic_validation_stats.csv   — Statistical results\n")
cat("  45b_scenic_regulon_targets.csv    — SCENIC target gene lists\n")
cat("\n45b complete.\n")
