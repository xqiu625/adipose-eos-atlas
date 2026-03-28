#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# Regenerate volcano plots from existing DEG CSVs (no model rerun needed)

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
deg_dir <- file.path(BASE, "deg")
fig_dir <- file.path(BASE, "figures")

# --- Volcano plot function ----------------------------------------------------
plot_volcano <- function(deg_df, title, filename, label_up, label_down,
                         fc_cutoff = 0.5, padj_cutoff = 0.05, top_n = 10) {
  deg_df <- deg_df %>%
    filter(!is.na(padj)) %>%
    mutate(
      sig = case_when(
        padj < padj_cutoff & log2FC >  fc_cutoff ~ "Up",
        padj < padj_cutoff & log2FC < -fc_cutoff ~ "Down",
        TRUE                                      ~ "NS"
      )
    )

  top_up <- deg_df %>% arrange(desc(log2FC)) %>% head(top_n)
  top_down <- deg_df %>% arrange(log2FC) %>% head(top_n)
  top_genes <- bind_rows(top_up, top_down) %>% distinct(gene, .keep_all = TRUE)

  x_range <- range(deg_df$log2FC, na.rm = TRUE)
  y_min <- min(-log10(deg_df$pval), na.rm = TRUE)

  p <- ggplot(deg_df, aes(x = log2FC, y = -log10(pval), color = sig)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text_repel(data = top_genes,
                    aes(label = gene), size = 2.5, max.overlaps = 30,
                    fontface = "italic", color = "black") +
    scale_color_manual(values = vol_colors) +
    geom_vline(xintercept = c(-fc_cutoff, fc_cutoff), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(padj_cutoff), linetype = "dashed", color = "grey40") +
    annotate("text", x = x_range[2] * 0.7, y = y_min, label = paste0("Up in ", label_up, " \u2192"),
             size = 3.5, fontface = "bold", color = "#E84A5F", hjust = 0.5) +
    annotate("text", x = x_range[1] * 0.7, y = y_min, label = paste0("\u2190 Up in ", label_down),
             size = 3.5, fontface = "bold", color = "#2166AC", hjust = 0.5) +
    labs(x = "log2 Fold Change",
         y = expression(-log[10](italic(p)-value)),
         title = title, color = "Significance") +
    theme_pub()

  cairo_pdf(file.path(fig_dir, filename), width = 8, height = 7)
  print(p)
  dev.off()
  cat(sprintf("  Saved: %s\n", filename))
}

# --- Replot from saved CSVs ---------------------------------------------------
comparisons <- list(
  list(csv = "deg_sex_effect.csv",         prefix = "07_volcano_sex_effect.pdf",
       desc = "WT_Female_HFD vs WT_Male_HFD (sex effect)", up = "Female", down = "Male"),
  list(csv = "deg_genotype_female.csv",    prefix = "07_volcano_genotype_female.pdf",
       desc = "KO_Female_HFD vs WT_Female_HFD (genotype effect, female)", up = "KO", down = "WT"),
  list(csv = "deg_genotype_male.csv",      prefix = "07_volcano_genotype_male.pdf",
       desc = "KO_Male_HFD vs WT_Male_HFD (genotype effect, male)", up = "KO", down = "WT"),
  list(csv = "deg_infection_female.csv",   prefix = "07_volcano_infection_female.pdf",
       desc = "Infected vs non-infected females", up = "Infected", down = "Non-infected"),

  list(csv = "deg_pseudobulk_sex_effect.csv",      prefix = "07c_volcano_pb_sex_effect.pdf",
       desc = "WT_Female_HFD vs WT_Male_HFD (sex effect)", up = "Female", down = "Male"),
  list(csv = "deg_pseudobulk_genotype_female.csv", prefix = "07c_volcano_pb_genotype_female.pdf",
       desc = "KO_Female_HFD vs WT_Female_HFD (genotype effect, female)", up = "KO", down = "WT"),
  list(csv = "deg_pseudobulk_genotype_male.csv",   prefix = "07c_volcano_pb_genotype_male.pdf",
       desc = "KO_Male_HFD vs WT_Male_HFD (genotype effect, male)", up = "KO", down = "WT"),
  list(csv = "deg_pseudobulk_infection_female.csv", prefix = "07c_volcano_pb_infection_female.pdf",
       desc = "Infected vs non-infected females", up = "Infected", down = "Non-infected")
)

for (comp in comparisons) {
  csv_path <- file.path(deg_dir, comp$csv)
  if (!file.exists(csv_path)) {
    cat(sprintf("  Skipping %s (not found)\n", comp$csv))
    next
  }
  cat(sprintf("Plotting %s...\n", comp$csv))
  deg <- read.csv(csv_path)
  # Standardize column names (MAST uses "pval", DESeq2 uses "pval" after rename)
  if (!"log2FC" %in% names(deg) && "log2FoldChange" %in% names(deg)) {
    deg <- deg %>% rename(log2FC = log2FoldChange, pval = pvalue)
  }
  plot_volcano(deg, title = comp$desc, filename = comp$prefix,
               label_up = comp$up, label_down = comp$down)
}

cat("\nDone! All volcano plots regenerated.\n")
