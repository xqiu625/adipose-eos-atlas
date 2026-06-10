#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 07c_deg_pseudobulk.R
# Pseudobulk DESeq2 DEG analysis for eosinophils (4 comparisons)
# Aggregates counts per sample, runs DESeq2, compares with MAST results
# SLURM: --mem=64G --cpus-per-task=4 --time=2:00:00 -p batch
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(ggplot2)
  library(dplyr)
  library(ggrepel)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
deg_dir <- file.path(BASE, "deg")
fig_dir <- file.path(BASE, "figures")
dir.create(deg_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# --- Helper: pseudobulk DESeq2 ------------------------------------------------
run_pseudobulk_deseq2 <- function(seu, group_col, test_level, ref_level,
                                   sample_col = "orig.ident") {
  # Subset to the two groups
  cells_use <- colnames(seu)[seu@meta.data[[group_col]] %in% c(test_level, ref_level)]
  seu_sub <- subset(seu, cells = cells_use)

  cat(sprintf("    Cells: %d test (%s), %d ref (%s)\n",
              sum(seu_sub@meta.data[[group_col]] == test_level), test_level,
              sum(seu_sub@meta.data[[group_col]] == ref_level), ref_level))

  # Aggregate raw counts per sample
  counts_mat <- GetAssayData(seu_sub, layer = "counts")
  meta <- seu_sub@meta.data
  meta$group <- meta[[group_col]]
  meta$sample <- meta[[sample_col]]

  samples <- unique(meta$sample)
  pb_list <- lapply(samples, function(s) {
    cells <- which(meta$sample == s)
    Matrix::rowSums(counts_mat[, cells, drop = FALSE])
  })
  pb_counts <- do.call(cbind, pb_list)
  colnames(pb_counts) <- samples

  # Sample-level metadata
  sample_group <- meta %>%
    select(sample, group) %>%
    distinct()
  col_data <- data.frame(
    sample = samples,
    group  = sample_group$group[match(samples, sample_group$sample)],
    row.names = samples
  )
  col_data$group <- factor(col_data$group, levels = c(ref_level, test_level))

  cat(sprintf("    Samples: %d test, %d ref\n",
              sum(col_data$group == test_level),
              sum(col_data$group == ref_level)))

  # Filter lowly expressed genes (keep genes with >= 10 counts total)
  pb_counts <- pb_counts[rowSums(pb_counts) >= 10, ]
  cat(sprintf("    Genes tested: %d\n", nrow(pb_counts)))

  # DESeq2
  cat("    Running DESeq2...\n")
  dds <- DESeqDataSetFromMatrix(countData = round(pb_counts),
                                colData   = col_data,
                                design    = ~ group)
  dds <- DESeq(dds)

  res <- results(dds, contrast = c("group", test_level, ref_level))
  res_df <- as.data.frame(res) %>%
    tibble::rownames_to_column("gene") %>%
    rename(log2FC = log2FoldChange, pval = pvalue, padj = padj) %>%
    arrange(padj)

  cat(sprintf("    Total genes with results: %d\n", sum(!is.na(res_df$padj))))
  return(res_df)
}

# --- Helper: volcano plot -----------------------------------------------------
plot_volcano <- function(deg_df, title, filename, label_up, label_down,
                         fc_cutoff = 0.25, padj_cutoff = 0.05, top_n = 10) {
  deg_df <- deg_df %>%
    filter(!is.na(padj)) %>%
    mutate(
      sig = case_when(
        padj < padj_cutoff & log2FC >  fc_cutoff ~ "Up",
        padj < padj_cutoff & log2FC < -fc_cutoff ~ "Down",
        TRUE                                      ~ "NS"
      )
    )

  # Always label top N genes in each direction by log2FC (regardless of sig)
  top_up <- deg_df %>% arrange(desc(log2FC)) %>% head(top_n)
  top_down <- deg_df %>% arrange(log2FC) %>% head(top_n)
  top_genes <- bind_rows(top_up, top_down) %>% distinct(gene, .keep_all = TRUE)

  # X-axis range for annotation placement
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
}

# --- Load eosinophil object ---------------------------------------------------
cat("Loading eos_subclustered.rds (decontX-cleaned, clustered)...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
cat(sprintf("  %d eosinophils, clusters: %s\n", ncol(eos),
            paste(sort(unique(Idents(eos))), collapse = ", ")))

# Remove C6 — too few cells (n=14)
n_c6 <- sum(Idents(eos) == "C6")
eos <- subset(eos, idents = "C6", invert = TRUE)
cat(sprintf("  Removed C6: %d cells. Remaining: %d eos (C0-C5)\n", n_c6, ncol(eos)))

# --- Define comparisons -------------------------------------------------------
comparisons <- list(
  sex_effect = list(
    name      = "sex_effect",
    desc      = "WT_Female_HFD vs WT_Male_HFD (sex effect)",
    subset_conditions = c("WT_Female_HFD", "WT_Male_HFD"),
    group_col = "condition",
    test      = "WT_Female_HFD",
    ref       = "WT_Male_HFD",
    label_up  = "Female",
    label_down = "Male"
  ),
  genotype_female = list(
    name      = "genotype_female",
    desc      = "KO_Female_HFD vs WT_Female_HFD (genotype effect, female)",
    subset_conditions = c("KO_Female_HFD", "WT_Female_HFD"),
    group_col = "condition",
    test      = "KO_Female_HFD",
    ref       = "WT_Female_HFD",
    label_up  = "KO",
    label_down = "WT"
  ),
  genotype_male = list(
    name      = "genotype_male",
    desc      = "KO_Male_HFD vs WT_Male_HFD (genotype effect, male)",
    subset_conditions = c("KO_Male_HFD", "WT_Male_HFD"),
    group_col = "condition",
    test      = "KO_Male_HFD",
    ref       = "WT_Male_HFD",
    label_up  = "KO",
    label_down = "WT"
  ),
  infection_female = list(
    name      = "infection_female",
    desc      = "Infected vs non-infected females",
    subset_conditions = c("WT_Female_HFDInfected", "KO_Female_HFDInfected",
                          "WT_Female_HFD", "KO_Female_HFD"),
    group_col = "infection_status",
    test      = "Infected",
    ref       = "NonInfected",
    label_up  = "Infected",
    label_down = "Non-infected"
  )
)

# --- Run pseudobulk DESeq2 for each comparison --------------------------------
cat("\n============================================================\n")
cat("Pseudobulk DESeq2 DEG Analysis\n")
cat("============================================================\n")

for (comp in comparisons) {
  cat(sprintf("\n=== Comparison: %s ===\n", comp$desc))

  # Subset to relevant conditions
  eos_sub <- subset(eos, condition %in% comp$subset_conditions)

  # For infection comparison, create a binary group variable
  if (comp$name == "infection_female") {
    eos_sub$infection_status <- ifelse(
      grepl("Infected", eos_sub$condition), "Infected", "NonInfected"
    )
  }

  cat(sprintf("  Total cells in comparison: %d\n", ncol(eos_sub)))

  tryCatch({
    deg_res <- run_pseudobulk_deseq2(eos_sub,
                                      group_col = comp$group_col,
                                      test_level = comp$test,
                                      ref_level  = comp$ref)

    # Report at multiple thresholds
    for (fc in c(0.5, 0.25)) {
      for (pcut in c(0.05, 0.1)) {
        sig_res <- deg_res %>%
          filter(!is.na(padj), padj < pcut, abs(log2FC) > fc)
        cat(sprintf("  DEGs (padj < %.2f, |log2FC| > %.2f): %d (Up: %d, Down: %d)\n",
                    pcut, fc, nrow(sig_res),
                    sum(sig_res$log2FC > 0), sum(sig_res$log2FC < 0)))
      }
    }

    # Save full results
    out_file <- file.path(deg_dir, sprintf("deg_pseudobulk_%s.csv", comp$name))
    write.csv(deg_res, out_file, row.names = FALSE)
    cat(sprintf("  Saved: %s\n", out_file))

    # Volcano plots at standard threshold
    plot_volcano(deg_res,
                 title      = comp$desc,
                 filename   = sprintf("07c_volcano_pb_%s.pdf", comp$name),
                 label_up   = comp$label_up,
                 label_down = comp$label_down)
    cat("  Volcano plot saved.\n")

  }, error = function(e) {
    cat(sprintf("  ERROR in comparison %s: %s\n", comp$name, e$message))
  })
}

# ==============================================================================
# Multi-factor DESeq2: infection females — ~ genotype + infection_status
# Extracts both the genotype effect (controlling for infection) and
# the infection effect (controlling for genotype)
# ==============================================================================
cat("\n============================================================\n")
cat("Multi-factor DESeq2: females ~ genotype + infection_status\n")
cat("============================================================\n")

eos_fem <- subset(eos, condition %in% c("WT_Female_HFD", "KO_Female_HFD",
                                         "WT_Female_HFDInfected", "KO_Female_HFDInfected"))
eos_fem$infection_status <- ifelse(grepl("Infected", eos_fem$condition), "Infected", "NonInfected")
eos_fem$genotype <- ifelse(grepl("^KO", eos_fem$condition), "KO", "WT")

cat(sprintf("  Total female cells: %d\n", ncol(eos_fem)))
cat(sprintf("  WT Naive: %d, WT Infected: %d, KO Naive: %d, KO Infected: %d\n",
            sum(eos_fem$condition == "WT_Female_HFD"),
            sum(eos_fem$condition == "WT_Female_HFDInfected"),
            sum(eos_fem$condition == "KO_Female_HFD"),
            sum(eos_fem$condition == "KO_Female_HFDInfected")))

# Aggregate per sample
counts_mat <- GetAssayData(eos_fem, layer = "counts")
meta_fem <- eos_fem@meta.data
meta_fem$sample <- meta_fem$orig.ident
samples <- unique(meta_fem$sample)

pb_list <- lapply(samples, function(s) {
  cells <- which(meta_fem$sample == s)
  Matrix::rowSums(counts_mat[, cells, drop = FALSE])
})
pb_counts <- do.call(cbind, pb_list)
colnames(pb_counts) <- samples

# Sample metadata
col_data <- meta_fem %>%
  select(sample, genotype, infection_status) %>%
  distinct()
rownames(col_data) <- col_data$sample
col_data <- col_data[samples, ]
col_data$genotype <- factor(col_data$genotype, levels = c("WT", "KO"))
col_data$infection_status <- factor(col_data$infection_status, levels = c("NonInfected", "Infected"))

cat(sprintf("  Samples per group:\n"))
print(table(col_data$genotype, col_data$infection_status))

pb_counts <- pb_counts[rowSums(pb_counts) >= 10, ]
cat(sprintf("  Genes tested: %d\n", nrow(pb_counts)))

# Fit multi-factor model
cat("  Running DESeq2 (~ genotype + infection_status)...\n")
dds_mf <- DESeqDataSetFromMatrix(countData = round(pb_counts),
                                  colData   = col_data,
                                  design    = ~ genotype + infection_status)
dds_mf <- DESeq(dds_mf)

# --- Extract genotype effect (controlling for infection) ---
cat("\n--- Multi-factor: Genotype effect (KO vs WT, controlling for infection) ---\n")
res_geno <- results(dds_mf, contrast = c("genotype", "KO", "WT"))
deg_geno <- as.data.frame(res_geno) %>%
  tibble::rownames_to_column("gene") %>%
  rename(log2FC = log2FoldChange, pval = pvalue, padj = padj) %>%
  arrange(padj)

for (fc in c(0.5, 0.25)) {
  for (pcut in c(0.05, 0.1)) {
    sig <- deg_geno %>% filter(!is.na(padj), padj < pcut, abs(log2FC) > fc)
    cat(sprintf("  DEGs (padj < %.2f, |log2FC| > %.2f): %d (Up KO: %d, Up WT: %d)\n",
                pcut, fc, nrow(sig), sum(sig$log2FC > 0), sum(sig$log2FC < 0)))
  }
}
write.csv(deg_geno, file.path(deg_dir, "deg_pseudobulk_multifactor_genotype.csv"), row.names = FALSE)
cat("  Saved: deg_pseudobulk_multifactor_genotype.csv\n")

plot_volcano(deg_geno,
             title    = "Multi-factor: Genotype (KO vs WT, ctrl for infection)",
             filename = "07c_volcano_pb_multifactor_genotype.pdf",
             label_up = "KO", label_down = "WT")

# --- Extract infection effect (controlling for genotype) ---
cat("\n--- Multi-factor: Infection effect (Infected vs Naive, controlling for genotype) ---\n")
res_inf <- results(dds_mf, contrast = c("infection_status", "Infected", "NonInfected"))
deg_inf <- as.data.frame(res_inf) %>%
  tibble::rownames_to_column("gene") %>%
  rename(log2FC = log2FoldChange, pval = pvalue, padj = padj) %>%
  arrange(padj)

for (fc in c(0.5, 0.25)) {
  for (pcut in c(0.05, 0.1)) {
    sig <- deg_inf %>% filter(!is.na(padj), padj < pcut, abs(log2FC) > fc)
    cat(sprintf("  DEGs (padj < %.2f, |log2FC| > %.2f): %d (Up Inf: %d, Up Naive: %d)\n",
                pcut, fc, nrow(sig), sum(sig$log2FC > 0), sum(sig$log2FC < 0)))
  }
}
write.csv(deg_inf, file.path(deg_dir, "deg_pseudobulk_multifactor_infection.csv"), row.names = FALSE)
cat("  Saved: deg_pseudobulk_multifactor_infection.csv\n")

plot_volcano(deg_inf,
             title    = "Multi-factor: Infection (Inf vs Naive, ctrl for genotype)",
             filename = "07c_volcano_pb_multifactor_infection.pdf",
             label_up = "Infected", label_down = "Naive")

rm(eos_fem, dds_mf); gc()

# --- Compare with MAST results -----------------------------------------------
cat("\n============================================================\n")
cat("Comparison: Pseudobulk DESeq2 vs MAST mixed-effects\n")
cat("============================================================\n")

for (comp in comparisons) {
  mast_file <- file.path(deg_dir, sprintf("deg_%s.csv", comp$name))
  pb_file   <- file.path(deg_dir, sprintf("deg_pseudobulk_%s.csv", comp$name))

  if (file.exists(mast_file) && file.exists(pb_file)) {
    mast <- read.csv(mast_file)
    pb   <- read.csv(pb_file)

    mast_sig <- mast %>% filter(padj < 0.05, abs(log2FC) > 0.5)
    pb_sig   <- pb %>% filter(!is.na(padj), padj < 0.05, abs(log2FC) > 0.5)

    overlap <- intersect(mast_sig$gene, pb_sig$gene)

    cat(sprintf("\n%s:\n", comp$desc))
    cat(sprintf("  MAST sig DEGs:       %d\n", nrow(mast_sig)))
    cat(sprintf("  Pseudobulk sig DEGs: %d\n", nrow(pb_sig)))
    cat(sprintf("  Overlap:             %d\n", length(overlap)))
    if (length(overlap) > 0) {
      cat(sprintf("  Shared genes: %s\n", paste(overlap, collapse = ", ")))
    }
  }
}

cat("\nDone! 07c_deg_pseudobulk.R completed successfully.\n")
