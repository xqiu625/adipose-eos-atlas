#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

## ADHOC: Pseudobulk DESeq2 for KO Female vs WT Female, excluding C6
## C6 excluded — too few cells (n=14)
## SLURM: --mem=64G --cpus-per-task=4 --time=1:00:00 -p batch

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
fig_dir <- file.path(BASE, "figures", "adhoc")
deg_dir <- file.path(BASE, "deg")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

cat("========================================\n")
cat("ADHOC: Pseudobulk DESeq2 — KO vs WT Female (no C6)\n")
cat("========================================\n\n")

# ── Load eos object ──────────────────────────────────────────────────────────
cat("Loading eos_subclustered.rds...\n")
eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
cat(sprintf("  %d genes x %d cells\n", nrow(eos), ncol(eos)))
cat(sprintf("  Clusters: %s\n", paste(sort(unique(Idents(eos))), collapse = ", ")))

# ── Remove C6 ───────────────────────────────────────────────────────────────
n_c6 <- sum(Idents(eos) == "C6")
eos <- subset(eos, idents = "C6", invert = TRUE)
cat(sprintf("  Removed C6: %d cells\n", n_c6))
cat(sprintf("  Remaining: %d cells in clusters %s\n",
            ncol(eos), paste(sort(unique(Idents(eos))), collapse = ", ")))

# ── Subset to KO Female + WT Female ─────────────────────────────────────────
eos_sub <- subset(eos, condition %in% c("KO_Female_HFD", "WT_Female_HFD"))
cat(sprintf("  KO_Female_HFD: %d cells\n", sum(eos_sub$condition == "KO_Female_HFD")))
cat(sprintf("  WT_Female_HFD: %d cells\n", sum(eos_sub$condition == "WT_Female_HFD")))

# ── Pseudobulk aggregation ──────────────────────────────────────────────────
cat("\nAggregating counts per sample...\n")
counts_mat <- GetAssayData(eos_sub, layer = "counts")
meta <- eos_sub@meta.data
meta$sample <- meta$orig.ident

samples <- unique(meta$sample)
pb_list <- lapply(samples, function(s) {
  cells <- which(meta$sample == s)
  Matrix::rowSums(counts_mat[, cells, drop = FALSE])
})
pb_counts <- do.call(cbind, pb_list)
colnames(pb_counts) <- samples

# Sample-level metadata
sample_group <- meta %>%
  select(sample, condition) %>%
  distinct()
col_data <- data.frame(
  sample = samples,
  group  = sample_group$condition[match(samples, sample_group$sample)],
  row.names = samples
)
col_data$group <- factor(col_data$group, levels = c("WT_Female_HFD", "KO_Female_HFD"))

cat(sprintf("  Samples: %d KO, %d WT\n",
            sum(col_data$group == "KO_Female_HFD"),
            sum(col_data$group == "WT_Female_HFD")))

# Filter lowly expressed genes
pb_counts <- pb_counts[rowSums(pb_counts) >= 10, ]
cat(sprintf("  Genes tested: %d\n", nrow(pb_counts)))

# ── DESeq2 ──────────────────────────────────────────────────────────────────
cat("\nRunning DESeq2...\n")
dds <- DESeqDataSetFromMatrix(countData = round(pb_counts),
                              colData   = col_data,
                              design    = ~ group)
dds <- DESeq(dds)

res <- results(dds, contrast = c("group", "KO_Female_HFD", "WT_Female_HFD"))
deg <- as.data.frame(res) %>%
  tibble::rownames_to_column("gene") %>%
  rename(log2FC = log2FoldChange, pval = pvalue, padj = padj) %>%
  arrange(padj)

# ── Report ──────────────────────────────────────────────────────────────────
cat("\n── DEG summary (no C6) ──\n")
for (fc in c(0.5, 0.25)) {
  for (pcut in c(0.05, 0.1)) {
    sig <- deg %>% filter(!is.na(padj), padj < pcut, abs(log2FC) > fc)
    cat(sprintf("  padj < %.2f, |log2FC| > %.2f: %d DEGs (Up in KO: %d, Up in WT: %d)\n",
                pcut, fc, nrow(sig), sum(sig$log2FC > 0), sum(sig$log2FC < 0)))
  }
}

# Top genes
cat("\n── Top 20 genes by padj ──\n")
print(head(deg, 20), row.names = FALSE)

# Check if suspect genes are still significant
cat("\n── Suspect genes status (no C6) ──\n")
suspects <- c("Cldn5", "Rasip1", "Shank3", "Cd163", "Mafb",
              "H2-Eb1", "H2-Ab1", "H2-Aa", "Cd74", "C1qa")
for (g in suspects) {
  row <- deg %>% filter(gene == g)
  if (nrow(row) > 0) {
    cat(sprintf("  %-10s  log2FC=%6.3f  padj=%s\n",
                g, row$log2FC[1],
                ifelse(is.na(row$padj[1]), "NA",
                       formatC(row$padj[1], format = "e", digits = 2))))
  } else {
    cat(sprintf("  %-10s  not tested\n", g))
  }
}

# ── Save CSV ────────────────────────────────────────────────────────────────
out_csv <- file.path(deg_dir, "deg_pseudobulk_genotype_female_noC6.csv")
write.csv(deg, out_csv, row.names = FALSE)
cat(sprintf("\nSaved: %s\n", out_csv))

# ── Compare with original (with C6) ────────────────────────────────────────
orig_file <- file.path(deg_dir, "deg_pseudobulk_genotype_female.csv")
if (file.exists(orig_file)) {
  orig <- read.csv(orig_file)
  cat("\n── Comparison: with C6 vs without C6 ──\n")

  orig_sig <- orig %>% filter(!is.na(padj), padj < 0.05, abs(log2FC) > 0.25)
  new_sig  <- deg %>% filter(!is.na(padj), padj < 0.05, abs(log2FC) > 0.25)

  cat(sprintf("  With C6:    %d DEGs\n", nrow(orig_sig)))
  cat(sprintf("  Without C6: %d DEGs\n", nrow(new_sig)))

  lost   <- setdiff(orig_sig$gene, new_sig$gene)
  gained <- setdiff(new_sig$gene, orig_sig$gene)
  kept   <- intersect(orig_sig$gene, new_sig$gene)

  if (length(kept) > 0)   cat(sprintf("  Kept:   %s\n", paste(kept, collapse = ", ")))
  if (length(lost) > 0)   cat(sprintf("  Lost:   %s\n", paste(lost, collapse = ", ")))
  if (length(gained) > 0) cat(sprintf("  Gained: %s\n", paste(gained, collapse = ", ")))
}

# ── Volcano plot ────────────────────────────────────────────────────────────
cat("\nGenerating volcano plot...\n")

deg_plot <- deg %>%
  filter(!is.na(padj)) %>%
  mutate(sig = case_when(
    padj < 0.05 & log2FC >  0.25 ~ "Up",
    padj < 0.05 & log2FC < -0.25 ~ "Down",
    TRUE ~ "NS"
  ))

top_up   <- deg_plot %>% arrange(desc(log2FC)) %>% head(10)
top_down <- deg_plot %>% arrange(log2FC) %>% head(10)
top_genes <- bind_rows(top_up, top_down) %>% distinct(gene, .keep_all = TRUE)

x_range <- range(deg_plot$log2FC, na.rm = TRUE)

p <- ggplot(deg_plot, aes(x = log2FC, y = -log10(pval), color = sig)) +
  geom_point(size = 0.8, alpha = 0.7) +
  geom_text_repel(data = top_genes, aes(label = gene),
                  size = 2.5, max.overlaps = 30,
                  fontface = "italic", color = "black") +
  scale_color_manual(values = vol_colors) +
  geom_vline(xintercept = c(-0.25, 0.25), linetype = "dashed", color = "grey40") +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
  annotate("text", x = x_range[2] * 0.7, y = 0,
           label = "Up in KO \u2192", size = 3.5, fontface = "bold",
           color = "#E84A5F", hjust = 0.5) +
  annotate("text", x = x_range[1] * 0.7, y = 0,
           label = "\u2190 Up in WT", size = 3.5, fontface = "bold",
           color = "#2166AC", hjust = 0.5) +
  labs(x = "log2 Fold Change", y = expression(-log[10](italic(p)-value)),
       title = "KO Female vs WT Female (C6 removed)",
       subtitle = "Pseudobulk DESeq2 | Eos clusters C0-C5 only",
       color = "Significance") +
  theme_pub()

cairo_pdf(file.path(fig_dir, "adhoc_volcano_genotype_female_noC6.pdf"),
          width = 8, height = 7)
print(p)
dev.off()
cat("  Saved volcano plot.\n")

cat("\nDone.\n")
