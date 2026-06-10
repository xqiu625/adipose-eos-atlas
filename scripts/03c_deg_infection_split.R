#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 03c_deg_infection_split.R
# Infection DEG comparisons split by genotype
#
# Comparisons: WT Inf vs KO Inf, WT Naive vs WT Inf, KO Naive vs KO Inf
# Method: pseudobulk DESeq2 (appropriate for n=2-3 samples per group)
#
# Input:  harmony_tier1_clustered.rds
# Output: DEG CSVs + volcano + heatmap + GO enrichment
# SLURM:  batch, 128G, 4 CPUs, 12h
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(DESeq2)
  library(ggplot2)
  library(ggrepel)
  library(dplyr)
  library(pheatmap)
  library(clusterProfiler)
  library(org.Mm.eg.db)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

# --- Paths ---
PNAS    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision"
out_dir <- file.path(PNAS, "revision/deg_infection_split")
fig_dir <- file.path(PNAS, "revision/figures/fig1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

# --- Load Tier 1 eos ---
cat("Loading Harmony Tier 1 object...\n")
seu <- readRDS(file.path(PNAS, "revision/harmony_tier1_clustered.rds"))
if (inherits(seu[["RNA"]], "Assay5")) seu[["RNA"]] <- JoinLayers(seu[["RNA"]])

# Ensure metadata
seu$sample_id <- seu$orig.ident
if (!"condition" %in% colnames(seu@meta.data)) {
  seu$genotype <- ifelse(grepl("^WT_", seu$orig.ident), "WT", "KO")
  seu$Sex <- ifelse(grepl("_Male_", seu$orig.ident), "Male", "Female")
  seu$infection_status <- ifelse(grepl("Infected", seu$orig.ident),
                                  "Infected", "NonInfected")
  seu$condition <- paste(seu$genotype, seu$Sex,
                         ifelse(seu$infection_status == "Infected",
                                "HFDInfected", "HFD"), sep = "_")
}

cat(sprintf("Tier 1 cells: %d\n", ncol(seu)))
cat("Condition distribution:\n")
print(table(seu$condition))

# --- Pseudobulk DESeq2 ---
run_pb_deseq2 <- function(seu_sub, group_col, test_level, ref_level,
                           sample_col = "orig.ident") {
  mat <- GetAssayData(seu_sub, layer = "counts")
  meta <- seu_sub@meta.data
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

# --- Volcano plot ---
make_volcano_pb <- function(deg_df, title, pdf_path, label_up, label_down) {
  deg_df <- deg_df %>%
    filter(!is.na(padj)) %>%
    mutate(sig = case_when(
      padj < 0.05 & log2FC >  0.5 ~ "Up",
      padj < 0.05 & log2FC < -0.5 ~ "Down",
      TRUE ~ "NS"
    ))

  top_up <- head(deg_df %>% arrange(desc(log2FC)), 10)
  top_dn <- head(deg_df %>% arrange(log2FC), 10)
  top_genes <- bind_rows(top_up, top_dn) %>% distinct(gene, .keep_all = TRUE)

  x_range <- range(deg_df$log2FC, na.rm = TRUE)

  p <- ggplot(deg_df, aes(x = log2FC, y = -log10(pval), color = sig)) +
    geom_point(size = 0.8, alpha = 0.7) +
    geom_text_repel(data = top_genes, aes(label = gene),
                    size = 2.5, max.overlaps = 30,
                    fontface = "italic", color = "black") +
    scale_color_manual(values = vol_colors) +
    geom_vline(xintercept = c(-0.5, 0.5), linetype = "dashed", color = "grey40") +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "grey40") +
    annotate("text", x = x_range[2] * 0.7, y = 0,
             label = paste0("Up in ", label_up), size = 3.5,
             fontface = "bold", color = "#E84A5F") +
    annotate("text", x = x_range[1] * 0.7, y = 0,
             label = paste0("Up in ", label_down), size = 3.5,
             fontface = "bold", color = "#2166AC") +
    labs(title = title, x = "log2 Fold Change",
         y = expression(-log[10](italic(p)-value))) +
    theme_pub()

  save_pdf(p, pdf_path, 8, 7)
}

# --- DEG heatmap (pseudobulk, Hu-style) ---
make_heatmap_pb <- function(deg_df, seu_obj, title, pdf_path, max_genes = 50) {
  sig <- deg_df %>%
    filter(!is.na(padj), padj < 0.05, abs(log2FC) > 0.5) %>%
    arrange(padj)

  if (nrow(sig) == 0) {
    sig <- deg_df %>%
      filter(!is.na(padj), padj < 0.1, abs(log2FC) > 0.25) %>%
      arrange(padj)
  }
  if (nrow(sig) == 0) { cat("  SKIP heatmap: no DEGs\n"); return() }

  genes_use <- head(sig$gene, max_genes)
  mat <- GetAssayData(seu_obj, layer = "data")[genes_use, , drop = FALSE]
  samples <- unique(seu_obj$orig.ident)

  pb_mat <- sapply(samples, function(s) {
    cells <- which(seu_obj$orig.ident == s)
    rowMeans(as.matrix(mat[, cells, drop = FALSE]))
  })
  pb_scaled <- t(scale(t(pb_mat)))
  pb_scaled[pb_scaled > 2.5] <- 2.5
  pb_scaled[pb_scaled < -2.5] <- -2.5

  sample_cond <- sapply(samples, function(s)
    unique(seu_obj$condition[seu_obj$orig.ident == s])[1])
  anno_col <- data.frame(Condition = sample_cond, row.names = samples)
  present <- intersect(names(cond_colors), unique(sample_cond))
  anno_colors <- list(Condition = cond_colors[present])

  # Order by direction then by condition
  up_genes <- sig$gene[sig$log2FC > 0]
  dn_genes <- sig$gene[sig$log2FC < 0]
  row_order <- c(intersect(up_genes, genes_use), intersect(dn_genes, genes_use))
  pb_scaled <- pb_scaled[row_order, order(sample_cond)]
  anno_col <- anno_col[order(sample_cond), , drop = FALSE]

  row_anno <- data.frame(
    Direction = ifelse(rownames(pb_scaled) %in% up_genes, "Up", "Down"),
    row.names = rownames(pb_scaled)
  )
  anno_colors$Direction <- c("Up" = "#B2182B", "Down" = "#2166AC")

  h <- max(6, length(genes_use) * 0.25 + 3)
  w <- max(7, ncol(pb_scaled) * 0.8 + 4)
  cairo_pdf(pdf_path, width = w, height = h)
  pheatmap(pb_scaled, cluster_rows = FALSE, cluster_cols = FALSE,
           annotation_col = anno_col, annotation_row = row_anno,
           annotation_colors = anno_colors,
           color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
           fontsize_row = 7, fontsize_col = 8, main = title, border_color = NA)
  dev.off()
  cat(sprintf("  Saved heatmap: %s\n", pdf_path))
}

# --- GO enrichment ---
run_go <- function(deg_df, prefix) {
  for (dir in c("up", "down")) {
    genes <- if (dir == "up") {
      deg_df$gene[!is.na(deg_df$padj) & deg_df$padj < 0.05 & deg_df$log2FC > 0.5]
    } else {
      deg_df$gene[!is.na(deg_df$padj) & deg_df$padj < 0.05 & deg_df$log2FC < -0.5]
    }
    if (length(genes) < 3) {
      cat(sprintf("  GO %s: %d genes, skipping\n", dir, length(genes)))
      next
    }
    ego <- tryCatch(
      enrichGO(gene = genes, OrgDb = org.Mm.eg.db, keyType = "SYMBOL",
               ont = "BP", pAdjustMethod = "BH", pvalueCutoff = 0.05),
      error = function(e) NULL)
    if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
      write.csv(as.data.frame(ego),
                file.path(out_dir, sprintf("%s_GO_BP_%s.csv", prefix, dir)),
                row.names = FALSE)
      p <- dotplot(ego, showCategory = 15) + ggtitle(sprintf("%s — %s", prefix, dir))
      save_pdf(p, file.path(out_dir, sprintf("%s_GO_%s_dotplot.pdf", prefix, dir)), 10, 8)
      cat(sprintf("  GO %s: %d terms\n", dir, nrow(as.data.frame(ego))))
    }
  }
}

# =============================================================================
# Define 3 new infection comparisons (Rev 1 Comment 11)
# =============================================================================
comparisons <- list(
  list(
    name = "WT_Inf_vs_WT_Naive",
    desc = "WT Infected vs WT Naive (females only)",
    test_cond = "WT_Female_HFDInfected",
    ref_cond  = "WT_Female_HFD",
    label_up  = "WT Infected",
    label_down = "WT Naive"
  ),
  list(
    name = "KO_Inf_vs_KO_Naive",
    desc = "KO Infected vs KO Naive (females only)",
    test_cond = "KO_Female_HFDInfected",
    ref_cond  = "KO_Female_HFD",
    label_up  = "KO Infected",
    label_down = "KO Naive"
  ),
  list(
    name = "WT_Inf_vs_KO_Inf",
    desc = "WT Infected vs KO Infected (females only)",
    test_cond = "WT_Female_HFDInfected",
    ref_cond  = "KO_Female_HFDInfected",
    label_up  = "WT Infected",
    label_down = "KO Infected"
  )
)

# =============================================================================
# Run analyses
# =============================================================================
for (comp in comparisons) {
  cat(sprintf("\n========================================\n"))
  cat(sprintf("Comparison: %s\n", comp$desc))
  cat(sprintf("========================================\n"))

  test_cells <- colnames(seu)[seu$condition == comp$test_cond]
  ref_cells  <- colnames(seu)[seu$condition == comp$ref_cond]
  cat(sprintf("  Test: %d cells (%s)\n", length(test_cells), comp$test_cond))
  cat(sprintf("  Ref:  %d cells (%s)\n", length(ref_cells), comp$ref_cond))

  seu_sub <- subset(seu, condition %in% c(comp$test_cond, comp$ref_cond))

  # --- Pseudobulk DESeq2 ---
  cat("\n  --- Pseudobulk DESeq2 ---\n")
  tryCatch({
    pb_res <- run_pb_deseq2(seu_sub, "condition", comp$test_cond, comp$ref_cond)
    n_sig <- sum(!is.na(pb_res$padj) & pb_res$padj < 0.05 & abs(pb_res$log2FC) > 0.5)
    cat(sprintf("  Pseudobulk DEGs (padj<0.05, |log2FC|>0.5): %d\n", n_sig))

    write.csv(pb_res, file.path(out_dir, sprintf("deg_pb_%s.csv", comp$name)),
              row.names = FALSE)

    make_volcano_pb(pb_res, comp$desc,
                    file.path(fig_dir, sprintf("fig1_inf_%s_volcano_pb.pdf", comp$name)),
                    comp$label_up, comp$label_down)

    make_heatmap_pb(pb_res, seu_sub, paste("Pseudobulk —", comp$desc),
                    file.path(fig_dir, sprintf("fig1_inf_%s_heatmap_pb.pdf", comp$name)))

    run_go(pb_res, sprintf("pb_%s", comp$name))

  }, error = function(e) {
    cat(sprintf("  ERROR pseudobulk: %s\n", e$message))
  })

}

# =============================================================================
# Also: Multi-factor model (genotype + infection as covariates)
# =============================================================================
cat("\n========================================\n")
cat("Multi-factor: ~ genotype + infection_status (females only)\n")
cat("========================================\n")

fem <- subset(seu, Sex == "Female")
fem$infection_status <- ifelse(grepl("Infected", fem$condition), "Infected", "Naive")
fem$genotype <- ifelse(grepl("^WT_", fem$orig.ident), "WT", "KO")

cat(sprintf("  Female cells: %d\n", ncol(fem)))
cat("  Design matrix:\n")
print(table(fem$genotype, fem$infection_status))

tryCatch({
  mat <- GetAssayData(fem, layer = "counts")
  samples <- unique(fem$orig.ident)

  pb_list <- lapply(samples, function(s) {
    cells <- which(fem$orig.ident == s)
    Matrix::rowSums(mat[, cells, drop = FALSE])
  })
  pb_counts <- do.call(cbind, pb_list)
  colnames(pb_counts) <- samples
  pb_counts <- pb_counts[rowSums(pb_counts) >= 10, ]

  col_data <- data.frame(
    sample = samples,
    genotype = sapply(samples, function(s) unique(fem$genotype[fem$orig.ident == s])[1]),
    infection = sapply(samples, function(s) unique(fem$infection_status[fem$orig.ident == s])[1]),
    row.names = samples
  )
  col_data$genotype <- factor(col_data$genotype, levels = c("WT", "KO"))
  col_data$infection <- factor(col_data$infection, levels = c("Naive", "Infected"))

  cat(sprintf("  Pseudobulk: %d genes, %d samples\n", nrow(pb_counts), ncol(pb_counts)))
  cat("  Sample design:\n")
  print(table(col_data$genotype, col_data$infection))

  dds <- DESeqDataSetFromMatrix(countData = round(pb_counts),
                                colData = col_data,
                                design = ~ genotype + infection)
  dds <- DESeq(dds)

  # Infection effect (controlling for genotype)
  res_inf <- as.data.frame(results(dds, contrast = c("infection", "Infected", "Naive")))
  res_inf$gene <- rownames(res_inf)
  res_inf <- res_inf %>% rename(log2FC = log2FoldChange, pval = pvalue) %>% arrange(padj)
  n_sig <- sum(!is.na(res_inf$padj) & res_inf$padj < 0.05 & abs(res_inf$log2FC) > 0.5)
  cat(sprintf("  Infection effect (controlling genotype): %d DEGs\n", n_sig))
  write.csv(res_inf, file.path(out_dir, "deg_pb_multifactor_infection.csv"), row.names = FALSE)

  make_volcano_pb(res_inf, "Infection effect (genotype-adjusted)",
                  file.path(fig_dir, "fig1_inf_multifactor_infection_volcano.pdf"),
                  "Infected", "Naive")

  # Genotype effect (controlling for infection)
  res_geno <- as.data.frame(results(dds, contrast = c("genotype", "KO", "WT")))
  res_geno$gene <- rownames(res_geno)
  res_geno <- res_geno %>% rename(log2FC = log2FoldChange, pval = pvalue) %>% arrange(padj)
  n_sig <- sum(!is.na(res_geno$padj) & res_geno$padj < 0.05 & abs(res_geno$log2FC) > 0.5)
  cat(sprintf("  Genotype effect (controlling infection): %d DEGs\n", n_sig))
  write.csv(res_geno, file.path(out_dir, "deg_pb_multifactor_genotype.csv"), row.names = FALSE)

  make_volcano_pb(res_geno, "Genotype effect (infection-adjusted)",
                  file.path(fig_dir, "fig1_inf_multifactor_genotype_volcano.pdf"),
                  "KO", "WT")

  # GO on infection effect
  run_go(res_inf, "multifactor_infection")

}, error = function(e) {
  cat(sprintf("  ERROR multi-factor: %s\n", e$message))
})

cat(sprintf("\nDone! Results saved to:\n  %s\n  %s\n", out_dir, fig_dir))
