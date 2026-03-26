#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))
library(Seurat); library(Matrix); library(ggplot2); library(dplyr)

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"

cat("=== Post-DecontX annotated object ===\n")
obj <- readRDS(file.path(BASE, "all_cells_annotated.rds"))

for (g in c("Igha", "Retnla")) {
  cat(sprintf("\n=== %s ===\n", g))
  if (g %in% rownames(obj)) {
    vals <- GetAssayData(obj, layer = "counts")[g, ]
    cat(sprintf("Total counts: %.1f, Cells expressing: %.1f%%\n",
                sum(vals), 100 * mean(vals > 0)))
    ct <- obj$celltype
    for (t in sort(unique(ct))) {
      cells <- vals[ct == t]
      cat(sprintf("  %-20s n=%5d mean=%.3f %%expr=%.1f%%\n",
                  t, length(cells), mean(cells), 100 * mean(cells > 0)))
    }
  } else {
    cat(g, "NOT in gene list\n")
  }
}

cat("\n=== Before vs After DecontX (sample 2, WT_Female_HFD_1) ===\n")
decon <- readRDS(file.path(BASE, "decontx/decontx_sample_2.rds"))
orig <- Read10X("/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data/RR0104/outs/per_sample_outs/2/count/sample_filtered_feature_bc_matrix")
corr <- decon$counts
common <- intersect(colnames(orig), colnames(corr))

for (g in c("Igha", "Retnla")) {
  if (g %in% rownames(orig) && g %in% rownames(corr)) {
    b <- mean(orig[g, common])
    a <- mean(corr[g, common])
    cat(sprintf("%s: before=%.4f after=%.4f reduction=%.1f%%\n",
                g, b, a, 100 * (b - a) / max(b, 1e-9)))
  } else {
    cat(g, "not found in one of the matrices\n")
  }
}

# ==============================================================================
# Gene expression bar plots by cell type (SiglecF, Igha, Retnla)
# ==============================================================================
fig_dir <- file.path(BASE, "figures")
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

immune_types <- c("Eosinophils", "Macrophages", "Dendritic cells",
                  "B cells", "NK cells", "Neutrophils", "Mast cells")

plot_gene_by_celltype <- function(obj, gene_name, display_name, filename) {
  cat(sprintf("\n=== %s expression by cell type ===\n", display_name))

  if (!gene_name %in% rownames(obj)) {
    cat(sprintf("  %s NOT in gene list\n", gene_name))
    return(invisible(NULL))
  }

  vals <- GetAssayData(obj, layer = "counts")[gene_name, ]
  ct <- obj$celltype

  df <- data.frame(celltype = ct, expr = vals) %>%
    group_by(celltype) %>%
    summarise(
      n = n(),
      pct_expr = 100 * mean(expr > 0),
      mean_expr = mean(expr),
      .groups = "drop"
    ) %>%
    arrange(desc(pct_expr))

  cat("  Per cell type:\n")
  for (i in seq_len(nrow(df))) {
    r <- df[i, ]
    cat(sprintf("  %-25s n=%5d  %%expr=%5.1f%%  mean=%.3f\n",
                r$celltype, r$n, r$pct_expr, r$mean_expr))
  }

  df$category <- ifelse(df$celltype %in% immune_types,
                        "Immune", "Stromal/Epithelial")
  df$celltype <- factor(df$celltype, levels = rev(df$celltype))

  p <- ggplot(df, aes(x = pct_expr, y = celltype, fill = category)) +
    geom_col(width = 0.7) +
    scale_fill_manual(values = c("Immune" = "#E02020",
                                 "Stromal/Epithelial" = "#1E88C7"),
                      name = "Cell Category") +
    labs(x = sprintf("%% Cells Expressing %s", display_name),
         y = "Cell Type",
         title = sprintf("Percentage of Cells Expressing %s by Cell Type",
                         display_name)) +
    theme_classic(base_size = 14) +
    theme(
      plot.title = element_text(face = "bold", size = 16),
      axis.text.y = element_text(size = 12),
      axis.text.x = element_text(size = 11),
      axis.title = element_text(size = 13, face = "bold"),
      legend.position = c(0.75, 0.3),
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_line(color = "grey90")
    )

  pdf(file.path(fig_dir, filename), width = 10, height = 7)
  print(p)
  dev.off()
  cat(sprintf("  Saved: %s\n", filename))
}

# --- Bar plot faceted by condition ---
plot_gene_by_celltype_condition <- function(obj, gene_name, display_name, filename) {
  cat(sprintf("\n=== %s expression by cell type x condition ===\n", display_name))

  if (!gene_name %in% rownames(obj)) {
    cat(sprintf("  %s NOT in gene list\n", gene_name))
    return(invisible(NULL))
  }

  vals <- GetAssayData(obj, layer = "counts")[gene_name, ]

  df <- data.frame(celltype = obj$celltype, condition = obj$condition,
                   expr = vals) %>%
    group_by(condition, celltype) %>%
    summarise(
      n = n(),
      pct_expr = 100 * mean(expr > 0),
      .groups = "drop"
    )

  # Use overall ordering (by mean pct across conditions)
  ct_order <- df %>%
    group_by(celltype) %>%
    summarise(mean_pct = mean(pct_expr), .groups = "drop") %>%
    arrange(desc(mean_pct)) %>%
    pull(celltype)

  df$category <- ifelse(df$celltype %in% immune_types,
                        "Immune", "Stromal/Epithelial")
  df$celltype <- factor(df$celltype, levels = rev(ct_order))

  p <- ggplot(df, aes(x = pct_expr, y = celltype, fill = category)) +
    geom_col(width = 0.7) +
    facet_wrap(~ condition, nrow = 1) +
    scale_fill_manual(values = c("Immune" = "#E02020",
                                 "Stromal/Epithelial" = "#1E88C7"),
                      name = "Cell Category") +
    labs(x = sprintf("%% Cells Expressing %s", display_name),
         y = "Cell Type",
         title = sprintf("%s Expression by Cell Type and Condition", display_name)) +
    theme_classic(base_size = 11) +
    theme(
      plot.title = element_text(face = "bold", size = 14),
      axis.text.y = element_text(size = 9),
      axis.text.x = element_text(size = 8),
      axis.title = element_text(size = 11, face = "bold"),
      strip.text = element_text(face = "bold", size = 8),
      legend.position = "bottom",
      legend.title = element_text(face = "bold"),
      panel.grid.major.x = element_line(color = "grey90")
    )

  pdf(file.path(fig_dir, filename), width = 18, height = 7)
  print(p)
  dev.off()
  cat(sprintf("  Saved: %s\n", filename))
}

# Find SiglecF gene name
siglec_gene <- grep("^Siglec[fF]", rownames(obj), value = TRUE)
if (length(siglec_gene) == 0) siglec_gene <- intersect(c("Siglecf", "SiglecF"), rownames(obj))
cat(sprintf("SiglecF gene name(s) found: %s\n", paste(siglec_gene, collapse = ", ")))

if (length(siglec_gene) >= 1) {
  plot_gene_by_celltype(obj, siglec_gene[1], "SiglecF", "adhoc_siglecf_by_celltype.pdf")
}
plot_gene_by_celltype(obj, "Igha", "Igha", "adhoc_igha_by_celltype.pdf")
plot_gene_by_celltype(obj, "Retnla", "Retnla", "adhoc_retnla_by_celltype.pdf")

# Per-condition plots
if (length(siglec_gene) >= 1) {
  plot_gene_by_celltype_condition(obj, siglec_gene[1], "SiglecF", "adhoc_siglecf_by_celltype_condition.pdf")
}
plot_gene_by_celltype_condition(obj, "Igha", "Igha", "adhoc_igha_by_celltype_condition.pdf")
plot_gene_by_celltype_condition(obj, "Retnla", "Retnla", "adhoc_retnla_by_celltype_condition.pdf")

cat("\nDone.\n")
