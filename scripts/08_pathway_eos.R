#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=2:00:00
#SBATCH --job-name=pathway_eos
#SBATCH --output=pathway_eos_%j.out

# ==============================================================================
# 08_pathway_eos.R
# GO and KEGG pathway enrichment on eosinophil DEG results
# ==============================================================================

suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Mm.eg.db)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
deg_dir <- file.path(BASE, "deg")
pw_dir  <- file.path(BASE, "pathway")
fig_dir <- file.path(BASE, "figures")
dir.create(pw_dir, showWarnings = FALSE, recursive = TRUE)
dir.create(fig_dir, showWarnings = FALSE, recursive = TRUE)

# --- Define comparisons to process -------------------------------------------
comparison_names <- c("sex_effect", "genotype_female", "genotype_male",
                      "infection_female")

# --- Helper: convert gene symbols to Entrez IDs ------------------------------
symbols_to_entrez <- function(genes) {
  mapped <- bitr(genes, fromType = "SYMBOL", toType = "ENTREZID",
                 OrgDb = org.Mm.eg.db)
  return(mapped)
}

# --- Process each comparison --------------------------------------------------
go_results_all <- list()
kegg_results_all <- list()

for (comp_name in comparison_names) {
  deg_file <- file.path(deg_dir, sprintf("deg_%s.csv", comp_name))

  if (!file.exists(deg_file)) {
    cat(sprintf("Skipping %s: DEG file not found.\n", comp_name))
    next
  }

  cat(sprintf("\n=== Processing: %s ===\n", comp_name))
  deg <- read.csv(deg_file)

  # Filter significant DEGs
  sig_deg <- deg %>% filter(padj < 0.05, abs(log2FC) > 0.5)
  cat(sprintf("  Significant DEGs: %d\n", nrow(sig_deg)))

  if (nrow(sig_deg) < 5) {
    cat("  Too few significant DEGs for enrichment, skipping.\n")
    next
  }

  # Convert gene symbols to Entrez IDs
  cat("  Converting gene symbols to Entrez IDs...\n")
  gene_map <- symbols_to_entrez(sig_deg$gene)
  cat(sprintf("  Mapped %d / %d genes to Entrez IDs\n",
              nrow(gene_map), nrow(sig_deg)))

  if (nrow(gene_map) < 3) {
    cat("  Too few mapped genes, skipping.\n")
    next
  }

  entrez_ids <- gene_map$ENTREZID

  # Also prepare background (all tested genes)
  bg_map <- symbols_to_entrez(deg$gene)
  bg_entrez <- bg_map$ENTREZID

  # --- GO enrichment (BP, MF, CC) ---
  for (ont in c("BP", "MF", "CC")) {
    cat(sprintf("  Running GO enrichment (%s)...\n", ont))
    tryCatch({
      ego <- enrichGO(
        gene          = entrez_ids,
        universe      = bg_entrez,
        OrgDb         = org.Mm.eg.db,
        ont           = ont,
        pAdjustMethod = "BH",
        pvalueCutoff  = 0.05,
        qvalueCutoff  = 0.2,
        readable      = TRUE
      )

      if (!is.null(ego) && nrow(as.data.frame(ego)) > 0) {
        go_df <- as.data.frame(ego)
        cat(sprintf("    %d significant GO %s terms\n", nrow(go_df), ont))

        write.csv(go_df,
                  file.path(pw_dir, sprintf("go_%s_%s.csv", ont, comp_name)),
                  row.names = FALSE)

        # Store for combined plot
        go_df$comparison <- comp_name
        go_df$ontology   <- ont
        go_results_all[[paste(comp_name, ont, sep = "_")]] <- go_df
      } else {
        cat(sprintf("    No significant GO %s terms.\n", ont))
      }
    }, error = function(e) {
      cat(sprintf("    GO %s error: %s\n", ont, e$message))
    })
  }

  # --- KEGG enrichment ---
  cat("  Running KEGG enrichment...\n")
  tryCatch({
    ekegg <- enrichKEGG(
      gene          = entrez_ids,
      organism      = "mmu",
      universe      = bg_entrez,
      pAdjustMethod = "BH",
      pvalueCutoff  = 0.05,
      qvalueCutoff  = 0.2
    )

    if (!is.null(ekegg) && nrow(as.data.frame(ekegg)) > 0) {
      kegg_df <- as.data.frame(ekegg)
      cat(sprintf("    %d significant KEGG pathways\n", nrow(kegg_df)))

      write.csv(kegg_df,
                file.path(pw_dir, sprintf("kegg_%s.csv", comp_name)),
                row.names = FALSE)

      kegg_df$comparison <- comp_name
      kegg_results_all[[comp_name]] <- kegg_df
    } else {
      cat("    No significant KEGG pathways.\n")
    }
  }, error = function(e) {
    cat(sprintf("    KEGG error: %s\n", e$message))
  })

  # --- Per-comparison dot plots (top 10 terms) ---
  cat("  Generating per-comparison dot plots...\n")

  # GO BP dot plot (most commonly useful)
  go_bp_key <- paste(comp_name, "BP", sep = "_")
  if (go_bp_key %in% names(go_results_all)) {
    go_bp <- go_results_all[[go_bp_key]] %>%
      arrange(p.adjust) %>%
      head(10)

    if (nrow(go_bp) > 0) {
      pdf(file.path(fig_dir, sprintf("08_go_bp_dotplot_%s.pdf", comp_name)),
          width = 10, height = 6)
      print(
        ggplot(go_bp, aes(x = Count, y = reorder(Description, -p.adjust),
                          fill = -log10(p.adjust))) +
          geom_bar(stat = "identity") +
          scale_fill_gradient(low = "#FDDBC7", high = "#B2182B") +
          labs(x = "Gene Count", y = "",
               title = sprintf("GO BP — %s (Top 10)", comp_name),
               fill = expression(-log[10](italic(p)[adj]))) +
          theme_pub() +
          theme(axis.text.y = element_text(size = 9))
      )
      dev.off()
    }
  }

  # KEGG dot plot
  if (comp_name %in% names(kegg_results_all)) {
    kegg_top <- kegg_results_all[[comp_name]] %>%
      arrange(p.adjust) %>%
      head(10)

    if (nrow(kegg_top) > 0) {
      pdf(file.path(fig_dir, sprintf("08_kegg_dotplot_%s.pdf", comp_name)),
          width = 10, height = 6)
      print(
        ggplot(kegg_top, aes(x = Count, y = reorder(Description, -p.adjust),
                             fill = -log10(p.adjust))) +
          geom_bar(stat = "identity") +
          scale_fill_gradient(low = "#FDDBC7", high = "#B2182B") +
          labs(x = "Gene Count", y = "",
               title = sprintf("KEGG — %s (Top 10)", comp_name),
               fill = expression(-log[10](italic(p)[adj]))) +
          theme_pub() +
          theme(axis.text.y = element_text(size = 9))
      )
      dev.off()
    }
  }
}

# --- Combined figure: GO BP across all comparisons ----------------------------
cat("\nGenerating combined GO BP dot plot...\n")
go_bp_combined <- bind_rows(
  go_results_all[grep("_BP$", names(go_results_all))]
)

if (nrow(go_bp_combined) > 0) {
  # Take top 5 per comparison
  go_bp_top <- go_bp_combined %>%
    group_by(comparison) %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) %>%
    ungroup()

  pdf(file.path(fig_dir, "08_go_bp_combined_dotplot.pdf"), width = 14, height = 10)
  print(
    ggplot(go_bp_top, aes(x = comparison, y = reorder(Description, -p.adjust),
                          color = -log10(p.adjust), size = Count)) +
      geom_point() +
      scale_color_gradient(low = "#FDDBC7", high = "#B2182B") +
      labs(x = "Comparison", y = "",
           title = "GO Biological Process — All Comparisons (Top 5 each)",
           color = expression(-log[10](italic(p)[adj]))) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 8))
  )
  dev.off()
}

# --- Combined KEGG figure -----------------------------------------------------
cat("Generating combined KEGG dot plot...\n")
kegg_combined <- bind_rows(kegg_results_all)

if (nrow(kegg_combined) > 0) {
  kegg_top <- kegg_combined %>%
    group_by(comparison) %>%
    arrange(p.adjust) %>%
    slice_head(n = 5) %>%
    ungroup()

  pdf(file.path(fig_dir, "08_kegg_combined_dotplot.pdf"), width = 14, height = 10)
  print(
    ggplot(kegg_top, aes(x = comparison, y = reorder(Description, -p.adjust),
                         color = -log10(p.adjust), size = Count)) +
      geom_point() +
      scale_color_gradient(low = "#FDDBC7", high = "#B2182B") +
      labs(x = "Comparison", y = "",
           title = "KEGG Pathways — All Comparisons (Top 5 each)",
           color = expression(-log[10](italic(p)[adj]))) +
      theme_pub() +
      theme(axis.text.x = element_text(angle = 45, hjust = 1),
            axis.text.y = element_text(size = 8))
  )
  dev.off()
}

cat("\nDone! 08_pathway_eos.R completed successfully.\n")
