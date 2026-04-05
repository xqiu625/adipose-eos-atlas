#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 26_gsea_custom_butterfly.R
# Custom GSEA butterfly plots with Rebecca's selected pathways
# SLURM: batch, 32G, 2 CPUs, 30min
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(stringr)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

PNAS    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision"
out_dir <- file.path(PNAS, "revision/figures/fig1")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# 5.5cm x 7.5cm
w_in <- 5.5 / 2.54
h_in <- 7.5 / 2.54

# --- Helper: custom butterfly from selected terms ---
make_custom_butterfly <- function(gsea_csv, selected_left, selected_right,
                                  label_left, label_right, prefix, out_dir) {
  gsea_df <- read.csv(gsea_csv)

  frames <- list()

  # Left side (NES < 0 = UP in label_left)
  if (length(selected_left) > 0) {
    matches <- gsea_df[sapply(gsea_df$Description, function(d)
      any(sapply(selected_left, function(s) grepl(s, d, ignore.case = TRUE)))), ]
    matches <- matches[matches$NES < 0 | nrow(matches) == 0, ]
    if (nrow(matches) > 0) {
      # Order by match order in selected_left
      match_order <- sapply(selected_left, function(s)
        which(sapply(matches$Description, function(d) grepl(s, d, ignore.case = TRUE)))[1])
      match_order <- match_order[!is.na(match_order)]
      matches <- matches[match_order, ]
      frames[["left"]] <- data.frame(
        Description = as.character(matches$Description),
        neg_log10p = -log10(matches$pvalue),
        plot_x = -(-log10(matches$pvalue)), side = "left",
        stringsAsFactors = FALSE)
    }
  }

  # Right side (NES > 0 = UP in label_right)
  if (length(selected_right) == 1 && selected_right[1] == "TOP5") {
    # Auto-select top 5 by p-value
    right_df <- gsea_df[gsea_df$NES > 0, ]
    if (nrow(right_df) > 0) {
      right_df <- head(right_df[order(right_df$pvalue), ], 5)
      frames[["right"]] <- data.frame(
        Description = as.character(right_df$Description),
        neg_log10p = -log10(right_df$pvalue),
        plot_x = -log10(right_df$pvalue), side = "right",
        stringsAsFactors = FALSE)
    }
  } else if (length(selected_right) > 0) {
    matches <- gsea_df[sapply(gsea_df$Description, function(d)
      any(sapply(selected_right, function(s) grepl(s, d, ignore.case = TRUE)))), ]
    matches <- matches[matches$NES > 0 | nrow(matches) == 0, ]
    if (nrow(matches) > 0) {
      match_order <- sapply(selected_right, function(s)
        which(sapply(matches$Description, function(d) grepl(s, d, ignore.case = TRUE)))[1])
      match_order <- match_order[!is.na(match_order)]
      matches <- matches[match_order, ]
      frames[["right"]] <- data.frame(
        Description = as.character(matches$Description),
        neg_log10p = -log10(matches$pvalue),
        plot_x = -log10(matches$pvalue), side = "right",
        stringsAsFactors = FALSE)
    }
  }

  # Left side: also support "TOP5" auto-select
  if (length(selected_left) == 1 && selected_left[1] == "TOP5") {
    left_df <- gsea_df[gsea_df$NES < 0, ]
    if (nrow(left_df) > 0) {
      left_df <- head(left_df[order(left_df$pvalue), ], 5)
      frames[["left"]] <- data.frame(
        Description = as.character(left_df$Description),
        neg_log10p = -log10(left_df$pvalue),
        plot_x = -(-log10(left_df$pvalue)), side = "left",
        stringsAsFactors = FALSE)
    }
  }

  if (length(frames) == 0) {
    cat(sprintf("  WARNING: no matching terms for %s\n", prefix))
    return()
  }

  plot_df <- do.call(rbind, frames)
  plot_df$label <- str_wrap(plot_df$Description, width = 25)

  left_desc  <- if ("left" %in% names(frames)) rev(frames[["left"]]$Description) else character(0)
  right_desc <- if ("right" %in% names(frames)) frames[["right"]]$Description else character(0)
  plot_df$Description <- factor(plot_df$Description, levels = c(left_desc, right_desc))
  plot_df$text_x     <- ifelse(plot_df$side == "left", -0.15, 0.15)
  plot_df$text_hjust <- ifelse(plot_df$side == "left", 1, 0)
  x_max <- max(abs(plot_df$plot_x)) * 1.15
  side_colors <- c("left" = "#2166AC", "right" = "#B2182B")

  p <- ggplot(plot_df, aes(x = plot_x, y = Description, fill = side)) +
    geom_col(width = 0.85, alpha = 0.85) +
    geom_text(aes(x = text_x, label = label, hjust = text_hjust),
              size = 1.8, lineheight = 0.8, color = "black", show.legend = FALSE) +
    scale_fill_manual(values = side_colors, guide = "none") +
    scale_x_continuous(labels = abs, limits = c(-x_max, x_max)) +
    geom_vline(xintercept = 0, color = "black", linewidth = 0.3) +
    labs(x = expression(-Log[10] ~ "(P-value)"), y = "", title = "GSEA GO BP") +
    annotate("text", x = -x_max * 0.5, y = nrow(plot_df) + 0.7,
             label = paste("UP in", label_left), hjust = 0.5,
             size = 2, fontface = "bold", color = "#2166AC") +
    annotate("text", x = x_max * 0.5, y = nrow(plot_df) + 0.7,
             label = paste("UP in", label_right), hjust = 0.5,
             size = 2, fontface = "bold", color = "#B2182B") +
    theme_minimal(base_size = 7) +
    theme(
      axis.text.y = element_blank(), axis.ticks.y = element_blank(),
      axis.title.x = element_text(size = 6),
      plot.title = element_text(face = "bold", hjust = 0.5, size = 7),
      panel.grid.major.y = element_blank(), panel.grid.minor = element_blank(),
      panel.border = element_rect(color = "black", fill = NA, linewidth = 0.5),
      plot.margin = margin(t = 15, r = 5, b = 5, l = 5),
      plot.background = element_rect(fill = "white", color = NA)
    ) +
    coord_cartesian(clip = "off")

  save_pdf(p, file.path(out_dir, sprintf("%s_GSEA_BP_butterfly.pdf", prefix)), w_in, h_in)
}

# =============================================================================
# Rebecca's selected pathways (2026-04-04 email)
# =============================================================================

BASE_DEG <- file.path(PNAS, "atlas/deg")

# --- Fig 1E: WT Female vs WT Male ---
cat("=== Fig 1E GSEA (custom) ===\n")
make_custom_butterfly(
  file.path(out_dir, "fig1E_GSEA_BP.csv"),
  selected_left = c(  # UP in WT Male (NES < 0)
    "long-chain fatty-acyl-CoA metabolic process",
    "regulation of T-helper 1 type immune response",
    "cholesterol biosynthetic process",
    "inflammatory response to wounding",
    "platelet activation"
  ),
  selected_right = c(  # UP in WT Female (NES > 0)
    "regulation of T cell mediated cytotoxicity",
    "regulation of leukocyte mediated cytotoxicity",
    "miRNA-mediated post-transcriptional gene silencing",
    "RNA splicing",
    "mRNA processing"
  ),
  label_left = "WT Male", label_right = "WT Female",
  prefix = "fig1E", out_dir = out_dir
)

# --- Fig 1G: KO Male vs WT Male ---
cat("=== Fig 1G GSEA (custom) ===\n")
make_custom_butterfly(
  file.path(out_dir, "fig1G_GSEA_BP.csv"),
  selected_left = c(  # UP in WT Male (NES < 0)
    "cellular response to misfolded protein",
    "cholesterol biosynthetic process",
    "negative regulation of carbohydrate metabolic process",
    "sodium ion transmembrane transport"
  ),
  selected_right = "TOP5",  # Keep original top 5 KO Male terms
  label_left = "WT Male", label_right = "KO Male",
  prefix = "fig1G", out_dir = out_dir
)

# --- Multifactor Infection (ctrl for genotype) ---
cat("=== Multifactor Infection GSEA (custom) ===\n")
make_custom_butterfly(
  file.path(out_dir, "multifactor_infection_GSEA_BP.csv"),
  selected_left = c(  # UP in Naive (NES < 0)
    "plasma lipoprotein particle remodeling",
    "complement activation, classical pathway",
    "antigen processing and presentation.*MHC class II",
    "apoptotic cell clearance",
    "positive chemotaxis"
  ),
  selected_right = "TOP5",  # Keep original top Infected terms
  label_left = "Naive", label_right = "Infected",
  prefix = "multifactor_infection", out_dir = out_dir
)

cat("\nDone! Custom GSEA butterfly plots saved.\n")
