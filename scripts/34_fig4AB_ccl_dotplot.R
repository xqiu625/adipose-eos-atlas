#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 34_fig4AB_ccl_dotplot.R
# Fig 4A: Top 10 CCL interactions Mac -> Eos (dot plot per condition)
# Fig 4B: Top 10 CCL interactions Eos -> Mac (dot plot per condition)
#
# Dot size + color = interaction strength (viridis, yellow=strong)
# X-axis = conditions, Y-axis = sender->receiver + L-R pair
#
# NOTE: Bypasses CellChat/igraph (raw deserialization)
#
# Input:  $BASE/cellchat/cellchat_<condition>.rds (6 conditions)
# Output: revision/figures/fig4/fig4A_ccl_mac_to_eos_dotplot.pdf
#         revision/figures/fig4/fig4B_ccl_eos_to_mac_dotplot.pdf
#         revision/figures/fig4/fig4AB_ccl_dotplots_combined.pdf
# SLURM:  batch, 32G, 4 CPUs, 1h
# ==============================================================================

suppressPackageStartupMessages({
  library(ggplot2)
  library(dplyr)
  library(cowplot)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

BASE    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
cc_dir  <- file.path(BASE, "cellchat")
PNAS    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision"
out_dir <- file.path(PNAS, "revision", "figures", "fig4")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

condition_order <- c("WT_Female_HFD", "WT_Male_HFD", "KO_Male_HFD",
                     "KO_Female_HFD", "KO_Female_HFDInfected",
                     "WT_Female_HFDInfected")

# =============================================================================
# 1. Load CellChat objects (raw deserialization, bypass igraph)
# =============================================================================
cat("Loading CellChat objects...\n")

extract_ccl_from_rds <- function(rds_path) {
  con <- gzfile(rds_path, "rb")
  cc_raw <- .Internal(unserialize(con, NULL))
  close(con)

  prob_arr <- cc_raw@net$prob
  lrsig    <- cc_raw@LR$LRsig

  ccl_idx <- which(lrsig$pathway_name == "CCL")
  if (length(ccl_idx) == 0) return(NULL)

  ccl_lrsig <- lrsig[ccl_idx, ]
  interaction_names <- dimnames(prob_arr)[[3]]
  source_names <- dimnames(prob_arr)[[1]]
  target_names <- dimnames(prob_arr)[[2]]

  rows <- list()
  for (i in seq_len(nrow(ccl_lrsig))) {
    int_name <- ccl_lrsig$interaction_name[i]
    ligand   <- ccl_lrsig$ligand[i]
    receptor <- ccl_lrsig$receptor[i]

    int_idx <- which(interaction_names == int_name)
    if (length(int_idx) == 0 && "interaction_name_2" %in% colnames(ccl_lrsig)) {
      int_idx <- which(interaction_names == ccl_lrsig$interaction_name_2[i])
    }
    if (length(int_idx) == 0) next

    prob_slice <- prob_arr[, , int_idx]
    nz <- which(prob_slice > 0, arr.ind = TRUE)
    if (nrow(nz) == 0) next

    for (j in seq_len(nrow(nz))) {
      rows <- c(rows, list(data.frame(
        source = source_names[nz[j, 1]],
        target = target_names[nz[j, 2]],
        ligand = ligand, receptor = receptor,
        prob = prob_slice[nz[j, 1], nz[j, 2]],
        stringsAsFactors = FALSE
      )))
    }
  }
  rm(cc_raw); gc(verbose = FALSE)
  if (length(rows) == 0) return(NULL)
  do.call(rbind, rows)
}

all_comm <- list()
for (cond in condition_order) {
  f <- file.path(cc_dir, paste0("cellchat_", cond, ".rds"))
  if (!file.exists(f)) next
  cat(sprintf("  %s ...\n", cond))
  ccl_df <- tryCatch(extract_ccl_from_rds(f),
                      error = function(e) { cat(sprintf("  ERROR: %s\n", e$message)); NULL })
  if (!is.null(ccl_df) && nrow(ccl_df) > 0) {
    ccl_df$condition <- cond
    all_comm <- c(all_comm, list(ccl_df))
    cat(sprintf("    %d CCL pairs\n", nrow(ccl_df)))
  }
}
comm_df <- do.call(rbind, all_comm)
cat(sprintf("\nTotal CCL interactions: %d\n", nrow(comm_df)))

# =============================================================================
# 2. Prepare labels
# =============================================================================
eos_types <- paste0("Eos_C", 0:5)
mac_types <- c("Mac_M1", "Mac_M2_Func", "Mac_M2_Surf")

ct_short <- c(
  "Eos_C0" = "EosC0", "Eos_C1" = "EosC1", "Eos_C2" = "EosC2",
  "Eos_C3" = "EosC3", "Eos_C4" = "EosC4", "Eos_C5" = "EosC5",
  "Mac_M1" = "M1", "Mac_M2_Func" = "M2 Funct.", "Mac_M2_Surf" = "M2 Surf"
)

# L-R pair in uppercase (CCL6-CCR3)
comm_df$lr_pair <- paste0(toupper(comm_df$ligand), "-", toupper(comm_df$receptor))

# Unique interaction ID
comm_df$int_id <- paste0(comm_df$source, "|", comm_df$target, "|", comm_df$lr_pair)

# Two-line display label: "M1 -> EosC3\nCCL6-CCR3"
comm_df$label <- paste0(ct_short[comm_df$source], " -> ",
                         ct_short[comm_df$target], "\n",
                         comm_df$lr_pair)

# Condition display labels with sex symbols
cond_display <- c(
  "WT_Female_HFD"          = "WT\u2640",
  "WT_Male_HFD"            = "WT\u2642",
  "KO_Male_HFD"            = "KO\u2642",
  "KO_Female_HFD"          = "KO\u2640",
  "KO_Female_HFDInfected"  = "KO\u2640\nInf",
  "WT_Female_HFDInfected"  = "WT\u2640\nInf"
)

# =============================================================================
# 3. Dot plot helper
# =============================================================================
make_dotplot <- function(df, title, top_n = 10) {
  # Get unique label per int_id
  label_map <- df %>%
    select(int_id, label) %>%
    distinct() %>%
    group_by(int_id) %>%
    summarise(label = first(label), .groups = "drop")

  # Per-condition probability (take max if duplicate L-R in same condition)
  df_slim <- df %>%
    group_by(int_id, condition) %>%
    summarise(prob = max(prob, na.rm = TRUE), .groups = "drop")

  # Complete grid: all interactions x all conditions (fill missing = 0)
  full_grid <- expand.grid(int_id = unique(df_slim$int_id),
                           condition = condition_order,
                           stringsAsFactors = FALSE)
  full_df <- full_grid %>%
    left_join(df_slim, by = c("int_id", "condition")) %>%
    mutate(prob = ifelse(is.na(prob), 0, prob))

  # Rank by mean probability across ALL 6 conditions
  ranking <- full_df %>%
    group_by(int_id) %>%
    summarise(mean_prob = mean(prob), .groups = "drop") %>%
    left_join(label_map, by = "int_id") %>%
    arrange(desc(mean_prob)) %>%
    head(top_n)

  cat(sprintf("  Top %d interactions for %s:\n", top_n, title))
  for (i in seq_len(nrow(ranking))) {
    cat(sprintf("    %d. %s (mean=%.4f)\n", i,
                gsub("\n", " / ", ranking$label[i]), ranking$mean_prob[i]))
  }

  # Filter to top N, exclude zero-prob dots
  dot_df <- full_df %>%
    filter(int_id %in% ranking$int_id) %>%
    left_join(label_map, by = "int_id") %>%
    filter(prob > 0)

  # Factor levels (strongest at top)
  dot_df$label <- factor(dot_df$label, levels = rev(ranking$label))
  dot_df$condition <- factor(dot_df$condition, levels = condition_order)

  # Global max for consistent scale
  max_prob <- max(dot_df$prob, na.rm = TRUE)
  size_breaks <- pretty(c(0, max_prob), n = 5)
  size_breaks <- size_breaks[size_breaks > 0]

  p <- ggplot(dot_df, aes(x = condition, y = label, size = prob, color = prob)) +
    geom_point() +
    scale_size_continuous(range = c(1, 9), name = "Strength",
                          breaks = size_breaks,
                          limits = c(0, max_prob * 1.05)) +
    scale_color_viridis_c(option = "viridis", name = "Strength",
                          breaks = size_breaks,
                          limits = c(0, max_prob * 1.05)) +
    scale_x_discrete(labels = cond_display[condition_order]) +
    labs(title = title, x = NULL, y = NULL) +
    theme_pub(base_size = 11) +
    theme(
      axis.text.y     = element_text(size = 9, lineheight = 0.85),
      axis.text.x     = element_text(size = 9, lineheight = 0.85),
      plot.title       = element_text(hjust = 0.5, face = "bold", size = 12),
      panel.grid.major = element_line(color = "grey92", linewidth = 0.3),
      legend.position  = "right"
    ) +
    guides(color = guide_colorbar(order = 1, barheight = unit(3, "cm"),
                                   barwidth = unit(0.4, "cm")),
           size  = guide_legend(order = 2, override.aes = list(color = "black")))

  p
}

# =============================================================================
# 4. Generate plots
# =============================================================================

# --- Fig 4A: Mac -> Eos ---
cat("\n--- Fig 4A: Macrophage -> Eosinophil ---\n")
mac_eos <- comm_df %>% filter(source %in% mac_types & target %in% eos_types)
cat(sprintf("  Mac->Eos interactions: %d\n", nrow(mac_eos)))
p_a <- make_dotplot(mac_eos, "Macrophage -> Eosinophil")

pdf_a <- file.path(out_dir, "fig4A_ccl_mac_to_eos_dotplot.pdf")
cairo_pdf(pdf_a, width = 7, height = 6)
print(p_a)
dev.off()
cat(sprintf("  Saved: %s\n", pdf_a))

# --- Fig 4B: Eos -> Mac ---
cat("\n--- Fig 4B: Eosinophil -> Macrophage ---\n")
eos_mac <- comm_df %>% filter(source %in% eos_types & target %in% mac_types)
cat(sprintf("  Eos->Mac interactions: %d\n", nrow(eos_mac)))
p_b <- make_dotplot(eos_mac, "Eosinophil -> Macrophage")

pdf_b <- file.path(out_dir, "fig4B_ccl_eos_to_mac_dotplot.pdf")
cairo_pdf(pdf_b, width = 7, height = 6)
print(p_b)
dev.off()
cat(sprintf("  Saved: %s\n", pdf_b))

# --- Combined ---
cat("\nCombining panels...\n")
combined <- plot_grid(p_a, p_b, labels = c("A", "B"), label_size = 14,
                      nrow = 1, rel_widths = c(1, 1), align = "h", axis = "tb")

pdf_ab <- file.path(out_dir, "fig4AB_ccl_dotplots_combined.pdf")
cairo_pdf(pdf_ab, width = 14, height = 6)
print(combined)
dev.off()
cat(sprintf("  Saved: %s\n", pdf_ab))

cat("\n34 complete.\n")
