#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))
if (nzchar(Sys.getenv("RERUN_R_LIB"))) {
  .rerun_r_lib_resolved <- normalizePath(.libPaths()[1], mustWork = FALSE)
  .rerun_r_lib_expected <- normalizePath(Sys.getenv("RERUN_R_LIB"), mustWork = FALSE)
  if (!identical(.rerun_r_lib_resolved, .rerun_r_lib_expected)) {
    stop(sprintf(
      "RERUN_R_LIB audit guard: .libPaths()[1] resolved to '%s' but RERUN_R_LIB='%s' -- the shared/legacy library would take precedence over the pinned isolated library. Refusing to proceed.",
      .rerun_r_lib_resolved, .rerun_r_lib_expected))
  }
  cat("=== Pinned R environment audit guard ===\n")
  cat("R.version.string:", R.version.string, "\n")
  cat("RERUN_R_LIB:", Sys.getenv("RERUN_R_LIB"), "\n")
  cat("CMDSTAN:", Sys.getenv("CMDSTAN"), "\n")
  cat(".libPaths():\n")
  print(.libPaths())
  cat("=========================================\n")
}

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
  library(grid)
})

BASE    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
cc_dir  <- file.path(BASE, "cellchat")
fig_dir <- file.path(BASE, "cellchat", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

condition_order <- c("WT_Female_HFD", "WT_Male_HFD", "KO_Male_HFD",
                     "KO_Female_HFD", "KO_Female_HFDInfected", "WT_Female_HFDInfected")

cat("Loading CellChat objects...\n")
cellchat_objects <- list()
for (cond in condition_order) {
  f <- file.path(cc_dir, paste0("cellchat_", cond, ".rds"))
  if (file.exists(f)) {
    cellchat_objects[[cond]] <- readRDS(f)
    cat(sprintf("  Loaded %s\n", cond))
  }
}
cat(sprintf("Loaded %d conditions\n", length(cellchat_objects)))

# Column order: C3 first, C6 excluded, then C5 C0 C1 C2 C4, then Mac
ct_order <- c("Eos_C3", "Eos_C5", "Eos_C0", "Eos_C1", "Eos_C2", "Eos_C4",
              "Mac_M1", "Mac_M2_Func", "Mac_M2_Surf")
ct_display <- c("Eos_C3" = "C3", "Eos_C5" = "C5",
                "Eos_C0" = "C0", "Eos_C1" = "C1", "Eos_C2" = "C2", "Eos_C4" = "C4",
                "Mac_M1" = "M1", "Mac_M2_Func" = "M2\nFunct.", "Mac_M2_Surf" = "M2\nSurf.")

cat("\nGenerating Fig 3D: CCL Ligand Senders & Receivers...\n")

sender_all   <- list()
receiver_all <- list()

for (cond in names(cellchat_objects)) {
  cc <- cellchat_objects[[cond]]

  ccl_comm <- tryCatch(
    subsetCommunication(cc, signaling = "CCL"),
    error = function(e) { cat(sprintf("  %s: %s\n", cond, e$message)); NULL }
  )
  if (is.null(ccl_comm) || nrow(ccl_comm) == 0) {
    cat(sprintf("  %s: no CCL communications\n", cond)); next
  }
  cat(sprintf("  %s: %d CCL pairs\n", cond, nrow(ccl_comm)))

  for (i in seq_len(nrow(ccl_comm))) {
    # Capitalize ligand: e.g. CCL6 -> Ccl6
    raw_lig   <- ccl_comm$ligand[i]
    ligand    <- paste0(toupper(substr(raw_lig, 1, 1)),
                        tolower(substr(raw_lig, 2, nchar(raw_lig))))
    source_ct <- ccl_comm$source[i]
    target_ct <- ccl_comm$target[i]
    prob_val  <- ccl_comm$prob[i]

    if (source_ct %in% ct_order) {
      sender_all <- c(sender_all, list(data.frame(
        ligand = ligand, celltype = source_ct, value = prob_val,
        condition = cond, stringsAsFactors = FALSE)))
    }
    if (target_ct %in% ct_order) {
      receiver_all <- c(receiver_all, list(data.frame(
        ligand = ligand, celltype = target_ct, value = prob_val,
        condition = cond, stringsAsFactors = FALSE)))
    }
  }
}

if (length(sender_all) == 0) stop("No CCL data found")

sender_df   <- do.call(rbind, sender_all)
receiver_df <- do.call(rbind, receiver_all)

# Average across conditions
avg_sender <- sender_df %>%
  group_by(ligand, celltype) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
avg_receiver <- receiver_df %>%
  group_by(ligand, celltype) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

# Top ligands by total signal (up to 6-8)
top_ligands_s <- avg_sender %>% group_by(ligand) %>%
  summarise(total = sum(value)) %>% arrange(desc(total)) %>%
  filter(total > 0) %>% head(6) %>% pull(ligand)
top_ligands_r <- avg_receiver %>% group_by(ligand) %>%
  summarise(total = sum(value)) %>% arrange(desc(total)) %>%
  filter(total > 0) %>% head(6) %>% pull(ligand)

cat(sprintf("  Top sender ligands:   %s\n", paste(top_ligands_s, collapse = ", ")))
cat(sprintf("  Top receiver ligands: %s\n", paste(top_ligands_r, collapse = ", ")))

# Filter and factor
avg_sender   <- avg_sender   %>% filter(ligand %in% top_ligands_s)
avg_receiver <- avg_receiver %>% filter(ligand %in% top_ligands_r)

# Ensure all cell type x ligand combinations exist (fill missing with 0)
full_sender <- expand.grid(ligand   = top_ligands_s,
                           celltype = ct_order, stringsAsFactors = FALSE)
full_receiver <- expand.grid(ligand   = top_ligands_r,
                             celltype = ct_order, stringsAsFactors = FALSE)
avg_sender   <- full_sender   %>% left_join(avg_sender,   by = c("ligand","celltype")) %>%
  mutate(value = replace(value, is.na(value), 0))
avg_receiver <- full_receiver %>% left_join(avg_receiver, by = c("ligand","celltype")) %>%
  mutate(value = replace(value, is.na(value), 0))

# Factor levels
avg_sender$ligand   <- factor(avg_sender$ligand,   levels = rev(top_ligands_s))
avg_sender$celltype <- factor(avg_sender$celltype, levels = ct_order)
avg_receiver$ligand   <- factor(avg_receiver$ligand,   levels = rev(top_ligands_r))
avg_receiver$celltype <- factor(avg_receiver$celltype, levels = ct_order)

ct_labels_vec <- ct_display[ct_order]

# Color scale limits
max_s <- max(avg_sender$value,   na.rm = TRUE)
max_r <- max(avg_receiver$value, na.rm = TRUE)

# ---- Helper: build one heatmap panel ----------------------------------------
make_heatmap <- function(df, title, fill_high, max_val, show_y = TRUE) {
  p <- ggplot(df, aes(x = celltype, y = ligand, fill = value)) +
    geom_tile(color = "grey70", linewidth = 0.4) +
    scale_fill_gradient(low = "white", high = fill_high,
                        limits = c(0, max_val),
                        breaks = pretty(c(0, max_val), n = 5),
                        name = "") +
    scale_x_discrete(labels = ct_labels_vec, expand = c(0, 0)) +
    scale_y_discrete(position = "right", expand = c(0, 0)) +
    labs(title = title, x = "", y = "") +
    theme_minimal(base_size = 11) +
    theme(
      axis.text.x     = element_text(size = 9, color = "black"),
      axis.text.y     = element_text(size = 10, color = "black",
                                     face = "italic", hjust = 0),
      plot.title      = element_text(size = 12, face = "bold", hjust = 0.5),
      panel.grid      = element_blank(),
      panel.border    = element_blank(),
      axis.ticks      = element_blank(),
      legend.position      = "right",
      legend.key.height    = unit(1.8, "cm"),
      legend.key.width     = unit(0.35, "cm"),
      legend.text          = element_text(size = 9),
      plot.margin     = margin(t = 5, r = 5, b = 30, l = 5)
    )

  if (!show_y) {
    p <- p + theme(axis.text.y = element_blank())
  }
  p
}

p_sender   <- make_heatmap(avg_sender,   "CCL Ligand Senders",   "#B2182B", max_s)
p_receiver <- make_heatmap(avg_receiver, "CCL Ligand Receivers", "#2166AC", max_r)

# ---- Add group underline annotations ----------------------------------------
n_eos <- 6  # C3 C5 C0 C1 C2 C4
n_mac <- 3  # M1 M2F M2S

add_group_labels <- function(p, n_eos, n_mac) {
  p +
    # Underline under Eos columns
    annotate("segment",
             x = 0.5, xend = n_eos + 0.5,
             y = 0.4, yend = 0.4,
             color = "black", linewidth = 0.6,
             xmin = -Inf, xmax = Inf, ymin = -Inf, ymax = Inf) +
    # "Eosinophils" label
    annotate("text",
             x = (n_eos + 1) / 2, y = 0.05,
             label = "Eosinophils", size = 3.2, fontface = "plain",
             hjust = 0.5, vjust = 1) +
    # "Funct." and "Surf." labels under M2 subtypes
    annotate("text", x = n_eos + 2, y = 0.05,
             label = "Funct.", size = 3, hjust = 0.5, vjust = 1) +
    annotate("text", x = n_eos + 3, y = 0.05,
             label = "Surf.",  size = 3, hjust = 0.5, vjust = 1) +
    coord_cartesian(clip = "off")
}

p_sender   <- add_group_labels(p_sender,   n_eos, n_mac)
p_receiver <- add_group_labels(p_receiver, n_eos, n_mac)

# ---- Combine and save -------------------------------------------------------
final <- plot_grid(p_sender, p_receiver, nrow = 1, rel_widths = c(1, 1),
                   align = "h", axis = "tb")

out_pdf <- file.path(fig_dir, "fig3D_ccl_ligand_sender_receiver.pdf")
pdf(out_pdf, width = 14, height = 4.5)
print(final)
dev.off()
cat(sprintf("  Saved %s\n", out_pdf))
cat("\n10c complete.\n")
