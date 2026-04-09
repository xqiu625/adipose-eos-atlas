#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  library(cowplot)
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

# Cell type display order
ct_order <- c("Eos_C3", "Eos_C5", "Eos_C0", "Eos_C1", "Eos_C2", "Eos_C4",
              "Mac_M1", "Mac_M2_Func", "Mac_M2_Surf")
ct_display <- c("Eos_C3" = "C3", "Eos_C5" = "C5", "Eos_C0" = "C0",
                "Eos_C1" = "C1", "Eos_C2" = "C2", "Eos_C4" = "C4",
                "Mac_M1" = "M1", "Mac_M2_Func" = "M2\nFunct.",
                "Mac_M2_Surf" = "M2\nSurf.")

cat("\nGenerating Fig 3D: CCL Ligand Senders & Receivers...\n")

sender_all   <- list()
receiver_all <- list()

for (cond in names(cellchat_objects)) {
  cc <- cellchat_objects[[cond]]

  # Use standard CellChat API to extract CCL communications
  ccl_comm <- tryCatch(
    subsetCommunication(cc, signaling = "CCL"),
    error = function(e) { cat(sprintf("  %s: subsetCommunication failed: %s\n", cond, e$message)); NULL }
  )

  if (is.null(ccl_comm) || nrow(ccl_comm) == 0) {
    cat(sprintf("  %s: no CCL communications found\n", cond))
    next
  }
  cat(sprintf("  %s: %d CCL L-R pairs\n", cond, nrow(ccl_comm)))

  for (i in seq_len(nrow(ccl_comm))) {
    ligand    <- ccl_comm$ligand[i]
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

if (length(sender_all) == 0) stop("No CCL sender data found in any condition")

sender_df   <- do.call(rbind, sender_all)
receiver_df <- do.call(rbind, receiver_all)

# Average across conditions
avg_sender <- sender_df %>%
  group_by(ligand, celltype) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")
avg_receiver <- receiver_df %>%
  group_by(ligand, celltype) %>%
  summarise(value = mean(value, na.rm = TRUE), .groups = "drop")

# Top ligands by total signal
top_ligands_s <- avg_sender %>% group_by(ligand) %>%
  summarise(total = sum(value)) %>% arrange(desc(total)) %>%
  filter(total > 0) %>% head(7) %>% pull(ligand)
top_ligands_r <- avg_receiver %>% group_by(ligand) %>%
  summarise(total = sum(value)) %>% arrange(desc(total)) %>%
  filter(total > 0) %>% head(7) %>% pull(ligand)

cat(sprintf("  Top sender ligands: %s\n",   paste(top_ligands_s, collapse = ", ")))
cat(sprintf("  Top receiver ligands: %s\n", paste(top_ligands_r, collapse = ", ")))

avg_sender   <- avg_sender   %>% filter(ligand %in% top_ligands_s)
avg_receiver <- avg_receiver %>% filter(ligand %in% top_ligands_r)

avg_sender$ligand   <- factor(avg_sender$ligand,   levels = rev(top_ligands_s))
avg_sender$celltype <- factor(avg_sender$celltype, levels = ct_order)
avg_receiver$ligand   <- factor(avg_receiver$ligand,   levels = rev(top_ligands_r))
avg_receiver$celltype <- factor(avg_receiver$celltype, levels = ct_order)

ct_labels_vec <- ct_display[ct_order]

p_sender <- ggplot(avg_sender, aes(x = celltype, y = ligand, fill = value)) +
  geom_tile(color = "grey70", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "#B2182B", name = "Signal\nStrength") +
  scale_x_discrete(labels = ct_labels_vec) +
  labs(title = "CCL Ligand Senders", x = "", y = "") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 10, face = "italic"),
        plot.title  = element_text(size = 13, face = "bold", hjust = 0.5),
        panel.grid  = element_blank(),
        legend.key.height = unit(1.5, "cm"))

p_receiver <- ggplot(avg_receiver, aes(x = celltype, y = ligand, fill = value)) +
  geom_tile(color = "grey70", linewidth = 0.3) +
  scale_fill_gradient(low = "white", high = "#2166AC", name = "Signal\nStrength") +
  scale_x_discrete(labels = ct_labels_vec) +
  labs(title = "CCL Ligand Receivers", x = "", y = "") +
  theme_minimal(base_size = 11) +
  theme(axis.text.x = element_text(size = 9),
        axis.text.y = element_text(size = 10, face = "italic"),
        plot.title  = element_text(size = 13, face = "bold", hjust = 0.5),
        panel.grid  = element_blank(),
        legend.key.height = unit(1.5, "cm"))

out_pdf <- file.path(fig_dir, "fig3D_ccl_ligand_sender_receiver.pdf")
pdf(out_pdf, width = 16, height = 5)
print(plot_grid(p_sender, p_receiver, nrow = 1))
dev.off()
cat(sprintf("  Saved %s\n", out_pdf))
cat("\n10c complete.\n")
