#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 10b_cellchat_manuscript_figures_fixed_v2.R
# Fixed: added circlize, robust loading, cleaner chord diagram
# ==============================================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(ggplot2)
  library(dplyr)
  library(patchwork)
  library(pheatmap)
  library(circlize)          # ← This was missing!
  library(RColorBrewer)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
cc_dir  <- file.path(BASE, "cellchat")
fig_dir <- file.path(BASE, "cellchat", "figures")
dir.create(fig_dir, recursive = TRUE, showWarnings = FALSE)

selected_pathways <- c("APP", "CCL", "FN1", "CD45", "THBS")
condition_order <- c("WT_Female_HFD", "WT_Male_HFD", "KO_Female_HFD", "KO_Male_HFD",
                     "WT_Female_HFDInfected", "KO_Female_HFDInfected")

# ==============================================================================
# 1. Load CellChat objects
# ==============================================================================
cat("Loading CellChat objects...\n")

cellchat_objects <- list()
cc_files <- list.files(cc_dir, pattern = "^cellchat_.*\\.rds$", full.names = TRUE)

for (f in cc_files) {
  cond <- gsub("^cellchat_|\\.rds$", "", basename(f))
  if (cond %in% condition_order) {
    cat(sprintf("  Loading %s\n", cond))
    cc <- readRDS(f)
    try(cc <- netAnalysis_computeCentrality(cc), silent = TRUE)
    cellchat_objects[[cond]] <- cc
  }
}

cat(sprintf("Loaded %d conditions\n", length(cellchat_objects)))

# ==============================================================================
# Helper: Chord Group
# ==============================================================================
chord_group <- function(x) {
  dplyr::case_when(
    grepl("^Eos_", x) ~ "Eosinophils",
    x == "Mac_M1"     ~ "M1_Mac",
    x == "Mac_M2_Func"~ "M2_Functional",
    x == "Mac_M2_Surf"~ "M2_Surface",
    TRUE              ~ NA_character_
  )
}

# ==============================================================================
# Fig 3A: CCL Chord Diagram (baseline = WT_Female_HFD)
# ==============================================================================
cat("\nGenerating Fig 3A: CCL Chord Diagram...\n")

baseline_cond <- "WT_Female_HFD"
if (!is.null(cellchat_objects[[baseline_cond]])) {
  cc <- cellchat_objects[[baseline_cond]]
  
  if ("CCL" %in% cc@netP$pathways) {
    prob_mat <- cc@netP$prob[,, "CCL"]
    
    groups <- chord_group(rownames(prob_mat))
    keep <- !is.na(groups)
    
    if (sum(keep) > 0) {
      prob_sub <- prob_mat[keep, keep]
      group_labels <- groups[keep]
      
      # Aggregate to 4 groups
      unique_groups <- c("Eosinophils", "M1_Mac", "M2_Functional", "M2_Surface")
      agg_mat <- matrix(0, nrow=4, ncol=4, dimnames = list(unique_groups, unique_groups))
      
      for (i in seq_along(group_labels)) {
        for (j in seq_along(group_labels)) {
          agg_mat[group_labels[i], group_labels[j]] <- 
            agg_mat[group_labels[i], group_labels[j]] + prob_sub[i, j]
        }
      }
      
      # Draw chord diagram
      pdf(file.path(fig_dir, "fig3A_ccl_chord_diagram.pdf"), width = 9, height = 9)
      circos.clear()
      circos.par(start.degree = 90, gap.after = c(rep(2, 3), 10))
      
      chordDiagram(agg_mat, 
                   grid.col = group_colors,   # from 00_theme.R
                   transparency = 0.25,
                   directional = 1,
                   direction.type = c("diffHeight", "arrows"),
                   link.arr.type = "big.arrow",
                   annotationTrack = "grid",
                   preAllocateTracks = list(track.height = 0.15))
      
      circos.track(track.index = 1, panel.fun = function(x, y) {
        circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                    facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.1)
      }, bg.border = NA)
      
      title(main = paste("CCL Signaling —", baseline_cond), cex.main = 1.4)
      circos.clear()
      dev.off()
      
      cat("  Saved fig3A_ccl_chord_diagram.pdf\n")
    }
  } else {
    cat("  CCL pathway not found in baseline CellChat object\n")
  }
} else {
  cat("  Baseline condition", baseline_cond, "not found\n")
}

cat("\n10b CellChat figures (v2) finished.\n")
cat("Check", fig_dir, "for output.\n")