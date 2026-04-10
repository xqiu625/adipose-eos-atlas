#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 35_fig3A_chord_revised.R
# Fig 3A revised: Chord diagram with enhanced arrow visibility (R1 Comment 18)
#
# Bypasses CellChat/igraph. Extracts @netP$prob via raw deserialization.
# Arrows made larger/bolder per reviewer request.
#
# Input:  $BASE/cellchat/cellchat_<condition>.rds (6 conditions)
# Output: revision/figures/fig3/fig3A_chord_all_pathways_revised.pdf
# SLURM:  batch, 32G, 4 CPUs, 1h
# ==============================================================================

suppressPackageStartupMessages({
  library(circlize)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

BASE    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
cc_dir  <- file.path(BASE, "cellchat")
PNAS    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision"
out_dir <- file.path(PNAS, "revision", "figures", "fig3")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

condition_order <- c("WT_Female_HFD", "WT_Male_HFD", "KO_Male_HFD",
                     "KO_Female_HFD", "KO_Female_HFDInfected",
                     "WT_Female_HFDInfected")

# --- Group mapping ---
chord_group <- function(x) {
  ifelse(grepl("^Eos_", x), "Eosinophils",
  ifelse(x == "Mac_M1", "M1 Mac",
  ifelse(x == "Mac_M2_Func", "M2 Functional",
  ifelse(x == "Mac_M2_Surf", "M2 Surface", NA_character_))))
}

group_colors <- c("Eosinophils"    = "#777B7E",
                  "M1 Mac"         = "#cb416b",
                  "M2 Functional"  = "#8FC50E",
                  "M2 Surface"     = "#5633AA")

unique_groups <- names(group_colors)

# =============================================================================
# 1. Load and aggregate pathway-level probabilities
# =============================================================================
cat("Loading CellChat objects (raw deserialization)...\n")

agg_total <- matrix(0, nrow = 4, ncol = 4,
                    dimnames = list(unique_groups, unique_groups))
n_loaded <- 0

for (cond in condition_order) {
  f <- file.path(cc_dir, paste0("cellchat_", cond, ".rds"))
  if (!file.exists(f)) next
  cat(sprintf("  %s ...\n", cond))

  con <- gzfile(f, "rb")
  cc_raw <- .Internal(unserialize(con, NULL))
  close(con)

  # Extract pathway-level prob: @netP$prob is a 3D array [source, target, pathway]
  netP_prob <- cc_raw@netP$prob
  pathways  <- dimnames(netP_prob)[[3]]
  ct_names  <- dimnames(netP_prob)[[1]]
  groups    <- chord_group(ct_names)

  cat(sprintf("    %d pathways, %d cell types\n", length(pathways), length(ct_names)))

  for (pw in pathways) {
    prob_mat <- netP_prob[, , pw]
    keep <- !is.na(groups)
    for (i in which(keep)) {
      for (j in which(keep)) {
        agg_total[groups[i], groups[j]] <-
          agg_total[groups[i], groups[j]] + prob_mat[i, j]
      }
    }
  }

  n_loaded <- n_loaded + 1
  rm(cc_raw); gc(verbose = FALSE)
}

agg_total <- agg_total / n_loaded
cat(sprintf("\nAggregated %d conditions. Matrix:\n", n_loaded))
print(round(agg_total, 4))

# =============================================================================
# 2. Chord diagram with enhanced arrows (Reviewer Comment 18)
# =============================================================================
cat("\nGenerating chord diagram with enhanced arrows...\n")

pdf_out <- file.path(out_dir, "fig3A_chord_all_pathways_revised.pdf")
cairo_pdf(pdf_out, width = 9, height = 9)

circos.clear()
circos.par(start.degree = 90, gap.after = c(rep(2, 3), 10))

chordDiagram(agg_total,
             grid.col = group_colors,
             transparency = 0.15,
             directional = 1,
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             link.arr.length = 0.06,
             link.arr.width = 0.04,
             link.lwd = 1.5,
             link.sort = TRUE,
             link.decreasing = TRUE,
             annotationTrack = "grid",
             preAllocateTracks = list(track.height = 0.15))

# Add sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5), cex = 1.2,
              font = 2)
}, bg.border = NA)

circos.clear()
dev.off()
cat(sprintf("  Saved: %s\n", pdf_out))

cat("\n35 complete.\n")
