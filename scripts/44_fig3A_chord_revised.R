#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 44_fig3A_chord_revised.R
# Revised Fig 3A: cleaner chord diagram
#
# Design:
#   Sectors = cell states (Eos, M1 Mac, M2 Func, M2 Surf)
#   Links = pathway communication (APP, CCL, FN1, CD45, THBS)
#   Link color = pathway identity
#   Link width = signaling strength (probability)
#   Directional arrows show sender -> receiver
#   Pathway names labeled directly on links
#
# Also generates: per-pathway individual panels (for supplement)
#
# Input:  $BASE/cellchat/cellchat_<condition>.rds
# Output: revision/figures/fig3/fig3A_chord_v2.pdf
# SLURM:  batch, 32G, 4 CPUs, 30min
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

# --- Cell state mapping ---
chord_group <- function(x) {
  ifelse(grepl("^Eos_", x), "Eos",
  ifelse(x == "Mac_M1", "M1 Mac",
  ifelse(x == "Mac_M2_Func", "M2 Func",
  ifelse(x == "Mac_M2_Surf", "M2 Surf", NA_character_))))
}

cell_colors <- c("Eos"     = "#777B7E",
                 "M1 Mac"  = "#cb416b",
                 "M2 Func" = "#8FC50E",
                 "M2 Surf" = "#5633AA")

FOCUS_PATHWAYS <- c("APP", "CCL", "FN1", "CD45", "THBS")

pathway_colors <- c(
  "APP"  = "#E41A1C",
  "CCL"  = "#377EB8",
  "FN1"  = "#4DAF4A",
  "CD45" = "#FF7F00",
  "THBS" = "#984EA3"
)

# =============================================================================
# 1. Extract pathway probabilities (same as script 40)
# =============================================================================
cat("Loading CellChat objects...\n")

cell_states <- names(cell_colors)

# per_condition[[key]][[cond]] holds one summed value per (pathway, broad
# source group, broad target group) per condition. A condition in which that
# pathway/group-pair never appears (pathway undetected entirely, or every
# underlying fine-subtype edge was zero) is filled with 0 below rather than
# omitted -- omitting it previously meant mean() divided only by the
# conditions where SOME signal existed, inflating rare/condition-specific
# links relative to ones present in every condition.
per_condition <- list()
loaded_conditions <- character()

for (cond in condition_order) {
  f <- file.path(cc_dir, paste0("cellchat_", cond, ".rds"))
  if (!file.exists(f)) next
  cat(sprintf("  %s ...\n", cond))
  loaded_conditions <- c(loaded_conditions, cond)

  con <- gzfile(f, "rb")
  cc_raw <- .Internal(unserialize(con, NULL))
  close(con)

  netP_prob <- cc_raw@netP$prob
  pathways  <- dimnames(netP_prob)[[3]]
  ct_names  <- dimnames(netP_prob)[[1]]
  groups    <- chord_group(ct_names)
  keep <- !is.na(groups)

  # Sum fine-subtype probabilities into one value per (pathway, broad source,
  # broad target) for this condition; pathways absent from this condition
  # contribute no keys here and are zero-filled below.
  cond_totals <- list()
  for (pw in FOCUS_PATHWAYS) {
    if (!(pw %in% pathways)) next
    prob_mat <- netP_prob[, , pw]
    for (i in which(keep)) {
      for (j in which(keep)) {
        key <- paste(pw, groups[i], groups[j], sep = "|||")
        cond_totals[[key]] <- (if (is.null(cond_totals[[key]])) 0 else cond_totals[[key]]) + prob_mat[i, j]
      }
    }
  }

  for (key in names(cond_totals)) {
    if (is.null(per_condition[[key]])) per_condition[[key]] <- list()
    per_condition[[key]][[cond]] <- cond_totals[[key]]
  }

  rm(cc_raw); gc(verbose = FALSE)
}

n_cond <- length(loaded_conditions)
cat(sprintf("  Loaded %d/%d conditions\n", n_cond, length(condition_order)))

# Average across conditions, zero-filling conditions where a key never appeared
pw_df <- do.call(rbind, lapply(names(per_condition), function(key) {
  parts <- strsplit(key, "\\|\\|\\|")[[1]]
  vals <- vapply(loaded_conditions, function(cond) {
    v <- per_condition[[key]][[cond]]
    if (is.null(v)) 0 else v
  }, numeric(1))
  data.frame(pathway = parts[1], source = parts[2], target = parts[3],
             prob = mean(vals), stringsAsFactors = FALSE)
}))

# Filter to focus pathways, cross-type only
pw_focus <- pw_df[pw_df$pathway %in% FOCUS_PATHWAYS, ]
is_cross <- (pw_focus$source == "Eos" & pw_focus$target != "Eos") |
            (pw_focus$source != "Eos" & pw_focus$target == "Eos")
pw_cross <- pw_focus[is_cross, ]

# Also include Mac-Mac interactions for completeness
pw_mac <- pw_focus[!is_cross & pw_focus$source != "Eos", ]
pw_all <- rbind(pw_cross, pw_mac)

cat(sprintf("\nFocus pathway interactions: %d cross-type, %d mac-mac\n",
            nrow(pw_cross), nrow(pw_mac)))

# Print summary
for (pw in FOCUS_PATHWAYS) {
  sub <- pw_all[pw_all$pathway == pw, ]
  cat(sprintf("\n%s (%d links):\n", pw, nrow(sub)))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("  %s -> %s : %.4f\n", sub$source[i], sub$target[i], sub$prob[i]))
  }
}

# =============================================================================
# 2. Version A: Combined chord — all pathways, pathway labels on links
# =============================================================================
cat("\n\nGenerating Version A: Combined chord with pathway labels...\n")

chord_df <- data.frame(
  from    = pw_cross$source,
  to      = pw_cross$target,
  value   = pw_cross$prob,
  pathway = pw_cross$pathway,
  stringsAsFactors = FALSE
)

link_cols <- adjustcolor(pathway_colors[chord_df$pathway], alpha.f = 0.5)
sector_order <- c("Eos", "M1 Mac", "M2 Func", "M2 Surf")

pdf_out <- file.path(out_dir, "fig3A_chord_v2.pdf")
cairo_pdf(pdf_out, width = 9, height = 9)

par(mar = c(2, 2, 3, 2))
circos.clear()
circos.par(start.degree = 90, gap.after = c(12, 6, 6, 12),
           track.margin = c(0.02, 0.02))

chordDiagram(chord_df[, c("from", "to", "value")],
             order = sector_order,
             grid.col = cell_colors,
             col = link_cols,
             transparency = 0,
             directional = 1,
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             link.arr.length = 0.08,
             link.arr.width = 0.05,
             link.lwd = 2,
             link.sort = TRUE,
             link.decreasing = TRUE,
             annotationTrack = c("grid"),
             preAllocateTracks = list(
               list(track.height = 0.22)
             ),
             scale = TRUE)

# Sector labels — bold, inside track
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector <- CELL_META$sector.index
  xlim   <- CELL_META$xlim
  ylim   <- CELL_META$ylim
  circos.text(mean(xlim), ylim[1] + 0.3, sector,
              facing = "bending.inside", niceFacing = TRUE,
              cex = 1.4, font = 2, col = "black")
}, bg.border = NA)

circos.clear()

# Title
title(main = "Eosinophil-Macrophage Signaling Network",
      cex.main = 1.3, font.main = 2, line = 0.5)

# Pathway legend — left side
legend("bottomleft",
       legend = names(pathway_colors),
       fill = pathway_colors,
       border = "grey30",
       bty = "n", cex = 1.1,
       title = expression(bold("Signaling Pathway")),
       title.adj = 0,
       inset = c(0.02, 0.02))

# Cell state legend — right side
legend("bottomright",
       legend = names(cell_colors),
       fill = cell_colors,
       border = "grey30",
       bty = "n", cex = 1.1,
       title = expression(bold("Cell State")),
       title.adj = 0,
       inset = c(0.02, 0.02))

# Arrow direction note
mtext("Arrow: sender \u2192 receiver", side = 1, line = 0.5, cex = 0.9, col = "grey40")

dev.off()
cat(sprintf("  Saved: %s\n", pdf_out))

# =============================================================================
# 3. Version B: Bipartite chord — pathways as sectors alongside cell types
#    Layout: cell states + pathways as sectors
# =============================================================================
cat("\nGenerating Version B: Bipartite chord (cell types + pathways as sectors)...\n")

# For each pathway-cell interaction, create two links:
#   sender_cell -> pathway  (weight = prob)
#   pathway -> receiver_cell (weight = prob)
# This shows: which cell types use which pathways to communicate

bipartite_links <- data.frame(
  from = character(), to = character(), value = numeric(),
  stringsAsFactors = FALSE
)

for (i in seq_len(nrow(pw_cross))) {
  pw <- pw_cross$pathway[i]
  src <- pw_cross$source[i]
  tgt <- pw_cross$target[i]
  prob <- pw_cross$prob[i]

  # sender -> pathway
  bipartite_links <- rbind(bipartite_links,
    data.frame(from = src, to = pw, value = prob, stringsAsFactors = FALSE))
  # pathway -> receiver
  bipartite_links <- rbind(bipartite_links,
    data.frame(from = pw, to = tgt, value = prob, stringsAsFactors = FALSE))
}

# Aggregate duplicate links
bipartite_agg <- aggregate(value ~ from + to, data = bipartite_links, FUN = sum)

# Determine which sectors actually appear in the data
all_sectors_in_data <- unique(c(bipartite_agg$from, bipartite_agg$to))

# Full desired order: cell types on left, pathways on right
bp_sector_full <- c("Eos", "M1 Mac", "M2 Func", "M2 Surf",
                    "APP", "CCL", "FN1", "CD45", "THBS")
bp_sector_order <- bp_sector_full[bp_sector_full %in% all_sectors_in_data]

# Colors: cell type colors + pathway colors (only present ones)
bp_colors <- c(cell_colors, pathway_colors)
bp_colors <- bp_colors[names(bp_colors) %in% bp_sector_order]

# Gap sizes: bigger gap at cell-type/pathway boundary
n_sectors <- length(bp_sector_order)
bp_gaps <- rep(4, n_sectors)
# Find boundary between cell types and pathways
cell_idx <- which(bp_sector_order %in% names(cell_colors))
pw_idx   <- which(bp_sector_order %in% FOCUS_PATHWAYS)
if (length(cell_idx) > 0 && length(pw_idx) > 0) {
  bp_gaps[max(cell_idx)] <- 15
  bp_gaps[max(pw_idx)]   <- 15
}

cat(sprintf("  Bipartite sectors present: %s\n", paste(bp_sector_order, collapse = ", ")))

# Link colors: by pathway (determine from "from" or "to")
bp_link_cols <- sapply(seq_len(nrow(bipartite_agg)), function(i) {
  fr <- bipartite_agg$from[i]
  to <- bipartite_agg$to[i]
  if (fr %in% FOCUS_PATHWAYS) return(adjustcolor(pathway_colors[fr], alpha.f = 0.45))
  if (to %in% FOCUS_PATHWAYS) return(adjustcolor(pathway_colors[to], alpha.f = 0.45))
  return(adjustcolor("grey", alpha.f = 0.3))
})

pdf_out2 <- file.path(out_dir, "fig3A_chord_v2_bipartite.pdf")
cairo_pdf(pdf_out2, width = 10, height = 10)

par(mar = c(2, 2, 3, 2))
circos.clear()
circos.par(start.degree = 90, gap.after = bp_gaps,
           track.margin = c(0.02, 0.02))

chordDiagram(bipartite_agg,
             order = bp_sector_order,
             grid.col = bp_colors,
             col = bp_link_cols,
             transparency = 0,
             directional = 1,
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             link.arr.length = 0.06,
             link.arr.width = 0.04,
             link.lwd = 1.5,
             link.sort = TRUE,
             link.decreasing = TRUE,
             annotationTrack = "grid",
             preAllocateTracks = list(
               list(track.height = 0.2)
             ))

# Sector labels
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector <- CELL_META$sector.index
  xlim   <- CELL_META$xlim
  ylim   <- CELL_META$ylim

  # Bold for cell types, italic for pathways
  if (sector %in% names(cell_colors)) {
    circos.text(mean(xlim), ylim[1] + 0.3, sector,
                facing = "bending.inside", niceFacing = TRUE,
                cex = 1.2, font = 2, col = "black")
  } else {
    circos.text(mean(xlim), ylim[1] + 0.3, sector,
                facing = "bending.inside", niceFacing = TRUE,
                cex = 1.1, font = 4, col = pathway_colors[sector])
  }
}, bg.border = NA)

circos.clear()

title(main = "Cell State \u2013 Pathway Signaling Architecture",
      cex.main = 1.3, font.main = 2, line = 0.5)

# Group labels
mtext("Cell States", side = 3, line = -1.5, adj = 0.15, cex = 0.9,
      font = 2, col = "grey40")
mtext("Pathways", side = 3, line = -1.5, adj = 0.85, cex = 0.9,
      font = 2, col = "grey40")

# Direction note
mtext("Arrow: sender \u2192 pathway \u2192 receiver", side = 1, line = 0.5,
      cex = 0.9, col = "grey40")

dev.off()
cat(sprintf("  Saved: %s\n", pdf_out2))

# =============================================================================
# 4. Version C: Faceted — one small chord per pathway (2x3 grid)
# =============================================================================
cat("\nGenerating Version C: Per-pathway faceted chord (2x3 grid)...\n")

pdf_out3 <- file.path(out_dir, "fig3A_chord_v2_perpathway.pdf")
cairo_pdf(pdf_out3, width = 12, height = 8)

par(mfrow = c(2, 3), mar = c(1, 1, 3, 1))

for (pw in FOCUS_PATHWAYS) {
  pw_sub <- pw_cross[pw_cross$pathway == pw, ]

  if (nrow(pw_sub) == 0) {
    plot.new()
    title(main = paste0(pw, "\n(no interactions)"), cex.main = 1.1, font.main = 2)
    next
  }

  pw_chord <- data.frame(
    from  = pw_sub$source,
    to    = pw_sub$target,
    value = pw_sub$prob,
    stringsAsFactors = FALSE
  )

  pw_link_col <- rep(adjustcolor(pathway_colors[pw], alpha.f = 0.5), nrow(pw_chord))

  # Only include sectors that appear in this pathway's data
  pw_sectors_present <- unique(c(pw_chord$from, pw_chord$to))
  pw_sector_order <- sector_order[sector_order %in% pw_sectors_present]
  n_pw_sectors <- length(pw_sector_order)
  # Gaps: larger gap between Eos and last Mac
  pw_gaps <- rep(5, n_pw_sectors)
  if (n_pw_sectors > 1) {
    pw_gaps[1] <- 10
    pw_gaps[n_pw_sectors] <- 10
  }

  circos.clear()
  circos.par(start.degree = 90, gap.after = pw_gaps,
             track.margin = c(0.01, 0.01),
             canvas.xlim = c(-1.2, 1.2), canvas.ylim = c(-1.2, 1.2))

  chordDiagram(pw_chord,
               order = pw_sector_order,
               grid.col = cell_colors[pw_sector_order],
               col = pw_link_col,
               transparency = 0,
               directional = 1,
               direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow",
               link.arr.length = 0.08,
               link.arr.width = 0.06,
               link.lwd = 2.5,
               link.sort = TRUE,
               link.decreasing = TRUE,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.18))

  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1] + 0.2,
                CELL_META$sector.index,
                facing = "bending.inside", niceFacing = TRUE,
                cex = 0.9, font = 2)
  }, bg.border = NA)

  circos.clear()

  # Pathway name as panel title, colored
  title(main = pw, cex.main = 1.5, font.main = 2,
        col.main = pathway_colors[pw])

  # Add interaction details below title
  details <- paste(sapply(seq_len(nrow(pw_sub)), function(i) {
    sprintf("%s\u2192%s", pw_sub$source[i], pw_sub$target[i])
  }), collapse = "  ")
  mtext(details, side = 1, line = -0.5, cex = 0.6, col = "grey40")
}

# 6th panel: legend
plot.new()
legend("center",
       legend = c(names(cell_colors), "", names(pathway_colors)),
       fill = c(cell_colors, NA, pathway_colors),
       border = c(rep("grey30", 4), NA, rep("grey30", 5)),
       bty = "n", cex = 1.1,
       title = expression(bold("Legend")))

dev.off()
cat(sprintf("  Saved: %s\n", pdf_out3))

# =============================================================================
# 5. Summary
# =============================================================================
cat("\n========================================\n")
cat("Summary of outputs:\n")
cat("========================================\n")
cat(sprintf("  fig3A_chord_v2.pdf           — Combined chord, pathway-colored links\n"))
cat(sprintf("  fig3A_chord_v2_bipartite.pdf — Bipartite: cell states + pathways as sectors\n"))
cat(sprintf("  fig3A_chord_v2_perpathway.pdf — Per-pathway faceted (2x3 grid)\n"))

cat(sprintf("\nTotal cross-type interactions: %d\n", nrow(pw_cross)))
for (pw in FOCUS_PATHWAYS) {
  n <- sum(pw_cross$pathway == pw)
  total_prob <- sum(pw_cross$prob[pw_cross$pathway == pw])
  cat(sprintf("  %s: %d links, total strength = %.4f\n", pw, n, total_prob))
}

cat("\n44 complete.\n")
