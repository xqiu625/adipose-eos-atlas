#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 44b_fig3A_chord_cellchat_style.R
# Fig 3A: MDPI/CellChat-style chord diagram
#
# Design (matching Cells 12 00045 g007 reference):
#   - Outer arcs = cell types (Eos, M1 Mac, M2 Func, M2 Surf)
#   - Each cell type arc subdivided by pathway (APP, CCL, FN1, CD45, THBS)
#   - Pathway names labeled on the subdivisions
#   - Links connect cell type-pathway subsections
#   - Link color = sender cell type
#   - Link width = signaling strength
#
# SLURM: batch, 32G, 4 CPUs, 30min
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

# =============================================================================
# 1. Extract pathway probabilities
# =============================================================================
cat("Loading CellChat objects...\n")

cell_states <- names(cell_colors)
all_pathway_data <- list()

for (cond in condition_order) {
  f <- file.path(cc_dir, paste0("cellchat_", cond, ".rds"))
  if (!file.exists(f)) next
  cat(sprintf("  %s ...\n", cond))

  con <- gzfile(f, "rb")
  cc_raw <- .Internal(unserialize(con, NULL))
  close(con)

  netP_prob <- cc_raw@netP$prob
  pathways  <- dimnames(netP_prob)[[3]]
  ct_names  <- dimnames(netP_prob)[[1]]
  groups    <- chord_group(ct_names)

  for (pw in intersect(FOCUS_PATHWAYS, pathways)) {
    prob_mat <- netP_prob[, , pw]
    keep <- !is.na(groups)
    for (i in which(keep)) {
      for (j in which(keep)) {
        val <- prob_mat[i, j]
        if (val > 0) {
          key <- paste(pw, groups[i], groups[j], sep = "|||")
          if (is.null(all_pathway_data[[key]])) {
            all_pathway_data[[key]] <- list(
              pathway = pw, source = groups[i], target = groups[j], vals = c()
            )
          }
          all_pathway_data[[key]]$vals <- c(all_pathway_data[[key]]$vals, val)
        }
      }
    }
  }
  rm(cc_raw); gc(verbose = FALSE)
}

# Average across conditions
pw_df <- do.call(rbind, lapply(all_pathway_data, function(x) {
  data.frame(pathway = x$pathway, source = x$source, target = x$target,
             prob = mean(x$vals), stringsAsFactors = FALSE)
}))
cat(sprintf("\nTotal interactions: %d\n", nrow(pw_df)))

# Print summary
for (pw in FOCUS_PATHWAYS) {
  sub <- pw_df[pw_df$pathway == pw, ]
  cat(sprintf("\n%s (%d links):\n", pw, nrow(sub)))
  for (i in seq_len(nrow(sub))) {
    cat(sprintf("  %s -> %s : %.4f\n", sub$source[i], sub$target[i], sub$prob[i]))
  }
}

# =============================================================================
# 2. Build compound sectors: cell_type.pathway
# =============================================================================
cat("\nBuilding compound chord data...\n")

# Create sector names: "CellType.Pathway"
# Each link goes from source_celltype.pathway to target_celltype.pathway
chord_links <- data.frame(
  from = paste0(pw_df$source, ".", pw_df$pathway),
  to   = paste0(pw_df$target, ".", pw_df$pathway),
  value = pw_df$prob,
  source_cell = pw_df$source,
  pathway = pw_df$pathway,
  stringsAsFactors = FALSE
)

# Determine which compound sectors actually exist
all_compound_sectors <- unique(c(chord_links$from, chord_links$to))

# Group by cell type for ordering
sector_groups <- sub("\\..*", "", all_compound_sectors)
names(sector_groups) <- all_compound_sectors

# Order: within each cell type, order pathways consistently
cell_order <- c("Eos", "M1 Mac", "M2 Func", "M2 Surf")
ordered_sectors <- c()
for (ct in cell_order) {
  ct_sectors <- all_compound_sectors[sector_groups == ct]
  # Order by pathway
  pw_order <- sapply(ct_sectors, function(s) match(sub(".*\\.", "", s), FOCUS_PATHWAYS))
  ordered_sectors <- c(ordered_sectors, ct_sectors[order(pw_order)])
}

cat(sprintf("  Compound sectors: %d\n", length(ordered_sectors)))
cat(sprintf("  Sector order: %s\n", paste(ordered_sectors, collapse = ", ")))

# Grid colors: each subsector inherits parent cell type color
grid_colors <- setNames(
  cell_colors[sub("\\..*", "", ordered_sectors)],
  ordered_sectors
)

# Link colors: by sender cell type (semi-transparent)
link_colors <- sapply(seq_len(nrow(chord_links)), function(i) {
  ct <- chord_links$source_cell[i]
  adjustcolor(cell_colors[ct], alpha.f = 0.4)
})

# Gap sizes: large gaps between cell types, small within
gap_sizes <- rep(1, length(ordered_sectors))
# Find boundaries between cell types
for (i in seq_len(length(ordered_sectors) - 1)) {
  ct_curr <- sub("\\..*", "", ordered_sectors[i])
  ct_next <- sub("\\..*", "", ordered_sectors[i + 1])
  if (ct_curr != ct_next) {
    gap_sizes[i] <- 8
  }
}
gap_sizes[length(ordered_sectors)] <- 8  # last to first

# =============================================================================
# 3. Plot
# =============================================================================
cat("\nGenerating chord diagram...\n")

pdf_out <- file.path(out_dir, "fig3A_chord_cellchat_style.pdf")
cairo_pdf(pdf_out, width = 11, height = 11)

par(mar = c(2, 2, 3, 2))
circos.clear()
circos.par(
  start.degree = 90,
  gap.after = gap_sizes,
  track.margin = c(0.01, 0.01),
  cell.padding = c(0, 0, 0, 0)
)

chordDiagram(
  chord_links[, c("from", "to", "value")],
  order = ordered_sectors,
  grid.col = grid_colors,
  col = link_colors,
  transparency = 0,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow",
  link.arr.length = 0.05,
  link.arr.width = 0.03,
  link.lwd = 1.2,
  link.sort = TRUE,
  link.decreasing = TRUE,
  annotationTrack = "grid",
  preAllocateTracks = list(
    list(track.height = 0.06),  # track 1: pathway labels
    list(track.height = 0.08)   # track 2: cell type arcs
  )
)

# Track 1 (inner): pathway labels on each subsector
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector <- CELL_META$sector.index
  xlim   <- CELL_META$xlim
  ylim   <- CELL_META$ylim

  # Extract pathway name
  pw_name <- sub(".*\\.", "", sector)

  circos.text(mean(xlim), mean(ylim), pw_name,
              facing = "clockwise", niceFacing = TRUE,
              cex = 0.7, font = 3, adj = c(0, 0.5))
}, bg.border = NA)

# Track 2 (outer): cell type colored arcs + bold labels
# We'll add cell type labels manually using highlight
circos.clear()

# Re-draw without clearing to add highlights...
# Actually, let's use a different approach: draw everything in one pass

circos.clear()
circos.par(
  start.degree = 90,
  gap.after = gap_sizes,
  track.margin = c(0.01, 0.01),
  cell.padding = c(0, 0, 0, 0)
)

chordDiagram(
  chord_links[, c("from", "to", "value")],
  order = ordered_sectors,
  grid.col = grid_colors,
  col = link_colors,
  transparency = 0,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow",
  link.arr.length = 0.05,
  link.arr.width = 0.03,
  link.lwd = 1.2,
  link.sort = TRUE,
  link.decreasing = TRUE,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.12)
)

# Pathway labels on each subsector
circos.track(track.index = 1, panel.fun = function(x, y) {
  sector <- CELL_META$sector.index
  xlim   <- CELL_META$xlim
  ylim   <- CELL_META$ylim

  pw_name <- sub(".*\\.", "", sector)

  circos.text(mean(xlim), ylim[1] + 0.15, pw_name,
              facing = "clockwise", niceFacing = TRUE,
              cex = 0.65, font = 3, adj = c(0, 0.5))
}, bg.border = NA)

# Add cell type group labels using highlight.sector
for (ct in cell_order) {
  ct_sectors <- ordered_sectors[sub("\\..*", "", ordered_sectors) == ct]
  if (length(ct_sectors) > 0) {
    highlight.sector(
      sector.index = ct_sectors,
      track.index = 1,
      col = adjustcolor(cell_colors[ct], alpha.f = 0.15),
      border = cell_colors[ct],
      lwd = 2,
      text = ct,
      text.col = cell_colors[ct],
      cex = 1.3,
      font = 2,
      facing = "bending.outside",
      niceFacing = TRUE
    )
  }
}

circos.clear()

# Title
title(main = "Eosinophil-Macrophage Signaling Network",
      cex.main = 1.4, font.main = 2, line = 0)
mtext("Pathway-resolved ligand-receptor communication",
      side = 3, line = -0.8, cex = 0.9, col = "grey40")

# Legend: cell state
legend("bottomleft",
       legend = names(cell_colors),
       fill = cell_colors,
       border = "grey30",
       bty = "n", cex = 1.0,
       title = expression(bold("Cell State")),
       title.adj = 0,
       inset = c(0.02, 0.02))

# Legend: link direction
mtext("Link color = sender cell type; arrow = direction",
      side = 1, line = 0.3, cex = 0.8, col = "grey40")

dev.off()
cat(sprintf("  Saved: %s\n", pdf_out))

# =============================================================================
# 4. Also generate a simplified version without self-loops
# =============================================================================
cat("\nGenerating version without self-loops...\n")

# Remove within-cell-type interactions
chord_cross <- chord_links[chord_links$source_cell != sub("\\..*", "", chord_links$to), ]

# Rebuild compound sectors
all_cross_sectors <- unique(c(chord_cross$from, chord_cross$to))
ordered_cross <- ordered_sectors[ordered_sectors %in% all_cross_sectors]

grid_colors_cross <- grid_colors[ordered_cross]
gap_cross <- rep(1, length(ordered_cross))
for (i in seq_len(length(ordered_cross) - 1)) {
  ct_curr <- sub("\\..*", "", ordered_cross[i])
  ct_next <- sub("\\..*", "", ordered_cross[i + 1])
  if (ct_curr != ct_next) gap_cross[i] <- 8
}
gap_cross[length(ordered_cross)] <- 8

link_colors_cross <- sapply(seq_len(nrow(chord_cross)), function(i) {
  ct <- chord_cross$source_cell[i]
  adjustcolor(cell_colors[ct], alpha.f = 0.45)
})

pdf_out2 <- file.path(out_dir, "fig3A_chord_cellchat_style_crossonly.pdf")
cairo_pdf(pdf_out2, width = 11, height = 11)

par(mar = c(2, 2, 3, 2))
circos.clear()
circos.par(start.degree = 90, gap.after = gap_cross,
           track.margin = c(0.01, 0.01), cell.padding = c(0, 0, 0, 0))

chordDiagram(
  chord_cross[, c("from", "to", "value")],
  order = ordered_cross,
  grid.col = grid_colors_cross,
  col = link_colors_cross,
  transparency = 0,
  directional = 1,
  direction.type = c("diffHeight", "arrows"),
  link.arr.type = "big.arrow",
  link.arr.length = 0.05,
  link.arr.width = 0.03,
  link.lwd = 1.5,
  link.sort = TRUE,
  link.decreasing = TRUE,
  annotationTrack = "grid",
  preAllocateTracks = list(track.height = 0.12)
)

circos.track(track.index = 1, panel.fun = function(x, y) {
  sector <- CELL_META$sector.index
  pw_name <- sub(".*\\.", "", sector)
  circos.text(mean(CELL_META$xlim), CELL_META$ylim[1] + 0.15, pw_name,
              facing = "clockwise", niceFacing = TRUE,
              cex = 0.7, font = 3, adj = c(0, 0.5))
}, bg.border = NA)

for (ct in cell_order) {
  ct_sectors <- ordered_cross[sub("\\..*", "", ordered_cross) == ct]
  if (length(ct_sectors) > 0) {
    highlight.sector(
      sector.index = ct_sectors,
      track.index = 1,
      col = adjustcolor(cell_colors[ct], alpha.f = 0.15),
      border = cell_colors[ct],
      lwd = 2,
      text = ct,
      text.col = cell_colors[ct],
      cex = 1.3,
      font = 2,
      facing = "bending.outside",
      niceFacing = TRUE
    )
  }
}

circos.clear()

title(main = "Eosinophil-Macrophage Signaling Network",
      cex.main = 1.4, font.main = 2, line = 0)
mtext("Cross-type communication only (Eos \u2194 Mac)",
      side = 3, line = -0.8, cex = 0.9, col = "grey40")

legend("bottomleft",
       legend = names(cell_colors),
       fill = cell_colors,
       border = "grey30",
       bty = "n", cex = 1.0,
       title = expression(bold("Cell State")),
       title.adj = 0,
       inset = c(0.02, 0.02))

mtext("Link color = sender; arrow = direction",
      side = 1, line = 0.3, cex = 0.8, col = "grey40")

dev.off()
cat(sprintf("  Saved: %s\n", pdf_out2))

# =============================================================================
# 5. Summary
# =============================================================================
cat("\n========================================\n")
cat("Summary:\n")
cat(sprintf("  fig3A_chord_cellchat_style.pdf          — All interactions\n"))
cat(sprintf("  fig3A_chord_cellchat_style_crossonly.pdf — Cross-type only (Eos<->Mac)\n"))
cat(sprintf("  Total links: %d (all), %d (cross-type)\n", nrow(chord_links), nrow(chord_cross)))
cat("\n44b complete.\n")
