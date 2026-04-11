#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 40_fig3A_pathway_chord.R
# Fig 3A revised: MDPI-style chord diagram
# (Reviewer 2 Comment 9 + Rebecca's request)
#
# Sectors = cell states (Eos, M1 Mac, M2 Functional, M2 Surface)
# Sector arcs colored by cell state
# Links = pathway-level communication, width ~ signal strength
# Link colors = pathway identity
# Pathway names labeled on left margin
# Focus pathways: APP, CCL, FN1, CD45, THBS (+ others if significant)
#
# Input:  $BASE/cellchat/cellchat_<condition>.rds (6 conditions)
# Output: revision/figures/fig3/fig3A_pathway_chord.pdf
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

# Focus pathways (Rebecca's selection from Fig 3B)
FOCUS_PATHWAYS <- c("APP", "CCL", "FN1", "CD45", "THBS")

# Pathway link colors — distinct per pathway
pathway_colors <- c(
  "APP"  = "#E41A1C",
  "CCL"  = "#377EB8",
  "FN1"  = "#4DAF4A",
  "CD45" = "#FF7F00",
  "THBS" = "#984EA3"
)

# =============================================================================
# 1. Extract per-pathway probabilities across all conditions
# =============================================================================
cat("Loading CellChat objects and extracting pathway-level probabilities...\n")

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

  cat(sprintf("    %d pathways, %d cell types\n", length(pathways), length(ct_names)))
  cat(sprintf("    Available focus pathways: %s\n",
              paste(intersect(FOCUS_PATHWAYS, pathways), collapse = ", ")))

  for (pw in pathways) {
    prob_mat <- netP_prob[, , pw]
    keep <- !is.na(groups)

    for (i in which(keep)) {
      for (j in which(keep)) {
        val <- prob_mat[i, j]
        if (val > 0) {
          key <- paste(pw, groups[i], groups[j], sep = "|||")
          if (is.null(all_pathway_data[[key]])) {
            all_pathway_data[[key]] <- list(
              pathway = pw, source = groups[i], target = groups[j],
              vals = c()
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
cat(sprintf("\nTotal pathway-cell interactions: %d\n", nrow(pw_df)))

# Filter to focus pathways
pw_focus <- pw_df[pw_df$pathway %in% FOCUS_PATHWAYS, ]
cat(sprintf("Focus pathway interactions: %d\n", nrow(pw_focus)))

# Also find top non-focus pathways by total probability
pw_other <- pw_df[!pw_df$pathway %in% FOCUS_PATHWAYS, ]
pw_other_total <- aggregate(prob ~ pathway, data = pw_other, FUN = sum)
pw_other_total <- pw_other_total[order(-pw_other_total$prob), ]
cat("\nTop non-focus pathways:\n")
print(head(pw_other_total, 10))

# Print focus pathway summary
cat("\nFocus pathway summary:\n")
for (pw in FOCUS_PATHWAYS) {
  sub <- pw_focus[pw_focus$pathway == pw, ]
  if (nrow(sub) > 0) {
    cat(sprintf("  %s: %d interactions, total prob = %.4f\n",
                pw, nrow(sub), sum(sub$prob)))
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("    %s -> %s : %.4f\n", sub$source[i], sub$target[i], sub$prob[i]))
    }
  } else {
    cat(sprintf("  %s: NOT FOUND in any condition\n", pw))
  }
}

# Keep only Eos <-> Mac interactions (exclude within-group)
is_cross <- (pw_focus$source == "Eos" & pw_focus$target != "Eos") |
            (pw_focus$source != "Eos" & pw_focus$target == "Eos")
pw_cross <- pw_focus[is_cross, ]
cat(sprintf("\nCross-type (Eos<->Mac) focus interactions: %d\n", nrow(pw_cross)))

# =============================================================================
# 2. Build chord diagram data
# =============================================================================
cat("\nBuilding chord diagram...\n")

# Create a unique link for each pathway x direction
# Each link: from = source, to = target, value = prob, colored by pathway
chord_df <- data.frame(
  from    = pw_cross$source,
  to      = pw_cross$target,
  value   = pw_cross$prob,
  pathway = pw_cross$pathway,
  stringsAsFactors = FALSE
)

# Link colors by pathway
link_cols <- adjustcolor(pathway_colors[chord_df$pathway], alpha.f = 0.55)

# Sector order: Eos first, then Mac subtypes
sector_order <- c("Eos", "M1 Mac", "M2 Func", "M2 Surf")
gap_sizes <- c(8, 4, 4, 8)

# =============================================================================
# 3. Plot
# =============================================================================
pdf_out <- file.path(out_dir, "fig3A_pathway_chord.pdf")
cairo_pdf(pdf_out, width = 10, height = 10)

circos.clear()
circos.par(start.degree = 90, gap.after = gap_sizes,
           track.margin = c(0.01, 0.01))

chordDiagram(chord_df[, c("from", "to", "value")],
             order = sector_order,
             grid.col = cell_colors,
             col = link_cols,
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
             preAllocateTracks = list(track.height = 0.15))

# Sector labels: cell state names, bold
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
              cex = 1.3, font = 2)
}, bg.border = NA)

circos.clear()

# Legend: pathway colors
legend("bottomleft",
       legend = names(pathway_colors),
       fill = pathway_colors,
       border = "black",
       bty = "n", cex = 0.9,
       title = "Signaling Pathway", title.adj = 0)

# Legend: cell state colors
legend("bottomright",
       legend = names(cell_colors),
       fill = cell_colors,
       border = "black",
       bty = "n", cex = 0.9,
       title = "Cell State", title.adj = 0)

dev.off()
cat(sprintf("  Saved: %s\n", pdf_out))

# =============================================================================
# 4. Also generate per-pathway individual chord diagrams
# =============================================================================
cat("\nGenerating per-pathway chord diagrams...\n")

for (pw in FOCUS_PATHWAYS) {
  pw_sub <- pw_cross[pw_cross$pathway == pw, ]
  if (nrow(pw_sub) == 0) {
    cat(sprintf("  %s: no cross-type interactions, skipping\n", pw))
    next
  }

  pw_link_cols <- rep(adjustcolor(pathway_colors[pw], alpha.f = 0.55), nrow(pw_sub))

  pw_pdf <- file.path(out_dir, paste0("fig3A_pathway_chord_", pw, ".pdf"))
  cairo_pdf(pw_pdf, width = 8, height = 8)

  circos.clear()
  circos.par(start.degree = 90, gap.after = gap_sizes,
             track.margin = c(0.01, 0.01))

  chordDiagram(pw_sub[, c("source", "target", "prob")],
               order = sector_order,
               grid.col = cell_colors,
               col = pw_link_cols,
               transparency = 0,
               directional = 1,
               direction.type = c("diffHeight", "arrows"),
               link.arr.type = "big.arrow",
               link.arr.length = 0.06,
               link.arr.width = 0.04,
               link.lwd = 2,
               link.sort = TRUE,
               link.decreasing = TRUE,
               annotationTrack = "grid",
               preAllocateTracks = list(track.height = 0.15))

  circos.track(track.index = 1, panel.fun = function(x, y) {
    circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
                facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
                cex = 1.3, font = 2)
  }, bg.border = NA)

  circos.clear()

  title(main = paste0(pw, " Signaling"), font.main = 2, cex.main = 1.4)

  dev.off()
  cat(sprintf("  Saved: %s\n", pw_pdf))
}

# =============================================================================
# 5. Summary
# =============================================================================
cat("\n========== Summary ==========\n")
cat(sprintf("Focus pathways: %s\n", paste(FOCUS_PATHWAYS, collapse = ", ")))
cat(sprintf("Cross-type links: %d\n", nrow(pw_cross)))
cat(sprintf("Cell states: %s\n", paste(sector_order, collapse = ", ")))

cat("\nSignaling directions per pathway:\n")
for (pw in FOCUS_PATHWAYS) {
  sub <- pw_cross[pw_cross$pathway == pw, ]
  if (nrow(sub) > 0) {
    for (i in seq_len(nrow(sub))) {
      cat(sprintf("  %s: %s -> %s (%.4f)\n", pw, sub$source[i], sub$target[i], sub$prob[i]))
    }
  }
}

cat(sprintf("\nOutput: %s\n", out_dir))
cat("40 complete.\n")
