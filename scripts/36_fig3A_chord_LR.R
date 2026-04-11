#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 36_fig3A_chord_LR.R
# Fig 3A alternative: Chord diagram with individual CCL L-R pair labels
# (Reviewer Comment 9: "label Fig 3A with specific ligand/receptor pairs")
#
# Sectors = individual CCL ligands & receptors (gene names on circumference)
# Sector arcs = proportional cell type group coloring (shows all Mac subtypes)
# Links = L-R pair interactions, width proportional to communication probability
# Link colors = direction (warm = Eos->Mac, cool = Mac->Eos)
#
# Input:  $BASE/cellchat/cellchat_<condition>.rds (6 conditions)
# Output: revision/figures/fig3/fig3A_chord_LR_pairs.pdf
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

# Direction colors for links (2-color: warm/cool)
dir_colors <- c("Eos_to_Mac" = "#E31A1C",   # warm red
                "Mac_to_Eos" = "#2171B5")    # cool blue

# =============================================================================
# 1. Extract CCL L-R pair interactions from all conditions
# =============================================================================
cat("Loading CellChat objects (raw deserialization)...\n")

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
# 2. Map cell types to groups and standardize molecule names
# =============================================================================
comm_df$source_group <- chord_group(comm_df$source)
comm_df$target_group <- chord_group(comm_df$target)

# Keep only Eos <-> Mac interactions (exclude within-group)
comm_df <- comm_df[!is.na(comm_df$source_group) & !is.na(comm_df$target_group), ]
is_eos_mac <- (comm_df$source_group == "Eosinophils" & comm_df$target_group != "Eosinophils") |
              (comm_df$source_group != "Eosinophils" & comm_df$target_group == "Eosinophils")
comm_df <- comm_df[is_eos_mac, ]

# Capitalize molecule names: Ccl6, Ccr2, etc.
cap_name <- function(x) {
  paste0(toupper(substr(x, 1, 1)), tolower(substr(x, 2, nchar(x))))
}
comm_df$ligand_cap <- cap_name(comm_df$ligand)
comm_df$receptor_cap <- cap_name(comm_df$receptor)

# Direction label (2-direction: Eos_to_Mac / Mac_to_Eos)
comm_df$direction <- ifelse(comm_df$source_group == "Eosinophils",
                            "Eos_to_Mac", "Mac_to_Eos")

cat("Ligands:", paste(sort(unique(comm_df$ligand_cap)), collapse = ", "), "\n")
cat("Receptors:", paste(sort(unique(comm_df$receptor_cap)), collapse = ", "), "\n")

# =============================================================================
# 3. Aggregate: average across conditions, then sum by direction per L-R pair
# =============================================================================

# First average across conditions for each (source, target, ligand, receptor)
agg_cond <- aggregate(prob ~ source + target + ligand_cap + receptor_cap +
                        source_group + target_group + direction,
                      data = comm_df, FUN = mean)

# Per L-R pair per direction: sum across all cell type pairs
agg_lr <- aggregate(prob ~ ligand_cap + receptor_cap + direction,
                    data = agg_cond, FUN = sum)

cat("\nL-R pair strengths by direction:\n")
agg_lr <- agg_lr[order(-agg_lr$prob), ]
print(agg_lr)

# Total strength per L-R pair (for ranking)
lr_total <- aggregate(prob ~ ligand_cap + receptor_cap, data = agg_lr, FUN = sum)
lr_total <- lr_total[order(-lr_total$prob), ]
cat("\nTotal L-R pair strengths:\n")
print(lr_total)

# =============================================================================
# 4. Proportional cell-type composition per molecule sector
# =============================================================================

ligands  <- sort(unique(agg_lr$ligand_cap))
receptors <- sort(unique(agg_lr$receptor_cap))
all_molecules <- c(ligands, receptors)

# Ligands: fraction of total signal sent by each cell type group
lig_grp <- aggregate(prob ~ ligand_cap + source_group, data = agg_cond, FUN = sum)
mol_props <- list()

for (mol in ligands) {
  sub <- lig_grp[lig_grp$ligand_cap == mol, ]
  total <- sum(sub$prob)
  all_g <- data.frame(group = names(group_colors), stringsAsFactors = FALSE)
  sub <- merge(all_g, sub, by.x = "group", by.y = "source_group", all.x = TRUE)
  sub$prob[is.na(sub$prob)] <- 0
  sub$frac <- if (total > 0) sub$prob / total else 0
  sub <- sub[order(match(sub$group, names(group_colors))), ]
  mol_props[[mol]] <- sub
}

# Receptors: fraction of total signal received by each cell type group
rec_grp <- aggregate(prob ~ receptor_cap + target_group, data = agg_cond, FUN = sum)

for (mol in receptors) {
  sub <- rec_grp[rec_grp$receptor_cap == mol, ]
  total <- sum(sub$prob)
  all_g <- data.frame(group = names(group_colors), stringsAsFactors = FALSE)
  sub <- merge(all_g, sub, by.x = "group", by.y = "target_group", all.x = TRUE)
  sub$prob[is.na(sub$prob)] <- 0
  sub$frac <- if (total > 0) sub$prob / total else 0
  sub <- sub[order(match(sub$group, names(group_colors))), ]
  mol_props[[mol]] <- sub
}

cat("\nProportional cell-type composition per molecule:\n")
for (mol in all_molecules) {
  props <- mol_props[[mol]]
  nz <- props[props$frac > 0, ]
  props_str <- paste(sprintf("%s=%.0f%%", nz$group, nz$frac * 100), collapse = ", ")
  cat(sprintf("  %s: %s\n", mol, props_str))
}

# Dominant group color (fallback for grid.col)
sector_colors <- character(length(all_molecules))
names(sector_colors) <- all_molecules
for (mol in all_molecules) {
  props <- mol_props[[mol]]
  sector_colors[mol] <- group_colors[props$group[which.max(props$frac)]]
}

# =============================================================================
# 5. Build chord diagram data frame with directional coloring
# =============================================================================

chord_df <- data.frame(
  from      = agg_lr$ligand_cap,
  to        = agg_lr$receptor_cap,
  value     = agg_lr$prob,
  direction = agg_lr$direction,
  stringsAsFactors = FALSE
)

# Link colors by direction (2-color: red/blue)
link_cols <- ifelse(chord_df$direction == "Eos_to_Mac",
                    adjustcolor(dir_colors["Eos_to_Mac"], alpha.f = 0.55),
                    adjustcolor(dir_colors["Mac_to_Eos"], alpha.f = 0.55))

# Sector order: ligands on left, receptors on right, big gap between groups
sector_order <- c(ligands, receptors)
n_lig <- length(ligands)
n_rec <- length(receptors)

gap_sizes <- rep(3, length(sector_order))
# Big gaps to separate ligands from receptors
if (n_lig > 0) gap_sizes[n_lig] <- 20
if (n_rec > 0) gap_sizes[n_lig + n_rec] <- 20

cat(sprintf("\n%d ligand sectors, %d receptor sectors\n", n_lig, n_rec))

# =============================================================================
# 6. Plot chord diagram with proportional sector arcs
# =============================================================================
cat("\nGenerating chord diagram...\n")

pdf_out <- file.path(out_dir, "fig3A_chord_LR_pairs.pdf")
cairo_pdf(pdf_out, width = 10, height = 10)

circos.clear()
circos.par(start.degree = 90, gap.after = gap_sizes,
           track.margin = c(0.01, 0.01))

chordDiagram(chord_df[, c("from", "to", "value")],
             order = sector_order,
             grid.col = sector_colors,
             col = link_cols,
             transparency = 0,
             directional = 1,
             direction.type = c("diffHeight", "arrows"),
             link.arr.type = "big.arrow",
             link.arr.length = 0.05,
             link.arr.width = 0.03,
             link.lwd = 1.5,
             link.sort = TRUE,
             link.decreasing = TRUE,
             annotationTrack = NULL,
             preAllocateTracks = list(
               list(track.height = 0.18),
               list(track.height = 0.06)
             ))

# Track 2 (inner): proportional cell-type colors for each sector arc
circos.track(track.index = 2, panel.fun = function(x, y) {
  sector <- CELL_META$sector.index
  xlim <- CELL_META$xlim
  ylim <- CELL_META$ylim

  props <- mol_props[[sector]]
  x_start <- xlim[1]
  x_range <- xlim[2] - xlim[1]
  for (i in seq_len(nrow(props))) {
    if (props$frac[i] == 0) next
    x_end <- x_start + x_range * props$frac[i]
    circos.rect(x_start, ylim[1], x_end, ylim[2],
                col = group_colors[props$group[i]], border = NA)
    x_start <- x_end
  }
}, bg.border = "grey60")

# Track 1 (outer): gene name labels in italic
circos.track(track.index = 1, panel.fun = function(x, y) {
  circos.text(CELL_META$xcenter, CELL_META$ylim[1], CELL_META$sector.index,
              facing = "clockwise", niceFacing = TRUE, adj = c(0, 0.5),
              cex = 1.2, font = 3)
}, bg.border = NA)

# Group labels: "Ligands" and "Receptors"
if (n_lig > 1) {
  highlight.sector(ligands, track.index = 1,
                   col = NA, border = NA,
                   text = "Ligands", text.col = "grey30",
                   cex = 1.3, facing = "bending.inside",
                   niceFacing = TRUE, padding = c(0, 0, 0, 0))
}
if (n_rec > 1) {
  highlight.sector(receptors, track.index = 1,
                   col = NA, border = NA,
                   text = "Receptors", text.col = "grey30",
                   cex = 1.3, facing = "bending.inside",
                   niceFacing = TRUE, padding = c(0, 0, 0, 0))
}

circos.clear()

# Legend: direction + cell type group
legend("bottomleft",
       legend = c("Eos \u2192 Mac signaling",
                  "Mac \u2192 Eos signaling",
                  "", names(group_colors)),
       fill = c(dir_colors["Eos_to_Mac"], dir_colors["Mac_to_Eos"],
                NA, group_colors),
       border = c("black", "black", NA, rep("black", length(group_colors))),
       bty = "n", cex = 0.85,
       title = "CCL Signaling", title.adj = 0)

dev.off()
cat(sprintf("  Saved: %s\n", pdf_out))

# =============================================================================
# 7. Summary
# =============================================================================
cat("\n========== Summary ==========\n")
cat(sprintf("Sectors: %d ligands + %d receptors = %d total\n",
            n_lig, n_rec, length(all_molecules)))
cat(sprintf("Links: %d directional L-R pair connections\n", nrow(chord_df)))

cat("\nTop L-R pairs:\n")
for (i in seq_len(min(10, nrow(lr_total)))) {
  lig <- lr_total$ligand_cap[i]
  rec <- lr_total$receptor_cap[i]
  cat(sprintf("  %s -> %s : %.4f\n", lig, rec, lr_total$prob[i]))
}

cat(sprintf("\nOutput: %s\n", pdf_out))
cat("36 complete.\n")
