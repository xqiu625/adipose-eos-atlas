#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=6:00:00
#SBATCH --job-name=cellchat
#SBATCH --output=cellchat_%j.out

# =============================================================================
# 10_cellchat.R
# CellChat analysis of eosinophil-macrophage communication across 6 conditions
# Input:  $BASE/eos_mac_combined.rds  (celltype: Eos_C0-C6, Mac_M1, Mac_M2_Surf, Mac_M2_Func)
# Output: Per-condition CellChat objects, merged object, interaction counts, PDFs
# =============================================================================

suppressPackageStartupMessages({
  library(CellChat)
  library(Seurat)
  library(patchwork)
  library(ggplot2)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
source(file.path(SCRIPTS, "00_theme.R"))

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
outdir <- file.path(BASE, "cellchat")
dir.create(outdir, recursive = TRUE, showWarnings = FALSE)

conditions <- c("WT_Female_HFD", "WT_Male_HFD", "KO_Female_HFD",
                "KO_Male_HFD", "WT_Female_HFDInfected", "KO_Female_HFDInfected")

# --- Load data ---------------------------------------------------------------
cat("Loading eos_mac_combined.rds ...\n")
obj <- tryCatch(
  readRDS(file.path(BASE, "eos_mac_combined.rds")),
  error = function(e) { cat("ERROR: Cannot load input:", e$message, "\n"); quit(status = 1) }
)
cat("  Cells:", ncol(obj), " Celltypes:", paste(unique(obj$celltype), collapse = ", "), "\n")

# --- Helper: classify cell as eos or mac -------------------------------------
is_eos <- function(ct) grepl("^Eos_", ct)
is_mac <- function(ct) grepl("^Mac_", ct)

# --- Run CellChat per condition ----------------------------------------------
cellchat_list <- list()
interaction_summary <- data.frame(condition = character(),
                                  n_eos = integer(),
                                  n_mac = integer(),
                                  n_interactions = integer(),
                                  n_pathways = integer(),
                                  stringsAsFactors = FALSE)

for (cond in conditions) {
  cat("\n========== Condition:", cond, "==========\n")

  tryCatch({
    # Subset cells for this condition
    cells_cond <- colnames(obj)[obj$condition == cond]
    if (length(cells_cond) == 0) {
      cat("  No cells for condition", cond, "- skipping\n")
      next
    }
    sub <- subset(obj, cells = cells_cond)

    n_eos <- sum(is_eos(sub$celltype))
    n_mac <- sum(is_mac(sub$celltype))
    cat("  Eos:", n_eos, " Mac:", n_mac, "\n")

    if (n_eos < 50 || n_mac < 50) {
      cat("  Fewer than 50 eos or 50 mac - skipping\n")
      next
    }

    # Prepare CellChat input
    data_input <- GetAssayData(sub, slot = "data")  # normalized
    meta <- data.frame(labels = sub$celltype, row.names = colnames(sub))

    cc <- createCellChat(object = data_input, meta = meta, group.by = "labels")
    cat("  CellChat object created\n")

    # Set ligand-receptor database
    CellChatDB <- CellChatDB.mouse
    cc@DB <- CellChatDB
    cat("  Database set (CellChatDB.mouse)\n")

    # Preprocessing
    cc <- subsetData(cc)
    cat("  Identifying overexpressed genes ...\n")
    cc <- identifyOverExpressedGenes(cc)
    cc <- identifyOverExpressedInteractions(cc)
    cat("  Overexpressed interactions identified\n")

    # Communication inference
    cat("  Computing communication probabilities (triMean) ...\n")
    cc <- computeCommunProb(cc, type = "triMean")
    cc <- filterCommunication(cc, min.cells = 10)
    cat("  Computing pathway-level probabilities ...\n")
    cc <- computeCommunProbPathway(cc)
    cc <- aggregateNet(cc)
    cat("  Computing network centrality ...\n")
    cc <- netAnalysis_computeCentrality(cc)

    # Record summary
    n_int <- nrow(subsetCommunication(cc))
    pathways_detected <- cc@netP$pathways
    n_path <- length(pathways_detected)
    cat("  Interactions:", n_int, " Pathways:", n_path, "\n")

    interaction_summary <- rbind(interaction_summary,
      data.frame(condition = cond, n_eos = n_eos, n_mac = n_mac,
                 n_interactions = n_int, n_pathways = n_path,
                 stringsAsFactors = FALSE))

    # Save individual object
    saveRDS(cc, file.path(outdir, paste0("cellchat_", cond, ".rds")))
    cat("  Saved cellchat_", cond, ".rds\n")

    cellchat_list[[cond]] <- cc

    # --- Per-condition figures ------------------------------------------------
    pdf_file <- file.path(outdir, paste0("cellchat_", cond, ".pdf"))
    pdf(pdf_file, width = 10, height = 8)

    # Interaction network circle plot
    tryCatch({
      netVisual_circle(cc@net$count,
                       vertex.weight = table(cc@idents),
                       weight.scale = TRUE,
                       title.name = paste0("Interaction count - ", cond))
    }, error = function(e) cat("  Warning: circle plot failed:", e$message, "\n"))

    # Outgoing signaling heatmap
    tryCatch({
      netAnalysis_signalingRole_heatmap(cc, pattern = "outgoing",
                                        title = paste0("Outgoing - ", cond))
    }, error = function(e) cat("  Warning: outgoing heatmap failed:", e$message, "\n"))

    # Incoming signaling heatmap
    tryCatch({
      netAnalysis_signalingRole_heatmap(cc, pattern = "incoming",
                                        title = paste0("Incoming - ", cond))
    }, error = function(e) cat("  Warning: incoming heatmap failed:", e$message, "\n"))

    # CCL pathway chord diagram (if detected)
    if ("CCL" %in% pathways_detected) {
      tryCatch({
        netVisual_chord_gene(cc, signaling = "CCL",
                             title.name = paste0("CCL pathway - ", cond))
      }, error = function(e) cat("  Warning: CCL chord failed:", e$message, "\n"))
    } else {
      cat("  CCL pathway not detected for", cond, "\n")
    }

    dev.off()
    cat("  Saved", pdf_file, "\n")

  }, error = function(e) {
    cat("  ERROR processing", cond, ":", e$message, "\n")
  })
}

# --- Merge CellChat objects ---------------------------------------------------
cat("\n========== Merging CellChat objects ==========\n")
if (length(cellchat_list) >= 2) {
  tryCatch({
    merged <- mergeCellChat(cellchat_list, add.names = names(cellchat_list))
    saveRDS(merged, file.path(outdir, "cellchat_merged.rds"))
    cat("Saved cellchat_merged.rds with", length(cellchat_list), "conditions\n")
  }, error = function(e) {
    cat("ERROR merging CellChat objects:", e$message, "\n")
  })
} else {
  cat("Fewer than 2 conditions passed filters - skipping merge\n")
}

# --- Save interaction summary -------------------------------------------------
write.csv(interaction_summary, file.path(outdir, "interaction_counts.csv"),
          row.names = FALSE)
cat("\nSaved interaction_counts.csv\n")

cat("\n========== CellChat analysis complete ==========\n")
cat("Output directory:", outdir, "\n")
sessionInfo()
