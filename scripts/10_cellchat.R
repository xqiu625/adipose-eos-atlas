#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

#SBATCH --partition=batch
#SBATCH --mem=64G
#SBATCH --cpus-per-task=4
#SBATCH --time=10:00:00
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
cat("  Cells:", ncol(obj), " Celltypes:", paste(sort(unique(obj$detailed_celltype)), collapse = ", "), "\n")

# Seurat v5: join split layers so GetAssayData works
cat("Joining layers...\n")
obj <- JoinLayers(obj)

# --- Helper: classify cell as eos or mac -------------------------------------
is_eos <- function(ct) grepl("^Eos_", ct)
is_mac <- function(ct) grepl("^Mac_", ct)

# --- Run CellChat per condition ----------------------------------------------
# NOTE: each object below pools ALL biological replicates of a condition into
# one network, so any pathway/network-level comparison drawn from these
# objects (including the merged object below) is a descriptive comparison of
# two point estimates, not a test across replicates. A genuine replicate-level
# statistical comparison (per-sample CellChat + Wilcoxon test) is added further
# down and should be what's cited for "does communication differ" claims; use
# these pooled networks only for which pathways/genes are involved.
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

    n_eos <- sum(is_eos(sub$detailed_celltype))
    n_mac <- sum(is_mac(sub$detailed_celltype))
    cat("  Eos:", n_eos, " Mac:", n_mac, "\n")

    if (n_eos < 50 || n_mac < 50) {
      cat("  Fewer than 50 eos or 50 mac - skipping\n")
      next
    }

    # Prepare CellChat input
    data_input <- GetAssayData(sub, layer = "data")  # normalized
    meta <- data.frame(labels = sub$detailed_celltype, row.names = colnames(sub))

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

# --- Replicate-level statistical test (per-sample CellChat) -------------------
# Runs CellChat once per biological replicate (orig.ident) instead of pooling
# by condition, then compares conditions on the resulting per-sample summary
# with a Wilcoxon rank-sum test -- giving actual n=samples-per-condition
# inference for "does overall Eos<->Mac communication strength differ between
# conditions," which the pooled per-condition networks above cannot provide.
cat("\n========== Replicate-level CellChat (per sample) ==========\n")

samples_all <- unique(obj$orig.ident)
sample_summary <- data.frame(sample = character(), condition = character(),
                             n_eos = integer(), n_mac = integer(),
                             total_interactions = integer(),
                             total_strength = numeric(),
                             stringsAsFactors = FALSE)

for (s in samples_all) {
  cells_s <- colnames(obj)[obj$orig.ident == s]
  sub_s <- subset(obj, cells = cells_s)
  cond_s <- as.character(unique(sub_s$condition))[1]

  n_eos_s <- sum(is_eos(sub_s$detailed_celltype))
  n_mac_s <- sum(is_mac(sub_s$detailed_celltype))

  if (n_eos_s < 10 || n_mac_s < 10) {
    cat(sprintf("  %s (%s): fewer than 10 eos or 10 mac (eos=%d, mac=%d) - skipping\n",
                s, cond_s, n_eos_s, n_mac_s))
    next
  }

  cc_s <- tryCatch({
    data_input_s <- GetAssayData(sub_s, layer = "data")
    meta_s <- data.frame(labels = sub_s$detailed_celltype, row.names = colnames(sub_s))
    cc_tmp <- createCellChat(object = data_input_s, meta = meta_s, group.by = "labels")
    cc_tmp@DB <- CellChatDB.mouse
    cc_tmp <- subsetData(cc_tmp)
    cc_tmp <- identifyOverExpressedGenes(cc_tmp)
    cc_tmp <- identifyOverExpressedInteractions(cc_tmp)
    cc_tmp <- computeCommunProb(cc_tmp, type = "triMean")
    cc_tmp <- filterCommunication(cc_tmp, min.cells = 10)
    cc_tmp <- aggregateNet(cc_tmp)
    cc_tmp
  }, error = function(e) {
    cat(sprintf("  %s (%s): ERROR - %s\n", s, cond_s, e$message))
    NULL
  })

  if (is.null(cc_s)) next

  total_int <- sum(cc_s@net$count)
  total_str <- sum(cc_s@net$weight)
  cat(sprintf("  %s (%s): eos=%d, mac=%d, interactions=%d, strength=%.4f\n",
              s, cond_s, n_eos_s, n_mac_s, total_int, total_str))

  sample_summary <- rbind(sample_summary, data.frame(
    sample = s, condition = cond_s, n_eos = n_eos_s, n_mac = n_mac_s,
    total_interactions = total_int, total_strength = total_str,
    stringsAsFactors = FALSE
  ))
}

write.csv(sample_summary, file.path(outdir, "cellchat_replicate_summary.csv"), row.names = FALSE)
cat("Saved cellchat_replicate_summary.csv\n")

# Wilcoxon rank-sum test: each condition vs WT_Female_HFD (reference)
ref_cond <- "WT_Female_HFD"
replicate_tests <- data.frame(condition = character(), n_ref = integer(),
                              n_test = integer(), p_value_count = numeric(),
                              p_value_strength = numeric(),
                              stringsAsFactors = FALSE)

if (nrow(sample_summary) > 0 && ref_cond %in% sample_summary$condition) {
  ref_vals <- sample_summary[sample_summary$condition == ref_cond, ]
  for (cond in setdiff(unique(sample_summary$condition), ref_cond)) {
    test_vals <- sample_summary[sample_summary$condition == cond, ]
    if (nrow(ref_vals) < 2 || nrow(test_vals) < 2) {
      cat(sprintf("  %s vs %s: fewer than 2 replicates on one side - skipping test\n",
                  cond, ref_cond))
      next
    }
    p_count <- tryCatch(wilcox.test(test_vals$total_interactions, ref_vals$total_interactions)$p.value,
                        error = function(e) NA)
    p_strength <- tryCatch(wilcox.test(test_vals$total_strength, ref_vals$total_strength)$p.value,
                           error = function(e) NA)
    replicate_tests <- rbind(replicate_tests, data.frame(
      condition = cond, n_ref = nrow(ref_vals), n_test = nrow(test_vals),
      p_value_count = p_count, p_value_strength = p_strength,
      stringsAsFactors = FALSE
    ))
  }
} else {
  cat(sprintf("  Reference condition %s has no replicate summary - skipping tests\n", ref_cond))
}

if (nrow(replicate_tests) > 0) {
  replicate_tests$fdr_count    <- p.adjust(replicate_tests$p_value_count, method = "BH")
  replicate_tests$fdr_strength <- p.adjust(replicate_tests$p_value_strength, method = "BH")
}
write.csv(replicate_tests, file.path(outdir, "cellchat_replicate_tests.csv"), row.names = FALSE)
cat("Saved cellchat_replicate_tests.csv (Wilcoxon, each condition vs WT_Female_HFD, BH-corrected)\n")

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
