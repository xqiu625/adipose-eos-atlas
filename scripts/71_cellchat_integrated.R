#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 71_cellchat_integrated.R   (INTEGRATION-based CellChat — Fig 3)
# Eosinophil-macrophage communication on the ComBat-integrated clusters.
# Combine integrated eosinophils (i0-i5) with macrophage subtypes (M1/M2 Surf/
# M2 Func), run CellChat per condition, focusing on the CCL pathway (the
# manuscript's CCL6-CCR3/CCR2 finding). Reports per-condition CCL communication
# so we can compare to the published (un-integrated) result.
# Inputs : $REDO/eos_combat010.rds (i0-i5), mac_subtyped.rds (mac_subtype)
# Output : $REDO/cellchat/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(CellChat); library(dplyr) })
options(stringsAsFactors = FALSE)

REDO <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/revision/redo"
BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out  <- file.path(REDO, "cellchat"); dir.create(out, recursive=TRUE, showWarnings=FALSE)

# ── combine integrated eos (i0-i5) + macrophage subtypes ─────────────────────
eos <- readRDS(file.path(REDO, "eos_combat010.rds")); DefaultAssay(eos) <- "RNA"
eos$cc_celltype <- as.character(eos$new_cluster)            # i0..i5
mac <- readRDS(file.path(BASE, "mac_subtyped.rds"))
if (inherits(mac[["RNA"]],"Assay5")) mac[["RNA"]] <- JoinLayers(mac[["RNA"]])
mac$cc_celltype <- as.character(mac$mac_subtype)            # Mac_M1 / Mac_M2_Surf / Mac_M2_Func
# common genes; align condition label
genes <- intersect(rownames(eos), rownames(mac))
eos_c <- subset(eos, features=genes); mac_c <- subset(mac, features=genes)
mac_c$condition <- sub("_\\d+$", "", mac_c$orig.ident)
eos_c$condition <- sub("_\\d+$", "", eos_c$orig.ident)
comb <- merge(eos_c, mac_c)
if (inherits(comb[["RNA"]],"Assay5")) comb[["RNA"]] <- JoinLayers(comb[["RNA"]])
comb <- NormalizeData(comb, verbose=FALSE)
comb$cc_celltype <- factor(comb$cc_celltype)
cat("Combined object:", ncol(comb), "cells; groups:\n"); print(table(comb$cc_celltype))
cat("Conditions:\n"); print(table(comb$condition))

# ── CellChat per condition ───────────────────────────────────────────────────
run_cc <- function(sub) {
  data <- GetAssayData(sub, layer="data")
  grp  <- droplevels(factor(sub$cc_celltype))          # drop empty groups per condition
  keepg <- names(which(table(grp) >= 10))              # CellChat needs >=10 cells/group
  cells <- colnames(sub)[as.character(grp) %in% keepg]
  data <- data[, cells]; grp <- droplevels(grp[cells])
  meta <- data.frame(group = grp, row.names = cells)
  cc <- createCellChat(object = data, meta = meta, group.by = "group")
  cc@DB <- CellChatDB.mouse
  cc <- subsetData(cc)
  cc <- identifyOverExpressedGenes(cc)
  cc <- identifyOverExpressedInteractions(cc)
  cc <- computeCommunProb(cc, type = "triMean")
  cc <- filterCommunication(cc, min.cells = 10)
  cc <- computeCommunProbPathway(cc)
  cc <- aggregateNet(cc)
  cc
}
cc_list <- list()
for (cond in names(table(comb$condition))) {
  sub <- subset(comb, subset = condition == cond)
  if (ncol(sub) < 50 || length(unique(sub$cc_celltype)) < 3) { cat("skip", cond, "\n"); next }
  cat("\n=== CellChat:", cond, "(", ncol(sub), "cells ) ===\n")
  cc <- tryCatch(run_cc(sub), error=function(e){cat("  error:",conditionMessage(e),"\n");NULL})
  if (is.null(cc)) next
  cc_list[[cond]] <- cc
  paths <- cc@netP$pathways
  cat("  pathways:", paste(head(paths,20), collapse=", "), "\n")
  if ("CCL" %in% paths) {
    df <- subsetCommunication(cc, signaling = "CCL")
    write.csv(df, file.path(out, paste0("CCL_", cond, ".csv")), row.names=FALSE)
    cat("  CCL interactions:", nrow(df), "\n")
  }
}
saveRDS(cc_list, file.path(out, "cellchat_integrated_list.rds"))
cat("\nSaved cellchat_integrated_list.rds (", length(cc_list), "conditions ). Done 71.\n")
