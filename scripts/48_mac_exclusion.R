#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))
if (nzchar(Sys.getenv("RERUN_R_LIB"))) {
  .rerun_r_lib_resolved <- normalizePath(.libPaths()[1], mustWork = FALSE)
  .rerun_r_lib_expected <- normalizePath(Sys.getenv("RERUN_R_LIB"), mustWork = FALSE)
  if (!identical(.rerun_r_lib_resolved, .rerun_r_lib_expected)) {
    stop(sprintf(
      "RERUN_R_LIB audit guard: .libPaths()[1] resolved to '%s' but RERUN_R_LIB='%s' -- the shared/legacy library would take precedence over the pinned isolated library. Refusing to proceed.",
      .rerun_r_lib_resolved, .rerun_r_lib_expected))
  }
  cat("=== Pinned R environment audit guard ===\n")
  cat("R.version.string:", R.version.string, "\n")
  cat("RERUN_R_LIB:", Sys.getenv("RERUN_R_LIB"), "\n")
  cat("CMDSTAN:", Sys.getenv("CMDSTAN"), "\n")
  cat(".libPaths():\n")
  print(.libPaths())
  cat("=========================================\n")
}

# ==============================================================================
# 48_mac_exclusion.R   (Reviewer #1, Concern #2 — complements Rebecca's FACS)
# Show the eosinophil clusters are NOT misclassified SiglecF+ macrophages:
#   - eos clusters are negative for macrophage markers (Mertk/Fcgr1/Adgre1/Csf1r/
#     Lyz2/C1qa/Cd68) and positive for eosinophil markers (Siglecf/Il5ra/Ccr3/Epx/
#     Prg2/Ear genes); the macrophages in the SAME data form a separate, marker+
#     cluster — a misclassified macrophage would land there, but none do.
#   - quantify cells co-expressing Siglecf AND a macrophage marker.
# Input : $BASE/all_cells_annotated.rds (eos + mac + mast + others)
# Output: $BASE/revision/figures/cluster_robustness/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr) })

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out  <- file.path(BASE, "revision", "figures", "cluster_robustness")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

obj_path <- file.path(BASE, "all_cells_annotated.rds")
if (!file.exists(obj_path)) obj_path <- file.path(BASE, "eos_mac_combined.rds")
cat("Loading", obj_path, "\n")
o <- readRDS(obj_path)
if (inherits(o[["RNA"]], "Assay5")) o[["RNA"]] <- JoinLayers(o[["RNA"]])

# pick the cell-type column
ctcol <- intersect(c("celltype","broad_celltype","cell_type","annotation"), colnames(o@meta.data))[1]
if (is.na(ctcol)) { o$celltype <- as.character(Idents(o)); ctcol <- "celltype" }
Idents(o) <- o[[ctcol]][,1]
cat("Cell-type column:", ctcol, "\n"); print(table(Idents(o)))

mac_markers <- intersect(c("Mertk","Fcgr1","Adgre1","Csf1r","Lyz2","C1qa","Cd68","Itgam"), rownames(o))
eos_markers <- intersect(c("Siglecf","Il5ra","Ccr3","Epx","Prg2","Ear2","Ear6","Il4","Ccr3"), rownames(o))
cat("Mac markers:", paste(mac_markers, collapse=","), "\nEos markers:", paste(eos_markers, collapse=","), "\n")

# ── Dot plot: markers x cell type ────────────────────────────────────────────
p <- DotPlot(o, features = c(eos_markers, mac_markers)) +
  RotatedAxis() +
  labs(title = "Eosinophil vs macrophage markers across cell types",
       subtitle = "Eosinophil clusters: eos-marker+ / mac-marker-negative") +
  theme(axis.text.x = element_text(size = 9))
ggsave(file.path(out, "mac_exclusion_dotplot.pdf"), p, width = 11, height = 6)

# ── Quantify mac-marker expression within eosinophils ────────────────────────
expr <- GetAssayData(o, layer = "data")
is_eos <- grepl("eos|Eos", as.character(Idents(o)))
tab <- lapply(c(eos_markers, mac_markers), function(g) {
  v <- expr[g, ]
  data.frame(gene = g,
             pct_pos_eos = round(100*mean(v[is_eos] > 0), 1),
             mean_eos = round(mean(v[is_eos]), 3),
             pct_pos_nonEos = round(100*mean(v[!is_eos] > 0), 1),
             mean_nonEos = round(mean(v[!is_eos]), 3))
}) %>% bind_rows()
write.csv(tab, file.path(out, "marker_expression_eos_vs_other.csv"), row.names = FALSE)
cat("\n=== marker expression: eosinophils vs other cells ===\n"); print(tab)

# ── SiglecF+ macrophage co-expression check ──────────────────────────────────
if ("Siglecf" %in% rownames(o) && length(mac_markers) > 0) {
  sig <- expr["Siglecf", ] > 0
  macpos <- Matrix::colSums(expr[mac_markers, , drop = FALSE] > 0) > 0
  co <- sig & macpos
  ct <- as.character(Idents(o))
  cocat <- table(ct[co])
  cat(sprintf("\nCells co-expressing Siglecf AND >=1 mac marker: %d (%.1f%% of all)\n",
              sum(co), 100*mean(co)))
  cat("  of those, by cell type (should be macrophages, NOT eos clusters):\n"); print(cocat)
  write.csv(data.frame(celltype = names(cocat), n_siglecf_macmarker = as.integer(cocat)),
            file.path(out, "siglecf_macmarker_coexpr_by_celltype.csv"), row.names = FALSE)
}

cat("\nDone 48_mac_exclusion.R. Outputs in", out, "\n")
