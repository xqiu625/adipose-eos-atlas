#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 48b_mac_exclusion_v2.R   (Reviewer #1, Concern #2 — corrected framing)
# Lead with the DISCRIMINATING macrophage markers MerTK (Mertk) and CD64 (Fcgr1):
# murine eosinophils also express F4/80 (Adgre1) and CD11b (Itgam), so those are
# NOT valid for excluding macrophages — we say so explicitly. SiglecF+ co-expression
# is recomputed using ONLY MerTK/CD64.
# Input : $BASE/all_cells_annotated.rds
# Output: $BASE/revision/figures/cluster_robustness/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(ggplot2); library(dplyr) })

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out  <- file.path(BASE, "revision", "figures", "cluster_robustness")
dir.create(out, recursive = TRUE, showWarnings = FALSE)

p <- file.path(BASE, "all_cells_annotated.rds")
o <- readRDS(p)
if (inherits(o[["RNA"]], "Assay5")) o[["RNA"]] <- JoinLayers(o[["RNA"]])
ctcol <- intersect(c("celltype","broad_celltype","cell_type","annotation"), colnames(o@meta.data))[1]
Idents(o) <- o[[ctcol]][,1]

disc_mac <- intersect(c("Mertk","Fcgr1"), rownames(o))               # discriminating
nondisc  <- intersect(c("Adgre1","Itgam","Csf1r","Lyz2","C1qa","Cd68"), rownames(o)) # shared/eos-expressed
eos_mk   <- intersect(c("Siglecf","Il5ra","Ccr3","Il4"), rownames(o))

# put cell types in a sensible order with Eosinophils/Macrophages first if present
lev <- levels(Idents(o)); ord <- c(grep("eos|Eos", lev, value=TRUE), grep("mac|Mac", lev, value=TRUE))
ord <- c(ord, setdiff(lev, ord)); Idents(o) <- factor(Idents(o), levels = ord)

# ── Panel A: discriminating markers ──────────────────────────────────────────
pA <- DotPlot(o, features = c(eos_mk, disc_mac)) + RotatedAxis() +
  labs(title = "Eosinophil markers + DISCRIMINATING macrophage markers (MerTK, CD64)",
       subtitle = "Eosinophil clusters are MerTK- / CD64- ; macrophages are MerTK+ / CD64+")
ggsave(file.path(out, "mac_exclusion_v2_discriminating.pdf"), pA, width = 9, height = 6)

# ── Panel B: non-discriminating markers (shown, labelled as such) ─────────────
pB <- DotPlot(o, features = nondisc) + RotatedAxis() +
  labs(title = "Markers NOT valid for macrophage exclusion in mouse",
       subtitle = "F4/80 (Adgre1) and CD11b (Itgam) are expressed by murine eosinophils too")
ggsave(file.path(out, "mac_exclusion_v2_nondiscriminating.pdf"), pB, width = 9, height = 5)

# ── Strict co-expression: SiglecF+ AND (MerTK or CD64) ───────────────────────
expr <- GetAssayData(o, layer = "data")
sig  <- expr["Siglecf", ] > 0
strict <- Matrix::colSums(expr[disc_mac, , drop = FALSE] > 0) > 0
co <- sig & strict
ct <- as.character(Idents(o))
cat(sprintf("\nCells co-expressing SiglecF AND (MerTK or CD64): %d (%.2f%% of all)\n",
            sum(co), 100*mean(co)))
tab <- as.data.frame(table(celltype = ct[co]))
# fraction WITHIN each cell type
n_by_ct <- table(ct)
tab$pct_within_celltype <- round(100*tab$Freq / as.integer(n_by_ct[as.character(tab$celltype)]), 2)
write.csv(tab, file.path(out, "siglecf_strict_macmarker_by_celltype.csv"), row.names = FALSE)
cat("\n=== SiglecF+ & (MerTK|CD64)+ cells, by cell type (eos should be ~0%) ===\n"); print(tab)

# ── Discriminating-marker expression table ───────────────────────────────────
is_eos <- grepl("eos|Eos", ct)
tab2 <- lapply(c(eos_mk, disc_mac, nondisc), function(g) {
  v <- expr[g, ]
  data.frame(gene = g, role = ifelse(g %in% disc_mac, "mac_discriminating",
                              ifelse(g %in% nondisc, "mac_nondiscriminating", "eosinophil")),
             pct_pos_eos = round(100*mean(v[is_eos]>0),1),
             pct_pos_nonEos = round(100*mean(v[!is_eos]>0),1))
}) %>% bind_rows()
write.csv(tab2, file.path(out, "marker_roles_eos_vs_other.csv"), row.names = FALSE)
cat("\n=== marker roles (eos vs other) ===\n"); print(tab2)

cat("\nDone 48b_mac_exclusion_v2.R\n")
