#!/usr/bin/env Rscript
lib_path <- "/bigdata/godziklab/shared/Xinru/R"
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 39_gurtner_intestinal_comparison.R
# Expand cross-dataset comparison to include Gurtner intestinal eosinophils
# (Reviewer 2 point 7: "compare intestinal eosinophils from Gurtner")
#
# Panels:
#   A: Integrated UMAP colored by tissue source (BM, Blood, Colon, SI, Adipose)
#   B: Marker dotplot across all tissue groups
#   C: Spearman correlation heatmap (Qiu C0-C5 vs BM/Blood/Colon/SI)
#   D: Dendrogram (Hu-paper style) showing cluster relationships
#
# Input:  $BASE/integration/qiu_gurtner_integrated.rds
# Output: revision/figures/supplement/supp_gurtner_intestinal_*.pdf
# SLURM:  batch, 64G, 4 CPUs, 2h
# ==============================================================================

suppressPackageStartupMessages({
  library(Seurat)
  library(ggplot2)
  library(dplyr)
  library(pheatmap)
  library(patchwork)
  library(dendextend)
})

SCRIPTS <- dirname(sub("--file=", "", grep("--file=", commandArgs(trailingOnly = FALSE), value = TRUE)[1]))
if (is.na(SCRIPTS) || SCRIPTS == "") {
  SCRIPTS <- "/bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-subtype-atlas/scripts"
}
source(file.path(SCRIPTS, "00_theme.R"))

BASE    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
PNAS    <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision"
out_dir <- file.path(PNAS, "revision", "figures", "supplement")
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

# =============================================================================
# 1. Load Gurtner integrated object
# =============================================================================
cat("Loading Qiu-Gurtner integrated object...\n")
combined <- readRDS(file.path(BASE, "integration", "qiu_gurtner_integrated.rds"))
if (inherits(combined[["RNA"]], "Assay5")) combined[["RNA"]] <- JoinLayers(combined[["RNA"]])
cat(sprintf("  %d cells\n", ncol(combined)))

# Check available organs
cat("\nAvailable organs in Gurtner data:\n")
print(table(combined$organ))

# =============================================================================
# 2. Create tissue groups including intestinal
# =============================================================================
cat("\nCreating tissue groups (including intestinal)...\n")

# Gurtner cluster label map (from 30_cross_dataset_similarity.R)
gurtner_label_map <- c(
  "4"  = "Precursors", "3"  = "Immature", "16" = "Cycling",
  "2"  = "Circulating", "1"  = "A-Eos", "5"  = "B-Eos"
)

combined$tissue_group <- case_when(
  combined$dataset == "Qiu_2025_Adipose" ~ paste0("Adipose_", combined$original_cluster),
  combined$organ == "Bone Marrow"         ~ "BM",
  combined$organ == "Blood"               ~ "Blood",
  combined$organ == "Colon"               ~ "Colon",
  combined$organ == "Small Intestine"     ~ "SI",
  combined$organ == "small intestine"     ~ "SI",
  grepl("[Cc]olon", combined$organ)       ~ "Colon",
  grepl("[Ss]mall|SI", combined$organ)    ~ "SI",
  TRUE                                    ~ paste0("Other_", combined$organ)
)

cat("  Group sizes:\n")
print(sort(table(combined$tissue_group), decreasing = TRUE))

# Define groups to keep
main_groups <- c(paste0("Adipose_C", 0:5), "BM", "Blood", "Colon", "SI")
has_groups <- main_groups[main_groups %in% unique(combined$tissue_group)]
cat(sprintf("  Available main groups: %s\n", paste(has_groups, collapse = ", ")))

# If no Colon/SI, check for other intestinal labels
if (!"Colon" %in% has_groups && !"SI" %in% has_groups) {
  cat("\n  WARNING: No intestinal eosinophils found in integrated object.\n")
  cat("  Available tissue_groups:\n")
  print(table(combined$tissue_group))
  cat("  Proceeding with available groups...\n")
}

main_cells <- combined$tissue_group %in% has_groups
cat(sprintf("  Keeping %d cells in %d groups\n", sum(main_cells), length(has_groups)))
sub <- subset(combined, cells = colnames(combined)[main_cells])
sub$tissue_group <- factor(sub$tissue_group, levels = has_groups)
Idents(sub) <- "tissue_group"

# =============================================================================
# 3. Colors
# =============================================================================
group_colors <- c(
  eos_colors[paste0("C", 0:5)],
  "BM"    = "#1D3557",
  "Blood" = "#457B9D",
  "Colon" = "#E6AB02",
  "SI"    = "#A6761D"
)
names(group_colors)[1:6] <- paste0("Adipose_C", 0:5)
group_colors <- group_colors[names(group_colors) %in% has_groups]

# =============================================================================
# 4. Panel A: UMAP by tissue source
# =============================================================================
cat("\nPanel A: UMAP by tissue source...\n")

combined$source_broad <- case_when(
  combined$dataset == "Qiu_2025_Adipose" ~ "Adipose (Qiu)",
  combined$organ == "Bone Marrow"         ~ "BM (Gurtner)",
  combined$organ == "Blood"               ~ "Blood (Gurtner)",
  combined$organ == "Colon"               ~ "Colon (Gurtner)",
  grepl("[Cc]olon", combined$organ)       ~ "Colon (Gurtner)",
  combined$organ == "Small Intestine"     ~ "SI (Gurtner)",
  grepl("[Ss]mall|SI", combined$organ)    ~ "SI (Gurtner)",
  TRUE                                    ~ paste0(combined$organ, " (Gurtner)")
)

source_colors <- c(
  "Adipose (Qiu)"   = "#E63946",
  "BM (Gurtner)"    = "#1D3557",
  "Blood (Gurtner)" = "#457B9D",
  "Colon (Gurtner)" = "#E6AB02",
  "SI (Gurtner)"    = "#A6761D"
)

p_a <- DimPlot(combined, group.by = "source_broad", pt.size = 0.3,
               cols = source_colors, shuffle = TRUE) +
  ggtitle("Integrated UMAP: Adipose vs BM/Blood/Intestinal") +
  theme_umap()

save_pdf(p_a, file.path(out_dir, "supp_gurtner_intestinal_A_umap.pdf"), 9, 7)

# =============================================================================
# 5. Panel B: Marker dotplot
# =============================================================================
cat("Panel B: Marker dotplot...\n")

marker_genes <- c(
  # Pan-eos
  "Siglecf", "Epx", "Prg2", "Ccr3", "Il5ra",
  # Maturation
  "Cd63", "Cd9", "Itgam", "Itgb2",
  # Cytokines / effectors
  "Il4", "Il13", "Ccl6", "Retnla",
  # MHC-II (C0)
  "H2-Ab1", "H2-Aa", "Cd74",
  # Migration (C1)
  "Sell", "Cxcr4",
  # Stress (C2)
  "Atf6", "Ddit3",
  # Tissue resident / intestinal
  "Cd44", "Cd22", "Clec4a4",
  # Progenitor
  "Kit", "Cd34",
  # Intestinal markers (from Hu et al.)
  "Alox15", "Retnlg", "Mmp9"
)
marker_genes <- marker_genes[marker_genes %in% rownames(sub)]
cat(sprintf("  %d markers found\n", length(marker_genes)))

p_b <- DotPlot(sub, features = marker_genes, scale = TRUE,
               col.min = -1.5, col.max = 2.5) +
  scale_color_gradientn(
    colors = c("#3B4CC0", "#7092D5", "#AAC7FD", "#DDDDDD",
               "#F7B799", "#E8744A", "#B40426"),
    name = "Scaled\nExpression") +
  scale_size_continuous(range = c(0.5, 6), name = "Percent\nExpressed") +
  labs(x = "", y = "",
       title = "Eosinophil Markers: Adipose vs BM/Blood/Intestinal") +
  theme_pub() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1, size = 8, face = "italic"),
        axis.text.y = element_text(size = 9))

save_pdf(p_b, file.path(out_dir, "supp_gurtner_intestinal_B_dotplot.pdf"),
         max(14, length(marker_genes) * 0.4), 7)

# =============================================================================
# 6. Panel C: Spearman correlation heatmap (Harmony embeddings)
# =============================================================================
cat("Panel C: Spearman correlation heatmap...\n")

# Use Harmony embeddings
harmony_emb <- Embeddings(sub, "harmony")
cat(sprintf("  Harmony: %d cells x %d dims\n", nrow(harmony_emb), ncol(harmony_emb)))

groups <- as.character(Idents(sub))
unique_groups <- levels(Idents(sub))
unique_groups <- unique_groups[unique_groups %in% unique(groups)]

avg_harmony <- do.call(rbind, lapply(unique_groups, function(g) {
  colMeans(harmony_emb[groups == g, , drop = FALSE])
}))
rownames(avg_harmony) <- unique_groups

# Spearman correlation
cor_mat <- cor(t(avg_harmony), method = "spearman")

# Annotation
ann_row <- data.frame(
  Source = ifelse(grepl("^Adipose", rownames(cor_mat)), "Qiu Adipose",
           ifelse(grepl("BM", rownames(cor_mat)), "Gurtner BM",
           ifelse(grepl("Blood", rownames(cor_mat)), "Gurtner Blood",
           ifelse(grepl("Colon", rownames(cor_mat)), "Gurtner Colon",
           ifelse(grepl("SI", rownames(cor_mat)), "Gurtner SI", "Other"))))),
  row.names = rownames(cor_mat)
)
ann_colors <- list(
  Source = c("Qiu Adipose" = "#E63946", "Gurtner BM" = "#1D3557",
             "Gurtner Blood" = "#457B9D", "Gurtner Colon" = "#E6AB02",
             "Gurtner SI" = "#A6761D")
)

# Display labels
display_labels <- gsub("Adipose_", "", rownames(cor_mat))

n_groups <- ncol(cor_mat)
pdf(file.path(out_dir, "supp_gurtner_intestinal_C_heatmap.pdf"),
    width = max(10, n_groups * 0.8), height = max(8, n_groups * 0.7))
pheatmap(cor_mat,
         color = colorRampPalette(c("#2166AC", "white", "#B2182B"))(100),
         breaks = seq(min(cor_mat, na.rm = TRUE), 1, length.out = 101),
         cluster_rows = TRUE, cluster_cols = TRUE,
         annotation_row = ann_row, annotation_col = ann_row,
         annotation_colors = ann_colors,
         display_numbers = TRUE, number_format = "%.2f",
         fontsize_number = max(5, 10 - n_groups * 0.3),
         fontsize_row = 10, fontsize_col = 10,
         labels_row = display_labels, labels_col = display_labels,
         main = "Cross-Tissue Eosinophil Similarity (Spearman)",
         border_color = NA)
dev.off()
cat("  Saved: supp_gurtner_intestinal_C_heatmap.pdf\n")

# Save correlation CSV
write.csv(cor_mat, file.path(out_dir, "gurtner_intestinal_spearman.csv"))

# =============================================================================
# 7. Panel D: Dendrogram
# =============================================================================
cat("Panel D: Dendrogram...\n")

dist_mat <- as.dist(1 - cor_mat)
hc <- hclust(dist_mat, method = "average")

dend <- as.dendrogram(hc)
orig_labels <- labels(dend)
labels(dend) <- gsub("Adipose_", "", orig_labels)

# Color labels by source
get_label_color <- function(x) {
  ifelse(grepl("^Adipose_C[0-5]$", x), "#E63946",
  ifelse(grepl("^BM$", x), "#1D3557",
  ifelse(grepl("^Blood$", x), "#457B9D",
  ifelse(grepl("^Colon$", x), "#E6AB02",
  ifelse(grepl("^SI$", x), "#A6761D", "black")))))
}
labels_colors(dend) <- sapply(orig_labels, get_label_color)
labels_cex(dend) <- 1.1

pdf(file.path(out_dir, "supp_gurtner_intestinal_D_dendrogram.pdf"), width = 12, height = 7)
par(mar = c(8, 4, 4, 2))
plot(dend, main = "Eosinophil Cross-Tissue Similarity",
     ylab = "1 - Spearman correlation", las = 2)
legend("topright",
       legend = c("Qiu Adipose (C0-C5)", "Gurtner BM",
                  "Gurtner Blood", "Gurtner Colon", "Gurtner SI"),
       text.col = c("#E63946", "#1D3557", "#457B9D", "#E6AB02", "#A6761D"),
       bty = "n", cex = 1.0)
dev.off()
cat("  Saved: supp_gurtner_intestinal_D_dendrogram.pdf\n")

# =============================================================================
# 8. Summary
# =============================================================================
cat("\n========== Summary ==========\n")
cat(sprintf("Total cells: %d\n", ncol(sub)))
cat(sprintf("Groups: %s\n", paste(has_groups, collapse = ", ")))

# Print closest tissue for each adipose cluster
cat("\nClosest tissue for each adipose cluster:\n")
for (cl in paste0("Adipose_C", 0:5)) {
  if (!cl %in% rownames(cor_mat)) next
  non_adipose <- setdiff(rownames(cor_mat), grep("^Adipose", rownames(cor_mat), value = TRUE))
  if (length(non_adipose) == 0) next
  cors <- cor_mat[cl, non_adipose]
  best <- names(which.max(cors))
  cat(sprintf("  %s -> %s (r = %.3f)\n", cl, best, cors[best]))
}

cat(sprintf("\nOutput: %s\n", out_dir))
cat("39 complete.\n")
