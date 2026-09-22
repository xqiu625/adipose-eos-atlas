#!/usr/bin/env Rscript
lib_path <- Sys.getenv("RERUN_R_LIB", "/bigdata/godziklab/shared/Xinru/R")
.libPaths(c(lib_path, .libPaths()))

# ==============================================================================
# 60_eos_vs_mac_dendrogram.R   (Reviewer #1 — Meera's requested panel)
# Similarity gene-tree / hierarchical-clustering dendrogram of the six eosinophil
# clusters (C0-C5) vs the macrophage subtypes (Mac_M1, Mac_M2_Surf, Mac_M2_Func).
# If C0-C5 form one branch cleanly separated from the macrophage subtypes, that is
# a direct, intuitive demonstration that none of the eosinophil clusters are
# macrophages. Pseudobulk per group on shared HVGs -> Spearman correlation ->
# hierarchical clustering (Ward on 1-correlation).
# Input : eos_subclustered.rds (C0-C5), mac_subtyped.rds (mac_subtype)
# Output: $BASE/revision/figures/eos_vs_mac/
# ==============================================================================
suppressPackageStartupMessages({ library(Seurat); library(Matrix) })

BASE <- "/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas"
out  <- file.path(BASE, "revision", "figures", "eos_vs_mac"); dir.create(out, recursive=TRUE, showWarnings=FALSE)

eos <- readRDS(file.path(BASE, "eos_subclustered.rds"))
if ("C6" %in% levels(Idents(eos))) { eos <- subset(eos, idents="C6", invert=TRUE); Idents(eos)<-droplevels(Idents(eos)) }
if (inherits(eos[["RNA"]],"Assay5")) eos[["RNA"]] <- JoinLayers(eos[["RNA"]])
eos$grp <- paste0("Eos_", as.character(Idents(eos)))

mac <- readRDS(file.path(BASE, "mac_subtyped.rds"))
if (inherits(mac[["RNA"]],"Assay5")) mac[["RNA"]] <- JoinLayers(mac[["RNA"]])
mac$grp <- as.character(mac$mac_subtype)   # Mac_M1, Mac_M2_Func, Mac_M2_Surf

common <- intersect(rownames(eos), rownames(mac))
cat("common genes:", length(common), "\n")

# shared HVGs across the combined set (so the tree is driven by informative genes)
ce <- GetAssayData(eos, layer="counts")[common, ]
cm <- GetAssayData(mac, layer="counts")[common, ]
comb <- CreateSeuratObject(cbind(ce, cm))
comb$grp <- c(eos$grp, mac$grp)
comb <- NormalizeData(comb, verbose=FALSE)
comb <- FindVariableFeatures(comb, nfeatures=2000, verbose=FALSE)
hvg <- VariableFeatures(comb)
data <- GetAssayData(comb, layer="data")[hvg, ]

groups <- sort(unique(comb$grp))
pb <- sapply(groups, function(g) rowMeans(data[, comb$grp == g, drop=FALSE]))
cormat <- cor(pb, method="spearman")
write.csv(round(cormat,3), file.path(out, "eos_mac_pseudobulk_correlation.csv"))
cat("\n=== pseudobulk Spearman correlation (groups) ===\n"); print(round(cormat,2))

# hierarchical clustering on 1 - correlation
d <- as.dist(1 - cormat)
hc <- hclust(d, method="ward.D2")
saveRDS(list(cormat=cormat, hc=hc, pb=pb), file.path(out, "eos_vs_mac_tree.rds"))

# is the eos/mac split the top split?
ct <- cutree(hc, k=2)
cat("\n=== top-level 2-cut membership (expect Eos in one cluster, Mac in the other) ===\n")
print(data.frame(group=names(ct), branch=as.integer(ct)))

# export newick-ish ordering + a base-R dendrogram PDF (figure styled later in Python)
pdf(file.path(out, "eos_vs_mac_dendrogram.pdf"), width=7, height=5)
plot(hc, main="Eosinophil clusters vs macrophage subtypes\n(pseudobulk, 1 - Spearman, Ward.D2)",
     xlab="", sub="", cex=0.9)
dev.off()
# also write linkage as a table for Python plotting
merge_df <- as.data.frame(hc$merge); merge_df$height <- hc$height
write.csv(merge_df, file.path(out, "eos_vs_mac_linkage.csv"), row.names=FALSE)
write.csv(data.frame(order=hc$order, label=hc$labels[hc$order]), file.path(out, "eos_vs_mac_leaforder.csv"), row.names=FALSE)
cat("\nDone 60_eos_vs_mac_dendrogram.R\n")
