# Methods — Adipose Eosinophil Subtype Atlas (DecontX Revision)

## Ambient RNA Removal (DecontX)
Ambient RNA contamination was removed per sample using DecontX (celda R package). For each of the 16 samples, raw (unfiltered) CellRanger counts were used as the background profile and filtered counts as the cell matrix. Cells with contamination scores exceeding 0.3 were excluded. Decontaminated counts replaced the original counts matrix for all downstream analyses.

## Quality Control and Doublet Removal
Doublets were identified using Scrublet (Python, via reticulate). QC filtering removed cells with contamination scores > 0.2. The 16 samples were merged into a single Seurat (v5) object for downstream processing.

## Normalization, Dimensionality Reduction, and Clustering
Gene expression was normalized using LogNormalize (scale factor 10,000). The top 2,000 highly variable genes were identified via variance-stabilizing transformation (vst). Principal component analysis (PCA) was performed on scaled data (30 PCs). A shared nearest neighbor (SNN) graph was constructed (dims = 1:20, k = 30), and UMAP was computed (n.neighbors = 50, min.dist = 0.1, spread = 1). Cells were clustered at multiple resolutions (0.5, 0.8, 1.0) using the Louvain algorithm; resolution 0.8 was used as the default.

## Cell Type Annotation
Cell types were assigned by label transfer from a previously annotated reference object (13 cell types, 97,061 cells) using Seurat's anchoring framework. Seven immune cell types were identified: eosinophils, macrophages, dendritic cells, B cells, NK cells, neutrophils, and mast cells.

## Eosinophil Subclustering and Resolution Justification
Eosinophils were subset and re-processed (normalization, 2,000 HVGs, PCA with 30 PCs, UMAP with dims 1:20). Clustering was performed across 12 resolutions (0.1-1.5). Resolution 0.4 was selected (7 clusters, C0-C6) based on four convergent lines of evidence:
1. **Clustree analysis**: 7 clusters were stable across resolutions 0.4-0.6
2. **Silhouette analysis**: local peak in mean silhouette score at resolution 0.4
3. **Subsampling stability**: adjusted Rand index (ARI) computed over 100 iterations of 80% subsampling showed high reproducibility
4. **sc-SHC significance testing** (Grabski et al., *Nat Methods* 2023): all clusters passed post-hoc significance testing

Cluster C6 was identified as too few cells and excluded from all downstream analyses, yielding 6 final eosinophil subtypes (C0-C5). Cluster markers were identified using FindAllMarkers (min.pct = 0.25, logfc.threshold = 0.25, Wilcoxon rank-sum test).

## Macrophage Subtype Classification
Macrophages were classified into M1, M2 surface, and M2 functional subtypes using Seurat's AddModuleScore (ctrl = 50 control genes). Module gene sets: M1 (*Itgax, Cd80, Tnf, Il6, Ccl2*), M2 surface (*Mrc1, Lyve1, Cd301b, Cd163*), M2 functional (*Chil3, Il10, Arg1, Tgfb1, Clec10a*). Each cell was assigned to the subtype with the highest module score.

## Differential Expression Analysis
### Single-cell level (MAST)
DEGs were identified using MAST with a hurdle model including a random intercept for sample: `zlm(~ group + cngeneson + (1|sample), method = "glmer")`. Genes expressed in <10% of cells were excluded. Significance thresholds: adjusted p-value < 0.05 (Benjamini-Hochberg), |log2FC| > 0.25. Four comparisons were performed: sex effect (WT Female vs Male HFD), genotype effect in females (KO vs WT Female HFD), genotype effect in males (KO vs WT Male HFD), and infection effect (HFD Infected vs Naive).

### Pseudobulk (DESeq2)
Raw counts were aggregated per sample per condition using pseudobulk summation. DESeq2 was applied with the same four comparisons. Concordance between MAST and DESeq2 results was assessed.

## Pathway Enrichment Analysis
Gene Ontology (GO) Biological Process and KEGG pathway enrichment was performed on cluster-specific markers using clusterProfiler (org.Mm.eg.db). Markers were filtered at adjusted p-value < 0.05 and avg_log2FC > 0.25, with a minimum of 5 genes per cluster.

## Cell-Cell Communication (CellChat)
Eosinophil-macrophage communication was analyzed using CellChat (v1.6.1) with the CellChatDB.mouse database. Communication probabilities were computed using the triMean method with a minimum of 10 cells per interaction. Analysis was performed separately for each of 6 experimental conditions (WT/KO x Male/Female x HFD/Infected). Per-condition CellChat objects were merged for comparative analysis of signaling pathway changes.

## Harmony Integration with Reference Datasets
### Two-way integrations
Qiu HFD adipose eosinophils (C0-C5) were independently integrated with (1) Hu et al. 2026 healthy adipose eosinophils and (2) Gurtner et al. 2022 blood/bone marrow eosinophils (GSE182001) using Harmony batch correction. For each integration: 2,000 variable features were selected, PCA was computed (30 PCs), Harmony was run with `group.by.vars = "dataset"`, and UMAP was generated on Harmony embeddings (dims 1:30). Integrated clustering used resolution 0.5.

### Three-way integration
All three datasets were merged and integrated using Harmony (group.by = dataset). The Hu dataset retained depot-level metadata (BAT, eWAT, scWAT). Gurtner clusters were mapped to their published nomenclature (Gurtner et al., *Nature* 2023) based on tissue distribution and validation against six canonical markers from their Figure 1d:

| Cluster | Tissue | Cd80 | Cd274 (PD-L1) | Retnla | Mki67 | Tgfb1 | Label |
|---------|--------|------|---------------|--------|-------|-------|-------|
| 4 | BM (94%) | 1.3% | 7.1% | 35.4% | **84.1%** | 18.0% | Eos Precursors |
| 16 | BM (91%) | 6.4% | 5.4% | 25.0% | **89.3%** | 25.7% | Cycling Precursors |
| 3 | BM (69%) | 2.8% | 7.5% | 17.7% | 8.0% | 5.0% | Immature Eos |
| 2 | Blood (64%) | 6.4% | 4.6% | **35.8%** | 2.2% | 1.3% | Circulating Eos |
| 1 | GI (64%) | **49.1%** | **47.2%** | 22.4% | 0.6% | 3.0% | A-Eos (active) |
| 5 | GI (68%) | 24.3% | 33.7% | 21.6% | 11.2% | **19.5%** | B-Eos (basal) |

Cluster 4 was identified as eosinophil precursors based on dominant Mki67 expression (84.1%), consistent with cycling progenitors. Cluster 2 was assigned as circulating eosinophils based on the highest Retnla expression (35.8%) and blood-dominant tissue distribution. Cluster 1 was confirmed as A-Eos (active) based on the highest Cd80 (49.1%) and Cd274/PD-L1 (47.2%) expression. Cluster 5 was assigned as B-Eos (basal) based on the highest Tgfb1 expression (19.5%) and GI-dominant tissue localization, though it showed intermediate Cd80/Cd274 levels suggesting a possible transitional phenotype. Cluster 3 was identified as immature eosinophils based on low expression of all subset-specific markers and BM-dominant localization. Cluster 16 (cycling precursors, Mki67 = 89.3%) was excluded from the main dendrogram analysis due to low cell count (n = 280). Data were exported for downstream trajectory analysis in Python (counts matrix, Harmony embeddings, UMAP coordinates, metadata).

## Trajectory Analysis
### CytoTRACE / CellRank
CytoTRACE scores were computed using the CytoTRACEKernel in CellRank (v2) on the three-way integrated data. Log-normalized expression (layer "X") was used directly. The transition matrix was computed with soft thresholding. Diffusion pseudotime (DPT) was calculated with 5 diffusion components, rooted at the Gurtner bone marrow cell with the highest CytoTRACE score. GPCCA macrostate analysis was performed with 20 Schur components and 5 macrostates.

### Monocle3
The three-way integrated Seurat object was converted to a Monocle3 cell_data_set. UMAP coordinates from the Harmony integration were transferred. The trajectory graph was learned using `learn_graph()` after preprocessing (30 dimensions). The root was set to the Gurtner bone marrow cell with the highest CytoTRACE score. Two figure sets were generated: (1) main figure with BM + Blood + Qiu C0-C5 only, (2) supplementary with Hu adipose depots included.

## Cross-Dataset Cluster Similarity
Spearman correlation was computed on average Harmony embeddings (30 dimensions) per group to assess cross-dataset similarity while controlling for batch effects. Groups included: Qiu C0-C5, Hu eWAT/scWAT/BAT, and Gurtner Precursors/Immature/Circulating/A-Eos/B-Eos. Hierarchical clustering (average linkage) on the distance matrix (1 - Spearman correlation) was used to generate a dendrogram. Harmony embeddings were used instead of raw expression to remove technical batch effects between datasets, enabling biologically meaningful cross-dataset comparisons.

## SCENIC Gene Regulatory Network Analysis
Single-cell regulatory network inference and clustering (pySCENIC) was used to identify transcription factor regulons. GRN inference was performed using GRNBoost2 (arboreto) on the DecontX-cleaned eosinophil expression matrix. Regulon pruning was performed using cisTarget with two ranking databases (mm10 refseq-r80, 10kb and 500bp windows) and motif annotations (motifs-v9-nr.mgi). AUCell was used to score regulon activity per cell. Regulon specificity scores (RSS) were computed per cluster using Jensen-Shannon divergence. Macrophage cells were scored with eosinophil-derived regulons for comparative analysis (Fig 5A). Condition-specific regulon activity was assessed using Kruskal-Wallis tests and pairwise comparisons across the 6 experimental conditions.

## Software Versions
- R 4.5.0, Seurat v5.4.0, Harmony v1.2.1, MAST, DESeq2, CellChat v1.6.1, monocle3 v1.4.26, clusterProfiler, celda (DecontX)
- Python 3.12, scanpy v1.10.3, cellrank v2, pySCENIC, Scrublet
- All analyses were performed on the UCR HPCC cluster using SLURM job scheduling
