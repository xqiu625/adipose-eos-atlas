# Methods — Adipose Eosinophil Cluster Atlas

This file mirrors the methods text in the manuscript. Numbers in brackets refer to scripts under `scripts/`.

## Ambient RNA Removal (DecontX)
Ambient RNA contamination was removed per sample using DecontX (celda R package). For each of 16 samples, raw (unfiltered) Cell Ranger counts served as the background profile and filtered counts as the cell matrix. Cells with contamination scores > 0.2 were excluded; decontaminated counts replaced raw counts for all downstream analyses. [`00_decontx_per_sample.R`, `21_ambient_rna_evidence.R`]

## Quality Control and Doublet Removal
Cells were retained with 200–5,000 detected features and ≤10% mitochondrial content. Doublets were identified using Scrublet (Python, via reticulate) with an expected doublet rate of 6%. The 16 samples were merged into a single Seurat v5 object (97,934 → 7,461 high-quality eosinophils across 16 samples). [`01_qc_doublet_filter.R`, `17_qc_supplement.R`]

## Normalization, Dimensionality Reduction, and Clustering
Gene expression was LogNormalized (scale factor 10,000); 2,000 highly variable features were selected (vst). PCA was performed on scaled data (30 PCs). SNN graphs were built on dims 1:20, and UMAP was computed (n.neighbors = 50, min.dist = 0.1). Cells were clustered using the Louvain algorithm.

## Cell Type Annotation
Cell types were assigned by marker-based annotation and label transfer using Seurat's anchoring framework. Seven immune cell types were identified (eosinophils, macrophages, dendritic cells, B cells, NK cells, neutrophils, mast cells). [`02_all_cells_umap.R`, `03_annotate_celltypes.R`]

## Eosinophil Subclustering — Re-processing and Resolution Justification
Eosinophils were subset and re-processed *de novo* (NormalizeData, FindVariableFeatures n = 2,000, ScaleData, RunPCA npcs = 30, RunUMAP dims = 1:20, FindNeighbors, FindClusters) to capture eosinophil-intrinsic variation. Resolution selection was assessed across 12 resolutions (0.1–1.5). **Resolution 0.4 was selected** based on four convergent lines of evidence:
1. **Clustree** stability — 7 clusters stable across resolutions 0.3–0.6.
2. **Silhouette** — mean silhouette plateaued at resolutions 0.3–0.5.
3. **Subsampling stability** — ARI over 100 iterations of 80% subsampling.
4. **sc-SHC** significance testing (Grabski et al., *Nat Methods* 2023).

One small cluster (n = 14 cells) was excluded due to its low cell number, yielding **six final eosinophil clusters (C0–C5)**. Cluster markers were identified using FindAllMarkers (min.pct = 0.25, logfc.threshold = 0.25, Wilcoxon). [`04_eos_subcluster.R`, `04b_rename_clusters.R`, `18_clustree_resolution.R`]

## Batch Effect Assessment and Harmony Sensitivity
Samples were multiplexed using the 10x Genomics Flex probe-based assay (4 samples per GEM reaction), which eliminates the primary technical batch sources at the experimental level (PMIDs 38612639, 35927335). Within-study Harmony integration was tested as a sensitivity analysis (Supplementary Fig. S7) but disrupted cluster structure (ARI = 0.28, mean Jaccard index = 0.54) because GEM reaction is partially confounded with biological condition; integration cannot distinguish these and removes biological signal along with technical noise (PMIDs 34949812, 31948481). Three independent analyses confirmed clusters are not batch-driven: (i) every cluster contains cells from all six experimental conditions, (ii) replicates of the same condition processed in different GEM reactions show concordant cluster proportions, (iii) per-cluster condition composition tracks biological variables rather than technical covariates. [`19_batch_effect_evidence.R`]

## Compositional Analysis (sccomp)
Cluster-proportion shifts across conditions were tested using Bayesian compositional analysis with `sccomp` v1.7 (Beta-Binomial framework, accounting for the sum-to-one constraint and sample-level pseudoreplication; PMID 37549298). Three contrasts: sex effect (WT Female vs WT Male HFD, n = 3 vs 3), genotype effect (KO vs WT females HFD, n = 3 vs 3), and infection effect (Infected vs Naive females, n = 4 vs 6). FDR-corrected posterior probabilities < 0.1 were considered significant. [`20_compositional_analysis.R`, `43_composition_biorep.R`]

## Macrophage Subtype Classification
Macrophages within the SiglecF⁺ fraction were classified using Seurat AddModuleScore (ctrl = 50) with curated gene sets: M1-like (*Itgax, Cd80, Tnf, Il6, Ccl2*), M2 Surface (*Mrc1, Lyve1, Cd301b, Cd163*), M2 Functional (*Chil3, Il10, Arg1, Tgfb1, Clec10a*). Each cell was assigned to the subtype with the highest module score. We refer to these as M1-like / M2-like states defined by module score. [`06_mac_subtype.R`]

## Differential Expression — Pseudobulk DESeq2
Condition-level differential expression was performed using a pseudobulk approach (Squair et al., *Nat Commun* 2021; PMID 34584091). Single-cell counts were aggregated per sample using `AggregateExpression()`, then tested with DESeq2. Multifactor designs (~genotype + infection) were used for the genotype-in-females and infection-in-females comparisons to control for confounding. Significance: adjusted p < 0.05 and |log₂FC| > 0.25.

Four comparisons:
- **Panel D** WT Female vs WT Male (n = 3 vs 3)
- **Panel E** KO Male vs WT Male (n = 3 vs 3)
- **Panel F** KO vs WT, all females, multifactor (n = 5 vs 5)
- **Panel G** Infected vs Naive, all females, multifactor (n = 4 vs 4)

Single-cell-level `FindMarkers` (Wilcoxon) was used only for identifying cluster-defining markers, not for condition-level inference. [`05_eos_pseudobulk.R`, `05b_pseudobulk_methods.R`, `07c_deg_pseudobulk.R`, `25_fig1_revision_v3.R`, `27_deg_excel_export.R`]

## Pathway Enrichment Analysis
GSEA was performed on DESeq2-ranked gene lists using `clusterProfiler::gseGO()` (Gene Ontology Biological Process; mouse `org.Mm.eg.db`). Genes were ranked by DESeq2 test statistic; enrichment was assessed with Benjamini–Hochberg correction. [`25_fig1_revision_v3.R`; cluster-marker over-representation in `08_pathway_eos_clusters.R`]

## Cell–Cell Communication (CellChat)
Eosinophil–macrophage communication was inferred using CellChat v1.6.1 with the mouse CellChatDB. Communication probabilities were computed using the triMean method with a minimum of 10 cells per interaction. CellChat objects were built separately per condition (WT/KO × Male/Female × HFD/Infected) and merged for comparative analysis. Five focal pathways were highlighted: CCL, FN1, APP, CD45, THBS. Network centrality (`netAnalysis_signalingRole_network()`, PMID 33597522) identified senders, receivers, mediators (betweenness centrality), and influencers (information centrality). [`10_cellchat.R`, `10b_cellchat_manuscript_figures.R`, `10c_cellchat_fig3d.R`, `33_fig3EF_cellchat_ccl_heatmap.R`, `44_fig3A_chord_revised.R`]

CellChat interaction probabilities represent computational predictions based on ligand–receptor co-expression, not experimentally validated signaling events.

## Spatial Analysis (RNAscope)
Crown-like structures from RNAscope confocal z-stacks (CCL6, CCR3, CCR2) were reconstructed in 3D using Imaris 11.0. Puncta and nuclei coordinates were exported for quantification. A modified G-cross function (PMID 29409671) was applied, taking the mean AUC at r' = 25 μm for each reference nucleus and normalizing via permutation. Two-way ANOVA with Tukey's HSD was applied (n = 4 biological replicates per group, two independent cohorts).

## Trajectory Analysis (Monocle 3)
The eosinophil dataset was integrated with reference bone-marrow and blood eosinophils from Gurtner et al. (PMID 36509106) using Harmony, and trajectory inference was performed using Monocle 3 v1.2.9. Trajectory alignment used `align_cds()`. The root was anchored to the Gurtner bone marrow reference population.

For the pseudotime gene-expression heatmap, pseudotime-dependent genes were identified via Monocle 3's `graph_test()` (Moran's I, q-value < 0.05) on the top 1,000 highly variable genes across pseudotime-ordered cells. Genes were grouped into six co-expression modules by hierarchical clustering of pseudotime dynamics; the top three by variance within each module were labeled as representative markers per maturation stage.

While pseudotime orders cells along a continuum, it does not constitute definitive evidence of a differentiation trajectory; this would require lineage tracing or temporal fate-mapping. [`11_monocle3.R`, `11b_monocle3_integrated.R`, `24_integrate_qiu_hu_gurtner.R`, `24b_cellrank_eos_maturation.py`, `24c_pseudotime_violin.py`]

## Cross-Dataset Integration and Similarity
Two cross-study integrations were performed using Harmony v1.2.1 (default parameters, top 30 PCs):
1. Our HFD adipose eosinophils with Hu et al. 2026 healthy adipose eosinophils (eWAT, scWAT, BAT) to assess cluster conservation across dietary states.
2. A three-study integration with Hu et al. 2026 (healthy adipose) and Gurtner et al. 2023 (GSE182001; bone marrow and blood) to map eosinophil maturation from bone marrow precursors through circulation to tissue-resident states.

Cross-dataset transcriptional similarity was assessed using Spearman correlation on Harmony embeddings (30 dims) per group; hierarchical clustering (average linkage on 1 − Spearman distance) generated the dendrogram. [`12_integrate_hu.R`, `13_integrate_gurtner.R`, `23_hu_mast_marker_validation.R`]

## Gurtner Cluster Labeling
Gurtner et al. clusters were mapped to their published nomenclature based on tissue distribution and canonical markers from their Fig 1d:

| Cluster | Tissue | Cd80 | Cd274 | Retnla | Mki67 | Tgfb1 | Label |
|---------|--------|------|-------|--------|-------|-------|-------|
| 4 | BM (94%) | 1.3% | 7.1% | 35.4% | **84.1%** | 18.0% | Eos Precursors |
| 16 | BM (91%) | 6.4% | 5.4% | 25.0% | **89.3%** | 25.7% | Cycling Precursors |
| 3 | BM (69%) | 2.8% | 7.5% | 17.7% | 8.0% | 5.0% | Immature Eos |
| 2 | Blood (64%) | 6.4% | 4.6% | **35.8%** | 2.2% | 1.3% | Circulating Eos |
| 1 | GI (64%) | **49.1%** | **47.2%** | 22.4% | 0.6% | 3.0% | A-Eos (active) |
| 5 | GI (68%) | 24.3% | 33.7% | 21.6% | 11.2% | **19.5%** | B-Eos (basal) |

Cluster 16 was excluded from the dendrogram due to low cell count (n = 280).

## SCENIC Gene Regulatory Network Analysis
pySCENIC was applied to infer transcription factor regulons and quantify regulatory activity at single-cell resolution. Co-expression networks were inferred with GRNBoost2 *independently per eosinophil cluster* (C0–C5) to amplify cluster-specific regulatory signal. Modules were refined using cisTarget (NES ≥ 3.0; mm10 motif database v9, 10kb and 500bp windows; motifs-v9-nr.mgi). Regulon activity was scored using AUCell (top 5% of expressed genes per cell).

Differential regulon activity between conditions was assessed using a random-subsampling strategy (200 cells per group; adapted from Xie et al. 2020, PMID 32719519). Regulons were flagged as "eosinophil-enriched" if maximum AUC across C0–C5 ≥ 0.05 and fold-change vs mean macrophage AUC ≥ 1.2. Sex-differential regulons were identified by comparing WT Female vs WT Male HFD-only samples (Wilcoxon p < 0.05); only TFs whose SCENIC-identified target gene expression confirmed the same sex-biased direction were retained. [`14_scenic_export.R`, `15_scenic_run.sh`, `15c_scenic_per_cluster.sh`, `16_scenic_visualize.py`, `16b–16e_*.py`, `22_scenic_by_condition.py`, `38_scenic_auc_violin.py`, `41–42*_hu_regulon_*`]

## Software Versions
- **R** 4.5.0 — Seurat 5.1.0, DESeq2, Harmony 1.2.1, CellChat 1.6.1, Monocle 3 1.2.9, clusterProfiler, celda (DecontX), sccomp 1.7
- **Python** 3.12 — scanpy 1.10.3, pyscenic, cellrank 2.0.7, scrublet
- Analyses run on the UCR HPCC cluster (SLURM); see `scripts/SLURM_COMMANDS.sh` for submission patterns.
