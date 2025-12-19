# Adipose Eosinophil and Macrophage Single-Cell Atlas

Single-cell RNA sequencing analysis of eosinophils and macrophages in adipose tissue under metabolic stress conditions, with focus on RELMα function in obesity.

**Publication:** Ruggiero-Ruff and Qiu et al. 2025. "An intra-adipose network of eosinophil subtypes is regulated by sex and RELMα and drives macrophage polarization in obesity-induced inflammation"

## Overview

This project investigates immune cell heterogeneity and cell-cell communication in adipose tissue using scRNA-seq. We characterize eosinophil and macrophage subpopulations across different metabolic conditions and genotypes, identifying seven eosinophil subtypes (C0-C6) that dynamically shift based on sex, RELMα expression, and helminth infection.

### Study Design

| Factor | Conditions |
|--------|------------|
| Genotype | Wild-Type (WT), RELMα knockout (KO) |
| Sex | Male, Female |
| Diet | High-Fat Diet (HFD) |
| Challenge | HFD + Helminth Infection |

## Repository Structure

```
├── scripts/
│   ├── 01_preprocessing/           # Data loading, QC, Seurat processing
│   ├── 02_cell_annotation/         # Cell type identification (UCE)
│   ├── 03_eosinophil_analysis/     # Subclustering, DEG, M1/M2 scoring
│   ├── 04_trajectory_analysis/     # scVelo, CellRank, pseudotime
│   ├── 05_cell_communication/      # CellChat, CellPhoneDB, CCL analysis
│   ├── 06_scenic_grn/              # SCENIC, gene regulatory networks
│   ├── 07_integration/             # Integration with Gurtner et al.
│   └── 08_figure_generation/       # Publication figure scripts
└── README.md
```

## Scripts Description

### 01_preprocessing/
- `01_local_data_seurat_preprocess.R` - Initial data loading and Seurat object creation
- `03_UMAP.R` - Dimensionality reduction and visualization
- `04_FindAllMarkers_cluster.R` - Cluster marker identification

### 02_cell_annotation/
- `05_UCE_annotate_by_mouse_kidney_norn.py` - Universal Cell Embedding annotation
- `06_label_UCE_identity.R` - Cell type labeling based on UCE results

### 03_eosinophil_analysis/
- `03_EOS_DEG.R` - Differential expression analysis between conditions
- `04_EOS_ClusterProfiler.R` - GO/KEGG pathway enrichment
- `06_EOS_subcluster_0.5.R` / `07_EOS_subcluster_0.8.R` - Eosinophil subclustering
- `08_GSVA.R` - Gene Set Variation Analysis
- `01_define_M1_M2_signatures.R` - M1/M2 macrophage signature definitions
- `02_score_cells_M1_M2.R` - Macrophage polarization scoring
- `01_Macrophage_Cluster_DEG.R` - Macrophage differential expression
- `01_pseudobulk_by_samplev3.R` - Pseudobulk analysis

### 04_trajectory_analysis/
- `02_EOS_scVelo_preparation.py` - Prepare data for RNA velocity
- `03_EOS_scVelo_dynamics.py` - RNA velocity analysis
- `01_EOS_scvelo_bold_arrow_1024.py` - Velocity stream visualization
- `02_CellRank2.py` - CellRank trajectory inference
- `04_peudotime_heatmap_1024.py` - Pseudotime heatmap generation

### 05_cell_communication/
- `02_cellchat_score_eos-mac_condition_1021.R` - CellChat analysis by condition
- `03_cellchat_score_eos-mac_eos-sub_1021.R` - CellChat with eosinophil subtypes
- `04.1_CCL_pathway_heatmaps.R` - CCL signaling pathway visualization
- `04.4_Ccl_selected_ligands_clean_v2_1113.R` - CCL ligand-receptor analysis
- `04.7.1-outgoing_and_incoming_signaling_condition_improved_v2_1113.R` - Signaling direction analysis
- `1128_CellChat_MultiPanel_Difference_Heatmap.R` - Multi-panel CellChat visualization
- `1211_CCL_Dotplots_TOP20_v2.R` - CCL interaction dotplots
- CellPhoneDB scripts (`03-11`) - Alternative cell-cell interaction analysis

### 06_scenic_grn/
- `1201_Subcluster_Regulon_Analysis.py` - Regulon activity analysis
- `1202_Eos_vs_Mac_Subcluster_Unique_Zscore.py` - Eosinophil-specific regulons
- `1203_Eos_GRN_Network.py` - Eosinophil gene regulatory network
- `1203_Mac_GRN_Network.py` - Macrophage gene regulatory network
- `SCENIC_Condition_Comparison_Zscore.py` - Condition-specific regulon comparison
- `SCENIC_Publication_Figures.py` - Publication-ready SCENIC figures
- `SCENIC_RELMa_Firework_GRN.py` - RELMα-dependent GRN visualization

### 07_integration/
- `01_process_Gurtner.R` - Process Gurtner et al. reference data
- `04_harmony_integration.R` - Harmony batch correction
- `05_combined_UMAP.R` - Integrated UMAP visualization

### 08_figure_generation/
- `Fig8C_UMAP_plot_by_clusters_scVelo_dynamics.py` - Figure 2G: RNA velocity
- `Fig8D_Eos_condition_cluster_composition.R` - Figure 2D: Cluster composition
- `Fig8E_top_5markers_Eos.R` - Figure 2B: Eosinophil markers
- `DEG_Summary_Statistics.R` - DEG statistics for manuscript
- `GO_Analysis_All_Terms_Macrophages.R` - GO enrichment figures

## Key Tools & Dependencies

### R
- Seurat (v5.1.0)
- clusterProfiler
- CellChat (v1.6.1)
- harmony (v1.2.1)
- monocle3 (v1.2.9)
- GSVA

### Python
- scanpy (v1.10.3)
- scvelo (v0.2.5)
- cellrank (v2.0.7)
- pyscenic
- cellphonedb

## Cluster Definitions

### Eosinophil Subclusters
| Cluster | Color | Key Features |
|---------|-------|--------------|
| C0 | #6CA6CD | MHC-II presenting, male-enriched |
| C1 | #E84A5F | Migratory (Il4, Itgb2), female-enriched |
| C2 | #FF6F1F | Stress response (Atf6), male-enriched |
| C3 | #CBAACB | Tissue-resident, CCL6 signaling hub |
| C4 | #B03A48 | Infection response, female-enriched |
| C5 | #5BA4A4 | Metabolically adapted |
| C6 | #2F3C7E | Self-renewing |

### Macrophage Subtypes
- M1 (Pro-inflammatory): Cd80, Il6, Ccl2
- M2 Surface: Mrc1, Lyve1, Cd301b
- M2 Functional: Arg1, Il10, Chil3

## Code Availability

All custom scripts and analysis code used in this study are publicly available at https://github.com/xqiu625/adipose-eos-subtype-atlas.

## Data Availability

- **Raw Data:** GEO accession GSE311941
- **Reference Data:** Gurtner et al. (GSE182001) - blood and bone marrow eosinophils

## Citation

Ruggiero-Ruff RE*, Qiu X*, He Y, Vance L, Coss D, Godzik A+, Nair MG+. An intra-adipose network of eosinophil subtypes is regulated by sex and RELMα and drives macrophage polarization in obesity-induced inflammation. 2025.

*Equal contribution; +Co-corresponding authors

## Contact

- Meera G Nair: meera.nair@medsch.ucr.edu
- Adam Godzik: adam.godzik@medsch.ucr.edu

Division of Biomedical Sciences, School of Medicine, University of California, Riverside
