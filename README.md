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

This repository uses **Snakemake** for reproducible workflow management.

```
├── config/
│   └── config.yaml                 # Pipeline configuration
├── workflow/
│   ├── Snakefile                   # Main workflow definition
│   ├── rules/
│   │   ├── preprocessing.smk       # Data loading, QC, Seurat
│   │   ├── cell_annotation.smk     # UCE cell type identification
│   │   ├── eosinophil_analysis.smk # Subclustering, DEG, M1/M2
│   │   ├── trajectory_analysis.smk # scVelo, CellRank, pseudotime
│   │   ├── cell_communication.smk  # CellChat, CellPhoneDB
│   │   ├── scenic_grn.smk          # SCENIC, GRN analysis
│   │   ├── integration.smk         # Harmony integration
│   │   └── figure_generation.smk   # Publication figures
│   └── scripts/
│       ├── 01_preprocessing/
│       ├── 02_cell_annotation/
│       ├── 03_eosinophil_analysis/
│       ├── 04_trajectory_analysis/
│       ├── 05_cell_communication/
│       ├── 06_scenic_grn/
│       ├── 07_integration/
│       └── 08_figure_generation/
└── README.md
```

## Quick Start

### Prerequisites

```bash
# Install Snakemake
conda install -c bioconda snakemake

# Or using mamba (recommended)
mamba install -c bioconda snakemake
```

### Running the Pipeline

```bash
# Run complete pipeline
snakemake --cores 8

# Run specific modules
snakemake preprocessing --cores 8
snakemake eosinophil --cores 8
snakemake trajectory --cores 8
snakemake communication --cores 8
snakemake scenic --cores 8
snakemake figures --cores 8

# Dry run to see planned jobs
snakemake -n

# Generate workflow DAG
snakemake --dag | dot -Tpdf > workflow_dag.pdf
```

## Workflow Modules

### 1. Preprocessing (`preprocessing.smk`)
- **01_local_data_seurat_preprocess.R** - Load 10x FLEX data, QC filtering
- **03_UMAP.R** - PCA and UMAP dimensionality reduction
- **04_FindAllMarkers_cluster.R** - Cluster marker identification

### 2. Cell Annotation (`cell_annotation.smk`)
- **05_UCE_annotate_by_mouse_kidney_norn.py** - Universal Cell Embedding
- **06_label_UCE_identity.R** - Final cell type assignments

### 3. Eosinophil Analysis (`eosinophil_analysis.smk`)
- **07_EOS_subcluster_0.8.R** - Eosinophil subclustering (C0-C6)
- **03_EOS_DEG.R** - Differential expression (sex, genotype, infection)
- **04_EOS_ClusterProfiler.R** - GO/KEGG pathway enrichment
- **08_GSVA.R** - Gene Set Variation Analysis
- **01_define_M1_M2_signatures.R** / **02_score_cells_M1_M2.R** - Macrophage polarization

### 4. Trajectory Analysis (`trajectory_analysis.smk`)
- **02_EOS_scVelo_preparation.py** - Prepare velocity data
- **03_EOS_scVelo_dynamics.py** - RNA velocity (dynamical model)
- **02_CellRank2.py** - CellRank trajectory inference
- **04_peudotime_heatmap_1024.py** - Pseudotime visualization

### 5. Cell Communication (`cell_communication.smk`)
- **02_cellchat_score_eos-mac_condition_1021.R** - CellChat by condition
- **03_cellchat_score_eos-mac_eos-sub_1021.R** - CellChat with eos subtypes
- **04.1_CCL_pathway_heatmaps.R** - CCL signaling analysis
- **04_run_cellphonedb_analysis.py** - CellPhoneDB alternative

### 6. SCENIC GRN (`scenic_grn.smk`)
- **1201_Subcluster_Regulon_Analysis.py** - Regulon activity
- **1203_Eos_GRN_Network.py** / **1203_Mac_GRN_Network.py** - GRN visualization
- **SCENIC_RELMa_Firework_GRN.py** - RELMα-dependent networks

### 7. Integration (`integration.smk`)
- **01_process_Gurtner.R** - Process Gurtner et al. reference
- **04_harmony_integration.R** - Harmony batch correction
- **05_combined_UMAP.R** - Integrated visualization

### 8. Figure Generation (`figure_generation.smk`)
- Publication-ready figures for manuscript

## Configuration

Edit `config/config.yaml` to customize:

```yaml
# Key parameters
preprocessing:
  min_features: 200
  max_mt_percent: 10

eosinophil:
  resolution: 0.8

deg:
  test: "MAST"
  logfc_threshold: 0.25

# Cluster colors
colors:
  eosinophil_clusters:
    C0: "#6CA6CD"
    C1: "#E84A5F"
    ...
```

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
