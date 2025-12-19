# Adipose Eosinophil and Macrophage Single-Cell Atlas

Single-cell RNA sequencing analysis of eosinophils and macrophages in adipose tissue under metabolic stress conditions, with focus on RELMα function in obesity.

## Overview

This project investigates immune cell heterogeneity and cell-cell communication in adipose tissue using scRNA-seq. We characterize eosinophil and macrophage subpopulations across different metabolic conditions and genotypes.

### Study Design

| Factor | Conditions |
|--------|------------|
| Genotype | Wild-Type (WT), RELMα knockout (KO) |
| Sex | Male, Female |
| Diet | High-Fat Diet (HFD) |
| Challenge | HFD + Infection |

## Repository Structure

```
├── scripts/
│   ├── R01_initial_analysis/     # Core analysis pipeline
│   │   ├── Eosinophil/           # Eosinophil-specific analyses
│   │   ├── R/                    # Integration and preprocessing
│   │   └── python/               # Python utilities
│   ├── R05-R21/                  # Progressive analyses and figures
│   └── ...
├── data/
│   └── public_datasets/          # Published reference datasets
├── results/                      # Output files and figures
└── Memo/                         # Analysis notes
```

## Analysis Pipeline

### 1. Preprocessing & Integration
- Quality control and filtering
- Batch correction (Harmony)
- Clustering and annotation (Seurat)

### 2. Cell Type Analysis

**Eosinophils (~10,700 cells)**
- 7 subclusters (Eos_C0–C6)
- Differential expression (MAST)
- RNA velocity (scVelo)
- Trajectory analysis (Monocle3, CellRank2)
- Gene regulatory networks (pySCENIC)

**Macrophages**
- Subtype classification (LAM, ProInflam, PVM, Metabolic, etc.)
- M1/M2 polarization scoring
- Trajectory and regulon analysis

### 3. Cell-Cell Communication
- CellChat pathway prediction
- CellPhoneDB analysis
- CCL ligand-receptor interactions
- Eosinophil-macrophage signaling networks

## Key Tools & Dependencies

### R
- Seurat (v4+)
- clusterProfiler
- monocle3
- CellChat
- harmony
- GSVA

### Python
- scanpy
- scvelo
- cellrank
- pyscenic
- cellphonedb

## Cluster Definitions

### Eosinophil Subclusters
| Cluster | Color | Proportion |
|---------|-------|------------|
| C0 | #6CA6CD | 28.9% |
| C1 | #E84A5F | 17.8% |
| C2 | #FF6F1F | 17.1% |
| C3 | #CBAACB | 12.5% |
| C4 | #B03A48 | 10.7% |
| C5 | #5BA4A4 | 10.5% |
| C6 | #2F3C7E | 2.5% |

### Macrophage Subtypes
- Mac_Cycling
- Mac_LAM (lipid-associated)
- Mac_Metabolic
- Mac_Mono
- Mac_ProInflam
- Mac_PVM (perivascular)

## Key Output Files

- `Eosinophil_enriched_genes.csv` - Eosinophil marker genes
- `CCL_Eosinophil_to_Macrophage_*.csv` - Interaction matrices
- SCENIC regulon files
- CellChat/CellPhoneDB results
- Publication-ready figures (PDF/SVG/EPS)

## Usage

### Running Eosinophil Analysis
```bash
cd scripts/R01_initial_analysis/Eosinophil
Rscript 02_EOS_proportion.R
Rscript 03_EOS_DEG.R
Rscript 04_EOS_ClusterProfiler.R
```

### Running Trajectory Analysis
```bash
python 05_scVelo/velocity_analysis.py
python scripts/R15_*/04_peudotime_heatmap_1024.py
```

### Running Cell-Cell Communication
```bash
Rscript 11_CellphoneDB/run_cellchat.R
```

## Code Availability

All custom scripts and analysis code used in this study are publicly available at https://github.com/xqiu625/adipose-eos-subtype-atlas.

## Data Availability

Primary scRNA-seq data: 10x Chromium 3' Gene Expression

Reference datasets:
- Felton et al. (adipose eosinophils)
- GSE248899
- Gurtner et al. (blood and bone marrow eosinophils)

## Citation

[Add citation when published]

## Contact

[Add contact information]
