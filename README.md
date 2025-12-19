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



[Add citation when published]

## Contact

[Add contact information]
