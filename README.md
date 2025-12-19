
Single-cell RNA sequencing analysis of eosinophils and macrophages in adipose tissue under metabolic stress conditions, with focus on RELMα function in obesity.

## Overview



### Study Design

| Factor | Conditions |
|--------|------------|
| Genotype | Wild-Type (WT), RELMα knockout (KO) |
| Sex | Male, Female |
| Diet | High-Fat Diet (HFD) |
| Challenge | HFD + Infection |

## Repository Structure

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
