# Upstream: Raw reads → count matrices (Cell Ranger)

This stage covers the raw-data processing that precedes the analysis pipeline:
converting raw 10x sequencing reads into per-sample gene-expression count
matrices. It is the documented origin of every downstream object in `scripts/`.

## Experimental design

Libraries were generated with **10x Genomics Chromium Fixed RNA Profiling
(FLEX)**, using the **Mouse Transcriptome Probe Set v1.0.1**
(`Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A`). The experiment
uses **4×4 probe-barcode multiplexing**: 16 samples are distributed across **4
GEM-reaction pools** (`RR0104`, `RR0508`, `RR0912`, `RR1316`), with **4 samples
per pool** distinguished by FLEX probe barcodes **BC001–BC004**.

The 16 samples span **WT/KO genotype × Male/Female sex × HFD**, with a subset of
females additionally **helminth-infected** (*N. brasiliensis*). RELMα is the
knockout (KO) gene of interest.

## Sample → pool → barcode → condition mapping

Reproduced exactly from `03_create_multiconfig_csv.sh` (the `[samples]` blocks):

| GEM pool | sample_id | probe barcode | description (condition) |
|----------|-----------|---------------|--------------------------|
| RR0104 | 1  | BC001 | WT Male HFD |
| RR0104 | 2  | BC002 | WT Female HFD |
| RR0104 | 3  | BC003 | KO Male HFD |
| RR0104 | 4  | BC004 | KO Female HFD |
| RR0508 | 5  | BC001 | WT Female HFD Infected |
| RR0508 | 6  | BC002 | WT Female HFD Infected |
| RR0508 | 7  | BC003 | KO Female HFD Infected |
| RR0508 | 8  | BC004 | KO Female HFD Infected |
| RR0912 | 9  | BC001 | WT Male HFD |
| RR0912 | 10 | BC002 | WT Male HFD |
| RR0912 | 11 | BC003 | WT Female HFD |
| RR0912 | 12 | BC004 | WT Female HFD |
| RR1316 | 13 | BC001 | KO Male HFD |
| RR1316 | 14 | BC002 | KO Male HFD |
| RR1316 | 15 | BC003 | KO Female HFD |
| RR1316 | 16 | BC004 | KO Female HFD |

## Scripts (run in order)

| # | Script | Purpose |
|---|--------|---------|
| 01 | `01_combine_fastqs.sh` | Concatenate per-lane (L007/L008) FASTQs into single per-pool R1/R2 files |
| 02 | `02_fastqc.sh` | SLURM FastQC (0.11.9) read-quality QC on the combined pool FASTQs |
| 03 | `03_create_multiconfig_csv.sh` | `cellranger multi` config CSVs for all 4 pools (reference, probe set, and the barcode→sample→condition mapping above) |
| 04 | `04_cellranger_multi.sh` | SLURM `cellranger multi` runs that demultiplex each pool into per-sample count matrices |

## Software / reference

- **Cell Ranger 8.0.1** (10x Genomics), `cellranger multi`
- **Reference genome:** `refdata-gex-mm10-2020-A` (mm10-2020-A)
- **Probe set:** `Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv`
- **FastQC 0.11.9**

## Data availability

Raw sequencing data are deposited at **GEO GSE311941**.

## Important

All paths in these scripts (FASTQ locations, reference/probe-set paths, config
CSV locations, SLURM partitions) are **specific to the UCR HPCC** environment in
which the published run was performed. They are recorded here to document the
exact run and **must be edited** to reproduce the analysis elsewhere.

## Next step

The per-sample count matrices produced by `04_cellranger_multi.sh` are the input
to the analysis pipeline. The analysis begins with Scrublet on the original
filtered count matrices in
[`../00a_doublet_detection_per_sample.py`](../00a_doublet_detection_per_sample.py),
followed by DecontX on predicted singlets in
[`../00_decontx_per_sample.R`](../00_decontx_per_sample.R). Probe-barcode
demultiplexing can resolve cross-sample multiplets, but it cannot identify all
same-sample or homotypic doublets; the computational calls therefore remain a
sensitivity analysis rather than experimental ground truth.
