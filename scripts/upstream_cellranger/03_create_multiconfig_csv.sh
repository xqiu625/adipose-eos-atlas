#!/usr/bin/bash

# 03_create_multiconfig_csv.sh — Upstream raw-data processing, step 3 of 4
# Documents the `cellranger multi` config CSVs for all 4 GEM-reaction pools
# (RR0104, RR0508, RR0912, RR1316) of the 10x FLEX 4x4 multiplexed design,
# mapping each FLEX probe barcode (BC001-BC004) to its sample and condition.
# Reference mm10-2020-A, probe set Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1.
# Paths are HPCC-specific to the GSE311941 run; edit to reproduce.

# https://www.10xgenomics.com/support/software/cell-ranger/latest/analysis/running-pipelines/cr-flex-multi-frp#mfrp-onebc

# Here is the CSV config template used in this study

[gene-expression]
reference,/path/to/transcriptome
probe-set,/path/to/probe-set.csv  #e.g., cellranger-x.y.z/probe_sets/Chromium_Human_Transcriptome_Probe_Set_v1.0.1_GRCh38-2020-A.csv
create-bam,false #do not generate BAM file

[libraries]
fastq_id,fastqs,feature_types
frp_gex,/path/to/fastqs,Gene Expression

[samples]
sample_id,probe_barcode_ids,description
sample1,BC001,Control
sample2,BC003,Treated


## RR1316
vi /bigdata/godziklab/shared/Xinru/302025/dataset/batch3/csv_config/RR1316_config.csv

[gene-expression]
reference,/bigdata/godziklab/shared/Xinru/yard/run_cellranger_count/refdata-gex-mm10-2020-A
probe-set,/bigdata/godziklab/shared/Xinru/yard/apps/cellranger-8.0.1/probe_sets/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv
create-bam,true

[libraries]
fastq_id,fastqs,feature_types
RR1316,/bigdata/godziklab/shared/Xinru/302025/dataset/batch3/raw_data,Gene Expression

[samples]
sample_id,probe_barcode_ids,description
13,BC001,KO Male HFD
14,BC002,KO Male HFD
15,BC003,KO Female HFD
16,BC004,KO Female HFD

## RR0104
vi /bigdata/godziklab/shared/Xinru/302025/dataset/batch3/csv_config/RR0104_config.csv

[gene-expression]
reference,/bigdata/godziklab/shared/Xinru/yard/run_cellranger_count/refdata-gex-mm10-2020-A
probe-set,/bigdata/godziklab/shared/Xinru/yard/apps/cellranger-8.0.1/probe_sets/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv
create-bam,true

[libraries]
fastq_id,fastqs,feature_types
RR0104,/bigdata/nairlab/shared/XQiu/EOS_project/dataset/raw_data,Gene Expression

[samples]
sample_id,probe_barcode_ids,description
1,BC001,WT Male HFD
2,BC002,WT Female HFD
3,BC003,KO Male HFD
4,BC004,KO Female HFD


## RR0508

vi /bigdata/godziklab/shared/Xinru/302025/dataset/batch3/csv_config/RR0508_config.csv

[gene-expression]
reference,/bigdata/godziklab/shared/Xinru/yard/run_cellranger_count/refdata-gex-mm10-2020-A
probe-set,/bigdata/godziklab/shared/Xinru/yard/apps/cellranger-8.0.1/probe_sets/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv
create-bam,true

[libraries]
fastq_id,fastqs,feature_types
RR0508,/bigdata/nairlab/shared/XQiu/EOS_project/dataset/raw_data,Gene Expression

[samples]
sample_id,probe_barcode_ids,description
5,BC001,WT Female HFD Infected
6,BC002,WT Female HFD Infected
7,BC003,KO Female HFD Infected
8,BC004,KO Female HFD Infected


## RR0912

[gene-expression]
reference,/bigdata/godziklab/shared/Xinru/yard/run_cellranger_count/refdata-gex-mm10-2020-A
probe-set,/bigdata/godziklab/shared/Xinru/yard/apps/cellranger-8.0.1/probe_sets/Chromium_Mouse_Transcriptome_Probe_Set_v1.0.1_mm10-2020-A.csv
create-bam,true

[libraries]
fastq_id,fastqs,feature_types
RR0912,/bigdata/nairlab/shared/XQiu/EOS_project/dataset/raw_data,Gene Expression

[samples]
sample_id,probe_barcode_ids,description
9,BC001,WT Male HFD
10,BC002,WT Male HFD
11,BC003,WT Female HFD
12,BC004,WT Female HFD
