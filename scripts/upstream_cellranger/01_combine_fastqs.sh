#!/usr/bin/bash

# 01_combine_fastqs.sh — Upstream raw-data processing, step 1 of 4
# Concatenates the per-lane (L007/L008) FASTQs for each of the 4 GEM-reaction
# pools (RR0104, RR0508, RR0912, RR1316) into single per-pool R1/R2 files,
# ready for FastQC (02) and `cellranger multi` (04).
# Paths/filenames are HPCC-specific to the GSE311941 run; edit to reproduce.

cat RR13-16_S3_L007_R1_001.fastq.gz RR13-16_S3_L008_R1_001.fastq.gz > RR1316_S3_L001_R1_001.fastq.gz
cat RR13-16_S3_L007_R2_001.fastq.gz RR13-16_S3_L008_R2_001.fastq.gz > RR1316_S3_L001_R2_001.fastq.gz

cat RR1-4_S1_L007_R1_001.fastq.gz RR1-4_S1_L008_R1_001.fastq.gz > RR0104_combined_R1_001.fastq.gz
cat RR1-4_S1_L007_R2_001.fastq.gz RR1-4_S1_L008_R2_001.fastq.gz > RR0104_combined_R2_001.fastq.gz


cat RR5-8_S4_L007_R1_001.fastq.gz RR5-8_S4_L008_R1_001.fastq.gz > RR0508_S4_L001_R1_001.fastq.gz
cat RR5-8_S4_L007_R2_001.fastq.gz RR5-8_S4_L008_R2_001.fastq.gz > RR0508_S4_L001_R2_001.fastq.gz


cat RR9-12_S2_L007_R1_001.fastq.gz RR9-12_S2_L008_R1_001.fastq.gz > RR0912_S2_L001_R1_001.fastq.gz
cat RR9-12_S2_L007_R2_001.fastq.gz RR9-12_S2_L008_R2_001.fastq.gz > RR0912_S2_L001_R2_001.fastq.gz


cat RR0104_combined_R1_001.fastq.gz RR001_S1_L001_R1_001.fastq.gz > RR0104_S1_L001_R1_001.fastq.gz
cat RR0104_combined_R2_001.fastq.gz RR001_S1_L001_R2_001.fastq.gz > RR0104_S1_L001_R2_001.fastq.gz
