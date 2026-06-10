#!/usr/bin/bash

# 02_fastqc.sh — Upstream raw-data processing, step 2 of 4
# SLURM job running FastQC (0.11.9) for read-quality QC on the combined R1/R2
# FASTQs of the 4 GEM-reaction pools produced by 01_combine_fastqs.sh.
# SLURM/module directives are HPCC-specific; edit to reproduce.

#SBATCH --ntasks=10
#SBATCH --mem=256gb
#SBATCH --time=2-00:00:00
#SBATCH -p batch
#SBATCH --job-name="fastqc"

module load fastqc
fastqc RR13-16_S3_combined_R1.fastq.gz
fastqc RR13-16_S3_combined_R2.fastq.gz
fastqc RR1-4_S1_combined_R1.fastq.gz
fastqc RR1-4_S1_combined_R2.fastq.gz
fastqc RR5-8_S4_combined_R1.fastq.gz
fastqc RR5-8_S4_combined_R2.fastq.gz
fastqc RR9-12_S2_combined_R1.fastq.gz
fastqc RR9-12_S2_combined_R2.fastq.gz

# fastqc/0.11.9
