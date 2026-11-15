#!/usr/bin/env bash
set -euo pipefail

# Run FastQC 0.11.9 on the four configured pool FASTQ pairs.
# Write reports under RERUN_ROOT.

#SBATCH --ntasks=10
#SBATCH --mem=256gb
#SBATCH --time=2-00:00:00
#SBATCH -p batch
#SBATCH --job-name="fastqc"

module load fastqc

: "${RERUN_ROOT:?RERUN_ROOT must be set, e.g. export RERUN_ROOT=/bigdata/nairlab/shared/XQiu/EOS_project/reruns/doublet_YYYYMMDD_HHMMSS}"
RR1316_DIR="/bigdata/godziklab/shared/Xinru/302025/dataset/batch3/raw_data"
SHARED_DIR="/bigdata/nairlab/shared/XQiu/EOS_project/dataset/raw_data"  # RR0104, RR0508, RR0912
OUT_DIR="$RERUN_ROOT/upstream_cellranger/fastqc"
mkdir -p "$OUT_DIR"

fastqc --outdir "$OUT_DIR" "$RR1316_DIR/RR1316_S1_L001_R1_001.fastq.gz"
fastqc --outdir "$OUT_DIR" "$RR1316_DIR/RR1316_S1_L001_R2_001.fastq.gz"
fastqc --outdir "$OUT_DIR" "$SHARED_DIR/RR0104_S1_L001_R1_001.fastq.gz"
fastqc --outdir "$OUT_DIR" "$SHARED_DIR/RR0104_S1_L001_R2_001.fastq.gz"
fastqc --outdir "$OUT_DIR" "$SHARED_DIR/RR0508_S4_L001_R1_001.fastq.gz"
fastqc --outdir "$OUT_DIR" "$SHARED_DIR/RR0508_S4_L001_R2_001.fastq.gz"
fastqc --outdir "$OUT_DIR" "$SHARED_DIR/RR0912_S2_L001_R1_001.fastq.gz"
fastqc --outdir "$OUT_DIR" "$SHARED_DIR/RR0912_S2_L001_R2_001.fastq.gz"

# fastqc/0.11.9
