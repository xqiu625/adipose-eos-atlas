#!/usr/bin/env bash
set -euo pipefail

# 02_fastqc.sh — Upstream raw-data processing, step 2 of 4
# SLURM job running FastQC (0.11.9) for read-quality QC on the combined R1/R2
# FASTQs of the 4 GEM-reaction pools produced by 01_combine_fastqs.sh.
# SLURM/module directives are HPCC-specific; edit to reproduce.
#
# Previously this listed filenames like RR13-16_S3_combined_R1.fastq.gz,
# which 01_combine_fastqs.sh never produced. Running FastQC on files that
# don't exist would just fail per-file with no files processed. This version
# targets the exact filenames/directories 01_combine_fastqs.sh actually
# establishes -- keep the two in sync. RR1316 is RR1316_S1_L001_* (not S3):
# per the GEO submission package's documented provenance, RR1316's raw data
# is RR001 (delivered under that name by the sequencing core, symlinked into
# place by 01_combine_fastqs.sh), not a combination of RR13-16_S3 lane files.
#
# ==== --outdir fix (2026-09-22, audit round 2) ===============================
# Without --outdir, FastQC writes each report next to its input FASTQ file,
# not under the versioned RERUN_ROOT -- a "legacy hard-coded root" failure
# mode caught on the first audited-rerun submission (job 29008756: confirmed
# via `scontrol show job` that WorkDir was the repo checkout, not RERUN_ROOT,
# and its reports landed next to the raw FASTQs). That run was already 7+
# hours in when caught and was NOT restarted (killing it would waste real
# compute for a pure bookkeeping fix); as of this fix the job is still
# RUNNING and no reports have been relocated yet. Its reports will be copied
# after successful completion into RERUN_ROOT using a checksum-verified
# transfer; the deviation and the transfer results will be recorded in
# RERUN_LOG.md once that actually happens. This version fixes the script
# itself so every subsequent submission writes directly under RERUN_ROOT.

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
