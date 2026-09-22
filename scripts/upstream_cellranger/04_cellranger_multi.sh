#!/usr/bin/env bash
set -euo pipefail

# 04_cellranger_multi.sh — Upstream raw-data processing, step 4 of 4
# Submits SLURM `cellranger multi` (Cell Ranger 8.0.1) jobs that demultiplex
# each of the 4 GEM-reaction pools (RR0104, RR0508, RR0912, RR1316) into
# per-sample count matrices, using the config CSVs from
# 03_create_multiconfig_csv.sh. Output count matrices feed the analysis
# pipeline starting at 00_decontx_per_sample.R.
# Both standard and highmem SLURM variants are recorded; paths are HPCC-specific.
#
# Previously this file concatenated 6 separate `#SBATCH ...` blocks followed
# by `cellranger multi` calls into one script. `sbatch` only honors #SBATCH
# directives that appear before the first executable line, so submitting this
# file as a single job (`sbatch 04_cellranger_multi.sh`) ran only the FIRST
# block's resource request and then executed every `cellranger multi` call
# sequentially inside that one allocation -- silently ignoring the highmem/
# longer-time requests for the later pools. It also reused --id=RR1316 and
# --id=RR0104 for both the standard and "highmem" runs of the same pool,
# which collide on the same pipestance output directory. This version submits
# one independent SLURM job per run via `sbatch`, each with its own resource
# request, and gives the highmem reruns a distinct --id/output directory.

CONFIG_DIR="/bigdata/godziklab/shared/Xinru/302025/dataset/batch3/csv_config"

submit() {
  local job_name="$1" ntasks="$2" mem="$3" time="$4" partition="$5"
  local run_id="$6" csv="$7"
  sbatch --job-name="$job_name" --ntasks="$ntasks" --mem="$mem" \
    --time="$time" -p "$partition" \
    --wrap="cellranger multi --id=$run_id --csv=$csv"
}

# RR1316
submit "RR1316"  10 256gb 4-00:00:00 batch    RR1316          "$CONFIG_DIR/RR1316_config.csv"
# RR1316 highmem rerun -- distinct --id so it writes to its own pipestance dir
submit "hRR1316" 16 512gb 4-00:00:00 highmem  RR1316_highmem  "$CONFIG_DIR/RR1316_config.csv"

# RR0104
submit "RR0104"  10 256gb 4-00:00:00 intel    RR0104          "$CONFIG_DIR/RR0104_config.csv"
# RR0104 highmem rerun -- distinct --id so it writes to its own pipestance dir
submit "hRR0104" 16 512gb 4-00:00:00 highmem  RR0104_highmem  "$CONFIG_DIR/RR0104_config.csv"

# RR0508
submit "RR0508"  10 256gb  4-00:00:00  intel  RR0508 "$CONFIG_DIR/RR0508_config.csv"

# RR0912
submit "RR0912"  16 256gb 10-00:00:00  intel  RR0912 "$CONFIG_DIR/RR0912_config.csv"

echo "Submitted 6 cellranger multi jobs."
