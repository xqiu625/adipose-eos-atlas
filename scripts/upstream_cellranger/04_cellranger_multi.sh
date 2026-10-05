#!/usr/bin/env bash
set -euo pipefail

# Submit one Cell Ranger 8.0.1 job per pool.
# Set RERUN_ROOT; jobs use separate output directories and explicit resources.

: "${RERUN_ROOT:?RERUN_ROOT must be set, e.g. export RERUN_ROOT=/bigdata/nairlab/shared/XQiu/EOS_project/reruns/doublet_YYYYMMDD_HHMMSS}"
CONFIG_DIR="/bigdata/godziklab/shared/Xinru/302025/dataset/batch3/csv_config"
OUT_DIR="$RERUN_ROOT/upstream_cellranger/cellranger_multi"
mkdir -p "$OUT_DIR"

# MD5-verified Cell Ranger 8.0.1 binary (tarball md5 1d08d8005b026fc1ca3c867bc67b2dbd,
# matches 10x Genomics' documented value for this release), pinned by absolute
# path so a batch job can never silently resolve a different version.
CELLRANGER_BIN="/bigdata/godziklab/shared/Xinru/yard/apps/cellranger-8.0.1/cellranger"

# Optional: set FASTQC_JOBID to the 02_fastqc.sh SLURM job ID to make every
# cellranger multi job wait on it (afterok) as a raw-read QC gate before
# committing multi-day compute.
DEP_ARGS=()
if [[ -n "${FASTQC_JOBID:-}" ]]; then
  DEP_ARGS=(--dependency="afterok:${FASTQC_JOBID}")
fi

submit() {
  local job_name="$1" cpus="$2" mem="$3" localmem_gb="$4" time="$5" partition="$6"
  local run_id="$7" csv="$8"
  local wrap_cmd
  wrap_cmd=$(cat <<WRAP
set -euo pipefail
VER=\$('$CELLRANGER_BIN' --version 2>&1)
echo "cellranger --version: \$VER"
if [[ "\$VER" != *cellranger-8.0.1* ]]; then
  echo "FATAL: pinned cellranger binary did not report version 8.0.1 -- refusing to run" >&2
  exit 1
fi
'$CELLRANGER_BIN' multi --id=$run_id --csv=$csv --localcores=\$SLURM_CPUS_PER_TASK --localmem=$localmem_gb
WRAP
)
  sbatch --job-name="$job_name" --ntasks=1 --cpus-per-task="$cpus" --mem="$mem" \
    --time="$time" -p "$partition" --chdir="$OUT_DIR" "${DEP_ARGS[@]}" \
    --wrap="$wrap_cmd"
}

# RR1316: 16 cores, 512 GB allocation, 480 GB for Cell Ranger.
submit "RR1316" 16 512gb 480 4-00:00:00 highmem RR1316 "$CONFIG_DIR/RR1316_config.csv"

# RR0104 -- highmem resources, same reasoning
submit "RR0104" 16 512gb 480 4-00:00:00 highmem RR0104 "$CONFIG_DIR/RR0104_config.csv"

# RR0508 -- 10 cores, 256gb SLURM allocation, 230gb to --localmem
submit "RR0508" 10 256gb 230 4-00:00:00 intel   RR0508 "$CONFIG_DIR/RR0508_config.csv"

# RR0912 -- 16 cores, 256gb SLURM allocation, 230gb to --localmem
submit "RR0912" 16 256gb 230 10-00:00:00 intel   RR0912 "$CONFIG_DIR/RR0912_config.csv"

echo "Submitted 4 cellranger multi jobs, output under $OUT_DIR, pinned to $CELLRANGER_BIN"
