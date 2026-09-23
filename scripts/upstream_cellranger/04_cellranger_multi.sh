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
# request.
#
# ==== RERUN_ROOT + redundant-job fix (2026-09-22) ============================
# Two fixes for the audited rerun:
# 1. `cellranger multi --id=X` writes its pipestance into the CURRENT
#    directory, so submitting via a bare --wrap from an arbitrary cwd landed
#    output in whatever directory sbatch happened to be invoked from -- a
#    "legacy hard-coded root" failure mode. This version requires RERUN_ROOT
#    to be set and --chdir's every job into
#    $RERUN_ROOT/upstream_cellranger/cellranger_multi so output always
#    resolves under the versioned rerun root.
# 2. Dropped the redundant standard-resource RR1316/RR0104 submissions. The
#    original script submitted BOTH a standard-resource run and a "highmem
#    rerun" of the same pool under a different --id; the "rerun" naming
#    implies the standard-resource attempt was insufficient (OOM), so
#    resubmitting it here would just spend ~4 days of batch-partition time on
#    a run likely to fail the same way. This version submits exactly one
#    (the higher-resource) job per pool -- 4 jobs total.
#
# ==== CPU allocation + binary pinning fix (2026-09-22, audit round 2) =======
# Two more fixes, applied before these jobs' first real submission (caught
# while still PENDING behind FastQC, no compute wasted):
# 3. `--ntasks=N` requests N MPI tasks, but each job runs exactly one
#    `cellranger multi` process -- Cell Ranger parallelizes internally via
#    threads, not via SLURM tasks. Changed to `--ntasks=1 --cpus-per-task=N`,
#    and `cellranger multi` is now called with explicit `--localcores` (tied
#    to $SLURM_CPUS_PER_TASK) and `--localmem` (a few GB under the SLURM
#    --mem request, so Cell Ranger's own resource autodetection can't
#    overcommit relative to what SLURM actually granted).
# 4. Jobs previously called bare `cellranger multi`, relying on PATH/module
#    state at run time to resolve to the login-node-verified 8.0.1 build --
#    not guaranteed on a compute node. Now calls the MD5-verified 8.0.1
#    binary by absolute path, and the wrapped command runs
#    `cellranger --version` first and hard-fails (nonzero exit, before
#    `multi` ever runs) if it doesn't report cellranger-8.0.1.

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

# RR1316 -- highmem resources (standard-resource attempt is documented as
# needing a rerun; submit only the resource level known to be sufficient).
# 16 cores, 512gb SLURM allocation, 480gb passed to --localmem (headroom for
# OS/SLURM overhead).
submit "RR1316" 16 512gb 480 4-00:00:00 highmem RR1316 "$CONFIG_DIR/RR1316_config.csv"

# RR0104 -- highmem resources, same reasoning
submit "RR0104" 16 512gb 480 4-00:00:00 highmem RR0104 "$CONFIG_DIR/RR0104_config.csv"

# RR0508 -- 10 cores, 256gb SLURM allocation, 230gb to --localmem
submit "RR0508" 10 256gb 230 4-00:00:00 intel   RR0508 "$CONFIG_DIR/RR0508_config.csv"

# RR0912 -- 16 cores, 256gb SLURM allocation, 230gb to --localmem
submit "RR0912" 16 256gb 230 10-00:00:00 intel   RR0912 "$CONFIG_DIR/RR0912_config.csv"

echo "Submitted 4 cellranger multi jobs, output under $OUT_DIR, pinned to $CELLRANGER_BIN"
