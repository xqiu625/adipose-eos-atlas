#!/usr/bin/env bash
# One audited stage per invocation; dry run unless --submit is explicit.
set -euo pipefail
if [[ $# -lt 3 || $# -gt 4 ]]; then
  echo "Usage: bash $0 scrublet|decontx INPUT_MANIFEST REVIEWED_GATE_JSON [--submit]" >&2
  exit 2
fi
stage="$1"
manifest="$2"
gate="$3"
mode="${4:---dry-run}"
case "$stage" in scrublet|decontx) ;; *) echo 'Invalid stage' >&2; exit 2;; esac
case "$mode" in --submit|--dry-run) ;; *) echo 'Invalid mode' >&2; exit 2;; esac
: "${RERUN_ROOT:?Set the versioned rerun root}"
: "${EOS_ATLAS_BASE:?Set the analysis base within RERUN_ROOT}"
: "${PYTHON_EXECUTABLE:?Set the absolute Python executable in the pinned analysis environment}"
[[ "$PYTHON_EXECUTABLE" == /* && -x "$PYTHON_EXECUTABLE" ]]
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
partition="${PREPROCESS_PARTITION:-batch}"
cpus="${PREPROCESS_CPUS:-2}"
mem="${PREPROCESS_MEM:-64G}"
limit="${PREPROCESS_TIME:-04:00:00}"
if [[ "$stage" == scrublet ]]; then
  concurrency="${PREPROCESS_CONCURRENCY:-8}"
else
  concurrency="${PREPROCESS_CONCURRENCY:-4}"
  : "${RERUN_R_LIB:?Set the pinned R library}"
  : "${CMDSTAN:?Set the pinned CmdStan directory}"
fi
[[ "$concurrency" =~ ^([1-9]|1[0-6])$ && "$cpus" =~ ^[1-9][0-9]*$ ]] || {
  echo 'Concurrency must be 1–16; CPUs must be a positive integer' >&2; exit 2;
}
prepare=("$PYTHON_EXECUTABLE" "$script_dir/parallel_preprocessing.py" prepare
  --stage "$stage" --base "$EOS_ATLAS_BASE" --manifest "$manifest" --gate "$gate"
  --expected-rate "${SCRUBLET_EXPECTED_RATE:-0.06}" --seed "${SCRUBLET_SEED:-42}")
if [[ "$mode" != --submit ]]; then
  printf 'Preparation command (not executed): '; printf '%q ' "${prepare[@]}"; printf '\n'
  printf 'Plan: %s, array=1-16%%%s, cpus/task=%s, mem/task=%s, time=%s, partition=%s\n' \
    "$stage" "$concurrency" "$cpus" "$mem" "$limit" "$partition"
  echo 'Then submit one afterok collector; stop for scientific QC before the next stage.'
  echo 'No files written, inputs/gates not yet validated, no jobs submitted.'
  exit 0
fi
command -v sbatch >/dev/null
state="$("${prepare[@]}")"
base="$(cd -- "$EOS_ATLAS_BASE" && pwd -P)"
runner="$state/run_preprocessing_task.sh"
array_id="$(sbatch --parsable --export=ALL --job-name="${stage}_samples" \
  --partition="$partition" --ntasks=1 --cpus-per-task="$cpus" --mem="$mem" --time="$limit" \
  --array="1-16%${concurrency}" --chdir="$base" \
  --output="$state/logs/%A_%a.out" --error="$state/logs/%A_%a.err" \
  "$runner" "$stage" "$base" worker "$PYTHON_EXECUTABLE")"
array_id="${array_id%%;*}"
[[ "$array_id" =~ ^[0-9]+$ ]]
printf 'array\t%s\n' "$array_id" > "$state/job_ids.tsv"
echo "Array submitted: $array_id (submission only, not PASS)"
collector_id="$(sbatch --parsable --export=ALL --job-name="${stage}_collect" \
  --dependency="afterok:${array_id}" --partition="$partition" --ntasks=1 \
  --cpus-per-task=1 --mem=8G --time=01:00:00 --chdir="$base" \
  --output="$state/logs/collect_%j.out" --error="$state/logs/collect_%j.err" \
  "$runner" "$stage" "$base" collect "$PYTHON_EXECUTABLE")"
collector_id="${collector_id%%;*}"
[[ "$collector_id" =~ ^[0-9]+$ ]]
printf 'collector\t%s\n' "$collector_id" >> "$state/job_ids.tsv"
echo "Collector submitted: $collector_id; record sacct and QC checkpoint after completion."
echo 'No dependent scientific stage was submitted.'
