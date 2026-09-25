#!/usr/bin/env bash
# Submitted only by submit_preprocessing_array.sh, from its frozen snapshot.
set -euo pipefail
stage="$1"
base="$2"
action="$3"
python_executable="$4"
script_dir="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
threads="${SLURM_CPUS_PER_TASK:?Missing SLURM_CPUS_PER_TASK}"
export OMP_NUM_THREADS="$threads" OPENBLAS_NUM_THREADS="$threads"
export MKL_NUM_THREADS="$threads" NUMEXPR_NUM_THREADS="$threads"
export EOS_ATLAS_BASE="$base"
if [[ "$stage" == decontx && "$action" == worker ]]; then
  : "${RERUN_R_LIB:?Set the validated pinned R library}"
  : "${CMDSTAN:?Set the validated CmdStan directory}"
  test -d "$RERUN_R_LIB"
  test -d "$CMDSTAN"
  # sbatch --export=ALL inherits login-shell modules (including R/4.5.2).
  # Clear them before loading the audited compiler/R pair; retain hard checks.
  module purge
  module load gcc/12.2.0
  module load R/4.5.0
  module list
  Rscript --vanilla -e 'stopifnot(as.character(getRversion()) == "4.5.0"); sessionInfo()'
fi
args=("$script_dir/parallel_preprocessing.py" "$action" --stage "$stage" --base "$base")
if [[ "$action" == worker ]]; then
  args+=(--sample-id "${SLURM_ARRAY_TASK_ID:?Missing array sample ID}")
fi
exec "$python_executable" "${args[@]}"
