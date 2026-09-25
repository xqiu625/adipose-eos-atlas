# Parallel preprocessing on HPCC

The audited entry point is `scripts/submit_preprocessing_array.sh`. It submits
one scientific stage and its collector, not the entire downstream pipeline.
Scrublet and DecontX now accept `--sample-id 1..16` and `--input-manifest PATH`.
Without these flags, their existing sequential interfaces remain available.
Scientific thresholds and Scrublet's per-sample seed are unchanged. DecontX's
existing function defaults are unchanged; verify serial/array equivalence in
the pinned HPCC environment before production use.

## Prerequisites and input identity

Read `CLAUDE.md` and `note.md`. Finish the four Cell Ranger runs, the web-summary
QC gate, RR1316 exact/functional reproduction, and the code/input audit first.
The GEO correction is a separate administrative task, deferred until **after
2026-10-07** (earliest 2026-10-08); it does not block a validated analysis branch.
Do not upload or contact GEO through this workflow.

Use a new analysis base inside the existing `RERUN_ROOT`. Pass a reviewed TSV
with exactly 16 rows and these required columns:

```text
sample_id	pool	filtered_dir	raw_dir	qc_status	provenance_status
```

`filtered_dir` and `raw_dir` are absolute paths to the Cell Ranger filtered and
raw matrix directories, each containing nonempty `matrix.mtx.gz`,
`barcodes.tsv.gz`, and `features.tsv.gz`. Both status columns must be `PASS`.
Also required are six preflight SHA-256 columns, with lowercase hex values:
`filtered_matrix_sha256`, `filtered_barcodes_sha256`,
`filtered_features_sha256`, `raw_matrix_sha256`, `raw_barcodes_sha256`,
`raw_features_sha256`. Generate these during the input audit; workers compare
actual file hashes against them on compute nodes before analysis. Additional
columns are allowed. Resolved matrix
directories must be within `RERUN_ROOT`; symlinks to historical outputs are
rejected. Samples 13–16 must resolve through `RR1316_L007L008`, explicitly
excluding the erroneous RR001-as-RR1316 branch. A directory name alone is not
proof of identity: the reviewed manifest and reproduction checkpoint supply
that evidence.

The gate JSON must reference an actual completed audit, never a placeholder:

```json
{
  "stage": "cellranger_qc_and_00a_audit",
  "status": "PASS",
  "manifest_sha256": "SHA256_OF_THE_EXACT_TSV",
  "script_sha256": "SHA256_OF_00a_doublet_detection_per_sample.py",
  "checkpoint": "PATH_OR_ID_OF_COMPLETED_CHECKPOINT",
  "reviewed_by": "ACTUAL_REVIEWER"
}
```

The launcher validates the schema and hash bindings; it cannot authenticate a
reviewer's judgment or check the truth of a manually supplied PASS. A
PASS_WITH_WARNING requires investigator approval and a recorded release before
providing a PASS gate for execution.

## Submit Scrublet

```bash
export RERUN_ROOT=/absolute/path/to/current/rerun
export EOS_ATLAS_BASE="$RERUN_ROOT/analysis_parallel_attempt_01"
export PYTHON_EXECUTABLE=/absolute/path/to/pinned/conda/env/bin/python
export RERUN_R_LIB="$RERUN_ROOT/r_env/renv_library"
export CMDSTAN="$RERUN_ROOT/r_env/cmdstan/cmdstan-2.40.0"

# Prints plan only. Does not create files, validate gates, or submit jobs.
bash scripts/submit_preprocessing_array.sh scrublet /absolute/manifest.tsv /absolute/00a_gate.json

# Only after the actual gate is approved:
bash scripts/submit_preprocessing_array.sh scrublet /absolute/manifest.tsv /absolute/00a_gate.json --submit
```

Defaults are 2 CPUs / 64 GB / 4 hours per task, `batch`, array `1-16%8`.
These are starting resource requests, not measured runtime predictions.
`PREPROCESS_CPUS`, `PREPROCESS_MEM`, `PREPROCESS_TIME`, `PREPROCESS_PARTITION`,
and `PREPROCESS_CONCURRENCY` override them. The account's limits, queue and
shared filesystem throughput may require lower concurrency. Thread counts are
bounded by the allocated CPUs. `SCRUBLET_EXPECTED_RATE` (default 0.06) and
`SCRUBLET_SEED` (42) are frozen at preparation; sensitivity runs need separate
analysis bases. Verify that the chosen Python executable loads the pinned
packages correctly in batch; use an already activated environment if its
activation scripts supply required runtime settings.

Preparation refuses an existing stage output directory, snapshots code,
manifest and gate under `parallel/scrublet/`, and records Git commit plus dirty
status. Each worker exclusively claims its sample, hashes its actual matrix
inputs, and writes a success receipt only after its subprocess exits zero and
all outputs exist. Logs use `%A_%a`; receipt filenames include the sample ID.
Array tasks never write shared aggregate outputs.

One collector uses `afterok:<array_job_id>` and verifies all 16 receipts,
output hashes, sample identities, duplicate barcodes, finite scores/thresholds,
call counts, package versions and analysis parameters. Its final status is
`COLLECTED_PENDING_SCIENTIFIC_QC`, **not PASS**. It writes the same shared
`doublets/scrublet_summary.csv` and run-parameters filename used downstream.

Record `sacct` for the array, its tasks and collector; inspect warnings and the
scientific metrics in `note.md`. If the collector submission fails, the array
may still be running: its ID has already been saved in `job_ids.tsv`. Inspect
and resubmit just the collector; do not invoke preparation again.

## Submit DecontX after Scrublet review

Use the same analysis base and manifest. Supply a new gate with:

- `stage`: `scrublet_qc_and_decontx_audit`;
- `script_sha256`: SHA-256 of `00_decontx_per_sample.R`;
- `collection_sha256`: SHA-256 of `parallel/scrublet/collection.json`;
- the same manifest hash, reviewed PASS, checkpoint and reviewer fields.

```bash
bash scripts/submit_preprocessing_array.sh decontx /absolute/manifest.tsv /absolute/00_gate.json
bash scripts/submit_preprocessing_array.sh decontx /absolute/manifest.tsv /absolute/00_gate.json --submit
```

Default concurrency is 4, with the same initial per-task resources. The worker
loads gcc/12.2.0 and R/4.5.0, checks the R version, requires `RERUN_R_LIB` and
`CMDSTAN`, then executes the existing library-resolution guard. Module conflicts
stop the job. Validate the pinned environment in a real batch smoke test.

DecontX requires the completed, unchanged Scrublet collection. Each task writes
one RDS and one sample scores CSV. The collector checks that every sample's
barcodes equal its predicted singlets and contamination values are finite and
in [0,1], then writes `decontx/decontx_all_scores.csv`. `01_qc_doublet_filter.R`
can consume these outputs using the same `EOS_ATLAS_BASE`, after the required
checkpoint. No QC/UMAP/annotation jobs are automatically submitted.

## Failure, testing and performance

Workers and collectors use persistent exclusive lock directories. Duplicate
tasks, missing samples, changed outputs and incomplete collections fail closed.
A failed sample can leave partial files; preserve them for diagnosis. This
initial implementation deliberately does not auto-resume or overwrite them.
Use a fresh attempt base and reviewed rerun plan; do not remove a lock to retry
over existing output. For a fresh DecontX-only attempt, any reuse of validated
Scrublet artifacts must be explicitly arranged and recorded; the launcher does
not copy them automatically.

Local checks:

```bash
python3 -m unittest discover -s tests -p 'test_parallel_preprocessing.py' -v
bash -n scripts/submit_preprocessing_array.sh scripts/run_preprocessing_task.sh
Rscript -e 'invisible(parse(file="scripts/00_decontx_per_sample.R"))'
```

Before production, compare a representative sample's sequential and single-task
outputs under the pinned HPCC packages. Repeat with a high-cell-count sample to
measure memory. The local scheduling tests use synthetic files and do not
establish scientific equivalence or real resource requirements.

Use `sacct` task/step records for `Elapsed`, `TotalCPU`, `AllocCPUS`, `MaxRSS`,
state and exit code. CPU efficiency is approximately TotalCPU divided by
(Elapsed × AllocCPUS), with durations converted to seconds and without double
counting parent/step records. Choose subsequent concurrency and memory from
these measurements and the account limits. A `--time` request is not a runtime
estimate. Do not restart active Cell Ranger jobs simply to change CPU requests.

After annotation, 04 and 06 can run concurrently. After their prerequisite
checkpoints, independent DEG, trajectory, communication and SCENIC branches can
overlap. Audit each script's real input/output paths and shared writes before
submitting: the old SLURM catalog is not an executable scheduler. This change
parallelizes preprocessing; it does not certify every later script for
concurrent execution.
