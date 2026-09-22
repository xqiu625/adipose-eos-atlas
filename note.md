# Doublet-aware rerun note

## Instructions for Claude Code on the UCR HPCC

Treat this section as the execution contract. Claude Code must read this entire
file before changing or submitting anything. It must work from the repository
checkout on the HPCC and must not infer success from `sbatch` accepting a job.
Success means that every required job finished with exit code 0, its expected
outputs exist and are nonempty, and the validation gates below passed.

### Non-negotiable safety and provenance rules

1. This requested rerun starts with the raw-data/configuration preflight and
   includes fresh Cell Ranger outputs in the new versioned rerun root. Reusing
   existing Cell Ranger matrices requires explicit investigator approval; record
   the decision, evidence, and exact reused paths in the audit log.
2. Do not overwrite, rename, or delete the existing analysis directory. The old
   results are the comparison baseline. Do not use `rm`, `git reset`, or a broad
   recursive copy/move command.
3. Do not mix old and new RDS, CSV, H5AD, Loom, CellChat, SCENIC, or figure
   artifacts. Create a new versioned analysis root for this rerun.
4. Do not submit `scripts/SLURM_COMMANDS.sh` as one shell script. It is a command
   catalog; its blocks do not provide all required cross-block dependencies.
5. Submit downstream jobs only with `afterok` dependencies or after explicitly
   verifying the prerequisite job and its outputs.
6. Never silently skip a failed sample or failed method. Stop the dependent
   pipeline and document the failure.
7. Do not change a threshold because it gives a preferred biological result.
   Record all threshold/rate choices before examining downstream conclusions.
8. Do not edit manuscript numbers until the entire comparison report has been
   generated and reviewed by the investigators.

### 1. Establish the run and preserve state

From the HPCC login node, locate the checkout and create a dedicated branch.
Use the actual checkout path if it differs from the example.

```bash
cd /bigdata/nairlab/shared/XQiu/EOS_project/adipose-eos-atlas
git status --short
git rev-parse HEAD
git switch -c doublet-aware-rerun
```

If the checkout is still named `adipose-eos-subtype-atlas`, use that directory
and record the discrepancy. Do not discard a dirty worktree. First save
`git status --short`, `git diff --stat`, and the commit hash in the run log.

Create a unique run identifier and directories. Do not call the variable
`HOME`.

```bash
RUN_ID="doublet_$(date +%Y%m%d_%H%M%S)"
RERUN_ROOT="/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/reruns/${RUN_ID}"
RUN_LOGS="/bigdata/godziklab/shared/Xinru/302025/06_LOGS/${RUN_ID}"
mkdir -p "$RERUN_ROOT" "$RUN_LOGS"
printf '%s\n' "$RUN_ID" > "$RERUN_ROOT/RUN_ID.txt"
git rev-parse HEAD > "$RERUN_ROOT/GIT_COMMIT.txt"
git status --short > "$RERUN_ROOT/GIT_STATUS_AT_START.txt"
```

Before running, Claude Code must make the selected pipeline scripts accept a
versioned base directory, preferably:

```text
EOS_ATLAS_BASE=/.../reruns/<RUN_ID>
EOS_DATA_ROOT=/.../dataset/processed_data
```

The current downstream scripts contain hard-coded `BASE` paths. Therefore,
merely exporting these variables is not sufficient until each script being run
has been updated to read them, for example in R with
`Sys.getenv("EOS_ATLAS_BASE", unset=<current default>)` and in Python with
`os.environ.get("EOS_ATLAS_BASE", <current default>)`. Make this configuration
change on the rerun branch, parse/compile every changed file, and record the
diff. Do not perform a blind replacement of unrelated external-reference paths.

### 2. Preflight gate

Claude Code must run and save the output of all of these checks before the first
job is submitted:

```bash
DATA_ROOT="/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data"
for gem in RR0104 RR0508 RR0912 RR1316; do
  test -d "$DATA_ROOT/$gem/outs/per_sample_outs" || exit 1
done

source activate /bigdata/godziklab/shared/Xinru/cellrank_env
python -c 'import scanpy, scrublet, numpy, pandas, scipy; print("Python imports OK")'
module load R/4.5.0
Rscript -e 'stopifnot(requireNamespace("Seurat"), requireNamespace("celda"), requireNamespace("Matrix")); cat("R imports OK\n")'

python -m py_compile scripts/00a_doublet_detection_per_sample.py
Rscript -e 'invisible(parse(file="scripts/00_decontx_per_sample.R")); invisible(parse(file="scripts/01_qc_doublet_filter.R"))'
bash -n scripts/SLURM_COMMANDS.sh
```

Also verify, rather than assume, all 16 directories
`<GEM>/outs/per_sample_outs/<sample_id>/count/sample_filtered_feature_bc_matrix`
and corresponding `sample_raw_feature_bc_matrix` directories. Compare the
sample/GEM/probe-barcode mapping against the four Cell Ranger `multi` config
files and the web summaries. Save a 16-row verification table in
`$RERUN_ROOT/preflight_sample_map.tsv`.

Record available disk space and stop if the new root cannot safely hold a full
second analysis. Record the Python and R session/package versions in the run
directory.

### 3. Configure the primary doublet run

Use `expected_doublet_rate=0.06` for the prespecified primary sensitivity run,
but do not call it experimentally established. The script must point to the
original Cell Ranger matrices and the new output root:

```bash
export EOS_ATLAS_BASE="$RERUN_ROOT"
export EOS_DATA_ROOT="$DATA_ROOT"
SCRIPTS="$(pwd)/scripts"
CONDA_ENV="/bigdata/godziklab/shared/Xinru/cellrank_env"
```

Update the two R preprocessing scripts to honor those environment variables
before submitting them. Confirm from their logs that they print the new base
path. If a log prints the legacy primary-analysis path, cancel dependent jobs
and fix the configuration.

### 4. Submit the preprocessing chain with SLURM dependencies

Use parsable job IDs and `afterok`. Adjust partition/time only if required by
current HPCC policy; record any change.

```bash
J00A=$(sbatch --parsable --job-name=00a_doublets --partition=batch \
  --mem=128G --cpus-per-task=8 --time=04:00:00 \
  --output="$RUN_LOGS/00a_doublets_%j.out" \
  --error="$RUN_LOGS/00a_doublets_%j.err" \
  --export=ALL,EOS_ATLAS_BASE="$RERUN_ROOT",EOS_DATA_ROOT="$DATA_ROOT" \
  --wrap="source activate $CONDA_ENV && python $SCRIPTS/00a_doublet_detection_per_sample.py --data-root $DATA_ROOT --output-dir $RERUN_ROOT/doublets --expected-rate 0.06 --seed 42")

J00=$(sbatch --parsable --dependency=afterok:$J00A \
  --job-name=00_decontx --partition=batch \
  --mem=128G --cpus-per-task=8 --time=04:00:00 \
  --output="$RUN_LOGS/00_decontx_%j.out" \
  --error="$RUN_LOGS/00_decontx_%j.err" \
  --export=ALL,EOS_ATLAS_BASE="$RERUN_ROOT",EOS_DATA_ROOT="$DATA_ROOT" \
  --wrap="module load R/4.5.0; Rscript $SCRIPTS/00_decontx_per_sample.R")

J01=$(sbatch --parsable --dependency=afterok:$J00 \
  --job-name=01_qc --partition=batch \
  --mem=128G --cpus-per-task=8 --time=04:00:00 \
  --output="$RUN_LOGS/01_qc_%j.out" \
  --error="$RUN_LOGS/01_qc_%j.err" \
  --export=ALL,EOS_ATLAS_BASE="$RERUN_ROOT",EOS_DATA_ROOT="$DATA_ROOT" \
  --wrap="module load R/4.5.0; Rscript $SCRIPTS/01_qc_doublet_filter.R")

printf '00a\t%s\n00\t%s\n01\t%s\n' "$J00A" "$J00" "$J01" \
  > "$RERUN_ROOT/preprocessing_job_ids.tsv"
```

Monitor with `squeue` and then `sacct`. Do not submit the all-cell analysis until
all three jobs show `COMPLETED` with exit code `0:0`.

```bash
squeue -j "$J00A,$J00,$J01"
sacct -j "$J00A,$J00,$J01" --format=JobID,JobName,State,ExitCode,Elapsed,MaxRSS
```

### 5. Preprocessing validation gate

Claude Code must stop before downstream analysis unless all conditions below
are true:

- `doublets/scrublet_summary.csv` has exactly 16 data rows and sample IDs 1-16.
- All 16 `scrublet_sample_<ID>.csv` files exist and are nonempty.
- No score or call is missing or non-finite.
- The sum of singlet and predicted-doublet calls equals the Cell Ranger filtered
  cell count for every sample.
- All 16 `decontx/decontx_sample_<ID>.rds` files exist and are nonempty.
- `decontx/decontx_all_scores.csv` contains all 16 samples.
- `seurat_qc_filtered.rds` exists and is nonempty.
- Logs contain no `ERROR`, traceback, killed process, missing file, or legacy
  output-root warning.

Generate `$RERUN_ROOT/preprocessing_validation.tsv` with one row per sample and
the following columns: sample, GEM, condition, Cell Ranger cells, predicted
doublets, doublet percentage, Scrublet threshold, retained singlets, DecontX
cells, cells above contamination 0.2, and final QC cells. Generate score
histograms per sample and a combined plot stratified by GEM and condition.

Flag for investigator review rather than automatically rejecting a sample if:

- the predicted-doublet percentage is below 0.5% or above 15%;
- one condition or GEM has a conspicuously different removal rate;
- the automatic threshold lies in a dense unimodal region; or
- eosinophil/macrophage marker hybrids are not enriched among predicted calls.

### 6. Submit the core object rebuild

After the preprocessing gate passes, submit this strict chain:

```text
01 QC
  -> 02 all-cell clustering
  -> 03 annotation
  -> 04 eosinophil subclustering
     (includes manuscript cluster naming and C6 removal)

03 annotation
  -> 06 macrophage subtyping

04 + 06
  -> 09 combined eosinophil/macrophage object
```

Use `afterok` dependencies. `04_eos_subcluster.R` now performs manuscript-label
assignment and C6 removal and must finish before any script consumes
`eos_subclustered.rds`. `04b_rename_clusters.R` is retired and intentionally
aborts; do not run it. Do not use the older manuscript object as an automatic
label source without documenting exactly where it enters. In particular,
review the previous-object dependency in `03_annotate_celltypes.R` and the
barcode-based correspondence and Hungarian assignment in
`04_eos_subcluster.R`; confirm that all new barcodes are accounted for, no cell
is silently dropped, and ambiguous assignments are recorded for investigator
review.

At the end of this phase validate:

- seven intended broad immune cell types are present;
- every sample and condition remains represented;
- C0-C5 exist after C6 exclusion;
- cluster labels contain no `NA`;
- `mac_subtyped.rds` has the three intended macrophage states; and
- cell totals reconcile across every parent/subset object.

### 7. Rerun downstream analyses in dependency waves

Use the commands in `scripts/SLURM_COMMANDS.sh` as templates, but add the new
base-directory export and explicit `afterok` dependencies. Run these waves in
order and validate each wave before continuing:

1. **Core eosinophil results after 04:** `05`, `05b`, `07c`, `08`, `17`, `18`,
   `19`, `20`, `21`, `25`, `27`, `31`, `43`.
2. **Macrophage/communication after 04 + 06 + 09:** `10`, `10b`, `10c`, `33`,
   `34`, `44`.
3. **External integration/trajectory after 04:** `11`, `11b`, `12`, `13`,
   `13b`, `23`, `23b`, `24`, `24b`, `24c`, `28`, `29`, `30`, `32`, `39`.
4. **SCENIC after 04 and 06 where required:** `14`, `14b`, `14c`, `14d`, `15`,
   `15c`, `16`, `16b`, `16c`, `16d`, `16e`, `22`, `38`, `41`-`42h`.
5. **Integration/robustness:** `37`, `45`-`61`, followed by `62a`-`75` only if
   that revision pipeline is declared part of the final analysis.

Before running `20_compositional_analysis.R`, verify the two implemented safety
conditions: the sex contrast is WT female versus WT male only, and all missing
sample-by-cluster combinations are explicitly completed with count zero. Record
the relevant log output and input-table checks in the stage audit checkpoint.

Do not assume that every numbered script is canonical merely because it exists.
For each manuscript panel, use `FIGURES.md`, record the exact producing script
and input object, and mark exploratory alternatives as such.

### 8. Required sensitivity runs

After the prespecified 0.06 run succeeds, rerun only the preprocessing and core
conclusion path at two additional plausible expected rates selected from Cell
Ranger loading information. If no better experimental estimate is available,
use 0.04 and 0.08 and label them sensitivity settings, not corrected truth.
Each rate must have its own output root. Compare, at minimum, C0-C5 recovery,
cluster proportions, primary pseudobulk contrasts, eosinophil/macrophage hybrid
markers, and central CCL CellChat conclusions.

### 9. Failure handling

If a job fails, Claude Code must:

1. stop or cancel jobs that depend on invalid outputs;
2. capture `sacct`, stdout, stderr, command, commit, and environment information;
3. identify whether the failure is code, data, memory/time, environment, or
   scheduler related;
4. make the smallest justified correction on the rerun branch;
5. rerun syntax checks and the failed job; and
6. record the correction in `$RERUN_ROOT/RERUN_LOG.md`.

Do not convert errors into warnings, skip samples, or reuse an output from a
previous run to make the pipeline continue.

### 10. Final deliverables Claude Code must produce

Do not report completion until the following exist:

- `RERUN_LOG.md`: chronological commands, job IDs, failures, corrections, and
  validation decisions.
- `software_versions.txt`: module, R session, Python, Conda, and package versions.
- `preflight_sample_map.tsv` and `preprocessing_validation.tsv`.
- `artifact_manifest.tsv`: relative path, producing script, producing job ID,
  input object(s), size, modification time, and SHA-256 for every final object,
  table, and manuscript figure.
- `old_vs_doublet_aware_summary.md`: old and new cell counts, clusters,
  composition, DE, CellChat, trajectory, and SCENIC comparisons.
- `manuscript_changes_required.md`: every number, method statement, figure, and
  conclusion that must change, plus conclusions that were reproduced.
- A final `sacct` table proving all reported jobs completed successfully.

The final message to the investigators must distinguish: completed and
validated; completed with a scientific warning; failed; not run; and blocked by
missing data or investigator choice. It must never describe queued/running jobs
as completed.

## Why this rerun is required

The manuscript states that Scrublet doublet removal was part of preprocessing,
but the provenance files needed to show that the reported objects were produced
with that step are not present. The previous repository implementation also ran
Scrublet after DecontX, although Scrublet is intended to receive an unnormalized
count matrix. The pipeline has therefore been changed to run Scrublet on each
sample's original Cell Ranger filtered matrix before DecontX.

This rerun is a sensitivity analysis unless original logs establish that the
published analysis used the same procedure. Do not describe the new outputs as
the original analysis. Retain both the pre-rerun results and the new results,
and report whether conclusions change.

The 10x Flex probe barcodes help resolve multiplets containing cells from
different samples in the same GEM. They do not eliminate all same-sample or
homotypic doublets. Scrublet calls are predictions and should be accompanied by
marker-based inspection and per-sample summaries.

## Configuration to verify before running

1. Confirm the 16 sample-to-GEM assignments in
   `scripts/00a_doublet_detection_per_sample.py` and
   `scripts/00_decontx_per_sample.R` against the Cell Ranger `multi` configs.
2. Confirm that `data_root` points to the original Cell Ranger 8.0.1 outputs,
   not a previously corrected matrix.
3. Review each Cell Ranger web summary and the `multiplexing_analysis` outputs.
   Archive the recovered-cell counts, barcode-overlap information, and any
   multiplet-related metrics for every GEM reaction.
4. Justify `--expected-rate 0.06` using the recovered-cell/loading information.
   Also run at least one lower and one higher plausible rate as a sensitivity
   check; do not choose a rate based on which result best preserves clusters.
5. Use a frozen Python environment and record versions of Python, Scanpy,
   Scrublet, NumPy, pandas, and SciPy. The repository still needs a lock file.

## Required rerun order

Start from raw-data/configuration validation, then generate fresh Cell Ranger
outputs in the versioned rerun root. Do not skip any step below. If the
investigator explicitly approves reuse of existing Cell Ranger matrices, record
that scope deviation and its justification before continuing.

1. `scripts/00a_doublet_detection_per_sample.py`
   - Input: each sample's uncorrected
     `sample_filtered_feature_bc_matrix`.
   - Output: `doublets/scrublet_sample_<ID>.csv`,
     `doublets/scrublet_summary.csv`, and recorded run parameters.
   - Stop and investigate if any sample is absent, has an implausible call rate,
     or lacks a clear score distribution.
2. `scripts/00_decontx_per_sample.R`
   - Requires calls for all 16 samples.
   - Removes predicted doublets before fitting DecontX.
   - Uses the unfiltered Cell Ranger matrix as ambient background.
   - Applies the manuscript-consistent contamination threshold of 0.2.
3. `scripts/01_qc_doublet_filter.R`
   - Applies feature/mitochondrial/contamination QC and merges the 16 samples.
   - It no longer reruns Scrublet after DecontX.
4. `scripts/02_all_cells_umap.R`
5. `scripts/03_annotate_celltypes.R`
6. `scripts/04_eos_subcluster.R` (includes manuscript-label assignment and C6
   removal; do not run the retired `scripts/04b_rename_clusters.R`)
7. `scripts/06_mac_subtype.R`

Use `scripts/SLURM_COMMANDS.sh` for the updated initial submission order.

## Downstream code that must be rerun

Because doublet removal changes the cells and count aggregates, all results
derived from the all-cell, eosinophil, or macrophage objects must be regenerated.
At minimum rerun:

- QC and cluster validation: `17`, `18`, `19`, `21`, `37`, `47`-`53`, `59`,
  and `59b`.
- Pseudobulk differential expression and enrichment: `05`, `05b`, `07c`, `08`,
  `25`, and `27`.
- Eosinophil composition: `20` and `43`. Before interpreting these outputs,
  separately fix the known WT-only sex-contrast and missing-zero-combination
  issues in `20_compositional_analysis.R`.
- Eosinophil/macrophage objects and CellChat: `06`, `09`, `10`, `10b`, `10c`,
  `33`, `34`, and `44`.
- Trajectory and external integration: `11`, `11b`, `12`, `13`, `13b`, `23`,
  `23b`, `24`, `24b`, `24c`, `28`, `29`, `30`, `31`, `32`, `39`, and `41`-`42h`.
- SCENIC: `14`, `14b`, `14c`, `14d`, `15`, `15c`, `16`-`16e`, `22`, and `38`.
- Integration and robustness analyses: `45`-`61` and, if they are retained as
  part of the final analysis, the revision-pipeline `62a`-`75` scripts.
- Every manuscript figure and numerical statement based on cell counts,
  proportions, markers, DEGs, pathways, pseudotime, CellChat, or regulons.

Do not reuse old downstream RDS, H5AD, Loom, CSV, or CellChat objects. Write the
rerun to a new versioned output directory or archive the old directory before
starting, so that old and new artifacts cannot be mixed.

## Validation required after rerunning

Create a table for every sample containing:

- Cell Ranger filtered-cell count.
- Scrublet score threshold and expected rate.
- Predicted doublet count and percentage.
- Counts after DecontX contamination filtering and standard QC.
- Eosinophil, macrophage, and other-cell counts after annotation.

Then compare old versus doublet-aware results for:

- all-cell and eosinophil UMAPs;
- C0-C5 size, markers, and assignment stability;
- per-sample cluster proportions and sccomp conclusions;
- pseudobulk DEG direction and significance;
- eosinophil-macrophage hybrid-marker cells;
- CellChat CCL interactions;
- trajectory ordering; and
- SCENIC regulon rankings.

Inspect predicted doublets by sample, GEM reaction, condition, original cluster,
`nCount_RNA`, `nFeature_RNA`, and canonical eosinophil/macrophage markers. A
large condition-specific removal rate must be investigated because it can create
or erase apparent compositional effects.

## Manuscript rule

If this rerun is retrospective, write that explicitly. Report the original
preprocessing accurately, describe doublet removal as a sensitivity analysis,
and state whether each principal conclusion was reproduced. Only state that
doublets were removed in the primary pipeline if archived provenance proves that
the primary reported objects actually came from that pipeline.
