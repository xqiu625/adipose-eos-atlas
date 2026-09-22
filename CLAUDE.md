# Claude Code: HPCC rerun mandate

Read [`note.md`](note.md) completely before editing code, submitting a job, or
reusing any artifact. It is the execution contract for the audited HPCC rerun.

Start at the raw-data preflight and provenance capture. Because the requested
scope is a from-the-beginning rerun, validate the FASTQs, checksums, Cell Ranger
multi configurations, reference, probe set, sample/probe-barcode map, and the
existing Cell Ranger web summaries. Do not overwrite the historical Cell Ranger
or analysis outputs. Submit a fresh Cell Ranger run only after the raw-input
preflight passes and the new versioned output root has been recorded. If the
investigator elects to reuse the already validated Cell Ranger matrices, record
that explicit decision and begin computation at `00a`; do not silently treat
matrix reuse as a raw-data rerun.

After every stage, stop and write an audit checkpoint before submitting any
dependent stage. Each checkpoint must contain:

- stage name and producing script;
- exact command and SLURM job ID;
- Git commit and dirty-worktree status;
- input paths and SHA-256 hashes where practical;
- expected versus observed output files, sizes, and row/cell counts;
- `sacct` state, exit code, elapsed time, and maximum RSS;
- scan of stdout/stderr for errors, warnings, tracebacks, killed jobs, and
  legacy output paths;
- scientific validation metrics required by `note.md`;
- one status: `PASS`, `PASS_WITH_WARNING`, `FAIL`, `BLOCKED`, or `NOT_RUN`;
- the reason and next authorized action.

Append the checkpoint to `$RERUN_ROOT/RERUN_LOG.md` and write a machine-readable
row to `$RERUN_ROOT/audit_checkpoints.tsv`. Print the same concise checkpoint in
the Claude Code response so the investigator receives feedback at every step.
Only `PASS` permits automatic submission of dependents. `PASS_WITH_WARNING`
requires explicit investigator approval. `FAIL` or `BLOCKED` stops the affected
branch.

The HPCC `$RERUN_ROOT` must mirror the relative directory layout of the local
destination `/Users/xinruqiu/Desktop/302025/adipose-eos-atlas-results-v2` (for
example, `04_eos_subcluster/results`, `04_eos_subcluster/figures`, and the stage
logs). After each `PASS`, report the exact remote relative paths ready for
transfer and add them to `$RERUN_ROOT/transfer_ready.tsv` with stage, path, byte
size, file count, and checkpoint status. Do not claim that HPCC pushed files to
the Mac: the transfer is initiated locally with the v2 directory's
`pull_from_hpcc.sh` over SSH/SCP. Never transfer failed or unvalidated partial
outputs as final results.

Important current pipeline facts:

- `scripts/04_eos_subcluster.R` performs manuscript cluster assignment and C6
  removal itself.
- `scripts/04b_rename_clusters.R` is retired and intentionally aborts. Never
  run it and never make downstream work depend on it.
- The WT-only sex contrast and missing sample-by-cluster zero combinations are
  already implemented in `scripts/20_compositional_analysis.R`; validate them
  during the rerun rather than reapplying an obsolete fix.
- Most downstream scripts still contain a hard-coded legacy analysis root.
  They must be changed and syntax-checked to honor the versioned rerun root
  before execution. A log containing the legacy output root is a failure.
- Do not claim completion from successful `sbatch` submission. Completion
  requires `COMPLETED`, exit code `0:0`, nonempty expected outputs, and a passed
  scientific validation gate.

Audit update for 2026-09-22: the stale-variable `NameError` in the final summary
of `scripts/42e_regulon_unique_analysis.py` has been corrected and the script
passes Python compilation. Before the later SCENIC wave, add explicit
empty-result guards before sorting `fc_data`/`diff_data` in
`42f_regulon_venn_and_comparison.py` and `42g_regulon_dotplots.py`, respectively,
then compile and smoke-test those paths with small synthetic tables. Treat any
partially written outputs from a failed run as invalid.
