#!/usr/bin/env python3
"""Audited array workers and collectors (standard library only).

This controls execution, not scientific QC. A collector's success does not
authorize the next stage: a separate reviewed gate is required by prepare.
"""
import argparse
import csv
import hashlib
import json
import math
import os
import re
from pathlib import Path
import shutil
import subprocess
import sys

STAGES = {"scrublet": "00a_doublet_detection_per_sample.py",
          "decontx": "00_decontx_per_sample.R"}
# Must match BACKGROUND_POLICY_ALLOWED in 00_decontx_per_sample.R exactly
# (Checkpoint 23/26). Frozen into run.json at prepare time and injected into
# the worker subprocess environment, overriding any ambient shell value, so
# the audit trail always records the actual reviewed policy that ran.
BACKGROUND_POLICIES = ("historical", "exclude_all_original_called_cells")
POOLS = ["RR0104"] * 4 + ["RR0508"] * 4 + ["RR0912"] * 4 + ["RR1316"] * 4
MATRIX_FILES = {"matrix": "matrix.mtx.gz", "barcodes": "barcodes.tsv.gz",
                "features": "features.tsv.gz"}


def digest(path):
    h = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(8 * 1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def read_json(path):
    return json.loads(Path(path).read_text())


def write_json(path, value):
    # Exclusive creation prevents reusing a successful task or aggregation.
    with Path(path).open("x") as handle:
        json.dump(value, handle, indent=2)
        handle.write("\n")


def table(path):
    with Path(path).open(newline="") as handle:
        return list(csv.DictReader(handle))


def manifest(path):
    with Path(path).open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    if len(rows) != 16 or sorted(int(r["sample_id"]) for r in rows) != list(range(1, 17)):
        raise ValueError("Manifest must contain each sample 1–16 exactly once")
    rows.sort(key=lambda r: int(r["sample_id"]))
    for row in rows:
        sid = int(row["sample_id"])
        if row["pool"] != POOLS[sid - 1]:
            raise ValueError(f"Wrong pool for sample {sid}")
        if row.get("qc_status") != "PASS" or row.get("provenance_status") != "PASS":
            raise ValueError(f"Sample {sid} has not passed input QC/provenance")
        for key in ("filtered_dir", "raw_dir"):
            directory = Path(row[key])
            if not directory.is_absolute() or not directory.is_dir():
                raise ValueError(f"Missing absolute matrix directory: {directory}")
            for kind, name in MATRIX_FILES.items():
                if not (directory / name).is_file() or (directory / name).stat().st_size == 0:
                    raise ValueError(f"Missing/empty input: {directory / name}")
                hash_column = f"{key.removesuffix('_dir')}_{kind}_sha256"
                if not re.fullmatch(r"[0-9a-f]{64}", row.get(hash_column, "")):
                    raise ValueError(f"Missing/invalid preflight checksum: {hash_column}, sample {sid}")
    return rows


def state_dir(base, stage):
    return base / "parallel" / stage


def configuration(base, stage):
    state = state_dir(base, stage)
    config = read_json(state / "run.json")
    for name, expected in config["snapshot_hashes"].items():
        if digest(state / name) != expected:
            raise ValueError(f"Frozen run snapshot changed: {name}")
    if str(base.resolve()) != config["base"]:
        raise ValueError("Run directory moved; refusing ambiguous provenance")
    return state, config


def validate_receipt(base, stage, sid):
    state, config = configuration(base, stage)
    receipt = read_json(state / f"sample_{sid}.json")
    if (receipt["sample_id"] != sid or receipt["stage"] != stage or
            receipt["run_sha256"] != digest(state / "run.json")):
        raise ValueError(f"Receipt belongs to a different run/sample: {sid}")
    for relative, expected in receipt["outputs"].items():
        path = base / relative
        if not path.is_file() or path.stat().st_size == 0 or digest(path) != expected:
            raise ValueError(f"Missing or changed task output: {path}")
    return receipt


def verify_collection(base, stage):
    state, _ = configuration(base, stage)
    report = read_json(state / "collection.json")
    if report["samples"] != list(range(1, 17)) or report["run_sha256"] != digest(state / "run.json"):
        raise ValueError("Invalid collection identity")
    for relative, expected in report["outputs"].items():
        if digest(base / relative) != expected:
            raise ValueError(f"Changed aggregate: {relative}")
    for sid in range(1, 17):
        receipt = validate_receipt(base, stage, sid)
        if digest(state / f"sample_{sid}.json") != report["receipts"][str(sid)]:
            raise ValueError(f"Changed receipt: {sid}")
    return report


def prepare(args):
    base = args.base.resolve()
    root = Path(os.environ["RERUN_ROOT"]).resolve()
    if base != root and root not in base.parents:
        raise ValueError("Analysis base must be inside RERUN_ROOT")
    if root == Path("/"):
        raise ValueError("RERUN_ROOT cannot be /")
    if not 0 < args.expected_rate < 1:
        raise ValueError("Expected rate must be between zero and one")
    if args.stage == "decontx" and args.background_policy not in BACKGROUND_POLICIES:
        raise ValueError(
            f"--background-policy is required for decontx and must be one of "
            f"{BACKGROUND_POLICIES}, got {args.background_policy!r}")
    rows = manifest(args.manifest)
    gate = read_json(args.gate)
    required_stage = "cellranger_qc_and_00a_audit" if args.stage == "scrublet" else "scrublet_qc_and_decontx_audit"
    source = Path(__file__).resolve().parent
    code_hash = digest(source / STAGES[args.stage])
    if (gate.get("status") != "PASS" or gate.get("stage") != required_stage or
            gate.get("manifest_sha256") != digest(args.manifest) or
            gate.get("script_sha256") != code_hash or not gate.get("checkpoint") or
            not gate.get("reviewed_by")):
        raise ValueError("Gate must be reviewed PASS and bind this manifest and analysis script")
    if args.stage == "decontx":
        verify_collection(base, "scrublet")
        previous = state_dir(base, "scrublet")
        if read_json(previous / "run.json")["manifest_sha256"] != digest(args.manifest):
            raise ValueError("DecontX and Scrublet input manifests differ")
        if gate.get("collection_sha256") != digest(previous / "collection.json"):
            raise ValueError("DecontX gate must bind the reviewed Scrublet collection")
    for row in rows:
        for key in ("filtered_dir", "raw_dir"):
            for name in MATRIX_FILES.values():
                resolved = (Path(row[key]) / name).resolve()
                if root not in resolved.parents:
                    raise ValueError(f"Input is not a fresh output inside RERUN_ROOT: {resolved}")
                if row["pool"] == "RR1316" and "RR1316_L007L008" not in resolved.parts:
                    raise ValueError("Samples 13–16 must use validated RR1316_L007L008 outputs")
    output = base / ("doublets" if args.stage == "scrublet" else "decontx")
    if output.exists():
        raise FileExistsError(f"Use a fresh analysis base; output already exists: {output}")
    state = state_dir(base, args.stage)
    state.mkdir(parents=True, exist_ok=False)
    output.mkdir()
    (state / "logs").mkdir()
    snapshots = {"input_manifest.tsv": args.manifest, "gate.json": args.gate,
                 STAGES[args.stage]: source / STAGES[args.stage],
                 "parallel_preprocessing.py": Path(__file__),
                 "run_preprocessing_task.sh": source / "run_preprocessing_task.sh"}
    for name, path in snapshots.items():
        shutil.copyfile(path, state / name)
    revision = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=source, text=True).strip()
    dirty = subprocess.check_output(["git", "status", "--short"], cwd=source, text=True)
    write_json(state / "run.json", {
        "stage": args.stage, "base": str(base), "git_commit": revision,
        "git_status": dirty, "manifest_sha256": digest(args.manifest),
        "snapshot_hashes": {name: digest(state / name) for name in snapshots},
        "resolved_inputs": {str(Path(row[key]) / name): str((Path(row[key]) / name).resolve())
                            for row in rows for key in ("filtered_dir", "raw_dir")
                            for name in MATRIX_FILES.values()},
        "expected_rate": args.expected_rate, "seed": args.seed,
        "background_policy": args.background_policy if args.stage == "decontx" else None,
    })
    print(state)


def worker(args):
    base, sid, stage = args.base.resolve(), args.sample_id, args.stage
    state, config = configuration(base, stage)
    # An exclusive persistent claim prevents duplicate submissions and partial reuse.
    (state / f"sample_{sid}.lock").mkdir()
    rows = manifest(state / "input_manifest.tsv")
    row = rows[sid - 1]
    inputs = {}
    for key in (("filtered_dir",) if stage == "scrublet" else ("filtered_dir", "raw_dir")):
        for kind, name in MATRIX_FILES.items():
            path = Path(row[key]) / name
            if str(path.resolve()) != config["resolved_inputs"][str(path)]:
                raise ValueError(f"Input symlink target changed since preflight: {path}")
            actual = digest(path)
            if actual != row[f"{key.removesuffix('_dir')}_{kind}_sha256"]:
                raise ValueError(f"Input checksum differs from reviewed manifest: {path}")
            inputs[str(path.resolve())] = actual
    if stage == "decontx":
        verify_collection(base, "scrublet")
        inputs[str(base / "doublets" / f"scrublet_sample_{sid}.csv")] = digest(
            base / "doublets" / f"scrublet_sample_{sid}.csv")
    command = ([sys.executable] if stage == "scrublet" else ["Rscript"]) + [
        str(state / STAGES[stage]), "--sample-id", str(sid),
        "--input-manifest", str(state / "input_manifest.tsv")]
    env = os.environ.copy()
    env["EOS_ATLAS_BASE"] = str(base)
    if stage == "scrublet":
        command += ["--output-dir", str(base / "doublets"), "--expected-rate",
                    str(config["expected_rate"]), "--seed", str(config["seed"])]
        names = [f"doublets/scrublet_sample_{sid}.csv", f"doublets/scrublet_summary_sample_{sid}.csv",
                 f"doublets/scrublet_run_parameters_sample_{sid}.json"]
    else:
        if config.get("background_policy") not in BACKGROUND_POLICIES:
            raise ValueError(f"Frozen run.json has no valid background_policy: {config.get('background_policy')!r}")
        # Overrides any ambient shell value: the policy that actually runs is
        # always exactly what was reviewed and frozen at prepare time.
        env["DECONTX_BACKGROUND_POLICY"] = config["background_policy"]
        names = [f"decontx/decontx_sample_{sid}.rds", f"decontx/decontx_scores_sample_{sid}.csv"]
    if any((base / name).exists() for name in names):
        raise FileExistsError("Task outputs already exist; use a new attempt")
    print(json.dumps({"command": command, "inputs": inputs}), flush=True)
    subprocess.run(command, env=env, check=True)
    if any(not (base / name).is_file() or (base / name).stat().st_size == 0 for name in names):
        raise ValueError("Task returned without all expected nonempty outputs")
    write_json(state / f"sample_{sid}.json", {
        "stage": stage, "sample_id": sid, "run_sha256": digest(state / "run.json"),
        "inputs": inputs, "command": command, "job_id": os.environ.get("SLURM_JOB_ID"),
        "array_job_id": os.environ.get("SLURM_ARRAY_JOB_ID"),
        "outputs": {name: digest(base / name) for name in names},
    })


def collect(args):
    base, stage = args.base.resolve(), args.stage
    state, config = configuration(base, stage)
    (state / "collect.lock").mkdir()
    combined = []
    packages = []
    receipts = {}
    for sid in range(1, 17):
        validate_receipt(base, stage, sid)
        receipts[str(sid)] = digest(state / f"sample_{sid}.json")
        name = (f"doublets/scrublet_sample_{sid}.csv" if stage == "scrublet"
                else f"decontx/decontx_scores_sample_{sid}.csv")
        rows = table(base / name)
        if not rows or any(int(r["sample_id"]) != sid for r in rows):
            raise ValueError(f"Empty/wrong sample table: {sid}")
        if len({r["barcode"] for r in rows}) != len(rows):
            raise ValueError(f"Duplicate barcode within sample {sid}")
        if any(not math.isfinite(float(r["scrublet_score"])) for r in rows):
            raise ValueError(f"Invalid Scrublet score: {sid}")
        if stage == "scrublet":
            if any(r["scrublet_doublet"].lower() not in ("true", "false") for r in rows):
                raise ValueError(f"Invalid doublet calls: {sid}")
            summaries = table(base / f"doublets/scrublet_summary_sample_{sid}.csv")
            if len(summaries) != 1 or int(summaries[0]["sample_id"]) != sid:
                raise ValueError(f"Invalid summary identity: {sid}")
            summary = summaries[0]
            n_doublets = sum(r["scrublet_doublet"].lower() == "true" for r in rows)
            if int(summary["n_cells"]) != len(rows) or int(summary["n_predicted_doublets"]) != n_doublets:
                raise ValueError(f"Summary/calls disagree: {sid}")
            if not math.isfinite(float(summary["threshold"])):
                raise ValueError(f"Missing/nonfinite threshold: {sid}")
            record = read_json(base / f"doublets/scrublet_run_parameters_sample_{sid}.json")
            if (record["parameters"]["seed"] != config["seed"] or
                    record["parameters"]["expected_rate"] != config["expected_rate"]):
                raise ValueError(f"Mixed analysis parameters: {sid}")
            packages.append((record["python"], record["packages"]))
            combined.append(summary)
        else:
            calls = table(base / f"doublets/scrublet_sample_{sid}.csv")
            expected = {r["barcode"] for r in calls if r["scrublet_doublet"].lower() == "false"}
            if {r["barcode"] for r in rows} != expected:
                raise ValueError(f"DecontX cells do not equal predicted singlets: {sid}")
            if any(not math.isfinite(float(r["contamination"])) or
                   not 0 <= float(r["contamination"]) <= 1 for r in rows):
                raise ValueError(f"Invalid contamination values: {sid}")
            combined.extend(rows)
    if packages and any(p != packages[0] for p in packages):
        raise ValueError("Mixed Python/package versions across tasks")
    relative = ("doublets/scrublet_summary.csv" if stage == "scrublet"
                else "decontx/decontx_all_scores.csv")
    with (base / relative).open("x", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(combined[0]))
        writer.writeheader()
        writer.writerows(combined)
    outputs = {relative: digest(base / relative)}
    if stage == "scrublet":
        record_path = base / "doublets/scrublet_run_parameters.json"
        write_json(record_path, {"parameters": config, "python": packages[0][0],
                                 "packages": packages[0][1]})
        outputs[str(record_path.relative_to(base))] = digest(record_path)
    write_json(state / "collection.json", {
        "status": "COLLECTED_PENDING_SCIENTIFIC_QC", "samples": list(range(1, 17)),
        "run_sha256": digest(state / "run.json"), "receipts": receipts, "outputs": outputs,
    })
    print("16 tasks collected. Record sacct and scientific QC checkpoint before the next stage.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("prepare", "worker", "collect", "verify"))
    parser.add_argument("--stage", required=True, choices=STAGES)
    parser.add_argument("--base", required=True, type=Path)
    parser.add_argument("--manifest", type=Path)
    parser.add_argument("--gate", type=Path)
    parser.add_argument("--sample-id", type=int, choices=range(1, 17))
    parser.add_argument("--expected-rate", type=float, default=0.06)
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--background-policy", choices=BACKGROUND_POLICIES)
    args = parser.parse_args()
    if args.action == "prepare" and (not args.manifest or not args.gate):
        parser.error("prepare needs --manifest and --gate")
    if args.action == "prepare" and args.stage == "decontx" and args.background_policy is None:
        parser.error("prepare needs --background-policy for stage=decontx")
    if args.action == "worker" and args.sample_id is None:
        parser.error("worker needs --sample-id")
    if args.action == "verify":
        verify_collection(args.base.resolve(), args.stage)
    else:
        globals()[args.action](args)


if __name__ == "__main__":
    main()
