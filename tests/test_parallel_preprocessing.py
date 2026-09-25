"""Orchestration tests; biological calculations are mocked, not validated here."""
import argparse
import csv
import importlib.util
import contextlib
import io
import json
import os
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest
from unittest.mock import patch

SCRIPTS = Path(__file__).resolve().parents[1] / "scripts"
spec = importlib.util.spec_from_file_location("parallel", SCRIPTS / "parallel_preprocessing.py")
p = importlib.util.module_from_spec(spec)
spec.loader.exec_module(p)


def csv_write(path, rows, delimiter=","):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter=delimiter)
        writer.writeheader()
        writer.writerows(rows)


class ParallelTests(unittest.TestCase):
    def setUp(self):
        silence = contextlib.redirect_stdout(io.StringIO())
        silence.__enter__()
        self.addCleanup(silence.__exit__, None, None, None)
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name).resolve()
        self.base = self.root / "analysis"
        self.env = patch.dict(os.environ, {"RERUN_ROOT": str(self.root)})
        self.env.start()
        self.addCleanup(self.env.stop)
        self.rows = []
        for sid in range(1, 17):
            pool = p.POOLS[sid - 1]
            run = "RR1316_L007L008" if pool == "RR1316" else pool
            row = {"sample_id": sid, "pool": pool, "qc_status": "PASS", "provenance_status": "PASS"}
            for key in ("filtered_dir", "raw_dir"):
                path = self.root / "cellranger" / run / str(sid) / key
                path.mkdir(parents=True)
                for name in ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"):
                    (path / name).write_bytes(f"synthetic fixture {sid} {name}".encode())
                row[key] = str(path)
                for kind, name in p.MATRIX_FILES.items():
                    row[f"{key.removesuffix('_dir')}_{kind}_sha256"] = p.digest(path / name)
            self.rows.append(row)
        self.manifest = self.root / "manifest.tsv"
        csv_write(self.manifest, self.rows, "\t")
        self.captured_envs = []

    def args(self, stage="scrublet", sample_id=None, background_policy="historical"):
        return argparse.Namespace(base=self.base, stage=stage, sample_id=sample_id,
                                  manifest=self.manifest, gate=self.root / f"{stage}_gate.json",
                                  expected_rate=0.06, seed=42,
                                  background_policy=background_policy if stage == "decontx" else None)

    def prepare(self, stage="scrublet", execute=True, background_policy="historical"):
        args = self.args(stage, background_policy=background_policy)
        gate = {"stage": "cellranger_qc_and_00a_audit" if stage == "scrublet" else "scrublet_qc_and_decontx_audit",
                "status": "PASS", "manifest_sha256": p.digest(self.manifest),
                "script_sha256": p.digest(SCRIPTS / p.STAGES[stage]),
                "checkpoint": "synthetic-test-only", "reviewed_by": "unit-test"}
        if stage == "decontx":
            gate["collection_sha256"] = p.digest(p.state_dir(self.base, "scrublet") / "collection.json")
        args.gate.write_text(json.dumps(gate))
        if execute:
            p.prepare(args)

    def fake_analysis(self, command, env, check):
        self.captured_envs.append(env)
        sid = int(command[command.index("--sample-id") + 1])
        stage = "scrublet" if "--output-dir" in command else "decontx"
        if stage == "scrublet":
            out = self.base / "doublets"
            rows = [{"barcode": "AA", "sample_id": sid, "scrublet_score": 0.01, "scrublet_doublet": False},
                    {"barcode": "BB", "sample_id": sid, "scrublet_score": 0.8, "scrublet_doublet": True}]
            csv_write(out / f"scrublet_sample_{sid}.csv", rows)
            csv_write(out / f"scrublet_summary_sample_{sid}.csv",
                      [{"sample_id": sid, "n_cells": 2, "n_predicted_doublets": 1, "threshold": 0.3}])
            (out / f"scrublet_run_parameters_sample_{sid}.json").write_text(json.dumps({
                "python": "test", "packages": {"scrublet": "test"},
                "parameters": {"seed": 42, "expected_rate": 0.06}}))
        else:
            out = self.base / "decontx"
            (out / f"decontx_sample_{sid}.rds").write_bytes(b"mock RDS, not biological validation")
            csv_write(out / f"decontx_scores_sample_{sid}.csv", [{
                "barcode": "AA", "sample_id": sid, "condition": "test",
                "contamination": 0.1, "scrublet_score": 0.01}])

    def worker(self, sid, stage="scrublet"):
        with patch.object(p.subprocess, "run", side_effect=self.fake_analysis):
            p.worker(self.args(stage, sid))

    def complete_scrublet(self):
        self.prepare()
        for sid in range(1, 17):
            self.worker(sid)
        p.collect(self.args())

    def test_two_stages_collect_all_samples(self):
        self.complete_scrublet()
        self.assertEqual(len(p.table(self.base / "doublets/scrublet_summary.csv")), 16)
        self.prepare("decontx")
        for sid in range(1, 17):
            self.worker(sid, "decontx")
        p.collect(self.args("decontx"))
        report = p.verify_collection(self.base, "decontx")
        self.assertEqual(report["status"], "COLLECTED_PENDING_SCIENTIFIC_QC")
        self.assertEqual(len(p.table(self.base / "decontx/decontx_all_scores.csv")), 16)

    def test_single_sample_never_writes_aggregate_or_other_sample(self):
        self.prepare()
        self.worker(14)
        self.assertFalse((self.base / "doublets/scrublet_summary.csv").exists())
        self.assertFalse((self.base / "doublets/scrublet_sample_13.csv").exists())
        self.assertEqual(p.validate_receipt(self.base, "scrublet", 14)["sample_id"], 14)

    def test_duplicate_worker_rejected_without_overwrite(self):
        self.prepare()
        self.worker(1)
        before = p.digest(self.base / "doublets/scrublet_sample_1.csv")
        with self.assertRaises(FileExistsError):
            self.worker(1)
        self.assertEqual(before, p.digest(self.base / "doublets/scrublet_sample_1.csv"))

    def test_missing_task_blocks_collection(self):
        self.prepare()
        self.worker(1)
        with self.assertRaises(FileNotFoundError):
            p.collect(self.args())
        self.assertFalse((self.base / "doublets/scrublet_summary.csv").exists())

    def test_failed_analysis_has_no_success_receipt(self):
        self.prepare()
        with patch.object(p.subprocess, "run", side_effect=subprocess.CalledProcessError(1, "analysis")):
            with self.assertRaises(subprocess.CalledProcessError):
                p.worker(self.args(sample_id=1))
        self.assertFalse((p.state_dir(self.base, "scrublet") / "sample_1.json").exists())

    def test_changed_input_fails_before_analysis(self):
        self.prepare()
        (Path(self.rows[0]["filtered_dir"]) / "matrix.mtx.gz").write_bytes(b"replaced")
        with patch.object(p.subprocess, "run") as run:
            with self.assertRaisesRegex(ValueError, "checksum differs"):
                p.worker(self.args(sample_id=1))
            run.assert_not_called()

    def test_background_policy_required_for_decontx_prepare(self):
        self.complete_scrublet()
        with self.assertRaisesRegex(ValueError, "background-policy"):
            self.prepare("decontx", background_policy=None)

    def test_background_policy_rejected_if_invalid(self):
        self.complete_scrublet()
        with self.assertRaisesRegex(ValueError, "background-policy"):
            self.prepare("decontx", background_policy="bogus_policy")

    def test_background_policy_frozen_into_run_json(self):
        self.complete_scrublet()
        self.prepare("decontx", background_policy="exclude_all_original_called_cells")
        run_json = p.read_json(p.state_dir(self.base, "decontx") / "run.json")
        self.assertEqual(run_json["background_policy"], "exclude_all_original_called_cells")

    def test_scrublet_run_json_has_no_background_policy(self):
        self.prepare("scrublet")
        run_json = p.read_json(p.state_dir(self.base, "scrublet") / "run.json")
        self.assertIsNone(run_json["background_policy"])

    def test_worker_injects_frozen_background_policy_overriding_ambient_env(self):
        self.complete_scrublet()
        self.prepare("decontx", background_policy="exclude_all_original_called_cells")
        self.captured_envs = []
        with patch.dict(os.environ, {"DECONTX_BACKGROUND_POLICY": "historical"}):
            self.worker(1, "decontx")
        self.assertEqual(len(self.captured_envs), 1)
        self.assertEqual(self.captured_envs[0]["DECONTX_BACKGROUND_POLICY"],
                          "exclude_all_original_called_cells")

    def test_worker_sets_policy_even_with_no_ambient_value(self):
        self.complete_scrublet()
        self.prepare("decontx", background_policy="historical")
        self.captured_envs = []
        with patch.dict(os.environ, {}, clear=False):
            os.environ.pop("DECONTX_BACKGROUND_POLICY", None)
            self.worker(1, "decontx")
        self.assertEqual(self.captured_envs[0]["DECONTX_BACKGROUND_POLICY"], "historical")

    def test_tampered_output_blocks_decontx(self):
        self.complete_scrublet()
        with (self.base / "doublets/scrublet_sample_14.csv").open("a") as handle:
            handle.write("changed")
        with self.assertRaisesRegex(ValueError, "changed task output"):
            self.prepare("decontx")

    def test_changed_snapshot_rejected(self):
        self.prepare()
        (p.state_dir(self.base, "scrublet") / "input_manifest.tsv").write_text("changed")
        with self.assertRaisesRegex(ValueError, "snapshot changed"):
            self.worker(1)

    def test_wrong_rr1316_branch_rejected(self):
        for row in self.rows[12:]:
            for key in ("filtered_dir", "raw_dir"):
                wrong = self.root / "wrong_RR1316" / str(row["sample_id"]) / key
                wrong.parent.mkdir(parents=True, exist_ok=True)
                wrong.mkdir()
                for name in ("matrix.mtx.gz", "barcodes.tsv.gz", "features.tsv.gz"):
                    (wrong / name).write_bytes(b"wrong input")
                row[key] = str(wrong)
        csv_write(self.manifest, self.rows, "\t")
        with self.assertRaisesRegex(ValueError, "RR1316_L007L008"):
            self.prepare()

    def test_manifest_mapping_and_status_rejected(self):
        self.rows[0]["pool"] = "RR1316"
        csv_write(self.manifest, self.rows, "\t")
        with self.assertRaisesRegex(ValueError, "Wrong pool"):
            p.manifest(self.manifest)
        self.rows[0]["pool"] = "RR0104"
        self.rows[0]["qc_status"] = "RUNNING"
        csv_write(self.manifest, self.rows, "\t")
        with self.assertRaisesRegex(ValueError, "not passed"):
            p.manifest(self.manifest)

    def test_unreviewed_gate_and_reprepare_rejected(self):
        self.prepare()
        with self.assertRaises(FileExistsError):
            p.prepare(self.args())
        args = self.args()
        gate = p.read_json(args.gate)
        gate["status"] = "BLOCKED"
        args.gate.write_text(json.dumps(gate))
        with self.assertRaisesRegex(ValueError, "Gate must"):
            p.prepare(args)

    def test_launcher_dry_run_creates_no_analysis(self):
        env = dict(os.environ, EOS_ATLAS_BASE=str(self.base), PYTHON_EXECUTABLE=sys.executable)
        result = subprocess.run(["bash", str(SCRIPTS / "submit_preprocessing_array.sh"),
                                 "scrublet", str(self.manifest), str(self.root / "not-yet-reviewed.json")],
                                env=env, text=True, capture_output=True, check=True)
        self.assertIn("array=1-16%8", result.stdout)
        self.assertIn("no jobs submitted", result.stdout)
        self.assertFalse(self.base.exists())

    def test_launcher_submits_only_array_and_dependent_collector(self):
        self.prepare(execute=False)
        fake_bin = self.root / "bin"
        fake_bin.mkdir()
        fake = fake_bin / "sbatch"
        fake.write_text("#!" + sys.executable + "\n" +
                        "import json, os, pathlib, sys\n"
                        "p = pathlib.Path(os.environ['FAKE_SBATCH_LOG'])\n"
                        "n = len(p.read_text().splitlines()) if p.exists() else 0\n"
                        "with p.open('a') as f: f.write(json.dumps(sys.argv[1:]) + '\\n')\n"
                        "print(100 + n)\n")
        fake.chmod(0o755)
        log = self.root / "sbatch.jsonl"
        env = dict(os.environ, EOS_ATLAS_BASE=str(self.base), PYTHON_EXECUTABLE=sys.executable,
                   PATH=str(fake_bin) + os.pathsep + os.environ["PATH"], FAKE_SBATCH_LOG=str(log))
        subprocess.run(["bash", str(SCRIPTS / "submit_preprocessing_array.sh"), "scrublet",
                        str(self.manifest), str(self.args().gate), "--submit"],
                       env=env, capture_output=True, text=True, check=True)
        calls = [json.loads(line) for line in log.read_text().splitlines()]
        self.assertEqual(len(calls), 2)
        self.assertIn("--array=1-16%8", calls[0])
        self.assertIn("--dependency=afterok:100", calls[1])
        self.assertIn("worker", calls[0])
        self.assertIn("collect", calls[1])
        self.assertEqual((p.state_dir(self.base, "scrublet") / "job_ids.tsv").read_text(),
                         "array\t100\ncollector\t101\n")


if __name__ == "__main__":
    unittest.main()
