#!/usr/bin/env python3
"""Run Scrublet on uncorrected Cell Ranger filtered counts, one sample at a time.

This stage deliberately precedes DecontX. It writes calls only; it does not
modify the Cell Ranger matrices. The next stage, 00_decontx_per_sample.R,
requires these calls and excludes predicted doublets before DecontX fitting.
"""

from __future__ import annotations

import argparse
import csv
import json
import platform
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

SAMPLES = [
    (1, "WT_Male_HFD_1", "RR0104"),
    (2, "WT_Female_HFD_1", "RR0104"),
    (3, "KO_Male_HFD_1", "RR0104"),
    (4, "KO_Female_HFD_1", "RR0104"),
    (5, "WT_Female_HFDInfected_1", "RR0508"),
    (6, "WT_Female_HFDInfected_2", "RR0508"),
    (7, "KO_Female_HFDInfected_1", "RR0508"),
    (8, "KO_Female_HFDInfected_2", "RR0508"),
    (9, "WT_Male_HFD_2", "RR0912"),
    (10, "WT_Male_HFD_3", "RR0912"),
    (11, "WT_Female_HFD_2", "RR0912"),
    (12, "WT_Female_HFD_3", "RR0912"),
    (13, "KO_Male_HFD_2", "RR1316"),
    (14, "KO_Male_HFD_3", "RR1316"),
    (15, "KO_Female_HFD_2", "RR1316"),
    (16, "KO_Female_HFD_3", "RR1316"),
]


def arguments() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--data-root",
        default="/bigdata/nairlab/shared/XQiu/EOS_project/dataset/processed_data",
        help="Directory containing <GEM>/outs/per_sample_outs/<sample>/count",
    )
    parser.add_argument(
        "--output-dir",
        default="/bigdata/nairlab/shared/XQiu/EOS_project/PNAS_revision/atlas/doublets",
    )
    parser.add_argument(
        "--expected-rate",
        type=float,
        default=0.06,
        help="Sensitivity parameter; justify against recovered-cell loading metrics",
    )
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--sample-id", type=int, choices=range(1, 17),
                        help="Process only this sample; do not write shared summaries")
    parser.add_argument("--input-manifest", type=Path,
                        help="TSV with sample_id, pool, filtered_dir, raw_dir")
    return parser.parse_args()


def main() -> None:
    args = arguments()
    import numpy as np
    import pandas as pd
    import scanpy as sc
    import scrublet as scr

    if not 0 < args.expected_rate < 1:
        raise ValueError("--expected-rate must be between zero and one")
    manifest = None
    if args.input_manifest:
        with args.input_manifest.open() as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        manifest = {int(row["sample_id"]): row for row in rows}
        if len(rows) != 16 or set(manifest) != set(range(1, 17)):
            raise ValueError("Input manifest must contain exactly samples 1–16")
        for sid, _, pool in SAMPLES:
            if manifest[sid]["pool"] != pool:
                raise ValueError(f"Wrong pool in input manifest for sample {sid}")
    data_root = Path(args.data_root)
    output_dir = Path(args.output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    summaries: list[dict[str, object]] = []

    for sample_id, condition, gem in SAMPLES:
        if args.sample_id is not None and sample_id != args.sample_id:
            continue
        matrix_dir = (
            data_root
            / gem
            / "outs"
            / "per_sample_outs"
            / str(sample_id)
            / "count"
            / "sample_filtered_feature_bc_matrix"
        )
        if manifest is not None:
            matrix_dir = Path(manifest[sample_id]["filtered_dir"])
        if not matrix_dir.is_dir():
            raise FileNotFoundError(f"Missing Cell Ranger matrix: {matrix_dir}")

        print(f"[{sample_id:02d}/16] {condition}: {matrix_dir}", flush=True)
        adata = sc.read_10x_mtx(matrix_dir, var_names="gene_symbols", cache=False)
        adata.var_names_make_unique()
        counts = adata.X.tocsc() if hasattr(adata.X, "tocsc") else adata.X

        detector = scr.Scrublet(
            counts,
            expected_doublet_rate=args.expected_rate,
            random_state=args.seed,
        )
        scores, calls = detector.scrub_doublets()
        if len(scores) != adata.n_obs or len(calls) != adata.n_obs:
            raise RuntimeError(f"Scrublet returned the wrong cell count for sample {sample_id}")
        if not np.isfinite(scores).all():
            raise RuntimeError(f"Non-finite Scrublet scores for sample {sample_id}")

        call_table = pd.DataFrame(
            {
                "barcode": adata.obs_names.astype(str),
                "sample_id": sample_id,
                "condition": condition,
                "gem_reaction": gem,
                "scrublet_score": scores,
                "scrublet_doublet": calls.astype(bool),
            }
        )
        call_path = output_dir / f"scrublet_sample_{sample_id}.csv"
        if args.sample_id is not None and call_path.exists():
            raise FileExistsError(f"Refusing to overwrite sample output: {call_path}")
        call_table.to_csv(call_path, index=False)

        threshold = getattr(detector, "threshold_", None)
        summaries.append(
            {
                "sample_id": sample_id,
                "condition": condition,
                "gem_reaction": gem,
                "n_cells": int(adata.n_obs),
                "n_predicted_doublets": int(calls.sum()),
                "predicted_doublet_fraction": float(calls.mean()),
                "expected_doublet_rate": args.expected_rate,
                "threshold": None if threshold is None else float(threshold),
                "seed": args.seed,
                "input_matrix": str(matrix_dir),
            }
        )
        print(
            f"  retained={int((~calls).sum())}; predicted_doublets={int(calls.sum())} "
            f"({100 * calls.mean():.2f}%); threshold={threshold}",
            flush=True,
        )

    summary = pd.DataFrame(summaries)
    expected = {args.sample_id} if args.sample_id is not None else set(range(1, 17))
    if set(summary["sample_id"]) != expected:
        raise RuntimeError("Doublet analysis did not complete for all selected samples")
    summary_name = (f"scrublet_summary_sample_{args.sample_id}.csv"
                    if args.sample_id is not None else "scrublet_summary.csv")
    summary.to_csv(output_dir / summary_name, index=False)
    package_versions = {}
    for package in ("scanpy", "scrublet", "numpy", "pandas", "scipy"):
        try:
            package_versions[package] = version(package)
        except PackageNotFoundError:
            package_versions[package] = "unknown"
    run_record = {
        "parameters": {k: str(v) if isinstance(v, Path) else v for k, v in vars(args).items()},
        "python": platform.python_version(),
        "packages": package_versions,
    }
    record_name = (f"scrublet_run_parameters_sample_{args.sample_id}.json"
                   if args.sample_id is not None else "scrublet_run_parameters.json")
    (output_dir / record_name).write_text(
        json.dumps(run_record, indent=2) + "\n"
    )
    print(f"Completed {len(expected)} sample(s). Outputs: {output_dir}")


if __name__ == "__main__":
    main()
