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
    parser.add_argument("--diagnostics", action="store_true",
                        help="Save observed/simulated scores and plots; no parameter changes")
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
        if args.diagnostics:
            # Diagnostics are produced AFTER fitting; they do not consume model RNG
            # or modify scores, calls, preprocessing, or the automatic threshold.
            simulated = np.asarray(detector.doublet_scores_sim_)
            if threshold is None or not np.isfinite(threshold) or not np.isfinite(simulated).all():
                raise RuntimeError("Missing/nonfinite automatic threshold or simulated scores")
            prefix = output_dir / f"scrublet_diagnostics_sample_{sample_id}"
            extra_paths = [prefix.with_suffix(ext) for ext in (".npz", ".json", ".png")]
            if any(path.exists() for path in extra_paths):
                raise FileExistsError("Refusing to overwrite diagnostic outputs")
            np.savez_compressed(extra_paths[0], observed_scores=scores,
                                simulated_scores=simulated, barcodes=adata.obs_names.to_numpy(dtype=str))
            diagnostic = {
                "sample_id": sample_id, "scope": "DIAGNOSTIC_ONLY_NO_CELL_REMOVAL",
                "n_observed": int(len(scores)), "n_simulated": int(len(simulated)),
                "threshold": float(threshold), "expected_rate": args.expected_rate,
                "seed": args.seed, "predicted_doublet_fraction": float(calls.mean()),
                "detectable_doublet_fraction": float(detector.detectable_doublet_fraction_),
                "overall_doublet_rate_estimate": float(detector.overall_doublet_rate_),
                "rate_review_flag": bool(calls.mean() < 0.005 or calls.mean() > 0.15),
                "observed_quantiles": np.quantile(scores, [0, .5, .95, .99, 1]).tolist(),
                "simulated_quantiles": np.quantile(simulated, [0, .5, .95, .99, 1]).tolist(),
            }
            with extra_paths[1].open("x") as handle:
                json.dump(diagnostic, handle, indent=2, allow_nan=False)
                handle.write("\n")
            import matplotlib
            matplotlib.use("Agg")
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(figsize=(7, 4))
            edges = np.linspace(0, max(float(max(scores)), float(max(simulated)), float(threshold)) * 1.05, 61)
            ax.hist(scores, bins=edges, density=True, histtype="step", label="Observed cells")
            ax.hist(simulated, bins=edges, density=True, histtype="step", label="Simulated doublets")
            ax.axvline(threshold, color="black", linestyle="--", label=f"Auto threshold {threshold:.3f}")
            ax.set(xlabel="Scrublet score", ylabel="Density", title=f"Sample {sample_id}: diagnostic only; no cells removed")
            ax.legend()
            fig.tight_layout()
            fig.savefig(extra_paths[2], dpi=160)
            plt.close(fig)
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
