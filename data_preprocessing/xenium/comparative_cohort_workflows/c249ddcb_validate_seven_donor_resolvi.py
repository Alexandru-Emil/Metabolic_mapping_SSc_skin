#!/usr/bin/env python3
"""Validate a completed seven-donor ResolVI run without loading dense layers."""

from __future__ import annotations

import argparse
import hashlib
import json
import re
from pathlib import Path

import anndata as ad
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-h5ad", required=True, type=Path)
    parser.add_argument("--run-dir", required=True, type=Path)
    parser.add_argument("--training-log", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--chunk-size", type=int, default=256)
    return parser.parse_args()


def sha256_file(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def layer_block(layer: object, start: int, stop: int) -> np.ndarray | sp.spmatrix:
    block = layer[start:stop]
    if hasattr(block, "to_memory"):
        block = block.to_memory()
    return block


def block_values(block: np.ndarray | sp.spmatrix) -> np.ndarray:
    if sp.issparse(block):
        return block.data
    return np.asarray(block)


def summarize_layer(
    layer: object,
    n_obs: int,
    chunk_size: int,
    expect_probability_rows: bool = False,
) -> dict[str, float | int | bool]:
    total = 0.0
    nonzero = 0
    minimum = np.inf
    maximum = -np.inf
    finite = True
    row_sum_minimum = np.inf
    row_sum_maximum = -np.inf
    for start in range(0, n_obs, chunk_size):
        block = layer_block(layer, start, min(start + chunk_size, n_obs))
        values = block_values(block)
        finite = finite and bool(np.isfinite(values).all())
        if values.size:
            minimum = min(minimum, float(values.min()))
            maximum = max(maximum, float(values.max()))
        if sp.issparse(block):
            total += float(block.sum())
            nonzero += int(block.nnz)
            row_sums = np.asarray(block.sum(axis=1)).ravel()
        else:
            array = np.asarray(block)
            total += float(array.sum())
            nonzero += int(np.count_nonzero(array))
            row_sums = array.sum(axis=1)
        row_sum_minimum = min(row_sum_minimum, float(row_sums.min()))
        row_sum_maximum = max(row_sum_maximum, float(row_sums.max()))

    probability_rows_valid = True
    if expect_probability_rows:
        probability_rows_valid = (
            abs(row_sum_minimum - 1.0) <= 1e-4
            and abs(row_sum_maximum - 1.0) <= 1e-4
        )
    return {
        "finite": finite,
        "minimum": minimum,
        "maximum": maximum,
        "sum": total,
        "nonzero": nonzero,
        "row_sum_minimum": row_sum_minimum,
        "row_sum_maximum": row_sum_maximum,
        "probability_rows_valid": probability_rows_valid,
    }


def parse_elbo(log_path: Path) -> pd.DataFrame:
    text = log_path.read_text(errors="replace").replace("\r", "\n")
    pattern = re.compile(
        r"Epoch\s+(\d+)/\d+:[^\n]*elbo_train=([0-9.eE+-]+)"
    )
    values: dict[int, float] = {}
    for epoch, value in pattern.findall(text):
        values[int(epoch)] = float(value)
    return pd.DataFrame(
        {"epoch": list(values.keys()), "elbo_train": list(values.values())}
    ).sort_values("epoch")


def main() -> None:
    args = parse_args()
    input_path = args.input_h5ad.resolve()
    run_dir = args.run_dir.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    model_h5ad_path = run_dir / "seven_donor_resolvi.h5ad"
    manifest_path = run_dir / "run_manifest.json"
    model_path = run_dir / "model" / "model.pt"

    for path in [
        input_path,
        model_h5ad_path,
        manifest_path,
        model_path,
        args.training_log,
    ]:
        if not path.is_file():
            raise FileNotFoundError(path)

    manifest = json.loads(manifest_path.read_text())
    if manifest.get("status") != "complete":
        raise ValueError("The training manifest is not complete.")
    expected_config = {
        "sample_key": "SampleId",
        "slide_key": "SlideId",
        "counts_layer": "counts",
        "spatial_key": "X_spatial",
        "n_latent": 30,
        "max_epochs": 100,
        "batch_size": 512,
        "seed": 20260729,
        "include_corrected_px_rate": True,
        "posterior_samples": 3,
        "smoke_test_cells_per_donor": 0,
    }
    config = manifest.get("config", {})
    config_mismatches = {
        key: {"expected": value, "observed": config.get(key)}
        for key, value in expected_config.items()
        if config.get(key) != value
    }
    expected_donor_set = {
        "Validation1",
        "Validation2",
        "ExcludedValidation",
        "Validation3",
        "Validation4",
        "Validation5",
        "Validation6",
    }
    if set(config.get("donors", [])) != expected_donor_set:
        config_mismatches["donors"] = {
            "expected": sorted(expected_donor_set),
            "observed": config.get("donors"),
        }
    if config_mismatches:
        raise ValueError(f"Unexpected production configuration: {config_mismatches}")
    if sha256_file(input_path) != manifest["config"]["input_sha256"]:
        raise ValueError("Input H5AD checksum differs from the training manifest.")
    output_sha256 = sha256_file(model_h5ad_path)
    if output_sha256 != manifest["output_h5ad_sha256"]:
        raise ValueError("Output H5AD checksum differs from the training manifest.")

    source = ad.read_h5ad(input_path, backed="r")
    result = ad.read_h5ad(model_h5ad_path, backed="r")
    if source.shape != result.shape or source.shape != (21271, 5099):
        raise ValueError(f"Unexpected input/output shapes: {source.shape}, {result.shape}")
    if not source.obs_names.equals(result.obs_names):
        raise ValueError("Input and output cell identifiers differ.")
    if not source.var_names.equals(result.var_names):
        raise ValueError("Input and output gene identifiers differ.")

    required_layers = {"counts", "scvi_norm_expr", "px_rate_q50"}
    missing_layers = sorted(required_layers - set(result.layers.keys()))
    if missing_layers:
        raise KeyError(f"Missing output layers: {missing_layers}")
    for key in ["X_spatial", "X_resolvi"]:
        if key not in result.obsm:
            raise KeyError(f"Missing output representation: {key}")
    required_scvi_obs = {"_scvi_batch", "_scvi_labels"}
    missing_scvi_obs = sorted(required_scvi_obs - set(result.obs.columns))
    if missing_scvi_obs:
        raise KeyError(f"Missing scvi registry columns: {missing_scvi_obs}")
    n_scvi_batches = int(result.obs["_scvi_batch"].nunique())
    n_scvi_labels = int(result.obs["_scvi_labels"].nunique())
    if n_scvi_batches != 7:
        raise ValueError(f"Expected seven scvi batches, observed {n_scvi_batches}.")
    if n_scvi_labels != 1:
        raise ValueError(
            f"Expected one unsupervised scvi label level, observed {n_scvi_labels}."
        )

    input_counts = summarize_layer(
        source.layers["counts"], source.n_obs, args.chunk_size
    )
    output_counts = summarize_layer(
        result.layers["counts"], result.n_obs, args.chunk_size
    )
    if input_counts["sum"] != output_counts["sum"]:
        raise ValueError("Input and output count sums differ.")
    if input_counts["nonzero"] != output_counts["nonzero"]:
        raise ValueError("Input and output count nonzero totals differ.")
    if input_counts["minimum"] < 0:
        raise ValueError("Negative raw counts detected.")

    normalized = summarize_layer(
        result.layers["scvi_norm_expr"],
        result.n_obs,
        args.chunk_size,
        expect_probability_rows=True,
    )
    corrected = summarize_layer(
        result.layers["px_rate_q50"],
        result.n_obs,
        args.chunk_size,
    )
    if not normalized["finite"] or not normalized["probability_rows_valid"]:
        raise ValueError("scvi_norm_expr failed finite/probability-row checks.")
    if not corrected["finite"] or corrected["minimum"] < 0:
        raise ValueError("px_rate_q50 failed finite/non-negative checks.")

    source_spatial = np.asarray(source.obsm["X_spatial"])
    result_spatial = np.asarray(result.obsm["X_spatial"])
    if not np.array_equal(source_spatial, result_spatial):
        raise ValueError("Spatial coordinates changed between input and output.")
    latent = np.asarray(result.obsm["X_resolvi"])
    latent_sd = latent.std(axis=0)
    if latent.shape != (21271, 30) or not np.isfinite(latent).all():
        raise ValueError("ResolVI latent failed shape/finite checks.")
    if np.any(latent_sd <= 1e-6):
        raise ValueError("At least one ResolVI latent dimension is collapsed.")

    donor_counts = result.obs["SampleId"].astype(str).value_counts(sort=False)
    expected_donors = {
        "Validation1": 2872,
        "Validation2": 3934,
        "ExcludedValidation": 1880,
        "Validation3": 1775,
        "Validation4": 3413,
        "Validation5": 4174,
        "Validation6": 3223,
    }
    if donor_counts.to_dict() != expected_donors:
        raise ValueError(f"Unexpected donor counts: {donor_counts.to_dict()}")

    elbo = parse_elbo(args.training_log.resolve())
    elbo.to_csv(output_dir / "elbo_train_from_console.csv", index=False)
    if (
        len(elbo) < 99
        or int(elbo["epoch"].max()) != 100
        or not np.isfinite(elbo["elbo_train"]).all()
        or float(elbo.iloc[-1]["elbo_train"])
        >= float(elbo.iloc[0]["elbo_train"])
    ):
        raise ValueError(
            "Console log does not contain a complete, finite, improving "
            "100-epoch ELBO trajectory."
        )
    figure, axis = plt.subplots(figsize=(6.0, 3.8))
    axis.plot(elbo["epoch"], elbo["elbo_train"], color="#2C6E9B", linewidth=1.8)
    axis.set(xlabel="Epoch", ylabel="Training ELBO", title="ResolVI training trajectory")
    axis.grid(axis="y", color="#DDDDDD", linewidth=0.6)
    figure.tight_layout()
    figure.savefig(output_dir / "resolvi_training_elbo.pdf")
    figure.savefig(output_dir / "resolvi_training_elbo.png", dpi=300)
    plt.close(figure)

    validation = {
        "status": "passed",
        "input_h5ad": str(input_path),
        "input_h5ad_sha256": sha256_file(input_path),
        "output_h5ad": str(model_h5ad_path),
        "output_h5ad_sha256": output_sha256,
        "model_pt": str(model_path),
        "model_pt_sha256": sha256_file(model_path),
        "shape": list(result.shape),
        "validated_config": expected_config,
        "n_scvi_batches": n_scvi_batches,
        "n_scvi_label_levels": n_scvi_labels,
        "donor_counts": expected_donors,
        "input_counts": input_counts,
        "output_counts": output_counts,
        "scvi_norm_expr": normalized,
        "px_rate_q50": corrected,
        "latent_sd_minimum": float(latent_sd.min()),
        "latent_sd_maximum": float(latent_sd.max()),
        "n_elbo_epochs_parsed": int(len(elbo)),
        "elbo_first": float(elbo.iloc[0]["elbo_train"]),
        "elbo_last": float(elbo.iloc[-1]["elbo_train"]),
    }
    (output_dir / "resolvi_validation.json").write_text(
        json.dumps(validation, indent=2, sort_keys=True) + "\n"
    )
    source.file.close()
    result.file.close()
    print(json.dumps(validation, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
