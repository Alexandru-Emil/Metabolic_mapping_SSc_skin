#!/usr/bin/env python3
"""Train and export an unsupervised ResolVI model for the seven-donor cohort."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import os
import platform
import sys
import time
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
import scipy
import scipy.sparse as sp
import scvi
import torch
from scvi.external import RESOLVI


DEFAULT_DONORS = [
    "Validation1",
    "Validation2",
    "ExcludedValidation",
    "Validation3",
    "Validation4",
    "Validation5",
    "Validation6",
]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "Fit an unsupervised ResolVI model on selected Xenium tissue sections "
            "and export the latent and normalized expression."
        )
    )
    parser.add_argument("--input-h5ad", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--donors", nargs="+", default=DEFAULT_DONORS)
    parser.add_argument("--sample-key", default="SampleId")
    parser.add_argument("--slide-key", default="SlideId")
    parser.add_argument("--counts-layer", default="counts")
    parser.add_argument("--spatial-key", default="X_spatial")
    parser.add_argument("--n-latent", type=int, default=30)
    parser.add_argument("--max-epochs", type=int, default=100)
    parser.add_argument("--batch-size", type=int, default=512)
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--accelerator", default="cpu")
    parser.add_argument(
        "--device",
        default="1",
        help="Pyro/scvi device argument. Use 1 for one CPU device or 0 for GPU 0.",
    )
    parser.add_argument("--num-workers", type=int, default=0)
    parser.add_argument(
        "--include-corrected-px-rate",
        action="store_true",
        help="Also generate posterior-median corrected expression as px_rate_q50.",
    )
    parser.add_argument("--posterior-samples", type=int, default=3)
    parser.add_argument(
        "--smoke-test-cells-per-donor",
        type=int,
        default=0,
        help="Use a deterministic subset per donor. Zero uses all selected cells.",
    )
    parser.add_argument("--force", action="store_true")
    return parser.parse_args()


def sha256_file(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def stable_hash(value: object) -> str:
    payload = json.dumps(value, sort_keys=True, separators=(",", ":")).encode()
    return hashlib.sha256(payload).hexdigest()


def portable_config(config: dict[str, object]) -> dict[str, object]:
    """Exclude only the machine-specific path; input identity uses SHA-256."""
    return {
        key: value
        for key, value in config.items()
        if key != "input_h5ad"
    }


def write_json(path: Path, value: object) -> None:
    temporary = path.with_suffix(path.suffix + ".tmp")
    temporary.write_text(json.dumps(value, indent=2, sort_keys=True) + "\n")
    os.replace(temporary, path)


def package_versions() -> dict[str, str | bool]:
    return {
        "python": platform.python_version(),
        "platform": platform.platform(),
        "anndata": importlib.metadata.version("anndata"),
        "numpy": np.__version__,
        "pandas": pd.__version__,
        "scanpy": importlib.metadata.version("scanpy"),
        "scipy": scipy.__version__,
        "scvi_tools": str(scvi.__version__),
        "torch": str(torch.__version__),
        "cuda_available": torch.cuda.is_available(),
    }


def resolve_device(value: str) -> int | str:
    return int(value) if value.isdigit() else value


def validate_and_subset(adata: ad.AnnData, args: argparse.Namespace) -> ad.AnnData:
    required_obs = [args.sample_key, args.slide_key]
    missing_obs = [key for key in required_obs if key not in adata.obs]
    if missing_obs:
        raise KeyError(f"Missing required obs columns: {missing_obs}")
    if args.counts_layer not in adata.layers:
        raise KeyError(f"Missing counts layer: {args.counts_layer}")
    if args.spatial_key not in adata.obsm:
        raise KeyError(f"Missing spatial coordinates: {args.spatial_key}")
    if not adata.obs_names.is_unique:
        raise ValueError("Cell identifiers are not unique; refusing to alter them.")

    observed_donors = set(adata.obs[args.sample_key].astype(str))
    missing_donors = sorted(set(args.donors) - observed_donors)
    if missing_donors:
        raise ValueError(f"Requested donors absent from input: {missing_donors}")

    selected = adata.obs[args.sample_key].astype(str).isin(args.donors).to_numpy()
    adata = adata[selected].copy()

    if args.smoke_test_cells_per_donor > 0:
        rng = np.random.default_rng(args.seed)
        sample_values = adata.obs[args.sample_key].astype(str).to_numpy()
        sampled_indices = []
        for donor in args.donors:
            donor_indices = np.flatnonzero(sample_values == donor)
            n_cells = min(args.smoke_test_cells_per_donor, donor_indices.size)
            sampled_indices.append(rng.choice(donor_indices, n_cells, replace=False))
        adata = adata[np.concatenate(sampled_indices)].copy()

    adata.obs[args.sample_key] = (
        adata.obs[args.sample_key].astype(str).astype("category")
    )
    adata.obs[args.slide_key] = (
        adata.obs[args.slide_key].astype(str).astype("category")
    )

    coords = np.asarray(adata.obsm[args.spatial_key])
    if coords.ndim != 2 or coords.shape[1] < 2:
        raise ValueError(f"{args.spatial_key} must have at least two columns.")
    if not np.isfinite(coords).all():
        raise ValueError(f"{args.spatial_key} contains non-finite coordinates.")

    counts = adata.layers[args.counts_layer]
    values = counts.data if sp.issparse(counts) else np.asarray(counts)
    if values.size == 0 or np.nanmin(values) < 0:
        raise ValueError("Counts must be a non-empty, non-negative matrix.")
    if not np.allclose(values, np.rint(values)):
        raise ValueError("Counts layer contains non-integer values.")

    counts_by_donor = (
        adata.obs[args.sample_key]
        .value_counts(sort=False)
        .rename_axis(args.sample_key)
        .rename("n_cells")
        .reset_index()
    )
    if (counts_by_donor["n_cells"] < 11).any():
        raise ValueError("Every tissue section needs at least 11 cells for neighbors.")
    adata.uns["seven_donor_cell_counts"] = counts_by_donor.to_dict("list")
    return adata


def save_training_history(model: RESOLVI, output_dir: Path) -> None:
    history_dir = output_dir / "training_history"
    history_dir.mkdir(exist_ok=True)
    for metric, values in model.history.items():
        if hasattr(values, "to_csv"):
            values.to_csv(history_dir / f"{metric}.csv", index=True)


def main() -> None:
    args = parse_args()
    input_path = args.input_h5ad.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)

    if not input_path.is_file():
        raise FileNotFoundError(input_path)

    input_sha256 = sha256_file(input_path)
    config = {
        "input_h5ad": str(input_path),
        "input_sha256": input_sha256,
        "donors": args.donors,
        "sample_key": args.sample_key,
        "slide_key": args.slide_key,
        "counts_layer": args.counts_layer,
        "spatial_key": args.spatial_key,
        "n_latent": args.n_latent,
        "max_epochs": args.max_epochs,
        "batch_size": args.batch_size,
        "seed": args.seed,
        "accelerator": args.accelerator,
        "device": args.device,
        "num_workers": args.num_workers,
        "include_corrected_px_rate": args.include_corrected_px_rate,
        "posterior_samples": args.posterior_samples,
        "smoke_test_cells_per_donor": args.smoke_test_cells_per_donor,
    }
    config_hash = stable_hash(portable_config(config))
    manifest_path = output_dir / "run_manifest.json"
    output_h5ad = output_dir / "seven_donor_resolvi.h5ad"
    model_dir = output_dir / "model"

    if manifest_path.exists() and output_h5ad.exists() and not args.force:
        previous = json.loads(manifest_path.read_text())
        previous_config_hash = stable_hash(
            portable_config(previous.get("config", {}))
        )
        if (
            previous.get("status") == "complete"
            and previous_config_hash == config_hash
        ):
            print(f"Cache hit: {output_h5ad}")
            return
        raise FileExistsError(
            f"{output_dir} contains a different or incomplete run. Use --force "
            "or choose a new output directory."
        )

    manifest = {
        "status": "running",
        "config_hash": config_hash,
        "config": config,
        "versions": package_versions(),
        "started_unix": time.time(),
    }
    write_json(manifest_path, manifest)

    scvi.settings.seed = args.seed
    torch.set_float32_matmul_precision("high")

    print(f"Reading {input_path}", flush=True)
    adata = ad.read_h5ad(input_path)
    adata = validate_and_subset(adata, args)
    print(
        adata.obs[args.sample_key].value_counts(sort=False).to_string(),
        flush=True,
    )

    # SampleId is the tissue-section batch. SlideId retains the acquisition-slide
    # covariate without duplicating SampleId through the ROI column.
    RESOLVI.setup_anndata(
        adata,
        layer=args.counts_layer,
        batch_key=args.sample_key,
        categorical_covariate_keys=[args.slide_key],
        prepare_data_kwargs={"spatial_rep": args.spatial_key},
    )
    model = RESOLVI(
        adata,
        n_latent=args.n_latent,
        semisupervised=False,
    )

    model.train(
        max_epochs=args.max_epochs,
        batch_size=args.batch_size,
        accelerator=args.accelerator,
        device=resolve_device(args.device),
        datasplitter_kwargs={"num_workers": args.num_workers},
        enable_checkpointing=False,
        logger=False,
    )

    adata.obsm["X_resolvi"] = model.get_latent_representation(
        batch_size=args.batch_size
    ).astype(np.float32)
    adata.layers["scvi_norm_expr"] = model.get_normalized_expression(
        batch_size=args.batch_size,
        return_numpy=True,
    ).astype(np.float32)

    if args.include_corrected_px_rate:
        posterior = model.sample_posterior(
            model=model.module.model_corrected,
            return_sites=["px_rate"],
            summary_fun={"post_sample_q50": np.median},
            num_samples=args.posterior_samples,
            summary_frequency=30,
            accelerator=args.accelerator,
            device=resolve_device(args.device),
            batch_size=args.batch_size,
        )
        adata.layers["px_rate_q50"] = posterior["post_sample_q50"][
            "px_rate"
        ].astype(np.float32)

    adata.uns["resolvi_run_config"] = config
    adata.uns["resolvi_package_versions"] = package_versions()
    save_training_history(model, output_dir)
    model.save(model_dir, overwrite=True, save_anndata=False)
    adata.write_h5ad(output_h5ad, compression="gzip", compression_opts=4)

    manifest.update(
        {
            "status": "complete",
            "finished_unix": time.time(),
            "n_cells": adata.n_obs,
            "n_genes": adata.n_vars,
            "output_h5ad": str(output_h5ad),
            "output_h5ad_sha256": sha256_file(output_h5ad),
            "model_dir": str(model_dir),
        }
    )
    write_json(manifest_path, manifest)
    print(f"Completed: {output_h5ad}", flush=True)


if __name__ == "__main__":
    try:
        main()
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr, flush=True)
        raise
