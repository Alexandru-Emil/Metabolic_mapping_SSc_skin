"""Original unsupervised six-donor ResolVI run (17922984), with portable I/O.
Yi-Nan Li, source commit c0b9ce60ca19e8ab2fdeee66c4c7ea9e18713cf0:
https://github.com/ynanli/ssc-metabolic/blob/c0b9ce60ca19e8ab2fdeee66c4c7ea9e18713cf0/scripts/resolvi.py
"""
import argparse
from pathlib import Path


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--model-dir", type=Path)
    parser.add_argument("--batch-key", default="SampleId")
    parser.add_argument("--allow-package-changes", action="store_true")
    args = parser.parse_args()
    model_dir = args.model_dir or args.output.with_suffix(".model")
    if args.output.exists() or model_dir.exists():
        parser.error("Choose new output and model paths; existing results will not be overwritten.")

    import importlib.metadata
    import numpy as np
    import scanpy as sc
    import scvi
    import torch
    from scipy import sparse

    reference = {"scvi-tools": "1.4.1", "torch": "2.4.1", "scanpy": "1.11.5"}
    if not args.allow_package_changes:
        for package, expected in reference.items():
            version = importlib.metadata.version(package).split("+")[0]
            if version != expected:
                raise ValueError(f"Use {package} {expected} (original run); found {version}.")
    torch.set_float32_matmul_precision("high")
    scvi.settings.seed = 0
    adata = sc.read_h5ad(args.input)
    if "counts" not in adata.layers or "X_spatial" not in adata.obsm:
        raise ValueError("Input requires layers['counts'] and obsm['X_spatial'].")
    if args.batch_key not in adata.obs:
        raise ValueError(f"Missing batch column {args.batch_key}.")
    samples = adata.obs[args.batch_key]
    if samples.isna().any() or samples.astype(str).eq("").any() or samples.nunique() != 6:
        raise ValueError("Input must contain exactly the original six donors.")
    if not adata.obs_names.is_unique or not adata.var_names.is_unique:
        raise ValueError("Cell and gene identifiers must be unique.")
    values = adata.layers["counts"]
    values = values.data if sparse.issparse(values) else np.asarray(values)
    if not np.isfinite(values).all() or (values < 0).any() or not np.equal(values, np.floor(values)).all():
        raise ValueError("The counts layer must contain finite nonnegative raw integer counts.")
    coords = np.asarray(adata.obsm["X_spatial"])
    if coords.shape != (adata.n_obs, 2) or not np.isfinite(coords).all():
        raise ValueError("X_spatial must contain finite x/y coordinates for every cell.")
    adata.obs[args.batch_key] = samples.astype("category")
    scvi.external.RESOLVI.setup_anndata(
        adata, layer="counts", batch_key=args.batch_key,
        categorical_covariate_keys=None, labels_key=None,
    )
    model = scvi.external.RESOLVI(
        adata, n_layers=2, n_latent=30, gene_likelihood="zinb",
        dispersion="gene-batch", semisupervised=False,
    )
    model.train(max_epochs=400)
    adata.obsm["X_resolvi"] = model.get_latent_representation()
    sc.pp.neighbors(adata, use_rep="X_resolvi")
    sc.tl.umap(adata)
    adata.layers["scvi_norm_expr"] = model.get_normalized_expression(n_samples=5)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    adata.write_h5ad(args.output)
    model_dir.parent.mkdir(parents=True, exist_ok=True)
    model.save(str(model_dir), overwrite=False)
    print(f"Saved {args.output} and {model_dir}")


if __name__ == "__main__":
    main()
