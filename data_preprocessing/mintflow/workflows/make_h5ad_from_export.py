# Scientific workflow derived from Revision/Nan annotation/Data/RMD files/MintFlow/00_scripts_and_reproducibility/make_h5ad_from_export.py
import os
from pathlib import Path
def project_path(relative):
    return str(Path(os.environ.get("METABOLIC_INPUT_DIR", "inputs")) / relative)
#!/usr/bin/env python

import argparse
import os

import anndata as ad
import pandas as pd
from scipy.io import mmread


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundle-dir", required=True)
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    counts = mmread(os.path.join(args.bundle_dir, "counts_genes_by_cells.mtx")).tocsr()
    obs = pd.read_csv(os.path.join(args.bundle_dir, "obs.csv"), dtype=str)
    var = pd.read_csv(os.path.join(args.bundle_dir, "var.csv"), dtype=str)

    keep_obs = [
        "cell_id",
        "metfiblabel",
        "metEClabel",
        "combined_met_label",
        "mintflow_slice_id",
        "mintflow_batch",
        "x_centroid",
        "y_centroid",
        "KEY",
        "sample_id",
        "lv2_anno_new",
    ]
    obs = obs[[col for col in keep_obs if col in obs.columns]].copy()
    obs.index = obs["cell_id"].astype(str)
    var.index = var["gene_id"].astype(str)

    for col in ["x_centroid", "y_centroid"]:
        obs[col] = pd.to_numeric(obs[col])

    adata = ad.AnnData(X=counts.T.tocsr(), obs=obs, var=var)
    adata.obsm["spatial"] = obs[["x_centroid", "y_centroid"]].to_numpy()
    os.makedirs(os.path.dirname(args.out), exist_ok=True)
    adata.write_h5ad(args.out)
    print(adata)
    print(obs["metfiblabel"].value_counts(dropna=False))


if __name__ == "__main__":
    main()
