#!/usr/bin/env python

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scipy.io
import scipy.sparse as sp


OBS_CELLTYPE = "celltype_for_MintFlow"
OBS_SLICE = "TissueSectionID_for_MintFlow"
OBS_BATCH = "batchID_for_MintFlow"
OBS_X = "x_centroid"
OBS_Y = "y_centroid"


def read_csv_with_index(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, index_col=0)


def convert_section(section_dir: Path, out_file: Path) -> None:
    obs = read_csv_with_index(section_dir / "obs.csv")
    var = read_csv_with_index(section_dir / "var.csv")
    x = scipy.io.mmread(section_dir / "counts_cells_by_genes.mtx")
    x = sp.csr_matrix(x)
    if x.data.size and (x.data.min() < 0 or not np.allclose(x.data, np.round(x.data))):
        raise ValueError(f"{section_dir}: counts do not look like raw non-negative integers")
    x.data = np.round(x.data).astype(np.int64)

    required = [OBS_CELLTYPE, OBS_SLICE, OBS_BATCH, OBS_X, OBS_Y]
    missing = [col for col in required if col not in obs.columns]
    if missing:
        raise ValueError(f"{section_dir}: missing required obs columns: {', '.join(missing)}")

    adata = ad.AnnData(X=x, obs=obs, var=var)
    adata.obs_names = adata.obs_names.astype(str)
    adata.var_names = adata.var_names.astype(str)
    adata.var_names_make_unique()
    adata.obsm["spatial"] = adata.obs[[OBS_X, OBS_Y]].to_numpy(dtype=float)
    for col in [OBS_CELLTYPE, OBS_SLICE, OBS_BATCH]:
        adata.obs[col] = adata.obs[col].astype("category")
    adata.write_h5ad(out_file)


def main() -> None:
    parser = argparse.ArgumentParser(description="Convert SFE Matrix Market exports into MintFlow h5ad files.")
    parser.add_argument("--mtx-dir", default="mintflow_mtx")
    parser.add_argument("--out-dir", default="mintflow_h5ad")
    args = parser.parse_args()

    mtx_dir = Path(args.mtx_dir)
    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    section_dirs = sorted(path for path in mtx_dir.glob("tissue_section_*") if path.is_dir())
    if not section_dirs:
        raise ValueError(f"No tissue_section_* directories found in {mtx_dir}")

    h5ad_files: list[str] = []
    for i, section_dir in enumerate(section_dirs, start=1):
        out_file = out_dir / f"tissue_section_{i:03d}.h5ad"
        convert_section(section_dir, out_file)
        h5ad_files.append(str(out_file.resolve()).replace("\\", "/"))

    (out_dir / "h5ad_files.txt").write_text("\n".join(h5ad_files) + "\n", encoding="utf-8")
    print(f"Wrote {len(h5ad_files)} h5ad file(s) to {out_dir}")


if __name__ == "__main__":
    main()
