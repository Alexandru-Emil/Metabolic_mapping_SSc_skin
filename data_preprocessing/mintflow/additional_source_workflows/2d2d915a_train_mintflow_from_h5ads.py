#!/usr/bin/env python

from __future__ import annotations

import argparse
import pickle
from pathlib import Path

import mintflow
import numpy as np
import scanpy as sc
import scipy.sparse as sp
import torch
from tqdm.auto import tqdm


OBS_CELLTYPE = "celltype_for_MintFlow"
OBS_SLICE = "TissueSectionID_for_MintFlow"
OBS_BATCH = "batchID_for_MintFlow"
OBS_X = "x_centroid"
OBS_Y = "y_centroid"


def check_and_rewrite_h5ad(path: str) -> None:
    adata = sc.read_h5ad(path)
    for col in [OBS_CELLTYPE, OBS_SLICE, OBS_BATCH, OBS_X, OBS_Y]:
        if col not in adata.obs:
            raise ValueError(f"{path}: missing required adata.obs column {col!r}")

    values = adata.X.data if sp.issparse(adata.X) else np.asarray(adata.X)
    if values.size and (values.min() < 0 or not np.allclose(values, np.round(values))):
        raise ValueError(f"{path}: adata.X does not look like raw non-negative integer counts")

    if sp.issparse(adata.X):
        adata.X.data = np.round(adata.X.data).astype(np.int64)
    else:
        adata.X = np.round(adata.X).astype(np.int64)

    adata.obsm["spatial"] = adata.obs[[OBS_X, OBS_Y]].to_numpy()
    for col in [OBS_CELLTYPE, OBS_SLICE, OBS_BATCH]:
        adata.obs[col] = adata.obs[col].astype("category")
    adata.write_h5ad(path)


def configure_tissues(config: dict, h5ads: list[str], width_window: int, n_neighs: int, train: bool) -> None:
    loader_key = "config_dataloader_train" if train else "config_dataloader_test"
    for i, path in enumerate(h5ads, start=1):
        tissue = config["list_tissue"][f"anndata{i}"]
        tissue["file"] = path
        tissue["obskey_cell_type"] = OBS_CELLTYPE
        tissue["obskey_sliceid_to_checkUnique"] = OBS_SLICE
        tissue["obskey_x"] = OBS_X
        tissue["obskey_y"] = OBS_Y
        tissue["obskey_biological_batch_key"] = OBS_BATCH
        tissue[loader_key]["width_window"] = width_window
        tissue["config_neighbourhood_graph"] = {
            "n_neighs": n_neighs,
            "set_diag": "False",
            "delaunay": "False",
        }


def main() -> None:
    parser = argparse.ArgumentParser(description="Train MintFlow from SFE-exported Xenium h5ad files.")
    parser.add_argument("--h5ad-list", default="mintflow_h5ad/h5ad_files.txt")
    parser.add_argument("--out-dir", default="mintflow_outputs")
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--width-window", type=int, default=300)
    parser.add_argument("--n-neighs", type=int, default=10)
    parser.add_argument("--wandb", action="store_true")
    args = parser.parse_args()

    h5ads = [line.strip() for line in Path(args.h5ad_list).read_text().splitlines() if line.strip()]
    if not h5ads:
        raise ValueError("No h5ad files found in --h5ad-list")
    for path in h5ads:
        check_and_rewrite_h5ad(path)

    out_dir = Path(args.out_dir)
    out_dir.mkdir(parents=True, exist_ok=True)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")

    config_data_train, config_data_evaluation, config_model, config_training = mintflow.get_default_configurations(
        num_tissue_sections_training=len(h5ads),
        num_tissue_sections_evaluation=len(h5ads),
    )
    configure_tissues(config_data_train, h5ads, args.width_window, args.n_neighs, train=True)
    configure_tissues(config_data_evaluation, h5ads, args.width_window, args.n_neighs, train=False)

    config_training["num_training_epochs"] = args.epochs
    config_training["flag_use_GPU"] = "True" if torch.cuda.is_available() else "False"
    config_training["flag_enable_wandb"] = "True" if args.wandb else "False"
    config_training["wandb_project_name"] = "MintFlow"
    config_training["wandb_run_name"] = "MintFlow_Xenium_SFE"

    config_data_train = mintflow.verify_and_postprocess_config_data_train(config_data_train)
    config_data_evaluation = mintflow.verify_and_postprocess_config_data_evaluation(config_data_evaluation)
    config_model = mintflow.verify_and_postprocess_config_model(config_model, num_tissue_sections=len(config_data_train))
    config_training = mintflow.verify_and_postprocess_config_training(config_training)

    dict_all4_configs = {
        "config_data_train": config_data_train,
        "config_data_evaluation": config_data_evaluation,
        "config_model": config_model,
        "config_training": config_training,
    }
    data_mintflow = mintflow.setup_data(dict_all4_configs=dict_all4_configs)
    model = mintflow.setup_model(dict_all4_configs=dict_all4_configs, data_mintflow=data_mintflow)
    trainer = mintflow.Trainer(dict_all4_configs=dict_all4_configs, model=model, data_mintflow=data_mintflow)

    for epoch in tqdm(range(config_training["num_training_epochs"]), desc="Training epoch"):
        trainer.train_one_epoch()
        predictions = mintflow.predict(
            device=device,
            dict_all4_configs=dict_all4_configs,
            data_mintflow=data_mintflow,
            model=model,
            evalulate_on_sections="all",
        )
        with (out_dir / f"predictions_epoch_{epoch}.pkl").open("wb") as handle:
            pickle.dump(predictions, handle)

        df_eval = mintflow.evaluate_by_known_signalling_genes(
            device=device,
            dict_all4_configs=dict_all4_configs,
            data_mintflow=data_mintflow,
            model=model,
            evalulate_on_sections="all",
            optional_list_colvaltype_toadd=[["training_epoch", epoch, "category"]],
        )
        df_eval.to_pickle(out_dir / f"df_evaluation_result_epoch_{epoch}.pkl")
        mintflow.dump_checkpoint(
            model=model,
            data_mintflow=data_mintflow,
            dict_all4_configs=dict_all4_configs,
            path_dump=out_dir / f"checkpoint_epoch_{epoch}.pt",
        )


if __name__ == "__main__":
    main()
