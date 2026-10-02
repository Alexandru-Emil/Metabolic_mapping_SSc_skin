# Scientific workflow derived from Revision/Nan annotation/Data/RMD files/MintFlow/00_scripts_and_reproducibility/train_mintflow_xenium.py
import os
from pathlib import Path
def project_path(relative):
    return str(Path(os.environ.get("METABOLIC_INPUT_DIR", "inputs")) / relative)
#!/usr/bin/env python

import argparse
import os
import pickle

import mintflow
import torch
from tqdm.autonotebook import tqdm


def configure_tissue(config, path_anndata, dataloader_key, width_window):
    tissue = config["list_tissue"]["anndata1"]
    tissue["file"] = path_anndata
    tissue["obskey_cell_type"] = "metfiblabel"
    tissue["obskey_sliceid_to_checkUnique"] = "mintflow_slice_id"
    tissue["obskey_x"] = "x_centroid"
    tissue["obskey_y"] = "y_centroid"
    tissue["obskey_biological_batch_key"] = "mintflow_batch"
    tissue[dataloader_key]["width_window"] = width_window
    tissue["config_neighbourhood_graph"] = {
        "n_neighs": 5,
        "set_diag": "False",
        "delaunay": "False",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--h5ad", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--obskey-cell-type", default="metfiblabel")
    parser.add_argument("--epochs", type=int, default=20)
    parser.add_argument("--width-window", type=int, default=100)
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)
    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    print(f"Using device: {device}")

    config_data_train, config_data_evaluation, config_model, config_training = mintflow.get_default_configurations(
        num_tissue_sections_training=1,
        num_tissue_sections_evaluation=1,
    )

    configure_tissue(config_data_train, args.h5ad, "config_dataloader_train", args.width_window)
    configure_tissue(config_data_evaluation, args.h5ad, "config_dataloader_test", args.width_window)
    config_data_train["list_tissue"]["anndata1"]["obskey_cell_type"] = args.obskey_cell_type
    config_data_evaluation["list_tissue"]["anndata1"]["obskey_cell_type"] = args.obskey_cell_type

    config_training["num_training_epochs"] = args.epochs
    config_training["flag_use_GPU"] = "True" if torch.cuda.is_available() else "False"
    config_training["flag_enable_wandb"] = "False"
    config_training["wandb_project_name"] = "MintFlow"
    config_training["wandb_run_name"] = "MintFlow_Xenium_SFE_09032026_metfiblabel"

    config_data_train = mintflow.verify_and_postprocess_config_data_train(config_data_train)
    config_data_evaluation = mintflow.verify_and_postprocess_config_data_evaluation(config_data_evaluation)
    config_model = mintflow.verify_and_postprocess_config_model(
        config_model,
        num_tissue_sections=len(config_data_train),
    )
    config_training = mintflow.verify_and_postprocess_config_training(config_training)

    dict_all4_configs = {
        "config_data_train": config_data_train,
        "config_data_evaluation": config_data_evaluation,
        "config_model": config_model,
        "config_training": config_training,
    }

    data_mintflow = mintflow.setup_data(dict_all4_configs=dict_all4_configs)
    model = mintflow.setup_model(dict_all4_configs=dict_all4_configs, data_mintflow=data_mintflow)
    trainer = mintflow.Trainer(
        dict_all4_configs=dict_all4_configs,
        model=model,
        data_mintflow=data_mintflow,
    )

    for index_epoch in tqdm(range(config_training["num_training_epochs"]), desc="Training epoch"):
        trainer.train_one_epoch()

        predictions = mintflow.predict(
            device=device,
            dict_all4_configs=dict_all4_configs,
            data_mintflow=data_mintflow,
            model=model,
            evalulate_on_sections="all",
        )
        with open(os.path.join(args.outdir, f"predictions_epoch_{index_epoch}.pkl"), "wb") as handle:
            pickle.dump(predictions, handle)

        try:
            df_evaluation_result = mintflow.evaluate_by_known_signalling_genes(
                device=device,
                dict_all4_configs=dict_all4_configs,
                data_mintflow=data_mintflow,
                model=model,
                evalulate_on_sections="all",
                optional_list_colvaltype_toadd=[["training_epoch", index_epoch, "category"]],
            )
            df_evaluation_result.to_pickle(
                os.path.join(args.outdir, f"df_evaluation_result_epoch_{index_epoch}.pkl")
            )
        except Exception as exc:
            print(f"Evaluation skipped for epoch {index_epoch}: {exc}")

        mintflow.dump_checkpoint(
            model=model,
            data_mintflow=data_mintflow,
            dict_all4_configs=dict_all4_configs,
            path_dump=os.path.join(args.outdir, f"checkpoint_epoch_{index_epoch}.pt"),
        )

    mintflow.dump_checkpoint(
        model=model,
        data_mintflow=data_mintflow,
        dict_all4_configs=dict_all4_configs,
        path_dump=os.path.join(args.outdir, "checkpoint_final.pt"),
    )
    print(f"Finished training. Outputs written to: {args.outdir}")


if __name__ == "__main__":
    main()
