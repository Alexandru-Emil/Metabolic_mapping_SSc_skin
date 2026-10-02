#!/usr/bin/env python3
"""Cluster and annotate the joint seven-donor ResolVI latent representation."""

from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from sklearn.metrics import (
    accuracy_score,
    adjusted_rand_score,
    balanced_accuracy_score,
    f1_score,
    normalized_mutual_info_score,
)
from sklearn.neighbors import KNeighborsClassifier


MAJOR_LEVELS = [
    "epithelial",
    "fibroblast",
    "pericyte",
    "endothelial",
    "muscle",
    "telocyte",
    "Schwann cell",
    "myeloid cell",
    "lymphocyte",
]
FIB_LEVELS = [
    "papillary_Fib",
    "COMP_Fib",
    "COL8A1_Fib",
    "PI16_Fib",
    "CXCL12_Fib",
    "CCL19_Fib",
    "NGFR_Fib",
    "COCH_Fib",
]
EC_LEVELS = [
    "cycling_EC",
    "ACKR1_EC",
    "EPC",
    "ACTA2_EC",
    "HEY1_EC",
    "lymphatic_EC",
]
RESOLUTIONS = [0.2, 0.4, 0.6, 0.8, 1.0, 1.2, 1.5, 2.0]
KNN_NEIGHBORS = [3, 5, 10, 15, 25, 40, 60]
KNN_WEIGHTS = ["uniform", "distance"]


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser()
    parser.add_argument("--model-h5ad", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--seed", type=int, default=20260729)
    parser.add_argument("--query-donor", default="ExcludedValidation")
    return parser.parse_args()


def sha256_file(path: Path, chunk_size: int = 16 * 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        while chunk := handle.read(chunk_size):
            digest.update(chunk)
    return digest.hexdigest()


def clean_labels(values: pd.Series) -> pd.Series:
    labels = values.astype("string")
    return labels.mask(labels.isna() | labels.isin(["", "NA", "<NA>", "nan"]))


def natural_cluster_key(value: str) -> tuple[int, float | str]:
    try:
        return (0, float(value))
    except ValueError:
        return (1, value)


def majority_cluster_map(
    clusters: pd.Series,
    labels: pd.Series,
    reference_mask: np.ndarray,
    donors: pd.Series,
) -> tuple[pd.Series, pd.DataFrame]:
    cluster_values = clusters.astype(str)
    label_values = clean_labels(labels)
    donor_values = donors.astype(str)
    reference_mask = (
        np.asarray(reference_mask, dtype=bool) & label_values.notna().to_numpy()
    )

    mapping: dict[str, str] = {}
    records: list[dict[str, object]] = []
    cluster_names = sorted(cluster_values.unique(), key=natural_cluster_key)
    for cluster in cluster_names:
        all_index = cluster_values.eq(cluster).to_numpy()
        reference_index = all_index & reference_mask
        counts = label_values[reference_index].value_counts()
        if not counts.empty:
            top_count = int(counts.max())
            top_label = sorted(counts[counts.eq(top_count)].index.astype(str))[0]
            mapping[cluster] = top_label
        else:
            top_count = 0
            top_label = pd.NA

        records.append(
            {
                "cluster": cluster,
                "assigned_label": top_label,
                "n_all": int(all_index.sum()),
                "n_reference": int(reference_index.sum()),
                "n_ssc230": int((all_index & donor_values.eq("ExcludedValidation")).sum()),
                "n_donors": int(donor_values[all_index].nunique()),
                "n_reference_donors": int(
                    donor_values[reference_index].nunique()
                ),
                "purity": (
                    top_count / int(reference_index.sum())
                    if reference_index.sum()
                    else np.nan
                ),
            }
        )

    predictions = cluster_values.map(mapping).astype("string")
    return predictions, pd.DataFrame.from_records(records)


def evaluate_resolution(
    clusters: pd.Series,
    labels: pd.Series,
    donors: pd.Series,
    query_donor: str,
    target_clusters: int,
) -> tuple[dict[str, object], pd.DataFrame]:
    label_values = clean_labels(labels)
    donor_values = donors.astype(str)
    reference_mask = (
        donor_values.ne(query_donor) & label_values.notna()
    ).to_numpy()
    predictions, cluster_stats = majority_cluster_map(
        clusters,
        label_values,
        reference_mask,
        donor_values,
    )

    query_clusters = set(clusters[donor_values.eq(query_donor)].astype(str))
    query_stats = cluster_stats[cluster_stats["cluster"].isin(query_clusters)]
    accuracy = float(
        np.mean(
            predictions[reference_mask].astype(str).to_numpy()
            == label_values[reference_mask].astype(str).to_numpy()
        )
    )
    ari = float(
        adjusted_rand_score(
            label_values[reference_mask].astype(str),
            clusters[reference_mask].astype(str),
        )
    )
    nmi = float(
        normalized_mutual_info_score(
            label_values[reference_mask].astype(str),
            clusters[reference_mask].astype(str),
        )
    )
    n_clusters = int(clusters.nunique())
    mean_donors = (
        float(query_stats["n_donors"].mean()) if len(query_stats) else 0.0
    )
    min_reference = (
        int(query_stats["n_reference"].min()) if len(query_stats) else 0
    )
    unsupported = (
        int(query_stats["n_reference"].eq(0).sum()) if len(query_stats) else 0
    )
    selection_score = (
        accuracy
        + 0.05 * ari
        + 0.03 * min(mean_donors / 7.0, 1.0)
        - 0.003 * abs(n_clusters - target_clusters)
        - 0.05 * unsupported
    )
    metrics = {
        "n_clusters": n_clusters,
        "reference_label_accuracy": accuracy,
        "adjusted_rand_index": ari,
        "normalized_mutual_information": nmi,
        "mean_donors_in_ssc230_clusters": mean_donors,
        "min_reference_cells_in_ssc230_clusters": min_reference,
        "unsupported_ssc230_clusters": unsupported,
        "selection_score": selection_score,
        "eligible": (
            unsupported == 0 and min_reference >= 20 and mean_donors >= 4
        ),
    }
    return metrics, cluster_stats


def label_missing_by_knn(
    latent: np.ndarray,
    labels: pd.Series,
    reference_mask: np.ndarray,
    predictions: pd.Series,
    k: int = 25,
) -> pd.Series:
    missing = predictions.isna().to_numpy()
    if not missing.any():
        return predictions

    reference_labels = clean_labels(labels)
    training = (
        np.asarray(reference_mask, dtype=bool)
        & reference_labels.notna().to_numpy()
    )
    if not training.any():
        raise ValueError("No labeled reference cells are available for kNN.")

    classifier = KNeighborsClassifier(
        n_neighbors=min(k, int(training.sum())),
        weights="uniform",
        metric="euclidean",
    )
    classifier.fit(latent[training], reference_labels[training].astype(str))
    output = predictions.copy()
    output.iloc[np.flatnonzero(missing)] = classifier.predict(latent[missing])
    return output


def run_grid(
    adata: ad.AnnData,
    labels: pd.Series,
    donors: pd.Series,
    prefix: str,
    target_clusters: int,
    query_donor: str,
    seed: int,
) -> tuple[pd.Series, pd.Series, pd.DataFrame, pd.DataFrame, float]:
    metrics_records: list[dict[str, object]] = []
    stats_by_resolution: dict[float, pd.DataFrame] = {}

    for resolution in RESOLUTIONS:
        key = f"{prefix}_res_{str(resolution).replace('.', '_')}"
        sc.tl.leiden(
            adata,
            resolution=resolution,
            key_added=key,
            neighbors_key=f"{prefix}_neighbors",
            random_state=seed,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
        metrics, cluster_stats = evaluate_resolution(
            adata.obs[key],
            labels,
            donors,
            query_donor=query_donor,
            target_clusters=target_clusters,
        )
        metrics_records.append(
            {"resolution": resolution, "cluster_col": key, **metrics}
        )
        stats_by_resolution[resolution] = cluster_stats

    metrics_frame = pd.DataFrame.from_records(metrics_records)
    candidates = metrics_frame[metrics_frame["eligible"]].copy()
    if candidates.empty:
        candidates = metrics_frame.copy()
    candidates["cluster_distance"] = (
        candidates["n_clusters"] - target_clusters
    ).abs()
    selected_row = candidates.sort_values(
        ["selection_score", "cluster_distance", "resolution"],
        ascending=[False, True, True],
    ).iloc[0]
    selected_resolution = float(selected_row["resolution"])
    selected_key = str(selected_row["cluster_col"])
    selected_clusters = adata.obs[selected_key].astype(str)
    reference_mask = (
        donors.astype(str).ne(query_donor) & clean_labels(labels).notna()
    ).to_numpy()
    predictions, cluster_stats = majority_cluster_map(
        selected_clusters,
        labels,
        reference_mask,
        donors,
    )
    predictions = label_missing_by_knn(
        np.asarray(adata.obsm["X_resolvi"]),
        labels,
        reference_mask,
        predictions,
    )
    return (
        predictions,
        selected_clusters,
        metrics_frame,
        cluster_stats,
        selected_resolution,
    )


def run_lineage(
    full: ad.AnnData,
    major_predictions: pd.Series,
    reference_subtypes: pd.Series,
    lineage: str,
    subtype_levels: list[str],
    target_clusters: int,
    query_donor: str,
    seed: int,
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.Series]:
    keep = major_predictions.eq(lineage).to_numpy()
    lineage_data = full[keep].copy()
    lineage_data.obs["reference_subtype"] = clean_labels(
        reference_subtypes[keep]
    ).to_numpy()
    lineage_labels = clean_labels(lineage_data.obs["reference_subtype"]).mask(
        ~clean_labels(lineage_data.obs["reference_subtype"]).isin(
            subtype_levels
        )
    )
    donors = lineage_data.obs["SampleId"].astype(str)

    prefix = "fibroblast" if lineage == "fibroblast" else "endothelial"
    sc.pp.neighbors(
        lineage_data,
        n_neighbors=10,
        use_rep="X_resolvi",
        key_added=f"{prefix}_neighbors",
        random_state=seed,
    )
    sc.tl.umap(
        lineage_data,
        neighbors_key=f"{prefix}_neighbors",
        random_state=seed,
        min_dist=0.3,
    )
    (
        subtype_predictions,
        selected_clusters,
        metrics,
        cluster_stats,
        selected_resolution,
    ) = run_grid(
        lineage_data,
        lineage_labels,
        donors,
        prefix=prefix,
        target_clusters=target_clusters,
        query_donor=query_donor,
        seed=seed,
    )
    subtype_predictions = subtype_predictions.mask(
        ~subtype_predictions.isin(subtype_levels)
    )
    result = pd.DataFrame(
        {
            "cell": lineage_data.obs_names,
            "SampleId": donors.to_numpy(),
            "reference_subtype": lineage_labels.to_numpy(),
            "cluster_majority_subtype": subtype_predictions.to_numpy(),
            "joint_cluster": selected_clusters.to_numpy(),
            "UMAP1": lineage_data.obsm["X_umap"][:, 0],
            "UMAP2": lineage_data.obsm["X_umap"][:, 1],
        }
    )
    metrics.insert(0, "selected", metrics["resolution"].eq(selected_resolution))
    metrics.to_csv(
        output_dir / f"{prefix}_joint_resolution_grid.csv", index=False
    )
    cluster_stats.to_csv(
        output_dir / f"{prefix}_joint_cluster_label_map.csv", index=False
    )
    result.to_csv(
        output_dir / f"{prefix}_joint_labels_and_umap.csv", index=False
    )
    full_predictions = pd.Series(
        pd.NA, index=full.obs_names, dtype="string"
    )
    full_predictions.loc[lineage_data.obs_names] = subtype_predictions.to_numpy()
    return result, full_predictions


def metric_record(
    truth: np.ndarray,
    predicted: np.ndarray,
) -> dict[str, float]:
    return {
        "accuracy": float(accuracy_score(truth, predicted)),
        "balanced_accuracy": float(balanced_accuracy_score(truth, predicted)),
        "macro_f1": float(f1_score(truth, predicted, average="macro")),
    }


def select_lodo_knn(
    latent: np.ndarray,
    labels: pd.Series,
    donors: pd.Series,
    reference_mask: np.ndarray,
    allowed_levels: list[str],
    prefix: str,
    output_dir: Path,
) -> tuple[dict[str, object], pd.Series]:
    label_values = clean_labels(labels)
    donor_values = donors.astype(str)
    eligible = (
        np.asarray(reference_mask, dtype=bool)
        & label_values.isin(allowed_levels).to_numpy()
    )
    reference_donors = sorted(donor_values[eligible].unique())
    if len(reference_donors) < 3:
        raise ValueError(f"{prefix}: fewer than three reference donors.")

    grid_records: list[dict[str, object]] = []
    predictions_by_config: dict[tuple[str, int], pd.Series] = {}
    folds_by_config: dict[tuple[str, int], pd.DataFrame] = {}
    for weights in KNN_WEIGHTS:
        for k in KNN_NEIGHBORS:
            cross_validated = pd.Series(
                pd.NA, index=labels.index, dtype="string"
            )
            fold_records: list[dict[str, object]] = []
            for held_out_donor in reference_donors:
                testing = eligible & donor_values.eq(held_out_donor).to_numpy()
                training = eligible & donor_values.ne(held_out_donor).to_numpy()
                classifier = KNeighborsClassifier(
                    n_neighbors=min(k, int(training.sum())),
                    weights=weights,
                    metric="euclidean",
                    n_jobs=-1,
                )
                classifier.fit(
                    latent[training],
                    label_values[training].astype(str),
                )
                fold_prediction = classifier.predict(latent[testing])
                cross_validated.iloc[np.flatnonzero(testing)] = fold_prediction
                fold_truth = label_values[testing].astype(str).to_numpy()
                fold_records.append(
                    {
                        "held_out_donor": held_out_donor,
                        "n_cells": int(testing.sum()),
                        **metric_record(fold_truth, fold_prediction),
                    }
                )

            truth = label_values[eligible].astype(str).to_numpy()
            predicted = cross_validated[eligible].astype(str).to_numpy()
            key = (weights, k)
            grid_records.append(
                {
                    "weights": weights,
                    "k": k,
                    "n_cells": int(eligible.sum()),
                    "n_donors": len(reference_donors),
                    **metric_record(truth, predicted),
                }
            )
            predictions_by_config[key] = cross_validated
            folds_by_config[key] = pd.DataFrame.from_records(fold_records)

    grid = pd.DataFrame.from_records(grid_records)
    selected_row = grid.sort_values(
        ["balanced_accuracy", "macro_f1", "accuracy", "k", "weights"],
        ascending=[False, False, False, True, True],
    ).iloc[0]
    selected_key = (str(selected_row["weights"]), int(selected_row["k"]))
    grid.insert(
        0,
        "selected",
        grid["weights"].eq(selected_key[0]) & grid["k"].eq(selected_key[1]),
    )
    grid.to_csv(output_dir / f"{prefix}_knn_lodo_grid.csv", index=False)
    folds_by_config[selected_key].to_csv(
        output_dir / f"{prefix}_knn_lodo_fold_metrics.csv",
        index=False,
    )

    selected_predictions = predictions_by_config[selected_key]
    prediction_table = pd.DataFrame(
        {
            "cell": labels.index[eligible],
            "SampleId": donor_values[eligible].to_numpy(),
            "observed_label": label_values[eligible].to_numpy(),
            "lodo_predicted_label": selected_predictions[eligible].to_numpy(),
        }
    )
    prediction_table.to_csv(
        output_dir / f"{prefix}_knn_lodo_predictions.csv",
        index=False,
    )
    confusion = (
        prediction_table.groupby(
            ["observed_label", "lodo_predicted_label"],
            dropna=False,
        )
        .size()
        .rename("n_cells")
        .reset_index()
    )
    confusion.to_csv(
        output_dir / f"{prefix}_knn_lodo_confusion.csv",
        index=False,
    )
    selected = {
        "weights": selected_key[0],
        "k": selected_key[1],
        "n_cells": int(selected_row["n_cells"]),
        "n_donors": int(selected_row["n_donors"]),
        "accuracy": float(selected_row["accuracy"]),
        "balanced_accuracy": float(selected_row["balanced_accuracy"]),
        "macro_f1": float(selected_row["macro_f1"]),
    }
    return selected, selected_predictions


def predict_query_knn(
    latent: np.ndarray,
    labels: pd.Series,
    reference_mask: np.ndarray,
    query_target_mask: np.ndarray,
    allowed_levels: list[str],
    config: dict[str, object],
) -> tuple[pd.Series, pd.Series]:
    label_values = clean_labels(labels)
    training = (
        np.asarray(reference_mask, dtype=bool)
        & label_values.isin(allowed_levels).to_numpy()
    )
    target = np.asarray(query_target_mask, dtype=bool)
    classifier = KNeighborsClassifier(
        n_neighbors=min(int(config["k"]), int(training.sum())),
        weights=str(config["weights"]),
        metric="euclidean",
        n_jobs=-1,
    )
    classifier.fit(latent[training], label_values[training].astype(str))
    predictions = pd.Series(pd.NA, index=labels.index, dtype="string")
    confidence = pd.Series(np.nan, index=labels.index, dtype=float)
    if target.any():
        probabilities = classifier.predict_proba(latent[target])
        predictions.iloc[np.flatnonzero(target)] = classifier.classes_[
            probabilities.argmax(axis=1)
        ]
        confidence.iloc[np.flatnonzero(target)] = probabilities.max(axis=1)
    return predictions, confidence


def agreement(
    predicted: pd.Series,
    reference: pd.Series,
    mask: np.ndarray,
    missing_as_mismatch: bool = False,
) -> tuple[float, int]:
    reference_values = clean_labels(reference)
    comparable = np.asarray(mask, dtype=bool) & reference_values.notna().to_numpy()
    if not missing_as_mismatch:
        comparable &= predicted.notna().to_numpy()
    if not comparable.any():
        return np.nan, 0
    predicted_values = predicted.fillna("__UNASSIGNED__").astype(str)
    value = float(
        np.mean(
            predicted_values[comparable].to_numpy()
            == reference_values[comparable].astype(str).to_numpy()
        )
    )
    return value, int(comparable.sum())


def main() -> None:
    args = parse_args()
    model_h5ad = args.model_h5ad.resolve()
    output_dir = args.output_dir.resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    if not model_h5ad.is_file():
        raise FileNotFoundError(model_h5ad)

    print(f"Reading model metadata and latent: {model_h5ad}", flush=True)
    source = ad.read_h5ad(model_h5ad, backed="r")
    required_obs = {"SampleId", "lv1_anno_final", "lv2_anno_final"}
    missing_obs = sorted(required_obs - set(source.obs.columns))
    if missing_obs:
        raise KeyError(f"Missing model obs columns: {missing_obs}")
    if "X_resolvi" not in source.obsm:
        raise KeyError("X_resolvi is missing from the model H5AD.")

    obs = source.obs.copy()
    latent = np.asarray(source.obsm["X_resolvi"]).astype(np.float32)
    obs_names = source.obs_names.astype(str).copy()
    source.file.close()
    if latent.shape != (len(obs_names), 30) or not np.isfinite(latent).all():
        raise ValueError(f"Unexpected ResolVI latent shape/content: {latent.shape}")

    joint = ad.AnnData(obs=obs)
    joint.obs_names = obs_names
    joint.obsm["X_resolvi"] = latent
    donors = joint.obs["SampleId"].astype(str)
    reference_major = clean_labels(joint.obs["lv1_anno_final"])
    reference_subtype = clean_labels(joint.obs["lv2_anno_final"])
    reference_mask = donors.ne(args.query_donor).to_numpy()
    query_mask = donors.eq(args.query_donor).to_numpy()

    print("Computing joint ResolVI neighbors, UMAP, and major clusters", flush=True)
    sc.pp.neighbors(
        joint,
        n_neighbors=20,
        use_rep="X_resolvi",
        key_added="major_neighbors",
        random_state=args.seed,
    )
    sc.tl.umap(
        joint,
        neighbors_key="major_neighbors",
        random_state=args.seed,
        min_dist=0.3,
    )
    (
        major_cluster_predictions,
        major_clusters,
        major_metrics,
        major_cluster_stats,
        major_resolution,
    ) = run_grid(
        joint,
        reference_major,
        donors,
        prefix="major",
        target_clusters=20,
        query_donor=args.query_donor,
        seed=args.seed,
    )
    major_cluster_predictions = major_cluster_predictions.mask(
        ~major_cluster_predictions.isin(MAJOR_LEVELS)
    )

    major_metrics.insert(
        0, "selected", major_metrics["resolution"].eq(major_resolution)
    )
    major_metrics.to_csv(
        output_dir / "major_joint_resolution_grid.csv", index=False
    )
    major_cluster_stats.to_csv(
        output_dir / "major_joint_cluster_label_map.csv", index=False
    )

    print("Tuning sample-aware major-label mapping in the joint latent", flush=True)
    major_knn_config, major_lodo_predictions = select_lodo_knn(
        latent,
        reference_major,
        donors,
        reference_mask,
        MAJOR_LEVELS,
        prefix="major",
        output_dir=output_dir,
    )
    major_query_predictions, major_query_confidence = predict_query_knn(
        latent,
        reference_major,
        reference_mask,
        query_mask,
        MAJOR_LEVELS,
        major_knn_config,
    )
    final_major = reference_major.copy()
    final_major.iloc[np.flatnonzero(query_mask)] = (
        major_query_predictions[query_mask].to_numpy()
    )
    final_major = final_major.mask(~final_major.isin(MAJOR_LEVELS))
    if final_major.isna().any():
        raise ValueError("Final major annotation contains unresolved cells.")

    major_result = pd.DataFrame(
        {
            "cell": joint.obs_names,
            "SampleId": donors.to_numpy(),
            "reference_lv1": reference_major.to_numpy(),
            "joint_lv1": final_major.to_numpy(),
            "cluster_majority_lv1": major_cluster_predictions.to_numpy(),
            "mapping_confidence": major_query_confidence.to_numpy(),
            "joint_cluster": major_clusters.to_numpy(),
            "UMAP1": joint.obsm["X_umap"][:, 0],
            "UMAP2": joint.obsm["X_umap"][:, 1],
        }
    )
    major_result.to_csv(
        output_dir / "major_joint_labels_and_umap.csv", index=False
    )

    print("Computing joint ResolVI fibroblast subclusters", flush=True)
    fibroblast_result, fibroblast_cluster_predictions = run_lineage(
        joint,
        final_major,
        reference_subtype,
        lineage="fibroblast",
        subtype_levels=FIB_LEVELS,
        target_clusters=9,
        query_donor=args.query_donor,
        seed=args.seed,
        output_dir=output_dir,
    )
    print("Computing joint ResolVI endothelial subclusters", flush=True)
    endothelial_result, endothelial_cluster_predictions = run_lineage(
        joint,
        final_major,
        reference_subtype,
        lineage="endothelial",
        subtype_levels=EC_LEVELS,
        target_clusters=7,
        query_donor=args.query_donor,
        seed=args.seed,
        output_dir=output_dir,
    )

    print(
        "Tuning sample-aware fibroblast and endothelial subtype mappings",
        flush=True,
    )
    fib_knn_config, fib_lodo_predictions = select_lodo_knn(
        latent,
        reference_subtype,
        donors,
        reference_mask,
        FIB_LEVELS,
        prefix="fibroblast_subtype",
        output_dir=output_dir,
    )
    ec_knn_config, ec_lodo_predictions = select_lodo_knn(
        latent,
        reference_subtype,
        donors,
        reference_mask,
        EC_LEVELS,
        prefix="endothelial_subtype",
        output_dir=output_dir,
    )
    fib_query_mask = query_mask & final_major.eq("fibroblast").to_numpy()
    ec_query_mask = query_mask & final_major.eq("endothelial").to_numpy()
    fib_query_predictions, fib_query_confidence = predict_query_knn(
        latent,
        reference_subtype,
        reference_mask,
        fib_query_mask,
        FIB_LEVELS,
        fib_knn_config,
    )
    ec_query_predictions, ec_query_confidence = predict_query_knn(
        latent,
        reference_subtype,
        reference_mask,
        ec_query_mask,
        EC_LEVELS,
        ec_knn_config,
    )

    final_subtype = pd.Series(pd.NA, index=joint.obs_names, dtype="string")
    valid_reference_subtype = (
        reference_mask
        & reference_subtype.isin(FIB_LEVELS + EC_LEVELS).to_numpy()
    )
    final_subtype.iloc[np.flatnonzero(valid_reference_subtype)] = (
        reference_subtype[valid_reference_subtype].to_numpy()
    )
    final_subtype.iloc[np.flatnonzero(fib_query_mask)] = (
        fib_query_predictions[fib_query_mask].to_numpy()
    )
    final_subtype.iloc[np.flatnonzero(ec_query_mask)] = (
        ec_query_predictions[ec_query_mask].to_numpy()
    )
    subtype_query_confidence = fib_query_confidence.combine_first(
        ec_query_confidence
    )
    cluster_majority_subtype = fibroblast_cluster_predictions.combine_first(
        endothelial_cluster_predictions
    )

    major_cv_metrics = metric_record(
        reference_major[reference_mask].astype(str).to_numpy(),
        major_lodo_predictions[reference_mask].astype(str).to_numpy(),
    )
    fib_cv_mask = reference_mask & reference_subtype.isin(FIB_LEVELS).to_numpy()
    fib_cv_metrics = metric_record(
        reference_subtype[fib_cv_mask].astype(str).to_numpy(),
        fib_lodo_predictions[fib_cv_mask].astype(str).to_numpy(),
    )
    ec_cv_mask = reference_mask & reference_subtype.isin(EC_LEVELS).to_numpy()
    ec_cv_metrics = metric_record(
        reference_subtype[ec_cv_mask].astype(str).to_numpy(),
        ec_lodo_predictions[ec_cv_mask].astype(str).to_numpy(),
    )
    major_query_agreement, major_query_n = agreement(
        final_major, reference_major, query_mask
    )
    subtype_query_agreement, subtype_query_n = agreement(
        final_subtype,
        reference_subtype,
        query_mask,
        missing_as_mismatch=True,
    )
    major_cluster_agreement, major_cluster_n = agreement(
        major_cluster_predictions,
        reference_major,
        reference_mask,
        missing_as_mismatch=True,
    )

    summary = pd.DataFrame(
        [
            {
                "metric": "Six-donor leave-one-donor-out major-label accuracy in joint ResolVI latent",
                "value": major_cv_metrics["accuracy"],
                "unit": "proportion",
                "n_cells": int(reference_mask.sum()),
            },
            {
                "metric": "Six-donor leave-one-donor-out major-label balanced accuracy in joint ResolVI latent",
                "value": major_cv_metrics["balanced_accuracy"],
                "unit": "proportion",
                "n_cells": int(reference_mask.sum()),
            },
            {
                "metric": "Six-donor leave-one-donor-out fibroblast-subtype accuracy in joint ResolVI latent",
                "value": fib_cv_metrics["accuracy"],
                "unit": "proportion",
                "n_cells": int(fib_cv_mask.sum()),
            },
            {
                "metric": "Six-donor leave-one-donor-out endothelial-subtype accuracy in joint ResolVI latent",
                "value": ec_cv_metrics["accuracy"],
                "unit": "proportion",
                "n_cells": int(ec_cv_mask.sum()),
            },
            {
                "metric": "ExcludedValidation major-label agreement: joint ResolVI latent mapping versus prior RNA-anchor mapping",
                "value": major_query_agreement,
                "unit": "proportion",
                "n_cells": major_query_n,
            },
            {
                "metric": "ExcludedValidation Fib/EC subtype agreement: joint ResolVI latent mapping versus prior RNA-anchor mapping",
                "value": subtype_query_agreement,
                "unit": "proportion",
                "n_cells": subtype_query_n,
            },
            {
                "metric": "Six-donor major-label accuracy of cluster-majority sensitivity annotation",
                "value": major_cluster_agreement,
                "unit": "proportion",
                "n_cells": major_cluster_n,
            },
        ]
    )
    summary.to_csv(
        output_dir / "joint_resolvi_annotation_summary.csv", index=False
    )

    annotations = pd.DataFrame(
        {
            "cell": joint.obs_names,
            "SampleId": donors.to_numpy(),
            "reference_lv1": reference_major.to_numpy(),
            "reference_lv2": reference_subtype.to_numpy(),
            "joint_resolvi_lv1": final_major.to_numpy(),
            "joint_resolvi_lv2": final_subtype.to_numpy(),
            "joint_resolvi_lv1_mapping_confidence": (
                major_query_confidence.to_numpy()
            ),
            "joint_resolvi_lv2_mapping_confidence": (
                subtype_query_confidence.to_numpy()
            ),
            "cluster_majority_lv1_sensitivity": (
                major_cluster_predictions.to_numpy()
            ),
            "cluster_majority_lv2_sensitivity": (
                cluster_majority_subtype.to_numpy()
            ),
            "joint_resolvi_major_cluster": major_clusters.to_numpy(),
            "joint_resolvi_umap1": joint.obsm["X_umap"][:, 0],
            "joint_resolvi_umap2": joint.obsm["X_umap"][:, 1],
        }
    )
    annotations.to_csv(
        output_dir / "joint_resolvi_cell_annotations.csv.gz",
        index=False,
        compression="gzip",
    )
    latent_frame = pd.DataFrame(
        latent,
        columns=[f"X_resolvi_{index}" for index in range(1, 31)],
    )
    latent_frame.insert(0, "cell", joint.obs_names)
    latent_frame.to_csv(
        output_dir / "joint_resolvi_latent.csv.gz",
        index=False,
        compression="gzip",
    )

    manifest = {
        "model_h5ad": str(model_h5ad),
        "model_h5ad_sha256": sha256_file(model_h5ad),
        "n_cells": int(joint.n_obs),
        "n_latent": int(latent.shape[1]),
        "seed": args.seed,
        "final_annotation_strategy": (
            "preserve curated six-donor labels and map ExcludedValidation with "
            "leave-one-donor-out-tuned kNN in the joint ResolVI latent"
        ),
        "query_donor_excluded_from_classifier_training": args.query_donor,
        "six_donor_reference_labels_preserved": True,
        "major_knn": major_knn_config,
        "fibroblast_subtype_knn": fib_knn_config,
        "endothelial_subtype_knn": ec_knn_config,
        "joint_leiden_role": "sensitivity analysis and visualization only",
        "major_selected_resolution": major_resolution,
        "fibroblast_n_cells": int(len(fibroblast_result)),
        "endothelial_n_cells": int(len(endothelial_result)),
        "versions": {
            "anndata": importlib.metadata.version("anndata"),
            "scanpy": importlib.metadata.version("scanpy"),
            "scikit_learn": importlib.metadata.version("scikit-learn"),
            "igraph": importlib.metadata.version("igraph"),
            "leidenalg": importlib.metadata.version("leidenalg"),
        },
    }
    (output_dir / "joint_resolvi_annotation_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    print(summary.to_string(index=False), flush=True)
    print(f"Completed annotation outputs: {output_dir}", flush=True)


if __name__ == "__main__":
    main()
