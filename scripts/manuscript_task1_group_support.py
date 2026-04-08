#!/usr/bin/env python3
"""
Shared Task1 group-metric helpers for the live A1 builder and retained
historical bridge builders.

Status:
- support-only helper

Consumed by:
- manuscript support builders; no standalone canonical manuscript output role

Architecture:
- helper surface only; see `scripts/ARCHITECTURE.md`

These helpers keep the manuscript bridge scripts aligned with the existing
Task1 contracts without reintroducing representation as a shared-group key.
"""

from __future__ import annotations

from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

try:
    from s1_task1_internal_metrics import (
        EXPECTED_TASK1_SNAPSHOT,
        GLOBAL_SEED,
        ProjectedPathwayStore,
        cosine_similarity,
        deterministic_subsample_by_uid,
        energy_distance_biascorr,
        init_global_rng,
        load_lincs_inputs,
        load_scperturb_inputs,
        make_group_id,
        normalize_perturbation_type,
        normalize_text,
        pearson_corr,
    )
except ModuleNotFoundError:
    from scripts.s1_task1_internal_metrics import (
        EXPECTED_TASK1_SNAPSHOT,
        GLOBAL_SEED,
        ProjectedPathwayStore,
        cosine_similarity,
        deterministic_subsample_by_uid,
        energy_distance_biascorr,
        init_global_rng,
        load_lincs_inputs,
        load_scperturb_inputs,
        make_group_id,
        normalize_perturbation_type,
        normalize_text,
        pearson_corr,
    )


COMMON_REPRESENTATIONS = {"Gene", "Pathway"}
TRIPLET_COLUMNS = ["dataset", "cell_line", "target"]
DETAIL_COLUMNS = TRIPLET_COLUMNS + ["representation_detail"]
GROUP_METRICS = ["cosine", "pcc", "edist"]


def paired_datasets_from_cross_direction(dataset_or_direction: str) -> list[str]:
    if not isinstance(dataset_or_direction, str) or "_to_" not in dataset_or_direction:
        raise ValueError(f"Unexpected Task1 cross direction: {dataset_or_direction}")
    source_dataset, target_dataset = dataset_or_direction.split("_to_", 1)
    return [source_dataset, target_dataset]


def source_dataset_from_cross_direction(dataset_or_direction: str) -> str | None:
    if not isinstance(dataset_or_direction, str) or "_to_" not in dataset_or_direction:
        return None
    return dataset_or_direction.split("_to_", 1)[0]


def _iter_internal_group_rows(
    project_root: Path,
    allowed_triplets: pd.DataFrame | None = None,
) -> Iterable[dict[str, object]]:
    snapshot_root = (project_root / EXPECTED_TASK1_SNAPSHOT).resolve()
    rng = init_global_rng(GLOBAL_SEED)
    lincs_meta, lincs_gene_store, pathway_w, _, _ = load_lincs_inputs(snapshot_root)
    sc_meta_map, sc_gene_map, sc_pathway_map, _ = load_scperturb_inputs(snapshot_root)

    lincs_pert = lincs_meta["pert_type"].map(normalize_perturbation_type)
    lincs_base = pd.DataFrame(
        {
            "cell_line": normalize_text(lincs_meta["cell_line"]),
            "target": normalize_text(lincs_meta["target"]),
            "canonical_query_uid": np.arange(len(lincs_meta), dtype=np.int64).astype(str),
            "delta_row_index": np.arange(len(lincs_meta), dtype=np.int64),
            "delta_valid_bool": True,
            "perturbation_type": lincs_pert,
        }
    )
    lincs_base = lincs_base.loc[lincs_base["perturbation_type"].eq("Genetic")].copy()
    lincs_meta_min = lincs_base[
        ["cell_line", "target", "canonical_query_uid", "delta_row_index", "delta_valid_bool"]
    ].reset_index(drop=True)
    lincs_pathway_store = ProjectedPathwayStore(lincs_gene_store, pathway_w)

    sc_src = sc_meta_map["Genetic"]
    sc_meta_min = pd.DataFrame(
        {
            "cell_line": normalize_text(sc_src["cell_std"]),
            "target": normalize_text(sc_src["target_std"]),
            "canonical_query_uid": sc_src["delta_row_idx"].astype(np.int64).astype(str),
            "delta_row_index": sc_src["delta_row_idx"].astype(np.int64),
            "delta_valid_bool": True,
        }
    )

    cohorts = [
        ("LINCS", "Gene", lincs_meta_min.copy(), lincs_gene_store),
        ("LINCS", "Pathway", lincs_meta_min.copy(), lincs_pathway_store),
        ("scPerturb", "Gene", sc_meta_min.copy(), sc_gene_map["Genetic"]),
        ("scPerturb", "Pathway", sc_meta_min.copy(), sc_pathway_map["Genetic"]),
    ]

    allowed_by_dataset: dict[str, pd.DataFrame] = {}
    if allowed_triplets is not None and not allowed_triplets.empty:
        allowed_frame = allowed_triplets[TRIPLET_COLUMNS].drop_duplicates().copy()
        for dataset, subset in allowed_frame.groupby("dataset", sort=False):
            allowed_by_dataset[str(dataset)] = subset[["cell_line", "target"]].drop_duplicates().copy()

    for dataset, representation_detail, meta, vector_store in cohorts:
        if allowed_by_dataset:
            allowed_subset = allowed_by_dataset.get(dataset)
            if allowed_subset is None or allowed_subset.empty:
                continue
            meta = meta.merge(allowed_subset, on=["cell_line", "target"], how="inner")
            if meta.empty:
                continue

        meta["delta_valid_bool"] = meta["delta_valid_bool"].fillna(False).astype(bool)
        meta_work = meta.loc[meta["delta_valid_bool"]].copy().reset_index(drop=True)
        if meta_work.empty:
            continue

        meta_work["group_id"] = [
            make_group_id(cell_line=row.cell_line, target=row.target)
            for row in meta_work.itertuples(index=False)
        ]

        group_to_positions = {
            group_id: pd.Index(sorted(idx)).to_numpy(dtype="int64")
            for group_id, idx in meta_work.groupby("group_id", sort=False).groups.items()
        }

        for group_id in sorted(group_to_positions):
            pos_all = group_to_positions[group_id]
            row_idx = meta_work.loc[pos_all, "delta_row_index"].to_numpy(dtype="int64")
            vectors = vector_store.load_rows(row_idx)
            finite_mask = np.isfinite(vectors).all(axis=1)
            if not finite_mask.any():
                continue

            pos_valid = pos_all[finite_mask]
            vec_valid = vectors[finite_mask]
            uid_valid = meta_work.loc[pos_valid, "canonical_query_uid"].astype(str).to_numpy()
            first_row = meta_work.loc[pos_valid[0]]
            n_group = int(vec_valid.shape[0])

            row: dict[str, object] = {
                "dataset": dataset,
                "cell_line": str(first_row["cell_line"]),
                "target": str(first_row["target"]),
                "representation_detail": representation_detail,
                "group_id": group_id,
                "n_total": n_group,
                "n_A": pd.NA,
                "n_B": pd.NA,
                "n_A_sub": pd.NA,
                "n_B_sub": pd.NA,
                "cosine": pd.NA,
                "pcc": pd.NA,
                "edist": pd.NA,
            }

            if n_group < 4:
                yield row
                continue

            sort_idx = uid_valid.argsort(kind="mergesort")
            vec_sorted = vec_valid[sort_idx]
            uid_sorted = uid_valid[sort_idx]

            perm = rng.permutation(n_group)
            half = n_group // 2
            a_idx = perm[:half]
            b_idx = perm[half:]

            a_vec = vec_sorted[a_idx]
            b_vec = vec_sorted[b_idx]
            a_uid = uid_sorted[a_idx]
            b_uid = uid_sorted[b_idx]

            centroid_a = a_vec.mean(axis=0)
            centroid_b = b_vec.mean(axis=0)

            row["n_A"] = int(len(a_vec))
            row["n_B"] = int(len(b_vec))
            row["cosine"] = float(cosine_similarity(centroid_a, centroid_b))
            row["pcc"] = float(pearson_corr(centroid_a, centroid_b))

            if len(a_vec) >= 2 and len(b_vec) >= 2:
                a_sub, n_a_sub = deterministic_subsample_by_uid(a_vec, a_uid, 256, rng)
                b_sub, n_b_sub = deterministic_subsample_by_uid(b_vec, b_uid, 256, rng)
                row["n_A_sub"] = int(n_a_sub)
                row["n_B_sub"] = int(n_b_sub)
                if len(a_sub) >= 2 and len(b_sub) >= 2:
                    row["edist"] = float(energy_distance_biascorr(a_sub, b_sub))
            yield row


def compute_task1_internal_group_metrics(
    project_root: Path,
    allowed_triplets: pd.DataFrame | None = None,
) -> pd.DataFrame:
    frame = pd.DataFrame(_iter_internal_group_rows(project_root, allowed_triplets=allowed_triplets))
    if frame.empty:
        return pd.DataFrame(
            columns=[
                *DETAIL_COLUMNS,
                "group_id",
                "n_total",
                "n_A",
                "n_B",
                "n_A_sub",
                "n_B_sub",
                *GROUP_METRICS,
            ]
        )
    return frame.sort_values(DETAIL_COLUMNS, kind="mergesort").reset_index(drop=True)


def load_or_compute_task1_internal_group_metrics(
    project_root: Path,
    cache_path: Path,
    allowed_triplets: pd.DataFrame | None = None,
) -> pd.DataFrame:
    if cache_path.exists():
        return pd.read_parquet(cache_path)

    frame = compute_task1_internal_group_metrics(project_root, allowed_triplets=allowed_triplets)
    cache_path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_parquet(cache_path, index=False)
    return frame


def load_task1_cross_group_metrics(path: Path) -> pd.DataFrame:
    frame = pd.read_parquet(
        path,
        columns=[
            "dataset_or_direction",
            "perturbation_type",
            "representation",
            "cell_line",
            "target_token",
            "group_id",
            "n_L",
            "n_S",
            "n_L_sub",
            "n_S_sub",
            "cosine",
            "pcc",
            "edist",
        ],
    )
    frame = frame.loc[
        frame["perturbation_type"].eq("Genetic")
        & frame["representation"].isin(COMMON_REPRESENTATIONS)
    ].copy()

    rows: list[dict[str, object]] = []
    for row in frame.itertuples(index=False):
        for dataset in paired_datasets_from_cross_direction(row.dataset_or_direction):
            rows.append(
                {
                    "dataset": dataset,
                    "cell_line": row.cell_line,
                    "target": row.target_token,
                    "representation_detail": row.representation,
                    "group_id": row.group_id,
                    "task1_pair_direction": row.dataset_or_direction,
                    "n_L": row.n_L,
                    "n_S": row.n_S,
                    "n_L_sub": row.n_L_sub,
                    "n_S_sub": row.n_S_sub,
                    "cosine": row.cosine,
                    "pcc": row.pcc,
                    "edist": row.edist,
                }
            )

    if not rows:
        return pd.DataFrame(
            columns=[
                *DETAIL_COLUMNS,
                "group_id",
                "task1_pair_direction",
                "n_L",
                "n_S",
                "n_L_sub",
                "n_S_sub",
                *GROUP_METRICS,
            ]
        )

    return pd.DataFrame(rows).sort_values(DETAIL_COLUMNS, kind="mergesort").reset_index(drop=True)
