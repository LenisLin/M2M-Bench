# SCRIPT_HEADER_CONTRACT
# Script: scripts/build_task1_snapshot_v2_candidate.py
# Purpose: Build a candidate Task1 snapshot directly from upstream LINCS/scPerturb raw sources.
# Outputs:
#   - <output-root>/lincs/lincs-engine1-meta.csv
#   - <output-root>/lincs/lincs-engine1-gene-delta.pt
#   - <output-root>/lincs/lincs-gene-alignment.csv
#   - <output-root>/scperturb_delta/scperturb-chemical-delta-meta.csv
#   - <output-root>/scperturb_delta/scperturb-chemical-gene-delta.npy
#   - <output-root>/scperturb_delta/scperturb-chemical-pathway-delta.npy
#   - <output-root>/scperturb_delta/scperturb-crispr-delta-meta.csv
#   - <output-root>/scperturb_delta/scperturb-crispr-gene-delta.npy
#   - <output-root>/scperturb_delta/scperturb-crispr-pathway-delta.npy
#   - <output-root>/pathway/{hallmark-w-2477x50.npy,lincs-pathway-policy.json}
#   - <output-root>/fm_delta/staging_manifest.json
#   - <output-root>/snapshot_manifest.json

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import shutil
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Iterable, Mapping, Sequence

import h5py
import numpy as np
import pandas as pd
import torch
import yaml

try:
    from m2mbench.utils.raw_rebuild import (
        ScPerturbSource,
        build_alignment_index,
        build_scperturb_target_fields,
        build_standardization_audit_rows,
        choose_scperturb_cell_column,
        dataframe_from_records,
        discover_scperturb_sources,
        is_control_like,
        join_tokens,
        normalize_lincs_subtype,
        normalize_optional_text,
        normalize_required_text,
        normalize_scperturb_subtype,
        perturbation_class_for_subtype,
        standardize_cell_line,
    )
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from m2mbench.utils.raw_rebuild import (
        ScPerturbSource,
        build_alignment_index,
        build_scperturb_target_fields,
        build_standardization_audit_rows,
        choose_scperturb_cell_column,
        dataframe_from_records,
        discover_scperturb_sources,
        is_control_like,
        join_tokens,
        normalize_lincs_subtype,
        normalize_optional_text,
        normalize_required_text,
        normalize_scperturb_subtype,
        perturbation_class_for_subtype,
        standardize_cell_line,
    )

try:
    import anndata as ad  # type: ignore
except Exception:  # pragma: no cover - optional runtime import
    ad = None

try:
    from path_policy import (
        DEFAULT_OSMOSIS_LINCS_ROOT,
        DEFAULT_OSMOSIS_SCPERTURB_ROOT,
        DEFAULT_TASK1_SNAPSHOT_ROOT,
        DEFAULT_TASK1_SNAPSHOT_V2_CANDIDATE_ROOT,
    )
except ModuleNotFoundError:
    from scripts.path_policy import (
        DEFAULT_OSMOSIS_LINCS_ROOT,
        DEFAULT_OSMOSIS_SCPERTURB_ROOT,
        DEFAULT_TASK1_SNAPSHOT_ROOT,
        DEFAULT_TASK1_SNAPSHOT_V2_CANDIDATE_ROOT,
    )

CONFIG_PATH = Path("config/config.yaml")
GLOBAL_SEED = 619
MAX_CONTROLS = 50


@dataclass(frozen=True)
class ScPerturbBuildResult:
    meta: pd.DataFrame
    gene_delta: np.ndarray
    pathway_delta: np.ndarray
    schema_audit: pd.DataFrame
    attrition: pd.DataFrame


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build raw-driven Task1 snapshot v2 candidate")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--output-root", type=Path, default=DEFAULT_TASK1_SNAPSHOT_V2_CANDIDATE_ROOT)
    parser.add_argument("--seed", type=int, default=GLOBAL_SEED)
    return parser.parse_args()


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def resolve_config_path(project_root: Path, raw_path: str | Path) -> Path:
    path = Path(str(raw_path))
    if path.is_absolute():
        return path.resolve()
    return (project_root / path).resolve()


def compute_sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8") as handle:
        json.dump(payload, handle, indent=2, ensure_ascii=True)
        handle.write("\n")


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)


def copy_if_needed(src: Path, dst: Path) -> None:
    dst.parent.mkdir(parents=True, exist_ok=True)
    if dst.exists():
        return
    shutil.copy2(src, dst)


def ensure_required_columns(frame: pd.DataFrame, required: Sequence[str], name: str) -> None:
    missing = [column for column in required if column not in frame.columns]
    if missing:
        raise ValueError(f"Missing required columns for {name}: {missing}")


def load_task1_sidecar_assets(sidecar_root: Path) -> tuple[pd.DataFrame, np.ndarray, Path, Path]:
    lincs_meta = pd.read_csv(sidecar_root / "lincs/lincs-engine1-meta.csv")
    alignment_path = sidecar_root / "lincs/lincs-gene-alignment.csv"
    pathway_w_path = sidecar_root / "pathway/hallmark-w-2477x50.npy"
    pathway_policy_path = sidecar_root / "pathway/lincs-pathway-policy.json"
    alignment = pd.read_csv(alignment_path)
    ensure_required_columns(alignment, ["gene_symbol", "local_idx"], "lincs-gene-alignment.csv")
    alignment = alignment.sort_values("local_idx", kind="mergesort").reset_index(drop=True)
    pathway_w = np.load(pathway_w_path).astype(np.float32, copy=False)
    if pathway_w.shape != (2477, 50):
        raise ValueError(f"Unexpected pathway W shape: {pathway_w.shape}")
    return lincs_meta, pathway_w, alignment_path, pathway_policy_path


def build_lincs_matrix_row_lookup(raw_gctx_path: Path) -> tuple[h5py.File, h5py.Dataset, dict[str, int]]:
    handle = h5py.File(raw_gctx_path, "r")
    matrix = handle["0"]["DATA"]["0"]["matrix"]
    col_ids = handle["0"]["META"]["COL"]["id"][:]
    lookup = {
        raw_id.decode("utf-8"): int(idx)
        for idx, raw_id in enumerate(col_ids.tolist())
    }
    return handle, matrix, lookup


def build_lincs_gene_alignment_index(alignment: pd.DataFrame, raw_gene_info_path: Path) -> np.ndarray:
    raw_gene = pd.read_csv(raw_gene_info_path, sep="\t")
    ensure_required_columns(raw_gene, ["pr_gene_symbol"], "raw LINCS gene_info")
    return build_alignment_index(
        alignment["gene_symbol"].astype(str).tolist(),
        raw_gene["pr_gene_symbol"].astype(str).tolist(),
    )


def load_lincs_aligned_rows(
    *,
    matrix: h5py.Dataset,
    row_positions: np.ndarray,
    gene_index_map: np.ndarray,
) -> np.ndarray:
    out = np.zeros((int(row_positions.size), int(gene_index_map.size)), dtype=np.float32)
    valid_local = np.where(gene_index_map >= 0)[0].astype(np.int64, copy=False)
    if valid_local.size == 0 or row_positions.size == 0:
        return out
    valid_source = gene_index_map[valid_local].astype(np.int64, copy=False)
    order = np.argsort(row_positions, kind="mergesort")
    sorted_rows = row_positions[order].astype(np.int64, copy=False)
    block = np.asarray(matrix[sorted_rows.tolist(), :][:, valid_source.tolist()], dtype=np.float32)
    reverse = np.empty_like(order)
    reverse[order] = np.arange(order.size, dtype=np.int64)
    out[:, valid_local] = block[reverse]
    return out


def build_lincs_meta_and_gene_delta(
    *,
    raw_lincs_root: Path,
    sidecar_root: Path,
    output_root: Path,
) -> tuple[pd.DataFrame, Path, pd.DataFrame]:
    raw_gctx_path = raw_lincs_root / "GSE92742/GSE92742_Broad_LINCS_Level5_COMPZ.MODZ_n473647x12328.gctx"
    raw_gene_info_path = raw_lincs_root / "GSE92742/GSE92742_Broad_LINCS_gene_info.txt.gz"
    raw_sig_info_path = raw_lincs_root / "GSE92742/GSE92742_Broad_LINCS_sig_info.txt.gz"
    beta_siginfo_path = raw_lincs_root / "beta/siginfo_beta.txt"

    sidecar_meta, _pathway_w, alignment_path, pathway_policy_path = load_task1_sidecar_assets(sidecar_root)
    raw_sig_info = pd.read_csv(raw_sig_info_path, sep="\t")
    beta_siginfo = pd.read_csv(beta_siginfo_path, sep="\t")

    ensure_required_columns(sidecar_meta, ["sig_id", "unique_id", "cell_line", "pert_type", "target"], "LINCS sidecar meta")
    ensure_required_columns(raw_sig_info, ["sig_id", "cell_id", "pert_type"], "raw LINCS sig_info")
    ensure_required_columns(beta_siginfo, ["sig_id", "cell_iname", "cmap_name"], "beta siginfo")

    meta = sidecar_meta.merge(
        raw_sig_info[["sig_id", "cell_id", "pert_type", "pert_time", "pert_dose"]],
        on="sig_id",
        how="left",
        suffixes=("", "_raw"),
    ).merge(
        beta_siginfo[["sig_id", "cell_iname", "cmap_name"]],
        on="sig_id",
        how="left",
    )
    meta["cell_line_raw"] = meta["cell_iname"].fillna(meta["cell_id"]).fillna(meta["cell_line"]).astype(str)
    meta["cell_line_std"] = meta["cell_line"].map(standardize_cell_line)
    meta["perturbation_subtype"] = meta["pert_type"].map(normalize_lincs_subtype)
    meta["perturbation_class"] = meta["perturbation_subtype"].map(
        lambda value: perturbation_class_for_subtype(value) if isinstance(value, str) and value else "NA"
    )
    meta["target_raw"] = meta["target"].map(normalize_required_text)
    meta["target_tokens"] = meta["target_raw"].map(
        lambda value: join_tokens(tuple(token for token in value.replace("|", ";").split(";") if token))
    )

    gctx_handle, matrix, row_lookup = build_lincs_matrix_row_lookup(raw_gctx_path)
    try:
        meta["matrix_row_index"] = meta["sig_id"].map(row_lookup)
        meta = meta.loc[meta["matrix_row_index"].notna()].copy().reset_index(drop=True)
        meta["matrix_row_index"] = pd.to_numeric(meta["matrix_row_index"], errors="raise").astype(np.int64)
        alignment = pd.read_csv(alignment_path).sort_values("local_idx", kind="mergesort").reset_index(drop=True)
        gene_index_map = build_lincs_gene_alignment_index(alignment, raw_gene_info_path)
        gene_delta = load_lincs_aligned_rows(
            matrix=matrix,
            row_positions=meta["matrix_row_index"].to_numpy(dtype=np.int64, copy=False),
            gene_index_map=gene_index_map,
        )
    finally:
        gctx_handle.close()

    meta_out = meta[
        [
            "unique_id",
            "cell_line",
            "pert_type",
            "target",
            "dose_val",
            "time_val",
            "sig_id",
            "paired_control_id",
            "cell_line_raw",
            "cell_line_std",
            "perturbation_subtype",
            "perturbation_class",
            "target_raw",
            "target_tokens",
            "matrix_row_index",
        ]
    ].copy()

    lincs_dir = output_root / "lincs"
    lincs_dir.mkdir(parents=True, exist_ok=True)
    gene_pt_path = lincs_dir / "lincs-engine1-gene-delta.pt"
    torch.save({"y_delta_gene": torch.from_numpy(gene_delta)}, gene_pt_path)
    copy_if_needed(alignment_path, lincs_dir / "lincs-gene-alignment.csv")
    copy_if_needed(sidecar_root / "pathway/hallmark-w-2477x50.npy", output_root / "pathway/hallmark-w-2477x50.npy")
    copy_if_needed(pathway_policy_path, output_root / "pathway/lincs-pathway-policy.json")
    audit = build_standardization_audit_rows(
        pd.DataFrame(
            {
                "perturbation_subtype": meta_out["perturbation_subtype"],
                "cell_line_raw": meta_out["cell_line_raw"],
                "cell_line_std": meta_out["cell_line_std"],
                "target_raw": meta_out["target_raw"],
                "target_tokens": meta_out["target_tokens"].map(
                    lambda value: tuple(token for token in normalize_required_text(value).split(";") if token)
                ),
            }
        )
    )
    return meta_out, gene_pt_path, audit


def _to_dense(x: Any) -> np.ndarray:
    if hasattr(x, "toarray"):
        return np.asarray(x.toarray(), dtype=np.float32)
    return np.asarray(x, dtype=np.float32)


def select_controls(
    *,
    candidate_indices: np.ndarray,
    control_ids: pd.Series,
    seed: int,
) -> np.ndarray:
    if candidate_indices.size <= MAX_CONTROLS:
        return np.sort(candidate_indices.astype(np.int64, copy=False), kind="mergesort")
    order = np.argsort(control_ids.iloc[candidate_indices].astype(str).to_numpy(), kind="mergesort")
    sorted_candidates = candidate_indices[order].astype(np.int64, copy=False)
    rng = np.random.default_rng(seed)
    picked = rng.permutation(sorted_candidates.size)[:MAX_CONTROLS]
    return np.sort(sorted_candidates[picked], kind="mergesort")


def load_aligned_scperturb_rows(
    *,
    adata: Any,
    row_indices: Sequence[int],
    gene_index_map: np.ndarray,
) -> np.ndarray:
    row_indices_array = np.asarray(list(row_indices), dtype=np.int64)
    out = np.zeros((int(row_indices_array.size), int(gene_index_map.size)), dtype=np.float32)
    valid_local = np.where(gene_index_map >= 0)[0].astype(np.int64, copy=False)
    if row_indices_array.size == 0 or valid_local.size == 0:
        return out
    source_idx = gene_index_map[valid_local].astype(np.int64, copy=False)
    block = _to_dense(adata[row_indices_array.tolist(), source_idx.tolist()].X)
    out[:, valid_local] = block
    return out


def build_scperturb_subtype_bundle(
    *,
    subtype: str,
    sources: Sequence[ScPerturbSource],
    alignment: pd.DataFrame,
    pathway_w: np.ndarray,
    seed: int,
) -> ScPerturbBuildResult:
    if ad is None:
        raise RuntimeError("anndata is required to build scPerturb bundles from raw h5ad files.")

    template_genes = alignment["gene_symbol"].astype(str).tolist()
    delta_rows: list[np.ndarray] = []
    meta_records: list[dict[str, object]] = []
    audit_records: list[dict[str, object]] = []
    attrition_records: list[dict[str, object]] = []

    next_row_idx = 0
    for source in sources:
        obs = pd.read_csv(source.obs_path)
        cell_column = choose_scperturb_cell_column(obs.columns)
        if cell_column is None:
            attrition_records.append(
                {
                    "dataset_id": source.dataset_id,
                    "perturbation_subtype": subtype,
                    "reason": "missing_cell_line_column",
                    "n_rows": int(len(obs)),
                }
            )
            continue

        ensure_required_columns(obs, ["perturbation_type", "perturbation"], f"{source.obs_path.name}")
        obs = obs.copy()
        obs["dataset_id"] = source.dataset_id
        obs["raw_obs_idx"] = np.arange(len(obs), dtype=np.int64)
        obs["raw_cell_id"] = (
            obs["Unnamed: 0"].astype(str)
            if "Unnamed: 0" in obs.columns
            else obs["raw_obs_idx"].map(lambda idx: f"{source.dataset_id}::{int(idx)}").astype(str)
        )
        obs["cell_line_raw"] = obs[cell_column].map(normalize_required_text)
        obs["cell_line_std"] = obs["cell_line_raw"].map(standardize_cell_line)
        obs["subtype_norm"] = obs["perturbation_type"].map(normalize_scperturb_subtype)
        obs["control_like"] = obs.apply(
            lambda row: is_control_like(
                row.get("perturbation"),
                row.get("perturbation_raw"),
                row.get("target"),
            ),
            axis=1,
        )

        adata = ad.read_h5ad(source.h5ad_path, backed="r")
        try:
            gene_index_map = build_alignment_index(template_genes, list(adata.var_names))

            treated_mask = obs["subtype_norm"].eq(subtype) & (~obs["control_like"])
            treated = obs.loc[treated_mask].copy().sort_values("raw_cell_id", kind="mergesort")
            if treated.empty:
                continue

            control_ids = obs["raw_cell_id"]
            for treated_row in treated.itertuples(index=False):
                target_raw, target_tokens, target_reason = build_scperturb_target_fields(
                    subtype=subtype,
                    target_value=getattr(treated_row, "target", ""),
                    perturbation_value=getattr(treated_row, "perturbation", ""),
                )
                if target_reason is not None:
                    attrition_records.append(
                        {
                            "dataset_id": source.dataset_id,
                            "perturbation_subtype": subtype,
                            "reason": target_reason,
                            "treated_cell_id": treated_row.raw_cell_id,
                        }
                    )
                    continue

                control_mask = obs["control_like"].copy()
                control_mask &= obs["subtype_norm"].eq(subtype)
                control_mask &= obs["cell_line_std"].eq(treated_row.cell_line_std)
                if subtype == "Chemical" and "time" in obs.columns and not pd.isna(getattr(treated_row, "time", np.nan)):
                    control_mask &= pd.to_numeric(obs["time"], errors="coerce").eq(float(getattr(treated_row, "time")))
                candidate_controls = np.flatnonzero(control_mask.to_numpy())
                if candidate_controls.size == 0:
                    attrition_records.append(
                        {
                            "dataset_id": source.dataset_id,
                            "perturbation_subtype": subtype,
                            "reason": "no_controls_available",
                            "treated_cell_id": treated_row.raw_cell_id,
                        }
                    )
                    continue

                selected_controls = select_controls(
                    candidate_indices=candidate_controls,
                    control_ids=control_ids,
                    seed=seed + int(treated_row.raw_obs_idx),
                )
                treated_aligned = load_aligned_scperturb_rows(
                    adata=adata,
                    row_indices=[int(treated_row.raw_obs_idx)],
                    gene_index_map=gene_index_map,
                )[0]
                control_aligned = load_aligned_scperturb_rows(
                    adata=adata,
                    row_indices=selected_controls.tolist(),
                    gene_index_map=gene_index_map,
                )
                delta = treated_aligned - np.mean(control_aligned, axis=0, dtype=np.float32)
                delta_rows.append(np.asarray(delta, dtype=np.float32))

                meta_records.append(
                    {
                        "delta_row_idx": int(next_row_idx),
                        "dataset_std": source.dataset_id,
                        "cell_std": treated_row.cell_line_std,
                        "cell_line_raw": treated_row.cell_line_raw,
                        "target_std": join_tokens(target_tokens),
                        "target_raw": target_raw,
                        "target_tokens": join_tokens(target_tokens),
                        "pert_type": normalize_required_text(treated_row.perturbation_type),
                        "perturbation_subtype": subtype,
                        "perturbation_class": perturbation_class_for_subtype(subtype),
                        "dose_val": float(pd.to_numeric([getattr(treated_row, "dose_value", np.nan)], errors="coerce")[0]),
                        "time_val": float(pd.to_numeric([getattr(treated_row, "time", np.nan)], errors="coerce")[0]),
                        "raw_cell_id": treated_row.raw_cell_id,
                        "raw_obs_idx": int(treated_row.raw_obs_idx),
                        "n_controls_used": int(selected_controls.size),
                    }
                )
                audit_records.append(
                    {
                        "perturbation_subtype": subtype,
                        "cell_line_raw": treated_row.cell_line_raw,
                        "cell_line_std": treated_row.cell_line_std,
                        "target_raw": target_raw,
                        "target_tokens": target_tokens,
                    }
                )
                next_row_idx += 1
        finally:
            adata.file.close()

    meta = dataframe_from_records(meta_records)
    if meta.empty:
        gene_delta = np.empty((0, 2477), dtype=np.float32)
    else:
        gene_delta = np.vstack(delta_rows).astype(np.float32, copy=False)
    pathway_delta = np.asarray(gene_delta @ pathway_w, dtype=np.float32) if gene_delta.size else np.empty((0, 50), dtype=np.float32)
    return ScPerturbBuildResult(
        meta=meta,
        gene_delta=gene_delta,
        pathway_delta=pathway_delta,
        schema_audit=build_standardization_audit_rows(dataframe_from_records(audit_records)),
        attrition=dataframe_from_records(attrition_records),
    )


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()

    config_path = project_root / CONFIG_PATH
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    raw_lincs_root = resolve_config_path(
        project_root,
        config.get("paths", {}).get("osmosis_lincs_raw", DEFAULT_OSMOSIS_LINCS_ROOT),
    )
    raw_scperturb_root = resolve_config_path(
        project_root,
        config.get("paths", {}).get("osmosis_scperturb_raw", DEFAULT_OSMOSIS_SCPERTURB_ROOT),
    )
    sidecar_root = resolve_config_path(
        project_root,
        config.get("paths", {}).get("task1_snapshot", DEFAULT_TASK1_SNAPSHOT_ROOT),
    )
    output_root = resolve_config_path(project_root, args.output_root)
    output_root.mkdir(parents=True, exist_ok=True)

    lincs_meta, gene_pt_path, lincs_audit = build_lincs_meta_and_gene_delta(
        raw_lincs_root=raw_lincs_root,
        sidecar_root=sidecar_root,
        output_root=output_root,
    )
    write_csv(lincs_meta, output_root / "lincs/lincs-engine1-meta.csv")
    write_csv(lincs_audit, output_root / "lincs/lincs-standardization-audit.csv")

    _sidecar_lincs_meta, pathway_w, alignment_path, _policy_path = load_task1_sidecar_assets(sidecar_root)
    alignment = pd.read_csv(alignment_path).sort_values("local_idx", kind="mergesort").reset_index(drop=True)
    sc_sources = discover_scperturb_sources(raw_scperturb_root)
    chemical = build_scperturb_subtype_bundle(
        subtype="Chemical",
        sources=sc_sources,
        alignment=alignment,
        pathway_w=pathway_w,
        seed=args.seed,
    )
    crispr = build_scperturb_subtype_bundle(
        subtype="CRISPR",
        sources=sc_sources,
        alignment=alignment,
        pathway_w=pathway_w,
        seed=args.seed,
    )

    sc_dir = output_root / "scperturb_delta"
    sc_dir.mkdir(parents=True, exist_ok=True)
    write_csv(chemical.meta, sc_dir / "scperturb-chemical-delta-meta.csv")
    np.save(sc_dir / "scperturb-chemical-gene-delta.npy", chemical.gene_delta)
    np.save(sc_dir / "scperturb-chemical-pathway-delta.npy", chemical.pathway_delta)
    write_csv(chemical.schema_audit, sc_dir / "scperturb-chemical-standardization-audit.csv")
    write_csv(chemical.attrition, sc_dir / "scperturb-chemical-attrition.csv")

    write_csv(crispr.meta, sc_dir / "scperturb-crispr-delta-meta.csv")
    np.save(sc_dir / "scperturb-crispr-gene-delta.npy", crispr.gene_delta)
    np.save(sc_dir / "scperturb-crispr-pathway-delta.npy", crispr.pathway_delta)
    write_csv(crispr.schema_audit, sc_dir / "scperturb-crispr-standardization-audit.csv")
    write_csv(crispr.attrition, sc_dir / "scperturb-crispr-attrition.csv")

    fm_root = output_root / "fm_delta"
    fm_root.mkdir(parents=True, exist_ok=True)
    write_json(
        fm_root / "staging_manifest.json",
        {
            "status": "staged_not_materialized",
            "note": "FM rebuild is intentionally downstream of the raw Task1 snapshot candidate build.",
        },
    )
    (output_root / "cross_contract").mkdir(parents=True, exist_ok=True)

    write_json(
        output_root / "snapshot_manifest.json",
        {
            "snapshot_version": "task1_snapshot_v2_candidate",
            "created_at": utc_now_iso(),
            "seed": int(args.seed),
            "roots": {
                "raw_lincs_root": str(raw_lincs_root),
                "raw_scperturb_root": str(raw_scperturb_root),
                "task1_sidecar_root": str(sidecar_root),
                "output_root": str(output_root),
            },
            "datasets": {
                "LINCS": {
                    "n_rows": int(len(lincs_meta)),
                    "subtypes": sorted(lincs_meta["perturbation_subtype"].dropna().astype(str).unique().tolist()),
                    "gene_delta_path": str(gene_pt_path.relative_to(output_root)),
                },
                "scPerturb": {
                    "chemical_rows": int(len(chemical.meta)),
                    "crispr_rows": int(len(crispr.meta)),
                    "source_files": [source.dataset_id for source in sc_sources],
                },
            },
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
