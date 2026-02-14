#!/usr/bin/env python3
"""
Task0 -> Task1 context attrition audit for M2M-Bench.

Purpose
-------
Quantify where and why contexts/rows are dropped from Task0 unified metadata
to Task1 modality analysis inputs.

This script is intentionally written in a "data-analysis script" style:
  1) load
  2) derive stages
  3) diagnose drop reasons
  4) export tables + markdown summary
"""

from __future__ import annotations

import argparse
import json
from pathlib import Path
from typing import Dict, Iterable, List, Tuple

import numpy as np
import pandas as pd
import torch


# ---------------------------------------------------------------------
# Section 1. IO helpers
# ---------------------------------------------------------------------

VALID_SOURCES = {"LINCS", "scPerturb"}
VALID_MODALITIES = {"Chemical", "Genetic"}
CONTEXT_COLS = ["cell_std", "modality", "target_std"]


def _read_parquet_or_csv(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    try:
        return pd.read_parquet(path)
    except Exception as exc:
        fallback = path.with_suffix(".csv")
        if fallback.exists():
            return pd.read_csv(fallback)
        raise RuntimeError(
            f"Cannot read {path}. Missing parquet engine and no CSV fallback at {fallback}."
        ) from exc


def load_unified_meta(bundle_path: Path, unified_meta_path: Path) -> pd.DataFrame:
    """
    Prefer loading from Task0 bundle (torch) to avoid parquet dependency issues.
    Fallback to explicit unified_meta file if bundle path unavailable.
    """
    if bundle_path.exists():
        bundle = torch.load(bundle_path, map_location="cpu", weights_only=False)
        if isinstance(bundle, dict) and "unified_meta" in bundle:
            meta = bundle["unified_meta"]
            if not isinstance(meta, pd.DataFrame):
                raise TypeError("bundle['unified_meta'] is not a DataFrame.")
            return meta.copy()
    return _read_parquet_or_csv(unified_meta_path)


def _write_markdown(path: Path, lines: Iterable[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------
# Section 2. Stage builders
# ---------------------------------------------------------------------

def _is_missing_text(series: pd.Series) -> pd.Series:
    s = series.astype("string")
    return s.isna() | (s.str.strip() == "")


def build_task1_stage(meta: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Reproduce Task1 stage-1 index logic used in src/m2m_bench/task1/groupwise_gap.py:
      - keep only LINCS/scPerturb
      - keep only Chemical/Genetic
      - reset index and assign task1_row_id

    Returns:
      stage1_kept, stage0_with_drop_reason
    """
    need_cols = [
        "cell_std",
        "target_std",
        "modality",
        "source_db",
        "global_idx",
        "chunk_file",
        "chunk_idx",
        "dose_val",
        "time_val",
        "pert_type",
        "cond_id",
        "uid",
    ]
    missing = [c for c in need_cols if c not in meta.columns]
    if missing:
        raise ValueError(f"unified_meta missing columns: {missing}")

    df = meta[need_cols].copy()
    df["drop_reason_stage1"] = "kept"

    missing_core = (
        _is_missing_text(df["cell_std"])
        | _is_missing_text(df["target_std"])
        | _is_missing_text(df["modality"])
        | _is_missing_text(df["source_db"])
    )
    invalid_source = ~df["source_db"].isin(VALID_SOURCES)
    invalid_modality = ~df["modality"].isin(VALID_MODALITIES)

    # Priority: missing core > invalid source > invalid modality
    df.loc[missing_core, "drop_reason_stage1"] = "missing_core_metadata"
    df.loc[(~missing_core) & invalid_source, "drop_reason_stage1"] = "invalid_source_db"
    df.loc[(~missing_core) & (~invalid_source) & invalid_modality, "drop_reason_stage1"] = "invalid_modality"

    kept = df[df["drop_reason_stage1"] == "kept"].copy().reset_index(drop=True)
    kept["dose_val"] = pd.to_numeric(kept["dose_val"], errors="coerce")
    kept["time_val"] = pd.to_numeric(kept["time_val"], errors="coerce")
    kept["task1_row_id"] = np.arange(len(kept), dtype=np.int64)
    return kept, df


def build_overlap_contexts(stage1: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build context-level LINCS/sc overlap table at key=(cell, modality, target).
    Returns:
      context_counts_wide, candidate_rows
    """
    counts = (
        stage1.groupby(CONTEXT_COLS + ["source_db"], dropna=False)
        .size()
        .reset_index(name="n_rows")
    )
    wide = (
        counts.pivot_table(
            index=CONTEXT_COLS,
            columns="source_db",
            values="n_rows",
            fill_value=0,
            aggfunc="sum",
        )
        .reset_index()
    )
    if "LINCS" not in wide.columns:
        wide["LINCS"] = 0
    if "scPerturb" not in wide.columns:
        wide["scPerturb"] = 0

    wide["is_overlap_context"] = (wide["LINCS"] > 0) & (wide["scPerturb"] > 0)
    overlap_keys = wide[wide["is_overlap_context"]][CONTEXT_COLS].copy()
    candidates = stage1.merge(overlap_keys, on=CONTEXT_COLS, how="inner")
    return wide, candidates


# ---------------------------------------------------------------------
# Section 3. Drop-reason diagnostics
# ---------------------------------------------------------------------

def classify_one_sided_context_reasons(stage1: pd.DataFrame, context_wide: pd.DataFrame) -> pd.DataFrame:
    """
    For one-sided contexts, classify *why* counterpart does not exist.

    Reason categories:
      - cell_not_shared_in_other_source
      - target_not_shared_in_other_source
      - cell_and_target_absent_in_other_source
      - cell_target_combination_missing_in_other_source
    """
    one_sided = context_wide[~context_wide["is_overlap_context"]].copy()
    if one_sided.empty:
        return pd.DataFrame(columns=CONTEXT_COLS + ["only_source", "drop_reason_context"])

    one_sided["only_source"] = np.where(
        (one_sided["LINCS"] > 0) & (one_sided["scPerturb"] == 0),
        "LINCS",
        np.where(
            (one_sided["scPerturb"] > 0) & (one_sided["LINCS"] == 0),
            "scPerturb",
            "unknown",
        ),
    )
    one_sided["other_source"] = np.where(
        one_sided["only_source"] == "LINCS",
        "scPerturb",
        np.where(one_sided["only_source"] == "scPerturb", "LINCS", "unknown"),
    )

    # Precompute source+modality -> set(cell), set(target)
    lookup: Dict[Tuple[str, str], Dict[str, set]] = {}
    for (src, mod), grp in stage1.groupby(["source_db", "modality"], dropna=False, sort=False):
        lookup[(str(src), str(mod))] = {
            "cells": set(grp["cell_std"].astype(str).tolist()),
            "targets": set(grp["target_std"].astype(str).tolist()),
        }

    def _check_membership(row: pd.Series) -> Tuple[bool, bool]:
        key = (str(row["other_source"]), str(row["modality"]))
        info = lookup.get(key, {"cells": set(), "targets": set()})
        cell_exists = str(row["cell_std"]) in info["cells"]
        target_exists = str(row["target_std"]) in info["targets"]
        return cell_exists, target_exists

    checks = one_sided.apply(_check_membership, axis=1, result_type="expand")
    checks.columns = ["cell_exists_in_other_source", "target_exists_in_other_source"]
    one_sided = pd.concat([one_sided, checks], axis=1)

    cond_cell_only = (~one_sided["cell_exists_in_other_source"]) & one_sided["target_exists_in_other_source"]
    cond_target_only = one_sided["cell_exists_in_other_source"] & (~one_sided["target_exists_in_other_source"])
    cond_none = (~one_sided["cell_exists_in_other_source"]) & (~one_sided["target_exists_in_other_source"])

    one_sided["drop_reason_context"] = np.select(
        [cond_cell_only, cond_target_only, cond_none],
        [
            "cell_not_shared_in_other_source",
            "target_not_shared_in_other_source",
            "cell_and_target_absent_in_other_source",
        ],
        default="cell_target_combination_missing_in_other_source",
    )

    keep_cols = CONTEXT_COLS + ["only_source", "drop_reason_context"]
    return one_sided[keep_cols].copy()


def build_stage_summary(
    meta_raw: pd.DataFrame,
    stage1: pd.DataFrame,
    candidates: pd.DataFrame,
    matched_pairs: pd.DataFrame | None,
) -> pd.DataFrame:
    s0_context = meta_raw[CONTEXT_COLS].dropna().drop_duplicates()
    s1_context = stage1[CONTEXT_COLS].dropna().drop_duplicates()
    s2_context = candidates[CONTEXT_COLS].dropna().drop_duplicates()

    rows = [
        {"stage": "S0_unified_raw", "n_rows": int(len(meta_raw)), "n_contexts": int(len(s0_context))},
        {"stage": "S1_task1_eligible", "n_rows": int(len(stage1)), "n_contexts": int(len(s1_context))},
        {"stage": "S2_overlap_candidates", "n_rows": int(len(candidates)), "n_contexts": int(len(s2_context))},
    ]

    if matched_pairs is not None and not matched_pairs.empty:
        rows.append(
            {
                "stage": "S3_matched_pairs",
                "n_rows": int(len(matched_pairs)),
                "n_contexts": int(matched_pairs[CONTEXT_COLS].drop_duplicates().shape[0]),
            }
        )
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------
# Section 4. Main analysis
# ---------------------------------------------------------------------

def run_audit(
    bundle_path: Path,
    unified_meta_path: Path,
    matched_pairs_path: Path,
    output_dir: Path,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    analysis_dir = output_dir / "analysis"
    qc_dir = output_dir / "qc"
    analysis_dir.mkdir(parents=True, exist_ok=True)
    qc_dir.mkdir(parents=True, exist_ok=True)

    # 4.1 load base data
    meta_raw = load_unified_meta(bundle_path=bundle_path, unified_meta_path=unified_meta_path)
    stage1, stage0_reason = build_task1_stage(meta_raw)
    context_wide, candidates = build_overlap_contexts(stage1)

    matched_pairs = None
    if matched_pairs_path.exists():
        matched_pairs = _read_parquet_or_csv(matched_pairs_path)

    # 4.2 stage summary
    stage_summary = build_stage_summary(
        meta_raw=meta_raw, stage1=stage1, candidates=candidates, matched_pairs=matched_pairs
    )
    stage_summary.to_csv(analysis_dir / "stage_summary.csv", index=False)

    # 4.3 stage1 drop reasons
    stage1_drop_summary = (
        stage0_reason.groupby("drop_reason_stage1", dropna=False)
        .size()
        .reset_index(name="n_rows")
        .sort_values("n_rows", ascending=False)
    )
    stage1_drop_summary.to_csv(analysis_dir / "stage1_drop_reason_summary.csv", index=False)

    # 4.4 context overlap + one-sided reasons
    context_wide.to_csv(analysis_dir / "context_overlap_counts.csv", index=False)
    one_sided = classify_one_sided_context_reasons(stage1=stage1, context_wide=context_wide)
    one_sided.to_csv(analysis_dir / "context_drop_reasons_detail.csv", index=False)

    one_sided_summary = (
        one_sided.groupby(["modality", "only_source", "drop_reason_context"], dropna=False)
        .size()
        .reset_index(name="n_contexts")
        .sort_values(["modality", "n_contexts"], ascending=[True, False])
    )
    one_sided_summary.to_csv(analysis_dir / "context_drop_reasons_summary.csv", index=False)

    # 4.5 candidate usage in matched pairs
    if matched_pairs is not None and not matched_pairs.empty:
        required_match_cols = ["lincs_row_id", "sc_row_id", "match_type", "cell_std", "modality", "target_std"]
        missing_match_cols = [c for c in required_match_cols if c not in matched_pairs.columns]
        if missing_match_cols:
            raise ValueError(f"matched_pairs missing required columns: {missing_match_cols}")

        used_ids = np.unique(
            np.concatenate(
                [
                    matched_pairs["lincs_row_id"].to_numpy(dtype=np.int64),
                    matched_pairs["sc_row_id"].to_numpy(dtype=np.int64),
                ]
            )
        )

        cand_rows = candidates.copy()
        cand_rows["is_used_in_matched_pairs"] = cand_rows["task1_row_id"].isin(set(used_ids))
        cand_rows["drop_reason_stage3"] = np.where(
            cand_rows["is_used_in_matched_pairs"],
            "used",
            "count_imbalance_unmatched",
        )

        stage3_drop_summary = (
            cand_rows.groupby(["modality", "source_db", "drop_reason_stage3"], dropna=False)
            .size()
            .reset_index(name="n_rows")
            .sort_values(["modality", "source_db", "n_rows"], ascending=[True, True, False])
        )
        stage3_drop_summary.to_csv(analysis_dir / "stage3_drop_reason_summary.csv", index=False)

        # Match-type + protocol gap summary
        match_type_summary = (
            matched_pairs.groupby(["modality", "match_type"], dropna=False)
            .size()
            .reset_index(name="n_pairs")
        )
        match_type_summary.to_csv(analysis_dir / "match_type_summary.csv", index=False)

        for col in ["dose_val_lincs", "dose_val_sc", "time_val_lincs", "time_val_sc"]:
            if col in matched_pairs.columns:
                matched_pairs[col] = pd.to_numeric(matched_pairs[col], errors="coerce")

        if {"dose_val_lincs", "dose_val_sc", "time_val_lincs", "time_val_sc"}.issubset(matched_pairs.columns):
            matched_pairs["dose_logdiff"] = np.abs(
                np.log1p(np.clip(matched_pairs["dose_val_lincs"], a_min=0.0, a_max=None))
                - np.log1p(np.clip(matched_pairs["dose_val_sc"], a_min=0.0, a_max=None))
            )
            matched_pairs["time_absdiff"] = np.abs(matched_pairs["time_val_lincs"] - matched_pairs["time_val_sc"])
            protocol_gap_summary = (
                matched_pairs.groupby(["modality", "match_type"], dropna=False)[["dose_logdiff", "time_absdiff"]]
                .agg(["mean", "median", "std"])
                .reset_index()
            )
            protocol_gap_summary.columns = [
                "_".join([str(x) for x in col if str(x) != ""]).strip("_")
                for col in protocol_gap_summary.columns.to_flat_index()
            ]
            protocol_gap_summary.to_csv(analysis_dir / "protocol_gap_summary.csv", index=False)

    # 4.6 metadata missing diagnostics in stage1
    meta_missing = pd.DataFrame(
        {
            "field": ["dose_val", "time_val", "cond_id", "pert_type"],
            "missing_rate_stage1": [
                float(pd.isna(stage1["dose_val"]).mean()),
                float(pd.isna(stage1["time_val"]).mean()),
                float(_is_missing_text(stage1["cond_id"]).mean()),
                float(_is_missing_text(stage1["pert_type"]).mean()),
            ],
        }
    )
    meta_missing.to_csv(qc_dir / "metadata_missing_summary_stage1.csv", index=False)

    # 4.7 compact markdown report
    report_lines: List[str] = []
    report_lines.append("# Task0 → Task1 Attrition Audit")
    report_lines.append("")
    report_lines.append("## Stage Counts")
    for _, r in stage_summary.iterrows():
        report_lines.append(
            f"- {r['stage']}: rows={int(r['n_rows'])}, contexts={int(r['n_contexts'])}"
        )

    report_lines.append("")
    report_lines.append("## Stage1 Drop Reasons")
    for _, r in stage1_drop_summary.iterrows():
        report_lines.append(f"- {r['drop_reason_stage1']}: {int(r['n_rows'])} rows")

    if not one_sided_summary.empty:
        report_lines.append("")
        report_lines.append("## One-Sided Context Drop Reasons")
        for _, r in one_sided_summary.iterrows():
            report_lines.append(
                f"- modality={r['modality']}, only_source={r['only_source']}, "
                f"reason={r['drop_reason_context']}: {int(r['n_contexts'])} contexts"
            )

    _write_markdown(analysis_dir / "audit_report.md", report_lines)

    # 4.8 run manifest
    run_manifest = {
        "bundle_path": str(bundle_path),
        "unified_meta_path": str(unified_meta_path),
        "matched_pairs_path": str(matched_pairs_path),
        "summary": {
            "n_stage0_rows": int(len(meta_raw)),
            "n_stage1_rows": int(len(stage1)),
            "n_candidate_rows": int(len(candidates)),
            "n_context_overlap": int(context_wide["is_overlap_context"].sum()),
            "n_context_one_sided": int((~context_wide["is_overlap_context"]).sum()),
        },
        "outputs": {
            "analysis_dir": str(analysis_dir),
            "qc_dir": str(qc_dir),
        },
    }
    (output_dir / "run_manifest_task01_attrition.json").write_text(
        json.dumps(run_manifest, indent=2), encoding="utf-8"
    )


def main() -> None:
    parser = argparse.ArgumentParser("Task0->Task1 attrition audit")
    parser.add_argument(
        "--bundle-path",
        type=str,
        default="./outputs/task0_curated/bundle/m2m_task0_bundle.pt",
        help="Task0 bundle path (preferred for robust loading).",
    )
    parser.add_argument(
        "--unified-meta-path",
        type=str,
        default="./outputs/task0_curated/metadata/unified_meta.parquet",
        help="Fallback unified_meta path (parquet/csv).",
    )
    parser.add_argument(
        "--matched-pairs-path",
        type=str,
        default="./outputs/task1/data/m1_matched_pairs.csv",
        help="Task1 matched pairs table (csv/parquet).",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./outputs/task1_audit",
    )
    args = parser.parse_args()

    run_audit(
        bundle_path=Path(args.bundle_path),
        unified_meta_path=Path(args.unified_meta_path),
        matched_pairs_path=Path(args.matched_pairs_path),
        output_dir=Path(args.output_dir),
    )


if __name__ == "__main__":
    main()
