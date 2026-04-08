#!/usr/bin/env python3
"""
Build formal modifier result objects for Figure 3 questions F3.6-F3.8.

Status:
- support-only modifier builder

Manuscript role:
- materializes 7 formal result objects for time, dose, and target-multiplicity
  modifier analyses on Task2 C2G retrieval

Architecture:
- see scripts/ARCHITECTURE.md for script-family classification
- inputs: S5 task2_retrieval_per_query.parquet (C2G direction only)
- outputs: 7 CSVs to manuscript_active/analysis/

Current revision contract:
- descriptive summaries are target-resolved and use (dataset, cell_line, target,
  representation, modifier_level) rows
- dose/time formal stats are per-target Spearman correlations using query-level
  observations within (dataset, cell_line, representation, target) units
- explicit testability flags capture minimum-observation and modifier-variation
  coverage for each target-level correlation
"""

from __future__ import annotations

import argparse
import re
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import spearmanr, wilcoxon


NAS_RUNS_ROOT = Path("/mnt/NAS_21T/ProjectData/M2M/runs")
DEFAULT_TASK2_RETRIEVAL_PER_QUERY = (
    NAS_RUNS_ROOT
    / "s5_multisource_impl_verify_20260311_a"
    / "s5_task2_retrieval_multisource"
    / "task2_retrieval_per_query.parquet"
)
DEFAULT_ANALYSIS_ROOT = NAS_RUNS_ROOT / "manuscript_active" / "analysis"

RETRIEVAL_METRICS = ["mrr_corrected", "hit1_corrected", "hit5_corrected", "hit10_corrected"]
STRAT_KEYS = ["dataset", "cell_line", "representation"]
FORMAL_REPRESENTATIONS = {"Gene", "Pathway"}
MIN_TARGETS_PER_LEVEL = 5
MIN_TARGETS_PER_TEST = 5
MIN_OBS_PER_TARGET = 3
MIN_MODIFIER_UNIQUE = 2
MIN_METRIC_UNIQUE = 2

# Dose binning: log10 breaks (dose values assumed in uM).
DOSE_BINS_LOG10 = [-np.inf, -2.0, -1.0, 0.0, 1.0, np.inf]
DOSE_BIN_LABELS = ["<0.01", "0.01-0.1", "0.1-1.0", "1.0-10.0", ">10.0"]
DOSE_BIN_MIDPOINTS_LOG10: dict[str, float] = {
    "<0.01": -2.5,
    "0.01-0.1": -1.5,
    "0.1-1.0": -0.5,
    "1.0-10.0": 0.5,
    ">10.0": 1.5,
}

SUMMARY_METRIC_COLUMNS = [f"mean_{metric}" for metric in RETRIEVAL_METRICS]
MODIFIER_STATS_COLS = [
    "modifier_type",
    "dataset",
    "cell_line",
    "representation",
    "comparison",
    "metric",
    "test_name",
    "test_status",
    "effect_direction",
    "effect_size",
    "raw_p",
    "bh_q",
    "testable_bool",
    "underpowered_bool",
    "n_levels",
    "n_test_units",
    "modifier_low_level",
    "modifier_high_level",
    "notes",
]
SPEARMAN_STATS_COLS = [
    "dataset",
    "cell_line",
    "representation",
    "target",
    "modifier_type",
    "metric_name",
    "metric",
    "spearman_rho",
    "p_value",
    "bh_q",
    "n_obs",
    "testable_bool",
    "untestable_reason",
    "time_coverage_complete_bool",
    "dose_coverage_complete_bool",
    "time_variation_bool",
    "dose_variation_bool",
]


def parse_args() -> argparse.Namespace:
    p = argparse.ArgumentParser(description="Build F3.6-F3.8 modifier result objects.")
    p.add_argument("--parquet-path", type=Path, default=DEFAULT_TASK2_RETRIEVAL_PER_QUERY)
    p.add_argument("--analysis-root", type=Path, default=DEFAULT_ANALYSIS_ROOT)
    return p.parse_args()


def require_path(path: Path, label: str) -> None:
    if not path.exists():
        raise FileNotFoundError(f"{label} does not exist: {path}")


def bh_adjust(series: pd.Series) -> pd.Series:
    if series.empty:
        return pd.Series(dtype="float64")
    order = np.argsort(series.to_numpy(), kind="mergesort")
    ranked = series.to_numpy()[order]
    n = ranked.size
    adjusted = ranked * n / np.arange(1, n + 1, dtype="float64")
    adjusted = np.minimum.accumulate(adjusted[::-1])[::-1]
    adjusted = np.clip(adjusted, 0.0, 1.0)
    out = np.empty(n, dtype="float64")
    out[order] = adjusted
    return pd.Series(out, index=series.index, dtype="float64")


def _apply_family_bh(stats_df: pd.DataFrame) -> pd.DataFrame:
    out = stats_df.copy()
    out["bh_q"] = np.nan
    testable_mask = out["testable_bool"].fillna(False) & out["raw_p"].notna()
    if testable_mask.sum() == 0:
        return out
    all_families = out[["metric", "test_name"]].apply(tuple, axis=1)
    for family in all_families[testable_mask].unique():
        fam_mask = testable_mask & (all_families == family)
        out.loc[fam_mask, "bh_q"] = bh_adjust(out.loc[fam_mask, "raw_p"])
    return out


def _apply_global_bh(stats_df: pd.DataFrame) -> pd.DataFrame:
    return _apply_family_bh(stats_df)


def _read_parquet_with_target_columns(parquet_path: Path) -> pd.DataFrame:
    base_columns = [
        "dataset",
        "cell_line",
        "representation",
        "direction",
        "query_time",
        "query_dose_value",
        "query_n_targets",
        "mrr_corrected",
        "hit1_corrected",
        "hit5_corrected",
        "hit10_corrected",
    ]
    target_candidates = [
        ["query_target_tokens", "query_target_token"],
        ["query_target_tokens"],
        ["query_target_token"],
    ]
    last_error: Exception | None = None
    for target_columns in target_candidates:
        try:
            return pd.read_parquet(parquet_path, columns=base_columns + target_columns)
        except Exception as exc:  # pragma: no cover - only used against full parquet schemas
            last_error = exc
    if last_error is not None:
        raise last_error
    raise RuntimeError("No target-column read path succeeded.")


def load_c2g_parquet(parquet_path: Path) -> pd.DataFrame:
    require_path(parquet_path, "S5 task2 retrieval per-query parquet")
    df = _read_parquet_with_target_columns(parquet_path)
    df = df.loc[df["direction"].eq("C2G")].copy()
    if "query_target_tokens" not in df.columns:
        df["query_target_tokens"] = df.get("query_target_token", pd.Series(pd.NA, index=df.index))
    elif "query_target_token" in df.columns:
        df["query_target_tokens"] = df["query_target_tokens"].fillna(df["query_target_token"])
    return df


def split_target_tokens(value: object) -> list[str]:
    if pd.isna(value):
        return []
    text = str(value).strip()
    if text in {"", "NA", "nan"}:
        return []
    parts = [
        piece.strip()
        for piece in re.split(r"[;|]", text)
        if piece is not None
    ]
    return sorted({part for part in parts if part and part not in {"NA", "nan"}})


def explode_targets(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["target"] = out["query_target_tokens"].map(split_target_tokens)
    out = out.explode("target", ignore_index=True)
    out = out.loc[out["target"].notna() & out["target"].astype(str).ne("")].copy()
    return out


def _assign_dose_bin(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    log10_dose = pd.Series(np.nan, index=series.index, dtype="float64")
    valid = numeric > 0
    if valid.any():
        log10_dose.loc[valid] = np.log10(numeric.loc[valid])
    return pd.cut(
        log10_dose,
        bins=DOSE_BINS_LOG10,
        labels=DOSE_BIN_LABELS,
        right=True,
    ).astype("object")


def _assign_multiplicity_bin(series: pd.Series) -> pd.Series:
    numeric = pd.to_numeric(series, errors="coerce")
    out = pd.Series(pd.NA, index=series.index, dtype="object")
    out.loc[numeric.eq(1)] = "1"
    out.loc[numeric.eq(2)] = "2"
    out.loc[numeric.ge(3)] = "3+"
    return out


def _summarize_target_level(
    df: pd.DataFrame,
    *,
    level_column: str,
    modifier_type: str,
    modifier_level_unit: str,
) -> pd.DataFrame:
    group_keys = [*STRAT_KEYS, "target", level_column]
    out = (
        df.groupby(group_keys, dropna=False, as_index=False)
        .agg(
            n_queries=(level_column, "size"),
            mean_mrr_corrected=("mrr_corrected", "mean"),
            mean_hit1_corrected=("hit1_corrected", "mean"),
            mean_hit5_corrected=("hit5_corrected", "mean"),
            mean_hit10_corrected=("hit10_corrected", "mean"),
        )
    )
    if level_column != "modifier_level":
        out = out.rename(columns={level_column: "modifier_level"})
    out["modifier_type"] = modifier_type
    out["modifier_level_unit"] = modifier_level_unit
    if modifier_type == "dose":
        out["dose_bin_midpoint_log10"] = out["modifier_level"].map(DOSE_BIN_MIDPOINTS_LOG10)
    return out.sort_values([*STRAT_KEYS, "target", "modifier_level"], kind="mergesort").reset_index(drop=True)


def _filter_formal_representations(df: pd.DataFrame) -> pd.DataFrame:
    return df.loc[df["representation"].isin(FORMAL_REPRESENTATIONS)].copy()


def build_time_effect_summary(df: pd.DataFrame) -> pd.DataFrame:
    exploded = explode_targets(df)
    exploded = _filter_formal_representations(exploded)
    exploded["query_time"] = pd.to_numeric(exploded["query_time"], errors="coerce").round(1)
    exploded = exploded.loc[exploded["query_time"].notna()].copy()
    return _summarize_target_level(
        exploded,
        level_column="query_time",
        modifier_type="time",
        modifier_level_unit="hours",
    )


def build_dose_effect_summary(df: pd.DataFrame) -> pd.DataFrame:
    exploded = explode_targets(df)
    exploded = _filter_formal_representations(exploded)
    exploded["dose_bin"] = _assign_dose_bin(exploded["query_dose_value"])
    exploded = exploded.loc[exploded["dose_bin"].notna()].copy()
    return _summarize_target_level(
        exploded,
        level_column="dose_bin",
        modifier_type="dose",
        modifier_level_unit="log10_uM_bin_assumed",
    )


def build_multiplicity_summary(df: pd.DataFrame) -> pd.DataFrame:
    exploded = explode_targets(df)
    exploded = _filter_formal_representations(exploded)
    exploded["multiplicity_bin"] = _assign_multiplicity_bin(exploded["query_n_targets"])
    exploded = exploded.loc[exploded["multiplicity_bin"].notna()].copy()
    return _summarize_target_level(
        exploded,
        level_column="multiplicity_bin",
        modifier_type="target_multiplicity",
        modifier_level_unit="query_n_targets_bin",
    )


def _not_testable_row(
    *,
    base: dict[str, object],
    comparison: str,
    notes: str,
    low_level: object = pd.NA,
    high_level: object = pd.NA,
    n_test_units: int = 0,
) -> dict[str, object]:
    return {
        **base,
        "comparison": comparison,
        "test_name": "wilcoxon_signed_rank",
        "test_status": "not_tested",
        "effect_direction": pd.NA,
        "effect_size": np.nan,
        "raw_p": np.nan,
        "bh_q": np.nan,
        "underpowered_bool": True,
        "n_test_units": n_test_units,
        "modifier_low_level": low_level,
        "modifier_high_level": high_level,
        "notes": notes,
    }


def _safe_wilcoxon(values_a: pd.Series, values_b: pd.Series) -> tuple[float, float, str]:
    delta = pd.to_numeric(values_a, errors="coerce") - pd.to_numeric(values_b, errors="coerce")
    delta = delta.dropna()
    if delta.empty:
        return (np.nan, np.nan, "not_tested")
    nonzero = delta.loc[delta.ne(0)]
    if nonzero.empty:
        return (0.0, 1.0, "tested")
    _, p_value = wilcoxon(delta, zero_method="wilcox", alternative="two-sided", mode="auto")
    return (float(delta.median()), float(p_value), "tested")


def _build_extreme_level_stats(
    summary: pd.DataFrame,
    *,
    modifier_type: str,
    ordered_levels: list[object] | None = None,
) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if summary.empty:
        return pd.DataFrame(columns=MODIFIER_STATS_COLS)

    metric_map = {metric: f"mean_{metric}" for metric in RETRIEVAL_METRICS}
    for strat, grp in summary.groupby(STRAT_KEYS, sort=True, dropna=False):
        dataset, cell_line, representation = strat
        level_order = ordered_levels if ordered_levels is not None else sorted(grp["modifier_level"].dropna().unique().tolist())
        for metric, value_column in metric_map.items():
            metric_frame = grp.loc[grp[value_column].notna(), ["target", "modifier_level", value_column]].copy()
            level_counts = metric_frame.groupby("modifier_level", dropna=False)["target"].nunique().to_dict()
            valid_levels = [level for level in level_order if level_counts.get(level, 0) >= MIN_TARGETS_PER_LEVEL]
            base = {
                "modifier_type": modifier_type,
                "dataset": dataset,
                "cell_line": cell_line,
                "representation": representation,
                "metric": metric,
                "testable_bool": len(valid_levels) >= 2,
                "underpowered_bool": False,
                "n_levels": int(len(valid_levels)),
            }
            if len(valid_levels) < 2:
                rows.append(
                    _not_testable_row(
                        base=base,
                        comparison=f"{modifier_type}_extreme_levels",
                        notes=f"fewer than 2 modifier levels with >={MIN_TARGETS_PER_LEVEL} targets",
                    )
                )
                continue

            low_level = valid_levels[0]
            high_level = valid_levels[-1]
            paired = (
                metric_frame.loc[metric_frame["modifier_level"].isin([low_level, high_level])]
                .pivot_table(index="target", columns="modifier_level", values=value_column, aggfunc="mean")
            )
            needed_levels = [level for level in [low_level, high_level] if level in paired.columns]
            paired = paired.loc[:, needed_levels]
            if len(needed_levels) < 2:
                rows.append(
                    _not_testable_row(
                        base=base,
                        comparison=f"{modifier_type}_{low_level}_vs_{high_level}",
                        notes="no matched targets across the extreme modifier levels",
                        low_level=low_level,
                        high_level=high_level,
                    )
                )
                continue

            paired = paired.dropna()
            n_test_units = int(len(paired))
            if n_test_units < MIN_TARGETS_PER_TEST:
                rows.append(
                    _not_testable_row(
                        base=base,
                        comparison=f"{modifier_type}_{low_level}_vs_{high_level}",
                        notes=f"matched target units below threshold ({n_test_units} < {MIN_TARGETS_PER_TEST})",
                        low_level=low_level,
                        high_level=high_level,
                        n_test_units=n_test_units,
                    )
                )
                continue

            median_delta, raw_p, test_status = _safe_wilcoxon(paired[low_level], paired[high_level])
            effect_direction = (
                "level_a_gt_b"
                if median_delta > 0
                else ("level_b_gt_a" if median_delta < 0 else "tie_or_zero_delta")
            )
            rows.append(
                {
                    **base,
                    "comparison": f"{modifier_type}_{low_level}_vs_{high_level}",
                    "test_name": "wilcoxon_signed_rank",
                    "test_status": test_status,
                    "effect_direction": effect_direction,
                    "effect_size": median_delta,
                    "raw_p": raw_p,
                    "bh_q": np.nan,
                    "n_test_units": n_test_units,
                    "modifier_low_level": low_level,
                    "modifier_high_level": high_level,
                    "notes": (
                        "Matched target-resolved comparison across extreme modifier levels; "
                        f"n_test_units counts shared (cell_line, target) units ({n_test_units})."
                    ),
                }
            )
    out = pd.DataFrame(rows, columns=MODIFIER_STATS_COLS)
    return _apply_family_bh(out).sort_values([*STRAT_KEYS, "metric", "comparison"], kind="mergesort").reset_index(drop=True)


def build_time_effect_stats(df: pd.DataFrame) -> pd.DataFrame:
    return build_modifier_spearman_stats(df, modifier_type="time")


def build_dose_effect_stats(df: pd.DataFrame) -> pd.DataFrame:
    return build_modifier_spearman_stats(df, modifier_type="dose")


def build_multiplicity_stats(df: pd.DataFrame) -> pd.DataFrame:
    summary = build_multiplicity_summary(df)
    rows: list[dict[str, object]] = []
    metric_map = {metric: f"mean_{metric}" for metric in RETRIEVAL_METRICS}
    for strat, grp in summary.groupby(STRAT_KEYS, sort=True, dropna=False):
        dataset, cell_line, representation = strat
        for metric, value_column in metric_map.items():
            metric_frame = grp.loc[grp[value_column].notna(), ["target", "modifier_level", value_column]].copy()
            level_counts = metric_frame.groupby("modifier_level", dropna=False)["target"].nunique().to_dict()
            single_count = int(level_counts.get("1", 0))
            multi_count = int(
                metric_frame.loc[metric_frame["modifier_level"].isin(["2", "3+"]), "target"].nunique()
            )
            base = {
                "modifier_type": "target_multiplicity",
                "dataset": dataset,
                "cell_line": cell_line,
                "representation": representation,
                "metric": metric,
                "testable_bool": single_count >= MIN_TARGETS_PER_LEVEL and multi_count >= MIN_TARGETS_PER_LEVEL,
                "underpowered_bool": False,
                "n_levels": int(sum(level_counts.get(level, 0) >= MIN_TARGETS_PER_LEVEL for level in ["1", "2", "3+"])),
            }
            if not base["testable_bool"]:
                rows.append(
                    _not_testable_row(
                        base=base,
                        comparison="single_vs_multi",
                        notes=(
                            "single or multi strata did not reach the matched-target threshold "
                            f"({single_count} single, {multi_count} multi)"
                        ),
                        low_level="1",
                        high_level="multi",
                    )
                )
                continue

            pivot = metric_frame.pivot_table(index="target", columns="modifier_level", values=value_column, aggfunc="mean")
            multi_columns = [column for column in ["2", "3+"] if column in pivot.columns]
            if "1" not in pivot.columns or not multi_columns:
                rows.append(
                    _not_testable_row(
                        base=base,
                        comparison="single_vs_multi",
                        notes="single or multi columns missing after target pivot",
                        low_level="1",
                        high_level="multi",
                    )
                )
                continue

            pivot["multi"] = pivot[multi_columns].mean(axis=1, skipna=True)
            paired = pivot.loc[:, ["1", "multi"]].dropna()
            n_test_units = int(len(paired))
            if n_test_units < MIN_TARGETS_PER_TEST:
                rows.append(
                    _not_testable_row(
                        base=base,
                        comparison="single_vs_multi",
                        notes=f"matched target units below threshold ({n_test_units} < {MIN_TARGETS_PER_TEST})",
                        low_level="1",
                        high_level="multi",
                        n_test_units=n_test_units,
                    )
                )
                continue

            median_delta, raw_p, test_status = _safe_wilcoxon(paired["1"], paired["multi"])
            effect_direction = (
                "single_gt_multi"
                if median_delta > 0
                else ("multi_gt_single" if median_delta < 0 else "tie_or_zero_delta")
            )
            rows.append(
                {
                    **base,
                    "comparison": "single_vs_multi",
                    "test_name": "wilcoxon_signed_rank",
                    "test_status": test_status,
                    "effect_direction": effect_direction,
                    "effect_size": median_delta,
                    "raw_p": raw_p,
                    "bh_q": np.nan,
                    "n_test_units": n_test_units,
                    "modifier_low_level": "1",
                    "modifier_high_level": "multi",
                    "notes": (
                        "Matched target-resolved single-vs-multi comparison; "
                        f"n_test_units counts shared (cell_line, target) units ({n_test_units})."
                    ),
                }
            )
    out = pd.DataFrame(rows, columns=MODIFIER_STATS_COLS)
    return _apply_family_bh(out).sort_values([*STRAT_KEYS, "metric", "comparison"], kind="mergesort").reset_index(drop=True)


def _prepare_modifier_rows(df: pd.DataFrame, *, modifier_type: str) -> pd.DataFrame:
    exploded = explode_targets(df)
    exploded = _filter_formal_representations(exploded)
    if modifier_type == "time":
        modifier_values = pd.to_numeric(exploded["query_time"], errors="coerce").round(1)
    elif modifier_type == "dose":
        modifier_values = pd.to_numeric(exploded["query_dose_value"], errors="coerce")
    else:
        raise ValueError(f"Unsupported modifier_type: {modifier_type}")
    exploded = exploded.assign(modifier_value=modifier_values)
    coverage_complete = exploded.groupby([*STRAT_KEYS, "target"], dropna=False)["modifier_value"].transform(
        lambda s: s.notna().all()
    )
    exploded = exploded.assign(coverage_complete=coverage_complete)
    return exploded.loc[exploded["modifier_value"].notna()].copy()


def _testability_flags(
    modifier_series: pd.Series,
    metric_series: pd.Series,
    n_obs: int,
) -> tuple[bool, bool, bool, bool, str | pd.NA]:
    meets_min_obs = n_obs >= MIN_OBS_PER_TARGET
    modifier_unique = modifier_series.nunique(dropna=True)
    modifier_varies = modifier_unique >= MIN_MODIFIER_UNIQUE
    metric_unique = metric_series.nunique(dropna=True)
    metric_varies = metric_unique >= MIN_METRIC_UNIQUE
    testable = meets_min_obs and modifier_varies and metric_varies
    reasons: list[str] = []
    if not meets_min_obs:
        reasons.append(f"insufficient_observations_lt_{MIN_OBS_PER_TARGET}")
    if not modifier_varies:
        reasons.append(f"no_modifier_variation_lt_{MIN_MODIFIER_UNIQUE}_unique")
    if not metric_varies:
        reasons.append(f"no_metric_variation_lt_{MIN_METRIC_UNIQUE}_unique")
    return meets_min_obs, modifier_varies, metric_varies, testable, "; ".join(reasons) if reasons else pd.NA


def _modifier_flag_fields(
    *,
    modifier_type: str,
    coverage_complete: bool,
    variation_bool: bool,
) -> dict[str, object]:
    if modifier_type == "time":
        return {
            "time_coverage_complete_bool": coverage_complete,
            "dose_coverage_complete_bool": pd.NA,
            "time_variation_bool": variation_bool,
            "dose_variation_bool": pd.NA,
        }
    if modifier_type == "dose":
        return {
            "time_coverage_complete_bool": pd.NA,
            "dose_coverage_complete_bool": coverage_complete,
            "time_variation_bool": pd.NA,
            "dose_variation_bool": variation_bool,
        }
    raise ValueError(f"Unsupported modifier_type for flags: {modifier_type}")


def _build_spearman_stats(prepared: pd.DataFrame, *, modifier_type: str) -> pd.DataFrame:
    rows: list[dict[str, object]] = []
    if prepared.empty:
        return pd.DataFrame(columns=SPEARMAN_STATS_COLS)

    for strat, grp in prepared.groupby([*STRAT_KEYS, "target"], sort=True, dropna=False):
        dataset, cell_line, representation, target = strat
        for metric in RETRIEVAL_METRICS:
            metric_values = pd.to_numeric(grp[metric], errors="coerce")
            modifier_values = pd.to_numeric(grp["modifier_value"], errors="coerce")
            valid = modifier_values.notna() & metric_values.notna()
            n_obs = int(valid.sum())
            modifier_valid = modifier_values.loc[valid]
            metric_valid = metric_values.loc[valid]
            coverage_complete = bool(grp["coverage_complete"].all())
            meets_min_obs, modifier_varies, metric_varies, testable, reason = _testability_flags(
                modifier_valid,
                metric_valid,
                n_obs,
            )
            if testable:
                rho, p_value = spearmanr(modifier_valid, metric_values.loc[valid], nan_policy="omit")
                if np.isnan(rho) or np.isnan(p_value):
                    testable = False
                    reason = (
                        "spearman_undefined"
                        if pd.isna(reason)
                        else f"{reason}; spearman_undefined"
                    )
                    rho, p_value = (np.nan, np.nan)
            else:
                rho, p_value = (np.nan, np.nan)
            rows.append(
                {
                    "dataset": dataset,
                    "cell_line": cell_line,
                    "representation": representation,
                    "target": target,
                    "modifier_type": modifier_type,
                    "metric_name": metric,
                    "metric": metric,
                    "spearman_rho": rho,
                    "p_value": p_value,
                    "bh_q": np.nan,
                    "n_obs": n_obs,
                    "testable_bool": testable,
                    "untestable_reason": reason,
                    **_modifier_flag_fields(
                        modifier_type=modifier_type,
                        coverage_complete=coverage_complete,
                        variation_bool=modifier_varies,
                    ),
                }
            )
    return pd.DataFrame(rows, columns=SPEARMAN_STATS_COLS)


def _apply_spearman_bh(stats_df: pd.DataFrame) -> pd.DataFrame:
    out = stats_df.copy()
    out["bh_q"] = np.nan
    testable_mask = out["testable_bool"].fillna(False) & out["p_value"].notna()
    if testable_mask.sum() == 0:
        return out
    families = out[["representation", "metric_name", "modifier_type"]].apply(tuple, axis=1)
    for family in families[testable_mask].unique():
        fam_mask = testable_mask & (families == family)
        out.loc[fam_mask, "bh_q"] = bh_adjust(out.loc[fam_mask, "p_value"])
    return out


def build_modifier_spearman_stats(df: pd.DataFrame, *, modifier_type: str) -> pd.DataFrame:
    prepared = _prepare_modifier_rows(df, modifier_type=modifier_type)
    stats = _build_spearman_stats(prepared, modifier_type=modifier_type)
    stats = _apply_spearman_bh(stats)
    return stats.sort_values([*STRAT_KEYS, "target", "metric_name"], kind="mergesort").reset_index(drop=True)


def build_dose_time_spearman_stats(time_stats: pd.DataFrame, dose_stats: pd.DataFrame) -> pd.DataFrame:
    combined = pd.concat([time_stats, dose_stats], ignore_index=True)
    return combined.sort_values(
        ["modifier_type", *STRAT_KEYS, "target", "metric_name"],
        kind="mergesort",
    ).reset_index(drop=True)


def build_modifier_comparison_statistics(
    time_stats: pd.DataFrame,
    dose_stats: pd.DataFrame,
    multiplicity_stats: pd.DataFrame,
) -> pd.DataFrame:
    combined = pd.concat([time_stats, dose_stats, multiplicity_stats], ignore_index=True)
    combined["bh_q"] = np.nan
    combined = _apply_global_bh(combined)
    return combined.sort_values(["modifier_type", *STRAT_KEYS, "comparison", "metric"], kind="mergesort").reset_index(drop=True)


def main() -> int:
    args = parse_args()
    require_path(args.parquet_path, "S5 task2 retrieval parquet")
    args.analysis_root.mkdir(parents=True, exist_ok=True)

    print(f"Loading C2G rows from {args.parquet_path} ...")
    df = load_c2g_parquet(args.parquet_path)
    print(f"  C2G rows: {len(df):,}")

    print("Building time effect summary ...")
    time_summary = build_time_effect_summary(df)
    out_ts = args.analysis_root / "figure3_task2_c2g_time_effect_summary.csv"
    time_summary.to_csv(out_ts, index=False)
    print(f"  -> {out_ts} ({len(time_summary):,} rows)")

    print("Building time effect stats ...")
    time_stats = build_time_effect_stats(df)
    out_tst = args.analysis_root / "figure3_task2_c2g_time_effect_stats.csv"
    time_stats.to_csv(out_tst, index=False)
    print(f"  -> {out_tst} ({len(time_stats):,} rows)")

    print("Building dose effect summary ...")
    dose_summary = build_dose_effect_summary(df)
    out_ds = args.analysis_root / "figure3_task2_c2g_dose_effect_summary.csv"
    dose_summary.to_csv(out_ds, index=False)
    print(f"  -> {out_ds} ({len(dose_summary):,} rows)")

    print("Building dose effect stats ...")
    dose_stats = build_dose_effect_stats(df)
    out_dst = args.analysis_root / "figure3_task2_c2g_dose_effect_stats.csv"
    dose_stats.to_csv(out_dst, index=False)
    print(f"  -> {out_dst} ({len(dose_stats):,} rows)")

    print("Building unified dose/time Spearman stats ...")
    dose_time_stats = build_dose_time_spearman_stats(time_stats, dose_stats)
    out_dt = args.analysis_root / "figure3_task2_c2g_dose_time_spearman_stats.csv"
    dose_time_stats.to_csv(out_dt, index=False)
    print(f"  -> {out_dt} ({len(dose_time_stats):,} rows)")

    print("Building target multiplicity summary ...")
    mult_summary = build_multiplicity_summary(df)
    out_ms = args.analysis_root / "figure3_task2_c2g_target_multiplicity_summary.csv"
    mult_summary.to_csv(out_ms, index=False)
    print(f"  -> {out_ms} ({len(mult_summary):,} rows)")

    print("Building target multiplicity stats ...")
    mult_stats = build_multiplicity_stats(df)
    out_mst = args.analysis_root / "figure3_task2_c2g_target_multiplicity_stats.csv"
    mult_stats.to_csv(out_mst, index=False)
    print(f"  -> {out_mst} ({len(mult_stats):,} rows)")

    print("Done.")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
