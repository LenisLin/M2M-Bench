#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Chem2Gen-Bench — Task 3 Visualization Prep (Script 11)

Goal
----
Prepare R-ready CSVs for 3 Task3 figures, WITHOUT recomputing Task3 metrics.

Figure 1 (Scoreboard)
  - Unified long table for heatmap/table plotting, with panels:
      All + {Cleanest, Family, Promiscuous}
  - Sources: Script10 outputs (PerQuery/PerTarget) + tier annotation from Drug_meta.

Figure 2 (Paired deltas & lollipop-ready summary)
  (a) viz_stats_deltas_query.csv:
      - Query-level paired deltas:
          (A) Systema - Standard, within each Track & Direction & Query_ID
          (B) Track - Gene, within each View & Direction & Query_ID
      - Output: tidy long table (one row per paired delta observation per metric)

  (b) viz_fig2_lollipop.csv (NEW/UPDATED):
      - Lollipop-ready summary table with stderr:
          ComparePanel in { View(Systema-Standard), Standard(Track-Gene) }
      - Metrics (ONLY):
          Pairwise: CentroidCosine
          Pairwise: ED_score = -log10(max(-edist_mean, eps))   (higher is better)
          Retrieval: Success_CRISPR2Drug (Success_Score), Direction fixed to --fig2_direction (default CRISPR->Drug)
      - Pairing:
          Pairwise paired by Target
          Retrieval paired by Query_ID
      - Summary stats:
          MeanDiff, MedianDiff, StdErrMean, StdErrMedian (bootstrap), Wilcoxon p, -log10(p)

Figure 3 (17 targets; tier-stratified gains vs Gene baseline)
  - Target-level aggregation from Retrieval_PerQuery (mean over queries per target)
  - Pairwise_PerTarget already target-level
  - For each (View, Direction, Target), compute Track-vs-Gene deltas
  - Join Target_Tier from Drug_meta majority mapping.

Required inputs (from Script10; under --analysis_dir):
  - Task3_Pairwise_PerTarget.csv
  - Task3_Pairwise_Summary.csv
  - Task3_Retrieval_PerQuery.csv
  - Task3_Retrieval_Summary.csv

Optional inputs (for tier annotation; under --eval_dir):
  - Drug_meta.csv
  - CRISPR_meta.csv (not required)

Outputs (CSV; to --out_dir):
  - viz_scoreboard_long.csv
  - viz_scoreboard_wide.csv
  - viz_stats_deltas_query.csv
  - viz_target_gain.csv
  - annotated_pairwise_per_target.csv
  - annotated_retrieval_per_query.csv
  - viz_fig2_lollipop.csv

Notes
-----
- NaN-safe: missing pairs are dropped.
- Rank metrics:
    Rank_Improvement = BaselineRank - ComparedRank  (positive=better)
"""

from __future__ import annotations

import os
import argparse
import warnings
from typing import Dict, Tuple, List

import numpy as np
import pandas as pd


# ---------------------------
# Constants / helpers
# ---------------------------

TIER_ORDER = [
    "The Cleanest Hits",
    "The Family Hits",
    "The Promiscuous Hits",
    "Control",
    "Unknown",
]

TIER_CANONICAL = {
    "cleanest": "The Cleanest Hits",
    "family": "The Family Hits",
    "promiscuous": "The Promiscuous Hits",
    "control": "Control",
}

DEFAULT_PANELS = ["All", "The Cleanest Hits", "The Family Hits", "The Promiscuous Hits"]


def _ensure_dir(path: str) -> None:
    os.makedirs(path, exist_ok=True)


def _read_required_csv(path: str, name: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(f"[CRITICAL] Missing required input: {name} -> {path}")
    df = pd.read_csv(path)
    if df.empty:
        raise ValueError(f"[CRITICAL] Input is empty: {name} -> {path}")
    return df


def _read_optional_csv(path: str) -> pd.DataFrame:
    if not path or (not os.path.exists(path)):
        return pd.DataFrame()
    # Drug_meta / CRISPR_meta were written with index_col=0 in your pipeline
    try:
        return pd.read_csv(path, index_col=0)
    except Exception:
        return pd.read_csv(path)


def _must_have_cols(df: pd.DataFrame, cols: List[str], label: str) -> None:
    miss = [c for c in cols if c not in df.columns]
    if miss:
        raise ValueError(f"[{label}] Missing required columns: {miss}")


def _canonicalize_tier(x: object) -> str:
    if x is None or (isinstance(x, float) and np.isnan(x)):
        return "Unknown"
    s = str(x).strip()
    if s in TIER_ORDER:
        return s
    s_low = s.lower()
    for k, v in TIER_CANONICAL.items():
        if k in s_low:
            return v
    return s if s else "Unknown"


def _scale01_clip(
    series: pd.Series,
    q: Tuple[float, float] = (0.05, 0.95),
    higher_is_better: bool = True,
) -> pd.Series:
    """
    Quantile-clipped min-max scaling to [0,1] within a block.
    If higher_is_better=False, invert after scaling so larger always means better.
    """
    x = pd.to_numeric(series, errors="coerce")
    if x.notna().sum() == 0:
        return pd.Series([np.nan] * len(series), index=series.index)

    lo = x.quantile(q[0])
    hi = x.quantile(q[1])
    if not np.isfinite(lo) or not np.isfinite(hi) or lo == hi:
        scaled = pd.Series([0.5] * len(series), index=series.index, dtype=float)
    else:
        x2 = x.clip(lower=lo, upper=hi)
        scaled = (x2 - lo) / (hi - lo)

    if not higher_is_better:
        scaled = 1.0 - scaled
    return scaled.astype(float)


def _safe_wilcoxon_p(diff: np.ndarray) -> float:
    """
    Two-sided Wilcoxon signed-rank test p-value vs 0.
    Returns NaN if scipy unavailable or insufficient samples.
    """
    x = np.asarray(diff, dtype=float)
    x = x[np.isfinite(x)]
    if x.size < 3:
        return float("nan")
    try:
        from scipy.stats import wilcoxon
        res = wilcoxon(x, zero_method="wilcox", correction=False, alternative="two-sided", mode="auto")
        return float(res.pvalue)
    except Exception:
        return float("nan")


def _neglog10_p(p: float) -> float:
    if p is None or (isinstance(p, float) and (not np.isfinite(p))):
        return float("nan")
    p2 = float(p)
    if p2 <= 0:
        p2 = 1e-300
    return float(-np.log10(max(p2, 1e-300)))


def _stderr_mean(diff: np.ndarray) -> float:
    """Standard error of the mean."""
    x = np.asarray(diff, dtype=float)
    x = x[np.isfinite(x)]
    n = int(x.size)
    if n < 2:
        return float("nan")
    sd = float(np.std(x, ddof=1))
    return sd / float(np.sqrt(n))


def _stderr_median_bootstrap(
    diff: np.ndarray,
    rng: np.random.Generator,
    n_boot: int = 200,
    max_n: int = 5000,
) -> float:
    """
    Bootstrap SE of median.
    If n > max_n: subsample to max_n without replacement, bootstrap on that sample,
                 then scale SE by sqrt(n_sub / n) as a rough approximation.
    """
    x = np.asarray(diff, dtype=float)
    x = x[np.isfinite(x)]
    n = int(x.size)
    if n < 3:
        return float("nan")

    if n > int(max_n):
        n_sub = int(max_n)
        x0 = rng.choice(x, size=n_sub, replace=False)
        scale = float(np.sqrt(n_sub / n))
    else:
        x0 = x
        scale = 1.0

    m = int(x0.size)
    if m < 3:
        return float("nan")

    meds = np.empty(int(n_boot), dtype=float)
    for b in range(int(n_boot)):
        samp = rng.choice(x0, size=m, replace=True)
        meds[b] = float(np.median(samp))

    se_sub = float(np.std(meds, ddof=1)) if meds.size > 1 else float("nan")
    return se_sub * scale


def _edist_to_edscore(edist_mean: pd.Series, eps: float = 1e-12) -> pd.Series:
    """
    Convert "negated edist" (edist_mean typically negative; closer to 0 is better) to a log score:
        ED_score = -log10(max(-edist_mean, eps))
    Higher ED_score means better (smaller distance).
    If (-edist_mean) <= eps => NaN.
    """
    x = pd.to_numeric(edist_mean, errors="coerce")
    dist = -x
    dist = dist.where(dist > eps, np.nan)
    return (-np.log10(dist)).astype(float)


# ---------------------------
# Tier mapping
# ---------------------------

def build_target_tier_mapping(drug_meta: pd.DataFrame) -> pd.DataFrame:
    """
    Build a Target -> majority specificity tier mapping from Drug_meta.csv.

    Expected columns in drug_meta:
      - specificity_tier
      - clean_target_mapped (may contain 'A;B;C')
    """
    if drug_meta is None or drug_meta.empty:
        return pd.DataFrame(columns=["Target", "Target_Tier", "Tier_Counts"])

    needed = ["specificity_tier", "clean_target_mapped"]
    for c in needed:
        if c not in drug_meta.columns:
            warnings.warn(f"[TierMap] Drug meta missing column: {c}. Tier mapping will be unavailable.")
            return pd.DataFrame(columns=["Target", "Target_Tier", "Tier_Counts"])

    dm = drug_meta.copy()
    dm["specificity_tier"] = dm["specificity_tier"].map(_canonicalize_tier)

    targets = (
        dm["clean_target_mapped"]
        .astype(str)
        .fillna("")
        .str.split(";")
        .explode()
        .str.strip()
    )
    dm2 = dm.loc[targets.index].copy()
    dm2["Target"] = targets.values
    dm2 = dm2[(dm2["Target"].notna()) & (dm2["Target"].astype(str).str.len() > 0)]

    if dm2.empty:
        return pd.DataFrame(columns=["Target", "Target_Tier", "Tier_Counts"])

    ct = (
        dm2.groupby(["Target", "specificity_tier"], dropna=False)
        .size()
        .reset_index(name="n")
    )

    def _tier_rank(t: str) -> int:
        return TIER_ORDER.index(t) if t in TIER_ORDER else (len(TIER_ORDER) - 1)

    rows = []
    for tgt, sub in ct.groupby("Target"):
        sub2 = sub.copy()
        sub2["tier_rank"] = sub2["specificity_tier"].apply(_tier_rank)
        # majority; tie-break by higher specificity (lower rank)
        sub2 = sub2.sort_values(["n", "tier_rank"], ascending=[False, True])
        tier = str(sub2.iloc[0]["specificity_tier"])
        counts = ";".join([f"{r.specificity_tier}:{int(r.n)}" for r in sub2.itertuples(index=False)])
        rows.append({"Target": str(tgt), "Target_Tier": tier, "Tier_Counts": counts})

    return pd.DataFrame(rows)


def attach_tiers_to_retrieval_per_query(
    perq: pd.DataFrame,
    drug_meta: pd.DataFrame,
    target_tier_map: pd.DataFrame,
) -> pd.DataFrame:
    """
    Attach specificity_tier to each retrieval query.
      - Drug queries: tier from Drug_meta (Query_ID should match Drug_meta index)
      - CRISPR queries: tier inferred from True_Target via target majority tier (from drugs)
    """
    perq = perq.copy()
    _must_have_cols(perq, ["Query_Modality", "Query_ID", "True_Target"], "Task3_Retrieval_PerQuery")

    tgt2tier = {}
    if target_tier_map is not None and (not target_tier_map.empty):
        tgt2tier = dict(zip(target_tier_map["Target"].astype(str), target_tier_map["Target_Tier"].astype(str)))

    q2tier = {}
    if drug_meta is not None and (not drug_meta.empty) and ("specificity_tier" in drug_meta.columns):
        dm = drug_meta.copy()
        dm["specificity_tier"] = dm["specificity_tier"].map(_canonicalize_tier)
        dm.index = dm.index.astype(str)
        q2tier = dm["specificity_tier"].to_dict()

    def _infer(row) -> str:
        qm = str(row["Query_Modality"]).strip().lower()
        qid = str(row["Query_ID"])
        tgt = str(row["True_Target"])
        if qm == "drug":
            return _canonicalize_tier(q2tier.get(qid, "Unknown"))
        return _canonicalize_tier(tgt2tier.get(tgt, "Unknown"))

    perq["specificity_tier"] = perq.apply(_infer, axis=1).map(_canonicalize_tier)
    return perq


def attach_tiers_to_pairwise_per_target(per_tgt: pd.DataFrame, target_tier_map: pd.DataFrame) -> pd.DataFrame:
    per_tgt = per_tgt.copy()
    _must_have_cols(per_tgt, ["Target"], "Task3_Pairwise_PerTarget")

    tgt2tier = {}
    tgt2counts = {}
    if target_tier_map is not None and (not target_tier_map.empty):
        tgt2tier = dict(zip(target_tier_map["Target"].astype(str), target_tier_map["Target_Tier"].astype(str)))
        tgt2counts = dict(zip(target_tier_map["Target"].astype(str), target_tier_map["Tier_Counts"].astype(str)))

    per_tgt["specificity_tier"] = per_tgt["Target"].astype(str).map(
        lambda t: _canonicalize_tier(tgt2tier.get(t, "Unknown"))
    )
    per_tgt["tier_counts"] = per_tgt["Target"].astype(str).map(lambda t: tgt2counts.get(t, ""))
    return per_tgt


# ---------------------------
# Scoreboard (Figure 1)
# ---------------------------

def summarize_pairwise_by_panel(per_tgt: pd.DataFrame, panels: List[str]) -> pd.DataFrame:
    req = ["Track", "View", "Target", "centroid_cosine", "edist_mean", "specificity_tier"]
    _must_have_cols(per_tgt, req, "Pairwise_PerTarget(annotated)")

    df = per_tgt.copy()
    df["Panel"] = df["specificity_tier"].map(_canonicalize_tier)

    out_rows = []

    all_sum = (
        df.groupby(["Track", "View"], dropna=False)
        .agg(
            n_targets=("Target", "nunique"),
            mean_centroid_cos=("centroid_cosine", "mean"),
            median_centroid_cos=("centroid_cosine", "median"),
            mean_edist=("edist_mean", "mean"),
            median_edist=("edist_mean", "median"),
        )
        .reset_index()
    )
    all_sum["Panel"] = "All"
    out_rows.append(all_sum)

    tiers = [p for p in panels if p != "All"]
    for tier in tiers:
        sub = df[df["Panel"] == tier]
        if sub.empty:
            continue
        s = (
            sub.groupby(["Track", "View"], dropna=False)
            .agg(
                n_targets=("Target", "nunique"),
                mean_centroid_cos=("centroid_cosine", "mean"),
                median_centroid_cos=("centroid_cosine", "median"),
                mean_edist=("edist_mean", "mean"),
                median_edist=("edist_mean", "median"),
            )
            .reset_index()
        )
        s["Panel"] = tier
        out_rows.append(s)

    out = pd.concat(out_rows, axis=0, ignore_index=True)
    return out


def summarize_retrieval_by_panel(per_q: pd.DataFrame, panels: List[str]) -> pd.DataFrame:
    req = ["Track", "View", "Direction", "MRR", "Success_Score", "True_Rank", "specificity_tier", "Query_ID"]
    _must_have_cols(per_q, req, "Retrieval_PerQuery(annotated)")

    df = per_q.copy()
    df["Panel"] = df["specificity_tier"].map(_canonicalize_tier)

    out_rows = []

    all_sum = (
        df.groupby(["Track", "View", "Direction"], dropna=False)
        .agg(
            n_queries=("Query_ID", "count"),
            mean_mrr=("MRR", "mean"),
            mean_success=("Success_Score", "mean"),
            median_rank=("True_Rank", "median"),
        )
        .reset_index()
    )
    all_sum["Panel"] = "All"
    out_rows.append(all_sum)

    tiers = [p for p in panels if p != "All"]
    for tier in tiers:
        sub = df[df["Panel"] == tier]
        if sub.empty:
            continue
        s = (
            sub.groupby(["Track", "View", "Direction"], dropna=False)
            .agg(
                n_queries=("Query_ID", "count"),
                mean_mrr=("MRR", "mean"),
                mean_success=("Success_Score", "mean"),
                median_rank=("True_Rank", "median"),
            )
            .reset_index()
        )
        s["Panel"] = tier
        out_rows.append(s)

    out = pd.concat(out_rows, axis=0, ignore_index=True)
    return out


def build_scoreboard_long(pair_sum_panel: pd.DataFrame, ret_sum_panel: pd.DataFrame) -> pd.DataFrame:
    """
    Tidy long scoreboard table for R.
    Columns: MetricGroup, Metric, Track, View, Panel, Direction, Value, Scaled_Best, Rank_Best, N
    """
    pw = pair_sum_panel.copy()
    _must_have_cols(pw, ["Track", "View", "Panel", "n_targets", "mean_centroid_cos", "mean_edist"], "pair_sum_panel")

    pw_long = []
    pw_long.append(
        pw.rename(columns={"n_targets": "N", "mean_centroid_cos": "Value"})[
            ["Track", "View", "Panel", "N", "Value"]
        ].assign(MetricGroup="Pairwise", Direction="NA", Metric="Mean_CentroidCosine", HigherIsBetter=True)
    )
    # In Script10, edist values can be negative and "closer to 0" is better; treat higher as better.
    pw_long.append(
        pw.rename(columns={"n_targets": "N", "mean_edist": "Value"})[
            ["Track", "View", "Panel", "N", "Value"]
        ].assign(MetricGroup="Pairwise", Direction="NA", Metric="Mean_NegEDist", HigherIsBetter=True)
    )
    pw_long_df = pd.concat(pw_long, ignore_index=True)

    rt = ret_sum_panel.copy()
    _must_have_cols(
        rt,
        ["Track", "View", "Panel", "Direction", "n_queries", "mean_mrr", "mean_success", "median_rank"],
        "ret_sum_panel",
    )

    rt_long = []
    rt_long.append(
        rt.rename(columns={"n_queries": "N", "mean_mrr": "Value"})[
            ["Track", "View", "Panel", "Direction", "N", "Value"]
        ].assign(MetricGroup="Retrieval", Metric="Mean_MRR", HigherIsBetter=True)
    )
    rt_long.append(
        rt.rename(columns={"n_queries": "N", "mean_success": "Value"})[
            ["Track", "View", "Panel", "Direction", "N", "Value"]
        ].assign(MetricGroup="Retrieval", Metric="Mean_Success", HigherIsBetter=True)
    )
    rt_long.append(
        rt.rename(columns={"n_queries": "N", "median_rank": "Value"})[
            ["Track", "View", "Panel", "Direction", "N", "Value"]
        ].assign(MetricGroup="Retrieval", Metric="Median_Rank", HigherIsBetter=False)
    )
    rt_long_df = pd.concat(rt_long, ignore_index=True)

    pw_long_df["Direction"] = "NA"
    rt_long_df["Direction"] = rt_long_df["Direction"].astype(str)

    out = pd.concat([pw_long_df, rt_long_df], ignore_index=True)

    def _apply_scale(block: pd.DataFrame) -> pd.DataFrame:
        higher = bool(block["HigherIsBetter"].iloc[0])
        block = block.copy()
        block["Scaled_Best"] = _scale01_clip(block["Value"], higher_is_better=higher)
        block["Rank_Best"] = (-block["Scaled_Best"]).rank(method="min").astype("Int64")
        return block

    group_cols = ["MetricGroup", "Metric", "View", "Panel", "Direction"]
    out = out.groupby(group_cols, dropna=False, group_keys=False).apply(_apply_scale)

    out["Panel"] = pd.Categorical(out["Panel"], categories=DEFAULT_PANELS + ["Control", "Unknown"], ordered=True)
    out = out.sort_values(["MetricGroup", "Metric", "View", "Panel", "Direction", "Rank_Best", "Track"])
    return out


def build_scoreboard_wide(score_long: pd.DataFrame) -> pd.DataFrame:
    df = score_long.copy()

    def _colkey(r) -> str:
        if r["MetricGroup"] == "Pairwise":
            return f"{r['Panel']}__{r['Metric']}"
        return f"{r['Panel']}__{r['Direction']}__{r['Metric']}"

    df["ColKey"] = df.apply(_colkey, axis=1)

    wide_val = (
        df.pivot_table(index=["View", "Track"], columns="ColKey", values="Value", aggfunc="first")
        .reset_index()
    )
    wide_scaled = (
        df.pivot_table(index=["View", "Track"], columns="ColKey", values="Scaled_Best", aggfunc="first")
        .reset_index()
    )
    wide_scaled = wide_scaled.rename(
        columns={c: f"{c}__SCALED" for c in wide_scaled.columns if c not in ["View", "Track"]}
    )
    out = pd.merge(wide_val, wide_scaled, on=["View", "Track"], how="left")
    return out


# ---------------------------
# Figure 2: Query-level paired deltas (viz_stats_deltas_query.csv)
# ---------------------------

def _make_query_delta_long(
    df_joined: pd.DataFrame,
    metric_map: List[Tuple[str, str, bool]],
    contrast: str,
    group_cols_keep: List[str],
    col_suffix_a: str,
    col_suffix_b: str,
    rank_improvement: bool,
) -> pd.DataFrame:
    rows = []
    for col, mname, hib in metric_map:
        a = pd.to_numeric(df_joined[f"{col}{col_suffix_a}"], errors="coerce")
        b = pd.to_numeric(df_joined[f"{col}{col_suffix_b}"], errors="coerce")

        ok = a.notna() & b.notna()
        if ok.sum() == 0:
            continue

        if col == "True_Rank" and rank_improvement:
            # positive means improvement (baseline rank - compared rank)
            delta = (b[ok] - a[ok]).astype(float)
            raw_delta = (a[ok] - b[ok]).astype(float)
            mname2 = "Rank_Improvement"
            hib2 = True
        else:
            delta = (a[ok] - b[ok]).astype(float)
            raw_delta = delta
            mname2 = mname
            hib2 = hib

        base_part = df_joined.loc[ok, group_cols_keep].copy()
        base_part["Contrast"] = contrast
        base_part["Metric"] = mname2
        base_part["HigherIsBetter"] = hib2
        base_part["Delta"] = delta.values
        base_part["Delta_Raw"] = raw_delta.values
        base_part["Value_A"] = a[ok].values
        base_part["Value_B"] = b[ok].values
        rows.append(base_part)

    if not rows:
        return pd.DataFrame()
    out = pd.concat(rows, ignore_index=True)
    return out


def compute_query_deltas_systema_minus_standard(per_q: pd.DataFrame) -> pd.DataFrame:
    """
    Within each (Track, Direction, Query_ID), compute Systema - Standard deltas.
    For Rank, provide Rank_Improvement = Rank_Standard - Rank_Systema (positive=better).
    """
    df = per_q.copy()
    req = ["Track", "View", "Direction", "Query_ID", "True_Target", "specificity_tier",
           "MRR", "Success_Score", "True_Rank"]
    _must_have_cols(df, req, "Retrieval_PerQuery(annotated)")

    df["Query_ID"] = df["Query_ID"].astype(str)
    df["True_Target"] = df["True_Target"].astype(str)

    a = df[df["View"] == "Systema"].copy()
    b = df[df["View"] == "Standard"].copy()
    if a.empty or b.empty:
        warnings.warn("[Deltas] Missing Systema or Standard rows in retrieval per-query. Systema-Standard deltas will be empty.")
        return pd.DataFrame()

    key = ["Track", "Direction", "Query_ID"]
    a = a.set_index(key)
    b = b.set_index(key)
    joined = a.join(b, how="inner", lsuffix="_sys", rsuffix="_std").reset_index()

    metric_map = [
        ("MRR", "MRR", True),
        ("Success_Score", "Success", True),
        ("True_Rank", "Rank", False),  # will become Rank_Improvement
    ]
    keep_cols = ["Track", "Direction", "Query_ID", "True_Target_sys", "specificity_tier_sys"]
    out = _make_query_delta_long(
        joined,
        metric_map=metric_map,
        contrast="Systema_minus_Standard",
        group_cols_keep=keep_cols,
        col_suffix_a="_sys",
        col_suffix_b="_std",
        rank_improvement=True,  # Rank_Improvement = std - sys
    )
    if out.empty:
        return out

    out = out.rename(columns={"True_Target_sys": "True_Target", "specificity_tier_sys": "specificity_tier"})
    out["View"] = "NA"
    out["Baseline"] = "Standard"
    out["Compared"] = "Systema"
    return out


def compute_query_deltas_track_minus_gene(per_q: pd.DataFrame, baseline_track: str = "Gene") -> pd.DataFrame:
    """
    Within each (View, Direction, Query_ID), compute Track - Gene deltas.
    For Rank, provide Rank_Improvement = Rank_Gene - Rank_Track (positive=better).
    """
    df = per_q.copy()
    req = ["Track", "View", "Direction", "Query_ID", "True_Target", "specificity_tier",
           "MRR", "Success_Score", "True_Rank"]
    _must_have_cols(df, req, "Retrieval_PerQuery(annotated)")

    df["Query_ID"] = df["Query_ID"].astype(str)
    df["True_Target"] = df["True_Target"].astype(str)

    out_all = []
    for view in sorted(df["View"].dropna().unique()):
        for direction in sorted(df["Direction"].dropna().unique()):
            sub = df[(df["View"] == view) & (df["Direction"] == direction)]
            base = sub[sub["Track"] == baseline_track].copy()
            if base.empty:
                continue
            base = base.set_index("Query_ID")

            for track in sorted(sub["Track"].dropna().unique()):
                if track == baseline_track:
                    continue
                cur = sub[sub["Track"] == track].copy().set_index("Query_ID")
                joined = cur.join(base, how="inner", lsuffix="_t", rsuffix="_b").reset_index()
                if joined.empty:
                    continue

                metric_map = [
                    ("MRR", "MRR", True),
                    ("Success_Score", "Success", True),
                    ("True_Rank", "Rank", False),  # will become Rank_Improvement
                ]
                keep_cols = ["Query_ID", "True_Target_t", "specificity_tier_t"]
                out = _make_query_delta_long(
                    joined,
                    metric_map=metric_map,
                    contrast="Track_minus_Gene",
                    group_cols_keep=keep_cols,
                    col_suffix_a="_t",
                    col_suffix_b="_b",
                    rank_improvement=True,  # Rank_Improvement = gene - track
                )
                if out.empty:
                    continue

                out = out.rename(columns={"True_Target_t": "True_Target", "specificity_tier_t": "specificity_tier"})
                out["View"] = view
                out["Direction"] = direction
                out["Track"] = track
                out["Baseline"] = baseline_track
                out["Compared"] = track
                out_all.append(out)

    if not out_all:
        return pd.DataFrame()
    return pd.concat(out_all, ignore_index=True)


# ---------------------------
# Figure 2 (NEW): Lollipop-ready summary table (viz_fig2_lollipop.csv)
# ---------------------------

def _summarize_diff_block(
    diffs: np.ndarray,
    metric: str,
    metric_group: str,
    compare_panel: str,
    baseline: str,
    compared: str,
    track: str,
    view: str,
    direction: str,
    rng: np.random.Generator,
    boot_median: int,
    boot_max_n: int,
) -> Dict:
    diffs = np.asarray(diffs, dtype=float)
    diffs = diffs[np.isfinite(diffs)]
    if diffs.size < 3:
        return {}

    p = _safe_wilcoxon_p(diffs)
    return {
        "ComparePanel": compare_panel,
        "MetricGroup": metric_group,
        "Metric": metric,
        "Track": track,
        "View": view,
        "Direction": direction,
        "Baseline": baseline,
        "Compared": compared,
        "N": int(diffs.size),
        "MedianDiff": float(np.median(diffs)),
        "MeanDiff": float(np.mean(diffs)),
        "StdErrMean": _stderr_mean(diffs),
        "StdErrMedian": _stderr_median_bootstrap(diffs, rng=rng, n_boot=int(boot_median), max_n=int(boot_max_n)),
        "P_Wilcoxon": p,
        "NegLog10P": _neglog10_p(p),
    }


def build_fig2_lollipop_summary(
    pair_per_target_anno: pd.DataFrame,
    ret_per_query_anno: pd.DataFrame,
    baseline_track: str = "Gene",
    direction_keep: str = "CRISPR->Drug",
    eps_ed: float = 1e-12,
    boot_median: int = 200,
    boot_max_n: int = 5000,
    seed: int = 0,
) -> pd.DataFrame:
    """
    Build Fig2 lollipop summary table with stderr.

    Output columns:
      ComparePanel, MetricGroup, Metric, Track, View, Direction,
      Baseline, Compared, N, MedianDiff, MeanDiff, StdErrMean, StdErrMedian,
      P_Wilcoxon, NegLog10P

    ComparePanel:
      - View(Systema-Standard): Systema - Standard, within each Track
      - Standard(Track-Gene):   Track - Gene, within Standard view
    """
    rng = np.random.default_rng(int(seed))
    rows: List[Dict] = []

    # ---------- Pairwise (per-target; paired by Target) ----------
    pw = pair_per_target_anno.copy()
    _must_have_cols(pw, ["Track", "View", "Target", "centroid_cosine", "edist_mean"], "annotated_pairwise_per_target")
    pw["Target"] = pw["Target"].astype(str)

    pw["ED_score"] = _edist_to_edscore(pw["edist_mean"], eps=eps_ed)

    # Panel A: View(Systema-Standard)  (Systema - Standard) within Track, by Target pairing
    pw_sys = pw[pw["View"].astype(str) == "Systema"][["Track", "Target", "centroid_cosine", "ED_score"]].copy()
    pw_std = pw[pw["View"].astype(str) == "Standard"][["Track", "Target", "centroid_cosine", "ED_score"]].copy()

    if (not pw_sys.empty) and (not pw_std.empty):
        j = pw_sys.merge(pw_std, on=["Track", "Target"], how="inner", suffixes=("_sys", "_std"))
        for track in sorted(j["Track"].dropna().astype(str).unique()):
            sub = j[j["Track"].astype(str) == track]

            d_cos = (pd.to_numeric(sub["centroid_cosine_sys"], errors="coerce")
                     - pd.to_numeric(sub["centroid_cosine_std"], errors="coerce")).to_numpy(dtype=float)
            d_ed = (pd.to_numeric(sub["ED_score_sys"], errors="coerce")
                    - pd.to_numeric(sub["ED_score_std"], errors="coerce")).to_numpy(dtype=float)

            r1 = _summarize_diff_block(
                d_cos, "CentroidCosine", "Pairwise", "View(Systema-Standard)",
                "Standard", "Systema", track, "NA", "NA",
                rng, boot_median, boot_max_n
            )
            if r1:
                rows.append(r1)

            r2 = _summarize_diff_block(
                d_ed, "ED_score", "Pairwise", "View(Systema-Standard)",
                "Standard", "Systema", track, "NA", "NA",
                rng, boot_median, boot_max_n
            )
            if r2:
                rows.append(r2)
    else:
        warnings.warn("[Fig2] Pairwise: missing Systema or Standard; View(Systema-Standard) pairwise panel may be empty.")

    # Panel B: Standard(Track-Gene) (Track - Gene) within Standard view, by Target pairing
    pw_std_all = pw[pw["View"].astype(str) == "Standard"][["Track", "Target", "centroid_cosine", "ED_score"]].copy()
    base_pw = pw_std_all[pw_std_all["Track"].astype(str) == str(baseline_track)][["Target", "centroid_cosine", "ED_score"]].copy()
    if not base_pw.empty:
        base_pw = base_pw.rename(columns={"centroid_cosine": "cos_b", "ED_score": "ed_b"})
        for track in sorted(pw_std_all["Track"].dropna().astype(str).unique()):
            if track == str(baseline_track):
                continue
            cur = pw_std_all[pw_std_all["Track"].astype(str) == track][["Target", "centroid_cosine", "ED_score"]].copy()
            if cur.empty:
                continue
            cur = cur.rename(columns={"centroid_cosine": "cos_t", "ED_score": "ed_t"})
            j = cur.merge(base_pw, on="Target", how="inner")
            if j.empty:
                continue

            d_cos = (pd.to_numeric(j["cos_t"], errors="coerce") - pd.to_numeric(j["cos_b"], errors="coerce")).to_numpy(dtype=float)
            d_ed = (pd.to_numeric(j["ed_t"], errors="coerce") - pd.to_numeric(j["ed_b"], errors="coerce")).to_numpy(dtype=float)

            r1 = _summarize_diff_block(
                d_cos, "CentroidCosine", "Pairwise", "Standard(Track-Gene)",
                str(baseline_track), track, track, "Standard", "NA",
                rng, boot_median, boot_max_n
            )
            if r1:
                rows.append(r1)

            r2 = _summarize_diff_block(
                d_ed, "ED_score", "Pairwise", "Standard(Track-Gene)",
                str(baseline_track), track, track, "Standard", "NA",
                rng, boot_median, boot_max_n
            )
            if r2:
                rows.append(r2)
    else:
        warnings.warn(f"[Fig2] Pairwise: baseline track '{baseline_track}' missing in Standard view; Standard(Track-Gene) pairwise panel may be empty.")

    # ---------- Retrieval Success (per-query; paired by Query_ID; Direction fixed) ----------
    rq = ret_per_query_anno.copy()
    _must_have_cols(rq, ["Track", "View", "Direction", "Query_ID", "Success_Score"], "annotated_retrieval_per_query")
    rq["Query_ID"] = rq["Query_ID"].astype(str)
    rq["Direction"] = rq["Direction"].astype(str)
    rq["Track"] = rq["Track"].astype(str)
    rq["View"] = rq["View"].astype(str)

    rq = rq[rq["Direction"] == str(direction_keep)].copy()
    if rq.empty:
        warnings.warn(f"[Fig2] Retrieval: no rows for Direction='{direction_keep}'. Success rows will be empty.")
    else:
        # Panel A: View(Systema-Standard)  (Systema - Standard) within Track, paired by Query_ID
        rq_sys = rq[rq["View"] == "Systema"][["Track", "Query_ID", "Success_Score"]].copy()
        rq_std = rq[rq["View"] == "Standard"][["Track", "Query_ID", "Success_Score"]].copy()
        if (not rq_sys.empty) and (not rq_std.empty):
            j = rq_sys.merge(rq_std, on=["Track", "Query_ID"], how="inner", suffixes=("_sys", "_std"))
            for track in sorted(j["Track"].dropna().astype(str).unique()):
                sub = j[j["Track"].astype(str) == track]
                d = (pd.to_numeric(sub["Success_Score_sys"], errors="coerce")
                     - pd.to_numeric(sub["Success_Score_std"], errors="coerce")).to_numpy(dtype=float)

                r = _summarize_diff_block(
                    d, "Success_CRISPR2Drug", "Retrieval", "View(Systema-Standard)",
                    "Standard", "Systema", track, "NA", str(direction_keep),
                    rng, boot_median, boot_max_n
                )
                if r:
                    rows.append(r)
        else:
            warnings.warn("[Fig2] Retrieval: missing Systema or Standard; View(Systema-Standard) retrieval panel may be empty.")

        # Panel B: Standard(Track-Gene) (Track - Gene) within Standard view, paired by Query_ID
        rq_std_all = rq[rq["View"] == "Standard"][["Track", "Query_ID", "Success_Score"]].copy()
        base = rq_std_all[rq_std_all["Track"] == str(baseline_track)][["Query_ID", "Success_Score"]].copy()
        if not base.empty:
            base = base.rename(columns={"Success_Score": "succ_b"})
            for track in sorted(rq_std_all["Track"].dropna().astype(str).unique()):
                if track == str(baseline_track):
                    continue
                cur = rq_std_all[rq_std_all["Track"] == track][["Query_ID", "Success_Score"]].copy()
                if cur.empty:
                    continue
                cur = cur.rename(columns={"Success_Score": "succ_t"})
                j = cur.merge(base, on="Query_ID", how="inner")
                if j.empty:
                    continue
                d = (pd.to_numeric(j["succ_t"], errors="coerce") - pd.to_numeric(j["succ_b"], errors="coerce")).to_numpy(dtype=float)

                r = _summarize_diff_block(
                    d, "Success_CRISPR2Drug", "Retrieval", "Standard(Track-Gene)",
                    str(baseline_track), track, track, "Standard", str(direction_keep),
                    rng, boot_median, boot_max_n
                )
                if r:
                    rows.append(r)
        else:
            warnings.warn(f"[Fig2] Retrieval: baseline track '{baseline_track}' missing in Standard view; Standard(Track-Gene) retrieval panel may be empty.")

    out = pd.DataFrame(rows)
    if out.empty:
        return out

    # enforce column order (stable for R)
    col_order = [
        "ComparePanel", "MetricGroup", "Metric", "Track", "View", "Direction",
        "Baseline", "Compared", "N",
        "MedianDiff", "MeanDiff", "StdErrMean", "StdErrMedian",
        "P_Wilcoxon", "NegLog10P",
    ]
    for c in col_order:
        if c not in out.columns:
            out[c] = np.nan
    out = out[col_order]

    # stable ordering hints
    metric_order = ["CentroidCosine", "ED_score", "Success_CRISPR2Drug"]
    out["Metric"] = pd.Categorical(out["Metric"].astype(str), categories=metric_order, ordered=True)
    out["ComparePanel"] = pd.Categorical(
        out["ComparePanel"].astype(str),
        categories=["View(Systema-Standard)", "Standard(Track-Gene)"],
        ordered=True
    )

    out = out.sort_values(["ComparePanel", "Track", "MetricGroup", "Metric"])
    return out


# ---------------------------
# Figure 3: Target-level gains vs Gene
# ---------------------------

def compute_target_gain_table(
    per_q: pd.DataFrame,
    per_tgt: pd.DataFrame,
    target_tier_map: pd.DataFrame,
    baseline_track: str = "Gene",
) -> pd.DataFrame:
    """
    Target-level deltas vs Gene baseline.
    Retrieval aggregated per target (mean MRR, mean Success, median Rank + Rank_Improvement)
    Pairwise already target-level (centroid_cosine, edist_mean)
    """
    tgt2tier = {}
    if target_tier_map is not None and (not target_tier_map.empty):
        tgt2tier = dict(zip(target_tier_map["Target"].astype(str), target_tier_map["Target_Tier"].astype(str)))

    # --- Retrieval: aggregate to target-level ---
    rq = per_q.copy()
    _must_have_cols(
        rq, ["Track", "View", "Direction", "True_Target", "Query_ID", "MRR", "Success_Score", "True_Rank"],
        "Retrieval_PerQuery(annotated)"
    )
    rq["True_Target"] = rq["True_Target"].astype(str)

    ret_tgt = (
        rq.groupby(["Track", "View", "Direction", "True_Target"], dropna=False)
        .agg(
            n_queries=("Query_ID", "count"),
            mean_mrr=("MRR", "mean"),
            mean_success=("Success_Score", "mean"),
            median_rank=("True_Rank", "median"),
        )
        .reset_index()
        .rename(columns={"True_Target": "Target"})
    )
    ret_tgt["Target_Tier"] = ret_tgt["Target"].map(lambda t: _canonicalize_tier(tgt2tier.get(str(t), "Unknown")))

    # compute deltas vs Gene
    ret_rows = []
    for (view, direction, target), sub in ret_tgt.groupby(["View", "Direction", "Target"], dropna=False):
        base = sub[sub["Track"].astype(str) == str(baseline_track)]
        if base.empty:
            continue
        base = base.iloc[0]

        for _, r in sub.iterrows():
            track = str(r["Track"])
            if track == str(baseline_track):
                continue

            ret_rows.append({
                "MetricGroup": "Retrieval",
                "Metric": "Mean_MRR",
                "View": view,
                "Direction": direction,
                "Target": target,
                "Target_Tier": r["Target_Tier"],
                "Track": track,
                "BaselineTrack": str(baseline_track),
                "N": int(r["n_queries"]),
                "Value": float(r["mean_mrr"]),
                "BaselineValue": float(base["mean_mrr"]),
                "Delta": float(r["mean_mrr"] - base["mean_mrr"]),
                "Delta_Better": float(r["mean_mrr"] - base["mean_mrr"]),
                "HigherIsBetter": True,
            })
            ret_rows.append({
                "MetricGroup": "Retrieval",
                "Metric": "Mean_Success",
                "View": view,
                "Direction": direction,
                "Target": target,
                "Target_Tier": r["Target_Tier"],
                "Track": track,
                "BaselineTrack": str(baseline_track),
                "N": int(r["n_queries"]),
                "Value": float(r["mean_success"]),
                "BaselineValue": float(base["mean_success"]),
                "Delta": float(r["mean_success"] - base["mean_success"]),
                "Delta_Better": float(r["mean_success"] - base["mean_success"]),
                "HigherIsBetter": True,
            })
            # Rank: smaller better; improvement = base - track
            rank_impr = float(base["median_rank"] - r["median_rank"])
            ret_rows.append({
                "MetricGroup": "Retrieval",
                "Metric": "Rank_Improvement",
                "View": view,
                "Direction": direction,
                "Target": target,
                "Target_Tier": r["Target_Tier"],
                "Track": track,
                "BaselineTrack": str(baseline_track),
                "N": int(r["n_queries"]),
                "Value": float(r["median_rank"]),
                "BaselineValue": float(base["median_rank"]),
                "Delta": float(r["median_rank"] - base["median_rank"]),  # raw delta
                "Delta_Better": rank_impr,
                "HigherIsBetter": True,
            })

    ret_gain = pd.DataFrame(ret_rows)

    # --- Pairwise: already target-level ---
    pw = per_tgt.copy()
    _must_have_cols(pw, ["Track", "View", "Target", "centroid_cosine", "edist_mean"], "Pairwise_PerTarget(annotated)")
    pw["Target"] = pw["Target"].astype(str)
    pw["Target_Tier"] = pw["Target"].map(lambda t: _canonicalize_tier(tgt2tier.get(str(t), "Unknown")))

    pw_rows = []
    for (view, target), sub in pw.groupby(["View", "Target"], dropna=False):
        base = sub[sub["Track"].astype(str) == str(baseline_track)]
        if base.empty:
            continue
        base = base.iloc[0]

        for _, r in sub.iterrows():
            track = str(r["Track"])
            if track == str(baseline_track):
                continue

            d_cos = float(pd.to_numeric(r["centroid_cosine"], errors="coerce") - pd.to_numeric(base["centroid_cosine"], errors="coerce"))
            d_ed = float(pd.to_numeric(r["edist_mean"], errors="coerce") - pd.to_numeric(base["edist_mean"], errors="coerce"))

            pw_rows.append({
                "MetricGroup": "Pairwise",
                "Metric": "CentroidCosine",
                "View": view,
                "Direction": "NA",
                "Target": target,
                "Target_Tier": r["Target_Tier"],
                "Track": track,
                "BaselineTrack": str(baseline_track),
                "N": 0,
                "Value": float(pd.to_numeric(r["centroid_cosine"], errors="coerce")),
                "BaselineValue": float(pd.to_numeric(base["centroid_cosine"], errors="coerce")),
                "Delta": d_cos,
                "Delta_Better": d_cos,
                "HigherIsBetter": True,
            })
            pw_rows.append({
                "MetricGroup": "Pairwise",
                "Metric": "NegEDist",
                "View": view,
                "Direction": "NA",
                "Target": target,
                "Target_Tier": r["Target_Tier"],
                "Track": track,
                "BaselineTrack": str(baseline_track),
                "N": 0,
                "Value": float(pd.to_numeric(r["edist_mean"], errors="coerce")),
                "BaselineValue": float(pd.to_numeric(base["edist_mean"], errors="coerce")),
                "Delta": d_ed,
                "Delta_Better": d_ed,
                "HigherIsBetter": True,
            })

    pw_gain = pd.DataFrame(pw_rows)

    out = pd.concat([ret_gain, pw_gain], ignore_index=True)
    if out.empty:
        return out

    out["Target_Tier"] = out["Target_Tier"].map(_canonicalize_tier)
    out["Target_Tier"] = pd.Categorical(out["Target_Tier"], categories=TIER_ORDER, ordered=True)
    out = out.sort_values(["MetricGroup", "Metric", "View", "Direction", "Target_Tier", "Target", "Track"])
    return out


# ---------------------------
# Main
# ---------------------------

def main():
    ap = argparse.ArgumentParser("Chem2Gen Task3 — Visualization CSV preparation (Script 11)")

    ap.add_argument("--eval_dir", type=str,
                    default="/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets/Evaluation_Set_K562")
    ap.add_argument("--analysis_dir", type=str,
                    default="/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/Task3_Analysis")
    ap.add_argument("--out_dir", type=str,
                    default="/mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task3_Unified")

    ap.add_argument("--panels", type=str, default=",".join(DEFAULT_PANELS),
                    help="Comma-separated panels. Default: All + Cleanest + Family + Promiscuous")
    ap.add_argument("--baseline_track", type=str, default="Gene",
                    help="Baseline track for deltas (Track-vs-Gene).")

    # Fig2 lollipop controls
    ap.add_argument("--fig2_direction", type=str, default="CRISPR->Drug",
                    help="Direction used for retrieval Success lollipop (default: CRISPR->Drug).")
    ap.add_argument("--fig2_eps_ed", type=float, default=1e-12,
                    help="Epsilon for ED_score = -log10(max(-edist_mean, eps)).")
    ap.add_argument("--fig2_boot_median", type=int, default=200,
                    help="Bootstrap replicates for StdErrMedian in Fig2 lollipop summary.")
    ap.add_argument("--fig2_boot_max_n", type=int, default=5000,
                    help="Max N used for bootstrap median SE per group (subsample then scale back).")
    ap.add_argument("--fig2_seed", type=int, default=0,
                    help="Random seed for bootstrap SE.")

    args = ap.parse_args()

    panels = [p.strip() for p in args.panels.split(",") if p.strip()]
    if "All" not in panels:
        panels = ["All"] + panels

    _ensure_dir(args.out_dir)

    # ---- Required Script10 outputs ----
    f_pair_per_target = os.path.join(args.analysis_dir, "Task3_Pairwise_PerTarget.csv")
    f_pair_summary = os.path.join(args.analysis_dir, "Task3_Pairwise_Summary.csv")
    f_ret_per_query = os.path.join(args.analysis_dir, "Task3_Retrieval_PerQuery.csv")
    f_ret_summary = os.path.join(args.analysis_dir, "Task3_Retrieval_Summary.csv")

    pair_per_target = _read_required_csv(f_pair_per_target, "Task3_Pairwise_PerTarget.csv")
    _ = _read_required_csv(f_pair_summary, "Task3_Pairwise_Summary.csv")   # traceability only
    ret_per_query = _read_required_csv(f_ret_per_query, "Task3_Retrieval_PerQuery.csv")
    _ = _read_required_csv(f_ret_summary, "Task3_Retrieval_Summary.csv")   # traceability only

    _must_have_cols(pair_per_target, ["Track", "View", "Target", "centroid_cosine", "edist_mean"], "Task3_Pairwise_PerTarget")
    _must_have_cols(ret_per_query, ["Track", "View", "Direction", "Query_Modality", "Query_ID",
                                    "True_Target", "True_Rank", "MRR", "Success_Score"], "Task3_Retrieval_PerQuery")

    # ---- Optional meta for tier annotation ----
    drug_meta_path = os.path.join(args.eval_dir, "Drug_meta.csv")
    drug_meta = _read_optional_csv(drug_meta_path)

    if drug_meta.empty:
        warnings.warn(f"[Meta] Drug_meta not found or empty at: {drug_meta_path}. Tier panels will be 'Unknown'.")
    else:
        if "specificity_tier" in drug_meta.columns:
            drug_meta["specificity_tier"] = drug_meta["specificity_tier"].map(_canonicalize_tier)
        drug_meta.index = drug_meta.index.astype(str)

    target_tier_map = build_target_tier_mapping(drug_meta)

    # ---- Annotate per-target / per-query with tiers ----
    pair_per_target_anno = attach_tiers_to_pairwise_per_target(pair_per_target, target_tier_map)
    ret_per_query_anno = attach_tiers_to_retrieval_per_query(ret_per_query, drug_meta, target_tier_map)

    pair_per_target_anno.to_csv(os.path.join(args.out_dir, "annotated_pairwise_per_target.csv"), index=False)
    ret_per_query_anno.to_csv(os.path.join(args.out_dir, "annotated_retrieval_per_query.csv"), index=False)

    # =========================
    # Figure 1: Scoreboard
    # =========================
    pair_panel_sum = summarize_pairwise_by_panel(pair_per_target_anno, panels=panels)
    ret_panel_sum = summarize_retrieval_by_panel(ret_per_query_anno, panels=panels)

    scoreboard_long = build_scoreboard_long(pair_sum_panel=pair_panel_sum, ret_sum_panel=ret_panel_sum)
    scoreboard_wide = build_scoreboard_wide(scoreboard_long)

    scoreboard_long.to_csv(os.path.join(args.out_dir, "viz_scoreboard_long.csv"), index=False)
    scoreboard_wide.to_csv(os.path.join(args.out_dir, "viz_scoreboard_wide.csv"), index=False)

    # =========================
    # Figure 2a: Query-level deltas (existing)
    # =========================
    deltas_sys = compute_query_deltas_systema_minus_standard(ret_per_query_anno)
    deltas_track = compute_query_deltas_track_minus_gene(ret_per_query_anno, baseline_track=args.baseline_track)

    deltas_all = pd.concat([deltas_sys, deltas_track], ignore_index=True)
    if deltas_all.empty:
        warnings.warn("[viz_stats_deltas_query] No deltas produced (check if both views exist, and Gene baseline exists).")
        deltas_all.to_csv(os.path.join(args.out_dir, "viz_stats_deltas_query.csv"), index=False)
    else:
        deltas_all["specificity_tier"] = deltas_all["specificity_tier"].map(_canonicalize_tier)
        deltas_all["specificity_tier"] = pd.Categorical(deltas_all["specificity_tier"], categories=TIER_ORDER, ordered=True)
        deltas_all["Direction"] = deltas_all.get("Direction", "NA").astype(str)
        deltas_all = deltas_all.rename(columns={"True_Target": "Target"})
        deltas_all = deltas_all[
            ["Contrast", "Metric", "HigherIsBetter", "Delta", "Delta_Raw",
             "Value_A", "Value_B",
             "Track", "Baseline", "Compared",
             "View", "Direction", "Query_ID", "Target", "specificity_tier"]
        ]
        deltas_all.to_csv(os.path.join(args.out_dir, "viz_stats_deltas_query.csv"), index=False)

    # =========================
    # Figure 2b: Lollipop-ready summary (WITH StdErrMean/StdErrMedian)
    # =========================
    fig2 = build_fig2_lollipop_summary(
        pair_per_target_anno=pair_per_target_anno,
        ret_per_query_anno=ret_per_query_anno,
        baseline_track=args.baseline_track,
        direction_keep=args.fig2_direction,
        eps_ed=float(args.fig2_eps_ed),
        boot_median=int(args.fig2_boot_median),
        boot_max_n=int(args.fig2_boot_max_n),
        seed=int(args.fig2_seed),
    )
    fig2_path = os.path.join(args.out_dir, "viz_fig2_lollipop.csv")
    if fig2.empty:
        warnings.warn("[viz_fig2_lollipop] Empty. Check Systema/Standard exist, baseline exists, and fig2_direction exists.")
    fig2.to_csv(fig2_path, index=False)

    # =========================
    # Figure 3: Target-level gains vs Gene
    # =========================
    target_gain = compute_target_gain_table(
        per_q=ret_per_query_anno,
        per_tgt=pair_per_target_anno,
        target_tier_map=target_tier_map,
        baseline_track=args.baseline_track,
    )
    target_gain.to_csv(os.path.join(args.out_dir, "viz_target_gain.csv"), index=False)

    # ---- Done ----
    print("\n[Task3 Script11] Done. Wrote R-ready CSVs to:")
    print(f"  {args.out_dir}\n")
    print("Key outputs:")
    print(f"  - {os.path.join(args.out_dir, 'viz_scoreboard_long.csv')}")
    print(f"  - {os.path.join(args.out_dir, 'viz_stats_deltas_query.csv')}")
    print(f"  - {os.path.join(args.out_dir, 'viz_fig2_lollipop.csv')}")
    print(f"  - {os.path.join(args.out_dir, 'viz_target_gain.csv')}")
    print(f"  - {os.path.join(args.out_dir, 'annotated_retrieval_per_query.csv')}")
    print(f"  - {os.path.join(args.out_dir, 'annotated_pairwise_per_target.csv')}")


if __name__ == "__main__":
    main()
