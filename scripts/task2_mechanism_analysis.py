#!/usr/bin/env python3
"""
M2M-Bench Task2 (Mechanism Fidelity) - refactored analysis script.

Design goals
------------
1) Data-analysis script style (single script, explicit sequential steps).
2) High readability: clear section blocks, typed helper functions, concise comments.
3) No DomainType main comparison:
   - We keep source fields for diagnostics.
   - Main tests/labels do not use SameDomain vs CrossDomain as primary contrast.

Input requirements
------------------
- Pairwise wide table from legacy Task1 pairwise module.
- Retrieval per-query table from legacy Task1 multi-scenario retrieval module.

This script does not rebuild raw embeddings; it is an analysis-layer consolidation.
"""

from __future__ import annotations

import argparse
import json
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact, spearmanr, wilcoxon


# =====================================================================
# Section 0. Configuration
# =====================================================================


@dataclass
class Task2NoDomainConfig:
    pairwise_wide_path: str = "/mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task1_Metrics/Task1_Pairwise_Metrics_Wide.csv"
    retrieval_per_query_path: str = "/mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task1_Retrieval_MultiScenario/Task1_Retrieval_MultiScenario_PerQuery.csv"
    output_dir: str = "./outputs/task2_nodomain"
    min_query_samples_per_context: int = 3
    min_contexts_for_enrichment: int = 10
    topk_cases_per_class: int = 20
    hero_track: str = "Gene"
    hero_view: str = "Systema"
    domain_scope: str = "same_only"
    random_seed: int = 42

    @property
    def analysis_dir(self) -> str:
        return str(Path(self.output_dir) / "analysis")

    @property
    def run_manifest_path(self) -> str:
        return str(Path(self.output_dir) / "run_manifest_task2_nodomain.json")


# =====================================================================
# Section 1. Generic helpers
# =====================================================================


def _ensure_dir(path: str) -> None:
    Path(path).mkdir(parents=True, exist_ok=True)


def _safe_wilcoxon(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """Paired Wilcoxon with robust handling for small N / all ties."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if x.size < 2:
        return np.nan, np.nan
    d = x - y
    if np.allclose(d, 0):
        return 0.0, 1.0
    stat, pval = wilcoxon(x, y)
    return float(stat), float(pval)


def _bh_fdr(p_values: np.ndarray) -> np.ndarray:
    """Benjamini-Hochberg FDR; NaN p-values treated as 1."""
    p = np.asarray(p_values, dtype=float)
    if p.size == 0:
        return p
    p = np.where(np.isfinite(p), p, 1.0)
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    corr = ranked * (n / np.arange(1, n + 1))
    corr = np.minimum.accumulate(corr[::-1])[::-1]
    out = np.empty_like(corr)
    out[order] = np.clip(corr, 0, 1)
    return out


def _calc_stable_log2_or(a: int, b: int, c: int, d: int) -> float:
    """Haldane-Anscombe corrected log2 odds ratio."""
    num = (a + 0.5) * (d + 0.5)
    den = (b + 0.5) * (c + 0.5)
    return float(np.log2(num / den))


# =====================================================================
# Section 2. Input loading + standardization
# =====================================================================


def load_pairwise_wide(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)

    required = [
        "Track",
        "cell_std",
        "target_std",
        "source_db_chem",
        "source_db_gene",
        "cosine_std",
        "cosine_sys",
        "systema_gain_cosine",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Pairwise wide missing columns: {missing}")

    out = df.rename(
        columns={
            "cell_std": "Cell",
            "target_std": "Target",
            "source_db_chem": "Source_Chem",
            "source_db_gene": "Source_Gene",
        }
    )
    keep = [
        "Track",
        "Cell",
        "Target",
        "Source_Chem",
        "Source_Gene",
        "cosine_std",
        "cosine_sys",
        "systema_gain_cosine",
    ]
    return out[keep].copy()


def load_retrieval(path: Path) -> pd.DataFrame:
    if not path.exists():
        raise FileNotFoundError(path)
    df = pd.read_csv(path)

    required = [
        "Track",
        "Scenario",
        "View",
        "Direction",
        "LabelType",
        "Cell",
        "Target",
        "Source_Chem",
        "Source_Gene",
        "Dose",
        "Time",
        "CondID",
        "True_Rank",
        "MRR",
        "Success_Score",
        "N_Gallery",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Retrieval per-query missing columns: {missing}")

    out = df.copy()
    out["DomainType"] = np.where(
        out["Source_Chem"].astype(str) == out["Source_Gene"].astype(str),
        "SameDomain",
        "CrossDomain",
    )
    out["Dose_num"] = pd.to_numeric(out["Dose"], errors="coerce")
    out["Time_num"] = pd.to_numeric(out["Time"], errors="coerce")
    out["True_Rank"] = pd.to_numeric(out["True_Rank"], errors="coerce")
    out["N_Gallery"] = pd.to_numeric(out["N_Gallery"], errors="coerce")
    out["LogDose"] = np.where(out["Dose_num"] > 0, np.log10(out["Dose_num"]), np.nan)
    return out


# =====================================================================
# Section 3. Step1 - L1 tidy + context + tests (no domain contrast)
# =====================================================================


def attach_pairwise_scores(l1: pd.DataFrame, pairwise: pd.DataFrame) -> pd.DataFrame:
    """
    Attach pairwise cosine to each L1 instance by:
    (Track, Cell, Target, Source_Chem, Source_Gene, View)
    """
    keys = ["Track", "Cell", "Target", "Source_Chem", "Source_Gene"]
    grouped = pairwise.groupby(keys, as_index=False).agg(
        cosine_std=("cosine_std", "mean"),
        cosine_sys=("cosine_sys", "mean"),
        systema_gain_cosine=("systema_gain_cosine", "mean"),
    )

    std = grouped.copy()
    std["View"] = "Standard"
    std["Pairwise_Cosine"] = std["cosine_std"]
    std["Pairwise_Systema_Gain"] = 0.0

    sys = grouped.copy()
    sys["View"] = "Systema"
    sys["Pairwise_Cosine"] = sys["cosine_sys"]
    sys["Pairwise_Systema_Gain"] = sys["systema_gain_cosine"]

    attach = pd.concat([std, sys], ignore_index=True)[keys + ["View", "Pairwise_Cosine", "Pairwise_Systema_Gain"]]
    return l1.merge(attach, on=keys + ["View"], how="left", validate="m:1")


def step1_build_l1(cfg: Task2NoDomainConfig, retrieval: pd.DataFrame, pairwise: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Primary slice definition (unchanged from legacy):
      Scenario=A_L1Base, LabelType=Target, Direction=Chem->Gene
    """
    l1 = retrieval[
        (retrieval["Scenario"] == "A_L1Base")
        & (retrieval["LabelType"] == "Target")
        & (retrieval["Direction"] == "Chem->Gene")
    ].copy()
    if cfg.domain_scope == "same_only":
        l1 = l1[l1["DomainType"] == "SameDomain"].copy()
    if l1.empty:
        raise ValueError("No rows match A_L1Base + Target + Chem->Gene.")

    l1 = attach_pairwise_scores(l1=l1, pairwise=pairwise)

    # Instance-level tidy
    inst_cols = [
        "Track",
        "Scenario",
        "View",
        "Direction",
        "DomainType",
        "Cell",
        "Target",
        "Source_Chem",
        "Source_Gene",
        "Dose",
        "Time",
        "CondID",
        "LogDose",
        "True_Rank",
        "MRR",
        "Success_Score",
        "N_Gallery",
        "Pairwise_Cosine",
        "Pairwise_Systema_Gain",
    ]
    inst = l1[inst_cols].copy()

    # Context-level table (no DomainType key)
    ctx_key = ["Track", "View", "Cell", "Target", "Source_Chem", "Source_Gene"]
    ctx = inst.groupby(ctx_key, as_index=False).agg(
        N_Instances=("Success_Score", "size"),
        Mean_Success=("Success_Score", "mean"),
        Median_Success=("Success_Score", "median"),
        Max_Success=("Success_Score", "max"),
        Mean_MRR=("MRR", "mean"),
        Pairwise_Cosine=("Pairwise_Cosine", "mean"),
    )
    ctx = ctx[ctx["N_Instances"] >= cfg.min_query_samples_per_context].reset_index(drop=True)
    return inst, ctx


def step1_tests_no_domain(ctx: pd.DataFrame) -> pd.DataFrame:
    """
    Keep only mechanism-relevant paired tests:
      1) Standard vs Systema (paired)
      2) Pathway vs Gene (paired)
    """
    rows: List[Dict] = []

    # Test A: Systema vs Standard for each track
    for track in ["Gene", "Pathway"]:
        sub = ctx[ctx["Track"] == track].copy()
        piv = sub.pivot_table(
            index=["Cell", "Target", "Source_Chem", "Source_Gene"],
            columns="View",
            values="Mean_Success",
            aggfunc="mean",
        ).dropna()

        if {"Standard", "Systema"}.issubset(piv.columns) and len(piv) >= 2:
            _, pval = _safe_wilcoxon(piv["Systema"].to_numpy(), piv["Standard"].to_numpy())
            diff = piv["Systema"].to_numpy() - piv["Standard"].to_numpy()
            rows.append(
                {
                    "Comparison": f"Systema vs Standard [Track={track}]",
                    "Test": "Paired_Wilcoxon",
                    "N": int(len(piv)),
                    "P_Value": float(pval),
                    "Median_Diff": float(np.nanmedian(diff)),
                    "Mean_Diff": float(np.nanmean(diff)),
                    "Frac_Positive": float(np.mean(diff > 0)),
                }
            )

    # Test B: Pathway vs Gene for each view
    for view in ["Standard", "Systema"]:
        sub = ctx[ctx["View"] == view].copy()
        piv = sub.pivot_table(
            index=["Cell", "Target", "Source_Chem", "Source_Gene"],
            columns="Track",
            values="Mean_Success",
            aggfunc="mean",
        ).dropna()

        if {"Gene", "Pathway"}.issubset(piv.columns) and len(piv) >= 2:
            _, pval = _safe_wilcoxon(piv["Pathway"].to_numpy(), piv["Gene"].to_numpy())
            diff = piv["Pathway"].to_numpy() - piv["Gene"].to_numpy()
            rows.append(
                {
                    "Comparison": f"Pathway vs Gene [View={view}]",
                    "Test": "Paired_Wilcoxon",
                    "N": int(len(piv)),
                    "P_Value": float(pval),
                    "Median_Diff": float(np.nanmedian(diff)),
                    "Mean_Diff": float(np.nanmean(diff)),
                    "Frac_Positive": float(np.mean(diff > 0)),
                }
            )

    out = pd.DataFrame(rows)
    if not out.empty:
        out["FDR_BH"] = _bh_fdr(out["P_Value"].to_numpy())
    return out


# =====================================================================
# Section 4. Step2 - Context labels + diagnostics (no domain)
# =====================================================================


def _safe_tertile(series: pd.Series) -> pd.Series:
    """Tertile labels with fallback when data are insufficient."""
    s = pd.to_numeric(series, errors="coerce")
    if s.notna().sum() < 3:
        return pd.Series(["Mid"] * len(series), index=series.index)
    ranked = s.rank(method="first")
    try:
        return pd.qcut(ranked, 3, labels=["Low", "Mid", "High"])
    except Exception:
        return pd.Series(["Mid"] * len(series), index=series.index)


def _classify_context(mean_tier: str, peak_tier: str) -> str:
    if str(peak_tier) == "Low":
        return "Robust_Low"
    if str(peak_tier) == "High" and str(mean_tier) == "High":
        return "Robust_High"
    if str(peak_tier) == "High" and str(mean_tier) != "High":
        return "Protocol_Sensitive"
    return "Intermediate"


def step2_context_labeling_no_domain(
    cfg: Task2NoDomainConfig,
    inst_l1: pd.DataFrame,
    retrieval_all: pd.DataFrame,
) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Build labels by hero slice (Track/View/Direction only; no DomainType condition).
    """
    hero = inst_l1[
        (inst_l1["Track"] == cfg.hero_track)
        & (inst_l1["View"] == cfg.hero_view)
        & (inst_l1["Direction"] == "Chem->Gene")
    ].copy()
    if hero.empty:
        raise ValueError("Hero slice is empty.")

    labels = hero.groupby(["Cell", "Target"], as_index=False).agg(
        Mean_Success=("Success_Score", "mean"),
        Peak_Success=("Success_Score", "max"),
        N_Instances=("Success_Score", "size"),
        Context_Cosine=("Pairwise_Cosine", "mean"),
    )

    labels["Mean_Tier"] = _safe_tertile(labels["Mean_Success"])
    labels["Peak_Tier"] = _safe_tertile(labels["Peak_Success"])
    labels["Performance_Class"] = [
        _classify_context(m, p) for m, p in zip(labels["Mean_Tier"], labels["Peak_Tier"])
    ]

    # L2/L3 diagnostics aggregated by source pair (keep source visibility)
    diag_frames: List[pd.DataFrame] = []
    for scenario, prefix in [("B_CellDep", "L2_CellDep"), ("C_TgtSpec", "L3_TgtSpec")]:
        sub = retrieval_all[
            (retrieval_all["Scenario"] == scenario)
            & (retrieval_all["Track"] == cfg.hero_track)
            & (retrieval_all["Direction"] == "Chem->Gene")
        ].copy()
        if cfg.domain_scope == "same_only":
            sub = sub[sub["DomainType"] == "SameDomain"].copy()
        if sub.empty:
            continue

        g = sub.groupby(["Cell", "Target", "Source_Chem", "Source_Gene", "View"], as_index=False).agg(
            MeanSuccess=("Success_Score", "mean"),
            MeanMRR=("MRR", "mean"),
            N=("Success_Score", "size"),
        )
        g = g.rename(
            columns={
                "MeanSuccess": f"{prefix}_MeanSuccess",
                "MeanMRR": f"{prefix}_MeanMRR",
                "N": f"{prefix}_N",
            }
        )
        g["Diag"] = prefix
        diag_frames.append(g)

    diag_sp = pd.concat(diag_frames, ignore_index=True) if diag_frames else pd.DataFrame()
    return labels, diag_sp


# =====================================================================
# Section 5. Step4 + Step5
# =====================================================================


def step4_case_tracer(
    cfg: Task2NoDomainConfig,
    labels: pd.DataFrame,
    retrieval_all: pd.DataFrame,
) -> pd.DataFrame:
    """
    Select top-K contexts per class and trace across A/B/C.
    """
    if labels.empty:
        return pd.DataFrame()

    labels2 = labels.copy()
    labels2["SensitivityScore"] = labels2["Peak_Success"] - labels2["Mean_Success"]

    selected: List[pd.DataFrame] = []
    k = int(max(1, cfg.topk_cases_per_class))

    for cls in ["Robust_High", "Protocol_Sensitive", "Robust_Low"]:
        sub = labels2[labels2["Performance_Class"] == cls].copy()
        if sub.empty:
            continue
        if cls == "Robust_High":
            sub = sub.sort_values("Mean_Success", ascending=False)
        elif cls == "Protocol_Sensitive":
            sub = sub.sort_values("SensitivityScore", ascending=False)
        else:
            sub = sub.sort_values("Mean_Success", ascending=True)
        selected.append(sub.head(k)[["Cell", "Target", "Performance_Class"]])

    if not selected:
        return pd.DataFrame()
    keep = pd.concat(selected, ignore_index=True).drop_duplicates(["Cell", "Target"], keep="first")

    df = retrieval_all.copy()
    df = df[df["Scenario"].isin(["A_L1Base", "B_CellDep", "C_TgtSpec"])].copy()
    if cfg.domain_scope == "same_only":
        df = df[df["DomainType"] == "SameDomain"].copy()
    df = df.merge(keep, on=["Cell", "Target"], how="inner")
    if df.empty:
        return pd.DataFrame()

    # For each scenario, prefer LabelType=Target if present.
    if "LabelType" in df.columns:
        parts = []
        for sc in ["A_L1Base", "B_CellDep", "C_TgtSpec"]:
            s = df[df["Scenario"] == sc].copy()
            if s.empty:
                continue
            st = s[s["LabelType"] == "Target"].copy()
            parts.append(st if not st.empty else s)
        df = pd.concat(parts, ignore_index=True) if parts else df

    tracer = df.groupby(
        ["Cell", "Target", "Performance_Class", "Scenario", "Track", "View"],
        as_index=False,
    ).agg(
        MeanSuccess=("Success_Score", "mean"),
        MeanMRR=("MRR", "mean"),
        N=("Success_Score", "size"),
    )
    return tracer


def run_fisher_enrichment(
    labels: pd.DataFrame,
    group_col: str,
    class_col: str,
    positive_class: str,
    min_contexts: int,
) -> pd.DataFrame:
    """
    One-sided Fisher enrichment for class membership.
    """
    if labels.empty:
        return pd.DataFrame()

    total_pos = int((labels[class_col] == positive_class).sum())
    total_neg = int((labels[class_col] != positive_class).sum())
    total = total_pos + total_neg
    if total == 0 or total_pos == 0 or total_neg == 0:
        return pd.DataFrame()

    baseline = total_pos / total
    rows = []
    for val in labels[group_col].dropna().astype(str).unique():
        is_g = labels[group_col].astype(str) == str(val)
        n_ctx = int(is_g.sum())
        if n_ctx < min_contexts:
            continue
        a = int((is_g & (labels[class_col] == positive_class)).sum())
        c = int(n_ctx - a)
        b = int(total_pos - a)
        d = int(total_neg - c)
        odds, pval = fisher_exact([[a, b], [c, d]], alternative="greater")
        rows.append(
            {
                "Feature": str(val),
                "GroupCol": group_col,
                "Target_Class": positive_class,
                "N_Contexts": n_ctx,
                "N_Positive": a,
                "Baseline_Rate": baseline,
                "Rate_Group": a / n_ctx if n_ctx > 0 else np.nan,
                "Rate_Background": b / (b + d) if (b + d) > 0 else np.nan,
                "Rate_Diff": (a / n_ctx) - (b / (b + d)) if (n_ctx > 0 and (b + d) > 0) else np.nan,
                "Odds_Ratio_Raw": float(odds) if np.isfinite(odds) else np.nan,
                "Log2_OR_Stable": _calc_stable_log2_or(a, b, c, d),
                "P_Value": float(pval),
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["FDR_BH"] = _bh_fdr(out["P_Value"].to_numpy())
    out = out.sort_values(["FDR_BH", "P_Value"], ascending=[True, True]).reset_index(drop=True)
    return out


def protocol_correlations(
    cfg: Task2NoDomainConfig,
    inst_l1: pd.DataFrame,
    labels: pd.DataFrame,
    focus_class: str = "Protocol_Sensitive",
) -> pd.DataFrame:
    """
    Correlate Success_Score vs LogDose/Time within each (Cell,Target,SourcePair).
    Uses hero Track/View/Direction but no DomainType filter.
    """
    lb = labels[["Cell", "Target", "Performance_Class"]].copy()
    sub = inst_l1.merge(lb, on=["Cell", "Target"], how="inner")
    sub = sub[
        (sub["Performance_Class"] == focus_class)
        & (sub["Track"] == cfg.hero_track)
        & (sub["View"] == cfg.hero_view)
        & (sub["Direction"] == "Chem->Gene")
    ].copy()
    if cfg.domain_scope == "same_only":
        sub = sub[sub["DomainType"] == "SameDomain"].copy()
    if sub.empty:
        return pd.DataFrame()

    rows = []
    for (cell, tgt, src_c, src_g), g in sub.groupby(["Cell", "Target", "Source_Chem", "Source_Gene"], sort=False):
        if len(g) < 6:
            continue

        dose_corr, dose_p = np.nan, np.nan
        if pd.to_numeric(g["LogDose"], errors="coerce").nunique(dropna=True) >= 3:
            dose_corr, dose_p = spearmanr(g["LogDose"], g["Success_Score"], nan_policy="omit")

        time_corr, time_p = np.nan, np.nan
        time_num = pd.to_numeric(g["Time"], errors="coerce")
        if time_num.nunique(dropna=True) >= 3:
            time_corr, time_p = spearmanr(time_num, g["Success_Score"], nan_policy="omit")

        rows.append(
            {
                "Cell": cell,
                "Target": tgt,
                "Source_Chem": src_c,
                "Source_Gene": src_g,
                "N": int(len(g)),
                "Dose_Corr": float(dose_corr) if np.isfinite(dose_corr) else np.nan,
                "Dose_P": float(dose_p) if np.isfinite(dose_p) else np.nan,
                "Time_Corr": float(time_corr) if np.isfinite(time_corr) else np.nan,
                "Time_P": float(time_p) if np.isfinite(time_p) else np.nan,
            }
        )
    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["Dose_FDR"] = _bh_fdr(out["Dose_P"].to_numpy())
    out["Time_FDR"] = _bh_fdr(out["Time_P"].to_numpy())
    return out


# =====================================================================
# Section 6. Main runner
# =====================================================================


def run_task2_nodomain(cfg: Task2NoDomainConfig) -> None:
    _ensure_dir(cfg.analysis_dir)
    out = Path(cfg.analysis_dir)

    # Step 0: load
    pairwise = load_pairwise_wide(Path(cfg.pairwise_wide_path))
    retrieval = load_retrieval(Path(cfg.retrieval_per_query_path))

    # Step 1: build primary data + tests
    inst_l1, ctx_l1 = step1_build_l1(cfg=cfg, retrieval=retrieval, pairwise=pairwise)
    step1_tests = step1_tests_no_domain(ctx_l1)

    inst_l1.to_csv(out / "Step1_L1_Instance_Tidy_NoDomain.csv", index=False)
    ctx_l1.to_csv(out / "Step1_L1_Context_Aggregated_NoDomain.csv", index=False)
    step1_tests.to_csv(out / "Step1_Tests_NoDomain.csv", index=False)

    # Step 2: context labels + diagnostics
    labels, diag_sp = step2_context_labeling_no_domain(cfg=cfg, inst_l1=inst_l1, retrieval_all=retrieval)
    labels.to_csv(out / "Step2_Context_Labels_NoDomain.csv", index=False)
    if not diag_sp.empty:
        diag_sp.to_csv(out / "Step2_Diagnostics_SourcePairLevel_NoDomain.csv", index=False)

    # Step 4: case tracer
    tracer = step4_case_tracer(cfg=cfg, labels=labels, retrieval_all=retrieval)
    tracer.to_csv(out / "Step4_CaseStudy_Tracer_NoDomain.csv", index=False)

    # Step 5: enrichment + protocol
    enr_target_high = run_fisher_enrichment(
        labels=labels,
        group_col="Target",
        class_col="Performance_Class",
        positive_class="Robust_High",
        min_contexts=cfg.min_contexts_for_enrichment,
    )
    enr_target_low = run_fisher_enrichment(
        labels=labels,
        group_col="Target",
        class_col="Performance_Class",
        positive_class="Robust_Low",
        min_contexts=cfg.min_contexts_for_enrichment,
    )
    enr_target_ps = run_fisher_enrichment(
        labels=labels,
        group_col="Target",
        class_col="Performance_Class",
        positive_class="Protocol_Sensitive",
        min_contexts=cfg.min_contexts_for_enrichment,
    )
    enr_cell_high = run_fisher_enrichment(
        labels=labels,
        group_col="Cell",
        class_col="Performance_Class",
        positive_class="Robust_High",
        min_contexts=cfg.min_contexts_for_enrichment,
    )
    protocol = protocol_correlations(cfg=cfg, inst_l1=inst_l1, labels=labels, focus_class="Protocol_Sensitive")

    enr_target_high.to_csv(out / "Step5_Enrichment_Targets_RobustHigh_NoDomain.csv", index=False)
    enr_target_low.to_csv(out / "Step5_Enrichment_Targets_RobustLow_NoDomain.csv", index=False)
    enr_target_ps.to_csv(out / "Step5_Enrichment_Targets_ProtocolSensitive_NoDomain.csv", index=False)
    enr_cell_high.to_csv(out / "Step5_Enrichment_Cells_RobustHigh_NoDomain.csv", index=False)
    protocol.to_csv(out / "Step5_Protocol_Correlations_ProtocolSensitive_NoDomain.csv", index=False)

    run_manifest = {
        "config": asdict(cfg),
        "summary": {
            "n_pairwise_rows": int(len(pairwise)),
            "n_retrieval_rows": int(len(retrieval)),
            "n_step1_instance": int(len(inst_l1)),
            "n_step1_context": int(len(ctx_l1)),
            "n_step2_labels": int(len(labels)),
            "n_step4_tracer": int(len(tracer)),
        },
        "outputs": {"analysis_dir": str(out)},
    }
    Path(cfg.run_manifest_path).write_text(json.dumps(run_manifest, indent=2), encoding="utf-8")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("M2M-Bench Task2 mechanism analysis (No DomainType main contrast)")
    parser.add_argument("--pairwise-wide-path", type=str, default=Task2NoDomainConfig.pairwise_wide_path)
    parser.add_argument("--retrieval-per-query-path", type=str, default=Task2NoDomainConfig.retrieval_per_query_path)
    parser.add_argument("--output-dir", type=str, default=Task2NoDomainConfig.output_dir)
    parser.add_argument("--min-query-samples-per-context", type=int, default=Task2NoDomainConfig.min_query_samples_per_context)
    parser.add_argument("--min-contexts-for-enrichment", type=int, default=Task2NoDomainConfig.min_contexts_for_enrichment)
    parser.add_argument("--topk-cases-per-class", type=int, default=Task2NoDomainConfig.topk_cases_per_class)
    parser.add_argument("--hero-track", type=str, default=Task2NoDomainConfig.hero_track)
    parser.add_argument("--hero-view", type=str, default=Task2NoDomainConfig.hero_view)
    parser.add_argument(
        "--domain-scope",
        type=str,
        choices=["same_only", "all"],
        default=Task2NoDomainConfig.domain_scope,
        help="same_only: keep only SameDomain rows; all: keep both Same/Cross domain rows.",
    )
    parser.add_argument("--random-seed", type=int, default=Task2NoDomainConfig.random_seed)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    cfg = Task2NoDomainConfig(
        pairwise_wide_path=args.pairwise_wide_path,
        retrieval_per_query_path=args.retrieval_per_query_path,
        output_dir=args.output_dir,
        min_query_samples_per_context=args.min_query_samples_per_context,
        min_contexts_for_enrichment=args.min_contexts_for_enrichment,
        topk_cases_per_class=args.topk_cases_per_class,
        hero_track=args.hero_track,
        hero_view=args.hero_view,
        domain_scope=args.domain_scope,
        random_seed=args.random_seed,
    )
    run_task2_nodomain(cfg)


if __name__ == "__main__":
    main()
