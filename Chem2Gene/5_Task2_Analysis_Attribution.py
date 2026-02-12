#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Chem2Gen-Bench — Task 2 Unified Analysis (for Task1 MultiScenario v3 outputs)

Inputs:
  1) Task1_Pairwise_Metrics_Wide.csv
  2) Task1_Retrieval_MultiScenario_PerQuery.csv

Core principles:
  - Scenario A_L1Base is PRIMARY outcome for Task2.
  - "Domain" means source database (LINCS vs scPerturb).
  - For protocol analysis, keep instance-level (dose/time) on chemical side.
  - Keep BG correction consistency by using Task1 outputs (no re-building here).

[PATCH summary vs your previous version]
  - Step1: explicitly filter Direction == "Chem->Gene" to avoid silent mixing.
  - Tests:
      * paired Wilcoxon for Standard vs Systema and Gene vs Pathway
      * unpaired Mann–Whitney U for DomainType (SameDomain vs CrossDomain),
        using (Cell,Target) as samples by collapsing over source pairs.
  - Vectorize membership / key building (avoid apply/lambda on large tables).
  - Step2 labels: professional naming + explicit definitions.
  - Step2 diagnostics: include Source_Chem/Source_Gene/DomainType intermediate layer
    and export a source-pair-level diagnostics table.
  - bh_fdr: NaN p-values -> treated as 1.0
"""

import os
import numpy as np
import pandas as pd
from dataclasses import dataclass
from typing import Optional, Tuple, Dict, List

from tqdm import tqdm
from scipy.stats import wilcoxon, spearmanr, fisher_exact, mannwhitneyu


# -----------------------------
# Config
# -----------------------------
@dataclass
class Task2Config:
    BASE_DIR: str = "/mnt/NAS_21T/ProjectData/Chem2Gen"

    # Task1 outputs
    PAIRWISE_WIDE: str = "R_Vis_Ready/Task1_Metrics/Task1_Pairwise_Metrics_Wide.csv"
    RETRIEVAL_MULTI: str = "R_Vis_Ready/Task1_Retrieval_MultiScenario/Task1_Retrieval_MultiScenario_PerQuery.csv"

    # Optional: bundle path for UMAP (only if RUN_UMAP=True)
    BUNDLE_PT: str = "Benchmark_Datasets/Chem2Gen_Benchmark_MixedSources.pt"

    # Output folder for Task2
    OUT_DIR: str = "R_Vis_Ready/Task2_Unified"

    # Analysis knobs
    MIN_QUERY_SAMPLES_PER_CONTEXT: int = 3
    MIN_CONTEXTS_FOR_ENRICHMENT: int = 10

    # Hero slice for labeling (and default protocol correlations)
    HERO_TRACK: str = "Gene"
    HERO_VIEW: str = "Systema"
    HERO_DOMAIN: str = "SameDomain"

    # Case study
    TOPK_CASES: int = 20

    # Protocol correlations
    TIME_BINS: Dict[str, Tuple[float, float]] = None

    # Optional UMAP
    RUN_UMAP: bool = True
    UMAP_N_NEIGHBORS: int = 30
    UMAP_MIN_DIST: float = 0.3
    UMAP_RANDOM_STATE: int = 42
    UMAP_MAX_POINTS: int = 50000  # cap for speed/memory


def _default_time_bins():
    return {"Early": (0.0, 12.0), "Mid": (12.1, 36.0), "Late": (36.1, 1e9)}


# -----------------------------
# Stats helpers
# -----------------------------
def safe_wilcoxon(x: np.ndarray, y: np.ndarray) -> Tuple[float, float]:
    """Paired Wilcoxon signed-rank test robust to small n / all-ties."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    m = np.isfinite(x) & np.isfinite(y)
    x = x[m]
    y = y[m]
    if x.size < 2:
        return np.nan, np.nan
    diff = x - y
    if np.allclose(diff, 0):
        return 0.0, 1.0
    stat, p = wilcoxon(x, y)
    return float(stat), float(p)


def safe_mannwhitneyu(x: np.ndarray, y: np.ndarray, alternative: str = "two-sided") -> Tuple[float, float]:
    """Unpaired Mann–Whitney U test robust to small n."""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    x = x[np.isfinite(x)]
    y = y[np.isfinite(y)]
    if x.size < 2 or y.size < 2:
        return np.nan, np.nan
    stat, p = mannwhitneyu(x, y, alternative=alternative)
    return float(stat), float(p)


def rank_biserial_from_u(u_stat: float, n1: int, n2: int) -> float:
    """
    Rank-biserial correlation (equivalent to Cliff's delta up to sign convention),
    computed from Mann–Whitney U without O(N*M) memory.
    Range: [-1, 1].
    """
    if not np.isfinite(u_stat) or n1 <= 0 or n2 <= 0:
        return np.nan
    return float((2.0 * u_stat) / (n1 * n2) - 1.0)


def bh_fdr(p_values: np.ndarray) -> np.ndarray:
    """Benjamini–Hochberg FDR correction. NaN treated as 1.0."""
    p = np.asarray(p_values, dtype=float)
    if p.size == 0:
        return p
    p = np.where(np.isfinite(p), p, 1.0)  # [PATCH] NaN -> 1.0
    n = p.size
    order = np.argsort(p)
    ranked = p[order]
    corr = ranked * (n / np.arange(1, n + 1))
    corr = np.minimum.accumulate(corr[::-1])[::-1]
    out = np.empty_like(corr)
    out[order] = np.clip(corr, 0, 1)
    return out


# -----------------------------
# Formatting / feature helpers
# -----------------------------
def ensure_dir(path: str):
    os.makedirs(path, exist_ok=True)


def parse_log10_dose_vec(dose: np.ndarray) -> np.ndarray:
    dose = np.asarray(dose, dtype=float)
    out = np.full(dose.shape, np.nan, dtype=float)
    m = np.isfinite(dose) & (dose > 0)
    out[m] = np.log10(dose[m])
    return out


def bin_time_vec(time: np.ndarray, time_bins: Dict[str, Tuple[float, float]]) -> np.ndarray:
    t = np.asarray(time, dtype=float)
    out = np.full(t.shape, "Unknown", dtype=object)
    m = np.isfinite(t)
    if not m.any():
        return out
    for label, (lo, hi) in time_bins.items():
        sel = m & (t >= float(lo)) & (t <= float(hi))
        out[sel] = label
    # Anything finite but not assigned
    out[m & (out == "Unknown")] = "Other"
    return out


def format_num_key_vec(x: np.ndarray, ndigits: int = 6) -> np.ndarray:
    """
    Vectorized stable numeric key:
      finite -> formatted '%.{ndigits}f'
      non-finite -> 'nan'
    """
    arr = np.asarray(x, dtype=float)
    out = np.full(arr.shape, "nan", dtype=object)
    m = np.isfinite(arr)
    if m.any():
        out[m] = np.char.mod(f"%.{ndigits}f", np.round(arr[m], ndigits))
    return out


# -----------------------------
# IO + standardization
# -----------------------------
def load_pairwise_wide(cfg: Task2Config) -> pd.DataFrame:
    path = os.path.join(cfg.BASE_DIR, cfg.PAIRWISE_WIDE)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing pairwise file: {path}")
    df = pd.read_csv(path)

    required = [
        "Track", "cell_std", "target_std", "domain_type", "source_db_chem", "source_db_gene",
        "cosine_std", "cosine_sys", "systema_gain_cosine"
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Pairwise wide missing columns: {missing}\nFound: {list(df.columns)}")

    df = df.rename(columns={
        "cell_std": "Cell",
        "target_std": "Target",
        "domain_type": "DomainType",
        "source_db_chem": "Source_Chem",
        "source_db_gene": "Source_Gene",
    })

    keep = ["Track", "Cell", "Target", "DomainType", "Source_Chem", "Source_Gene",
            "cosine_std", "cosine_sys", "systema_gain_cosine"]
    return df[keep].copy()


def load_retrieval(cfg: Task2Config) -> pd.DataFrame:
    path = os.path.join(cfg.BASE_DIR, cfg.RETRIEVAL_MULTI)
    if not os.path.exists(path):
        raise FileNotFoundError(f"Missing retrieval file: {path}")
    df = pd.read_csv(path)

    required = [
        "Track", "Scenario", "View", "Direction", "LabelType", "Cell", "Target",
        "Source_Chem", "Source_Gene", "Dose", "Time", "CondID",
        "True_Rank", "MRR", "Success_Score", "N_Gallery"
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"Retrieval CSV missing columns: {missing}\nFound: {list(df.columns)}")

    df["DomainType"] = np.where(df["Source_Chem"] == df["Source_Gene"], "SameDomain", "CrossDomain")

    df["Dose_num"] = pd.to_numeric(df["Dose"], errors="coerce")
    df["Time_num"] = pd.to_numeric(df["Time"], errors="coerce")
    df["True_Rank"] = pd.to_numeric(df["True_Rank"], errors="coerce")
    df["N_Gallery"] = pd.to_numeric(df["N_Gallery"], errors="coerce")

    # Ensure Success_Score exists
    if "Success_Score" not in df.columns or df["Success_Score"].isna().all():
        n = df["N_Gallery"].fillna(1).clip(lower=1)
        r = df["True_Rank"].fillna(n)
        denom = (n - 1.0).replace(0, 1.0)
        df["Success_Score"] = 1.0 - (r - 1.0) / denom

    # Protocol features (vectorized)
    if cfg.TIME_BINS is None:
        cfg.TIME_BINS = _default_time_bins()

    df["LogDose"] = parse_log10_dose_vec(df["Dose_num"].to_numpy())
    df["TimeBin"] = bin_time_vec(df["Time_num"].to_numpy(), cfg.TIME_BINS)

    # Dedup-safe key (vectorized keys for Dose/Time)
    df["_DoseKey"] = format_num_key_vec(df["Dose_num"].to_numpy(), 6)
    df["_TimeKey"] = format_num_key_vec(df["Time_num"].to_numpy(), 6)

    key_cols = [
        "Track", "Scenario", "View", "Direction", "LabelType",
        "Cell", "Target", "Source_Chem", "Source_Gene", "CondID", "_DoseKey", "_TimeKey"
    ]
    df["_key"] = df[key_cols].astype(str).agg("|".join, axis=1)

    # keep best rank if duplicates exist
    df = df.sort_values(["True_Rank"], ascending=True).drop_duplicates("_key").reset_index(drop=True)

    return df


# -----------------------------
# Step 1: L1 tidy + attach pairwise
# -----------------------------
def attach_pairwise_context(df_l1: pd.DataFrame, df_pair: pd.DataFrame) -> pd.DataFrame:
    """
    Attach pairwise cosine to each L1 retrieval instance by:
      (Track, Cell, Target, DomainType, Source_Chem, Source_Gene, View)

    Pairwise wide has cosine_std/cosine_sys. We map by View:
      Standard -> cosine_std
      Systema  -> cosine_sys
    """
    gcols = ["Track", "Cell", "Target", "DomainType", "Source_Chem", "Source_Gene"]
    dfp = df_pair.groupby(gcols, as_index=False).agg({
        "cosine_std": "mean",
        "cosine_sys": "mean",
        "systema_gain_cosine": "mean"
    })

    std = dfp.copy()
    std["View"] = "Standard"
    std["Pairwise_Cosine"] = std["cosine_std"]
    std["Pairwise_Systema_Gain"] = 0.0  # kept as discussed (avoid changing downstream)

    sys = dfp.copy()
    sys["View"] = "Systema"
    sys["Pairwise_Cosine"] = sys["cosine_sys"]
    sys["Pairwise_Systema_Gain"] = sys["systema_gain_cosine"]

    ctx = pd.concat([std, sys], ignore_index=True)
    ctx = ctx[gcols + ["View", "Pairwise_Cosine", "Pairwise_Systema_Gain"]]

    out = df_l1.merge(ctx, on=gcols + ["View"], how="left", validate="m:1")
    return out


def step1_build_l1_tidy(cfg: Task2Config, df_retr: pd.DataFrame, df_pair: pd.DataFrame) -> Tuple[pd.DataFrame, pd.DataFrame]:
    """
    Primary outcome: Scenario A_L1Base / Target retrieval / Chem->Gene.
    Produces:
      - Instance-level tidy
      - Context-aggregated tidy at (Track,View,DomainType,Cell,Target,Source_Chem,Source_Gene)
    """
    # [PATCH] Explicit Chem->Gene filter to prevent silent mixing
    l1 = df_retr[
        (df_retr["Scenario"] == "A_L1Base") &
        (df_retr["LabelType"] == "Target") &
        (df_retr["Direction"] == "Chem->Gene")
    ].copy()

    if l1.empty:
        raise ValueError("No Scenario==A_L1Base & LabelType==Target & Direction==Chem->Gene rows found.")

    l1 = attach_pairwise_context(l1, df_pair)

    inst_cols = [
        "Track", "Scenario", "View", "Direction", "DomainType",
        "Cell", "Target", "Source_Chem", "Source_Gene",
        "Dose", "Time", "CondID", "LogDose", "TimeBin",
        "True_Rank", "MRR", "Success_Score", "N_Gallery",
        "Pairwise_Cosine", "Pairwise_Systema_Gain"
    ]
    inst = l1[inst_cols].copy()

    ctx_gcols = ["Track", "View", "DomainType", "Cell", "Target", "Source_Chem", "Source_Gene"]
    ctx = inst.groupby(ctx_gcols, as_index=False).agg(
        N_Instances=("Success_Score", "size"),
        Mean_Success=("Success_Score", "mean"),
        Median_Success=("Success_Score", "median"),
        Max_Success=("Success_Score", "max"),
        Mean_MRR=("MRR", "mean"),
        Pairwise_Cosine=("Pairwise_Cosine", "mean"),
    )

    # enforce minimum instances per context
    ctx = ctx[ctx["N_Instances"] >= cfg.MIN_QUERY_SAMPLES_PER_CONTEXT].reset_index(drop=True)
    return inst, ctx


# -----------------------------
# Step 1b: tests (paired & unpaired)
# -----------------------------
def paired_tests(ctx: pd.DataFrame) -> pd.DataFrame:
    """
    Tests on context-aggregated Mean_Success.

    - Standard vs Systema: paired Wilcoxon on same context
      index = (DomainType, Cell, Target, Source_Chem, Source_Gene)
    - Gene vs Pathway: paired Wilcoxon on same context
      index = (DomainType, Cell, Target, Source_Chem, Source_Gene)
    - SameDomain vs CrossDomain: unpaired Mann–Whitney U (U-test),
      using (Cell,Target) as samples by first collapsing over source-pairs.
    """
    rows = []

    # 1) Standard vs Systema (paired)
    for track in ["Gene", "Pathway"]:
        sub = ctx[ctx["Track"] == track].copy()
        piv = sub.pivot_table(
            index=["DomainType", "Cell", "Target", "Source_Chem", "Source_Gene"],
            columns="View",
            values="Mean_Success",
            aggfunc="mean"
        ).dropna()

        if {"Standard", "Systema"}.issubset(piv.columns) and len(piv) >= 2:
            stat, p = safe_wilcoxon(piv["Systema"].values, piv["Standard"].values)
            diff = piv["Systema"].values - piv["Standard"].values
            rows.append({
                "Comparison": f"Systema vs Standard (Mean_Success) [{track}]",
                "Test": "Paired_Wilcoxon",
                "N": int(len(piv)),
                "P_Value": p,
                "Median_Diff": float(np.nanmedian(diff)),
                "Mean_Diff": float(np.nanmean(diff)),
                "Frac_Positive": float(np.mean(diff > 0)),
            })

    # 2) Gene vs Pathway (paired) within view
    for view in ["Standard", "Systema"]:
        sub = ctx[ctx["View"] == view].copy()
        piv = sub.pivot_table(
            index=["DomainType", "Cell", "Target", "Source_Chem", "Source_Gene"],
            columns="Track",
            values="Mean_Success",
            aggfunc="mean"
        ).dropna()

        if {"Gene", "Pathway"}.issubset(piv.columns) and len(piv) >= 2:
            stat, p = safe_wilcoxon(piv["Pathway"].values, piv["Gene"].values)
            diff = piv["Pathway"].values - piv["Gene"].values
            rows.append({
                "Comparison": f"Pathway vs Gene (Mean_Success) [{view}]",
                "Test": "Paired_Wilcoxon",
                "N": int(len(piv)),
                "P_Value": p,
                "Median_Diff": float(np.nanmedian(diff)),
                "Mean_Diff": float(np.nanmean(diff)),
                "Frac_Positive": float(np.mean(diff > 0)),
            })

    # 3) SameDomain vs CrossDomain (unpaired U-test), using (Cell,Target) as samples
    # [PATCH] collapse across source-pairs first, so each (Cell,Target) contributes one value per DomainType
    collapsed = ctx.groupby(["Track", "View", "Cell", "Target", "DomainType"], as_index=False).agg(
        Mean_Success=("Mean_Success", "mean"),
        N_SourcePairs=("Mean_Success", "size"),
    )

    for track in ["Gene", "Pathway"]:
        for view in ["Standard", "Systema"]:
            sub = collapsed[(collapsed["Track"] == track) & (collapsed["View"] == view)]
            same = sub[sub["DomainType"] == "SameDomain"]["Mean_Success"].to_numpy()
            cross = sub[sub["DomainType"] == "CrossDomain"]["Mean_Success"].to_numpy()

            if len(same) >= 2 and len(cross) >= 2:
                u, p = safe_mannwhitneyu(same, cross, alternative="two-sided")
                rbc = rank_biserial_from_u(u, n1=len(same), n2=len(cross))
                rows.append({
                    "Comparison": f"SameDomain vs CrossDomain (Mean_Success) [{track}/{view}]",
                    "Test": "MannWhitneyU",
                    "N_Same": int(len(same)),
                    "N_Cross": int(len(cross)),
                    "P_Value": p,
                    "Median_Diff": float(np.nanmedian(same) - np.nanmedian(cross)),
                    "Mean_Diff": float(np.nanmean(same) - np.nanmean(cross)),
                    "RankBiserial": rbc,
                })

    out = pd.DataFrame(rows)
    if not out.empty:
        out["FDR_BH"] = bh_fdr(out["P_Value"].values)
    return out


# -----------------------------
# Step 2: Context labeling + diagnostics
# -----------------------------
def step2_context_labeling(cfg: Task2Config, inst_l1: pd.DataFrame, ctx_l1: pd.DataFrame, df_retr_all: pd.DataFrame):
    """
    Professional context labeling on (Cell,Target) from HERO slice.

    HERO slice defaults: Track=cfg.HERO_TRACK, View=cfg.HERO_VIEW, DomainType=cfg.HERO_DOMAIN, Direction=Chem->Gene.

    We compute per (Cell,Target):
      - Mean_Success, Peak_Success (=max), N_Instances, Context_Cosine

    Tiers:
      - Mean_Tier: tertiles of Mean_Success
      - Peak_Tier: tertiles of Peak_Success

    Performance_Class (explicit definitions):
      - Robust_High:       Peak_Tier=High AND Mean_Tier=High
      - Protocol_Sensitive:Peak_Tier=High AND Mean_Tier!=High
      - Robust_Low:        Peak_Tier=Low
      - Intermediate:      otherwise
    """
    hero = inst_l1[
        (inst_l1["Track"] == cfg.HERO_TRACK) &
        (inst_l1["View"] == cfg.HERO_VIEW) &
        (inst_l1["Direction"] == "Chem->Gene") &
        (inst_l1["DomainType"] == cfg.HERO_DOMAIN)
    ].copy()

    if hero.empty:
        hero = inst_l1[
            (inst_l1["Track"] == cfg.HERO_TRACK) &
            (inst_l1["View"] == cfg.HERO_VIEW) &
            (inst_l1["Direction"] == "Chem->Gene")
        ].copy()

    if hero.empty:
        hero = inst_l1[(inst_l1["View"] == cfg.HERO_VIEW) & (inst_l1["Direction"] == "Chem->Gene")].copy()

    if hero.empty:
        raise ValueError("Hero slice is empty; cannot build context labels.")

    sc = hero.groupby(["Cell", "Target"], as_index=False).agg(
        Mean_Success=("Success_Score", "mean"),
        Peak_Success=("Success_Score", "max"),
        N_Instances=("Success_Score", "size"),
        Context_Cosine=("Pairwise_Cosine", "mean"),
    )

    if len(sc) < 3:
        sc["Mean_Tier"] = "Mid"
        sc["Peak_Tier"] = "Mid"
    else:
        sc["Mean_Tier"] = pd.qcut(sc["Mean_Success"].rank(method="first"), 3, labels=["Low", "Mid", "High"])
        sc["Peak_Tier"] = pd.qcut(sc["Peak_Success"].rank(method="first"), 3, labels=["Low", "Mid", "High"])

    def classify(mean_t, peak_t):
        if str(peak_t) == "Low":
            return "Robust_Low"
        if str(peak_t) == "High" and str(mean_t) == "High":
            return "Robust_High"
        if str(peak_t) == "High" and str(mean_t) != "High":
            return "Protocol_Sensitive"
        return "Intermediate"

    sc["Performance_Class"] = [classify(m, p) for m, p in zip(sc["Mean_Tier"], sc["Peak_Tier"])]

    # ---- Diagnostics from L2/L3 (same retrieval file) ----
    # [PATCH] Keep Source_Chem/Source_Gene/DomainType in an intermediate layer, then collapse.
    diag_track = cfg.HERO_TRACK

    def _diag_sourcepair_level(df: pd.DataFrame, scenario: str, prefix: str) -> pd.DataFrame:
        sub = df[(df["Scenario"] == scenario) & (df["Track"] == diag_track) & (df["Direction"] == "Chem->Gene")].copy()
        if sub.empty:
            return pd.DataFrame()

        gcols = ["Cell", "Target", "DomainType", "Source_Chem", "Source_Gene", "View"]
        agg = sub.groupby(gcols, as_index=False).agg(
            MeanSuccess=("Success_Score", "mean"),
            MeanMRR=("MRR", "mean"),
            N=("Success_Score", "size"),
        )
        agg = agg.rename(columns={
            "MeanSuccess": f"{prefix}_MeanSuccess",
            "MeanMRR": f"{prefix}_MeanMRR",
            "N": f"{prefix}_N",
        })
        return agg

    l2_sp = _diag_sourcepair_level(df_retr_all, "B_CellDep", "L2_CellDep")
    l3_sp = _diag_sourcepair_level(df_retr_all, "C_TgtSpec", "L3_TgtSpec")

    # collapse source-pair-level -> (Cell,Target,DomainType,View)
    def _collapse_sp(sp: pd.DataFrame, prefix: str) -> pd.DataFrame:
        if sp.empty:
            return pd.DataFrame()
        gcols = ["Cell", "Target", "DomainType", "View"]
        collapsed = sp.groupby(gcols, as_index=False).agg(
            **{
                f"{prefix}_MeanSuccess": (f"{prefix}_MeanSuccess", "mean"),
                f"{prefix}_MeanMRR": (f"{prefix}_MeanMRR", "mean"),
                f"{prefix}_N": (f"{prefix}_N", "sum"),
                f"{prefix}_N_SourcePairs": ("Source_Chem", "size"),
            }
        )
        return collapsed

    l2_col = _collapse_sp(l2_sp, "L2_CellDep")
    l3_col = _collapse_sp(l3_sp, "L3_TgtSpec")

    def _pivot_diag(collapsed: pd.DataFrame, prefix: str) -> pd.DataFrame:
        if collapsed.empty:
            return pd.DataFrame(columns=["Cell", "Target"])
        # pivot to wide: columns like {prefix}_MeanSuccess__{DomainType}__{View}
        idx = ["Cell", "Target"]
        val_cols = [f"{prefix}_MeanSuccess", f"{prefix}_MeanMRR", f"{prefix}_N", f"{prefix}_N_SourcePairs"]
        piv = collapsed.pivot_table(index=idx, columns=["DomainType", "View"], values=val_cols, aggfunc="first")
        piv.columns = [f"{a}__{d}__{v}" for (a, d, v) in piv.columns]
        piv = piv.reset_index()
        return piv

    l2_piv = _pivot_diag(l2_col, "L2_CellDep")
    l3_piv = _pivot_diag(l3_col, "L3_TgtSpec")

    out = sc.merge(l2_piv, on=["Cell", "Target"], how="left").merge(l3_piv, on=["Cell", "Target"], how="left")

    # also return source-pair-level diagnostics for separate export
    diag_sp = None
    if (not l2_sp.empty) or (not l3_sp.empty):
        l2_sp2 = l2_sp.copy()
        l2_sp2["Diag"] = "L2_CellDep"
        l3_sp2 = l3_sp.copy()
        l3_sp2["Diag"] = "L3_TgtSpec"
        diag_sp = pd.concat([l2_sp2, l3_sp2], ignore_index=True)

    return out, diag_sp


# -----------------------------
# Step 4: Case study tracer across A/B/C
# -----------------------------
# def step4_case_study_tracer(cfg: Task2Config, labels: pd.DataFrame, df_retr_all: pd.DataFrame) -> pd.DataFrame:
#     """
#     Select top contexts from Robust_High (by Mean_Success), then trace across A/B/C.
#     """
#     good = labels[labels["Performance_Class"] == "Robust_High"].copy()
#     if good.empty:
#         return pd.DataFrame()

#     good = good.sort_values("Mean_Success", ascending=False).head(cfg.TOPK_CASES)
#     good["_ct_key"] = good["Cell"].astype(str) + "||" + good["Target"].astype(str)
#     keep = set(good["_ct_key"].values.tolist())

#     df = df_retr_all.copy()
#     df["_ct_key"] = df["Cell"].astype(str) + "||" + df["Target"].astype(str)
#     df = df[df["_ct_key"].isin(keep)].copy()

#     tracer = df.groupby(["Cell", "Target", "Scenario", "Track", "View"], as_index=False).agg(
#         MeanSuccess=("Success_Score", "mean"),
#         MeanMRR=("MRR", "mean"),
#         N=("Success_Score", "size"),
#     )
#     return tracer
def step4_case_study_tracer(cfg: Task2Config, labels: pd.DataFrame, df_retr_all: pd.DataFrame) -> pd.DataFrame:
    """
    Step 4: Case-study tracer across A/B/C.

    Requirements (per latest request):
      1) Keep case selection strategy consistent with original:
         - select top-K contexts *within each Performance_Class* using class-consistent ranking
         - Robust_High selection behavior remains identical to original (Mean_Success desc)
      2) Add Performance_Class column in tracer output (acceptable).
      3) Do NOT add extra hard filters (no MIN_QUERY filter; no extra Direction/Domain hard filtering).
      4) Fix Level2 disappearance:
         - prefer LabelType=='Target' PER SCENARIO if it exists; otherwise keep all rows for that scenario.

    Note:
      - K is cfg.TOPK_CASES per class. If you want exactly 10 per class to match R side,
        set TOPK_CASES=10 in config.
    """
    if labels is None or labels.empty:
        return pd.DataFrame()

    if "Performance_Class" not in labels.columns:
        raise ValueError("labels table missing required column: Performance_Class")

    k = int(cfg.TOPK_CASES) if cfg.TOPK_CASES is not None else 20
    k = max(1, k)

    lab = labels.copy()
    lab["Cell"] = lab["Cell"].astype(str)
    lab["Target"] = lab["Target"].astype(str)

    # Build a sensitivity score if possible (for Protocol_Sensitive ranking)
    if ("Peak_Success" in lab.columns) and ("Mean_Success" in lab.columns):
        lab["SensitivityScore"] = pd.to_numeric(lab["Peak_Success"], errors="coerce") - pd.to_numeric(lab["Mean_Success"], errors="coerce")
    else:
        lab["SensitivityScore"] = np.nan

    # ---- Case picking: consistent per-class TopK ----
    picks = []

    # (A) Robust_High: keep original behavior (Mean_Success desc)
    rh = lab[lab["Performance_Class"] == "Robust_High"].copy()
    if not rh.empty and "Mean_Success" in rh.columns:
        rh = rh.sort_values("Mean_Success", ascending=False).head(k)
        picks.append(rh[["Cell", "Target", "Performance_Class"]])

    # (B) Protocol_Sensitive: rank by SensitivityScore desc (fallback to Peak_Success desc if needed)
    ps = lab[lab["Performance_Class"] == "Protocol_Sensitive"].copy()
    if not ps.empty:
        if ps["SensitivityScore"].notna().any():
            ps = ps.sort_values("SensitivityScore", ascending=False).head(k)
        elif "Peak_Success" in ps.columns:
            ps = ps.sort_values("Peak_Success", ascending=False).head(k)
        elif "Mean_Success" in ps.columns:
            ps = ps.sort_values("Mean_Success", ascending=False).head(k)
        picks.append(ps[["Cell", "Target", "Performance_Class"]])

    # (C) Robust_Low: rank by Mean_Success asc (low examples)
    rl = lab[lab["Performance_Class"] == "Robust_Low"].copy()
    if not rl.empty and "Mean_Success" in rl.columns:
        rl = rl.sort_values("Mean_Success", ascending=True).head(k)
        picks.append(rl[["Cell", "Target", "Performance_Class"]])

    if not picks:
        return pd.DataFrame()

    pick_ctx = (
        pd.concat(picks, ignore_index=True)
        .drop_duplicates(["Cell", "Target"], keep="first")
        .reset_index(drop=True)
    )

    # ---- Filter retrieval to selected contexts only (no extra hard filters beyond this) ----
    df = df_retr_all.copy()
    df["Cell"] = df["Cell"].astype(str)
    df["Target"] = df["Target"].astype(str)

    # keep only A/B/C scenarios (this matches the intended tracer figure and typical retrieval schema)
    if "Scenario" in df.columns:
        df = df[df["Scenario"].isin(["A_L1Base", "B_CellDep", "C_TgtSpec"])].copy()

    df = df.merge(pick_ctx[["Cell", "Target"]], on=["Cell", "Target"], how="inner")
    if df.empty:
        return pd.DataFrame()

    # ---- Fix "Level2 disappears": LabelType preference per scenario (fallback if none) ----
    if "LabelType" in df.columns:
        parts = []
        for scen in ["A_L1Base", "B_CellDep", "C_TgtSpec"]:
            dsc = df[df["Scenario"] == scen].copy()
            if dsc.empty:
                continue
            dsc_t = dsc[dsc["LabelType"] == "Target"].copy()
            parts.append(dsc_t if not dsc_t.empty else dsc)
        if parts:
            df = pd.concat(parts, ignore_index=True)

    # attach Performance_Class (allowed by your requirement #2)
    df = df.merge(
        labels[["Cell", "Target", "Performance_Class"]].drop_duplicates(["Cell", "Target"]),
        on=["Cell", "Target"],
        how="left"
    )

    # ---- Aggregate tracer (now includes Performance_Class) ----
    tracer = df.groupby(
        ["Cell", "Target", "Performance_Class", "Scenario", "Track", "View"],
        as_index=False
    ).agg(
        MeanSuccess=("Success_Score", "mean"),
        MeanMRR=("MRR", "mean"),
        N=("Success_Score", "size"),
    )

    return tracer

# -----------------------------
# Step 5: Enrichment + protocol correlations
# -----------------------------
def calculate_stable_log2or(a, b, c, d) -> float:
    """Haldane–Anscombe correction for stable log2 odds ratio."""
    num = (a + 0.5) * (d + 0.5)
    den = (b + 0.5) * (c + 0.5)
    return float(np.log2(num / den))


def run_fisher_enrichment(df: pd.DataFrame,
                          group_col: str,
                          class_col: str,
                          positive_class: str,
                          min_contexts: int,
                          alternative: str = "greater") -> pd.DataFrame:
    """
    One-sided Fisher enrichment for membership in positive_class.
    Operates on labels table (one row per Cell,Target).
    """
    if df.empty:
        return pd.DataFrame()

    total_pos = int((df[class_col] == positive_class).sum())
    total_neg = int((df[class_col] != positive_class).sum())
    total_all = total_pos + total_neg
    if total_all == 0 or total_pos == 0 or total_neg == 0:
        return pd.DataFrame()

    baseline = total_pos / total_all
    results = []

    uniq = df[group_col].dropna().astype(str).unique()
    for g in tqdm(uniq, desc=f"Fisher {group_col} -> {positive_class}"):
        is_g = (df[group_col].astype(str) == str(g))
        n_ctx = int(is_g.sum())
        if n_ctx < min_contexts:
            continue

        a = int((is_g & (df[class_col] == positive_class)).sum())
        c = int(n_ctx - a)
        b = int(total_pos - a)
        d = int(total_neg - c)

        table = [[a, b], [c, d]]
        or_raw, p = fisher_exact(table, alternative=alternative)

        results.append({
            "Feature": str(g),
            "GroupCol": group_col,
            "Target_Class": positive_class,
            "Alternative": alternative,
            "N_Contexts": n_ctx,
            "N_Positive": a,
            "Baseline_Rate": baseline,
            "Rate_Group": a / n_ctx if n_ctx > 0 else np.nan,
            "Rate_Background": b / (b + d) if (b + d) > 0 else np.nan,
            "Rate_Diff": (a / n_ctx) - (b / (b + d)) if (n_ctx > 0 and (b + d) > 0) else np.nan,
            "Odds_Ratio_Raw": float(or_raw) if np.isfinite(or_raw) else np.nan,
            "Log2_OR_Stable": calculate_stable_log2or(a, b, c, d),
            "P_Value": float(p),
        })

    res = pd.DataFrame(results)
    if res.empty:
        return res

    res["FDR_BH"] = bh_fdr(res["P_Value"].values)
    res["NegLog10_P"] = -np.log10(np.where(np.isfinite(res["P_Value"].values), res["P_Value"].values, 1.0) + 1e-300)
    return res.sort_values(["FDR_BH", "P_Value", "Log2_OR_Stable"],
                           ascending=[True, True, False],
                           kind="mergesort").reset_index(drop=True)


def protocol_correlations_by_class(cfg: Task2Config, inst_l1: pd.DataFrame, labels: pd.DataFrame,
                                   class_col: str = "Performance_Class",
                                   focus_class: str = "Protocol_Sensitive") -> pd.DataFrame:
    """
    For a selected class (default: Protocol_Sensitive), compute within-(Cell,Target,SourcePair)
    Spearman correlations between Success_Score and LogDose / Time.

    [PATCH] Default restriction to HERO slice to avoid mixing View/Domain/Direction.
    """
    lb = labels[["Cell", "Target", class_col]].copy()
    sub = inst_l1.merge(lb, on=["Cell", "Target"], how="inner")

    # focus class
    sub = sub[sub[class_col] == focus_class].copy()
    if sub.empty:
        return pd.DataFrame()

    # [PATCH] Restrict to hero slice by default to avoid mixing interpretation
    sub = sub[
        (sub["Track"] == cfg.HERO_TRACK) &
        (sub["View"] == cfg.HERO_VIEW) &
        (sub["DomainType"] == cfg.HERO_DOMAIN) &
        (sub["Direction"] == "Chem->Gene")
    ].copy()

    if sub.empty:
        return pd.DataFrame()

    rows = []
    group_cols = ["Cell", "Target", "Source_Chem", "Source_Gene"]
    for (cell, tgt, sc, sg), g in tqdm(sub.groupby(group_cols), desc=f"Protocol corr [{focus_class}]"):
        if len(g) < 6:
            continue

        # Dose correlation: within most common TimeBin (prefer Mid if tied)
        dose_corr, dose_p = np.nan, np.nan
        tb_counts = g["TimeBin"].value_counts()
        if len(tb_counts) > 0:
            top_tb = "Mid" if ("Mid" in tb_counts and tb_counts["Mid"] == tb_counts.max()) else tb_counts.idxmax()
            gd = g[g["TimeBin"] == top_tb]
            if gd["LogDose"].nunique(dropna=True) >= 3:
                dose_corr, dose_p = spearmanr(gd["LogDose"], gd["Success_Score"], nan_policy="omit")

        # Time correlation: within top-half doses
        time_corr, time_p = np.nan, np.nan
        if g["LogDose"].notna().sum() >= 6:
            med = np.nanmedian(g["LogDose"].values)
            gt = g[g["LogDose"] >= med]
            if len(gt) >= 6 and pd.to_numeric(gt["Time"], errors="coerce").nunique(dropna=True) >= 3:
                time_corr, time_p = spearmanr(pd.to_numeric(gt["Time"], errors="coerce"),
                                              gt["Success_Score"], nan_policy="omit")

        rows.append({
            "Cell": cell,
            "Target": tgt,
            "Source_Chem": sc,
            "Source_Gene": sg,
            "N": int(len(g)),
            "Dose_Corr": float(dose_corr) if np.isfinite(dose_corr) else np.nan,
            "Dose_P": float(dose_p) if np.isfinite(dose_p) else np.nan,
            "Time_Corr": float(time_corr) if np.isfinite(time_corr) else np.nan,
            "Time_P": float(time_p) if np.isfinite(time_p) else np.nan,
        })

    out = pd.DataFrame(rows)
    if out.empty:
        return out
    out["Dose_FDR"] = bh_fdr(out["Dose_P"].values)
    out["Time_FDR"] = bh_fdr(out["Time_P"].values)
    return out


# -----------------------------
# Optional Step 3: UMAP (kept as-is; skips if bundle missing)
# -----------------------------
def step3_dual_umap(cfg: Task2Config, labels: pd.DataFrame) -> Optional[pd.DataFrame]:
    try:
        import torch
        import umap  # type: ignore
    except Exception as e:
        print(f"[WARN] Step3 UMAP skipped (missing deps): {e}")
        return None

    bundle_path = os.path.join(cfg.BASE_DIR, cfg.BUNDLE_PT)
    if not os.path.exists(bundle_path):
        print(f"[WARN] Step3 UMAP skipped (missing bundle): {bundle_path}")
        return None

    print(">>> [Step 3] Running UMAP on centroids: (Cell,Target,PerturbType) (Standard)...")
    bundle = torch.load(bundle_path, map_location="cpu")

    if "Level1" not in bundle or "pairs_meta" not in bundle["Level1"]:
        print("[WARN] Step3 UMAP skipped: bundle missing Level1/pairs_meta.")
        return None

    lvl1 = bundle["Level1"]
    meta = lvl1["pairs_meta"].copy().reset_index(drop=True)

    if "chem_tensors" not in lvl1 or "gene_tensors" not in lvl1:
        print("[WARN] Step3 UMAP skipped: bundle missing Level1/chem_tensors or Level1/gene_tensors.")
        return None

    chem_tensors = lvl1["chem_tensors"]
    gene_tensors = lvl1["gene_tensors"]

    req = ["uid_chem", "uid_gene", "cell_std", "target_std"]
    miss = [c for c in req if c not in meta.columns]
    if miss:
        print(f"[WARN] Step3 UMAP skipped: meta missing required columns: {miss}")
        return None

    m = meta.copy()
    m["Cell"] = m["cell_std"].astype(str)
    m["Target"] = m["target_std"].astype(str)

    lb = labels[["Cell", "Target", "Performance_Class"]].copy()
    lb["Cell"] = lb["Cell"].astype(str)
    lb["Target"] = lb["Target"].astype(str)

    def centroids_chunked(meta_df, X_all, group_cols, dedup_col, perturb_type, track, chunk_size=50000):
        meta_u = meta_df.drop_duplicates(dedup_col, keep="first").copy()
        meta_u = meta_u.reset_index(drop=False)

        if meta_u.empty:
            return pd.DataFrame(), np.zeros((0, 1), dtype=np.float32)

        key_str = meta_u[group_cols].astype(str).agg("|".join, axis=1).values
        codes, uniq = pd.factorize(key_str, sort=True)
        K = len(uniq)
        D = int(X_all.shape[1])

        sums = np.zeros((K, D), dtype=np.float64)
        cnts = np.zeros(K, dtype=np.int64)

        orig_idx = meta_u["index"].values.astype(int)
        n = len(orig_idx)

        for s in range(0, n, chunk_size):
            e = min(n, s + chunk_size)
            idx_chunk = orig_idx[s:e]
            code_chunk = codes[s:e]
            X_chunk = X_all[torch.as_tensor(idx_chunk, dtype=torch.long)].numpy().astype(np.float32)
            np.add.at(sums, code_chunk, X_chunk)
            cnts += np.bincount(code_chunk, minlength=K).astype(np.int64)

        cents = (sums / np.maximum(cnts, 1)[:, None]).astype(np.float32)

        meta_u["_code"] = codes
        rep = meta_u.groupby("_code", as_index=False).first()

        out_keys = rep[group_cols].copy()
        out_keys["N_Samples"] = cnts.astype(int)
        out_keys["PerturbType"] = perturb_type
        out_keys["Track"] = track
        out_keys["View"] = "Standard"
        return out_keys.reset_index(drop=True), cents

    tracks = [("Gene", "y_gene"), ("Pathway", "y_path")]
    all_rows = []

    for track_name, tensor_key in tracks:
        if tensor_key not in chem_tensors or tensor_key not in gene_tensors:
            print(f"[WARN] Step3 UMAP: missing {tensor_key} for {track_name}, skipping.")
            continue

        Xc_all = chem_tensors[tensor_key]
        Xg_all = gene_tensors[tensor_key]

        if len(Xc_all) != len(m) or len(Xg_all) != len(m):
            raise ValueError(f"CRITICAL: Tensor/meta misalignment for {track_name}/{tensor_key}.")

        group_cols = ["Cell", "Target"]

        chem_df, chem_cent = centroids_chunked(m, Xc_all, group_cols, "uid_chem", "Chemical", track_name)
        gene_df, gene_cent = centroids_chunked(m, Xg_all, group_cols, "uid_gene", "Genetic", track_name)

        if len(chem_df) == 0 or len(gene_df) == 0:
            continue

        chem_df = chem_df.merge(lb, on=["Cell", "Target"], how="left")
        gene_df = gene_df.merge(lb, on=["Cell", "Target"], how="left")

        X_umap = np.vstack([chem_cent, gene_cent]).astype(np.float32)
        reducer = umap.UMAP(
            n_neighbors=cfg.UMAP_N_NEIGHBORS,
            min_dist=cfg.UMAP_MIN_DIST,
            random_state=cfg.UMAP_RANDOM_STATE
        )
        emb = reducer.fit_transform(X_umap)

        chem_df["UMAP1"] = emb[:len(chem_cent), 0]
        chem_df["UMAP2"] = emb[:len(chem_cent), 1]
        gene_df["UMAP1"] = emb[len(chem_cent):, 0]
        gene_df["UMAP2"] = emb[len(chem_cent):, 1]

        all_rows.append(chem_df)
        all_rows.append(gene_df)

    if not all_rows:
        print("[WARN] Step3 UMAP produced no rows.")
        return None

    out = pd.concat(all_rows, ignore_index=True)
    front = ["UMAP1", "UMAP2", "Track", "View", "PerturbType", "Cell", "Target", "N_Samples", "Performance_Class"]
    cols = [c for c in front if c in out.columns] + [c for c in out.columns if c not in front]
    return out[cols].copy()


# -----------------------------
# Main
# -----------------------------
def main():
    cfg = Task2Config()
    cfg.TIME_BINS = _default_time_bins()

    out_dir = os.path.join(cfg.BASE_DIR, cfg.OUT_DIR)
    ensure_dir(out_dir)

    print(">>> Loading Task1 outputs...")
    df_pair = load_pairwise_wide(cfg)
    df_retr = load_retrieval(cfg)

    print(f"    Pairwise rows : {len(df_pair)}")
    print(f"    Retrieval rows: {len(df_retr)}")

    # Step 1
    print("\n>>> [Step 1] Build L1 tidy + context aggregated + tests")
    inst_l1, ctx_l1 = step1_build_l1_tidy(cfg, df_retr, df_pair)

    inst_path = os.path.join(out_dir, "Step1_L1_Instance_Tidy.csv")
    ctx_path = os.path.join(out_dir, "Step1_L1_Context_Aggregated.csv")
    inst_l1.to_csv(inst_path, index=False)
    ctx_l1.to_csv(ctx_path, index=False)
    print(f"    Saved: {inst_path} ({len(inst_l1)})")
    print(f"    Saved: {ctx_path} ({len(ctx_l1)})")

    tests = paired_tests(ctx_l1)
    tests_path = os.path.join(out_dir, "Step1_Tests.csv")
    tests.to_csv(tests_path, index=False)
    print(f"    Saved: {tests_path} ({len(tests)})")

    # Step 2
    print("\n>>> [Step 2] Context labeling + diagnostics")
    labels, diag_sp = step2_context_labeling(cfg, inst_l1, ctx_l1, df_retr)
    labels_path = os.path.join(out_dir, "Step2_Context_Labels.csv")
    labels.to_csv(labels_path, index=False)
    print(f"    Saved: {labels_path} ({len(labels)})")
    print("    Performance_Class counts:", labels["Performance_Class"].value_counts(dropna=False).to_dict())

    if diag_sp is not None and not diag_sp.empty:
        diag_path = os.path.join(out_dir, "Step2_Diagnostics_SourcePairLevel.csv")
        diag_sp.to_csv(diag_path, index=False)
        print(f"    Saved: {diag_path} ({len(diag_sp)})")

    # Step 4
    print("\n>>> [Step 4] Case-study tracer (A/B/C)")
    tracer = step4_case_study_tracer(cfg, labels, df_retr)
    tracer_path = os.path.join(out_dir, "Step4_CaseStudy_Tracer.csv")
    tracer.to_csv(tracer_path, index=False)
    print(f"    Saved: {tracer_path} ({len(tracer)})")

    # Step 5: enrichment + protocol correlations
    print("\n>>> [Step 5] Enrichment + protocol correlations")

    enr_targets_high = run_fisher_enrichment(
        labels, group_col="Target", class_col="Performance_Class",
        positive_class="Robust_High",
        min_contexts=cfg.MIN_CONTEXTS_FOR_ENRICHMENT,
        alternative="greater"
    )
    enr_targets_low = run_fisher_enrichment(
        labels, group_col="Target", class_col="Performance_Class",
        positive_class="Robust_Low",
        min_contexts=cfg.MIN_CONTEXTS_FOR_ENRICHMENT,
        alternative="greater"
    )
    enr_targets_ps = run_fisher_enrichment(
        labels, group_col="Target", class_col="Performance_Class",
        positive_class="Protocol_Sensitive",
        min_contexts=cfg.MIN_CONTEXTS_FOR_ENRICHMENT,
        alternative="greater"
    )
    enr_cells_high = run_fisher_enrichment(
        labels, group_col="Cell", class_col="Performance_Class",
        positive_class="Robust_High",
        min_contexts=cfg.MIN_CONTEXTS_FOR_ENRICHMENT,
        alternative="greater"
    )

    enr_targets_high.to_csv(os.path.join(out_dir, "Step5_Enrichment_Targets_RobustHigh.csv"), index=False)
    enr_targets_low.to_csv(os.path.join(out_dir, "Step5_Enrichment_Targets_RobustLow.csv"), index=False)
    enr_targets_ps.to_csv(os.path.join(out_dir, "Step5_Enrichment_Targets_ProtocolSensitive.csv"), index=False)
    enr_cells_high.to_csv(os.path.join(out_dir, "Step5_Enrichment_Cells_RobustHigh.csv"), index=False)
    print(f"    Enrichment saved under {out_dir}")

    proto = protocol_correlations_by_class(cfg, inst_l1, labels,
                                          class_col="Performance_Class",
                                          focus_class="Protocol_Sensitive")
    proto_path = os.path.join(out_dir, "Step5_Protocol_Correlations_ProtocolSensitive.csv")
    proto.to_csv(proto_path, index=False)
    print(f"    Saved: {proto_path} ({len(proto)})")

    # Step 3 optional
    if cfg.RUN_UMAP:
        umap_df = step3_dual_umap(cfg, labels)
        if umap_df is not None and not umap_df.empty:
            umap_path = os.path.join(out_dir, "Step3_UMAP_GroupCentroids.csv")
            umap_df.to_csv(umap_path, index=False)
            print(f"    Saved: {umap_path} ({len(umap_df)})")

    print("\n✅ Task 2 unified analysis complete.")


if __name__ == "__main__":
    main()
