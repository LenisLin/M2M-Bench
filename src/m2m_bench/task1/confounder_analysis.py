from __future__ import annotations

import argparse
import json
import os
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd


@dataclass
class Task1ConfounderConfig:
    per_pair_path: str = "./outputs/task1/analysis/modality_gap_per_pair.csv"
    matched_pairs_path: str = "./outputs/task1/data/m1_matched_pairs.parquet"
    output_dir: str = "./outputs/task1/confounder"
    n_perm: int = 200
    min_group_n: int = 10
    max_groups_for_report: int = 100
    random_seed: int = 42
    make_qc_plots: bool = True

    @property
    def analysis_dir(self) -> str:
        return str(Path(self.output_dir) / "analysis")

    @property
    def qc_dir(self) -> str:
        return str(Path(self.output_dir) / "qc")

    @property
    def run_manifest_path(self) -> str:
        return str(Path(self.output_dir) / "run_manifest_task1_confounder.json")


def _atomic_json_save(obj: dict, path: str) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = str(p) + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, p)


def _ensure_dirs(cfg: Task1ConfounderConfig) -> None:
    Path(cfg.analysis_dir).mkdir(parents=True, exist_ok=True)
    Path(cfg.qc_dir).mkdir(parents=True, exist_ok=True)


def _load_table(path_like: str) -> pd.DataFrame:
    path = Path(path_like)
    if not path.exists():
        raise FileNotFoundError(path)
    if path.suffix.lower() == ".csv":
        return pd.read_csv(path)
    try:
        return pd.read_parquet(path)
    except Exception as exc:
        csv_fallback = path.with_suffix(".csv")
        if csv_fallback.exists():
            return pd.read_csv(csv_fallback)
        raise RuntimeError(
            f"Failed to read {path}. Parquet engine may be missing and CSV fallback not found at {csv_fallback}."
        ) from exc


def _safe_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _prepare_dataset(cfg: Task1ConfounderConfig) -> pd.DataFrame:
    per_pair = _load_table(cfg.per_pair_path).copy()
    need_per_pair = [
        "pair_id_task1",
        "cell_std",
        "modality",
        "target_std",
        "match_type",
        "match_score",
        "dose_val_lincs",
        "dose_val_sc",
        "time_val_lincs",
        "time_val_sc",
        "cosine_gene",
        "cosine_path",
        "l2_gene",
        "l2_path",
    ]
    missing = [c for c in need_per_pair if c not in per_pair.columns]
    if missing:
        raise ValueError(f"per_pair input missing required columns: {missing}")

    matched = _load_table(cfg.matched_pairs_path).copy()
    if "pair_id_task1" not in matched.columns:
        raise ValueError("matched_pairs input missing required column: pair_id_task1")

    keep_cols = [
        "pair_id_task1",
        "pert_type_lincs",
        "pert_type_sc",
        "cond_id_lincs",
        "cond_id_sc",
    ]
    keep_cols = [c for c in keep_cols if c in matched.columns]
    merged = per_pair.merge(matched[keep_cols], on="pair_id_task1", how="left")

    merged["dose_val_lincs"] = _safe_numeric(merged["dose_val_lincs"])
    merged["dose_val_sc"] = _safe_numeric(merged["dose_val_sc"])
    merged["time_val_lincs"] = _safe_numeric(merged["time_val_lincs"])
    merged["time_val_sc"] = _safe_numeric(merged["time_val_sc"])
    merged["cosine_gene"] = _safe_numeric(merged["cosine_gene"])
    merged["cosine_path"] = _safe_numeric(merged["cosine_path"])
    merged["l2_gene"] = _safe_numeric(merged["l2_gene"])
    merged["l2_path"] = _safe_numeric(merged["l2_path"])
    merged["match_score"] = _safe_numeric(merged["match_score"])

    merged["gap_gene"] = 1.0 - merged["cosine_gene"]
    merged["gap_path"] = 1.0 - merged["cosine_path"]
    merged["dose_logdiff"] = np.abs(np.log1p(np.clip(merged["dose_val_lincs"], a_min=0.0, a_max=None)) - np.log1p(np.clip(merged["dose_val_sc"], a_min=0.0, a_max=None)))
    merged["time_absdiff"] = np.abs(merged["time_val_lincs"] - merged["time_val_sc"])
    merged["is_exact_cond"] = (merged["match_type"].astype(str) == "exact_cond").astype(int)
    merged["context_key"] = (
        merged["cell_std"].astype(str) + "||" + merged["modality"].astype(str) + "||" + merged["target_std"].astype(str)
    )

    return merged


def _summarize_by_groups(df: pd.DataFrame, cfg: Task1ConfounderConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    main_rows: list[dict] = []
    for gcol in ["modality", "match_type"]:
        for gval, grp in df.groupby(gcol, sort=False):
            if len(grp) < cfg.min_group_n:
                continue
            main_rows.append(
                {
                    "group_type": gcol,
                    "group_value": str(gval),
                    "n_pairs": int(len(grp)),
                    "mean_gap_gene": float(np.nanmean(grp["gap_gene"].to_numpy())),
                    "mean_gap_path": float(np.nanmean(grp["gap_path"].to_numpy())),
                    "mean_cosine_gene": float(np.nanmean(grp["cosine_gene"].to_numpy())),
                    "mean_cosine_path": float(np.nanmean(grp["cosine_path"].to_numpy())),
                    "mean_l2_gene": float(np.nanmean(grp["l2_gene"].to_numpy())),
                    "mean_l2_path": float(np.nanmean(grp["l2_path"].to_numpy())),
                    "mean_dose_logdiff": float(np.nanmean(grp["dose_logdiff"].to_numpy())),
                    "mean_time_absdiff": float(np.nanmean(grp["time_absdiff"].to_numpy())),
                }
            )

    strata_rows: list[dict] = []
    for gcol in ["cell_std", "target_std"]:
        counts = df[gcol].value_counts().head(cfg.max_groups_for_report)
        keep = set(counts.index.astype(str).tolist())
        sub = df[df[gcol].astype(str).isin(keep)]
        for gval, grp in sub.groupby(gcol, sort=False):
            if len(grp) < cfg.min_group_n:
                continue
            strata_rows.append(
                {
                    "group_type": gcol,
                    "group_value": str(gval),
                    "n_pairs": int(len(grp)),
                    "mean_gap_gene": float(np.nanmean(grp["gap_gene"].to_numpy())),
                    "mean_gap_path": float(np.nanmean(grp["gap_path"].to_numpy())),
                    "mean_cosine_gene": float(np.nanmean(grp["cosine_gene"].to_numpy())),
                    "mean_cosine_path": float(np.nanmean(grp["cosine_path"].to_numpy())),
                }
            )

    return pd.DataFrame(main_rows), pd.DataFrame(strata_rows)


def _eta_squared(y: np.ndarray, groups: np.ndarray) -> float:
    valid = np.isfinite(y) & np.asarray(pd.notna(groups))
    yv = y[valid]
    gv = groups[valid]
    if len(yv) < 3:
        return float("nan")
    grand = float(np.mean(yv))
    sst = float(np.sum((yv - grand) ** 2))
    if sst <= 0:
        return float("nan")
    ssb = 0.0
    for g in np.unique(gv):
        yi = yv[gv == g]
        if len(yi) == 0:
            continue
        ssb += float(len(yi) * (np.mean(yi) - grand) ** 2)
    return float(ssb / sst)


def _perm_test_eta(y: np.ndarray, groups: np.ndarray, n_perm: int, rng: np.random.Generator) -> tuple[float, float]:
    obs = _eta_squared(y, groups)
    if not np.isfinite(obs):
        return float("nan"), float("nan")
    null = np.empty(n_perm, dtype=np.float64)
    for i in range(n_perm):
        g_perm = rng.permutation(groups)
        null[i] = _eta_squared(y, g_perm)
    p_ge = float((np.sum(null >= obs) + 1) / (n_perm + 1))
    return obs, p_ge


def _abs_spearman(y: np.ndarray, x: np.ndarray) -> float:
    valid = np.isfinite(y) & np.isfinite(x)
    if np.sum(valid) < 5:
        return float("nan")
    yv = y[valid]
    xv = x[valid]
    ry = pd.Series(yv).rank(method="average").to_numpy()
    rx = pd.Series(xv).rank(method="average").to_numpy()
    if np.std(ry) <= 0 or np.std(rx) <= 0:
        return float("nan")
    corr = np.corrcoef(ry, rx)[0, 1]
    return float(np.abs(corr))


def _perm_test_abs_spearman(y: np.ndarray, x: np.ndarray, n_perm: int, rng: np.random.Generator) -> tuple[float, float]:
    obs = _abs_spearman(y, x)
    if not np.isfinite(obs):
        return float("nan"), float("nan")
    null = np.empty(n_perm, dtype=np.float64)
    valid = np.isfinite(y) & np.isfinite(x)
    yv = y[valid]
    xv = x[valid]
    for i in range(n_perm):
        yp = rng.permutation(yv)
        null[i] = _abs_spearman(yp, xv)
    p_ge = float((np.sum(null >= obs) + 1) / (n_perm + 1))
    return obs, p_ge


def _run_factor_tests(df: pd.DataFrame, cfg: Task1ConfounderConfig, rng: np.random.Generator) -> pd.DataFrame:
    outcomes = ["gap_gene", "gap_path"]
    categorical_factors = ["modality", "match_type", "cell_std", "target_std"]
    numeric_factors = ["dose_logdiff", "time_absdiff", "match_score"]

    rows: list[dict] = []
    for outcome in outcomes:
        y = _safe_numeric(df[outcome]).to_numpy(dtype=float)

        for factor in categorical_factors:
            g = df[factor].astype(str).to_numpy()
            if len(np.unique(g)) < 2:
                continue
            obs, p_ge = _perm_test_eta(y, g, cfg.n_perm, rng)
            rows.append(
                {
                    "outcome": outcome,
                    "factor": factor,
                    "factor_type": "categorical",
                    "effect_stat": obs,
                    "effect_name": "eta_squared",
                    "p_ge_null": p_ge,
                }
            )

        for factor in numeric_factors:
            x = _safe_numeric(df[factor]).to_numpy(dtype=float)
            obs, p_ge = _perm_test_abs_spearman(y, x, cfg.n_perm, rng)
            rows.append(
                {
                    "outcome": outcome,
                    "factor": factor,
                    "factor_type": "continuous",
                    "effect_stat": obs,
                    "effect_name": "abs_spearman",
                    "p_ge_null": p_ge,
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out = out.sort_values(["outcome", "p_ge_null", "effect_stat"], ascending=[True, True, False]).reset_index(drop=True)
        out["p_ge_null_fdr"] = np.nan
        for outcome, idx in out.groupby("outcome").groups.items():
            p = out.loc[list(idx), "p_ge_null"].to_numpy(dtype=float)
            q = _safe_fdr_bh(p)
            out.loc[list(idx), "p_ge_null_fdr"] = q
    return out


def _safe_fdr_bh(pvals: np.ndarray) -> np.ndarray:
    p = np.array(pvals, dtype=float)
    n = len(p)
    order = np.argsort(np.nan_to_num(p, nan=1.0))
    ranked = p[order]
    q = np.full(n, np.nan, dtype=float)
    prev = 1.0
    for i in range(n - 1, -1, -1):
        rank = i + 1
        if np.isnan(ranked[i]):
            cur = np.nan
        else:
            cur = min(prev, ranked[i] * n / rank)
        q[order[i]] = cur
        if not np.isnan(cur):
            prev = cur
    return q


def _build_design_matrix(
    df: pd.DataFrame,
    categorical_factors: list[str],
    numeric_factors: list[str],
) -> tuple[np.ndarray, np.ndarray, dict[str, list[int]]]:
    x_parts = [np.ones((len(df), 1), dtype=np.float64)]
    term_names = ["Intercept"]
    factor_to_cols: dict[str, list[int]] = {"Intercept": [0]}
    col_cursor = 1

    for factor in numeric_factors:
        x = _safe_numeric(df[factor]).to_numpy(dtype=float)
        fill = np.nanmedian(x) if np.any(np.isfinite(x)) else 0.0
        x = np.where(np.isfinite(x), x, fill)
        x = x.reshape(-1, 1).astype(np.float64)
        x_parts.append(x)
        term_names.append(factor)
        factor_to_cols[factor] = [col_cursor]
        col_cursor += 1

    for factor in categorical_factors:
        d = pd.get_dummies(df[factor].astype(str), prefix=factor, drop_first=True, dtype=float)
        if d.shape[1] == 0:
            continue
        x = d.to_numpy(dtype=np.float64)
        x_parts.append(x)
        cols = list(range(col_cursor, col_cursor + x.shape[1]))
        factor_to_cols[factor] = cols
        term_names.extend(d.columns.tolist())
        col_cursor += x.shape[1]

    x_all = np.concatenate(x_parts, axis=1)
    return x_all, np.array(term_names, dtype=object), factor_to_cols


def _fit_ols(x: np.ndarray, y: np.ndarray) -> tuple[np.ndarray, np.ndarray, float, float]:
    beta, _, _, _ = np.linalg.lstsq(x, y, rcond=None)
    pred = x @ beta
    resid = y - pred
    rss = float(np.sum(resid**2))
    tss = float(np.sum((y - np.mean(y)) ** 2))
    r2 = float(1.0 - rss / tss) if tss > 0 else float("nan")
    return beta, pred, rss, r2


def _model_decomposition(df: pd.DataFrame, outcome: str) -> tuple[pd.DataFrame, pd.DataFrame]:
    work = df.copy()
    y = _safe_numeric(work[outcome]).to_numpy(dtype=float)
    valid = np.isfinite(y)
    work = work.loc[valid].reset_index(drop=True)
    y = y[valid]
    if len(work) < 20:
        return pd.DataFrame(), pd.DataFrame()

    categorical = ["modality", "match_type", "cell_std", "target_std"]
    numeric = ["dose_logdiff", "time_absdiff", "match_score"]
    x, terms, f2c = _build_design_matrix(work, categorical_factors=categorical, numeric_factors=numeric)
    beta, _, rss_full, r2_full = _fit_ols(x, y)

    coef = pd.DataFrame(
        {
            "outcome": outcome,
            "term": terms,
            "coefficient": beta,
        }
    )

    contrib_rows: list[dict] = []
    factors = [f for f in numeric + categorical if f in f2c]
    for factor in factors:
        drop_cols = set(f2c.get(factor, []))
        keep_cols = [i for i in range(x.shape[1]) if i not in drop_cols]
        if 0 not in keep_cols:
            keep_cols = [0] + keep_cols
        xr = x[:, keep_cols]
        _, _, rss_red, r2_red = _fit_ols(xr, y)
        delta_r2 = float(max(0.0, r2_full - r2_red)) if np.isfinite(r2_full) and np.isfinite(r2_red) else float("nan")
        contrib_rows.append(
            {
                "outcome": outcome,
                "factor": factor,
                "n_obs": int(len(y)),
                "r2_full": r2_full,
                "r2_reduced": r2_red,
                "delta_r2": delta_r2,
                "rss_full": rss_full,
                "rss_reduced": rss_red,
            }
        )

    contrib = pd.DataFrame(contrib_rows)
    if not contrib.empty:
        total = float(np.nansum(np.clip(contrib["delta_r2"].to_numpy(dtype=float), a_min=0.0, a_max=None)))
        if total > 0:
            contrib["delta_r2_share"] = contrib["delta_r2"] / total
        else:
            contrib["delta_r2_share"] = np.nan
        contrib = contrib.sort_values("delta_r2", ascending=False).reset_index(drop=True)
    return coef, contrib


def _build_context_labels(df: pd.DataFrame, cfg: Task1ConfounderConfig) -> pd.DataFrame:
    summary = (
        df.groupby(["cell_std", "modality", "target_std"], sort=False)
        .agg(
            n_pairs=("pair_id_task1", "count"),
            mean_gap_gene=("gap_gene", "mean"),
            mean_gap_path=("gap_path", "mean"),
            mean_cosine_gene=("cosine_gene", "mean"),
            mean_cosine_path=("cosine_path", "mean"),
            mean_dose_logdiff=("dose_logdiff", "mean"),
            mean_time_absdiff=("time_absdiff", "mean"),
        )
        .reset_index()
    )
    summary = summary[summary["n_pairs"] >= cfg.min_group_n].copy()
    if summary.empty:
        return summary

    summary["combined_gap"] = 0.5 * summary["mean_gap_gene"] + 0.5 * summary["mean_gap_path"]
    labels = []
    for modality, grp in summary.groupby("modality", sort=False):
        q33 = float(np.nanquantile(grp["combined_gap"].to_numpy(), 0.33))
        q67 = float(np.nanquantile(grp["combined_gap"].to_numpy(), 0.67))
        for _, row in grp.iterrows():
            value = float(row["combined_gap"])
            if value <= q33:
                label = "Robust"
            elif value <= q67:
                label = "Conditional"
            else:
                label = "Unstable"
            labels.append(
                {
                    "cell_std": row["cell_std"],
                    "modality": modality,
                    "target_std": row["target_std"],
                    "n_pairs": int(row["n_pairs"]),
                    "mean_gap_gene": float(row["mean_gap_gene"]),
                    "mean_gap_path": float(row["mean_gap_path"]),
                    "mean_cosine_gene": float(row["mean_cosine_gene"]),
                    "mean_cosine_path": float(row["mean_cosine_path"]),
                    "combined_gap": value,
                    "label_rule_q33": q33,
                    "label_rule_q67": q67,
                    "reliability_label": label,
                }
            )
    out = pd.DataFrame(labels).sort_values(["modality", "combined_gap"]).reset_index(drop=True)
    return out


def _write_qc_plots(df: pd.DataFrame, cfg: Task1ConfounderConfig) -> None:
    if not cfg.make_qc_plots:
        return
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(df["gap_gene"].dropna().to_numpy(), bins=50, alpha=0.7, label="gap_gene")
    ax.hist(df["gap_path"].dropna().to_numpy(), bins=50, alpha=0.5, label="gap_path")
    ax.set_title("Task1 confounder input: gap distributions")
    ax.set_xlabel("gap value")
    ax.set_ylabel("frequency")
    ax.legend()
    fig.tight_layout()
    fig.savefig(Path(cfg.qc_dir) / "gap_distribution_hist.png", dpi=150)
    plt.close(fig)

    sub = df.dropna(subset=["dose_logdiff", "gap_gene"]).copy()
    if len(sub) > 0:
        sub = sub.sample(n=min(2000, len(sub)), random_state=cfg.random_seed)
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(sub["dose_logdiff"], sub["gap_gene"], s=6, alpha=0.5)
        ax.set_title("Dose mismatch vs gene gap")
        ax.set_xlabel("dose_logdiff")
        ax.set_ylabel("gap_gene")
        fig.tight_layout()
        fig.savefig(Path(cfg.qc_dir) / "dose_mismatch_vs_gap_gene.png", dpi=150)
        plt.close(fig)


def run_task1_confounder(cfg: Task1ConfounderConfig) -> None:
    t0 = time.time()
    rng = np.random.default_rng(cfg.random_seed)
    _ensure_dirs(cfg)

    df = _prepare_dataset(cfg)
    df.to_csv(Path(cfg.qc_dir) / "confounder_input_snapshot.csv", index=False)

    group_summary, strata_summary = _summarize_by_groups(df, cfg)
    group_summary.to_csv(Path(cfg.analysis_dir) / "group_summary.csv", index=False)
    strata_summary.to_csv(Path(cfg.analysis_dir) / "strata_summary.csv", index=False)

    factor_tests = _run_factor_tests(df, cfg, rng=rng)
    factor_tests.to_csv(Path(cfg.analysis_dir) / "factor_permutation_tests.csv", index=False)

    coef_gene, contrib_gene = _model_decomposition(df, outcome="gap_gene")
    coef_path, contrib_path = _model_decomposition(df, outcome="gap_path")
    coef_gene.to_csv(Path(cfg.analysis_dir) / "ols_coefficients_gap_gene.csv", index=False)
    coef_path.to_csv(Path(cfg.analysis_dir) / "ols_coefficients_gap_path.csv", index=False)
    contrib_gene.to_csv(Path(cfg.analysis_dir) / "variance_contribution_gap_gene.csv", index=False)
    contrib_path.to_csv(Path(cfg.analysis_dir) / "variance_contribution_gap_path.csv", index=False)

    labels = _build_context_labels(df, cfg)
    labels.to_csv(Path(cfg.analysis_dir) / "context_reliability_labels.csv", index=False)

    example_rows = []
    if not labels.empty:
        for tag in ["Robust", "Unstable"]:
            sub = labels[labels["reliability_label"] == tag].copy()
            if len(sub) == 0:
                continue
            picks = sub.head(5) if tag == "Robust" else sub.tail(5)
            picks = picks.copy()
            picks["example_type"] = tag.lower()
            example_rows.append(picks)
    if example_rows:
        pd.concat(example_rows, axis=0, ignore_index=True).to_csv(
            Path(cfg.analysis_dir) / "context_reliability_examples.csv", index=False
        )

    _write_qc_plots(df, cfg)

    run_manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": {
            "n_pairs": int(len(df)),
            "n_group_summary_rows": int(len(group_summary)),
            "n_factor_tests": int(len(factor_tests)),
            "n_context_labels": int(len(labels)),
        },
        "outputs": {
            "analysis_dir": cfg.analysis_dir,
            "qc_dir": cfg.qc_dir,
        },
    }
    _atomic_json_save(run_manifest, cfg.run_manifest_path)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("M2M-Bench Task1 confounder analysis")
    parser.add_argument("--per-pair-path", type=str, default=Task1ConfounderConfig.per_pair_path)
    parser.add_argument("--matched-pairs-path", type=str, default=Task1ConfounderConfig.matched_pairs_path)
    parser.add_argument("--output-dir", type=str, default=Task1ConfounderConfig.output_dir)
    parser.add_argument("--n-perm", type=int, default=Task1ConfounderConfig.n_perm)
    parser.add_argument("--min-group-n", type=int, default=Task1ConfounderConfig.min_group_n)
    parser.add_argument("--max-groups-for-report", type=int, default=Task1ConfounderConfig.max_groups_for_report)
    parser.add_argument("--random-seed", type=int, default=Task1ConfounderConfig.random_seed)
    parser.add_argument("--make-qc-plots", action=argparse.BooleanOptionalAction, default=Task1ConfounderConfig.make_qc_plots)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    cfg = Task1ConfounderConfig(
        per_pair_path=args.per_pair_path,
        matched_pairs_path=args.matched_pairs_path,
        output_dir=args.output_dir,
        n_perm=args.n_perm,
        min_group_n=args.min_group_n,
        max_groups_for_report=args.max_groups_for_report,
        random_seed=args.random_seed,
        make_qc_plots=args.make_qc_plots,
    )
    run_task1_confounder(cfg)


if __name__ == "__main__":
    main()
