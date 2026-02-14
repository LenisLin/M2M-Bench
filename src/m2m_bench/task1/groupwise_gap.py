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
import torch
from scipy.spatial.distance import cdist, pdist


@dataclass
class Task1Config:
    task0_unified_meta_path: str = "./outputs/task0_curated/metadata/unified_meta.parquet"
    processed_dir: str = "/mnt/NAS_21T/ProjectData/OSMOSIS/processed"
    lincs_tensor_name: str = "LINCS_Engine1_TrainData.pt"
    output_dir: str = "./outputs/task1"
    max_pairs_total: Optional[int] = 200000
    min_group_n_for_stats: int = 20
    max_groups_per_type: int = 300
    n_perm: int = 200
    n_bootstrap: int = 200
    max_n_for_null: int = 1500
    max_n_for_edist: int = 800
    random_seed: int = 42
    make_qc_plots: bool = True
    enable_strict_protocol_subset: bool = True
    strict_chemical_max_dose_logdiff: Optional[float] = 0.5
    strict_chemical_max_time_absdiff: Optional[float] = 24.0

    @property
    def lincs_data_path(self) -> str:
        return str(Path(self.processed_dir) / "LINCS_Processed" / self.lincs_tensor_name)

    @property
    def sc_dir(self) -> str:
        return str(Path(self.processed_dir) / "scPerturb_Processed")

    @property
    def data_dir(self) -> str:
        return str(Path(self.output_dir) / "data")

    @property
    def tensor_dir(self) -> str:
        return str(Path(self.data_dir) / "m1_tensors")

    @property
    def qc_table_dir(self) -> str:
        return str(Path(self.output_dir) / "qc" / "tables")

    @property
    def qc_fig_dir(self) -> str:
        return str(Path(self.output_dir) / "qc" / "figures")

    @property
    def analysis_dir(self) -> str:
        return str(Path(self.output_dir) / "analysis")

    @property
    def strict_analysis_dir(self) -> str:
        return str(Path(self.analysis_dir) / "strict")

    @property
    def run_manifest_path(self) -> str:
        return str(Path(self.output_dir) / "run_manifest_task1.json")


def _atomic_json_save(obj: dict, path: str) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = str(p) + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, p)


def _safe_parquet_write(df: pd.DataFrame, path: Path) -> bool:
    try:
        df.to_parquet(path, index=False)
        return True
    except Exception as exc:
        print(f"[WARN] Failed to write parquet {path}: {exc}")
        return False


def _ensure_dirs(cfg: Task1Config) -> None:
    for path in [cfg.data_dir, cfg.tensor_dir, cfg.qc_table_dir, cfg.qc_fig_dir, cfg.analysis_dir, cfg.strict_analysis_dir]:
        Path(path).mkdir(parents=True, exist_ok=True)


def _normalize_numeric(series: pd.Series) -> pd.Series:
    return pd.to_numeric(series, errors="coerce")


def _load_task1_index(cfg: Task1Config) -> pd.DataFrame:
    path = Path(cfg.task0_unified_meta_path)
    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path)
    else:
        try:
            df = pd.read_parquet(path)
        except Exception as exc:
            csv_fallback = path.with_suffix(".csv")
            if csv_fallback.exists():
                df = pd.read_csv(csv_fallback)
            else:
                raise RuntimeError(
                    f"Failed to read task0 unified metadata from {path}. "
                    f"Parquet engine missing and CSV fallback not found at {csv_fallback}."
                ) from exc
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
    missing = [c for c in need_cols if c not in df.columns]
    if missing:
        raise ValueError(f"Task0 unified_meta missing required columns: {missing}")

    out = df[need_cols].copy()
    out = out[out["source_db"].isin(["LINCS", "scPerturb"])].copy()
    out = out[out["modality"].isin(["Chemical", "Genetic"])].copy()
    out["dose_val"] = _normalize_numeric(out["dose_val"])
    out["time_val"] = _normalize_numeric(out["time_val"])
    out = out.reset_index(drop=True)
    out["task1_row_id"] = np.arange(len(out), dtype=np.int64)
    return out


def _build_m1_candidates(task1_index: pd.DataFrame) -> pd.DataFrame:
    key_cols = ["cell_std", "modality", "target_std"]
    counts = task1_index.groupby(key_cols + ["source_db"]).size().reset_index(name="n")
    pivot = counts.pivot_table(index=key_cols, columns="source_db", values="n", fill_value=0).reset_index()
    if "LINCS" not in pivot.columns:
        pivot["LINCS"] = 0
    if "scPerturb" not in pivot.columns:
        pivot["scPerturb"] = 0
    valid_keys = pivot[(pivot["LINCS"] > 0) & (pivot["scPerturb"] > 0)][key_cols]
    out = task1_index.merge(valid_keys, on=key_cols, how="inner")
    out = out.reset_index(drop=True)
    return out


def _pair_cost(
    dose_a: float,
    time_a: float,
    dose_b_arr: np.ndarray,
    time_b_arr: np.ndarray,
) -> np.ndarray:
    if np.isnan(dose_a):
        dose_term = np.full_like(dose_b_arr, 5.0, dtype=np.float64)
    else:
        dose_term = np.where(
            np.isnan(dose_b_arr),
            5.0,
            np.abs(np.log1p(np.clip(dose_a, a_min=0.0, a_max=None)) - np.log1p(np.clip(dose_b_arr, a_min=0.0, a_max=None))),
        )

    if np.isnan(time_a):
        time_term = np.full_like(time_b_arr, 5.0, dtype=np.float64)
    else:
        time_term = np.where(np.isnan(time_b_arr), 5.0, np.abs(time_a - time_b_arr) / 24.0)

    return dose_term + time_term


def _match_group(group_df: pd.DataFrame, rng: np.random.Generator) -> list[dict]:
    key_vals = {
        "cell_std": str(group_df["cell_std"].iloc[0]),
        "modality": str(group_df["modality"].iloc[0]),
        "target_std": str(group_df["target_std"].iloc[0]),
    }

    lincs = group_df[group_df["source_db"] == "LINCS"].copy().reset_index(drop=True)
    sc = group_df[group_df["source_db"] == "scPerturb"].copy().reset_index(drop=True)
    if len(lincs) == 0 or len(sc) == 0:
        return []

    matched: list[dict] = []
    l_unused = np.ones(len(lincs), dtype=bool)
    s_unused = np.ones(len(sc), dtype=bool)

    lincs["cond_id"] = lincs["cond_id"].astype(str)
    sc["cond_id"] = sc["cond_id"].astype(str)
    common_cond = sorted(set(lincs["cond_id"]).intersection(set(sc["cond_id"])))
    for cond in common_cond:
        li = np.where(l_unused & (lincs["cond_id"].values == cond))[0]
        si = np.where(s_unused & (sc["cond_id"].values == cond))[0]
        if len(li) == 0 or len(si) == 0:
            continue
        li = li.copy()
        si = si.copy()
        rng.shuffle(li)
        rng.shuffle(si)
        n = min(len(li), len(si))
        for i in range(n):
            l_idx = int(li[i])
            s_idx = int(si[i])
            l_unused[l_idx] = False
            s_unused[s_idx] = False
            matched.append(
                {
                    **key_vals,
                    "lincs_row_id": int(lincs.loc[l_idx, "task1_row_id"]),
                    "sc_row_id": int(sc.loc[s_idx, "task1_row_id"]),
                    "match_type": "exact_cond",
                    "match_score": 0.0,
                }
            )

    l_left = np.where(l_unused)[0]
    s_left = np.where(s_unused)[0]
    if len(l_left) == 0 or len(s_left) == 0:
        return matched

    if len(l_left) <= len(s_left):
        small_side = "lincs"
        small_idx = l_left
        big_idx = s_left
    else:
        small_side = "sc"
        small_idx = s_left
        big_idx = l_left

    big_available = np.ones(len(big_idx), dtype=bool)
    if small_side == "lincs":
        dose_big = sc.iloc[big_idx]["dose_val"].to_numpy(dtype=float)
        time_big = sc.iloc[big_idx]["time_val"].to_numpy(dtype=float)
    else:
        dose_big = lincs.iloc[big_idx]["dose_val"].to_numpy(dtype=float)
        time_big = lincs.iloc[big_idx]["time_val"].to_numpy(dtype=float)

    for sid in small_idx:
        avail = np.where(big_available)[0]
        if len(avail) == 0:
            break
        if small_side == "lincs":
            dose_s = float(lincs.loc[sid, "dose_val"])
            time_s = float(lincs.loc[sid, "time_val"])
        else:
            dose_s = float(sc.loc[sid, "dose_val"])
            time_s = float(sc.loc[sid, "time_val"])

        costs = _pair_cost(dose_s, time_s, dose_big[avail], time_big[avail])
        j_local = int(np.argmin(costs))
        j = int(avail[j_local])
        big_available[j] = False
        score = float(costs[j_local])
        if small_side == "lincs":
            l_idx = int(sid)
            s_idx = int(big_idx[j])
        else:
            l_idx = int(big_idx[j])
            s_idx = int(sid)
        matched.append(
            {
                **key_vals,
                "lincs_row_id": int(lincs.loc[l_idx, "task1_row_id"]),
                "sc_row_id": int(sc.loc[s_idx, "task1_row_id"]),
                "match_type": "nearest_dose_time",
                "match_score": score,
            }
        )

    return matched


def _match_m1_pairs(m1_candidates: pd.DataFrame, cfg: Task1Config, rng: np.random.Generator) -> pd.DataFrame:
    rows: list[dict] = []
    key_cols = ["cell_std", "modality", "target_std"]
    for _, grp in m1_candidates.groupby(key_cols, sort=False):
        rows.extend(_match_group(grp, rng=rng))
    out = pd.DataFrame(rows)
    if out.empty:
        raise RuntimeError("No matched pairs built for Task 1.")

    out = out.sort_values(["cell_std", "modality", "target_std", "match_type", "match_score"]).reset_index(drop=True)
    out["pair_id_task1"] = np.arange(len(out), dtype=np.int64)

    if cfg.max_pairs_total is not None and len(out) > cfg.max_pairs_total:
        keep = np.sort(rng.choice(len(out), size=cfg.max_pairs_total, replace=False))
        out = out.iloc[keep].reset_index(drop=True)

    return out


def _fetch_embeddings_for_rows(meta_rows: pd.DataFrame, cfg: Task1Config) -> tuple[np.ndarray, np.ndarray]:
    n = len(meta_rows)
    y_gene = [None] * n
    y_path = [None] * n

    lincs_rows = meta_rows[meta_rows["source_db"] == "LINCS"]
    if len(lincs_rows) > 0:
        lincs_data = torch.load(cfg.lincs_data_path, map_location="cpu")
        idx = torch.as_tensor(lincs_rows["global_idx"].to_numpy(dtype=int), dtype=torch.long)
        gene_sub = lincs_data["y_delta_gene"].index_select(0, idx).detach().cpu().numpy().astype(np.float32)
        path_sub = lincs_data["y_delta_pathway"].index_select(0, idx).detach().cpu().numpy().astype(np.float32)
        for i, rid in enumerate(lincs_rows["task1_row_id"].to_numpy(dtype=int)):
            pos = int(meta_rows.index[meta_rows["task1_row_id"] == rid][0])
            y_gene[pos] = gene_sub[i]
            y_path[pos] = path_sub[i]

    sc_rows = meta_rows[meta_rows["source_db"] == "scPerturb"]
    if len(sc_rows) > 0:
        for chunk_file, grp in sc_rows.groupby("chunk_file", sort=False):
            chunk_path = str(Path(cfg.sc_dir) / str(chunk_file))
            if not os.path.exists(chunk_path):
                raise FileNotFoundError(chunk_path)
            chunk = torch.load(chunk_path, map_location="cpu")
            idx = torch.as_tensor(grp["chunk_idx"].to_numpy(dtype=int), dtype=torch.long)
            gene_sub = chunk["y_delta_gene"].index_select(0, idx).detach().cpu().numpy().astype(np.float32)
            path_sub = chunk["y_delta_pathway"].index_select(0, idx).detach().cpu().numpy().astype(np.float32)
            for i, rid in enumerate(grp["task1_row_id"].to_numpy(dtype=int)):
                pos = int(meta_rows.index[meta_rows["task1_row_id"] == rid][0])
                y_gene[pos] = gene_sub[i]
                y_path[pos] = path_sub[i]

    if any(v is None for v in y_gene) or any(v is None for v in y_path):
        raise RuntimeError("Failed to retrieve some embeddings for matched rows.")

    return np.stack(y_gene, axis=0), np.stack(y_path, axis=0)


def _pairwise_cosine(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    num = np.sum(a * b, axis=1)
    den = np.linalg.norm(a, axis=1) * np.linalg.norm(b, axis=1)
    den = np.where(den == 0, np.nan, den)
    return num / den


def _pairwise_l2(a: np.ndarray, b: np.ndarray) -> np.ndarray:
    return np.linalg.norm(a - b, axis=1)


def _add_protocol_gap_columns(per_pair: pd.DataFrame) -> pd.DataFrame:
    """Attach explicit dose/time mismatch columns used in strict-subset analysis."""
    out = per_pair.copy()
    out["dose_logdiff"] = np.abs(
        np.log1p(np.clip(pd.to_numeric(out["dose_val_lincs"], errors="coerce"), a_min=0.0, a_max=None))
        - np.log1p(np.clip(pd.to_numeric(out["dose_val_sc"], errors="coerce"), a_min=0.0, a_max=None))
    )
    out["time_absdiff"] = np.abs(
        pd.to_numeric(out["time_val_lincs"], errors="coerce") - pd.to_numeric(out["time_val_sc"], errors="coerce")
    )
    return out


def _strict_protocol_mask(per_pair: pd.DataFrame, cfg: Task1Config) -> pd.Series:
    """
    Build strict protocol subset mask.
    - Genetic: keep exact condition matches.
    - Chemical: keep nearest-dose-time pairs under configured dose/time thresholds.
    """
    is_genetic_exact = (per_pair["modality"].astype(str) == "Genetic") & (
        per_pair["match_type"].astype(str) == "exact_cond"
    )

    is_chemical_nearest = (per_pair["modality"].astype(str) == "Chemical") & (
        per_pair["match_type"].astype(str) == "nearest_dose_time"
    )
    is_chemical_strict = is_chemical_nearest
    if cfg.strict_chemical_max_dose_logdiff is not None:
        is_chemical_strict = is_chemical_strict & (
            pd.to_numeric(per_pair["dose_logdiff"], errors="coerce") <= float(cfg.strict_chemical_max_dose_logdiff)
        )
    if cfg.strict_chemical_max_time_absdiff is not None:
        is_chemical_strict = is_chemical_strict & (
            pd.to_numeric(per_pair["time_absdiff"], errors="coerce") <= float(cfg.strict_chemical_max_time_absdiff)
        )

    return (is_genetic_exact | is_chemical_strict).fillna(False)


def _multivariate_energy_distance(x: np.ndarray, y: np.ndarray) -> float:
    if len(x) < 2 or len(y) < 2:
        return float("nan")
    exy = float(cdist(x, y, metric="euclidean").mean())
    exx = float(pdist(x, metric="euclidean").mean())
    eyy = float(pdist(y, metric="euclidean").mean())
    return 2.0 * exy - exx - eyy


def _bootstrap_ci_mean(values: np.ndarray, n_bootstrap: int, rng: np.random.Generator) -> tuple[float, float]:
    clean = values[np.isfinite(values)]
    if len(clean) == 0:
        return float("nan"), float("nan")
    means = np.empty(n_bootstrap, dtype=np.float64)
    n = len(clean)
    for i in range(n_bootstrap):
        idx = rng.integers(0, n, size=n)
        means[i] = float(np.mean(clean[idx]))
    return float(np.quantile(means, 0.025)), float(np.quantile(means, 0.975))


def _perm_null_mean_cos(
    x: np.ndarray,
    y: np.ndarray,
    n_perm: int,
    rng: np.random.Generator,
) -> np.ndarray:
    out = np.empty(n_perm, dtype=np.float64)
    n = len(x)
    for i in range(n_perm):
        perm = rng.permutation(n)
        out[i] = float(np.nanmean(_pairwise_cosine(x, y[perm])))
    return out


def _perm_null_edist(
    x: np.ndarray,
    y: np.ndarray,
    n_perm: int,
    rng: np.random.Generator,
) -> np.ndarray:
    out = np.empty(n_perm, dtype=np.float64)
    n = len(x)
    z = np.concatenate([x, y], axis=0)
    for i in range(n_perm):
        perm = rng.permutation(2 * n)
        xa = z[perm[:n]]
        ya = z[perm[n:]]
        out[i] = _multivariate_energy_distance(xa, ya)
    return out


def _downsample_embeddings(x: np.ndarray, y: np.ndarray, max_n: int, rng: np.random.Generator) -> tuple[np.ndarray, np.ndarray]:
    n = len(x)
    if n <= max_n:
        return x, y
    idx = np.sort(rng.choice(n, size=max_n, replace=False))
    return x[idx], y[idx]


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


def _write_qc_tables(index_df: pd.DataFrame, matched_df: pd.DataFrame, cfg: Task1Config) -> None:
    row_counts = index_df.groupby(["source_db", "modality"]).size().reset_index(name="n_rows")
    row_counts.to_csv(Path(cfg.qc_table_dir) / "qc_row_counts.csv", index=False)

    m_counts = matched_df.groupby(["modality", "match_type"]).size().reset_index(name="n_pairs")
    m_counts.to_csv(Path(cfg.qc_table_dir) / "qc_match_counts.csv", index=False)

    quality = (
        matched_df.groupby(["modality", "match_type"])["match_score"]
        .agg(["count", "mean", "median", "std", "min", "max"])
        .reset_index()
    )
    quality.to_csv(Path(cfg.qc_table_dir) / "qc_match_quality.csv", index=False)

    missingness = (
        index_df[["source_db", "modality", "dose_val", "time_val"]]
        .assign(
            dose_missing=lambda x: x["dose_val"].isna().astype(int),
            time_missing=lambda x: x["time_val"].isna().astype(int),
        )
        .groupby(["source_db", "modality"])[["dose_missing", "time_missing"]]
        .mean()
        .reset_index()
    )
    missingness.to_csv(Path(cfg.qc_table_dir) / "qc_missingness.csv", index=False)


def _write_qc_plots(index_df: pd.DataFrame, matched_df: pd.DataFrame, cfg: Task1Config) -> None:
    if not cfg.make_qc_plots:
        return
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    counts = index_df.groupby(["source_db", "modality"]).size().reset_index(name="n")
    fig, ax = plt.subplots(figsize=(8, 4))
    labels = counts["source_db"] + "|" + counts["modality"]
    ax.bar(labels, counts["n"])
    ax.set_title("Task1 index counts by source and modality")
    ax.set_ylabel("rows")
    ax.tick_params(axis="x", rotation=45)
    fig.tight_layout()
    fig.savefig(Path(cfg.qc_fig_dir) / "qc_counts_source_modality.png", dpi=150)
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(8, 4))
    ax.hist(matched_df["match_score"].dropna().to_numpy(), bins=60)
    ax.set_title("Task1 match score distribution")
    ax.set_xlabel("match_score")
    ax.set_ylabel("frequency")
    fig.tight_layout()
    fig.savefig(Path(cfg.qc_fig_dir) / "qc_match_score_hist.png", dpi=150)
    plt.close(fig)

    if "dose_val_lincs" in matched_df.columns and "dose_val_sc" in matched_df.columns:
        sub = matched_df.dropna(subset=["dose_val_lincs", "dose_val_sc"]).copy()
        if len(sub) > 0:
            sub = sub.sample(n=min(2000, len(sub)), random_state=cfg.random_seed)
            fig, ax = plt.subplots(figsize=(6, 6))
            ax.scatter(sub["dose_val_lincs"], sub["dose_val_sc"], s=5, alpha=0.5)
            ax.set_title("Dose matching sample (LINCS vs scPerturb)")
            ax.set_xlabel("dose_val_lincs")
            ax.set_ylabel("dose_val_sc")
            fig.tight_layout()
            fig.savefig(Path(cfg.qc_fig_dir) / "qc_dose_scatter.png", dpi=150)
            plt.close(fig)


def _run_groupwise_stats(
    per_pair: pd.DataFrame,
    l_gene: np.ndarray,
    s_gene: np.ndarray,
    l_path: np.ndarray,
    s_path: np.ndarray,
    cfg: Task1Config,
    rng: np.random.Generator,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, pd.DataFrame]:
    summary_rows: list[dict] = []
    strata_rows: list[dict] = []
    null_rows: list[dict] = []

    group_specs = [("global", None), ("modality", "modality")]
    for gname, gcol in group_specs:
        if gcol is None:
            groups = [("__ALL__", per_pair.index.to_numpy())]
        else:
            groups = []
            for gv, idx in per_pair.groupby(gcol).groups.items():
                groups.append((str(gv), np.array(list(idx), dtype=int)))

        for gv, idx in groups:
            n = len(idx)
            if n < cfg.min_group_n_for_stats:
                continue

            xg = l_gene[idx]
            yg = s_gene[idx]
            xp = l_path[idx]
            yp = s_path[idx]

            xg_null, yg_null = _downsample_embeddings(xg, yg, cfg.max_n_for_null, rng)
            xp_null, yp_null = _downsample_embeddings(xp, yp, cfg.max_n_for_null, rng)
            xg_ed, yg_ed = _downsample_embeddings(xg, yg, cfg.max_n_for_edist, rng)
            xp_ed, yp_ed = _downsample_embeddings(xp, yp, cfg.max_n_for_edist, rng)

            obs_cos_gene = float(np.nanmean(per_pair.loc[idx, "cosine_gene"].to_numpy()))
            obs_cos_path = float(np.nanmean(per_pair.loc[idx, "cosine_path"].to_numpy()))
            obs_l2_gene = float(np.nanmean(per_pair.loc[idx, "l2_gene"].to_numpy()))
            obs_l2_path = float(np.nanmean(per_pair.loc[idx, "l2_path"].to_numpy()))
            obs_edist_gene = float(_multivariate_energy_distance(xg_ed, yg_ed))
            obs_edist_path = float(_multivariate_energy_distance(xp_ed, yp_ed))

            null_cos_gene = _perm_null_mean_cos(xg_null, yg_null, cfg.n_perm, rng)
            null_cos_path = _perm_null_mean_cos(xp_null, yp_null, cfg.n_perm, rng)
            null_edist_gene = _perm_null_edist(xg_ed, yg_ed, cfg.n_perm, rng)
            null_edist_path = _perm_null_edist(xp_ed, yp_ed, cfg.n_perm, rng)

            p_cos_gene = float((np.sum(null_cos_gene >= obs_cos_gene) + 1) / (len(null_cos_gene) + 1))
            p_cos_path = float((np.sum(null_cos_path >= obs_cos_path) + 1) / (len(null_cos_path) + 1))
            p_ed_gene = float((np.sum(null_edist_gene <= obs_edist_gene) + 1) / (len(null_edist_gene) + 1))
            p_ed_path = float((np.sum(null_edist_path <= obs_edist_path) + 1) / (len(null_edist_path) + 1))

            z_cos_gene = float((obs_cos_gene - np.mean(null_cos_gene)) / (np.std(null_cos_gene) + 1e-12))
            z_cos_path = float((obs_cos_path - np.mean(null_cos_path)) / (np.std(null_cos_path) + 1e-12))
            z_ed_gene = float((obs_edist_gene - np.mean(null_edist_gene)) / (np.std(null_edist_gene) + 1e-12))
            z_ed_path = float((obs_edist_path - np.mean(null_edist_path)) / (np.std(null_edist_path) + 1e-12))

            ci_cos_gene = _bootstrap_ci_mean(per_pair.loc[idx, "cosine_gene"].to_numpy(), cfg.n_bootstrap, rng)
            ci_cos_path = _bootstrap_ci_mean(per_pair.loc[idx, "cosine_path"].to_numpy(), cfg.n_bootstrap, rng)
            ci_l2_gene = _bootstrap_ci_mean(per_pair.loc[idx, "l2_gene"].to_numpy(), cfg.n_bootstrap, rng)
            ci_l2_path = _bootstrap_ci_mean(per_pair.loc[idx, "l2_path"].to_numpy(), cfg.n_bootstrap, rng)

            summary_rows.append(
                {
                    "group_type": gname,
                    "group_value": gv,
                    "n_pairs": int(n),
                    "mean_cosine_gene": obs_cos_gene,
                    "mean_cosine_path": obs_cos_path,
                    "mean_l2_gene": obs_l2_gene,
                    "mean_l2_path": obs_l2_path,
                    "edist_gene": obs_edist_gene,
                    "edist_path": obs_edist_path,
                    "ci95_cosine_gene_lo": ci_cos_gene[0],
                    "ci95_cosine_gene_hi": ci_cos_gene[1],
                    "ci95_cosine_path_lo": ci_cos_path[0],
                    "ci95_cosine_path_hi": ci_cos_path[1],
                    "ci95_l2_gene_lo": ci_l2_gene[0],
                    "ci95_l2_gene_hi": ci_l2_gene[1],
                    "ci95_l2_path_lo": ci_l2_path[0],
                    "ci95_l2_path_hi": ci_l2_path[1],
                    "p_cosine_gene": p_cos_gene,
                    "p_cosine_path": p_cos_path,
                    "p_edist_gene": p_ed_gene,
                    "p_edist_path": p_ed_path,
                    "z_cosine_gene": z_cos_gene,
                    "z_cosine_path": z_cos_path,
                    "z_edist_gene": z_ed_gene,
                    "z_edist_path": z_ed_path,
                }
            )

            null_rows.append(
                {
                    "group_type": gname,
                    "group_value": gv,
                    "metric": "null_cosine_gene",
                    "null_mean": float(np.mean(null_cos_gene)),
                    "null_std": float(np.std(null_cos_gene)),
                    "null_q05": float(np.quantile(null_cos_gene, 0.05)),
                    "null_q95": float(np.quantile(null_cos_gene, 0.95)),
                }
            )
            null_rows.append(
                {
                    "group_type": gname,
                    "group_value": gv,
                    "metric": "null_cosine_path",
                    "null_mean": float(np.mean(null_cos_path)),
                    "null_std": float(np.std(null_cos_path)),
                    "null_q05": float(np.quantile(null_cos_path, 0.05)),
                    "null_q95": float(np.quantile(null_cos_path, 0.95)),
                }
            )
            null_rows.append(
                {
                    "group_type": gname,
                    "group_value": gv,
                    "metric": "null_edist_gene",
                    "null_mean": float(np.mean(null_edist_gene)),
                    "null_std": float(np.std(null_edist_gene)),
                    "null_q05": float(np.quantile(null_edist_gene, 0.05)),
                    "null_q95": float(np.quantile(null_edist_gene, 0.95)),
                }
            )
            null_rows.append(
                {
                    "group_type": gname,
                    "group_value": gv,
                    "metric": "null_edist_path",
                    "null_mean": float(np.mean(null_edist_path)),
                    "null_std": float(np.std(null_edist_path)),
                    "null_q05": float(np.quantile(null_edist_path, 0.05)),
                    "null_q95": float(np.quantile(null_edist_path, 0.95)),
                }
            )

    for gtype in ["cell_std", "target_std"]:
        counts = per_pair[gtype].value_counts().head(cfg.max_groups_per_type)
        keep = set(counts.index.astype(str).tolist())
        for gv, grp in per_pair[per_pair[gtype].astype(str).isin(keep)].groupby(gtype):
            idx = grp.index.to_numpy()
            if len(idx) < cfg.min_group_n_for_stats:
                continue
            strata_rows.append(
                {
                    "group_type": gtype,
                    "group_value": str(gv),
                    "n_pairs": int(len(idx)),
                    "mean_cosine_gene": float(np.nanmean(grp["cosine_gene"].to_numpy())),
                    "mean_cosine_path": float(np.nanmean(grp["cosine_path"].to_numpy())),
                    "mean_l2_gene": float(np.nanmean(grp["l2_gene"].to_numpy())),
                    "mean_l2_path": float(np.nanmean(grp["l2_path"].to_numpy())),
                }
            )

    summary = pd.DataFrame(summary_rows)
    strata = pd.DataFrame(strata_rows)
    null_summary = pd.DataFrame(null_rows)

    if not summary.empty:
        for metric in ["p_cosine_gene", "p_cosine_path", "p_edist_gene", "p_edist_path"]:
            summary[f"{metric}_fdr"] = _safe_fdr_bh(summary[metric].to_numpy())

    significance = summary[
        [
            "group_type",
            "group_value",
            "n_pairs",
            "p_cosine_gene",
            "p_cosine_gene_fdr",
            "p_cosine_path",
            "p_cosine_path_fdr",
            "p_edist_gene",
            "p_edist_gene_fdr",
            "p_edist_path",
            "p_edist_path_fdr",
        ]
    ].copy()

    effect_size = summary[
        [
            "group_type",
            "group_value",
            "n_pairs",
            "mean_cosine_gene",
            "z_cosine_gene",
            "mean_cosine_path",
            "z_cosine_path",
            "edist_gene",
            "z_edist_gene",
            "edist_path",
            "z_edist_path",
            "mean_l2_gene",
            "mean_l2_path",
        ]
    ].copy()

    return summary, strata, significance, effect_size, null_summary


def run_task1(cfg: Task1Config) -> None:
    t0 = time.time()
    rng = np.random.default_rng(cfg.random_seed)
    _ensure_dirs(cfg)

    task1_index = _load_task1_index(cfg)
    _safe_parquet_write(task1_index, Path(cfg.data_dir) / "task1_index.parquet")

    m1_candidates = _build_m1_candidates(task1_index)
    _safe_parquet_write(m1_candidates, Path(cfg.data_dir) / "m1_candidates.parquet")
    m1_candidates.to_csv(Path(cfg.data_dir) / "m1_candidates.csv", index=False)

    matched = _match_m1_pairs(m1_candidates, cfg=cfg, rng=rng)
    index_cols = [
        "task1_row_id",
        "source_db",
        "dose_val",
        "time_val",
        "cond_id",
        "pert_type",
        "cell_std",
        "modality",
        "target_std",
        "global_idx",
        "chunk_file",
        "chunk_idx",
    ]
    lookup = task1_index[index_cols].copy()

    lmeta = lookup.rename(columns=lambda c: f"{c}_lincs" if c != "task1_row_id" else c)
    smeta = lookup.rename(columns=lambda c: f"{c}_sc" if c != "task1_row_id" else c)
    matched = matched.merge(lmeta, left_on="lincs_row_id", right_on="task1_row_id", how="left").drop(columns=["task1_row_id"])
    matched = matched.merge(smeta, left_on="sc_row_id", right_on="task1_row_id", how="left").drop(columns=["task1_row_id"])
    _safe_parquet_write(matched, Path(cfg.data_dir) / "m1_matched_pairs.parquet")
    matched.to_csv(Path(cfg.data_dir) / "m1_matched_pairs.csv", index=False)

    all_ids = np.unique(np.concatenate([matched["lincs_row_id"].to_numpy(), matched["sc_row_id"].to_numpy()]))
    meta_rows = task1_index[task1_index["task1_row_id"].isin(all_ids)].copy().reset_index(drop=True)
    y_gene, y_path = _fetch_embeddings_for_rows(meta_rows, cfg=cfg)

    row_to_pos = dict(zip(meta_rows["task1_row_id"].to_numpy(dtype=int), np.arange(len(meta_rows), dtype=int)))
    lidx = matched["lincs_row_id"].map(row_to_pos).to_numpy(dtype=int)
    sidx = matched["sc_row_id"].map(row_to_pos).to_numpy(dtype=int)
    l_gene = y_gene[lidx]
    s_gene = y_gene[sidx]
    l_path = y_path[lidx]
    s_path = y_path[sidx]

    np.save(Path(cfg.tensor_dir) / "m1_lincs_y_gene.npy", l_gene.astype(np.float32))
    np.save(Path(cfg.tensor_dir) / "m1_sc_y_gene.npy", s_gene.astype(np.float32))
    np.save(Path(cfg.tensor_dir) / "m1_lincs_y_path.npy", l_path.astype(np.float32))
    np.save(Path(cfg.tensor_dir) / "m1_sc_y_path.npy", s_path.astype(np.float32))

    per_pair = matched[
        [
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
        ]
    ].copy()
    per_pair = _add_protocol_gap_columns(per_pair)
    per_pair["cosine_gene"] = _pairwise_cosine(l_gene, s_gene)
    per_pair["cosine_path"] = _pairwise_cosine(l_path, s_path)
    per_pair["l2_gene"] = _pairwise_l2(l_gene, s_gene)
    per_pair["l2_path"] = _pairwise_l2(l_path, s_path)

    gene_mu = float(np.nanmean(per_pair["cosine_gene"].to_numpy()))
    gene_sd = float(np.nanstd(per_pair["cosine_gene"].to_numpy()) + 1e-12)
    path_mu = float(np.nanmean(per_pair["cosine_path"].to_numpy()))
    path_sd = float(np.nanstd(per_pair["cosine_path"].to_numpy()) + 1e-12)
    per_pair["cosine_gene_z"] = (per_pair["cosine_gene"] - gene_mu) / gene_sd
    per_pair["cosine_path_z"] = (per_pair["cosine_path"] - path_mu) / path_sd
    if cfg.enable_strict_protocol_subset:
        per_pair["is_strict_protocol_pair"] = _strict_protocol_mask(per_pair, cfg).astype(bool)
    else:
        per_pair["is_strict_protocol_pair"] = False
    per_pair.to_csv(Path(cfg.analysis_dir) / "modality_gap_per_pair.csv", index=False)

    _write_qc_tables(task1_index, matched, cfg)
    _write_qc_plots(task1_index, matched, cfg)

    summary, strata, significance, effect_size, null_summary = _run_groupwise_stats(
        per_pair=per_pair,
        l_gene=l_gene,
        s_gene=s_gene,
        l_path=l_path,
        s_path=s_path,
        cfg=cfg,
        rng=rng,
    )

    summary.to_csv(Path(cfg.analysis_dir) / "modality_gap_summary.csv", index=False)
    strata.to_csv(Path(cfg.analysis_dir) / "modality_gap_by_strata.csv", index=False)
    significance.to_csv(Path(cfg.analysis_dir) / "significance_results.csv", index=False)
    effect_size.to_csv(Path(cfg.analysis_dir) / "effect_size_table.csv", index=False)
    null_summary.to_csv(Path(cfg.analysis_dir) / "null_distribution_summary.csv", index=False)

    strict_summary = pd.DataFrame()
    if cfg.enable_strict_protocol_subset:
        strict_per_pair = per_pair[per_pair["is_strict_protocol_pair"]].copy().reset_index(drop=True)
        strict_per_pair.to_csv(Path(cfg.strict_analysis_dir) / "modality_gap_per_pair_strict.csv", index=False)
        strict_counts = (
            strict_per_pair.groupby(["modality", "match_type"], dropna=False).size().reset_index(name="n_pairs")
            if not strict_per_pair.empty
            else pd.DataFrame(columns=["modality", "match_type", "n_pairs"])
        )
        strict_counts.to_csv(Path(cfg.qc_table_dir) / "qc_strict_subset_counts.csv", index=False)

        if len(strict_per_pair) >= cfg.min_group_n_for_stats:
            strict_mask = per_pair["is_strict_protocol_pair"].to_numpy(dtype=bool)
            strict_l_gene = l_gene[strict_mask]
            strict_s_gene = s_gene[strict_mask]
            strict_l_path = l_path[strict_mask]
            strict_s_path = s_path[strict_mask]
            strict_summary, strict_strata, strict_significance, strict_effect_size, strict_null_summary = _run_groupwise_stats(
                per_pair=strict_per_pair,
                l_gene=strict_l_gene,
                s_gene=strict_s_gene,
                l_path=strict_l_path,
                s_path=strict_s_path,
                cfg=cfg,
                rng=rng,
            )
            strict_summary.to_csv(Path(cfg.strict_analysis_dir) / "modality_gap_summary_strict.csv", index=False)
            strict_strata.to_csv(Path(cfg.strict_analysis_dir) / "modality_gap_by_strata_strict.csv", index=False)
            strict_significance.to_csv(Path(cfg.strict_analysis_dir) / "significance_results_strict.csv", index=False)
            strict_effect_size.to_csv(Path(cfg.strict_analysis_dir) / "effect_size_table_strict.csv", index=False)
            strict_null_summary.to_csv(Path(cfg.strict_analysis_dir) / "null_distribution_summary_strict.csv", index=False)

    run_manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": {
            "n_task1_index": int(len(task1_index)),
            "n_m1_candidates": int(len(m1_candidates)),
            "n_m1_matched_pairs": int(len(matched)),
            "n_strict_protocol_pairs": int(per_pair["is_strict_protocol_pair"].sum()),
            "strict_stats_written": bool(not strict_summary.empty),
        },
        "outputs": {
            "data_dir": cfg.data_dir,
            "qc_table_dir": cfg.qc_table_dir,
            "qc_fig_dir": cfg.qc_fig_dir,
            "analysis_dir": cfg.analysis_dir,
            "strict_analysis_dir": cfg.strict_analysis_dir,
        },
    }
    _atomic_json_save(run_manifest, cfg.run_manifest_path)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("M2M-Bench Task1 group-wise modality gap analysis")
    parser.add_argument("--task0-unified-meta-path", type=str, default=Task1Config.task0_unified_meta_path)
    parser.add_argument("--processed-dir", type=str, default=Task1Config.processed_dir)
    parser.add_argument("--lincs-tensor-name", type=str, default=Task1Config.lincs_tensor_name)
    parser.add_argument("--output-dir", type=str, default=Task1Config.output_dir)
    parser.add_argument("--max-pairs-total", type=int, default=Task1Config.max_pairs_total)
    parser.add_argument("--min-group-n-for-stats", type=int, default=Task1Config.min_group_n_for_stats)
    parser.add_argument("--max-groups-per-type", type=int, default=Task1Config.max_groups_per_type)
    parser.add_argument("--n-perm", type=int, default=Task1Config.n_perm)
    parser.add_argument("--n-bootstrap", type=int, default=Task1Config.n_bootstrap)
    parser.add_argument("--max-n-for-null", type=int, default=Task1Config.max_n_for_null)
    parser.add_argument("--max-n-for-edist", type=int, default=Task1Config.max_n_for_edist)
    parser.add_argument("--random-seed", type=int, default=Task1Config.random_seed)
    parser.add_argument("--make-qc-plots", action=argparse.BooleanOptionalAction, default=Task1Config.make_qc_plots)
    parser.add_argument(
        "--enable-strict-protocol-subset",
        action=argparse.BooleanOptionalAction,
        default=Task1Config.enable_strict_protocol_subset,
        help="Write an additional strict protocol-matched subset analysis.",
    )
    parser.add_argument(
        "--strict-chemical-max-dose-logdiff",
        type=float,
        default=Task1Config.strict_chemical_max_dose_logdiff,
        help="Chemical strict threshold for abs(log1p dose diff).",
    )
    parser.add_argument(
        "--strict-chemical-max-time-absdiff",
        type=float,
        default=Task1Config.strict_chemical_max_time_absdiff,
        help="Chemical strict threshold for absolute time difference (hours).",
    )
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    cfg = Task1Config(
        task0_unified_meta_path=args.task0_unified_meta_path,
        processed_dir=args.processed_dir,
        lincs_tensor_name=args.lincs_tensor_name,
        output_dir=args.output_dir,
        max_pairs_total=args.max_pairs_total,
        min_group_n_for_stats=args.min_group_n_for_stats,
        max_groups_per_type=args.max_groups_per_type,
        n_perm=args.n_perm,
        n_bootstrap=args.n_bootstrap,
        max_n_for_null=args.max_n_for_null,
        max_n_for_edist=args.max_n_for_edist,
        random_seed=args.random_seed,
        make_qc_plots=args.make_qc_plots,
        enable_strict_protocol_subset=args.enable_strict_protocol_subset,
        strict_chemical_max_dose_logdiff=args.strict_chemical_max_dose_logdiff,
        strict_chemical_max_time_absdiff=args.strict_chemical_max_time_absdiff,
    )
    run_task1(cfg)


if __name__ == "__main__":
    main()
