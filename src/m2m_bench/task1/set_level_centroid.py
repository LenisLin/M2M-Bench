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

from .retrieval_instance import Task1RetrievalConfig, _fetch_embeddings_for_rows, _load_candidates


@dataclass
class Task1SetLevelConfig:
    m1_candidates_path: str = "./outputs/task1/data/m1_candidates.parquet"
    processed_dir: str = "/mnt/NAS_21T/ProjectData/OSMOSIS/processed"
    lincs_tensor_name: str = "LINCS_Engine1_TrainData.pt"
    output_dir: str = "./outputs/task1/set_level"
    per_pair_path: Optional[str] = "./outputs/task1/analysis/modality_gap_per_pair.csv"
    min_n_per_source: int = 5
    n_perm: int = 200
    n_bootstrap: int = 200
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
        return str(Path(self.output_dir) / "run_manifest_task1_set_level.json")


def _atomic_json_save(obj: dict, path: str) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = str(p) + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, p)


def _ensure_dirs(cfg: Task1SetLevelConfig) -> None:
    Path(cfg.analysis_dir).mkdir(parents=True, exist_ok=True)
    Path(cfg.qc_dir).mkdir(parents=True, exist_ok=True)


def _cosine(a: np.ndarray, b: np.ndarray) -> float:
    den = float(np.linalg.norm(a) * np.linalg.norm(b))
    if den <= 0:
        return float("nan")
    return float(np.dot(a, b) / den)


def _l2(a: np.ndarray, b: np.ndarray) -> float:
    return float(np.linalg.norm(a - b))


def _centroid_metrics(
    x_lincs: np.ndarray,
    x_sc: np.ndarray,
    n_perm: int,
    n_bootstrap: int,
    rng: np.random.Generator,
) -> dict:
    c_l = x_lincs.mean(axis=0)
    c_s = x_sc.mean(axis=0)
    obs_cos = _cosine(c_l, c_s)
    obs_l2 = _l2(c_l, c_s)

    # null by shuffled source assignment within same context (size-preserving)
    pool = np.concatenate([x_lincs, x_sc], axis=0)
    n_l = len(x_lincs)
    n_total = len(pool)
    null_cos = np.empty(n_perm, dtype=np.float64)
    for i in range(n_perm):
        perm = rng.permutation(n_total)
        a = pool[perm[:n_l]]
        b = pool[perm[n_l:]]
        null_cos[i] = _cosine(a.mean(axis=0), b.mean(axis=0))
    p_cos_ge_null = float((np.sum(null_cos >= obs_cos) + 1) / (n_perm + 1))
    z_cos = float((obs_cos - np.mean(null_cos)) / (np.std(null_cos) + 1e-12))

    # bootstrap CI for centroid cosine
    boot_cos = np.empty(n_bootstrap, dtype=np.float64)
    for i in range(n_bootstrap):
        idx_l = rng.integers(0, len(x_lincs), size=len(x_lincs))
        idx_s = rng.integers(0, len(x_sc), size=len(x_sc))
        bl = x_lincs[idx_l].mean(axis=0)
        bs = x_sc[idx_s].mean(axis=0)
        boot_cos[i] = _cosine(bl, bs)
    ci_lo = float(np.quantile(boot_cos, 0.025))
    ci_hi = float(np.quantile(boot_cos, 0.975))

    disp_l = float(np.mean(np.linalg.norm(x_lincs - c_l, axis=1)))
    disp_s = float(np.mean(np.linalg.norm(x_sc - c_s, axis=1)))

    return {
        "cosine_centroid": obs_cos,
        "l2_centroid": obs_l2,
        "p_cosine_ge_null": p_cos_ge_null,
        "z_cosine": z_cos,
        "ci95_cosine_lo": ci_lo,
        "ci95_cosine_hi": ci_hi,
        "within_disp_lincs": disp_l,
        "within_disp_sc": disp_s,
        "null_cosine_mean": float(np.mean(null_cos)),
        "null_cosine_std": float(np.std(null_cos)),
    }


def _summarize_modalities(context_df: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for modality, grp in context_df.groupby("modality", sort=False):
        w = grp["n_total"].to_numpy(dtype=float)
        w = w / np.sum(w)
        rows.append(
            {
                "group_type": "modality",
                "group_value": str(modality),
                "n_contexts": int(len(grp)),
                "mean_cosine_gene_unweighted": float(np.nanmean(grp["cosine_centroid_gene"].to_numpy())),
                "mean_cosine_path_unweighted": float(np.nanmean(grp["cosine_centroid_path"].to_numpy())),
                "mean_cosine_gene_weighted": float(np.nansum(grp["cosine_centroid_gene"].to_numpy() * w)),
                "mean_cosine_path_weighted": float(np.nansum(grp["cosine_centroid_path"].to_numpy() * w)),
                "mean_l2_gene_unweighted": float(np.nanmean(grp["l2_centroid_gene"].to_numpy())),
                "mean_l2_path_unweighted": float(np.nanmean(grp["l2_centroid_path"].to_numpy())),
            }
        )

    rows.append(
        {
            "group_type": "global",
            "group_value": "__ALL__",
            "n_contexts": int(len(context_df)),
            "mean_cosine_gene_unweighted": float(np.nanmean(context_df["cosine_centroid_gene"].to_numpy())),
            "mean_cosine_path_unweighted": float(np.nanmean(context_df["cosine_centroid_path"].to_numpy())),
            "mean_cosine_gene_weighted": float(np.nansum(context_df["cosine_centroid_gene"] * context_df["n_total"]) / np.sum(context_df["n_total"])),
            "mean_cosine_path_weighted": float(np.nansum(context_df["cosine_centroid_path"] * context_df["n_total"]) / np.sum(context_df["n_total"])),
            "mean_l2_gene_unweighted": float(np.nanmean(context_df["l2_centroid_gene"].to_numpy())),
            "mean_l2_path_unweighted": float(np.nanmean(context_df["l2_centroid_path"].to_numpy())),
        }
    )
    return pd.DataFrame(rows)


def _compare_with_pair_level(context_df: pd.DataFrame, per_pair_path: Optional[str], out_dir: Path) -> tuple[pd.DataFrame, dict]:
    if per_pair_path is None:
        return pd.DataFrame(), {}
    p = Path(per_pair_path)
    if not p.exists():
        return pd.DataFrame(), {}
    per_pair = pd.read_csv(p)
    need = ["cell_std", "modality", "target_std", "cosine_gene", "cosine_path"]
    if not all(c in per_pair.columns for c in need):
        return pd.DataFrame(), {}

    pair_ctx = (
        per_pair.groupby(["cell_std", "modality", "target_std"], sort=False)
        .agg(
            n_pairs=("cosine_gene", "count"),
            mean_pair_cosine_gene=("cosine_gene", "mean"),
            mean_pair_cosine_path=("cosine_path", "mean"),
        )
        .reset_index()
    )

    merged = context_df.merge(pair_ctx, on=["cell_std", "modality", "target_std"], how="left")
    merged["delta_gene_centroid_minus_pair"] = merged["cosine_centroid_gene"] - merged["mean_pair_cosine_gene"]
    merged["delta_path_centroid_minus_pair"] = merged["cosine_centroid_path"] - merged["mean_pair_cosine_path"]
    merged.to_csv(out_dir / "context_centroid_vs_pairwise.csv", index=False)

    stats = {
        "n_overlap_contexts": int(np.sum(np.isfinite(merged["mean_pair_cosine_gene"].to_numpy()))),
        "corr_gene": float(pd.Series(merged["cosine_centroid_gene"]).corr(pd.Series(merged["mean_pair_cosine_gene"]))),
        "corr_path": float(pd.Series(merged["cosine_centroid_path"]).corr(pd.Series(merged["mean_pair_cosine_path"]))),
        "mean_delta_gene": float(np.nanmean(merged["delta_gene_centroid_minus_pair"].to_numpy())),
        "mean_delta_path": float(np.nanmean(merged["delta_path_centroid_minus_pair"].to_numpy())),
    }
    return merged, stats


def _write_qc_plots(context_df: pd.DataFrame, compare_df: pd.DataFrame, cfg: Task1SetLevelConfig) -> None:
    if not cfg.make_qc_plots:
        return
    try:
        import matplotlib.pyplot as plt
    except Exception:
        return

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(context_df["cosine_centroid_gene"].dropna().to_numpy(), bins=40, alpha=0.7, label="centroid_gene")
    ax.hist(context_df["cosine_centroid_path"].dropna().to_numpy(), bins=40, alpha=0.5, label="centroid_path")
    ax.set_title("Set-level centroid cosine distributions")
    ax.set_xlabel("cosine")
    ax.set_ylabel("n_context")
    ax.legend()
    fig.tight_layout()
    fig.savefig(Path(cfg.qc_dir) / "centroid_cosine_hist.png", dpi=160)
    plt.close(fig)

    if not compare_df.empty:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(
            compare_df["mean_pair_cosine_gene"].to_numpy(),
            compare_df["cosine_centroid_gene"].to_numpy(),
            s=10,
            alpha=0.7,
            label="gene",
        )
        ax.scatter(
            compare_df["mean_pair_cosine_path"].to_numpy(),
            compare_df["cosine_centroid_path"].to_numpy(),
            s=10,
            alpha=0.7,
            label="path",
        )
        ax.set_xlabel("pairwise context mean cosine")
        ax.set_ylabel("centroid cosine")
        ax.set_title("Set-level centroid vs pairwise context means")
        ax.legend()
        fig.tight_layout()
        fig.savefig(Path(cfg.qc_dir) / "centroid_vs_pair_scatter.png", dpi=160)
        plt.close(fig)


def run_task1_set_level(cfg: Task1SetLevelConfig) -> None:
    t0 = time.time()
    rng = np.random.default_rng(cfg.random_seed)
    _ensure_dirs(cfg)

    retrieval_cfg = Task1RetrievalConfig(
        m1_candidates_path=cfg.m1_candidates_path,
        processed_dir=cfg.processed_dir,
        lincs_tensor_name=cfg.lincs_tensor_name,
        output_dir="./tmp_not_used",
    )
    candidate_source_path = retrieval_cfg.m1_candidates_path
    try:
        candidates = _load_candidates(retrieval_cfg)
    except Exception:
        fallback = Path("./outputs/task1/retrieval/qc/retrieval_candidates_snapshot.csv")
        if not fallback.exists():
            raise
        retrieval_cfg = Task1RetrievalConfig(
            m1_candidates_path=str(fallback),
            processed_dir=cfg.processed_dir,
            lincs_tensor_name=cfg.lincs_tensor_name,
            output_dir="./tmp_not_used",
        )
        candidates = _load_candidates(retrieval_cfg)
        candidate_source_path = str(fallback)
    y_gene, y_path = _fetch_embeddings_for_rows(candidates, retrieval_cfg)

    rows = []
    key_cols = ["cell_std", "modality", "target_std"]
    for keys, grp in candidates.groupby(key_cols, sort=False):
        cell_std, modality, target_std = keys
        l = grp[grp["source_db"] == "LINCS"]
        s = grp[grp["source_db"] == "scPerturb"]
        if len(l) < cfg.min_n_per_source or len(s) < cfg.min_n_per_source:
            continue

        lidx = l["retrieval_row_idx"].to_numpy(dtype=int)
        sidx = s["retrieval_row_idx"].to_numpy(dtype=int)
        gene_stats = _centroid_metrics(
            x_lincs=y_gene[lidx],
            x_sc=y_gene[sidx],
            n_perm=cfg.n_perm,
            n_bootstrap=cfg.n_bootstrap,
            rng=rng,
        )
        path_stats = _centroid_metrics(
            x_lincs=y_path[lidx],
            x_sc=y_path[sidx],
            n_perm=cfg.n_perm,
            n_bootstrap=cfg.n_bootstrap,
            rng=rng,
        )

        rows.append(
            {
                "cell_std": str(cell_std),
                "modality": str(modality),
                "target_std": str(target_std),
                "n_lincs": int(len(l)),
                "n_sc": int(len(s)),
                "n_total": int(len(l) + len(s)),
                "cosine_centroid_gene": gene_stats["cosine_centroid"],
                "l2_centroid_gene": gene_stats["l2_centroid"],
                "p_cosine_gene_ge_null": gene_stats["p_cosine_ge_null"],
                "z_cosine_gene": gene_stats["z_cosine"],
                "ci95_cosine_gene_lo": gene_stats["ci95_cosine_lo"],
                "ci95_cosine_gene_hi": gene_stats["ci95_cosine_hi"],
                "within_disp_gene_lincs": gene_stats["within_disp_lincs"],
                "within_disp_gene_sc": gene_stats["within_disp_sc"],
                "null_cosine_gene_mean": gene_stats["null_cosine_mean"],
                "null_cosine_gene_std": gene_stats["null_cosine_std"],
                "cosine_centroid_path": path_stats["cosine_centroid"],
                "l2_centroid_path": path_stats["l2_centroid"],
                "p_cosine_path_ge_null": path_stats["p_cosine_ge_null"],
                "z_cosine_path": path_stats["z_cosine"],
                "ci95_cosine_path_lo": path_stats["ci95_cosine_lo"],
                "ci95_cosine_path_hi": path_stats["ci95_cosine_hi"],
                "within_disp_path_lincs": path_stats["within_disp_lincs"],
                "within_disp_path_sc": path_stats["within_disp_sc"],
                "null_cosine_path_mean": path_stats["null_cosine_mean"],
                "null_cosine_path_std": path_stats["null_cosine_std"],
            }
        )

    context_df = pd.DataFrame(rows)
    if context_df.empty:
        raise RuntimeError("No contexts met min_n_per_source for set-level analysis.")
    context_df.to_csv(Path(cfg.analysis_dir) / "set_level_context_metrics.csv", index=False)

    summary = _summarize_modalities(context_df)
    summary.to_csv(Path(cfg.analysis_dir) / "set_level_summary.csv", index=False)

    compare_df, compare_stats = _compare_with_pair_level(
        context_df=context_df,
        per_pair_path=cfg.per_pair_path,
        out_dir=Path(cfg.analysis_dir),
    )
    if compare_stats:
        pd.DataFrame([compare_stats]).to_csv(
            Path(cfg.analysis_dir) / "set_level_vs_pairwise_stats.csv", index=False
        )

    # Simple reliability labels from set-level combined centroid gap
    rel = context_df.copy()
    rel["combined_gap_centroid"] = 1.0 - 0.5 * (
        rel["cosine_centroid_gene"].to_numpy() + rel["cosine_centroid_path"].to_numpy()
    )
    label_rows = []
    for modality, grp in rel.groupby("modality", sort=False):
        q33 = float(np.nanquantile(grp["combined_gap_centroid"].to_numpy(), 0.33))
        q67 = float(np.nanquantile(grp["combined_gap_centroid"].to_numpy(), 0.67))
        for _, row in grp.iterrows():
            val = float(row["combined_gap_centroid"])
            if val <= q33:
                lab = "Robust"
            elif val <= q67:
                lab = "Conditional"
            else:
                lab = "Unstable"
            label_rows.append(
                {
                    "cell_std": row["cell_std"],
                    "modality": row["modality"],
                    "target_std": row["target_std"],
                    "n_lincs": int(row["n_lincs"]),
                    "n_sc": int(row["n_sc"]),
                    "cosine_centroid_gene": float(row["cosine_centroid_gene"]),
                    "cosine_centroid_path": float(row["cosine_centroid_path"]),
                    "combined_gap_centroid": val,
                    "q33_rule": q33,
                    "q67_rule": q67,
                    "reliability_label": lab,
                }
            )
    pd.DataFrame(label_rows).to_csv(Path(cfg.analysis_dir) / "set_level_reliability_labels.csv", index=False)

    _write_qc_plots(context_df=context_df, compare_df=compare_df, cfg=cfg)

    run_manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": {
            "candidate_source_path": candidate_source_path,
            "n_candidates": int(len(candidates)),
            "n_contexts": int(len(context_df)),
            "n_contexts_chemical": int(np.sum(context_df["modality"] == "Chemical")),
            "n_contexts_genetic": int(np.sum(context_df["modality"] == "Genetic")),
        },
        "outputs": {
            "analysis_dir": cfg.analysis_dir,
            "qc_dir": cfg.qc_dir,
        },
    }
    _atomic_json_save(run_manifest, cfg.run_manifest_path)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("M2M-Bench Task1 set-level centroid analysis")
    parser.add_argument("--m1-candidates-path", type=str, default=Task1SetLevelConfig.m1_candidates_path)
    parser.add_argument("--processed-dir", type=str, default=Task1SetLevelConfig.processed_dir)
    parser.add_argument("--lincs-tensor-name", type=str, default=Task1SetLevelConfig.lincs_tensor_name)
    parser.add_argument("--output-dir", type=str, default=Task1SetLevelConfig.output_dir)
    parser.add_argument("--per-pair-path", type=str, default=Task1SetLevelConfig.per_pair_path)
    parser.add_argument("--min-n-per-source", type=int, default=Task1SetLevelConfig.min_n_per_source)
    parser.add_argument("--n-perm", type=int, default=Task1SetLevelConfig.n_perm)
    parser.add_argument("--n-bootstrap", type=int, default=Task1SetLevelConfig.n_bootstrap)
    parser.add_argument("--random-seed", type=int, default=Task1SetLevelConfig.random_seed)
    parser.add_argument("--make-qc-plots", action=argparse.BooleanOptionalAction, default=Task1SetLevelConfig.make_qc_plots)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    cfg = Task1SetLevelConfig(
        m1_candidates_path=args.m1_candidates_path,
        processed_dir=args.processed_dir,
        lincs_tensor_name=args.lincs_tensor_name,
        output_dir=args.output_dir,
        per_pair_path=args.per_pair_path,
        min_n_per_source=args.min_n_per_source,
        n_perm=args.n_perm,
        n_bootstrap=args.n_bootstrap,
        random_seed=args.random_seed,
        make_qc_plots=args.make_qc_plots,
    )
    run_task1_set_level(cfg)


if __name__ == "__main__":
    main()
