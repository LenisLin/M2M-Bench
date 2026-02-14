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


@dataclass
class Task1RetrievalConfig:
    m1_candidates_path: str = "./outputs/task1/data/m1_candidates.parquet"
    processed_dir: str = "/mnt/NAS_21T/ProjectData/OSMOSIS/processed"
    lincs_tensor_name: str = "LINCS_Engine1_TrainData.pt"
    output_dir: str = "./outputs/task1/retrieval"
    tracks: str = "gene,path"
    directions: str = "LINCS->scPerturb,scPerturb->LINCS"
    topk: str = "1,5,10"
    batch_size: int = 256
    n_perm: int = 200
    min_queries_per_group: int = 10
    random_seed: int = 42
    max_queries_per_group: Optional[int] = None
    max_gallery_per_group: Optional[int] = None
    run_balanced_eval: bool = True
    balanced_gallery_size: int = 256
    balanced_true_per_query: int = 1
    balanced_n_repeats: int = 50

    @property
    def lincs_data_path(self) -> str:
        return str(Path(self.processed_dir) / "LINCS_Processed" / self.lincs_tensor_name)

    @property
    def sc_dir(self) -> str:
        return str(Path(self.processed_dir) / "scPerturb_Processed")

    @property
    def analysis_dir(self) -> str:
        return str(Path(self.output_dir) / "analysis")

    @property
    def qc_dir(self) -> str:
        return str(Path(self.output_dir) / "qc")

    @property
    def run_manifest_path(self) -> str:
        return str(Path(self.output_dir) / "run_manifest_task1_retrieval.json")


def _atomic_json_save(obj: dict, path: str) -> None:
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    tmp = str(p) + ".tmp"
    with open(tmp, "w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp, p)


def _ensure_dirs(cfg: Task1RetrievalConfig) -> None:
    Path(cfg.analysis_dir).mkdir(parents=True, exist_ok=True)
    Path(cfg.qc_dir).mkdir(parents=True, exist_ok=True)


def _parse_list_str(x: str, cast=str) -> list:
    items = [v.strip() for v in str(x).split(",") if v.strip()]
    return [cast(v) for v in items]


def _parse_directions(x: str) -> list[tuple[str, str]]:
    items = _parse_list_str(x, str)
    out: list[tuple[str, str]] = []
    for item in items:
        if "->" not in item:
            raise ValueError(f"Invalid direction format: {item} (expected src->dst)")
        src, dst = item.split("->", 1)
        out.append((src.strip(), dst.strip()))
    return out


def _load_candidates(cfg: Task1RetrievalConfig) -> pd.DataFrame:
    path = Path(cfg.m1_candidates_path)
    if path.suffix.lower() == ".csv":
        df = pd.read_csv(path).copy()
    else:
        try:
            df = pd.read_parquet(path).copy()
        except Exception as exc:
            csv_fallback = path.with_suffix(".csv")
            if csv_fallback.exists():
                df = pd.read_csv(csv_fallback).copy()
            else:
                raise RuntimeError(
                    f"Failed to read candidates from {path}. "
                    f"Parquet engine may be missing and CSV fallback was not found at {csv_fallback}."
                ) from exc
    required = [
        "task1_row_id",
        "cell_std",
        "target_std",
        "modality",
        "source_db",
        "global_idx",
        "chunk_file",
        "chunk_idx",
        "dose_val",
        "time_val",
    ]
    missing = [c for c in required if c not in df.columns]
    if missing:
        raise ValueError(f"m1_candidates missing required columns: {missing}")

    df = df[df["source_db"].isin(["LINCS", "scPerturb"])].copy()
    df = df[df["modality"].isin(["Chemical", "Genetic"])].copy()
    df = df.reset_index(drop=True)
    df["retrieval_row_idx"] = np.arange(len(df), dtype=np.int64)
    df["label_key"] = (
        df["cell_std"].astype(str) + "||" + df["modality"].astype(str) + "||" + df["target_std"].astype(str)
    )
    return df


def _fetch_embeddings_for_rows(
    meta_rows: pd.DataFrame,
    cfg: Task1RetrievalConfig,
) -> tuple[np.ndarray, np.ndarray]:
    n = len(meta_rows)
    y_gene = [None] * n
    y_path = [None] * n

    lincs_rows = meta_rows[meta_rows["source_db"] == "LINCS"]
    if len(lincs_rows) > 0:
        lincs_data = torch.load(cfg.lincs_data_path, map_location="cpu")
        idx = torch.as_tensor(lincs_rows["global_idx"].to_numpy(dtype=int), dtype=torch.long)
        gene_sub = lincs_data["y_delta_gene"].index_select(0, idx).detach().cpu().numpy().astype(np.float32)
        path_sub = lincs_data["y_delta_pathway"].index_select(0, idx).detach().cpu().numpy().astype(np.float32)
        for i, ridx in enumerate(lincs_rows["retrieval_row_idx"].to_numpy(dtype=int)):
            pos = int(meta_rows.index[meta_rows["retrieval_row_idx"] == ridx][0])
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
            for i, ridx in enumerate(grp["retrieval_row_idx"].to_numpy(dtype=int)):
                pos = int(meta_rows.index[meta_rows["retrieval_row_idx"] == ridx][0])
                y_gene[pos] = gene_sub[i]
                y_path[pos] = path_sub[i]

    if any(v is None for v in y_gene) or any(v is None for v in y_path):
        raise RuntimeError("Failed to fetch all embeddings for retrieval candidates.")

    return np.stack(y_gene, axis=0), np.stack(y_path, axis=0)


def _cosine_sim(q: np.ndarray, g: np.ndarray) -> np.ndarray:
    qn = np.linalg.norm(q, axis=1, keepdims=True)
    gn = np.linalg.norm(g, axis=1, keepdims=True).T
    denom = np.clip(qn * gn, a_min=1e-12, a_max=None)
    return (q @ g.T) / denom


def _simulate_null_group(
    n_match_arr: np.ndarray,
    n_gallery: int,
    n_perm: int,
    rng: np.random.Generator,
    topk_values: list[int],
) -> dict:
    null_top1 = np.empty(n_perm, dtype=np.float64)
    null_mrr = np.empty(n_perm, dtype=np.float64)
    null_topk = {k: np.empty(n_perm, dtype=np.float64) for k in topk_values}

    for p in range(n_perm):
        ranks = np.empty(len(n_match_arr), dtype=np.int32)
        for i, n_match in enumerate(n_match_arr):
            pos = rng.choice(n_gallery, size=int(n_match), replace=False)
            ranks[i] = int(np.min(pos)) + 1
        mrr = np.mean(1.0 / ranks)
        top1 = np.mean(ranks <= 1)
        null_mrr[p] = float(mrr)
        null_top1[p] = float(top1)
        for k in topk_values:
            null_topk[k][p] = float(np.mean(ranks <= k))

    out = {
        "null_mrr_mean": float(np.mean(null_mrr)),
        "null_mrr_std": float(np.std(null_mrr)),
        "null_mrr_q05": float(np.quantile(null_mrr, 0.05)),
        "null_mrr_q95": float(np.quantile(null_mrr, 0.95)),
        "null_top1_mean": float(np.mean(null_top1)),
        "null_top1_std": float(np.std(null_top1)),
        "null_top1_q05": float(np.quantile(null_top1, 0.05)),
        "null_top1_q95": float(np.quantile(null_top1, 0.95)),
    }
    for k in topk_values:
        arr = null_topk[k]
        out[f"null_top{k}_mean"] = float(np.mean(arr))
        out[f"null_top{k}_std"] = float(np.std(arr))
        out[f"null_top{k}_q05"] = float(np.quantile(arr, 0.05))
        out[f"null_top{k}_q95"] = float(np.quantile(arr, 0.95))
    out["_null_arrays"] = {"mrr": null_mrr, "top1": null_top1, **{f"top{k}": null_topk[k] for k in topk_values}}
    return out


def _balanced_query_metrics(
    scores: np.ndarray,
    match_mask: np.ndarray,
    *,
    n_gallery_balanced: int,
    n_true_balanced: int,
    n_repeats: int,
    topk_values: list[int],
    rng: np.random.Generator,
) -> dict | None:
    """
    Evaluate one query under balanced candidate composition:
    fixed gallery size + fixed positive count via repeated subsampling.
    """
    pos_scores = scores[match_mask]
    neg_scores = scores[~match_mask]
    n_pos = int(len(pos_scores))
    n_neg = int(len(neg_scores))

    n_true = int(n_true_balanced)
    n_gallery = int(n_gallery_balanced)
    n_neg_needed = n_gallery - n_true
    if n_true <= 0 or n_gallery <= n_true:
        return None
    if n_pos < n_true or n_neg < n_neg_needed:
        return None

    ranks = np.empty(n_repeats, dtype=np.float64)
    topk_hits: dict[int, np.ndarray] = {k: np.empty(n_repeats, dtype=np.float64) for k in topk_values}
    for i in range(n_repeats):
        pos_idx = rng.choice(n_pos, size=n_true, replace=False)
        neg_idx = rng.choice(n_neg, size=n_neg_needed, replace=False)
        sampled_pos = pos_scores[pos_idx]
        sampled_neg = neg_scores[neg_idx]
        best_pos = float(np.max(sampled_pos))
        rank = 1 + int(np.sum(sampled_neg > best_pos))
        ranks[i] = float(rank)
        for k in topk_values:
            topk_hits[k][i] = float(rank <= k)

    out = {
        "balanced_valid": True,
        "balanced_n_gallery_eval": int(n_gallery),
        "balanced_n_true_eval": int(n_true),
        "true_rank_balanced_mean": float(np.mean(ranks)),
        "mrr_balanced": float(np.mean(1.0 / ranks)),
    }
    for k in topk_values:
        out[f"success_top{k}_balanced"] = float(np.mean(topk_hits[k]))
    return out


def run_task1_retrieval(cfg: Task1RetrievalConfig) -> None:
    t0 = time.time()
    rng = np.random.default_rng(cfg.random_seed)
    _ensure_dirs(cfg)

    tracks = [t.lower() for t in _parse_list_str(cfg.tracks, str)]
    topk_values = sorted(set(_parse_list_str(cfg.topk, int)))
    directions = _parse_directions(cfg.directions)

    candidates = _load_candidates(cfg)
    y_gene, y_path = _fetch_embeddings_for_rows(candidates, cfg=cfg)

    candidates.to_csv(Path(cfg.qc_dir) / "retrieval_candidates_snapshot.csv", index=False)
    pd.DataFrame(
        {
            "source_db": candidates["source_db"],
            "modality": candidates["modality"],
            "n_rows": 1,
        }
    ).groupby(["source_db", "modality"]).sum().reset_index().to_csv(
        Path(cfg.qc_dir) / "retrieval_input_counts.csv", index=False
    )

    emb_map = {"gene": y_gene, "path": y_path}
    per_query_rows: list[dict] = []

    for modality in sorted(candidates["modality"].unique()):
        sub_mod = candidates[candidates["modality"] == modality].copy()
        for q_src, g_src in directions:
            query_df = sub_mod[sub_mod["source_db"] == q_src].copy().reset_index(drop=True)
            gallery_df = sub_mod[sub_mod["source_db"] == g_src].copy().reset_index(drop=True)
            if len(query_df) == 0 or len(gallery_df) == 0:
                continue

            if cfg.max_queries_per_group is not None and len(query_df) > cfg.max_queries_per_group:
                pick = np.sort(rng.choice(len(query_df), size=cfg.max_queries_per_group, replace=False))
                query_df = query_df.iloc[pick].reset_index(drop=True)
            if cfg.max_gallery_per_group is not None and len(gallery_df) > cfg.max_gallery_per_group:
                pick = np.sort(rng.choice(len(gallery_df), size=cfg.max_gallery_per_group, replace=False))
                gallery_df = gallery_df.iloc[pick].reset_index(drop=True)

            q_idx = query_df["retrieval_row_idx"].to_numpy(dtype=int)
            g_idx = gallery_df["retrieval_row_idx"].to_numpy(dtype=int)
            q_labels = query_df["label_key"].to_numpy(dtype=str)
            g_labels = gallery_df["label_key"].to_numpy(dtype=str)

            for track in tracks:
                if track not in emb_map:
                    raise ValueError(f"Unknown track '{track}'. Supported: gene,path")

                q_emb = emb_map[track][q_idx]
                g_emb = emb_map[track][g_idx]

                for start in range(0, len(query_df), cfg.batch_size):
                    end = min(start + cfg.batch_size, len(query_df))
                    sim = _cosine_sim(q_emb[start:end], g_emb)
                    top1_idx = np.argmax(sim, axis=1)

                    for i in range(end - start):
                        q_row = query_df.iloc[start + i]
                        scores = sim[i]
                        label = q_labels[start + i]
                        match_mask = g_labels == label
                        n_match = int(np.sum(match_mask))
                        if n_match <= 0:
                            continue

                        best_match_score = float(np.max(scores[match_mask]))
                        true_rank = 1 + int(np.sum(scores > best_match_score))
                        top1_g_idx = int(top1_idx[i])
                        top1_label = str(g_labels[top1_g_idx])
                        top1_score = float(scores[top1_g_idx])
                        row = {
                            "track": track,
                            "direction": f"{q_src}->{g_src}",
                            "query_source": q_src,
                            "gallery_source": g_src,
                            "modality": str(modality),
                            "query_row_id": int(q_row["task1_row_id"]),
                            "query_cell_std": str(q_row["cell_std"]),
                            "query_target_std": str(q_row["target_std"]),
                            "query_label_key": label,
                            "query_dose_val": q_row["dose_val"],
                            "query_time_val": q_row["time_val"],
                            "n_gallery": int(len(gallery_df)),
                            "n_true_in_gallery": int(n_match),
                            "true_rank": int(true_rank),
                            "mrr": float(1.0 / true_rank),
                            "best_match_score": best_match_score,
                            "top1_label_key": top1_label,
                            "top1_score": top1_score,
                            "top1_is_true": bool(top1_label == label),
                        }
                        for k in topk_values:
                            row[f"success_top{k}"] = bool(true_rank <= k)

                        if cfg.run_balanced_eval:
                            balanced = _balanced_query_metrics(
                                scores=scores,
                                match_mask=match_mask,
                                n_gallery_balanced=cfg.balanced_gallery_size,
                                n_true_balanced=cfg.balanced_true_per_query,
                                n_repeats=cfg.balanced_n_repeats,
                                topk_values=topk_values,
                                rng=rng,
                            )
                            if balanced is None:
                                row["balanced_valid"] = False
                                row["balanced_n_gallery_eval"] = cfg.balanced_gallery_size
                                row["balanced_n_true_eval"] = cfg.balanced_true_per_query
                                row["true_rank_balanced_mean"] = np.nan
                                row["mrr_balanced"] = np.nan
                                for k in topk_values:
                                    row[f"success_top{k}_balanced"] = np.nan
                            else:
                                row.update(balanced)
                        per_query_rows.append(row)

    per_query = pd.DataFrame(per_query_rows)
    if per_query.empty:
        raise RuntimeError("No retrieval per-query results were produced.")

    per_query.to_csv(Path(cfg.analysis_dir) / "retrieval_per_query.csv", index=False)

    summary_rows: list[dict] = []
    null_rows: list[dict] = []
    group_cols = ["track", "direction", "modality"]
    for keys, grp in per_query.groupby(group_cols, sort=False):
        track, direction, modality = keys
        n_query = len(grp)
        if n_query < cfg.min_queries_per_group:
            continue
        balanced_valid = None
        if cfg.run_balanced_eval and "balanced_valid" in grp.columns:
            balanced_valid = grp["balanced_valid"].astype(bool)
        n_balanced = int(balanced_valid.sum()) if balanced_valid is not None else 0
        row = {
            "track": track,
            "direction": direction,
            "modality": modality,
            "n_query": int(n_query),
            "n_gallery": int(grp["n_gallery"].iloc[0]),
            "mean_true_rank": float(grp["true_rank"].mean()),
            "median_true_rank": float(grp["true_rank"].median()),
            "mrr": float(grp["mrr"].mean()),
            "top1": float(grp["success_top1"].mean()) if "success_top1" in grp.columns else float(np.mean(grp["true_rank"] <= 1)),
            "n_query_balanced": n_balanced,
        }
        for k in topk_values:
            col = f"success_top{k}"
            row[f"top{k}"] = float(grp[col].mean()) if col in grp.columns else float(np.mean(grp["true_rank"] <= k))

        if cfg.run_balanced_eval and balanced_valid is not None and n_balanced > 0:
            balanced_grp = grp[balanced_valid].copy()
            row["mrr_balanced"] = float(balanced_grp["mrr_balanced"].mean())
            row["top1_balanced"] = float(balanced_grp["success_top1_balanced"].mean())
            for k in topk_values:
                col = f"success_top{k}_balanced"
                if col in balanced_grp.columns:
                    row[f"top{k}_balanced"] = float(balanced_grp[col].mean())
        else:
            row["mrr_balanced"] = np.nan
            row["top1_balanced"] = np.nan
            for k in topk_values:
                row[f"top{k}_balanced"] = np.nan
        summary_rows.append(row)

        n_match_arr = grp["n_true_in_gallery"].to_numpy(dtype=int)
        n_gallery = int(grp["n_gallery"].iloc[0])
        null_stat = _simulate_null_group(
            n_match_arr=n_match_arr,
            n_gallery=n_gallery,
            n_perm=cfg.n_perm,
            rng=rng,
            topk_values=topk_values,
        )
        null_row = {
            "track": track,
            "direction": direction,
            "modality": modality,
            "n_query": int(n_query),
            "n_gallery": n_gallery,
            "obs_mrr": float(row["mrr"]),
            "obs_top1": float(row["top1"]),
            "p_mrr_ge_null": float((np.sum(null_stat["_null_arrays"]["mrr"] >= row["mrr"]) + 1) / (cfg.n_perm + 1)),
            "p_top1_ge_null": float((np.sum(null_stat["_null_arrays"]["top1"] >= row["top1"]) + 1) / (cfg.n_perm + 1)),
            "null_mrr_mean": null_stat["null_mrr_mean"],
            "null_mrr_std": null_stat["null_mrr_std"],
            "null_top1_mean": null_stat["null_top1_mean"],
            "null_top1_std": null_stat["null_top1_std"],
            "null_mrr_q05": null_stat["null_mrr_q05"],
            "null_mrr_q95": null_stat["null_mrr_q95"],
            "null_top1_q05": null_stat["null_top1_q05"],
            "null_top1_q95": null_stat["null_top1_q95"],
        }
        for k in topk_values:
            obs = float(row[f"top{k}"])
            null_arr = null_stat["_null_arrays"][f"top{k}"]
            null_row[f"obs_top{k}"] = obs
            null_row[f"p_top{k}_ge_null"] = float((np.sum(null_arr >= obs) + 1) / (cfg.n_perm + 1))
            null_row[f"null_top{k}_mean"] = null_stat[f"null_top{k}_mean"]
            null_row[f"null_top{k}_std"] = null_stat[f"null_top{k}_std"]

        if cfg.run_balanced_eval and n_balanced >= cfg.min_queries_per_group:
            n_true = int(cfg.balanced_true_per_query)
            n_gallery_bal = int(cfg.balanced_gallery_size)
            balanced_null = _simulate_null_group(
                n_match_arr=np.full(n_balanced, n_true, dtype=int),
                n_gallery=n_gallery_bal,
                n_perm=cfg.n_perm,
                rng=rng,
                topk_values=topk_values,
            )
            null_row["obs_mrr_balanced"] = float(row["mrr_balanced"])
            null_row["obs_top1_balanced"] = float(row["top1_balanced"])
            null_row["p_mrr_balanced_ge_null"] = float(
                (np.sum(balanced_null["_null_arrays"]["mrr"] >= row["mrr_balanced"]) + 1) / (cfg.n_perm + 1)
            )
            null_row["p_top1_balanced_ge_null"] = float(
                (np.sum(balanced_null["_null_arrays"]["top1"] >= row["top1_balanced"]) + 1) / (cfg.n_perm + 1)
            )
            null_row["null_mrr_balanced_mean"] = balanced_null["null_mrr_mean"]
            null_row["null_mrr_balanced_std"] = balanced_null["null_mrr_std"]
            null_row["null_top1_balanced_mean"] = balanced_null["null_top1_mean"]
            null_row["null_top1_balanced_std"] = balanced_null["null_top1_std"]
            for k in topk_values:
                obs = float(row.get(f"top{k}_balanced", np.nan))
                null_arr = balanced_null["_null_arrays"][f"top{k}"]
                null_row[f"obs_top{k}_balanced"] = obs
                null_row[f"p_top{k}_balanced_ge_null"] = float((np.sum(null_arr >= obs) + 1) / (cfg.n_perm + 1))
                null_row[f"null_top{k}_balanced_mean"] = balanced_null[f"null_top{k}_mean"]
                null_row[f"null_top{k}_balanced_std"] = balanced_null[f"null_top{k}_std"]
        null_rows.append(null_row)

    summary = pd.DataFrame(summary_rows)
    null_summary = pd.DataFrame(null_rows)
    summary.to_csv(Path(cfg.analysis_dir) / "retrieval_summary.csv", index=False)
    null_summary.to_csv(Path(cfg.analysis_dir) / "retrieval_null_summary.csv", index=False)

    # helpful examples
    examples = []
    for keys, grp in per_query.groupby(group_cols, sort=False):
        track, direction, modality = keys
        best = grp.nsmallest(3, "true_rank")
        worst = grp.nlargest(3, "true_rank")
        for tag, sub in [("best", best), ("worst", worst)]:
            tmp = sub.copy()
            tmp["example_type"] = tag
            tmp["track"] = track
            tmp["direction"] = direction
            tmp["modality"] = modality
            examples.append(tmp)
    if examples:
        example_df = pd.concat(examples, axis=0, ignore_index=True)
        example_df.to_csv(Path(cfg.analysis_dir) / "retrieval_examples.csv", index=False)

    run_manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": {
            "n_candidates": int(len(candidates)),
            "n_per_query_rows": int(len(per_query)),
            "n_summary_groups": int(len(summary)),
            "balanced_eval_enabled": bool(cfg.run_balanced_eval),
        },
        "outputs": {
            "analysis_dir": cfg.analysis_dir,
            "qc_dir": cfg.qc_dir,
        },
    }
    _atomic_json_save(run_manifest, cfg.run_manifest_path)


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("M2M-Bench Task1 instance-level retrieval analysis")
    parser.add_argument("--m1-candidates-path", type=str, default=Task1RetrievalConfig.m1_candidates_path)
    parser.add_argument("--processed-dir", type=str, default=Task1RetrievalConfig.processed_dir)
    parser.add_argument("--lincs-tensor-name", type=str, default=Task1RetrievalConfig.lincs_tensor_name)
    parser.add_argument("--output-dir", type=str, default=Task1RetrievalConfig.output_dir)
    parser.add_argument("--tracks", type=str, default=Task1RetrievalConfig.tracks)
    parser.add_argument("--directions", type=str, default=Task1RetrievalConfig.directions)
    parser.add_argument("--topk", type=str, default=Task1RetrievalConfig.topk)
    parser.add_argument("--batch-size", type=int, default=Task1RetrievalConfig.batch_size)
    parser.add_argument("--n-perm", type=int, default=Task1RetrievalConfig.n_perm)
    parser.add_argument("--min-queries-per-group", type=int, default=Task1RetrievalConfig.min_queries_per_group)
    parser.add_argument("--random-seed", type=int, default=Task1RetrievalConfig.random_seed)
    parser.add_argument("--max-queries-per-group", type=int, default=Task1RetrievalConfig.max_queries_per_group)
    parser.add_argument("--max-gallery-per-group", type=int, default=Task1RetrievalConfig.max_gallery_per_group)
    parser.add_argument(
        "--run-balanced-eval",
        action=argparse.BooleanOptionalAction,
        default=Task1RetrievalConfig.run_balanced_eval,
        help="Compute balanced candidate-space retrieval metrics.",
    )
    parser.add_argument(
        "--balanced-gallery-size",
        type=int,
        default=Task1RetrievalConfig.balanced_gallery_size,
        help="Balanced evaluation gallery size.",
    )
    parser.add_argument(
        "--balanced-true-per-query",
        type=int,
        default=Task1RetrievalConfig.balanced_true_per_query,
        help="Balanced evaluation positives per query.",
    )
    parser.add_argument(
        "--balanced-n-repeats",
        type=int,
        default=Task1RetrievalConfig.balanced_n_repeats,
        help="Repeated subsampling iterations for balanced evaluation.",
    )
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    cfg = Task1RetrievalConfig(
        m1_candidates_path=args.m1_candidates_path,
        processed_dir=args.processed_dir,
        lincs_tensor_name=args.lincs_tensor_name,
        output_dir=args.output_dir,
        tracks=args.tracks,
        directions=args.directions,
        topk=args.topk,
        batch_size=args.batch_size,
        n_perm=args.n_perm,
        min_queries_per_group=args.min_queries_per_group,
        random_seed=args.random_seed,
        max_queries_per_group=args.max_queries_per_group,
        max_gallery_per_group=args.max_gallery_per_group,
        run_balanced_eval=args.run_balanced_eval,
        balanced_gallery_size=args.balanced_gallery_size,
        balanced_true_per_query=args.balanced_true_per_query,
        balanced_n_repeats=args.balanced_n_repeats,
    )
    run_task1_retrieval(cfg)


if __name__ == "__main__":
    main()
