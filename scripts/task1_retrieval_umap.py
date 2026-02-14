#!/usr/bin/env python3

from __future__ import annotations

import argparse
from pathlib import Path

import numpy as np
import pandas as pd

from _script_bootstrap import setup_project_imports

setup_project_imports()

from m2m_bench.task1.retrieval_instance import (
    Task1RetrievalConfig,
    _fetch_embeddings_for_rows,
    _load_candidates,
)


def _make_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("Task1 retrieval embedding visualization (UMAP/PCA)")
    parser.add_argument(
        "--m1-candidates-path",
        type=str,
        default="./outputs/task1/retrieval/qc/retrieval_candidates_snapshot.csv",
        help="Candidates table used for retrieval (csv/parquet).",
    )
    parser.add_argument(
        "--processed-dir",
        type=str,
        default="/mnt/NAS_21T/ProjectData/OSMOSIS/processed",
    )
    parser.add_argument(
        "--lincs-tensor-name",
        type=str,
        default="LINCS_Engine1_TrainData.pt",
    )
    parser.add_argument(
        "--output-dir",
        type=str,
        default="./outputs/task1/retrieval/qc",
    )
    parser.add_argument(
        "--max-points",
        type=int,
        default=8000,
        help="Max points to plot per track.",
    )
    parser.add_argument(
        "--random-seed",
        type=int,
        default=42,
    )
    parser.add_argument(
        "--point-size",
        type=float,
        default=5.0,
    )
    return parser


def _get_2d_projection(emb: np.ndarray, random_seed: int) -> tuple[str, np.ndarray]:
    try:
        import umap.umap_ as umap

        reducer = umap.UMAP(
            n_neighbors=30,
            min_dist=0.15,
            metric="cosine",
            random_state=random_seed,
        )
        z = reducer.fit_transform(emb)
        return "UMAP", z
    except Exception:
        from sklearn.decomposition import PCA

        z = PCA(n_components=2, random_state=random_seed).fit_transform(emb)
        return "PCA", z


def _plot_scatter(
    df: pd.DataFrame,
    xcol: str,
    ycol: str,
    color_col: str,
    title: str,
    output_path: Path,
    point_size: float,
) -> None:
    try:
        import matplotlib.pyplot as plt
    except Exception as exc:
        raise RuntimeError("matplotlib is required to draw UMAP/PCA plots.") from exc

    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig, ax = plt.subplots(figsize=(7, 6))

    values = df[color_col].astype(str)
    uniq = sorted(values.unique().tolist())
    cmap = plt.get_cmap("tab20")
    for i, v in enumerate(uniq):
        sub = df[values == v]
        ax.scatter(
            sub[xcol].to_numpy(),
            sub[ycol].to_numpy(),
            s=point_size,
            alpha=0.6,
            label=v,
            color=cmap(i % 20),
            edgecolors="none",
        )

    ax.set_title(title)
    ax.set_xlabel(xcol)
    ax.set_ylabel(ycol)
    if len(uniq) <= 12:
        ax.legend(loc="best", fontsize=8, frameon=False)
    fig.tight_layout()
    fig.savefig(output_path, dpi=170)
    plt.close(fig)


def main() -> None:
    args = _make_parser().parse_args()
    rng = np.random.default_rng(args.random_seed)

    cfg = Task1RetrievalConfig(
        m1_candidates_path=args.m1_candidates_path,
        processed_dir=args.processed_dir,
        lincs_tensor_name=args.lincs_tensor_name,
        output_dir="./tmp_not_used",
    )

    candidates = _load_candidates(cfg)
    if len(candidates) > args.max_points:
        keep = np.sort(rng.choice(len(candidates), size=args.max_points, replace=False))
        candidates = candidates.iloc[keep].reset_index(drop=True)
        candidates["retrieval_row_idx"] = np.arange(len(candidates), dtype=np.int64)

    y_gene, y_path = _fetch_embeddings_for_rows(candidates, cfg)

    out_dir = Path(args.output_dir)
    out_dir.mkdir(parents=True, exist_ok=True)

    manifest_rows = []
    for track, emb in [("gene", y_gene), ("path", y_path)]:
        method, z = _get_2d_projection(emb, random_seed=args.random_seed)
        plot_df = candidates[["source_db", "modality", "cell_std", "target_std", "label_key"]].copy()
        plot_df[f"{track}_x"] = z[:, 0]
        plot_df[f"{track}_y"] = z[:, 1]

        plot_df.to_csv(out_dir / f"retrieval_{track}_2d_points.csv", index=False)

        _plot_scatter(
            plot_df,
            xcol=f"{track}_x",
            ycol=f"{track}_y",
            color_col="source_db",
            title=f"Task1 Retrieval {track} space ({method}) by source_db",
            output_path=out_dir / f"retrieval_{track}_{method.lower()}_by_source.png",
            point_size=args.point_size,
        )
        _plot_scatter(
            plot_df,
            xcol=f"{track}_x",
            ycol=f"{track}_y",
            color_col="modality",
            title=f"Task1 Retrieval {track} space ({method}) by modality",
            output_path=out_dir / f"retrieval_{track}_{method.lower()}_by_modality.png",
            point_size=args.point_size,
        )

        top_labels = (
            plot_df["label_key"].value_counts().head(8).index.astype(str).tolist()
        )
        sub = plot_df[plot_df["label_key"].isin(top_labels)].copy()
        if len(sub) > 0:
            _plot_scatter(
                sub,
                xcol=f"{track}_x",
                ycol=f"{track}_y",
                color_col="label_key",
                title=f"Task1 Retrieval {track} space ({method}) top labels",
                output_path=out_dir / f"retrieval_{track}_{method.lower()}_top_labels.png",
                point_size=max(args.point_size, 8.0),
            )

        manifest_rows.append(
            {
                "track": track,
                "projection_method": method,
                "n_points": int(len(plot_df)),
                "output_points_csv": str(out_dir / f"retrieval_{track}_2d_points.csv"),
            }
        )

    pd.DataFrame(manifest_rows).to_csv(out_dir / "retrieval_2d_manifest.csv", index=False)


if __name__ == "__main__":
    main()
