#!/usr/bin/env python3
"""
Task2 target-confidence and polypharmacology stratification analysis.

This script extends existing Task2 outputs with an explicit target-tier layer:
1) assign each target into confidence/polypharmacology tiers,
2) summarize performance classes by tier,
3) quantify enrichment with Fisher exact + BH-FDR correction.
"""

from __future__ import annotations

import argparse
import json
import re
import time
from dataclasses import asdict, dataclass
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.stats import fisher_exact


# ---------------------------------------------------------------------
# Section 0. Configuration
# ---------------------------------------------------------------------


@dataclass
class Task2TargetTierConfig:
    context_labels_path: str = "./outputs/task2_nodomain/analysis/Step2_Context_Labels_NoDomain.csv"
    output_dir: str = "./outputs/task2_nodomain/target_tier"
    target_tier_map_path: str | None = None
    min_contexts_per_tier: int = 10

    @property
    def analysis_dir(self) -> str:
        return str(Path(self.output_dir) / "analysis")

    @property
    def run_manifest_path(self) -> str:
        return str(Path(self.output_dir) / "run_manifest_task2_target_tier.json")


# ---------------------------------------------------------------------
# Section 1. Helpers
# ---------------------------------------------------------------------


def _ensure_dir(path_like: str) -> None:
    Path(path_like).mkdir(parents=True, exist_ok=True)


def _bh_fdr(p_values: np.ndarray) -> np.ndarray:
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


def _stable_log2_or(a: int, b: int, c: int, d: int) -> float:
    num = (a + 0.5) * (d + 0.5)
    den = (b + 0.5) * (c + 0.5)
    return float(np.log2(num / den))


# ---------------------------------------------------------------------
# Section 2. Target tier assignment
# ---------------------------------------------------------------------


HIGH_CONFIDENCE_TARGETS = {
    "BRAF",
    "EGFR",
    "MAP2K1",
    "MTOR",
    "PIK3CA",
    "PSMB1",
}

POLYPHARM_OR_AMBIG_TARGETS = {
    "HSP90AA1",
    "NR3C1",
    "PYGM",
}


def _normalize_target(x: str) -> str:
    return str(x).strip().upper()


def _heuristic_target_tier(target: str) -> tuple[str, str, str]:
    """
    Returns (confidence_tier, polypharm_tier, reason).
    """
    t = _normalize_target(target)
    if t in HIGH_CONFIDENCE_TARGETS:
        return "Tier_A_HighConfidence", "LowPolyRisk", "heuristic_high_confidence_list"
    if t in POLYPHARM_OR_AMBIG_TARGETS:
        return "Tier_C_AmbiguousOrPolypharm", "HighPolyRisk", "heuristic_polypharm_list"
    if re.search(r"[\|,;/]|\bAND\b", t):
        return "Tier_C_AmbiguousOrPolypharm", "HighPolyRisk", "multi_target_pattern"
    if t.startswith("SG"):
        return "Tier_B_ModerateOrUnknown", "UnknownPolyRisk", "sgRNA_label_target"
    return "Tier_B_ModerateOrUnknown", "UnknownPolyRisk", "default_fallback"


def _load_external_tier_map(path: str) -> pd.DataFrame:
    df = pd.read_csv(path)
    cols = {c.lower(): c for c in df.columns}
    target_col = cols.get("target") or cols.get("target_std")
    conf_col = cols.get("confidence_tier")
    poly_col = cols.get("polypharm_tier")
    if target_col is None or conf_col is None or poly_col is None:
        raise ValueError(
            "target_tier_map must contain columns: Target(or target_std), confidence_tier, polypharm_tier"
        )
    out = df[[target_col, conf_col, poly_col]].copy()
    out.columns = ["Target", "confidence_tier", "polypharm_tier"]
    out["Target"] = out["Target"].map(_normalize_target)
    out["mapping_source"] = "external_csv"
    return out.drop_duplicates(subset=["Target"], keep="first")


def assign_target_tiers(ctx: pd.DataFrame, cfg: Task2TargetTierConfig) -> tuple[pd.DataFrame, pd.DataFrame]:
    ctx = ctx.copy()
    ctx["Target_norm"] = ctx["Target"].map(_normalize_target)

    ext_map = None
    if cfg.target_tier_map_path:
        ext_map = _load_external_tier_map(cfg.target_tier_map_path)

    # Start with heuristics
    heur_rows = []
    for target in sorted(ctx["Target_norm"].dropna().unique().tolist()):
        conf, poly, reason = _heuristic_target_tier(target)
        heur_rows.append(
            {
                "Target": target,
                "confidence_tier": conf,
                "polypharm_tier": poly,
                "mapping_source": reason,
            }
        )
    mapping = pd.DataFrame(heur_rows)

    # External map overrides heuristics when provided
    if ext_map is not None and not ext_map.empty:
        mapping = mapping.set_index("Target")
        ext_map = ext_map.set_index("Target")
        mapping.update(ext_map[["confidence_tier", "polypharm_tier", "mapping_source"]])
        mapping = mapping.reset_index()

    out = ctx.merge(mapping, left_on="Target_norm", right_on="Target", how="left")
    out = out.drop(columns=["Target_norm", "Target_y"]).rename(columns={"Target_x": "Target"})
    return out, mapping


# ---------------------------------------------------------------------
# Section 3. Tier-level summaries and enrichment
# ---------------------------------------------------------------------


def build_tier_summary(ctx_tier: pd.DataFrame, min_contexts: int) -> pd.DataFrame:
    summary = (
        ctx_tier.groupby(["confidence_tier", "polypharm_tier"], as_index=False)
        .agg(
            n_contexts=("Target", "size"),
            n_targets=("Target", "nunique"),
            mean_success=("Mean_Success", "mean"),
            median_success=("Mean_Success", "median"),
            mean_peak_success=("Peak_Success", "mean"),
        )
        .sort_values(["n_contexts"], ascending=False)
    )

    class_mix = (
        ctx_tier.groupby(["confidence_tier", "polypharm_tier", "Performance_Class"], as_index=False)
        .size()
        .rename(columns={"size": "n_contexts_class"})
    )
    summary = summary.merge(class_mix, on=["confidence_tier", "polypharm_tier"], how="left")
    summary["class_fraction"] = summary["n_contexts_class"] / summary["n_contexts"]
    summary["meets_min_contexts"] = summary["n_contexts"] >= int(min_contexts)
    return summary


def build_tier_enrichment(ctx_tier: pd.DataFrame, min_contexts: int) -> pd.DataFrame:
    rows: list[dict] = []
    classes = sorted(ctx_tier["Performance_Class"].dropna().unique().tolist())
    tiers = sorted(ctx_tier["confidence_tier"].dropna().unique().tolist())
    for cls in classes:
        y = (ctx_tier["Performance_Class"] == cls).to_numpy()
        for tier in tiers:
            x = (ctx_tier["confidence_tier"] == tier).to_numpy()
            a = int(np.sum(x & y))
            b = int(np.sum(x & (~y)))
            c = int(np.sum((~x) & y))
            d = int(np.sum((~x) & (~y)))
            tier_n = int(a + b)
            if tier_n < int(min_contexts):
                continue
            _, pval = fisher_exact([[a, b], [c, d]], alternative="two-sided")
            rows.append(
                {
                    "Performance_Class": cls,
                    "confidence_tier": tier,
                    "tier_contexts": tier_n,
                    "class_in_tier": a,
                    "class_not_in_tier": b,
                    "class_in_other_tiers": c,
                    "class_not_in_other_tiers": d,
                    "fraction_in_tier": float(a / tier_n) if tier_n > 0 else np.nan,
                    "fraction_in_other_tiers": float(c / (c + d)) if (c + d) > 0 else np.nan,
                    "log2_odds_ratio": _stable_log2_or(a, b, c, d),
                    "p_value": float(pval),
                }
            )
    out = pd.DataFrame(rows)
    if not out.empty:
        out["fdr_bh"] = _bh_fdr(out["p_value"].to_numpy())
    return out.sort_values(["Performance_Class", "fdr_bh", "p_value"])


# ---------------------------------------------------------------------
# Section 4. Main
# ---------------------------------------------------------------------


def run_task2_target_tier(cfg: Task2TargetTierConfig) -> None:
    t0 = time.time()
    _ensure_dir(cfg.analysis_dir)

    ctx = pd.read_csv(cfg.context_labels_path)
    need = ["Cell", "Target", "Mean_Success", "Peak_Success", "N_Instances", "Performance_Class"]
    missing = [c for c in need if c not in ctx.columns]
    if missing:
        raise ValueError(f"Context label table missing required columns: {missing}")

    ctx_tier, mapping = assign_target_tiers(ctx=ctx, cfg=cfg)
    summary = build_tier_summary(ctx_tier=ctx_tier, min_contexts=cfg.min_contexts_per_tier)
    enrich = build_tier_enrichment(ctx_tier=ctx_tier, min_contexts=cfg.min_contexts_per_tier)

    out_dir = Path(cfg.analysis_dir)
    mapping.to_csv(out_dir / "task2_target_tier_map_used.csv", index=False)
    ctx_tier.to_csv(out_dir / "task2_context_with_target_tier.csv", index=False)
    summary.to_csv(out_dir / "task2_target_tier_summary.csv", index=False)
    enrich.to_csv(out_dir / "task2_target_tier_enrichment.csv", index=False)

    key_summary = pd.DataFrame(
        [
            {"metric": "n_context_rows", "value": int(len(ctx_tier))},
            {"metric": "n_unique_targets", "value": int(ctx_tier["Target"].nunique())},
            {"metric": "n_confidence_tiers", "value": int(ctx_tier["confidence_tier"].nunique())},
            {"metric": "n_polypharm_tiers", "value": int(ctx_tier["polypharm_tier"].nunique())},
            {"metric": "n_enrichment_tests", "value": int(len(enrich))},
        ]
    )
    key_summary.to_csv(out_dir / "task2_target_tier_key_summary.csv", index=False)

    manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": {
            "n_context_rows": int(len(ctx_tier)),
            "n_mapping_rows": int(len(mapping)),
            "n_summary_rows": int(len(summary)),
            "n_enrichment_rows": int(len(enrich)),
        },
    }
    Path(cfg.run_manifest_path).write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    print(f"[Task2 Target Tier] done. Outputs: {cfg.analysis_dir}")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("Task2 target tier / polypharmacology analysis")
    parser.add_argument("--context-labels-path", type=str, default=Task2TargetTierConfig.context_labels_path)
    parser.add_argument("--output-dir", type=str, default=Task2TargetTierConfig.output_dir)
    parser.add_argument("--target-tier-map-path", type=str, default=Task2TargetTierConfig.target_tier_map_path)
    parser.add_argument("--min-contexts-per-tier", type=int, default=Task2TargetTierConfig.min_contexts_per_tier)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    cfg = Task2TargetTierConfig(
        context_labels_path=args.context_labels_path,
        output_dir=args.output_dir,
        target_tier_map_path=args.target_tier_map_path,
        min_contexts_per_tier=args.min_contexts_per_tier,
    )
    run_task2_target_tier(cfg)


if __name__ == "__main__":
    main()
