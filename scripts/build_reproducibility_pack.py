#!/usr/bin/env python3
"""
Build reproducibility pack for M2M-Bench.

Outputs:
1) manifest inventory (what runs exist),
2) seed registry extracted from manifests,
3) environment snapshot (python/platform/git/core package versions),
4) one-command reproducibility path markdown.
"""

from __future__ import annotations

import argparse
import json
import platform
import subprocess
import sys
import time
from dataclasses import asdict, dataclass
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path
from typing import Any

import pandas as pd


@dataclass
class ReproPackConfig:
    output_dir: str = "./outputs/reproducibility"
    manifest_paths: list[str] | None = None

    @property
    def analysis_dir(self) -> str:
        return str(Path(self.output_dir) / "analysis")

    @property
    def run_manifest_path(self) -> str:
        return str(Path(self.output_dir) / "run_manifest_reproducibility_pack.json")


DEFAULT_MANIFEST_PATHS = [
    "./outputs/task0_curated/run_manifest_task0.json",
    "./outputs/task1/run_manifest_task1_pipeline.json",
    "./outputs/task1_audit/run_manifest_task01_attrition.json",
    "./outputs/task1_reviewer_fixes/run_manifest_task1_reviewer_addons.json",
    "./outputs/task1_reviewer_fixes/retrieval_sensitivity/run_manifest_task1_retrieval_sensitivity.json",
    "./outputs/task2_nodomain/run_manifest_task2_nodomain.json",
    "./outputs/task2_nodomain/target_tier/run_manifest_task2_target_tier.json",
    "./outputs/task3_audit/run_manifest_task3_audit.json",
    "./outputs/task3_meta/run_manifest_task3_meta.json",
    "./outputs/project_alignment_audit/run_manifest_project_alignment_audit.json",
]


def _safe_json_load(path: Path) -> dict[str, Any] | None:
    if not path.exists():
        return None
    try:
        return json.loads(path.read_text(encoding="utf-8"))
    except Exception:
        return None


def _extract_seed_rows(obj: Any, manifest_path: str, prefix: str = "") -> list[dict]:
    rows: list[dict] = []
    if isinstance(obj, dict):
        for k, v in obj.items():
            kp = f"{prefix}.{k}" if prefix else str(k)
            if "seed" in str(k).lower():
                rows.append(
                    {
                        "manifest_path": manifest_path,
                        "key_path": kp,
                        "value": v,
                    }
                )
            rows.extend(_extract_seed_rows(v, manifest_path=manifest_path, prefix=kp))
    elif isinstance(obj, list):
        for i, v in enumerate(obj):
            kp = f"{prefix}[{i}]"
            rows.extend(_extract_seed_rows(v, manifest_path=manifest_path, prefix=kp))
    return rows


def _pkg_version(name: str) -> str:
    try:
        return version(name)
    except PackageNotFoundError:
        return "not_installed"


def _git_info() -> dict[str, str]:
    info = {"git_commit": "unknown", "git_branch": "unknown", "git_status_short": "unknown"}
    try:
        info["git_commit"] = (
            subprocess.check_output(["git", "rev-parse", "HEAD"], text=True, stderr=subprocess.DEVNULL).strip()
        )
        info["git_branch"] = (
            subprocess.check_output(["git", "rev-parse", "--abbrev-ref", "HEAD"], text=True, stderr=subprocess.DEVNULL).strip()
        )
        info["git_status_short"] = (
            subprocess.check_output(["git", "status", "--short"], text=True, stderr=subprocess.DEVNULL).strip()
        )
    except Exception:
        pass
    return info


def _write_one_command_markdown(path: Path) -> None:
    text = """# M2M-Bench one-command reproducibility paths

> Note: data download/provisioning is external; commands assume input files already available.

## Task1 reviewer-critical outputs

```bash
python scripts/task1_reviewer_addons.py \\
  --m1-candidates-path ./outputs/task1/data/m1_candidates.csv \\
  --per-pair-path ./outputs/task1/analysis/modality_gap_per_pair.csv \\
  --retrieval-summary-path ./outputs/task1/retrieval/analysis/retrieval_summary.csv \\
  --retrieval-null-path ./outputs/task1/retrieval/analysis/retrieval_null_summary.csv \\
  --retrieval-per-query-path ./outputs/task1/retrieval/analysis/retrieval_per_query.csv \\
  --context-overlap-path ./outputs/task1_audit/analysis/context_overlap_counts.csv \\
  --processed-dir /mnt/NAS_21T/ProjectData/OSMOSIS/processed \\
  --output-dir ./outputs/task1_reviewer_fixes
```

## Task2 target-tier sensitivity layer

```bash
python scripts/task2_target_tier_analysis.py \\
  --context-labels-path ./outputs/task2_nodomain/analysis/Step2_Context_Labels_NoDomain.csv \\
  --output-dir ./outputs/task2_nodomain/target_tier
```

## Task3 model meta-analysis

```bash
python scripts/task3_fm_meta_analysis.py \\
  --scoreboard-long-path /mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Task3_Unified/viz_scoreboard_long.csv \\
  --output-dir ./outputs/task3_meta
```

## Build reproducibility pack (this script)

```bash
python scripts/build_reproducibility_pack.py --output-dir ./outputs/reproducibility
```
"""
    path.write_text(text, encoding="utf-8")


def run_repro_pack(cfg: ReproPackConfig) -> None:
    t0 = time.time()
    manifests = cfg.manifest_paths or DEFAULT_MANIFEST_PATHS
    out = Path(cfg.analysis_dir)
    out.mkdir(parents=True, exist_ok=True)

    inv_rows: list[dict] = []
    seed_rows: list[dict] = []
    for mp in manifests:
        p = Path(mp)
        obj = _safe_json_load(p)
        exists = obj is not None
        inv_rows.append(
            {
                "manifest_path": mp,
                "exists": bool(exists),
                "created_at": obj.get("created_at") if exists else None,
                "elapsed_seconds": obj.get("elapsed_seconds") if exists else None,
                "summary_json": json.dumps(obj.get("summary", {}), ensure_ascii=False) if exists else None,
            }
        )
        if exists:
            seed_rows.extend(_extract_seed_rows(obj, manifest_path=mp))

    inv = pd.DataFrame(inv_rows)
    seed = pd.DataFrame(seed_rows)
    inv.to_csv(out / "repro_manifest_inventory.csv", index=False)
    seed.to_csv(out / "repro_seed_registry.csv", index=False)

    env = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "python_version": sys.version,
        "platform": platform.platform(),
        "packages": {
            "numpy": _pkg_version("numpy"),
            "pandas": _pkg_version("pandas"),
            "scipy": _pkg_version("scipy"),
            "torch": _pkg_version("torch"),
            "matplotlib": _pkg_version("matplotlib"),
        },
        **_git_info(),
    }
    (out / "repro_environment.json").write_text(json.dumps(env, indent=2), encoding="utf-8")
    _write_one_command_markdown(out / "repro_one_command_paths.md")

    run_manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": {
            "n_manifest_candidates": int(len(manifests)),
            "n_manifest_found": int(inv["exists"].sum()) if not inv.empty else 0,
            "n_seed_entries": int(len(seed)),
        },
    }
    Path(cfg.run_manifest_path).write_text(json.dumps(run_manifest, indent=2), encoding="utf-8")
    print(f"[Repro Pack] done. Outputs: {cfg.analysis_dir}")


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser("Build reproducibility pack for M2M-Bench")
    p.add_argument("--output-dir", type=str, default=ReproPackConfig.output_dir)
    p.add_argument(
        "--manifest-path",
        action="append",
        default=None,
        help="Optional manifest path. Repeat this flag to provide multiple manifests.",
    )
    return p


def main() -> None:
    args = _build_parser().parse_args()
    cfg = ReproPackConfig(
        output_dir=args.output_dir,
        manifest_paths=args.manifest_path if args.manifest_path else None,
    )
    run_repro_pack(cfg)


if __name__ == "__main__":
    main()
