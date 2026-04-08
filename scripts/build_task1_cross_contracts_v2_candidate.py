from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Mapping, Sequence

import pandas as pd
import yaml

try:
    from m2mbench.utils.raw_rebuild import standardize_cell_line, standardize_token
except ModuleNotFoundError:
    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))
    from m2mbench.utils.raw_rebuild import standardize_cell_line, standardize_token

try:
    from path_policy import DEFAULT_TASK1_SNAPSHOT_V2_CANDIDATE_ROOT
except ModuleNotFoundError:
    from scripts.path_policy import DEFAULT_TASK1_SNAPSHOT_V2_CANDIDATE_ROOT

CONFIG_PATH = Path("config/config.yaml")


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(description="Build Task1 cross contracts v2 candidate")
    parser.add_argument("--project-root", type=Path, default=Path("."))
    parser.add_argument("--snapshot-root", type=Path, default=DEFAULT_TASK1_SNAPSHOT_V2_CANDIDATE_ROOT)
    return parser.parse_args()


def utc_now_iso() -> str:
    return datetime.now(timezone.utc).replace(microsecond=0).isoformat()


def resolve_config_path(project_root: Path, raw_path: str | Path) -> Path:
    path = Path(str(raw_path))
    if path.is_absolute():
        return path.resolve()
    return (project_root / path).resolve()


def write_json(path: Path, payload: Mapping[str, object]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, ensure_ascii=True) + "\n", encoding="utf-8")


def write_csv(frame: pd.DataFrame, path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    frame.to_csv(path, index=False, quoting=csv.QUOTE_MINIMAL)


def explode_membership(
    frame: pd.DataFrame,
    *,
    subtype: str,
    row_id_column: str,
) -> pd.DataFrame:
    out = frame.copy()
    out["target_token"] = out["target_tokens"].fillna("").astype(str).map(
        lambda text: tuple(token for token in text.split(";") if token)
    )
    if subtype == "CRISPR":
        out = out.loc[out["target_token"].map(len).eq(1)].copy()
    out = out.explode("target_token")
    out["target_token"] = out["target_token"].fillna("").astype(str).map(standardize_token)
    out["cell_line"] = out["cell_line_raw"].fillna(out["cell_line"]).astype(str).map(standardize_cell_line)
    out["perturbation_subtype"] = subtype
    out["row_id"] = pd.to_numeric(out[row_id_column], errors="raise").astype("int64")
    out = out.loc[out["target_token"].ne("") & out["cell_line"].ne("")].copy()
    return out[["perturbation_subtype", "cell_line", "target_token", "row_id"]].reset_index(drop=True)


def build_pair_contract(lincs_membership: pd.DataFrame, sc_membership: pd.DataFrame) -> tuple[pd.DataFrame, pd.DataFrame]:
    merged = lincs_membership.merge(
        sc_membership,
        on=["perturbation_subtype", "cell_line", "target_token"],
        how="inner",
        suffixes=("_lincs", "_sc"),
    )
    merged = merged.rename(columns={"row_id_lincs": "lincs_row_idx", "row_id_sc": "sc_delta_row_idx"})
    merged = merged.sort_values(
        ["perturbation_subtype", "cell_line", "target_token", "lincs_row_idx", "sc_delta_row_idx"],
        kind="mergesort",
    ).reset_index(drop=True)
    proof = (
        merged.groupby(["perturbation_subtype", "cell_line", "target_token"], sort=False)
        .agg(
            n_lincs_rows=("lincs_row_idx", "nunique"),
            n_sc_rows=("sc_delta_row_idx", "nunique"),
            n_pairs=("lincs_row_idx", "size"),
        )
        .reset_index()
        .sort_values(["perturbation_subtype", "cell_line", "target_token"], kind="mergesort")
        .reset_index(drop=True)
    )
    return merged, proof


def main() -> int:
    args = parse_args()
    project_root = args.project_root.resolve()
    config_path = project_root / CONFIG_PATH
    config = yaml.safe_load(config_path.read_text(encoding="utf-8"))
    snapshot_root = resolve_config_path(project_root, args.snapshot_root)
    cross_root = snapshot_root / "cross_contract"
    cross_root.mkdir(parents=True, exist_ok=True)

    lincs_meta = pd.read_csv(snapshot_root / "lincs/lincs-engine1-meta.csv").reset_index().rename(columns={"index": "snapshot_row_idx"})
    chemical_meta = pd.read_csv(snapshot_root / "scperturb_delta/scperturb-chemical-delta-meta.csv")
    crispr_meta = pd.read_csv(snapshot_root / "scperturb_delta/scperturb-crispr-delta-meta.csv")

    lincs_chemical = explode_membership(
        lincs_meta.loc[lincs_meta["perturbation_subtype"].astype(str).eq("Chemical")].copy(),
        subtype="Chemical",
        row_id_column="snapshot_row_idx",
    )
    lincs_crispr = explode_membership(
        lincs_meta.loc[lincs_meta["perturbation_subtype"].astype(str).eq("CRISPR")].copy().assign(cell_line_raw=lincs_meta["cell_line_raw"]),
        subtype="CRISPR",
        row_id_column="snapshot_row_idx",
    )
    sc_chemical = explode_membership(
        chemical_meta.assign(cell_line=chemical_meta["cell_std"], cell_line_raw=chemical_meta["cell_line_raw"]),
        subtype="Chemical",
        row_id_column="delta_row_idx",
    )
    sc_crispr = explode_membership(
        crispr_meta.assign(cell_line=crispr_meta["cell_std"], cell_line_raw=crispr_meta["cell_line_raw"]),
        subtype="CRISPR",
        row_id_column="delta_row_idx",
    )

    chemical_contract, chemical_proof = build_pair_contract(lincs_chemical, sc_chemical)
    crispr_contract, crispr_proof = build_pair_contract(lincs_crispr, sc_crispr)

    write_csv(chemical_contract, cross_root / "cross-pairs-chemical-contract.csv")
    write_csv(crispr_contract, cross_root / "cross-pairs-crispr-contract.csv")
    write_csv(
        pd.concat([chemical_proof, crispr_proof], ignore_index=True, sort=False),
        cross_root / "alignment_proof.csv",
    )
    write_csv(
        pd.DataFrame(
            [
                {"perturbation_subtype": "Chemical", "matched_keys": int(len(chemical_proof)), "n_pairs": int(len(chemical_contract))},
                {"perturbation_subtype": "CRISPR", "matched_keys": int(len(crispr_proof)), "n_pairs": int(len(crispr_contract))},
            ]
        ),
        cross_root / "attrition_summary.csv",
    )
    write_json(
        cross_root / "exact_standardization_rule_audit.json",
        {
            "created_at": utc_now_iso(),
            "snapshot_root": str(snapshot_root),
            "cross_key": ["perturbation_subtype", "cell_line", "target_token"],
            "rules": {
                "cell_line": "uppercase whitespace-normalized string",
                "target_token": "uppercase token, chemical membership by token explode, CRISPR single-token only",
                "excluded_subtypes": ["sh", "oe"],
            },
            "matched_keys": {
                "Chemical": int(len(chemical_proof)),
                "CRISPR": int(len(crispr_proof)),
            },
        },
    )
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
