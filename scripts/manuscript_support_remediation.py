from __future__ import annotations

from pathlib import Path

import pandas as pd
import pyarrow.parquet as pq


ROOT = Path(__file__).resolve().parents[1]
DOCS = ROOT / "docs"
RUNS = ROOT / "runs"
DATA = ROOT / "data"


def split_pipe(value: str) -> list[str]:
    return [item.strip() for item in str(value).split("|") if item.strip()]


def split_source_files(value: str) -> list[str]:
    return [item.strip() for item in str(value).split(" | ") if item.strip()]


def read_columns(rel_path: str) -> set[str]:
    path = ROOT / rel_path
    if path.suffix == ".csv":
        return set(pd.read_csv(path, nrows=0).columns)
    if path.suffix == ".parquet":
        return set(pq.read_schema(path).names)
    return set()


def manuscript_role(row: pd.Series) -> str:
    if row["claim_id"] == "C2_F2g_cross_chemical_policy_exclusion":
        return "policy_exclusion_note"
    if row["claim_id"] == "C4_F3g_fm_scope_policy":
        return "scope_disclosure_note"
    if row["task"] in {"Task2_group", "Task2_FM_scope"} and "mean_cosine_centroid" in row["metric"]:
        return "primary_group_evidence"
    if row["task"] == "Task1_cross" and "mean_cosine_centroid" in row["metric"]:
        return "primary_group_evidence"
    if row["task"] == "Task1_internal" and "mean_cosine_centroid" in row["metric"]:
        return "primary_group_evidence"
    return "supporting_retrieval_evidence"


def scope_restrictions(row: pd.Series) -> str:
    claim_id = row["claim_id"]
    if claim_id == "C2_F2g_cross_chemical_policy_exclusion":
        return "cross chemical excluded by frozen policy/support gating; numeric rationale not manuscript-safe"
    if claim_id.startswith("C2_F2") and row["task"] == "Task1_cross":
        return "genetic-only cross scope for retained metric rows; keep retrieval directions separate; cross alignment contract required"
    if claim_id.startswith("C2_F2") and row["task"] == "Task1_internal":
        return "internal modality-preserving scope; common Gene/Pathway comparison only for retained main-text Task1 rows"
    if claim_id in {"C3_F3a_group_common_scope"}:
        return "dataset/cell_line stratified; Gene/Pathway common scope only; no raw dataset pooling; edist remains informational-only"
    if claim_id in {"C3_F3b_retrieval_c2g_common_scope", "C3_F3c_retrieval_g2c_common_scope"}:
        return "dataset/cell_line stratified; Gene/Pathway common scope only; keep directions separate; no raw dataset pooling"
    if claim_id == "C4_F3e_fm_group_k562_local":
        return "scPerturb K562 local FM scope only; group-level evidence primary; no benchmark-wide FM generalization"
    if claim_id == "C4_F3g_fm_scope_policy":
        return "scPerturb K562 FM scope only; LINCS FM absence is not_applicable_scope rather than attrition"
    return row["stratification"]


def build_panel_support_ledger() -> int:
    claim_map = pd.read_csv(DOCS / "manuscript_claim_to_evidence_map.csv")
    ledger_rows: list[dict[str, str]] = []

    for _, row in claim_map.iterrows():
        if row["support_decision"] != "main_text_candidate":
            continue

        available_columns: set[str] = set()
        source_files = split_source_files(row["source_file"])
        for rel_path in source_files:
            available_columns |= read_columns(rel_path)

        denom_fields = split_pipe(row["denominator_fields"])
        missing_fields = [field for field in denom_fields if field not in available_columns]
        if missing_fields:
            denom_values_available = "partial_missing:" + "|".join(missing_fields)
        else:
            denom_values_available = "all_present"

        notes = str(row["notes"])
        if missing_fields:
            notes += " Missing denominator columns in source inspection: " + ",".join(missing_fields) + "."

        ledger_rows.append(
            {
                "claim_id": row["claim_id"],
                "figure_panel": row["figure_panel"],
                "source_stage": row["source_stage"],
                "source_file": row["source_file"],
                "metric": row["metric"],
                "representation": row["representation"],
                "stratification": row["stratification"],
                "denominator_fields": row["denominator_fields"],
                "denominator_values_available": denom_values_available,
                "caution_codes": row["caution_codes"],
                "scope_restrictions": scope_restrictions(row),
                "support_decision": row["support_decision"],
                "manuscript_role": manuscript_role(row),
                "notes": notes,
            }
        )

    ledger = pd.DataFrame(ledger_rows).sort_values(["figure_panel", "claim_id"])
    ledger.to_csv(DOCS / "manuscript_panel_support_ledger.csv", index=False)
    return len(ledger)


def sparse_note(n_queries: int, base: str) -> str:
    if n_queries < 100:
        return base + "; sparse_support"
    return base


def build_n_targets_sensitivity() -> int:
    query_path = (
        RUNS
        / "s5_multisource_impl_verify_20260311_a"
        / "s5_task2_retrieval_multisource"
        / "task2_retrieval_per_query.parquet"
    )
    df = pd.read_parquet(
        query_path,
        columns=[
            "dataset",
            "cell_line",
            "representation",
            "direction",
            "query_n_targets",
            "mrr_corrected",
            "hit1_corrected",
            "hit5_corrected",
            "hit10_corrected",
            "N_gallery",
            "m_pos",
        ],
    )
    df = df[df["direction"] == "C2G"].copy()

    out = (
        df.groupby(
            ["dataset", "cell_line", "representation", "direction", "query_n_targets"],
            dropna=False,
            as_index=False,
        )
        .agg(
            n_queries=("query_n_targets", "size"),
            mean_mrr_corrected=("mrr_corrected", "mean"),
            mean_hit1_corrected=("hit1_corrected", "mean"),
            mean_hit5_corrected=("hit5_corrected", "mean"),
            mean_hit10_corrected=("hit10_corrected", "mean"),
            mean_N_gallery=("N_gallery", "mean"),
            mean_m_pos=("m_pos", "mean"),
        )
        .sort_values(["dataset", "cell_line", "representation", "query_n_targets"])
    )

    base_note = "downstream manuscript sensitivity; exact query_n_targets; C2G only; dataset/cell_line explicit"
    out["sensitivity_scope_note"] = out["n_queries"].map(lambda n: sparse_note(int(n), base_note))
    out.to_csv(DOCS / "task2_c2g_query_n_targets_sensitivity.csv", index=False)
    return len(out)


def build_specificity_tier_sensitivity() -> int:
    query_path = (
        RUNS
        / "s5_multisource_impl_verify_20260311_a"
        / "s5_task2_retrieval_multisource"
        / "task2_retrieval_per_query.parquet"
    )
    meta_path = (
        DATA
        / "task2_snapshot_v2"
        / "scperturb_k562"
        / "derived"
        / "delta_meta.csv"
    )

    queries = pd.read_parquet(
        query_path,
        columns=[
            "dataset",
            "cell_line",
            "representation",
            "direction",
            "query_row_id",
            "mrr_corrected",
            "hit1_corrected",
            "hit5_corrected",
            "hit10_corrected",
            "N_gallery",
            "m_pos",
        ],
    )
    queries = queries[
        (queries["dataset"] == "scPerturb")
        & (queries["cell_line"] == "K562")
        & (queries["direction"] == "C2G")
    ].copy()

    meta = pd.read_csv(meta_path, usecols=["row_id", "specificity_tier"]).rename(
        columns={"row_id": "query_row_id"}
    )
    merged = queries.merge(meta, on="query_row_id", how="left")

    out = (
        merged.groupby(
            ["dataset", "cell_line", "representation", "direction", "specificity_tier"],
            dropna=False,
            as_index=False,
        )
        .agg(
            n_queries=("query_row_id", "size"),
            mean_mrr_corrected=("mrr_corrected", "mean"),
            mean_hit1_corrected=("hit1_corrected", "mean"),
            mean_hit5_corrected=("hit5_corrected", "mean"),
            mean_hit10_corrected=("hit10_corrected", "mean"),
            mean_N_gallery=("N_gallery", "mean"),
            mean_m_pos=("m_pos", "mean"),
        )
        .sort_values(["representation", "specificity_tier"], na_position="last")
    )

    out["specificity_tier"] = out["specificity_tier"].fillna("")
    base_note = "downstream manuscript sensitivity; local scPerturb K562 C2G only; not benchmark-wide; tier counts may be unbalanced"
    out["sensitivity_scope_note"] = out["n_queries"].map(lambda n: sparse_note(int(n), base_note))
    out.to_csv(DOCS / "task2_scperturb_k562_c2g_specificity_tier_sensitivity.csv", index=False)
    return len(out)


def main() -> None:
    ledger_rows = build_panel_support_ledger()
    n_targets_rows = build_n_targets_sensitivity()
    specificity_rows = build_specificity_tier_sensitivity()
    print(f"manuscript_panel_support_ledger.csv rows: {ledger_rows}")
    print(f"task2_c2g_query_n_targets_sensitivity.csv rows: {n_targets_rows}")
    print(f"task2_scperturb_k562_c2g_specificity_tier_sensitivity.csv rows: {specificity_rows}")


if __name__ == "__main__":
    main()
