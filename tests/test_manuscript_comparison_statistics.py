from __future__ import annotations

import json
import sys
from pathlib import Path

import pandas as pd
import pytest

import scripts.manuscript_comparison_statistics as mcs
from scripts.manuscript_comparison_statistics import (
    ComparisonSpec,
    LOW_N_THRESHOLD,
    OUTPUT_COLUMNS,
    apply_family_bh,
    bh_adjust,
    build_manifest_entry,
    build_task1_contextual_long,
    expand_question_rows,
    parse_args,
    run_paired_tests_from_long,
    spec1_questions,
    update_manifest_registration,
    validate_output_schema,
)


def make_spec() -> ComparisonSpec:
    return ComparisonSpec(
        comparison_spec_id="test_spec",
        comparison_family="test_family",
        figure_id="FTEST",
        primary_source_object="/tmp/source.csv",
        pairing_type="paired",
        pairing_key_definition="shared unit",
        question_id_resolver=lambda _: ("Q1",),
        default_analysis_family="group_concordance",
    )


def test_duplicate_unit_group_rejection() -> None:
    frame = pd.DataFrame(
        {
            "analysis_family": ["group_concordance", "group_concordance", "group_concordance"],
            "metric_name": ["cosine", "cosine", "cosine"],
            "dataset": ["D1", "D1", "D1"],
            "cell_line": ["C1", "C1", "C1"],
            "group_label": ["Gene", "Gene", "Pathway"],
            "metric_value": [0.1, 0.2, 0.3],
        }
    )
    with pytest.raises(ValueError, match="duplicate rows"):
        run_paired_tests_from_long(
            long_frame=frame,
            spec=make_spec(),
            fixed_columns=["analysis_family", "metric_name"],
            unit_columns=["dataset", "cell_line"],
            group_column="group_label",
            value_column="metric_value",
            group_pairs=[("Gene", "Pathway")],
            duplicate_label="test duplicate",
        )


def test_pairing_completeness_uses_only_complete_pairs() -> None:
    frame = pd.DataFrame(
        {
            "analysis_family": ["group_concordance"] * 5,
            "metric_name": ["cosine"] * 5,
            "dataset": ["D1"] * 5,
            "cell_line": ["C1", "C1", "C2", "C2", "C3"],
            "group_label": ["Gene", "Pathway", "Gene", "Pathway", "Gene"],
            "metric_value": [0.9, 0.7, 0.8, 0.4, 0.5],
        }
    )
    result = run_paired_tests_from_long(
        long_frame=frame,
        spec=make_spec(),
        fixed_columns=["analysis_family", "metric_name"],
        unit_columns=["dataset", "cell_line"],
        group_column="group_label",
        value_column="metric_value",
        group_pairs=[("Gene", "Pathway")],
        duplicate_label="pair completeness",
    )
    assert int(result.loc[0, "n_test_units"]) == 2


def test_low_n_handling() -> None:
    frame = pd.DataFrame(
        {
            "analysis_family": ["group_concordance"] * 4,
            "metric_name": ["cosine"] * 4,
            "dataset": ["D1"] * 4,
            "cell_line": ["C1", "C1", "C2", "C2"],
            "group_label": ["Gene", "Pathway", "Gene", "Pathway"],
            "metric_value": [0.8, 0.5, 0.9, 0.4],
        }
    )
    result = run_paired_tests_from_long(
        long_frame=frame,
        spec=make_spec(),
        fixed_columns=["analysis_family", "metric_name"],
        unit_columns=["dataset", "cell_line"],
        group_column="group_label",
        value_column="metric_value",
        group_pairs=[("Gene", "Pathway")],
        duplicate_label="low n",
    )
    assert LOW_N_THRESHOLD == 5
    assert result.loc[0, "test_status"] == "not_tested_low_n"
    assert pd.isna(result.loc[0, "raw_p"])
    assert pd.isna(result.loc[0, "bh_q"])


def test_target_level_pairing_counts_cell_line_target_units() -> None:
    spec = ComparisonSpec(
        comparison_spec_id="f3_task2_performance_target_level_representation",
        comparison_family="figure3_task2_performance_target_level_representation",
        figure_id="F3",
        primary_source_object="/tmp/support.csv",
        pairing_type="paired",
        pairing_key_definition="shared dataset + cell_line + target",
        question_id_resolver=lambda _: ("Q3E",),
        default_analysis_family="retrieval",
    )
    frame = pd.DataFrame(
        [
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target": "T1",
                "analysis_family": "retrieval",
                "direction": "C2G",
                "metric_name": "mrr_corrected",
                "group_label": "Gene",
                "metric_value": 0.42,
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target": "T1",
                "analysis_family": "retrieval",
                "direction": "C2G",
                "metric_name": "mrr_corrected",
                "group_label": "Pathway",
                "metric_value": 0.31,
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target": "T2",
                "analysis_family": "retrieval",
                "direction": "C2G",
                "metric_name": "mrr_corrected",
                "group_label": "Gene",
                "metric_value": 0.39,
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target": "T2",
                "analysis_family": "retrieval",
                "direction": "C2G",
                "metric_name": "mrr_corrected",
                "group_label": "Pathway",
                "metric_value": 0.24,
            },
        ]
    )

    result = run_paired_tests_from_long(
        long_frame=frame,
        spec=spec,
        fixed_columns=["dataset", "analysis_family", "direction", "metric_name"],
        unit_columns=["cell_line", "target"],
        group_column="group_label",
        value_column="metric_value",
        group_pairs=[("Gene", "Pathway")],
        duplicate_label="Task2 target-level performance units",
    )

    assert result.loc[0, "pairing_key_definition"] == "shared dataset + cell_line + target"
    assert int(result.loc[0, "n_test_units"]) == 2


def test_bh_correction_within_family() -> None:
    frame = pd.DataFrame(
        {
            "comparison_spec_id": ["a", "b", "c"],
            "comparison_family": ["fam1", "fam1", "fam2"],
            "figure_id": ["F", "F", "F"],
            "analysis_family": ["group_concordance"] * 3,
            "metric_name": ["m1", "m2", "m3"],
            "group_a": ["Gene"] * 3,
            "group_b": ["Pathway"] * 3,
            "pairing_type": ["paired"] * 3,
            "pairing_key_definition": ["shared unit"] * 3,
            "n_test_units": [5, 5, 5],
            "raw_p": [0.01, 0.04, 0.03],
            "bh_q": [float("nan")] * 3,
            "significant_bool": [pd.NA] * 3,
            "effect_direction": ["group_a_gt_group_b"] * 3,
            "median_delta": [0.1, 0.1, 0.1],
            "test_name": ["wilcoxon_signed_rank"] * 3,
            "test_status": ["tested"] * 3,
            "notes": [""] * 3,
            "primary_source_object": ["/tmp/source.csv"] * 3,
            "dataset_scope": ["ALL"] * 3,
            "direction_scope": ["ALL"] * 3,
            "question_ids": [("Q1",), ("Q2",), ("Q3",)],
        }
    )
    adjusted = apply_family_bh(frame)
    fam1_expected = bh_adjust(pd.Series([0.01, 0.04]))
    assert adjusted.loc[0, "bh_q"] == pytest.approx(fam1_expected.iloc[0])
    assert adjusted.loc[1, "bh_q"] == pytest.approx(fam1_expected.iloc[1])
    assert adjusted.loc[2, "bh_q"] == pytest.approx(0.03)


def test_schema_validation() -> None:
    frame = pd.DataFrame([{column: pd.NA for column in OUTPUT_COLUMNS}])
    validate_output_schema(frame)
    with pytest.raises(ValueError, match="missing required columns"):
        validate_output_schema(frame.drop(columns=["notes"]))


def test_manifest_registration_behavior(tmp_path: Path) -> None:
    manifest_path = tmp_path / "framework_analysis_manifest.json"
    output_path = tmp_path / "manuscript_comparison_statistics.csv"
    payload = {
        "generated_at": "2026-03-19T00:00:00+00:00",
        "output_root": str(tmp_path),
        "source_roots": {},
        "contract_notes": [],
        "canonical_objects": [],
        "historical_noncanonical_objects": [],
    }
    manifest_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")
    output_path.write_text("comparison_spec_id\nstub\n", encoding="utf-8")

    update_manifest_registration(manifest_path, output_path)

    updated = json.loads(manifest_path.read_text(encoding="utf-8"))
    matches = [
        entry
        for entry in updated["canonical_objects"]
        if str(output_path) in {str(path) for path in entry.get("output_files", [])}
    ]
    assert len(matches) == 1
    assert matches[0]["canonical_filename"] == "manuscript_comparison_statistics.csv"

    output_path.unlink()
    with pytest.raises(FileNotFoundError):
        update_manifest_registration(manifest_path, output_path)


def test_task1_contextual_builder_collapses_perturbation_type_rows(tmp_path: Path) -> None:
    path = tmp_path / "figure3_task1_internal_contextual_support_summary.csv"
    pd.DataFrame(
        [
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "target": "ABCB1",
                "perturbation_type": "Chemical",
                "representation": "Gene",
                "analysis_family": "retrieval",
                "metric_name": "mrr_corrected",
                "metric_value": 0.8,
            },
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "target": "ABCB1",
                "perturbation_type": "Genetic",
                "representation": "Gene",
                "analysis_family": "retrieval",
                "metric_name": "mrr_corrected",
                "metric_value": 0.4,
            },
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "target": "ABCB1",
                "perturbation_type": "Chemical",
                "representation": "Pathway",
                "analysis_family": "retrieval",
                "metric_name": "mrr_corrected",
                "metric_value": 0.3,
            },
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "target": "ABCB1",
                "perturbation_type": "Genetic",
                "representation": "Pathway",
                "analysis_family": "retrieval",
                "metric_name": "mrr_corrected",
                "metric_value": 0.1,
            },
        ]
    ).to_csv(path, index=False)

    out = build_task1_contextual_long(path)

    assert len(out) == 2
    gene_value = out.loc[out["group_label"].eq("Gene"), "metric_value"].iloc[0]
    pathway_value = out.loc[out["group_label"].eq("Pathway"), "metric_value"].iloc[0]
    assert gene_value == pytest.approx(0.6)
    assert pathway_value == pytest.approx(0.2)


def test_spec1_questions_are_dataset_specific() -> None:
    assert spec1_questions({"dataset": "LINCS", "analysis_family": "group_concordance"}) == (
        "C2_F2a_internal_group_lincs",
    )
    assert spec1_questions({"dataset": "scPerturb", "analysis_family": "retrieval"}) == (
        "C2_F2a_internal_retrieval_scperturb_common_scope",
    )


def test_question_expansion_preserves_values() -> None:
    frame = pd.DataFrame(
        [
            {
                "comparison_spec_id": "spec",
                "comparison_family": "family",
                "figure_id": "F2",
                "analysis_family": "group_concordance",
                "metric_name": "cosine",
                "group_a": "Gene",
                "group_b": "Pathway",
                "pairing_type": "paired",
                "pairing_key_definition": "shared dataset + cell_line",
                "n_test_units": 6,
                "raw_p": 0.01,
                "bh_q": 0.02,
                "significant_bool": True,
                "effect_direction": "group_a_gt_group_b",
                "median_delta": 0.3,
                "test_name": "wilcoxon_signed_rank",
                "test_status": "tested",
                "notes": "",
                "primary_source_object": "/tmp/source.csv",
                "dataset_scope": "ALL",
                "direction_scope": "ALL",
                "question_ids": ("Q1", "Q2"),
            }
        ]
    )
    expanded = expand_question_rows(frame)
    assert list(expanded["question_id"]) == ["Q1", "Q2"]
    assert expanded["bh_q"].tolist() == [0.02, 0.02]



def test_parse_args_excludes_task1_heavy_inputs(monkeypatch: pytest.MonkeyPatch) -> None:
    monkeypatch.setattr(sys, "argv", ["prog"])
    args = parse_args()
    assert not hasattr(args, "task1_snapshot_root")
    assert not hasattr(args, "task1_retrieval_per_query_path")
    assert hasattr(args, "task2_retrieval_per_query_path")


def test_manifest_entry_excludes_removed_task1_fm_sources() -> None:
    entry = build_manifest_entry(Path("/tmp/manuscript_comparison_statistics.csv"))
    questions = set(entry["manuscript_questions"])
    assert "C2_F2c_internal_group_scperturb_fm" not in questions
    assert "C2_F2d_internal_retrieval_scperturb_fm" not in questions
    assert "C2_F2i_internal_to_cross_degradation_summary" in questions
    assert "C4_F3f_fm_retrieval_k562_local" in questions

    sources = set(entry["source_inputs"])
    assert "/home/lenislin/Experiment/projects/M2M/data/task1_snapshot_v1" not in sources
    assert not any(path.endswith("task1_retrieval_per_query.parquet") for path in sources)
    assert any(path.endswith("task2_retrieval_per_query.parquet") for path in sources)
    assert "Task1-FM comparison surfaces are excluded" in entry["fm_scope"]


def test_build_unique_results_uses_only_reduced_scope_specs(monkeypatch: pytest.MonkeyPatch) -> None:
    paired_specs: list[str] = []
    fm_specs: list[str] = []

    def fake_run_paired_tests_from_long(*, spec: ComparisonSpec, **_: object) -> pd.DataFrame:
        paired_specs.append(spec.comparison_spec_id)
        return pd.DataFrame([{"comparison_spec_id": spec.comparison_spec_id}])

    def fake_expand_fm_vs_base_tests(*, spec: ComparisonSpec, **_: object) -> pd.DataFrame:
        fm_specs.append(spec.comparison_spec_id)
        return pd.DataFrame([{"comparison_spec_id": spec.comparison_spec_id}])

    def empty_frame(columns: list[str]) -> pd.DataFrame:
        return pd.DataFrame({column: pd.Series(dtype="object") for column in columns})

    monkeypatch.setattr(mcs, "run_paired_tests_from_long", fake_run_paired_tests_from_long)
    monkeypatch.setattr(mcs, "expand_fm_vs_base_tests", fake_expand_fm_vs_base_tests)
    monkeypatch.setattr(
        mcs,
        "build_task1_bridge_summary_long",
        lambda *_: empty_frame(["dataset", "analysis_family", "metric_name", "group_label", "cell_line", "metric_value"]),
    )
    monkeypatch.setattr(
        mcs,
        "build_task1_bridge_detail_long",
        lambda *_: empty_frame(["representation_detail", "analysis_family", "metric_name", "dataset", "cell_line", "target", "group_label", "metric_value"]),
    )
    monkeypatch.setattr(
        mcs,
        "build_task2_direction_long",
        lambda *_: empty_frame(["analysis_family", "representation_label", "metric_name", "dataset", "cell_line", "direction", "metric_value"]),
    )
    monkeypatch.setattr(
        mcs,
        "build_task2_performance_long",
        lambda *_: empty_frame(["analysis_family", "direction", "metric_name", "dataset", "cell_line", "group_label", "metric_value"]),
    )
    monkeypatch.setattr(
        mcs,
        "build_task2_target_level_performance_long",
        lambda **_: empty_frame(["dataset", "cell_line", "target", "analysis_family", "direction", "metric_name", "group_label", "metric_value"]),
    )
    monkeypatch.setattr(
        mcs,
        "build_task2_cell_line_pattern_long",
        lambda *_: empty_frame(["dataset", "analysis_family", "direction", "metric_name", "cell_line", "group_label", "metric_value"]),
    )
    monkeypatch.setattr(
        mcs,
        "build_task2_target_pattern_long",
        lambda *_: empty_frame(["dataset", "analysis_family", "direction", "metric_name", "target", "group_label", "metric_value"]),
    )
    monkeypatch.setattr(
        mcs,
        "build_task1_contextual_long",
        lambda *_: empty_frame(["analysis_family", "metric_name", "dataset", "cell_line", "target", "group_label", "metric_value"]),
    )
    monkeypatch.setattr(
        mcs,
        "build_task2_k562_c2g_query_long",
        lambda *_: empty_frame(["analysis_family", "direction", "metric_name", "dataset", "cell_line", "query_row_id", "group_label", "metric_value"]),
    )

    out = mcs.build_unique_results(
        task1_bridge_summary_path=Path("/tmp/task1_bridge_summary.csv"),
        task1_bridge_detail_path=Path("/tmp/task1_bridge_detail.csv"),
        task2_retrieval_per_query_path=Path("/tmp/task2_retrieval_per_query.parquet"),
        task2_group_root=Path("/tmp/task2_group_root"),
        analysis_root=Path("/tmp/analysis"),
    )

    assert paired_specs == [
        "f2_task1_internal_common_representation",
        "f2_task1_cross_common_representation",
        "f2_task1_internal_to_cross_degradation",
        "f3_task2_direction_retrieval",
        "f3_task2_common_representation",
        "f3_task2_performance_target_level_representation",
        "f3_task2_cell_line_pattern_representation",
        "f3_task2_target_pattern_representation",
        "f3_task1_contextual_common_representation",
    ]
    assert fm_specs == [
        "f3_task2_k562_fm_target_anchored",
        "f3_task2_k562_fm_c2g_query_anchored",
    ]
    assert out["comparison_spec_id"].tolist() == paired_specs + fm_specs
    assert "f2_task1_scperturb_internal_fm_group" not in out["comparison_spec_id"].tolist()
    assert "f2_task1_scperturb_internal_fm_retrieval" not in out["comparison_spec_id"].tolist()
    assert "f3_task1_contextual_fm_representation" not in out["comparison_spec_id"].tolist()
