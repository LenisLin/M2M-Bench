from __future__ import annotations

import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

import scripts.manuscript_plot_ready_tables as mpr
import scripts.manuscript_task1_high_concordance as mthc


def test_build_figure2_panel_2a_preserves_cross_eligibility_all_bin(tmp_path: Path) -> None:
    source_path = tmp_path / "figure2_task1_scope_summary.csv"
    output_path = tmp_path / "figure2_panel_2a_task1_scope.csv"
    pd.DataFrame(
        [
            {
                "scope": "cross",
                "dataset_or_direction": "cross",
                "perturbation_type": "Chemical",
                "representation": "ALL",
                "analysis_family": "cross_eligibility",
                "scope_status": "excluded_by_support_gate",
                "metric_names": pd.NA,
                "metric_count": pd.NA,
                "n_total": pd.NA,
                "n_valid": pd.NA,
                "n_excluded": pd.NA,
                "N_gallery_max": pd.NA,
                "cross_alignment_contract": pd.NA,
                "n_matched_keys": 0,
                "eligible_bool": False,
                "exclusion_reason": "matched_keys_lt5",
                "scope_note": "Cross chemical excluded by policy.",
                "fm_scope_note": "FM excluded in cross scope.",
            },
            {
                "scope": "internal",
                "dataset_or_direction": "scPerturb",
                "perturbation_type": "Genetic",
                "representation": "FM:geneformer",
                "analysis_family": "group_concordance",
                "scope_status": "materialized",
                "metric_names": "mean_cosine_centroid",
                "metric_count": 1,
                "n_total": 20,
                "n_valid": 20,
                "n_excluded": 0,
                "N_gallery_max": pd.NA,
                "cross_alignment_contract": pd.NA,
                "n_matched_keys": pd.NA,
                "eligible_bool": pd.NA,
                "exclusion_reason": pd.NA,
                "scope_note": "Internal FM row.",
                "fm_scope_note": "FM local only.",
            },
        ]
    ).to_csv(source_path, index=False)

    mpr.build_figure2_panel_2a(source_path=source_path, output_path=output_path)

    out = pd.read_csv(output_path)
    eligibility_row = out.loc[out["analysis_family"].eq("cross_eligibility")].iloc[0]
    fm_row = out.loc[out["representation_class"].eq("FM")].iloc[0]
    assert eligibility_row["representation_class"] == "ALL"
    assert eligibility_row["coverage_denominator"] == 1
    assert eligibility_row["coverage_fraction"] == pytest.approx(0.0)
    assert eligibility_row["count_annotation_type"] == "n_matched_keys"
    assert fm_row["representation_class"] == "FM"
    assert fm_row["coverage_denominator"] == 3
    assert fm_row["coverage_fraction"] == pytest.approx(1 / 3)
    assert fm_row["count_annotation"] == pytest.approx(20)
    assert fm_row["count_annotation_type"] == "n_valid"


def test_build_task1_performance_overview_filters_scope_and_adds_fill_value(tmp_path: Path) -> None:
    source_path = tmp_path / "figure2_task1_performance_structure.csv"
    output_path = tmp_path / "figure2_task1_performance_overview.csv"
    pd.DataFrame(
        [
            {
                "scope": "cross",
                "dataset_or_direction": "LINCS_to_scPerturb",
                "perturbation_type": "Genetic",
                "representation": "Gene",
                "analysis_family": "group_concordance",
                "metric_name": "mean_cosine_centroid",
                "value_variant": "group",
                "value": 0.2,
                "n_total": 10,
                "n_valid": 10,
                "n_excluded": 0,
                "N_gallery_max": pd.NA,
                "cross_alignment_contract": "x",
                "chance_check_available_bool": False,
                "chance_abs_delta": pd.NA,
                "fm_scope_note": "x",
            },
            {
                "scope": "cross",
                "dataset_or_direction": "LINCS_to_scPerturb",
                "perturbation_type": "Genetic",
                "representation": "Pathway",
                "analysis_family": "group_concordance",
                "metric_name": "mean_cosine_centroid",
                "value_variant": "group",
                "value": 0.8,
                "n_total": 10,
                "n_valid": 10,
                "n_excluded": 0,
                "N_gallery_max": pd.NA,
                "cross_alignment_contract": "x",
                "chance_check_available_bool": False,
                "chance_abs_delta": pd.NA,
                "fm_scope_note": "x",
            },
            {
                "scope": "cross",
                "dataset_or_direction": "LINCS_to_scPerturb",
                "perturbation_type": "Chemical",
                "representation": "Gene",
                "analysis_family": "group_concordance",
                "metric_name": "mean_cosine_centroid",
                "value_variant": "group",
                "value": 0.9,
                "n_total": 10,
                "n_valid": 10,
                "n_excluded": 0,
                "N_gallery_max": pd.NA,
                "cross_alignment_contract": "x",
                "chance_check_available_bool": False,
                "chance_abs_delta": pd.NA,
                "fm_scope_note": "x",
            },
            {
                "scope": "cross",
                "dataset_or_direction": "LINCS_to_scPerturb",
                "perturbation_type": "Genetic",
                "representation": "Gene",
                "analysis_family": "retrieval",
                "metric_name": "mrr",
                "value_variant": "retrieval_raw",
                "value": 0.9,
                "n_total": 10,
                "n_valid": 10,
                "n_excluded": 0,
                "N_gallery_max": 10,
                "cross_alignment_contract": "x",
                "chance_check_available_bool": True,
                "chance_abs_delta": 0.0,
                "fm_scope_note": "x",
            },
            {
                "scope": "cross",
                "dataset_or_direction": "LINCS_to_scPerturb",
                "perturbation_type": "Genetic",
                "representation": "Gene",
                "analysis_family": "retrieval",
                "metric_name": "mrr",
                "value_variant": "retrieval_corrected",
                "value": 0.4,
                "n_total": 10,
                "n_valid": 10,
                "n_excluded": 0,
                "N_gallery_max": 10,
                "cross_alignment_contract": "x",
                "chance_check_available_bool": True,
                "chance_abs_delta": 0.0,
                "fm_scope_note": "x",
            },
        ]
    ).to_csv(source_path, index=False)

    mpr.build_task1_performance_overview(
        source_path=source_path,
        scope_value="cross",
        output_path=output_path,
    )

    out = pd.read_csv(output_path)
    assert set(out["perturbation_type"]) == {"Genetic"}
    assert set(out["value_variant"]) == {"group", "retrieval_corrected"}
    group_rows = out.loc[out["analysis_family"].eq("group_concordance")].sort_values("representation")
    assert group_rows["fill_value"].tolist() == pytest.approx([0.5, 1.0])


def test_build_figure2_panel_2c_builds_raw_and_summary_rows(tmp_path: Path) -> None:
    bridge_detail_path = tmp_path / "task1_internal_vs_cross_group_bridge.csv"
    stats_path = tmp_path / "manuscript_comparison_statistics.csv"
    output_path = tmp_path / "figure2_panel_2c_gene_vs_pathway_matched_units.csv"
    pd.DataFrame(
        [
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "target": "T1",
                "metric_family": "group_concordance",
                "metric_name": "cosine",
                "representation_detail": "Gene",
                "internal_value": 0.8,
                "cross_value": 0.6,
            },
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "target": "T1",
                "metric_family": "group_concordance",
                "metric_name": "cosine",
                "representation_detail": "Pathway",
                "internal_value": 0.5,
                "cross_value": 0.4,
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target": "T2",
                "metric_family": "group_concordance",
                "metric_name": "cosine",
                "representation_detail": "Gene",
                "internal_value": 0.7,
                "cross_value": 0.55,
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target": "T2",
                "metric_family": "group_concordance",
                "metric_name": "cosine",
                "representation_detail": "Pathway",
                "internal_value": 0.45,
                "cross_value": 0.35,
            },
        ]
    ).to_csv(bridge_detail_path, index=False)
    pd.DataFrame(
        [
            {
                "comparison_family": "figure2_task1_internal_common_representation",
                "analysis_family": "group_concordance",
                "metric_name": "cosine",
                "dataset_scope": "LINCS",
                "bh_q": 0.02,
                "test_status": "tested",
                "effect_direction": "group_a_gt_group_b",
                "median_delta": 0.3,
                "n_test_units": 2,
            },
            {
                "comparison_family": "figure2_task1_internal_common_representation",
                "analysis_family": "group_concordance",
                "metric_name": "cosine",
                "dataset_scope": "scPerturb",
                "bh_q": 0.03,
                "test_status": "tested",
                "effect_direction": "group_a_gt_group_b",
                "median_delta": 0.2,
                "n_test_units": 2,
            },
            {
                "comparison_family": "figure2_task1_cross_common_representation",
                "analysis_family": "group_concordance",
                "metric_name": "cosine",
                "dataset_scope": "LINCS",
                "bh_q": 0.04,
                "test_status": "tested",
                "effect_direction": "group_a_gt_group_b",
                "median_delta": 0.1,
                "n_test_units": 1,
            },
            {
                "comparison_family": "figure2_task1_cross_common_representation",
                "analysis_family": "group_concordance",
                "metric_name": "cosine",
                "dataset_scope": "scPerturb",
                "bh_q": 0.04,
                "test_status": "tested",
                "effect_direction": "group_a_gt_group_b",
                "median_delta": 0.1,
                "n_test_units": 1,
            },
        ]
    ).to_csv(stats_path, index=False)

    mpr.build_figure2_panel_2c(
        bridge_detail_path=bridge_detail_path,
        comparison_stats_path=stats_path,
        output_path=output_path,
    )

    out = pd.read_csv(output_path)
    assert set(out["row_kind"]) == {"raw", "summary"}
    assert set(out["comparison_context"]) == {
        "Internal LINCS",
        "Internal scPerturb",
        "LINCS -> scPerturb",
        "scPerturb -> LINCS",
    }
    summary = out.loc[out["row_kind"].eq("summary")].copy()
    assert {"q10_value", "q25_value", "q50_value", "q75_value", "q90_value"}.issubset(summary.columns)
    assert summary["q50_value"].notna().all()
    assert summary["metric_value"].equals(summary["q50_value"])


def test_build_figure2_panel_2d_uses_comparison_stats_and_long_scope_schema(tmp_path: Path) -> None:
    source_path = tmp_path / "task1_internal_vs_cross_group_bridge.csv"
    stats_path = tmp_path / "manuscript_comparison_statistics.csv"
    output_path = tmp_path / "figure2_panel_2d_internal_to_cross_degradation.csv"
    pd.DataFrame(
        [
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "target": "T1",
                "metric_family": "group_concordance",
                "metric_name": "cosine",
                "representation_detail": "Gene",
                "internal_value": 0.8,
                "cross_value": 0.6,
            },
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "target": "T1",
                "metric_family": "group_concordance",
                "metric_name": "cosine",
                "representation_detail": "Pathway",
                "internal_value": 0.5,
                "cross_value": 0.4,
            },
            {
                "dataset": "LINCS",
                "cell_line": "MCF7",
                "target": "T2",
                "metric_family": "group_concordance",
                "metric_name": "cosine",
                "representation_detail": "Gene",
                "internal_value": 0.7,
                "cross_value": 0.55,
            },
            {
                "dataset": "LINCS",
                "cell_line": "MCF7",
                "target": "T2",
                "metric_family": "group_concordance",
                "metric_name": "cosine",
                "representation_detail": "Pathway",
                "internal_value": 0.45,
                "cross_value": 0.35,
            },
        ]
    ).to_csv(source_path, index=False)
    pd.DataFrame(
        [
            {
                "comparison_family": "figure2_task1_internal_to_cross_degradation",
                "dataset_scope": "LINCS",
                "representation_scope": "Gene",
                "analysis_family": "group_concordance",
                "metric_name": "cosine",
                "bh_q": 0.04,
                "test_status": "tested",
                "effect_direction": "group_b_gt_group_a",
                "median_delta": -0.15,
                "n_test_units": 2,
            },
            {
                "comparison_family": "figure2_task1_internal_to_cross_degradation",
                "dataset_scope": "LINCS",
                "representation_scope": "Pathway",
                "analysis_family": "group_concordance",
                "metric_name": "cosine",
                "bh_q": 0.06,
                "test_status": "tested",
                "effect_direction": "group_b_gt_group_a",
                "median_delta": -0.10,
                "n_test_units": 2,
            },
        ]
    ).to_csv(stats_path, index=False)

    mpr.build_figure2_panel_2d(
        source_path=source_path,
        comparison_stats_path=stats_path,
        output_path=output_path,
    )

    out = pd.read_csv(output_path)
    assert set(out["scope"]) == {"Internal", "Cross"}
    assert {"summary_mean", "summary_median", "q10_value", "q25_value", "q50_value", "q75_value", "q90_value"}.issubset(out.columns)
    gene_rows = out.loc[out["representation"].eq("Gene")].sort_values("scope", kind="mergesort")
    assert set(gene_rows["n_test_units"]) == {2}
    assert gene_rows["summary_median"].notna().all()


def test_build_figure3_panel_3a_collapses_fm_scope_rows(tmp_path: Path) -> None:
    source_path = tmp_path / "figure3_task2_scope_summary.csv"
    output_path = tmp_path / "figure3_panel_3a_task2_scope.csv"
    pd.DataFrame(
        [
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "representation": "Gene",
                "availability_status": "available",
                "availability_reason": "base",
                "n_targets_eligible": 10,
                "n_targets_ineligible": 1,
                "group_slice_materialized_bool": True,
                "retrieval_c2g_materialized_bool": True,
                "retrieval_g2c_materialized_bool": True,
                "retrieval_directions_present": "C2G|G2C",
                "n_targets_total": 11,
                "n_attrition_target_rows": 1,
                "group_attrition_rows": 1,
                "scope_note": "base",
                "fm_scope_note": "base",
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "representation": "geneformer",
                "availability_status": "available",
                "availability_reason": "fm1",
                "n_targets_eligible": 9,
                "n_targets_ineligible": 2,
                "group_slice_materialized_bool": True,
                "retrieval_c2g_materialized_bool": True,
                "retrieval_g2c_materialized_bool": True,
                "retrieval_directions_present": "C2G|G2C",
                "n_targets_total": 11,
                "n_attrition_target_rows": 2,
                "group_attrition_rows": 2,
                "scope_note": "fm",
                "fm_scope_note": "fm",
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "representation": "scgpt",
                "availability_status": "not_applicable_scope",
                "availability_reason": "fm2",
                "n_targets_eligible": 8,
                "n_targets_ineligible": 3,
                "group_slice_materialized_bool": False,
                "retrieval_c2g_materialized_bool": False,
                "retrieval_g2c_materialized_bool": False,
                "retrieval_directions_present": pd.NA,
                "n_targets_total": 11,
                "n_attrition_target_rows": 3,
                "group_attrition_rows": 3,
                "scope_note": "fm",
                "fm_scope_note": "fm",
            },
        ]
    ).to_csv(source_path, index=False)

    mpr.build_figure3_panel_3a(source_path=source_path, output_path=output_path)

    out = pd.read_csv(output_path)
    assert set(out["representation_class"]) == {"Gene", "Pathway", "FM"}
    fm_row = out.loc[out["representation_class"].eq("FM")].iloc[0]
    assert fm_row["status_group"] == "available"
    assert fm_row["eligible_target_count"] == 10
    assert fm_row["coverage_denominator"] == 3
    assert fm_row["coverage_fraction"] == pytest.approx(1 / 3)
    pathway_row = out.loc[out["representation_class"].eq("Pathway")].iloc[0]
    assert pathway_row["coverage_fraction"] == pytest.approx(0.0)
    assert pathway_row["count_annotation_type"] == "eligible_target_count"


def test_build_figure3_panel_3b_keeps_lincs_common_scope_and_scperturb_tier_rows(
    tmp_path: Path,
) -> None:
    performance_path = tmp_path / "figure3_task2_performance_structure.csv"
    direction_support_path = tmp_path / "figure3_task2_direction_support_summary.csv"
    retrieval_per_query_path = tmp_path / "task2_retrieval_per_query.parquet"
    drug_meta_path = tmp_path / "Drug_meta.csv"
    output_path = tmp_path / "figure3_panel_3b_c2g_performance_overview.csv"
    pd.DataFrame(
        [
            {
                "analysis_family": "group_concordance",
                "dataset": "LINCS",
                "cell_line": "A375",
                "direction": pd.NA,
                "representation": "Gene",
                "metric_name": "mean_cosine_centroid",
                "metric_value": 0.2,
                "n_targets_total": 10,
                "n_targets_metric_valid": 10,
                "performance_scope_note": "group",
                "fm_scope_note": "base",
            },
            {
                "analysis_family": "group_concordance",
                "dataset": "LINCS",
                "cell_line": "MCF7",
                "direction": pd.NA,
                "representation": "Pathway",
                "metric_name": "mean_cosine_centroid",
                "metric_value": 0.8,
                "n_targets_total": 10,
                "n_targets_metric_valid": 10,
                "performance_scope_note": "group",
                "fm_scope_note": "base",
            },
            {
                "analysis_family": "group_concordance",
                "dataset": "scPerturb",
                "cell_line": "K562",
                "direction": pd.NA,
                "representation": "Gene",
                "metric_name": "mean_cosine_centroid",
                "metric_value": 0.95,
                "n_targets_total": 10,
                "n_targets_metric_valid": 10,
                "performance_scope_note": "group",
                "fm_scope_note": "local",
            },
            {
                "analysis_family": "group_concordance",
                "dataset": "scPerturb",
                "cell_line": "K562",
                "direction": pd.NA,
                "representation": "Pathway",
                "metric_name": "mean_cosine_centroid",
                "metric_value": 0.85,
                "n_targets_total": 10,
                "n_targets_metric_valid": 10,
                "performance_scope_note": "group",
                "fm_scope_note": "local",
            },
            {
                "analysis_family": "group_concordance",
                "dataset": "LINCS",
                "cell_line": "A375",
                "direction": pd.NA,
                "representation": "geneformer",
                "metric_name": "mean_cosine_centroid",
                "metric_value": 0.9,
                "n_targets_total": 10,
                "n_targets_metric_valid": 10,
                "performance_scope_note": "group",
                "fm_scope_note": "fm",
            },
        ]
    ).to_csv(performance_path, index=False)
    pd.DataFrame(
        [
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "direction": "C2G",
                "representation": "Gene",
                "n_total": 100,
                "n_valid": 100,
                "n_excluded_missing_metric_or_mpos0": 0,
                "N_gallery_mean": 10,
                "N_gallery_max": 10,
                "m_pos_mean": 2,
                "m_pos_p50": 2,
                "m_pos_p90": 3,
                "mean_mrr_corrected": 0.1,
                "mean_hit1_corrected": 0.2,
                "mean_hit5_corrected": 0.3,
                "mean_hit10_corrected": 0.4,
                "fm_scope_note": "base",
            },
            {
                "dataset": "LINCS",
                "cell_line": "MCF7",
                "direction": "C2G",
                "representation": "Gene",
                "n_total": 100,
                "n_valid": 100,
                "n_excluded_missing_metric_or_mpos0": 0,
                "N_gallery_mean": 10,
                "N_gallery_max": 10,
                "m_pos_mean": 2,
                "m_pos_p50": 2,
                "m_pos_p90": 3,
                "mean_mrr_corrected": 0.9,
                "mean_hit1_corrected": 0.8,
                "mean_hit5_corrected": 0.7,
                "mean_hit10_corrected": 0.6,
                "fm_scope_note": "base",
            },
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "direction": "C2G",
                "representation": "Pathway",
                "n_total": 100,
                "n_valid": 100,
                "n_excluded_missing_metric_or_mpos0": 0,
                "N_gallery_mean": 10,
                "N_gallery_max": 10,
                "m_pos_mean": 2,
                "m_pos_p50": 2,
                "m_pos_p90": 3,
                "mean_mrr_corrected": 0.15,
                "mean_hit1_corrected": 0.25,
                "mean_hit5_corrected": 0.35,
                "mean_hit10_corrected": 0.45,
                "fm_scope_note": "base",
            },
            {
                "dataset": "LINCS",
                "cell_line": "MCF7",
                "direction": "C2G",
                "representation": "Pathway",
                "n_total": 100,
                "n_valid": 100,
                "n_excluded_missing_metric_or_mpos0": 0,
                "N_gallery_mean": 10,
                "N_gallery_max": 10,
                "m_pos_mean": 2,
                "m_pos_p50": 2,
                "m_pos_p90": 3,
                "mean_mrr_corrected": 0.75,
                "mean_hit1_corrected": 0.65,
                "mean_hit5_corrected": 0.55,
                "mean_hit10_corrected": 0.45,
                "fm_scope_note": "base",
            },
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "direction": "G2C",
                "representation": "Gene",
                "n_total": 100,
                "n_valid": 100,
                "n_excluded_missing_metric_or_mpos0": 0,
                "N_gallery_mean": 10,
                "N_gallery_max": 10,
                "m_pos_mean": 2,
                "m_pos_p50": 2,
                "m_pos_p90": 3,
                "mean_mrr_corrected": 0.5,
                "mean_hit1_corrected": 0.5,
                "mean_hit5_corrected": 0.5,
                "mean_hit10_corrected": 0.5,
                "fm_scope_note": "base",
            },
        ]
    ).to_csv(direction_support_path, index=False)
    pd.DataFrame(
        [
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "direction": "C2G",
                "representation": "Gene",
                "query_row_id": 10,
                "treated_cell_id": "drug-1",
                "query_specificity_tier": "The Cleanest Hits",
                "N_gallery": 12,
                "m_pos": 2,
                "mrr_corrected": 0.91,
                "hit1_corrected": 0.81,
                "hit5_corrected": 0.71,
                "hit10_corrected": 0.61,
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "direction": "C2G",
                "representation": "Pathway",
                "query_row_id": 11,
                "treated_cell_id": "drug-2",
                "query_specificity_tier": "The Family Hits",
                "N_gallery": 14,
                "m_pos": 4,
                "mrr_corrected": 0.73,
                "hit1_corrected": 0.63,
                "hit5_corrected": 0.53,
                "hit10_corrected": 0.43,
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "direction": "C2G",
                "representation": "Gene",
                "query_row_id": 12,
                "treated_cell_id": "control-1",
                "query_specificity_tier": "Control",
                "N_gallery": 9,
                "m_pos": 1,
                "mrr_corrected": 0.05,
                "hit1_corrected": 0.04,
                "hit5_corrected": 0.03,
                "hit10_corrected": 0.02,
            },
        ]
    ).to_parquet(retrieval_per_query_path, index=False)
    pd.DataFrame(
        [
            {
                "treated_cell_id": "drug-1",
                "benchmark_group": "Test",
                "specificity_tier": "The Cleanest Hits",
            },
            {
                "treated_cell_id": "drug-2",
                "benchmark_group": "Test",
                "specificity_tier": "The Family Hits",
            },
            {
                "treated_cell_id": "control-1",
                "benchmark_group": "Control",
                "specificity_tier": "Control",
            },
        ]
    ).to_csv(drug_meta_path, index=False)

    mpr.build_figure3_panel_3b(
        performance_path=performance_path,
        direction_support_path=direction_support_path,
        task2_retrieval_per_query_path=retrieval_per_query_path,
        drug_meta_path=drug_meta_path,
        output_path=output_path,
    )

    out = pd.read_csv(output_path)
    assert set(out["analysis_family"]) == {"group_concordance", "retrieval"}
    assert set(out["dataset"]) == {"LINCS", "scPerturb"}
    assert "scPerturb" not in set(out.loc[out["analysis_family"].eq("group_concordance"), "dataset"])
    assert set(out["representation"]) == {"Gene", "Pathway"}
    assert set(out["direction"].fillna("ALL")) == {"ALL", "C2G"}
    retrieval_rows = out.loc[out["analysis_family"].eq("retrieval") & out["metric_name"].eq("mean_mrr_corrected")]
    lincs_gene_retrieval = retrieval_rows.loc[
        retrieval_rows["dataset"].eq("LINCS") & retrieval_rows["representation"].eq("Gene")
    ]
    assert lincs_gene_retrieval["fill_value"].tolist() == pytest.approx([0.5, 1.0])
    scperturb_retrieval = retrieval_rows.loc[retrieval_rows["dataset"].eq("scPerturb")]
    assert set(scperturb_retrieval["specificity_tier"].dropna()) == {
        "The Cleanest Hits",
        "The Family Hits",
    }
    assert "Control" not in set(scperturb_retrieval["specificity_tier"].dropna())
    assert scperturb_retrieval["cell_line"].eq("K562").all()
    assert scperturb_retrieval["n_total"].notna().all()


def test_build_figure3_panel_3c_builds_lincs_only_ranked_suitability_panel(tmp_path: Path) -> None:
    source_path = tmp_path / "figure3_task2_cell_line_pattern_summary.csv"
    output_path = tmp_path / "figure3_panel_3c_cell_line_pattern.csv"
    rows: list[dict[str, object]] = []
    cell_lines = ["A375", "MCF7", "PC3", "HT29", "HEPG2", "BT20"]
    value_map = {
        "A375": 0.96,
        "MCF7": 0.91,
        "PC3": 0.86,
        "HT29": 0.20,
        "HEPG2": 0.16,
        "BT20": 0.12,
        "K562": 0.55,
    }
    for dataset, dataset_cell_lines in [("LINCS", cell_lines), ("scPerturb", ["K562"])]:
        for cell_line in dataset_cell_lines:
            for analysis_family, direction, metric_name in [
                ("group_concordance", pd.NA, "cosine_centroid"),
                ("retrieval", "C2G", "mrr_corrected"),
            ]:
                for representation, offset in [("Gene", 0.03), ("Pathway", 0.0)]:
                    value = value_map[cell_line] + offset
                    rows.append(
                        {
                            "dataset": dataset,
                            "cell_line": cell_line,
                            "analysis_family": analysis_family,
                            "direction": direction,
                            "representation": representation,
                            "metric_name": metric_name,
                            "n_targets": 12,
                            "n_source_rows": 12,
                            "mean_metric_value": value,
                            "median_metric_value": value,
                            "min_metric_value": value,
                            "max_metric_value": value,
                            "n_queries_total": 80,
                            "n_queries_mean": 20,
                            "n_chem_instances_used_total": 5,
                            "n_gen_instances_used_total": 5,
                            "n_chem_sub_total": 5,
                            "n_gen_sub_total": 5,
                            "pattern_summary_scope": "test",
                            "pattern_note": "test",
                            "fm_scope_note": "test",
                        }
                    )
        rows.append(
            {
                "dataset": dataset,
                "cell_line": dataset_cell_lines[0],
                "analysis_family": "retrieval",
                "direction": "G2C",
                "representation": "Gene",
                "metric_name": "mrr_corrected",
                "n_targets": 12,
                "n_source_rows": 12,
                "mean_metric_value": 0.99,
                "median_metric_value": 0.99,
                "min_metric_value": 0.99,
                "max_metric_value": 0.99,
                "n_queries_total": 80,
                "n_queries_mean": 20,
                "n_chem_instances_used_total": 5,
                "n_gen_instances_used_total": 5,
                "n_chem_sub_total": 5,
                "n_gen_sub_total": 5,
                "pattern_summary_scope": "exclude",
                "pattern_note": "exclude",
                "fm_scope_note": "exclude",
            }
        )
    rows.append(
        {
            "dataset": "LINCS",
            "cell_line": "A375",
            "analysis_family": "group_concordance",
            "direction": pd.NA,
            "representation": "geneformer",
            "metric_name": "cosine_centroid",
            "n_targets": 12,
            "n_source_rows": 12,
            "mean_metric_value": 0.99,
            "median_metric_value": 0.99,
            "min_metric_value": 0.99,
            "max_metric_value": 0.99,
            "n_queries_total": pd.NA,
            "n_queries_mean": pd.NA,
            "n_chem_instances_used_total": 5,
            "n_gen_instances_used_total": 5,
            "n_chem_sub_total": 5,
            "n_gen_sub_total": 5,
            "pattern_summary_scope": "exclude",
            "pattern_note": "exclude",
            "fm_scope_note": "exclude",
        }
    )
    pd.DataFrame(rows).to_csv(source_path, index=False)

    mpr.build_figure3_panel_3c(source_path=source_path, output_path=output_path)

    out = pd.read_csv(output_path)
    assert set(out["dataset"]) == {"LINCS"}
    assert len(out) == 6
    assert out["cell_line"].nunique() == 6
    assert out["row_rank"].tolist() == [1, 2, 3, 4, 5, 6]
    assert out["suitability_score"].is_monotonic_decreasing
    assert "support_n" in out.columns
    assert "selection_direction" not in out.columns
    assert "fill_value" not in out.columns


def test_build_figure3_panel_3e_uses_target_level_units_and_corrected_stats(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
) -> None:
    stats_path = tmp_path / "manuscript_comparison_statistics.csv"
    output_path = tmp_path / "figure3_panel_3e_gene_vs_pathway_paired.csv"
    support_long = pd.DataFrame(
        [
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "target": "T1",
                "representation": "Gene",
                "analysis_family": "group_concordance",
                "direction": pd.NA,
                "metric_name": "cosine_centroid",
                "metric_value": 0.30,
            },
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "target": "T1",
                "representation": "Pathway",
                "analysis_family": "group_concordance",
                "direction": pd.NA,
                "metric_name": "cosine_centroid",
                "metric_value": 0.20,
            },
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "target": "T2",
                "representation": "Gene",
                "analysis_family": "group_concordance",
                "direction": pd.NA,
                "metric_name": "cosine_centroid",
                "metric_value": 0.28,
            },
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "target": "T2",
                "representation": "Pathway",
                "analysis_family": "group_concordance",
                "direction": pd.NA,
                "metric_name": "cosine_centroid",
                "metric_value": 0.19,
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target": "TP53",
                "representation": "Gene",
                "analysis_family": "retrieval",
                "direction": "C2G",
                "metric_name": "mrr_corrected",
                "metric_value": 0.42,
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target": "TP53",
                "representation": "Pathway",
                "analysis_family": "retrieval",
                "direction": "C2G",
                "metric_name": "mrr_corrected",
                "metric_value": 0.31,
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target": "HDAC1",
                "representation": "Gene",
                "analysis_family": "retrieval",
                "direction": "C2G",
                "metric_name": "mrr_corrected",
                "metric_value": 0.39,
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target": "HDAC1",
                "representation": "Pathway",
                "analysis_family": "retrieval",
                "direction": "C2G",
                "metric_name": "mrr_corrected",
                "metric_value": 0.24,
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target": "TP53",
                "representation": "Gene",
                "analysis_family": "retrieval",
                "direction": "G2C",
                "metric_name": "mrr_corrected",
                "metric_value": 0.99,
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target": "TP53",
                "representation": "Pathway",
                "analysis_family": "retrieval",
                "direction": "G2C",
                "metric_name": "mrr_corrected",
                "metric_value": 0.01,
            },
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "target": "T1",
                "representation": "geneformer",
                "analysis_family": "group_concordance",
                "direction": pd.NA,
                "metric_name": "cosine_centroid",
                "metric_value": 0.90,
            },
        ]
    )
    monkeypatch.setattr(mpr, "build_task2_target_support_long", lambda *_args, **_kwargs: support_long.copy())

    pd.DataFrame(
        [
            {
                "comparison_family": "figure3_task2_performance_target_level_representation",
                "dataset_scope": "LINCS",
                "analysis_family": "group_concordance",
                "direction_scope": "ALL",
                "metric_name": "cosine_centroid",
                "bh_q": 0.02,
                "test_status": "tested",
                "effect_direction": "group_a_gt_group_b",
                "median_delta": 0.10,
                "n_test_units": 2,
            },
            {
                "comparison_family": "figure3_task2_performance_target_level_representation",
                "dataset_scope": "scPerturb",
                "analysis_family": "retrieval",
                "direction_scope": "C2G",
                "metric_name": "mrr_corrected",
                "bh_q": 0.04,
                "test_status": "tested",
                "effect_direction": "group_a_gt_group_b",
                "median_delta": 0.11,
                "n_test_units": 2,
            },
        ]
    ).to_csv(stats_path, index=False)

    mpr.build_figure3_panel_3e(
        task2_group_root=tmp_path / "unused_group_root",
        task2_retrieval_per_query_path=tmp_path / "unused_retrieval_per_query.parquet",
        comparison_stats_path=stats_path,
        output_path=output_path,
    )

    out = pd.read_csv(output_path)
    assert set(out["unit_source"]) == {"target_level_performance"}
    assert set(out["row_kind"]) == {"raw", "summary"}
    assert set(out["representation"]) == {"Gene", "Pathway"}
    assert set(out["dataset"]) == {"LINCS", "scPerturb"}
    assert set(out["direction"].fillna("ALL")) == {"ALL", "C2G"}
    assert "G2C" not in set(out["direction"].fillna("ALL"))
    summary = out.loc[out["row_kind"].eq("summary")].copy()
    assert sorted(summary["bh_q"].dropna().unique().tolist()) == pytest.approx([0.02, 0.04])
    assert summary["n_units"].notna().all()
    assert summary["ci_low"].notna().all()
    assert summary["ci_high"].notna().all()
    scperturb_retrieval = summary.loc[
        summary["dataset"].eq("scPerturb")
        & summary["analysis_family"].eq("retrieval")
        & summary["direction"].eq("C2G")
    ]
    assert set(scperturb_retrieval["n_units"]) == {2}


def test_build_figure2_panel_2e_selects_blockwise_dual_tail_examples(tmp_path: Path) -> None:
    source_path = tmp_path / "cell_line_high_concordance.csv"
    output_path = tmp_path / "plot_ready.csv"
    rows: list[dict[str, object]] = []
    cell_lines = ["A375", "MCF7", "PC3", "HT29", "HEPG2", "BT20"]
    high_lines = {"A375", "MCF7", "PC3"}
    for dataset in ["LINCS", "scPerturb"]:
        for perturbation_type in ["Chemical", "Genetic"]:
            for cell_line in cell_lines:
                for representation, positive_shift in [("Gene", 2), ("Pathway", 0)]:
                    if cell_line in high_lines:
                        n_positive = 28 + positive_shift
                        n_negative = 12 - min(positive_shift, 1)
                        background_positive = 8
                        background_negative = 32
                    else:
                        n_positive = 6 + positive_shift
                        n_negative = 34 - min(positive_shift, 1)
                        background_positive = 20
                        background_negative = 20
                    rows.append(
                        {
                            "scope": "internal",
                            "dataset_or_direction": dataset,
                            "perturbation_type": perturbation_type,
                            "representation": representation,
                            "cell_line": cell_line,
                            "n_queries_total": 40,
                            "n_positive": n_positive,
                            "n_negative": n_negative,
                            "n_background_positive": background_positive,
                            "n_background_negative": background_negative,
                            "success_rate": n_positive / 40,
                            "background_success_rate": background_positive / (background_positive + background_negative),
                            "odds_ratio": 1.0,
                            "raw_p": 0.5,
                            "bh_q": 0.5,
                            "significant_bool": False,
                            "success_definition": "x",
                            "source_path": "/tmp/source.csv",
                            "underpowered_strata_bool": False,
                        }
                    )
    rows.extend(
        [
            {
                "scope": "internal",
                "dataset_or_direction": "LINCS",
                "perturbation_type": "Chemical",
                "representation": representation,
                "cell_line": "INTRAHEPATIC CHOLANGIOCYTE ORGANOIDS",
                "n_queries_total": 40,
                "n_positive": 35,
                "n_negative": 5,
                "n_background_positive": 5,
                "n_background_negative": 35,
                "success_rate": 0.875,
                "background_success_rate": 0.125,
                "odds_ratio": 1.0,
                "raw_p": 0.01,
                "bh_q": 0.01,
                "significant_bool": True,
                "success_definition": "x",
                "source_path": "/tmp/exclude.csv",
                "underpowered_strata_bool": False,
            }
            for representation in ["Gene", "Pathway"]
        ]
    )
    rows.extend(
        [
            {
                "scope": "internal",
                "dataset_or_direction": "LINCS",
                "perturbation_type": "Chemical",
                "representation": representation,
                "cell_line": "SKMEL5",
                "n_queries_total": 10,
                "n_positive": 9,
                "n_negative": 1,
                "n_background_positive": 1,
                "n_background_negative": 9,
                "success_rate": 0.9,
                "background_success_rate": 0.1,
                "odds_ratio": 1.0,
                "raw_p": 0.01,
                "bh_q": 0.01,
                "significant_bool": True,
                "success_definition": "x",
                "source_path": "/tmp/underpowered.csv",
                "underpowered_strata_bool": True,
            }
            for representation in ["Gene", "Pathway"]
        ]
    )
    pd.DataFrame(rows).to_csv(source_path, index=False)

    mpr.build_figure2_panel_2e(source_path=source_path, output_path=output_path)

    out = pd.read_csv(output_path)
    assert set(out["scope"]) == {"internal"}
    assert set(out["col_block"]) == {"LINCS", "scPerturb"}
    assert set(out["row_block"]) == {"Chemical", "Genetic"}
    assert "INTRAHEPATIC CHOLANGIOCYTE ORGANOIDS" not in set(out["cell_line"])
    assert "SKMEL5" not in set(out["cell_line"])
    direction_counts = (
        out.groupby(["col_block", "row_block", "selection_direction"], dropna=False)["cell_line"]
        .nunique()
        .to_dict()
    )
    assert direction_counts == {
        ("LINCS", "Chemical", "high"): 3,
        ("LINCS", "Chemical", "low"): 3,
        ("LINCS", "Genetic", "high"): 3,
        ("LINCS", "Genetic", "low"): 3,
        ("scPerturb", "Chemical", "high"): 3,
        ("scPerturb", "Chemical", "low"): 3,
        ("scPerturb", "Genetic", "high"): 3,
        ("scPerturb", "Genetic", "low"): 3,
    }
    assert out["support_n"].notna().all()
    assert out.groupby(["col_block", "row_block"])["shared_flag"].apply(lambda values: values.astype(bool).sum()).ge(2).all()
    assert set(out["perturbation_type_label"]) == {"Chemical", "Genetic"}


def test_repair_target_high_concordance_summary_rejects_multi_target_rows(tmp_path: Path) -> None:
    source_path = tmp_path / "figure2_task1_target_high_concordance_summary.csv"
    output_path = tmp_path / "figure2_task1_target_high_concordance_summary_repaired.csv"
    pd.DataFrame(
        [
            {
                "scope": "internal",
                "dataset_or_direction": "scPerturb",
                "perturbation_type": "Chemical",
                "representation": "Gene",
                "target_token": "HDAC1_HDAC2",
                "n_queries_total": 8,
                "n_positive": 4,
                "n_negative": 4,
                "n_background_positive": 10,
                "n_background_negative": 10,
                "success_rate": 0.5,
                "background_success_rate": 0.5,
                "odds_ratio": 1.0,
                "raw_p": 1.0,
                "bh_q": 1.0,
                "significant_bool": False,
                "success_definition": "x",
                "source_path": "/tmp/a.csv",
                "underpowered_strata_bool": True,
            },
            {
                "scope": "internal",
                "dataset_or_direction": "scPerturb",
                "perturbation_type": "Chemical",
                "representation": "Gene",
                "target_token": "HDAC1|HDAC2",
                "n_queries_total": 12,
                "n_positive": 9,
                "n_negative": 3,
                "n_background_positive": 5,
                "n_background_negative": 15,
                "success_rate": 0.75,
                "background_success_rate": 0.25,
                "odds_ratio": 9.0,
                "raw_p": 0.01,
                "bh_q": 0.01,
                "significant_bool": True,
                "success_definition": "x",
                "source_path": "/tmp/b.csv",
                "underpowered_strata_bool": True,
            },
        ]
    ).to_csv(source_path, index=False)

    with pytest.raises(ValueError, match="repair the analysis-layer builder first"):
        mpr.repair_target_high_concordance_summary(source_path=source_path, output_path=output_path)


def test_task1_high_concordance_main_splits_atomic_targets_upstream(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> None:
    internal_path = tmp_path / "task1_internal.parquet"
    cross_path = tmp_path / "task1_cross.parquet"
    cell_line_output_path = tmp_path / "figure2_task1_cell_line_high_concordance_summary.csv"
    target_output_path = tmp_path / "figure2_task1_target_high_concordance_summary.csv"

    pd.DataFrame(
        [
            {
                "scope": "internal",
                "dataset_or_direction": "scPerturb",
                "perturbation_type": "Chemical",
                "representation": "Gene",
                "cell_line": "K562",
                "target_token": "HDAC1_HDAC2",
                "mrr_corrected": 0.5,
                "hit10_corrected": 1.0,
            },
            {
                "scope": "internal",
                "dataset_or_direction": "scPerturb",
                "perturbation_type": "Chemical",
                "representation": "Gene",
                "cell_line": "K562",
                "target_token": "HDAC1|HDAC2",
                "mrr_corrected": 0.0,
                "hit10_corrected": 0.0,
            },
            {
                "scope": "internal",
                "dataset_or_direction": "scPerturb",
                "perturbation_type": "Chemical",
                "representation": "Gene",
                "cell_line": "K562",
                "target_token": "ABL1",
                "mrr_corrected": 0.0,
                "hit10_corrected": 0.0,
            },
        ]
    ).to_parquet(internal_path, index=False)
    pd.DataFrame(columns=mthc.REQUIRED_COLUMNS).to_parquet(cross_path, index=False)

    monkeypatch.setattr(
        mthc,
        "parse_args",
        lambda: SimpleNamespace(
            project_root=tmp_path,
            internal_per_query_path=internal_path,
            cross_per_query_path=cross_path,
            cell_line_output_path=cell_line_output_path,
            target_output_path=target_output_path,
        ),
    )

    mthc.main()

    out = pd.read_csv(target_output_path).sort_values("target_token", kind="mergesort").reset_index(drop=True)
    assert out["target_token"].tolist() == ["ABL1", "HDAC1", "HDAC2"]
    assert not out["target_token"].str.contains(r"[|_]", regex=True).any()
    hdac1 = out.loc[out["target_token"].eq("HDAC1")].iloc[0]
    hdac2 = out.loc[out["target_token"].eq("HDAC2")].iloc[0]
    assert hdac1["n_queries_total"] == 2
    assert hdac1["n_positive"] == 1
    assert hdac1["n_background_positive"] == 0
    assert hdac1["n_background_negative"] == 1
    assert hdac2["n_queries_total"] == 2
    assert hdac2["n_positive"] == 1
    assert hdac2["n_background_positive"] == 0
    assert hdac2["n_background_negative"] == 1
    assert out["raw_p"].notna().all()
    assert out["bh_q"].notna().all()
    assert out["underpowered_strata_bool"].eq(True).all()


def test_build_figure3_panel_3f_tradeoff_joins_gene_baseline_q_values(tmp_path: Path) -> None:
    target_pattern_path = tmp_path / "figure3_task2_target_pattern_summary.csv"
    comparison_stats_path = tmp_path / "manuscript_comparison_statistics.csv"
    figure3_scope_path = tmp_path / "figure3_task2_scope_summary.csv"
    output_path = tmp_path / "figure3_panel_3f_fm_local_tradeoff.csv"

    pd.DataFrame(
        [
            {
                "dataset": "scPerturb",
                "target": "T1",
                "analysis_family": "group_concordance",
                "direction": pd.NA,
                "representation": "Gene",
                "metric_name": "cosine_centroid",
                "n_cell_lines": 1,
                "n_source_rows": 1,
                "mean_metric_value": 0.8,
                "n_queries_total": pd.NA,
                "fm_scope_note": "base",
            },
            {
                "dataset": "scPerturb",
                "target": "T2",
                "analysis_family": "group_concordance",
                "direction": pd.NA,
                "representation": "Gene",
                "metric_name": "cosine_centroid",
                "n_cell_lines": 1,
                "n_source_rows": 1,
                "mean_metric_value": 0.6,
                "n_queries_total": pd.NA,
                "fm_scope_note": "base",
            },
            {
                "dataset": "scPerturb",
                "target": "T1",
                "analysis_family": "group_concordance",
                "direction": pd.NA,
                "representation": "Pathway",
                "metric_name": "cosine_centroid",
                "n_cell_lines": 1,
                "n_source_rows": 1,
                "mean_metric_value": 0.5,
                "n_queries_total": pd.NA,
                "fm_scope_note": "base",
            },
            {
                "dataset": "scPerturb",
                "target": "T1",
                "analysis_family": "group_concordance",
                "direction": pd.NA,
                "representation": "geneformer",
                "metric_name": "cosine_centroid",
                "n_cell_lines": 1,
                "n_source_rows": 1,
                "mean_metric_value": 0.4,
                "n_queries_total": pd.NA,
                "fm_scope_note": "fm",
            },
            {
                "dataset": "scPerturb",
                "target": "T1",
                "analysis_family": "retrieval",
                "direction": "C2G",
                "representation": "Gene",
                "metric_name": "mrr_corrected",
                "n_cell_lines": 1,
                "n_source_rows": 1,
                "mean_metric_value": 0.7,
                "n_queries_total": 10,
                "fm_scope_note": "base",
            },
            {
                "dataset": "scPerturb",
                "target": "T2",
                "analysis_family": "retrieval",
                "direction": "C2G",
                "representation": "Gene",
                "metric_name": "mrr_corrected",
                "n_cell_lines": 1,
                "n_source_rows": 1,
                "mean_metric_value": 0.5,
                "n_queries_total": 12,
                "fm_scope_note": "base",
            },
            {
                "dataset": "scPerturb",
                "target": "T1",
                "analysis_family": "retrieval",
                "direction": "C2G",
                "representation": "Pathway",
                "metric_name": "mrr_corrected",
                "n_cell_lines": 1,
                "n_source_rows": 1,
                "mean_metric_value": 0.45,
                "n_queries_total": 10,
                "fm_scope_note": "base",
            },
            {
                "dataset": "scPerturb",
                "target": "T1",
                "analysis_family": "retrieval",
                "direction": "C2G",
                "representation": "geneformer",
                "metric_name": "mrr_corrected",
                "n_cell_lines": 1,
                "n_source_rows": 1,
                "mean_metric_value": 0.3,
                "n_queries_total": 10,
                "fm_scope_note": "fm",
            },
        ]
    ).to_csv(target_pattern_path, index=False)

    pd.DataFrame(
        [
            {"dataset": "scPerturb", "cell_line": "K562", "representation": "Gene", "availability_status": "available"},
            {"dataset": "scPerturb", "cell_line": "K562", "representation": "Pathway", "availability_status": "available"},
            {"dataset": "scPerturb", "cell_line": "K562", "representation": "geneformer", "availability_status": "available"},
        ]
    ).to_csv(figure3_scope_path, index=False)

    pd.DataFrame(
        [
            {
                "comparison_family": "figure3_task2_target_pattern_representation",
                "analysis_family": "group_concordance",
                "metric_name": "cosine_centroid",
                "group_a": "Gene",
                "group_b": "Pathway",
                "bh_q": 0.04,
                "test_status": "tested",
                "dataset_scope": "scPerturb",
                "direction_scope": "ALL",
            },
            {
                "comparison_family": "figure3_task2_target_pattern_representation",
                "analysis_family": "retrieval",
                "metric_name": "mrr_corrected",
                "group_a": "Gene",
                "group_b": "Pathway",
                "bh_q": 0.08,
                "test_status": "tested",
                "dataset_scope": "scPerturb",
                "direction_scope": "C2G",
            },
            {
                "comparison_family": "figure3_task2_k562_fm_target_anchored",
                "analysis_family": "group_concordance",
                "metric_name": "cosine_centroid",
                "group_a": "geneformer",
                "group_b": "Gene",
                "bh_q": 0.002,
                "test_status": "tested",
                "dataset_scope": "ALL",
                "direction_scope": "ALL",
            },
            {
                "comparison_family": "figure3_task2_k562_fm_c2g_query_anchored",
                "analysis_family": "retrieval",
                "metric_name": "mrr_corrected",
                "group_a": "geneformer",
                "group_b": "Gene",
                "bh_q": 0.003,
                "test_status": "tested",
                "dataset_scope": "ALL",
                "direction_scope": "C2G",
            },
        ]
    ).to_csv(comparison_stats_path, index=False)

    mpr.build_figure3_panel_3f_tradeoff(
        target_pattern_path=target_pattern_path,
        performance_path=tmp_path / "unused_performance.csv",
        direction_support_path=tmp_path / "unused_direction.csv",
        comparison_stats_path=comparison_stats_path,
        figure3_scope_path=figure3_scope_path,
        task2_retrieval_per_query_path=tmp_path / "unused.parquet",
        output_path=output_path,
    )

    out = pd.read_csv(output_path)
    assert set(out["representation"]) == {"Gene", "Pathway", "geneformer"}
    assert set(out["dataset"]) == {"scPerturb"}
    assert set(out["cell_line"]) == {"K562"}
    gene_row = out.loc[(out["representation"].eq("Gene")) & (out["metric_name"].eq("cosine_centroid"))].iloc[0]
    pathway_row = out.loc[(out["representation"].eq("Pathway")) & (out["metric_name"].eq("cosine_centroid"))].iloc[0]
    fm_row = out.loc[(out["representation"].eq("geneformer")) & (out["metric_name"].eq("mrr_corrected"))].iloc[0]
    assert gene_row["is_gene_baseline"] == True
    assert pd.isna(gene_row["bh_q_vs_gene"])
    assert pathway_row["bh_q_vs_gene"] == pytest.approx(0.04)
    assert fm_row["bh_q_vs_gene"] == pytest.approx(0.003)
    assert fm_row["comparison_source"] == "figure3_task2_k562_fm_c2g_query_anchored"
    assert gene_row["gene_reference_metric_value"] == pytest.approx(gene_row["metric_value"])
    assert pathway_row["gene_reference_metric_value"] == pytest.approx(gene_row["metric_value"])


def test_build_figure3_panel_3d_builds_ranked_target_suitability_panel(tmp_path: Path) -> None:
    source_path = tmp_path / "target_pattern.csv"
    output_path = tmp_path / "figure3_panel_3d.csv"
    rows: list[dict[str, object]] = []
    targets = ["T1", "T2", "T3", "T4", "T5", "T6"]
    value_map = {
        "T1": 0.95,
        "T2": 0.90,
        "T3": 0.85,
        "T4": 0.20,
        "T5": 0.15,
        "T6": 0.10,
    }
    for dataset in ["LINCS", "scPerturb"]:
        for target in targets:
            for analysis_family, direction, metric_name in [
                ("group_concordance", pd.NA, "cosine_centroid"),
                ("retrieval", "C2G", "mrr_corrected"),
            ]:
                for representation, offset in [("Gene", 0.03), ("Pathway", 0.0)]:
                    value = value_map[target] + offset
                    rows.append(
                        {
                            "dataset": dataset,
                            "target": target,
                            "analysis_family": analysis_family,
                            "direction": direction,
                            "representation": representation,
                            "metric_name": metric_name,
                            "n_cell_lines": 4,
                            "n_source_rows": 4,
                            "mean_metric_value": value,
                            "median_metric_value": value,
                            "min_metric_value": value,
                            "max_metric_value": value,
                            "n_queries_total": 40,
                            "n_queries_mean": 20,
                            "n_chem_instances_used_total": 2,
                            "n_gen_instances_used_total": 2,
                            "n_chem_sub_total": 2,
                            "n_gen_sub_total": 2,
                            "pattern_summary_scope": "test",
                            "pattern_note": "test",
                            "fm_scope_note": "test",
                        }
                    )
        rows.append(
            {
                "dataset": dataset,
                "target": "DROP_ME",
                "analysis_family": "retrieval",
                "direction": "G2C",
                "representation": "Gene",
                "metric_name": "mrr_corrected",
                "n_cell_lines": 4,
                "n_source_rows": 4,
                "mean_metric_value": 0.99,
                "median_metric_value": 0.99,
                "min_metric_value": 0.99,
                "max_metric_value": 0.99,
                "n_queries_total": 40,
                "n_queries_mean": 20,
                "n_chem_instances_used_total": 2,
                "n_gen_instances_used_total": 2,
                "n_chem_sub_total": 2,
                "n_gen_sub_total": 2,
                "pattern_summary_scope": "exclude",
                "pattern_note": "exclude",
                "fm_scope_note": "exclude",
            }
        )
    rows.append(
        {
            "dataset": "scPerturb",
            "target": "T1",
            "analysis_family": "group_concordance",
            "direction": pd.NA,
            "representation": "geneformer",
            "metric_name": "cosine_centroid",
            "n_cell_lines": 4,
            "n_source_rows": 4,
            "mean_metric_value": 0.99,
            "median_metric_value": 0.99,
            "min_metric_value": 0.99,
            "max_metric_value": 0.99,
            "n_queries_total": pd.NA,
            "n_queries_mean": pd.NA,
            "n_chem_instances_used_total": 2,
            "n_gen_instances_used_total": 2,
            "n_chem_sub_total": 2,
            "n_gen_sub_total": 2,
            "pattern_summary_scope": "exclude",
            "pattern_note": "exclude",
            "fm_scope_note": "exclude",
        }
    )
    pd.DataFrame(rows).to_csv(source_path, index=False)

    mpr.build_figure3_panel_3d(source_path=source_path, output_path=output_path)

    out = pd.read_csv(output_path)
    assert set(out["dataset"]) == {"LINCS", "scPerturb"}
    assert len(out) == 12
    assert out.groupby("dataset")["target"].nunique().to_dict() == {"LINCS": 6, "scPerturb": 6}
    assert out.groupby("dataset")["row_rank"].apply(list).to_dict() == {
        "LINCS": [1, 2, 3, 4, 5, 6],
        "scPerturb": [1, 2, 3, 4, 5, 6],
    }
    assert out.groupby("dataset")["suitability_score"].apply(lambda values: values.is_monotonic_decreasing).all()
    assert "selection_direction" not in out.columns
    assert out["shared_across_modalities_bool"].astype(bool).all()


def test_build_figure3_panel_3j_collapses_repeated_perturbation_type_rows(tmp_path: Path) -> None:
    source_path = tmp_path / "contextual_support.csv"
    output_path = tmp_path / "figure3_panel_3j.csv"
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
                "metric_value": 0.2,
            },
            {
                "dataset": "LINCS",
                "cell_line": "A375",
                "target": "ABCB1",
                "perturbation_type": "Genetic",
                "representation": "Pathway",
                "analysis_family": "retrieval",
                "metric_name": "mrr_corrected",
                "metric_value": 0.0,
            },
        ]
    ).to_csv(source_path, index=False)

    mpr.build_figure3_panel_3j(source_path=source_path, output_path=output_path)

    out = pd.read_csv(output_path)
    assert len(out) == 2
    gene_row = out.loc[out["representation"].eq("Gene")].iloc[0]
    pathway_row = out.loc[out["representation"].eq("Pathway")].iloc[0]
    assert gene_row["mean_metric_value"] == pytest.approx(0.6)
    assert pathway_row["mean_metric_value"] == pytest.approx(0.1)
    assert gene_row["n_triplets"] == 1
    assert set(out["reference_only_note"]) == {"background/reference only; not bridge/ceiling."}


def test_build_figure1_object_map_includes_support_and_historical_rows(tmp_path: Path) -> None:
    manifest_path = tmp_path / "framework_analysis_manifest.json"
    output_path = tmp_path / "figure1_panel_1e_result_object_map.csv"
    manifest = {
        "canonical_objects": [
            {
                "object_name": "Task1 scope summary",
                "canonical_filename": "figure2_task1_scope_summary.csv",
                "justification": "canonical scope summary",
            }
        ],
        "historical_noncanonical_objects": [
            {
                "object_name": "Legacy bridge",
                "path": "/tmp/task1_task2_group_bridge.csv",
                "reason_no_longer_canonical": "historical only",
                "disposition": "historical_only",
            }
        ],
    }
    manifest_path.write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")

    mpr.build_figure1_object_map(manifest_path=manifest_path, output_path=output_path)

    out = pd.read_csv(output_path)
    assert "canonical" in set(out["object_status"])
    assert "support_only" in set(out["object_status"])
    assert "historical_only" in set(out["object_status"])
    canonical = out.loc[out["canonical_filename"].eq("figure2_task1_scope_summary.csv")].iloc[0]
    assert canonical["panel_target"] == "2A"


def test_build_figure1_scope_matrix_contains_cross_chemical_policy_row(tmp_path: Path) -> None:
    output_path = tmp_path / "figure1_panel_1c_lawful_scope_matrix.csv"

    mpr.build_figure1_scope_matrix(output_path=output_path)

    out = pd.read_csv(output_path)
    row = out.loc[
        out["task"].eq("Task1")
        & out["scope"].eq("cross")
        & out["dataset_or_direction"].eq("LINCS_to_scPerturb")
        & out["perturbation_type"].eq("Chemical")
        & out["representation_class"].eq("Gene")
    ].iloc[0]
    assert row["status"] == "excluded_by_support_gate"
    assert "must not be framed as attrition" in row["scope_note"]


def test_spatial_revision_targets_follow_requested_order() -> None:
    assert mpr.SPATIAL_REVISION_PLOT_READY_OUTPUTS == [
        "figure2_panel_2a_task1_scope.csv",
        "figure2_panel_2b_internal_performance_overview.csv",
        "figure2_panel_2c_gene_vs_pathway_matched_units.csv",
        "figure2_panel_2d_internal_to_cross_degradation.csv",
        "figure2_panel_2e_cell_line_high_concordance_summary.csv",
        "figure2_panel_2f_target_high_concordance_summary.csv",
        "figure3_panel_3b_c2g_performance_overview.csv",
        "figure3_panel_3c_cell_line_pattern.csv",
        "figure3_panel_3d_target_pattern_summary.csv",
        "figure3_panel_3e_gene_vs_pathway_paired.csv",
        "figure3_panel_3f_fm_local_tradeoff.csv",
        "figure1_panel_1b_dataset_context_coverage.csv",
    ]


def test_build_figure1_coverage_keeps_leaf_level_support_columns(tmp_path: Path) -> None:
    task1_inventory_path = tmp_path / "task1_inventory.csv"
    task2_pairs_coverage_path = tmp_path / "task2_pairs_coverage.csv"
    output_path = tmp_path / "figure1_panel_1b_dataset_context_coverage.csv"

    pd.DataFrame(
        [
            {
                "dataset": "LINCS",
                "perturbation_type": "Chemical",
                "cell_line": "MCF7",
                "target_token": "HDAC1_HDAC2",
                "n_instances": 5,
            }
        ]
    ).to_csv(task1_inventory_path, index=False)
    pd.DataFrame(
        [
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target_token": "HDAC1|HDAC2",
                "n_chem_instances": 3,
                "n_gen_instances": 5,
                "is_eligible_bool": True,
            },
            {
                "dataset": "scPerturb",
                "cell_line": "K562",
                "target_token": "ABL1",
                "n_chem_instances": 10,
                "n_gen_instances": 10,
                "is_eligible_bool": False,
            },
        ]
    ).to_csv(task2_pairs_coverage_path, index=False)

    mpr.build_figure1_coverage(
        task1_inventory_path=task1_inventory_path,
        task2_pairs_coverage_path=task2_pairs_coverage_path,
        output_path=output_path,
    )

    out = pd.read_csv(output_path)
    assert len(out) == 3
    task1_row = out.loc[out["task"].eq("Task1")].iloc[0]
    task2_rows = out.loc[out["task"].eq("Task2")].sort_values("perturbation_type", kind="mergesort").reset_index(drop=True)
    assert task1_row["context_label"] == "MCF7"
    assert task1_row["leaf_weight"] == pytest.approx(5)
    assert task1_row["context_summary_label"] == "MCF7\n1 targets"
    assert task1_row["source_table"] == str(task1_inventory_path)
    assert task2_rows["leaf_weight"].tolist() == pytest.approx([3, 5])
    assert set(task2_rows["perturbation_type"]) == {"Chemical", "Genetic"}
    assert set(task2_rows["n_chem_instances"]) == {3}
    assert set(task2_rows["n_gen_instances"]) == {5}
    assert set(task2_rows["source_table"]) == {str(task2_pairs_coverage_path)}
    assert set(task2_rows["context_label"]) == {"K562"}
    assert set(task2_rows["context_note"]) == {"1 targets"}


def test_build_figure2_panel_2f_rejects_multi_target_rows_before_plotting(tmp_path: Path) -> None:
    source_path = tmp_path / "target_high_concordance.csv"
    output_path = tmp_path / "figure2_panel_2f_target_high_concordance_summary.csv"

    pd.DataFrame(
        [
            {
                "scope": "internal",
                "dataset_or_direction": "scPerturb",
                "perturbation_type": "Chemical",
                "representation": "Gene",
                "target_token": "HDAC1_HDAC2",
                "n_queries_total": 40,
                "n_positive": 20,
                "n_negative": 20,
                "n_background_positive": 10,
                "n_background_negative": 30,
                "success_rate": 0.5,
                "background_success_rate": 0.25,
                "odds_ratio": 3.0,
                "raw_p": 0.02,
                "bh_q": 0.03,
                "significant_bool": True,
                "success_definition": "x",
                "source_path": "/tmp/gene_a.csv",
                "underpowered_strata_bool": False,
            },
            {
                "scope": "internal",
                "dataset_or_direction": "scPerturb",
                "perturbation_type": "Chemical",
                "representation": "Gene",
                "target_token": "HDAC1|HDAC2",
                "n_queries_total": 20,
                "n_positive": 10,
                "n_negative": 10,
                "n_background_positive": 5,
                "n_background_negative": 15,
                "success_rate": 0.5,
                "background_success_rate": 0.25,
                "odds_ratio": 3.0,
                "raw_p": 0.02,
                "bh_q": 0.03,
                "significant_bool": True,
                "success_definition": "x",
                "source_path": "/tmp/gene_b.csv",
                "underpowered_strata_bool": False,
            },
            {
                "scope": "internal",
                "dataset_or_direction": "scPerturb",
                "perturbation_type": "Chemical",
                "representation": "Pathway",
                "target_token": "HDAC1_HDAC2",
                "n_queries_total": 30,
                "n_positive": 15,
                "n_negative": 15,
                "n_background_positive": 10,
                "n_background_negative": 20,
                "success_rate": 0.5,
                "background_success_rate": 0.3333333333,
                "odds_ratio": 2.0,
                "raw_p": 0.04,
                "bh_q": 0.05,
                "significant_bool": True,
                "success_definition": "x",
                "source_path": "/tmp/pathway.csv",
                "underpowered_strata_bool": False,
            },
            {
                "scope": "internal",
                "dataset_or_direction": "scPerturb",
                "perturbation_type": "Chemical",
                "representation": "Gene",
                "target_token": "ABL1",
                "n_queries_total": 25,
                "n_positive": 5,
                "n_negative": 20,
                "n_background_positive": 25,
                "n_background_negative": 25,
                "success_rate": 0.2,
                "background_success_rate": 0.5,
                "odds_ratio": 0.25,
                "raw_p": 0.2,
                "bh_q": 0.2,
                "significant_bool": False,
                "success_definition": "x",
                "source_path": "/tmp/abl1.csv",
                "underpowered_strata_bool": False,
            },
        ]
    ).to_csv(source_path, index=False)

    with pytest.raises(ValueError, match="atomic target_token rows"):
        mpr.build_figure2_panel_2f(source_path=source_path, output_path=output_path)


def write_minimal_task2_cell_line_pattern(path: Path) -> None:
    rows = []
    scores = {
        ("LINCS", "A375"): {"Gene": 0.95, "Pathway": 0.72},
        ("LINCS", "MCF7"): {"Gene": 0.66, "Pathway": 0.88},
        ("LINCS", "PC3"): {"Gene": 0.25, "Pathway": 0.18},
        ("scPerturb", "K562"): {"Gene": 0.91, "Pathway": 0.63},
        ("scPerturb", "RPE1"): {"Gene": 0.32, "Pathway": 0.54},
    }
    for (dataset, cell_line), rep_scores in scores.items():
        for representation, value in rep_scores.items():
            rows.append(
                {
                    "dataset": dataset,
                    "cell_line": cell_line,
                    "analysis_family": "group_concordance",
                    "direction": "ALL",
                    "representation": representation,
                    "metric_name": "cosine_centroid",
                    "n_targets": 12 if dataset == "LINCS" else 8,
                    "n_source_rows": 12 if dataset == "LINCS" else 8,
                    "mean_metric_value": value,
                    "pattern_summary_scope": "cell_line_over_canonical_target_aggregation",
                    "pattern_note": "test",
                    "fm_scope_note": "FM local only.",
                    "n_queries_total": 40 if dataset == "LINCS" else 16,
                }
            )
    pd.DataFrame(rows).to_csv(path, index=False)


def write_minimal_task2_target_pattern(path: Path) -> None:
    rows = []
    scores = {
        ("LINCS", "MTOR"): {"Gene": 0.99, "Pathway": 0.98},
        ("LINCS", "EGFR"): {"Gene": 0.97, "Pathway": 0.24},
        ("LINCS", "MAP2K1"): {"Gene": 0.31, "Pathway": 0.95},
        ("LINCS", "BRD4"): {"Gene": 0.94, "Pathway": 0.92},
        ("LINCS", "LINCSONLY"): {"Gene": 0.96, "Pathway": 0.93},
        ("LINCS", "LOWLINCS"): {"Gene": 0.12, "Pathway": 0.14},
        ("scPerturb", "MTOR"): {"Gene": 0.98, "Pathway": 0.97},
        ("scPerturb", "EGFR"): {"Gene": 0.95, "Pathway": 0.21},
        ("scPerturb", "MAP2K1"): {"Gene": 0.26, "Pathway": 0.94},
        ("scPerturb", "BRD4"): {"Gene": 0.18, "Pathway": 0.16},
        ("scPerturb", "SCPONLY"): {"Gene": 0.97, "Pathway": 0.91},
        ("scPerturb", "LOWSCP"): {"Gene": 0.11, "Pathway": 0.13},
    }
    for (dataset, target), rep_scores in scores.items():
        for representation, value in rep_scores.items():
            rows.append(
                {
                    "dataset": dataset,
                    "target": target,
                    "analysis_family": "group_concordance",
                    "direction": "ALL",
                    "representation": representation,
                    "metric_name": "cosine_centroid",
                    "n_cell_lines": 6 if dataset == "LINCS" else 3,
                    "n_source_rows": 6 if dataset == "LINCS" else 3,
                    "mean_metric_value": value,
                    "median_metric_value": value,
                    "min_metric_value": value,
                    "max_metric_value": value,
                    "n_queries_total": 60 if dataset == "LINCS" else 18,
                    "n_queries_mean": 60 if dataset == "LINCS" else 18,
                    "n_chem_instances_used_total": 20 if dataset == "LINCS" else 6,
                    "n_gen_instances_used_total": 20 if dataset == "LINCS" else 6,
                    "n_chem_sub_total": 10 if dataset == "LINCS" else 3,
                    "n_gen_sub_total": 10 if dataset == "LINCS" else 3,
                    "pattern_summary_scope": "target_over_cell_line_aggregation",
                    "pattern_note": "test",
                    "fm_scope_note": "FM local only.",
                }
            )
    pd.DataFrame(rows).to_csv(path, index=False)


def test_build_extended_figure9_support_diagnostic_combines_cell_line_and_target_surfaces(tmp_path: Path) -> None:
    cell_path = tmp_path / "figure3_task2_cell_line_pattern_summary.csv"
    target_path = tmp_path / "figure3_task2_target_pattern_summary.csv"
    output_path = tmp_path / "extended_figure9_support_vs_suitability.csv"
    write_minimal_task2_cell_line_pattern(cell_path)
    write_minimal_task2_target_pattern(target_path)

    mpr.build_extended_figure9_support_diagnostic(
        cell_line_source_path=cell_path,
        target_source_path=target_path,
        output_path=output_path,
    )

    out = pd.read_csv(output_path)
    assert set(out["surface_type"]) == {"cell_line", "target"}
    assert {"support_log10", "representation_gap", "rank_within_dataset_surface"} <= set(out.columns)
    assert out["support_log10"].notna().all()
    assert out["shared_label"].isin({"shared", "dataset_specific"}).all()


def test_build_extended_figure9_target_exemplars_selects_unique_target_categories(tmp_path: Path) -> None:
    target_path = tmp_path / "figure3_task2_target_pattern_summary.csv"
    output_path = tmp_path / "extended_figure9_target_exemplars.csv"
    write_minimal_task2_target_pattern(target_path)

    mpr.build_extended_figure9_target_exemplars(
        source_path=target_path,
        output_path=output_path,
    )

    out = pd.read_csv(output_path)
    selected_targets = out[["selection_category", "target"]].drop_duplicates()
    assert selected_targets["selection_category"].nunique() >= 5
    assert selected_targets["target"].nunique() == selected_targets.shape[0]
    assert set(out["representation"]) == {"Gene", "Pathway"}
    assert out["row_display_label"].str.contains(r"^\d+\. ", regex=True).all()
