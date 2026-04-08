from __future__ import annotations

import pandas as pd
import pytest

import scripts.manuscript_task2_c2g_modifiers as mtm


def make_query_frame() -> pd.DataFrame:
    rows: list[dict[str, object]] = []

    def add_row(
        *,
        target: str,
        representation: str,
        query_time: float,
        query_dose_value: float,
        mrr: float,
        dataset: str = "scPerturb",
        cell_line: str = "K562",
    ) -> None:
        rows.append(
            {
                "dataset": dataset,
                "cell_line": cell_line,
                "representation": representation,
                "direction": "C2G",
                "query_target_tokens": target,
                "query_time": query_time,
                "query_dose_value": query_dose_value,
                "query_n_targets": 1,
                "mrr_corrected": mrr,
                "hit1_corrected": mrr * 0.5,
                "hit5_corrected": mrr * 0.75,
                "hit10_corrected": mrr * 0.9,
            }
        )

    for time, dose, mrr in [(6.0, 0.1, 0.20), (6.0, 1.0, 0.25), (24.0, 10.0, 0.35)]:
        add_row(target="TP53", representation="Gene", query_time=time, query_dose_value=dose, mrr=mrr)

    for time, dose, mrr in [(6.0, 0.1, 0.30), (24.0, 1.0, 0.32), (48.0, 10.0, 0.50)]:
        add_row(target="HDAC1", representation="Pathway", query_time=time, query_dose_value=dose, mrr=mrr)

    for time, dose, mrr in [(6.0, 0.1, 0.40), (6.0, 1.0, 0.41), (6.0, 10.0, 0.42)]:
        add_row(target="NOVAR", representation="Gene", query_time=time, query_dose_value=dose, mrr=mrr)

    for time, dose, mrr in [(6.0, 0.1, 0.15), (24.0, 1.0, 0.16)]:
        add_row(target="TOOFEW", representation="Gene", query_time=time, query_dose_value=dose, mrr=mrr)

    for time, dose, mrr in [(6.0, 0.1, 0.22), (24.0, 1.0, 0.24), (48.0, 10.0, 0.26)]:
        add_row(
            target="FM1",
            representation="FM:geneformer",
            query_time=time,
            query_dose_value=dose,
            mrr=mrr,
        )

    return pd.DataFrame(rows)


def test_spearman_rows_for_gene_and_pathway_targets(monkeypatch) -> None:
    monkeypatch.setattr(mtm, "MIN_OBS_PER_TARGET", 3)
    monkeypatch.setattr(mtm, "MIN_MODIFIER_UNIQUE", 2)

    stats = mtm.build_time_effect_stats(make_query_frame())
    mrr_rows = stats.loc[stats["metric"].eq("mrr_corrected")]

    assert set(mrr_rows["representation"]) <= {"Gene", "Pathway"}
    assert set(mrr_rows["target"]) >= {"TP53", "HDAC1"}
    assert "FM1" not in set(mrr_rows["target"])


def test_n_obs_counts_observations_not_unique_levels(monkeypatch) -> None:
    monkeypatch.setattr(mtm, "MIN_OBS_PER_TARGET", 3)
    monkeypatch.setattr(mtm, "MIN_MODIFIER_UNIQUE", 2)

    stats = mtm.build_time_effect_stats(make_query_frame())
    row = stats.loc[
        (stats["target"] == "TP53")
        & (stats["metric"] == "mrr_corrected")
        & (stats["modifier_type"] == "time")
    ].iloc[0]

    assert int(row["n_obs"]) == 3
    assert bool(row["testable_bool"]) is True


def test_not_testable_rows_emit_reasons(monkeypatch) -> None:
    monkeypatch.setattr(mtm, "MIN_OBS_PER_TARGET", 3)
    monkeypatch.setattr(mtm, "MIN_MODIFIER_UNIQUE", 2)

    stats = mtm.build_time_effect_stats(make_query_frame())

    no_var = stats.loc[
        (stats["target"] == "NOVAR")
        & (stats["metric"] == "mrr_corrected")
        & (stats["modifier_type"] == "time")
    ].iloc[0]
    assert bool(no_var["testable_bool"]) is False
    assert "no_modifier_variation" in str(no_var["untestable_reason"])

    too_few = stats.loc[
        (stats["target"] == "TOOFEW")
        & (stats["metric"] == "mrr_corrected")
        & (stats["modifier_type"] == "time")
    ].iloc[0]
    assert bool(too_few["testable_bool"]) is False
    assert "insufficient_observations" in str(too_few["untestable_reason"])


def test_bh_correction_is_representation_aware() -> None:
    df = pd.DataFrame(
        [
            {
                "representation": "Gene",
                "metric": "mrr_corrected",
                "modifier_type": "time",
                "testable_bool": True,
                "p_value": 0.01,
            },
            {
                "representation": "Gene",
                "metric": "mrr_corrected",
                "modifier_type": "time",
                "testable_bool": True,
                "p_value": 0.04,
            },
            {
                "representation": "Pathway",
                "metric": "mrr_corrected",
                "modifier_type": "time",
                "testable_bool": True,
                "p_value": 0.01,
            },
            {
                "representation": "Pathway",
                "metric": "mrr_corrected",
                "modifier_type": "time",
                "testable_bool": True,
                "p_value": 0.04,
            },
        ]
    )

    out = mtm._apply_spearman_bh(df)
    gene_row = out.loc[(out["representation"] == "Gene") & (out["p_value"] == 0.01)].iloc[0]
    pathway_row = out.loc[(out["representation"] == "Pathway") & (out["p_value"] == 0.01)].iloc[0]

    assert gene_row["bh_q"] == pytest.approx(0.02)
    assert pathway_row["bh_q"] == pytest.approx(0.02)


def test_unified_table_contains_both_modifier_types(monkeypatch) -> None:
    monkeypatch.setattr(mtm, "MIN_OBS_PER_TARGET", 3)
    monkeypatch.setattr(mtm, "MIN_MODIFIER_UNIQUE", 2)

    df = make_query_frame()
    time_stats = mtm.build_time_effect_stats(df)
    dose_stats = mtm.build_dose_effect_stats(df)
    combined = mtm.build_dose_time_spearman_stats(time_stats, dose_stats)

    assert set(combined["modifier_type"]) == {"time", "dose"}
