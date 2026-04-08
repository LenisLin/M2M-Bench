from __future__ import annotations

import json

import numpy as np
import pandas as pd
import torch

import scripts.s3_build_task2_multisource_snapshot as corrected_s3
import scripts.s3_build_task2_snapshot as legacy_s3
import scripts.s4_task2_group_concordance_multisource as task2_group


def test_corrected_scperturb_tokenization_accepts_canonical_and_osmosis_delimiters() -> None:
    assert corrected_s3.tokenize_scperturb_target("JAK1;JAK2") == ("JAK1", "JAK2")
    assert corrected_s3.tokenize_scperturb_target("JAK1_JAK2") == ("JAK1", "JAK2")
    assert corrected_s3.tokenize_scperturb_target("JAK1_JAK2;JAK3") == ("JAK1", "JAK2", "JAK3")


def test_lincs_delta_meta_preserves_full_target_tokens_while_membership_stays_eligible() -> None:
    supported_rows = pd.DataFrame(
        [
            {
                "source_row_index": 0,
                "treated_cell_id": "chem-1",
                "cell_line_raw": "A549",
                "cell_line": "A549",
                "perturbation_class": "Chemical",
                "target_raw": "A|B",
                "target_tokens_all": ("A", "B"),
                "time": 24.0,
                "dose_value": 1.0,
                "pert_type_raw": "drug",
                "paired_control_id": "ctrl-1",
                "sig_id": "sig-1",
            },
            {
                "source_row_index": 1,
                "treated_cell_id": "gen-1",
                "cell_line_raw": "A549",
                "cell_line": "A549",
                "perturbation_class": "Genetic",
                "target_raw": "A",
                "target_tokens_all": ("A",),
                "time": 24.0,
                "dose_value": 1.0,
                "pert_type_raw": "crispr",
                "paired_control_id": "ctrl-2",
                "sig_id": "sig-2",
            },
        ]
    )

    delta_meta, membership, coverage, eligible_pairs, exclusion_summary, summary = corrected_s3.build_lincs_task2_tables(
        supported_rows=supported_rows,
        excluded_rows=pd.DataFrame(),
    )

    assert delta_meta.loc[delta_meta["treated_cell_id"] == "chem-1", "target_tokens"].iloc[0] == "A;B"
    chem_membership = membership.loc[membership["treated_cell_id"] == "chem-1", "target_token"].tolist()
    assert chem_membership == ["A"]

    eligible_row = coverage.loc[coverage["target_token"] == "A"].iloc[0]
    assert bool(eligible_row["is_eligible_bool"]) is True
    ineligible_row = coverage.loc[coverage["target_token"] == "B"].iloc[0]
    assert bool(ineligible_row["is_eligible_bool"]) is False

    assert eligible_pairs["target_token"].tolist() == ["A"]
    assert exclusion_summary.empty
    assert int(summary["n_rows_included"]) == 2


def test_legacy_scperturb_drug_tokenization_supports_underscore_without_splitting_crispr_like_ids() -> None:
    assert legacy_s3.tokenize_target("JAK1_JAK2_JAK3", delimiters=("_", ";")) == ("JAK1", "JAK2", "JAK3")
    assert legacy_s3.tokenize_target("TSAI_HEK293_4", delimiters=(";",)) == ("TSAI_HEK293_4",)


def test_corrected_lincs_cell_line_normalization_preserves_distinct_labels() -> None:
    assert corrected_s3.normalize_lincs_cell_line("1HAE") == "1HAE"
    assert corrected_s3.normalize_lincs_cell_line("HA1E") == "HA1E"


def test_corrected_lincs_pathway_subset_is_projected_from_gene_delta(tmp_path) -> None:
    gene_tensor = torch.zeros((2, 2477), dtype=torch.float32)
    gene_tensor[0, 0] = 1.0
    gene_tensor[0, 1] = 2.0
    gene_tensor[1, 0] = -1.0
    gene_tensor[1, 1] = 3.0

    w = np.zeros((2477, 50), dtype=np.float32)
    w[0, 0] = 10.0
    w[1, 0] = 1.0
    w[1, 1] = 5.0

    w_path = tmp_path / "hallmark-w-2477x50.npy"
    np.save(w_path, w)
    policy_path = tmp_path / "lincs-pathway-policy.json"
    policy_path.write_text(
        json.dumps(
            {
                "mode": "project_on_load",
                "W_sha256": corrected_s3.compute_sha256(w_path),
            }
        ),
        encoding="utf-8",
    )

    loaded_w = corrected_s3.load_lincs_pathway_projection_contract(
        w_path=w_path,
        policy_path=policy_path,
    )
    out_path = tmp_path / "pathway_delta.npy"
    stats = corrected_s3.materialize_lincs_subset_pathway_from_gene_projection(
        source_gene_tensor=gene_tensor,
        source_row_indices=np.array([1, 0], dtype=np.int64),
        pathway_w=loaded_w,
        output_path=out_path,
    )

    out = np.load(out_path)
    assert stats == (2, 50)
    np.testing.assert_allclose(out[0], gene_tensor[1].numpy() @ w)
    np.testing.assert_allclose(out[1], gene_tensor[0].numpy() @ w)


def test_corrected_lincs_multi_target_genetic_rows_are_excluded_pre_support() -> None:
    lincs_meta = pd.DataFrame(
        [
            {
                "unique_id": "oe-bad",
                "cell_line": "A549",
                "target": "SMAD7|CPS1|SUZ12|FAM5C",
                "pert_type": "oe",
                "time_val": 96.0,
                "dose_val": np.nan,
                "paired_control_id": "ctrl-oe",
                "sig_id": "sig-oe",
            },
            {
                "unique_id": "crispr-good",
                "cell_line": "A549",
                "target": "SMAD7",
                "pert_type": "crispr",
                "time_val": 96.0,
                "dose_val": np.nan,
                "paired_control_id": "ctrl-cr",
                "sig_id": "sig-cr",
            },
        ]
    )

    supported, excluded = corrected_s3.normalize_lincs_metadata(lincs_meta)

    assert supported["treated_cell_id"].tolist() == ["crispr-good"]
    excluded_row = excluded.loc[excluded["treated_cell_id"] == "oe-bad"].iloc[0]
    assert excluded_row["exclusion_reason"] == "non_single_gene_genetic_target"


def test_membership_rebuild_uses_eligible_tokens_for_lincs_chemical_rows() -> None:
    delta_meta = pd.DataFrame(
        [
            {
                "row_id": 0,
                "source_row_index": 10,
                "treated_cell_id": "chem-1",
                "cell_line_raw": "A549",
                "dataset": "LINCS",
                "cell_line": "A549",
                "perturbation_class": "Chemical",
                "target_raw": "A|B",
                "target_tokens": "A;B",
                "parsed_target_tokens": ("A", "B"),
            },
            {
                "row_id": 1,
                "source_row_index": 11,
                "treated_cell_id": "gen-1",
                "cell_line_raw": "A549",
                "dataset": "LINCS",
                "cell_line": "A549",
                "perturbation_class": "Genetic",
                "target_raw": "A",
                "target_tokens": "A",
                "parsed_target_tokens": ("A",),
            },
        ]
    )
    eligible_subset = pd.DataFrame(
        [
            {"dataset": "LINCS", "cell_line": "A549", "target_token": "A"},
        ]
    )

    membership = task2_group.build_membership_from_delta_meta("LINCS", delta_meta, eligible_subset)

    assert membership["target_token"].tolist() == ["A", "A"]
    assert membership["treated_cell_id"].tolist() == ["chem-1", "gen-1"]
