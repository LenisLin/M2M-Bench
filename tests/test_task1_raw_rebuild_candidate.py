from __future__ import annotations

import numpy as np
import pandas as pd
import pytest
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "src"))

from m2mbench.utils.raw_rebuild import (
    build_scperturb_target_fields,
    is_control_like,
    normalize_lincs_subtype,
    normalize_scperturb_subtype,
)
from scripts.build_task1_cross_contracts_v2_candidate import build_pair_contract, explode_membership
from scripts.build_task1_snapshot_v2_candidate import ScPerturbSource, build_scperturb_subtype_bundle


def test_subtype_normalization_and_control_detection() -> None:
    assert normalize_lincs_subtype("drug") == "Chemical"
    assert normalize_lincs_subtype("crispr") == "CRISPR"
    assert normalize_lincs_subtype("sh") == "sh"
    assert normalize_lincs_subtype("oe") == "oe"

    assert normalize_scperturb_subtype("drug") == "Chemical"
    assert normalize_scperturb_subtype("CRISPR-cas9") == "CRISPR"
    assert normalize_scperturb_subtype("CRISPRi") == "CRISPR"
    assert normalize_scperturb_subtype("cytokines") is None

    assert is_control_like("control")
    assert is_control_like("DMSO")
    assert not is_control_like("EGFR")


def test_scperturb_target_field_rules() -> None:
    target_raw, tokens, reason = build_scperturb_target_fields(
        subtype="CRISPR",
        target_value="",
        perturbation_value="STAT1",
    )
    assert reason is None
    assert target_raw == "STAT1"
    assert tokens == ("STAT1",)

    target_raw, tokens, reason = build_scperturb_target_fields(
        subtype="Chemical",
        target_value="HDAC1_HDAC2",
        perturbation_value="drug-x",
    )
    assert reason is None
    assert target_raw == "HDAC1_HDAC2"
    assert tokens == ("HDAC1", "HDAC2")


def test_cross_contract_builder_uses_token_membership() -> None:
    lincs_meta = pd.DataFrame(
        [
            {
                "snapshot_row_idx": 0,
                "perturbation_subtype": "Chemical",
                "cell_line": "A549",
                "cell_line_raw": "A549",
                "target_tokens": "A;B",
            },
            {
                "snapshot_row_idx": 1,
                "perturbation_subtype": "CRISPR",
                "cell_line": "A549",
                "cell_line_raw": "A549",
                "target_tokens": "A",
            },
            {
                "snapshot_row_idx": 2,
                "perturbation_subtype": "sh",
                "cell_line": "A549",
                "cell_line_raw": "A549",
                "target_tokens": "A",
            },
        ]
    )
    sc_chemical = pd.DataFrame(
        [
            {
                "delta_row_idx": 10,
                "cell_std": "A549",
                "cell_line_raw": "A549",
                "target_tokens": "B;C",
            }
        ]
    )
    sc_crispr = pd.DataFrame(
        [
            {
                "delta_row_idx": 20,
                "cell_std": "A549",
                "cell_line_raw": "A549",
                "target_tokens": "A",
            }
        ]
    )

    chemical_contract, chemical_proof = build_pair_contract(
        explode_membership(lincs_meta.loc[lincs_meta["perturbation_subtype"].eq("Chemical")], subtype="Chemical", row_id_column="snapshot_row_idx"),
        explode_membership(sc_chemical.assign(cell_line=sc_chemical["cell_std"]), subtype="Chemical", row_id_column="delta_row_idx"),
    )
    crispr_contract, crispr_proof = build_pair_contract(
        explode_membership(lincs_meta.loc[lincs_meta["perturbation_subtype"].eq("CRISPR")], subtype="CRISPR", row_id_column="snapshot_row_idx"),
        explode_membership(sc_crispr.assign(cell_line=sc_crispr["cell_std"]), subtype="CRISPR", row_id_column="delta_row_idx"),
    )

    assert chemical_contract[["target_token", "lincs_row_idx", "sc_delta_row_idx"]].to_dict("records") == [
        {"target_token": "B", "lincs_row_idx": 0, "sc_delta_row_idx": 10}
    ]
    assert int(chemical_proof["matched_pairs" if "matched_pairs" in chemical_proof.columns else "n_pairs"].iloc[0]) == 1
    assert crispr_contract[["target_token", "lincs_row_idx", "sc_delta_row_idx"]].to_dict("records") == [
        {"target_token": "A", "lincs_row_idx": 1, "sc_delta_row_idx": 20}
    ]
    assert int(crispr_proof["matched_pairs" if "matched_pairs" in crispr_proof.columns else "n_pairs"].iloc[0]) == 1


def test_scperturb_bundle_builds_crispr_fallback_and_chemical_tokens(tmp_path) -> None:
    ad = pytest.importorskip("anndata")

    obs = pd.DataFrame(
        [
            {
                "Unnamed: 0": "cell-control-cr",
                "cell_line": "K562",
                "perturbation_type": "CRISPR",
                "perturbation": "control",
                "target": "",
                "time": np.nan,
                "dose_value": np.nan,
            },
            {
                "Unnamed: 0": "cell-treated-cr",
                "cell_line": "K562",
                "perturbation_type": "CRISPR",
                "perturbation": "STAT1",
                "target": "",
                "time": np.nan,
                "dose_value": np.nan,
            },
            {
                "Unnamed: 0": "cell-control-drug",
                "cell_line": "K562",
                "perturbation_type": "drug",
                "perturbation": "control",
                "target": "control",
                "time": 24.0,
                "dose_value": 0.0,
            },
            {
                "Unnamed: 0": "cell-treated-drug",
                "cell_line": "K562",
                "perturbation_type": "drug",
                "perturbation": "Tubastatin",
                "target": "HDAC6_HDAC8",
                "time": 24.0,
                "dose_value": 1.0,
            },
        ]
    )
    x = np.array(
        [
            [0.0, 0.0, 0.0],
            [2.0, 0.0, 0.0],
            [0.0, 1.0, 1.0],
            [0.0, 3.0, 5.0],
        ],
        dtype=np.float32,
    )
    adata = ad.AnnData(
        X=x,
        obs=obs.copy(),
        var=pd.DataFrame(index=["STAT1", "HDAC6", "HDAC8"]),
    )
    obs_path = tmp_path / "Cleaned_demo_obs.csv"
    h5ad_path = tmp_path / "Cleaned_demo.h5ad"
    obs.to_csv(obs_path, index=False)
    adata.write_h5ad(h5ad_path)

    source = ScPerturbSource(dataset_id="demo", obs_path=obs_path, h5ad_path=h5ad_path)
    alignment = pd.DataFrame({"gene_symbol": ["STAT1", "HDAC6", "HDAC8"], "local_idx": [0, 1, 2]})
    pathway_w = np.array(
        [
            [1.0, 0.0],
            [0.0, 1.0],
            [1.0, 1.0],
        ],
        dtype=np.float32,
    )

    crispr = build_scperturb_subtype_bundle(
        subtype="CRISPR",
        sources=[source],
        alignment=alignment,
        pathway_w=pathway_w,
        seed=619,
    )
    chemical = build_scperturb_subtype_bundle(
        subtype="Chemical",
        sources=[source],
        alignment=alignment,
        pathway_w=pathway_w,
        seed=619,
    )

    assert crispr.meta["target_raw"].tolist() == ["STAT1"]
    assert crispr.meta["target_tokens"].tolist() == ["STAT1"]
    np.testing.assert_allclose(crispr.gene_delta[0], np.array([2.0, 0.0, 0.0], dtype=np.float32))

    assert chemical.meta["target_tokens"].tolist() == ["HDAC6;HDAC8"]
    np.testing.assert_allclose(chemical.gene_delta[0], np.array([0.0, 2.0, 4.0], dtype=np.float32))
    np.testing.assert_allclose(
        chemical.pathway_delta[0],
        np.array([4.0, 6.0], dtype=np.float32),
    )
