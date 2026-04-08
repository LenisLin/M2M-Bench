from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable, Iterator, Mapping, Sequence

import numpy as np
import pandas as pd

CONTROL_PATTERN = re.compile(r"(?:^|[\W_])(control|vehicle|dmso|non-target)(?:$|[\W_])", re.IGNORECASE)

TASK1_INTERNAL_LINCS_SUBTYPES = ("Chemical", "CRISPR", "sh")
TASK1_CROSS_SUBTYPES = ("Chemical", "CRISPR")
TASK2_ELIGIBLE_SUBTYPES = ("Chemical", "CRISPR")


@dataclass(frozen=True)
class ScPerturbSource:
    dataset_id: str
    obs_path: Path
    h5ad_path: Path


def normalize_required_text(value: object) -> str:
    if pd.isna(value):
        return ""
    return str(value).strip()


def normalize_optional_text(value: object, default: str = "NA") -> str:
    text = normalize_required_text(value)
    return text if text else default


def standardize_cell_line(value: object) -> str:
    text = normalize_required_text(value)
    if not text:
        return ""
    return re.sub(r"\s+", " ", text).upper()


def standardize_token(value: object) -> str:
    text = normalize_required_text(value)
    if not text:
        return ""
    return re.sub(r"\s+", " ", text).upper()


def split_target_tokens(raw_target: object, delimiters: Sequence[str]) -> tuple[str, ...]:
    chunks = [normalize_required_text(raw_target)]
    for delimiter in delimiters:
        next_chunks: list[str] = []
        for chunk in chunks:
            next_chunks.extend(chunk.split(delimiter))
        chunks = next_chunks
    tokens = sorted({standardize_token(chunk) for chunk in chunks if standardize_token(chunk)})
    return tuple(tokens)


def join_tokens(tokens: Sequence[str]) -> str:
    return ";".join(str(token) for token in tokens if str(token).strip())


def is_control_like(*values: object) -> bool:
    for value in values:
        text = normalize_required_text(value)
        if text and CONTROL_PATTERN.search(text):
            return True
    return False


def normalize_lincs_subtype(raw_value: object) -> str | None:
    raw = normalize_required_text(raw_value).lower()
    mapping = {
        "drug": "Chemical",
        "trt_cp": "Chemical",
        "crispr": "CRISPR",
        "trt_xpr": "CRISPR",
        "trt_sh": "sh",
        "sh": "sh",
        "oe": "oe",
        "trt_oe": "oe",
    }
    return mapping.get(raw)


def normalize_scperturb_subtype(raw_value: object) -> str | None:
    raw = normalize_required_text(raw_value).lower()
    if raw == "drug":
        return "Chemical"
    if raw.startswith("crispr"):
        return "CRISPR"
    return None


def perturbation_class_for_subtype(subtype: str) -> str:
    return "Chemical" if subtype == "Chemical" else "Genetic"


def choose_scperturb_cell_column(columns: Iterable[str]) -> str | None:
    seen = {str(column) for column in columns}
    if "cell_line" in seen:
        return "cell_line"
    if "celltype" in seen:
        return "celltype"
    return None


def discover_scperturb_sources(root: Path) -> list[ScPerturbSource]:
    sources: list[ScPerturbSource] = []
    for obs_path in sorted(root.glob("Cleaned_*_obs.csv")):
        h5ad_name = obs_path.name.replace("_obs.csv", ".h5ad")
        h5ad_path = obs_path.with_name(h5ad_name)
        if not h5ad_path.is_file():
            continue
        dataset_id = obs_path.name.removeprefix("Cleaned_").removesuffix("_obs.csv")
        sources.append(ScPerturbSource(dataset_id=dataset_id, obs_path=obs_path, h5ad_path=h5ad_path))
    return sources


def build_scperturb_target_fields(
    *,
    subtype: str,
    target_value: object,
    perturbation_value: object,
) -> tuple[str, tuple[str, ...], str | None]:
    if subtype == "CRISPR":
        target_raw = normalize_required_text(target_value) or normalize_required_text(perturbation_value)
        if not target_raw:
            return "", tuple(), "missing_target"
        if is_control_like(target_raw):
            return "", tuple(), "control_like_target"
        tokens = split_target_tokens(target_raw, delimiters=(";", "_", "|"))
        if len(tokens) != 1:
            return target_raw, tokens, "non_single_gene_target"
        return target_raw, tokens, None

    if subtype == "Chemical":
        target_raw = normalize_required_text(target_value)
        if not target_raw:
            return "", tuple(), "missing_target"
        if is_control_like(target_raw):
            return "", tuple(), "control_like_target"
        tokens = split_target_tokens(target_raw, delimiters=(";", "_", "|"))
        if not tokens:
            return target_raw, tuple(), "missing_target_tokens"
        return target_raw, tokens, None

    return "", tuple(), "unsupported_subtype"


def build_standardization_audit_rows(frame: pd.DataFrame) -> pd.DataFrame:
    if frame.empty:
        return pd.DataFrame(
            columns=[
                "perturbation_subtype",
                "cell_line_raw",
                "cell_line_std",
                "target_raw",
                "target_tokens",
            ]
        )
    out = frame.copy()
    out["target_tokens"] = out["target_tokens"].map(lambda tokens: join_tokens(tokens if isinstance(tokens, tuple) else ()))
    return out[
        [
            "perturbation_subtype",
            "cell_line_raw",
            "cell_line_std",
            "target_raw",
            "target_tokens",
        ]
    ].drop_duplicates().reset_index(drop=True)


def build_alignment_index(template_genes: Sequence[str], source_genes: Sequence[str]) -> np.ndarray:
    lookup = {standardize_token(gene): idx for idx, gene in enumerate(source_genes)}
    index = np.full((len(template_genes),), -1, dtype=np.int64)
    for local_idx, gene in enumerate(template_genes):
        mapped = lookup.get(standardize_token(gene))
        if mapped is not None:
            index[local_idx] = int(mapped)
    return index


def dataframe_from_records(records: Sequence[Mapping[str, object]]) -> pd.DataFrame:
    if not records:
        return pd.DataFrame()
    return pd.DataFrame.from_records(records)


def iter_supported_task1_subtypes(frame: pd.DataFrame, column: str = "perturbation_subtype") -> Iterator[str]:
    if column not in frame.columns:
        return iter(())
    values = frame[column].fillna("").astype(str)
    supported = [value for value in TASK1_INTERNAL_LINCS_SUBTYPES if values.eq(value).any()]
    return iter(supported)
