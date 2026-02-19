from __future__ import annotations

import pandas as pd


class DataContractError(ValueError):
    """Raised when a data contract is violated."""


def validate_for_r_export(df: pd.DataFrame, primary_key: str) -> None:
    if not isinstance(df, pd.DataFrame):
        raise DataContractError("Input must be a pandas DataFrame.")
    if df.shape[0] == 0:
        raise DataContractError("DataFrame is empty; refusing to export.")
    if primary_key not in df.columns:
        raise DataContractError(f"Primary key column '{primary_key}' missing.")
    if df[primary_key].isna().any():
        raise DataContractError(f"Primary key column '{primary_key}' contains NaN.")
    if not df[primary_key].is_unique:
        raise DataContractError(f"Primary key column '{primary_key}' is not unique.")
