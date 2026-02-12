# -*- coding: utf-8 -*-
"""
Chem2Gen-Bench — Split Bundle Loader (functions.py)

This module provides a drop-in replacement loader for the original
Chem2Gen_Benchmark_MixedSources.pt bundle, using:
  - Parquet metadata snapshots in Benchmark_Datasets/
  - Split Level1 tensors in Benchmark_Datasets/Level1_Tensors_Split/
  - A manifest JSON describing tensor files

It returns a dict "data" compatible with prior scripts:
  data["unified_meta"] (pd.DataFrame or None)
  data["Level1"]["pairs_meta"] (pd.DataFrame)
  data["Level1"]["chem_tensors"]["y_gene"] (torch.Tensor, lazily loaded)
  data["Level2"]["chem_pool"] ...
  data["Level3"]["gene_pool"] ...
"""

from __future__ import annotations

import os
import json
from collections import OrderedDict
from typing import Dict, Optional, Any, Tuple

import torch
import pandas as pd


# -----------------------------
# Small I/O helpers
# -----------------------------
def _pjoin(*parts: str) -> str:
    return os.path.join(*parts)


def _exists(p: Optional[str]) -> bool:
    return p is not None and os.path.exists(p)


def _read_parquet(path: str) -> pd.DataFrame:
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    return pd.read_parquet(path)


# -----------------------------
# Lazy tensor dict with LRU cache
# -----------------------------
class LazyTensorDict(dict):
    """
    A dict-like object that loads torch tensors lazily from file paths.
    It caches loaded tensors with an LRU policy to control RAM usage.

    Access pattern is identical to a normal dict:
      t = data["Level1"]["chem_tensors"]["y_gene"]
    """

    def __init__(
        self,
        key_to_path: Dict[str, str],
        map_location: str = "cpu",
        cache_max_items: int = 2,
        name: str = "LazyTensorDict",
    ):
        super().__init__()
        self._paths = dict(key_to_path)
        self._map_location = map_location
        self._cache_max_items = int(cache_max_items)
        self._cache: "OrderedDict[str, torch.Tensor]" = OrderedDict()
        self._name = str(name)

    def __contains__(self, key: object) -> bool:
        return key in self._paths

    def keys(self):
        return self._paths.keys()

    def get(self, key: str, default=None):
        if key not in self._paths:
            return default
        return self.__getitem__(key)

    def _evict_if_needed(self):
        while len(self._cache) > self._cache_max_items:
            k, v = self._cache.popitem(last=False)
            # free tensor reference
            del v

    def __getitem__(self, key: str) -> torch.Tensor:
        if key in self._cache:
            self._cache.move_to_end(key, last=True)
            return self._cache[key]

        if key not in self._paths:
            raise KeyError(f"{self._name}: unknown tensor key '{key}'. Available: {sorted(self._paths.keys())}")

        path = self._paths[key]
        if not os.path.exists(path):
            raise FileNotFoundError(f"{self._name}: missing tensor file for key '{key}': {path}")

        t = torch.load(path, map_location=self._map_location)

        # cache + evict
        self._cache[key] = t
        self._cache.move_to_end(key, last=True)
        self._evict_if_needed()
        return t

    def materialize_all(self) -> Dict[str, torch.Tensor]:
        """
        Force-load all tensors and return a normal dict.
        Warning: can be very large (multiple 6.3G tensors).
        """
        out = {}
        for k in self._paths.keys():
            out[k] = self[k]
        return out


# -----------------------------
# Manifest parsing
# -----------------------------
def _resolve_path(p: str, tensor_dir: str) -> str:
    p = str(p)
    return p if os.path.isabs(p) else os.path.join(tensor_dir, p)


def _parse_manifest_level1_tensors(
    manifest_path: str,
    tensor_dir: str,
) -> Tuple[Dict[str, str], Dict[str, str]]:
    """
    Parse your manifest JSON and return:
      chem_paths = {"x": ..., "y_gene": ..., "y_path": ...}
      gene_paths = {"x": ..., "y_gene": ..., "y_path": ...}

    Supports manifest keys like:
      "L1_chem:x", "L1_chem:y_gene", ...
    Also supports fallback by basename if needed.
    """
    with open(manifest_path, "r") as f:
        man = json.load(f)

    tensors = man.get("tensors", {})
    if not isinstance(tensors, dict) or len(tensors) == 0:
        raise ValueError(f"Manifest has no 'tensors' dict: {manifest_path}")

    # 1) First try: parse keys like "L1_chem:x"
    chem_map: Dict[str, str] = {}
    gene_map: Dict[str, str] = {}

    for k, info in tensors.items():
        if not isinstance(info, dict):
            continue
        p = info.get("path", None)
        if not p:
            continue
        p = _resolve_path(p, tensor_dir)

        ks = str(k)
        # expected: "L1_chem:x" or "L1_gene:y_gene"
        if ks.startswith("L1_") and ":" in ks:
            left, right = ks.split(":", 1)
            # left: "L1_chem" or "L1_gene"
            side = left.replace("L1_", "").strip().lower()
            tkey = right.strip()
            if side == "chem":
                chem_map[tkey] = p
            elif side == "gene":
                gene_map[tkey] = p

    # 2) If still missing, build basename index as a robust fallback
    if len(chem_map) == 0 or len(gene_map) == 0:
        by_base: Dict[str, str] = {}
        for _, info in tensors.items():
            if not isinstance(info, dict):
                continue
            p = info.get("path", None)
            if not p:
                continue
            p = _resolve_path(p, tensor_dir)
            by_base[os.path.basename(p)] = p

        want = {
            "chem": {"x": "L1_chem_x.pt", "y_gene": "L1_chem_y_gene.pt", "y_path": "L1_chem_y_path.pt"},
            "gene": {"x": "L1_gene_x.pt", "y_gene": "L1_gene_y_gene.pt", "y_path": "L1_gene_y_path.pt"},
        }
        for tkey, bn in want["chem"].items():
            if tkey not in chem_map and bn in by_base:
                chem_map[tkey] = by_base[bn]
        for tkey, bn in want["gene"].items():
            if tkey not in gene_map and bn in by_base:
                gene_map[tkey] = by_base[bn]

    # Validate
    need = {"x", "y_gene", "y_path"}
    miss_c = sorted(list(need - set(chem_map.keys())))
    miss_g = sorted(list(need - set(gene_map.keys())))
    if miss_c or miss_g:
        raise KeyError(
            "Manifest missing required Level1 tensor entries.\n"
            f"Missing chem keys: {miss_c}\n"
            f"Missing gene keys: {miss_g}\n"
            f"Manifest: {manifest_path}\n"
            f"Found chem keys: {sorted(chem_map.keys())}\n"
            f"Found gene keys: {sorted(gene_map.keys())}"
        )

    return chem_map, gene_map


# -----------------------------
# Main loader
# -----------------------------
def load_chem2gen_bundle_from_split(
    base_dir: str,
    tensor_subdir: str = "Level1_Tensors_Split",
    manifest_name: str = "manifest_level1_tensors.json",
    load_unified_meta: bool = True,
    map_location: str = "cpu",
    tensor_cache_max_items: int = 2,
    strict: bool = True,
) -> Dict[str, Any]:
    """
    Load Chem2Gen bundle from split artifacts under base_dir.

    Expected files under base_dir:
      - unified_meta.parquet (optional)
      - level1_pairs_meta.parquet
      - level1_pairs_unique.parquet
      - level1_pair_stats.parquet
      - level2_chem_pool.parquet / level2_gene_pool.parquet / level2_pair_stats.parquet
      - level3_chem_pool.parquet / level3_gene_pool.parquet / level3_pair_stats.parquet
      - Level1_Tensors_Split/*.pt + manifest_level1_tensors.json

    Returns data dict compatible with old torch bundle.
    """
    base_dir = os.path.abspath(base_dir)
    if not os.path.isdir(base_dir):
        raise NotADirectoryError(base_dir)

    # ---- metadata ----
    paths = {
        "unified_meta": _pjoin(base_dir, "unified_meta.parquet"),
        "l1_meta": _pjoin(base_dir, "level1_pairs_meta.parquet"),
        "l1_unique": _pjoin(base_dir, "level1_pairs_unique.parquet"),
        "l1_stats": _pjoin(base_dir, "level1_pair_stats.parquet"),
        "l2c": _pjoin(base_dir, "level2_chem_pool.parquet"),
        "l2g": _pjoin(base_dir, "level2_gene_pool.parquet"),
        "l2stats": _pjoin(base_dir, "level2_pair_stats.parquet"),
        "l3c": _pjoin(base_dir, "level3_chem_pool.parquet"),
        "l3g": _pjoin(base_dir, "level3_gene_pool.parquet"),
        "l3stats": _pjoin(base_dir, "level3_pair_stats.parquet"),
    }

    if strict:
        for k in ["l1_meta", "l1_unique", "l1_stats", "l2c", "l2g", "l2stats", "l3c", "l3g", "l3stats"]:
            if not os.path.exists(paths[k]):
                raise FileNotFoundError(f"Missing required parquet: {paths[k]}")

    unified_meta = _read_parquet(paths["unified_meta"]) if (load_unified_meta and _exists(paths["unified_meta"])) else None
    l1_meta = _read_parquet(paths["l1_meta"]) if _exists(paths["l1_meta"]) else pd.DataFrame()
    l1_unique = _read_parquet(paths["l1_unique"]) if _exists(paths["l1_unique"]) else pd.DataFrame()
    l1_stats = _read_parquet(paths["l1_stats"]) if _exists(paths["l1_stats"]) else pd.DataFrame()

    l2c = _read_parquet(paths["l2c"]) if _exists(paths["l2c"]) else pd.DataFrame()
    l2g = _read_parquet(paths["l2g"]) if _exists(paths["l2g"]) else pd.DataFrame()
    l2stats = _read_parquet(paths["l2stats"]) if _exists(paths["l2stats"]) else pd.DataFrame()

    l3c = _read_parquet(paths["l3c"]) if _exists(paths["l3c"]) else pd.DataFrame()
    l3g = _read_parquet(paths["l3g"]) if _exists(paths["l3g"]) else pd.DataFrame()
    l3stats = _read_parquet(paths["l3stats"]) if _exists(paths["l3stats"]) else pd.DataFrame()

    # derive targets/cells from Level1 unique if available
    l1_targets = sorted(l1_unique["target_std"].astype(str).unique().tolist()) if ("target_std" in l1_unique.columns) else []
    l1_cells = sorted(l1_unique["cell_std"].astype(str).unique().tolist()) if ("cell_std" in l1_unique.columns) else []

    # ---- tensors ----
    tensor_dir = _pjoin(base_dir, tensor_subdir)
    manifest_path = _pjoin(tensor_dir, manifest_name)

    chem_paths: Dict[str, str]
    gene_paths: Dict[str, str]

    if os.path.exists(manifest_path):
        chem_paths, gene_paths = _parse_manifest_level1_tensors(manifest_path, tensor_dir=tensor_dir)
    else:
        # fallback to conventional filenames
        chem_paths = {
            "x": _pjoin(tensor_dir, "L1_chem_x.pt"),
            "y_gene": _pjoin(tensor_dir, "L1_chem_y_gene.pt"),
            "y_path": _pjoin(tensor_dir, "L1_chem_y_path.pt"),
        }
        gene_paths = {
            "x": _pjoin(tensor_dir, "L1_gene_x.pt"),
            "y_gene": _pjoin(tensor_dir, "L1_gene_y_gene.pt"),
            "y_path": _pjoin(tensor_dir, "L1_gene_y_path.pt"),
        }

    # strict existence checks for tensors
    if strict:
        for nm, p in {**{f"chem:{k}": v for k, v in chem_paths.items()},
                      **{f"gene:{k}": v for k, v in gene_paths.items()}}.items():
            if not os.path.exists(p):
                raise FileNotFoundError(f"Missing tensor file for {nm}: {p}")

    # Wrap as LazyTensorDict; keeps the same interface as normal dict
    chem_tensors = LazyTensorDict(
        key_to_path=chem_paths,
        map_location=map_location,
        cache_max_items=tensor_cache_max_items,
        name="Level1.chem_tensors",
    )
    gene_tensors = LazyTensorDict(
        key_to_path=gene_paths,
        map_location=map_location,
        cache_max_items=tensor_cache_max_items,
        name="Level1.gene_tensors",
    )

    data = {
        "unified_meta": unified_meta,
        "Level1": {
            "pairs_meta": l1_meta,
            "pairs_unique": l1_unique,
            "pairs_stats": l1_stats,
            "chem_tensors": chem_tensors,
            "gene_tensors": gene_tensors,
        },
        "Level2": {
            "chem_pool": l2c,
            "gene_pool": l2g,
            "pair_stats": l2stats,
            "targets_from_level1": l1_targets,
        },
        "Level3": {
            "chem_pool": l3c,
            "gene_pool": l3g,
            "pair_stats": l3stats,
            "cells_from_level1": l1_cells,
        },
        "_bundle_type": "split",
        "_base_dir": base_dir,
        "_tensor_dir": tensor_dir,
        "_manifest": manifest_path if os.path.exists(manifest_path) else None,
    }

    return data


# -----------------------------
# Optional: lightweight sanity checks (no raw LINCS/scPerturb needed)
# -----------------------------
def quick_check_split_bundle(data: Dict[str, Any], strict: bool = True) -> None:
    """
    Quick checks to ensure:
      - Level1 meta rows match tensor first dimension
      - required columns exist
    It only loads one tensor per side (y_gene) to check shape.
    """
    l1 = data["Level1"]["pairs_meta"]
    if strict and (l1 is None or len(l1) == 0):
        raise ValueError("Level1 pairs_meta is empty.")

    # load one tensor to validate shape
    chem_yg = data["Level1"]["chem_tensors"]["y_gene"]
    gene_yg = data["Level1"]["gene_tensors"]["y_gene"]

    n = len(l1)
    if strict:
        if chem_yg.shape[0] != n or gene_yg.shape[0] != n:
            raise ValueError(f"Row mismatch: pairs_meta={n}, chem_y_gene={chem_yg.shape[0]}, gene_y_gene={gene_yg.shape[0]}")

    # basic column presence
    must_cols = ["cell_std", "target_std", "source_db_chem", "source_db_gene"]
    miss = [c for c in must_cols if c not in l1.columns]
    if strict and miss:
        raise KeyError(f"pairs_meta missing required columns: {miss}")


if __name__ == "__main__":
    # Example usage
    BASE = "/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets"
    data = load_chem2gen_bundle_from_split(BASE, load_unified_meta=True, tensor_cache_max_items=2, strict=True)
    print("Loaded bundle keys:", list(data.keys()))
    print("Level1 rows:", len(data["Level1"]["pairs_meta"]))
    # Trigger a lazy tensor load:
    t = data["Level1"]["chem_tensors"]["y_gene"]
    print("chem y_gene:", tuple(t.shape), t.dtype, t.device)
