from __future__ import annotations

import argparse
import glob
import hashlib
import json
import os
import time
from dataclasses import asdict, dataclass
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
import torch
from tqdm import tqdm


@dataclass
class Task0Config:
    processed_dir: str = "/mnt/NAS_21T/ProjectData/OSMOSIS/processed"
    output_dir: str = "./outputs/task0_curated"
    lincs_meta_name: str = "LINCS_Engine1_Meta.csv"
    lincs_tensor_name: str = "LINCS_Engine1_TrainData.pt"
    sc_meta_glob: str = "*Meta*.csv"
    save_bundle: bool = True
    save_parquet: bool = True
    save_unified_csv: bool = False
    save_tensors_separately: bool = True
    keep_tensors_in_bundle: bool = True
    write_tensor_manifest: bool = True
    ignore_genetic_dose_time_in_pair_meta: bool = True
    strict_load: bool = True
    max_pairs_per_group_level1: Optional[int] = None
    post_save_reload_check: bool = True
    reload_check_n: int = 2048
    enable_qc_gates: bool = True
    qc_seed: int = 42
    qc_zero_norm_thresh: float = 1e-12
    qc_max_zero_rate_lincs: float = 1e-3
    qc_spotcheck_n_lincs: int = 2000
    qc_spotcheck_max_mismatch_rate: float = 1e-4
    compute_sha256: bool = False

    @property
    def lincs_dir(self) -> str:
        return str(Path(self.processed_dir) / "LINCS_Processed")

    @property
    def sc_dir(self) -> str:
        return str(Path(self.processed_dir) / "scPerturb_Processed")

    @property
    def lincs_meta_path(self) -> str:
        return str(Path(self.lincs_dir) / self.lincs_meta_name)

    @property
    def lincs_data_path(self) -> str:
        return str(Path(self.lincs_dir) / self.lincs_tensor_name)

    @property
    def metadata_dir(self) -> str:
        return str(Path(self.output_dir) / "metadata")

    @property
    def bundle_dir(self) -> str:
        return str(Path(self.output_dir) / "bundle")

    @property
    def tensor_dir(self) -> str:
        return str(Path(self.output_dir) / "tensors" / "Level1_Tensors_Split")

    @property
    def bundle_path(self) -> str:
        return str(Path(self.bundle_dir) / "m2m_task0_bundle.pt")

    @property
    def tensor_manifest_path(self) -> str:
        return str(Path(self.tensor_dir) / "manifest_level1_tensors.json")

    @property
    def run_manifest_path(self) -> str:
        return str(Path(self.output_dir) / "run_manifest_task0.json")


def standardize_str(series: pd.Series) -> pd.Series:
    return series.astype(str).str.upper().str.strip()


def determine_lincs_modality(pert_type: str) -> str:
    pt = str(pert_type).lower()
    if pt in ["trt_cp", "drug"]:
        return "Chemical"
    if pt in ["trt_sh", "sh", "trt_oe", "oe", "trt_xpr", "crispr"]:
        return "Genetic"
    return "Other"


def determine_sc_modality(pert_type: str) -> str:
    pt = str(pert_type).lower()
    if any(x in pt for x in ["crispr", "cas9", "knockout", "ko", "kd", "activation", "crispri", "crispra"]):
        return "Genetic"
    if any(x in pt for x in ["drug", "compound", "small molecule", "treatment", "trt_cp"]):
        return "Chemical"
    return "Other"


def _is_multitarget(target_std: pd.Series) -> pd.Series:
    multi_sep = r"[\|,;/]|\bAND\b"
    return target_std.str.contains(multi_sep, regex=True, na=True)


def _bulk_flag(source_db: pd.Series) -> pd.Series:
    return source_db.eq("LINCS")


def make_condition_id(df: pd.DataFrame, dose_col: str = "dose_val", time_col: str = "time_val") -> pd.Series:
    dose = df[dose_col].astype(str).replace("nan", "NA")
    time_value = df[time_col].astype(str).replace("nan", "NA")
    return "DOSE=" + dose + "|TIME=" + time_value


def add_pair_id(df: pd.DataFrame) -> pd.DataFrame:
    out = df.copy()
    out["pair_id"] = out["cell_std"].astype(str) + ":" + out["target_std"].astype(str)
    return out


def atomic_torch_save(obj, path: str) -> None:
    tmp_path = path + ".tmp"
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(tmp_path, "wb") as handle:
        torch.save(obj, handle)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp_path, path)


def atomic_json_save(obj: dict, path: str) -> None:
    tmp_path = path + ".tmp"
    Path(path).parent.mkdir(parents=True, exist_ok=True)
    with open(tmp_path, "w", encoding="utf-8") as handle:
        json.dump(obj, handle, indent=2)
        handle.flush()
        os.fsync(handle.fileno())
    os.replace(tmp_path, path)


def file_sha256(path: str, chunk_size: int = 1024 * 1024) -> str:
    digest = hashlib.sha256()
    with open(path, "rb") as handle:
        while True:
            block = handle.read(chunk_size)
            if not block:
                break
            digest.update(block)
    return digest.hexdigest()


def tensor_basic_stats(tensor: torch.Tensor) -> dict:
    tensor_cpu = tensor.detach().cpu()
    return {
        "shape": list(tensor_cpu.shape),
        "dtype": str(tensor_cpu.dtype).replace("torch.", ""),
        "numel": int(tensor_cpu.numel()),
        "bytes": int(tensor_cpu.numel() * tensor_cpu.element_size()),
    }


def save_tensor_split(
    tensor: torch.Tensor,
    out_path: str,
    cfg: Task0Config,
    rng: np.random.Generator,
) -> dict:
    atomic_torch_save(tensor, out_path)
    file_size = os.path.getsize(out_path)

    sha256 = None
    if cfg.compute_sha256:
        sha256 = file_sha256(out_path)

    if cfg.post_save_reload_check and cfg.reload_check_n and tensor.ndim >= 1:
        tensor_reloaded = torch.load(out_path, map_location="cpu")
        n_rows = int(tensor.shape[0])
        n_check = int(min(cfg.reload_check_n, n_rows))
        if n_check > 0:
            idx = rng.choice(n_rows, size=n_check, replace=False)
            idx_t = torch.as_tensor(idx, dtype=torch.long)
            lhs = tensor.detach().cpu().index_select(0, idx_t)
            rhs = tensor_reloaded.detach().cpu().index_select(0, idx_t)
            if not torch.allclose(lhs, rhs):
                max_abs = float((lhs - rhs).abs().max().item())
                raise RuntimeError(f"[PostSaveCheck] Reload mismatch for {out_path}. max_abs_diff={max_abs}")

    return {"path": out_path, "file_size": int(file_size), "sha256": sha256}


def save_tensor_dict_split(
    tensors: dict,
    tag: str,
    out_dir: str,
    cfg: Task0Config,
    rng: np.random.Generator,
) -> dict:
    info = {}
    for key, value in tensors.items():
        if not torch.is_tensor(value):
            continue
        out_path = str(Path(out_dir) / f"{tag}_{key}.pt")
        meta = save_tensor_split(value, out_path, cfg=cfg, rng=rng)
        meta.update(tensor_basic_stats(value))
        info[f"{tag}:{key}"] = meta
        print(f"Saved split tensor: {out_path} | shape={meta['shape']} | dtype={meta['dtype']}")
    return info


def load_lincs_meta(cfg: Task0Config) -> pd.DataFrame:
    print("Loading LINCS metadata...")
    df = pd.read_csv(cfg.lincs_meta_path)
    df["global_idx"] = np.arange(len(df))
    df["cell_std"] = standardize_str(df["cell_line"])
    df["target_std"] = standardize_str(df["target"])
    df["modality"] = df["pert_type"].apply(determine_lincs_modality)

    mask = (
        df["modality"].isin(["Chemical", "Genetic"])
        & df["target"].notna()
        & (~_is_multitarget(df["target_std"]))
        & (df["cell_std"] != "UNKNOWN")
    )
    out = df.loc[mask].copy()
    out["source_db"] = "LINCS"
    out["chunk_file"] = cfg.lincs_tensor_name
    out["chunk_idx"] = -1
    return out


def infer_sc_data_filename(meta_filename: str) -> str:
    fname = os.path.basename(meta_filename)
    if "part" in fname:
        chunk_suffix = fname.split("Meta_")[1].replace(".csv", "")
        return f"scPerturb_Engine1_TrainData_{chunk_suffix}.pt"
    return fname.replace("_Meta.csv", ".pt")


def load_scperturb_meta(cfg: Task0Config) -> pd.DataFrame:
    print("Loading scPerturb metadata...")
    meta_files = sorted(glob.glob(str(Path(cfg.sc_dir) / cfg.sc_meta_glob)))
    if len(meta_files) == 0:
        raise FileNotFoundError(f"No scPerturb meta files found under: {cfg.sc_dir}")

    df_list: list[pd.DataFrame] = []
    for meta_path in tqdm(meta_files, desc="scPerturb meta files"):
        chunk_meta = pd.read_csv(meta_path)
        chunk_meta["chunk_file"] = infer_sc_data_filename(meta_path)
        chunk_meta["chunk_idx"] = np.arange(len(chunk_meta))
        df_list.append(chunk_meta)

    df = pd.concat(df_list, ignore_index=True)
    df["cell_std"] = standardize_str(df["cell_line"])
    df["target_std"] = standardize_str(df["target"])
    df["global_idx"] = np.arange(df.shape[0])
    df["modality"] = df["pert_type"].apply(determine_sc_modality)

    bad_targets = {"NON-TARGETING", "NONTARGETING", "NT", "CONTROL", "NEGATIVE"}
    mask = (
        df["modality"].isin(["Chemical", "Genetic"])
        & df["target"].notna()
        & (~df["target_std"].isin(bad_targets))
        & (df["cell_std"] != "UNKNOWN")
    )
    out = df.loc[mask].copy()
    out["source_db"] = "scPerturb"
    return out


def fetch_tensors_mixed(
    meta_df: pd.DataFrame,
    cfg: Task0Config,
    lincs_cache: Optional[dict] = None,
) -> tuple[dict, Optional[dict]]:
    required = {"source_db", "global_idx", "chunk_file", "chunk_idx"}
    missing_cols = required - set(meta_df.columns)
    if missing_cols:
        raise ValueError(f"meta_df missing required columns: {missing_cols}")

    retrieved = [None] * len(meta_df)
    cache = lincs_cache

    lincs_mask = meta_df["source_db"].eq("LINCS").values
    if lincs_mask.any():
        if cache is None:
            cache = torch.load(cfg.lincs_data_path, map_location="cpu")

        global_idxs = meta_df.loc[lincs_mask, "global_idx"].to_numpy()
        global_idxs_t = torch.as_tensor(global_idxs, dtype=torch.long)
        x_sub = cache["x_baseline"][global_idxs_t]
        yp_sub = cache["y_delta_pathway"][global_idxs_t]
        yg_sub = cache["y_delta_gene"][global_idxs_t]

        row_positions = np.flatnonzero(lincs_mask)
        for i, pos in enumerate(row_positions):
            retrieved[pos] = (x_sub[i], yp_sub[i], yg_sub[i])

    sc_mask = meta_df["source_db"].eq("scPerturb").values
    if sc_mask.any():
        sc_req = meta_df.loc[sc_mask].copy()
        for fname, group in tqdm(sc_req.groupby("chunk_file"), desc="scPerturb tensor chunks"):
            path = str(Path(cfg.sc_dir) / fname)
            if not os.path.exists(path):
                message = f"Missing scPerturb chunk: {path}"
                if cfg.strict_load:
                    raise FileNotFoundError(message)
                print("[WARN]", message)
                continue

            chunk = torch.load(path, map_location="cpu")
            local_idxs = group["chunk_idx"].to_numpy()
            local_idxs_t = torch.as_tensor(local_idxs, dtype=torch.long)
            x_sub = chunk["x_baseline"][local_idxs_t]
            yp_sub = chunk["y_delta_pathway"][local_idxs_t]
            yg_sub = chunk["y_delta_gene"][local_idxs_t]

            pos_in_meta = meta_df.index.get_indexer(group.index)
            for i, pos in enumerate(pos_in_meta):
                retrieved[pos] = (x_sub[i], yp_sub[i], yg_sub[i])

    missing = [idx for idx, value in enumerate(retrieved) if value is None]
    if missing:
        message = f"Tensor retrieval missing {len(missing)} rows (first 5 positions: {missing[:5]})"
        if cfg.strict_load:
            raise RuntimeError(message)
        print("[WARN]", message)

    x = torch.stack([value[0] for value in retrieved])
    y_path = torch.stack([value[1] for value in retrieved])
    y_gene = torch.stack([value[2] for value in retrieved])
    return {"x": x, "y_path": y_path, "y_gene": y_gene}, cache


def build_unified_pool(lincs_df: pd.DataFrame, sc_df: pd.DataFrame) -> pd.DataFrame:
    cols_to_keep = [
        "cell_std",
        "target_std",
        "modality",
        "source_db",
        "global_idx",
        "chunk_file",
        "chunk_idx",
        "pert_type",
        "dose_val",
        "time_val",
        "tissue",
        "disease",
    ]
    for column in cols_to_keep:
        if column not in lincs_df.columns:
            lincs_df[column] = np.nan
        if column not in sc_df.columns:
            sc_df[column] = np.nan

    unified = pd.concat([lincs_df[cols_to_keep], sc_df[cols_to_keep]], ignore_index=True)
    unified = unified.loc[unified["cell_std"] != "UNKNOWN"].copy()
    unified["uid"] = np.arange(len(unified), dtype=np.int64)
    unified["is_bulk"] = _bulk_flag(unified["source_db"])
    unified["cond_id"] = make_condition_id(unified)
    return add_pair_id(unified)


def meta_from_uids(unified_df: pd.DataFrame, uids: np.ndarray) -> pd.DataFrame:
    out = unified_df.set_index("uid").loc[uids, ["source_db", "global_idx", "chunk_file", "chunk_idx"]].copy()
    return out.reset_index(drop=True)


def build_level1_pairs_unique(chem_pool: pd.DataFrame, gene_pool: pd.DataFrame) -> pd.DataFrame:
    chem_pairs = chem_pool[["cell_std", "target_std"]].drop_duplicates()
    gene_pairs = gene_pool[["cell_std", "target_std"]].drop_duplicates()
    level1_pairs = chem_pairs.merge(gene_pairs, on=["cell_std", "target_std"], how="inner")
    return add_pair_id(level1_pairs)


def build_level1_instance_pairs(
    chem_pool: pd.DataFrame,
    gene_pool: pd.DataFrame,
    l1_pairs_unique: pd.DataFrame,
    cfg: Task0Config,
) -> pd.DataFrame:
    chem_l1 = chem_pool.merge(l1_pairs_unique[["cell_std", "target_std"]], on=["cell_std", "target_std"], how="inner")
    gene_l1 = gene_pool.merge(l1_pairs_unique[["cell_std", "target_std"]], on=["cell_std", "target_std"], how="inner")

    pairs = pd.merge(
        chem_l1,
        gene_l1,
        on=["cell_std", "target_std"],
        suffixes=("_chem", "_gene"),
        how="inner",
    )

    if cfg.max_pairs_per_group_level1 is not None:
        pairs["_grp"] = pairs["cell_std"].astype(str) + "||" + pairs["target_std"].astype(str)
        pairs = pairs.groupby("_grp", group_keys=False).head(cfg.max_pairs_per_group_level1).drop(columns=["_grp"])

    pairs["domain_combo"] = pairs["source_db_chem"].astype(str) + "->" + pairs["source_db_gene"].astype(str)
    pairs["bulk_combo"] = pairs["is_bulk_chem"].astype(str) + "->" + pairs["is_bulk_gene"].astype(str)
    pairs["chem_cond_id"] = pairs["cond_id_chem"]
    if cfg.ignore_genetic_dose_time_in_pair_meta:
        pairs["gene_cond_id"] = "IGNORED"
        pairs["dose_val_gene"] = np.nan
        pairs["time_val_gene"] = np.nan
    else:
        pairs["gene_cond_id"] = pairs["cond_id_gene"]

    keep = [
        "uid_chem",
        "uid_gene",
        "pair_id_chem",
        "pair_id_gene",
        "cell_std",
        "target_std",
        "modality_chem",
        "modality_gene",
        "source_db_chem",
        "source_db_gene",
        "is_bulk_chem",
        "is_bulk_gene",
        "domain_combo",
        "bulk_combo",
        "pert_type_chem",
        "pert_type_gene",
        "dose_val_chem",
        "time_val_chem",
        "dose_val_gene",
        "time_val_gene",
        "chem_cond_id",
        "gene_cond_id",
        "global_idx_chem",
        "chunk_file_chem",
        "chunk_idx_chem",
        "global_idx_gene",
        "chunk_file_gene",
        "chunk_idx_gene",
        "tissue_chem",
        "disease_chem",
        "tissue_gene",
        "disease_gene",
    ]
    keep = [column for column in keep if column in pairs.columns]
    out = pairs[keep].copy()
    return add_pair_id(out)


def pair_stats_from_pools(chem_pool: pd.DataFrame, gene_pool: pd.DataFrame) -> pd.DataFrame:
    chem = add_pair_id(chem_pool)
    gene = add_pair_id(gene_pool)

    c_stats = (
        chem.groupby(["pair_id", "cell_std", "target_std"])
        .agg(
            n_chem=("uid", "size"),
            chem_sources=("source_db", lambda s: "|".join(sorted(set(map(str, s))))),
            chem_bulk=("is_bulk", lambda s: int(np.any(s))),
        )
        .reset_index()
    )
    g_stats = (
        gene.groupby(["pair_id", "cell_std", "target_std"])
        .agg(
            n_gene=("uid", "size"),
            gene_sources=("source_db", lambda s: "|".join(sorted(set(map(str, s))))),
            gene_bulk=("is_bulk", lambda s: int(np.any(s))),
        )
        .reset_index()
    )

    out = pd.merge(c_stats, g_stats, on=["pair_id", "cell_std", "target_std"], how="outer")
    for column in ["n_chem", "n_gene"]:
        if column in out.columns:
            out[column] = out[column].fillna(0).astype(int)
    out["has_chem"] = out["n_chem"] > 0
    out["has_gene"] = out["n_gene"] > 0
    return out


def qc_zero_norm_rate_level1_side(
    level1_pairs: pd.DataFrame,
    tensors: dict,
    side: str,
    label: str,
    cfg: Task0Config,
) -> None:
    if not cfg.enable_qc_gates:
        return

    src_col = f"source_db_{side}"
    mask = level1_pairs[src_col].values == "LINCS"
    if mask.sum() == 0:
        print(f"[QC] {label}: no LINCS rows found; skip.")
        return

    y = tensors["y_gene"].detach().cpu()
    norms = torch.linalg.vector_norm(y[torch.as_tensor(np.where(mask)[0], dtype=torch.long)], dim=1)
    zero_rate = float((norms < cfg.qc_zero_norm_thresh).float().mean().item())
    print(f"[QC] {label}: LINCS zero-norm rate(y_gene)={zero_rate:.6f} (N={int(mask.sum())})")
    if zero_rate > cfg.qc_max_zero_rate_lincs:
        raise RuntimeError(
            f"[QC FAILED] {label}: zero-norm rate too high for LINCS y_gene: "
            f"{zero_rate:.6f} (threshold={cfg.qc_max_zero_rate_lincs})"
        )


def qc_spotcheck_lincs_vs_raw(
    level1_pairs: pd.DataFrame,
    tensors: dict,
    side: str,
    label: str,
    cfg: Task0Config,
    rng: np.random.Generator,
) -> None:
    if not cfg.enable_qc_gates:
        return

    src_col = f"source_db_{side}"
    gidx_col = f"global_idx_{side}"
    idx_all = np.where(level1_pairs[src_col].values == "LINCS")[0]
    if idx_all.size == 0:
        print(f"[QC] {label}: no LINCS rows found; skip spotcheck.")
        return

    n_check = int(min(cfg.qc_spotcheck_n_lincs, idx_all.size))
    pick = rng.choice(idx_all, size=n_check, replace=False)

    raw = torch.load(cfg.lincs_data_path, map_location="cpu")
    gidx = level1_pairs.loc[pick, gidx_col].astype(int).to_numpy()
    y_raw = raw["y_delta_gene"].index_select(0, torch.as_tensor(gidx, dtype=torch.long)).detach().cpu()
    y_new = tensors["y_gene"].index_select(0, torch.as_tensor(pick, dtype=torch.long)).detach().cpu()

    max_abs = (y_raw - y_new).abs().amax(dim=1).numpy()
    mismatch_rate = float((max_abs > 1e-6).mean())
    print(f"[QC] {label}: LINCS spotcheck mismatch_rate={mismatch_rate:.6f} (k={n_check})")

    if mismatch_rate > cfg.qc_spotcheck_max_mismatch_rate:
        worst = int(np.argmax(max_abs))
        raise RuntimeError(
            f"[QC FAILED] {label}: spotcheck mismatch too high: {mismatch_rate:.6f} "
            f"(threshold={cfg.qc_spotcheck_max_mismatch_rate}). "
            f"Worst max_abs_diff={float(max_abs[worst])}, row={int(pick[worst])}, global_idx={int(gidx[worst])}."
        )


def ensure_output_dirs(cfg: Task0Config) -> None:
    Path(cfg.metadata_dir).mkdir(parents=True, exist_ok=True)
    Path(cfg.bundle_dir).mkdir(parents=True, exist_ok=True)
    Path(cfg.tensor_dir).mkdir(parents=True, exist_ok=True)


def run_task0(cfg: Task0Config) -> None:
    t0 = time.time()
    rng = np.random.default_rng(cfg.qc_seed)
    ensure_output_dirs(cfg)

    lincs_df = load_lincs_meta(cfg)
    sc_df = load_scperturb_meta(cfg)

    unified_df = build_unified_pool(lincs_df, sc_df)
    print(f"Unified samples: {len(unified_df)}")
    print(unified_df.groupby(["source_db", "modality"]).size())

    chem_pool = unified_df[unified_df["modality"] == "Chemical"].copy()
    gene_pool = unified_df[unified_df["modality"] == "Genetic"].copy()

    level1_pairs_unique = build_level1_pairs_unique(chem_pool, gene_pool)
    level1_pairs = build_level1_instance_pairs(chem_pool, gene_pool, level1_pairs_unique, cfg)

    print(f"Level1 unique pairs: {len(level1_pairs_unique)}")
    print(f"Level1 instance pairs: {len(level1_pairs)}")

    lincs_cache = None
    chem_fetch = meta_from_uids(unified_df, level1_pairs["uid_chem"].to_numpy())
    chem_tensors, lincs_cache = fetch_tensors_mixed(chem_fetch, cfg, lincs_cache=lincs_cache)
    gene_fetch = meta_from_uids(unified_df, level1_pairs["uid_gene"].to_numpy())
    gene_tensors, lincs_cache = fetch_tensors_mixed(gene_fetch, cfg, lincs_cache=lincs_cache)

    qc_zero_norm_rate_level1_side(level1_pairs, chem_tensors, side="chem", label="Level1-Chem", cfg=cfg)
    qc_zero_norm_rate_level1_side(level1_pairs, gene_tensors, side="gene", label="Level1-Gene", cfg=cfg)
    qc_spotcheck_lincs_vs_raw(level1_pairs, chem_tensors, side="chem", label="Level1-Chem", cfg=cfg, rng=rng)
    qc_spotcheck_lincs_vs_raw(level1_pairs, gene_tensors, side="gene", label="Level1-Gene", cfg=cfg, rng=rng)

    split_manifest = {}
    if cfg.save_tensors_separately:
        split_manifest.update(save_tensor_dict_split(chem_tensors, "L1_chem", cfg.tensor_dir, cfg, rng))
        split_manifest.update(save_tensor_dict_split(gene_tensors, "L1_gene", cfg.tensor_dir, cfg, rng))

        if cfg.write_tensor_manifest:
            tensor_manifest = {
                "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
                "processed_dir": os.path.realpath(cfg.processed_dir),
                "lincs_data_path": os.path.realpath(cfg.lincs_data_path),
                "sc_dir": os.path.realpath(cfg.sc_dir),
                "tensor_out_dir": os.path.realpath(cfg.tensor_dir),
                "reload_check_n": int(cfg.reload_check_n) if cfg.post_save_reload_check else 0,
                "compute_sha256": bool(cfg.compute_sha256),
                "tensors": split_manifest,
            }
            atomic_json_save(tensor_manifest, cfg.tensor_manifest_path)
            print(f"Wrote tensor manifest: {cfg.tensor_manifest_path}")

    l1_targets = sorted(level1_pairs_unique["target_std"].astype(str).unique().tolist())
    l1_cells = sorted(level1_pairs_unique["cell_std"].astype(str).unique().tolist())

    level2_chem_pool = chem_pool[chem_pool["target_std"].isin(l1_targets)].copy()
    level2_gene_pool = gene_pool[gene_pool["target_std"].isin(l1_targets)].copy()
    level3_chem_pool = chem_pool[chem_pool["cell_std"].isin(l1_cells)].copy()
    level3_gene_pool = gene_pool[gene_pool["cell_std"].isin(l1_cells)].copy()

    level1_pair_stats = pair_stats_from_pools(
        chem_pool[chem_pool["pair_id"].isin(level1_pairs_unique["pair_id"])],
        gene_pool[gene_pool["pair_id"].isin(level1_pairs_unique["pair_id"])],
    )
    level2_pair_stats = pair_stats_from_pools(level2_chem_pool, level2_gene_pool)
    level3_pair_stats = pair_stats_from_pools(level3_chem_pool, level3_gene_pool)

    if cfg.save_bundle:
        output = {
            "unified_meta": unified_df,
            "Level1": {
                "pairs_meta": level1_pairs,
                "pairs_unique": level1_pairs_unique,
                "pairs_stats": level1_pair_stats,
                "chem_tensors": chem_tensors if cfg.keep_tensors_in_bundle else None,
                "gene_tensors": gene_tensors if cfg.keep_tensors_in_bundle else None,
            },
            "Level2": {
                "chem_pool": level2_chem_pool,
                "gene_pool": level2_gene_pool,
                "pair_stats": level2_pair_stats,
                "targets_from_level1": l1_targets,
            },
            "Level3": {
                "chem_pool": level3_chem_pool,
                "gene_pool": level3_gene_pool,
                "pair_stats": level3_pair_stats,
                "cells_from_level1": l1_cells,
            },
        }
        atomic_torch_save(output, cfg.bundle_path)
        print(f"Saved bundle: {cfg.bundle_path}")

    if cfg.save_parquet:
        metadata_dir = Path(cfg.metadata_dir)
        try:
            unified_df.to_parquet(metadata_dir / "unified_meta.parquet", index=False)
            level1_pairs.to_parquet(metadata_dir / "level1_pairs_meta.parquet", index=False)
            level1_pairs_unique.to_parquet(metadata_dir / "level1_pairs_unique.parquet", index=False)
            level1_pair_stats.to_parquet(metadata_dir / "level1_pair_stats.parquet", index=False)
            level2_chem_pool.to_parquet(metadata_dir / "level2_chem_pool.parquet", index=False)
            level2_gene_pool.to_parquet(metadata_dir / "level2_gene_pool.parquet", index=False)
            level2_pair_stats.to_parquet(metadata_dir / "level2_pair_stats.parquet", index=False)
            level3_chem_pool.to_parquet(metadata_dir / "level3_chem_pool.parquet", index=False)
            level3_gene_pool.to_parquet(metadata_dir / "level3_gene_pool.parquet", index=False)
            level3_pair_stats.to_parquet(metadata_dir / "level3_pair_stats.parquet", index=False)
            print(f"Saved parquet snapshots to: {cfg.metadata_dir}")
        except Exception as exc:
            print(f"[WARN] Failed to write parquet metadata: {exc}")
            print("[WARN] Proceeding without parquet outputs. Consider installing pyarrow/fastparquet.")

    if cfg.save_unified_csv:
        metadata_dir = Path(cfg.metadata_dir)
        metadata_dir.mkdir(parents=True, exist_ok=True)
        unified_df.to_csv(metadata_dir / "unified_meta.csv", index=False)
        print(f"Saved CSV snapshot: {metadata_dir / 'unified_meta.csv'}")

    run_manifest = {
        "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
        "elapsed_seconds": float(time.time() - t0),
        "config": asdict(cfg),
        "summary": {
            "n_unified": int(len(unified_df)),
            "n_level1_unique_pairs": int(len(level1_pairs_unique)),
            "n_level1_instance_pairs": int(len(level1_pairs)),
            "n_level2_chem": int(len(level2_chem_pool)),
            "n_level2_gene": int(len(level2_gene_pool)),
            "n_level3_chem": int(len(level3_chem_pool)),
            "n_level3_gene": int(len(level3_gene_pool)),
        },
        "outputs": {
            "bundle_path": cfg.bundle_path if cfg.save_bundle else None,
            "metadata_dir": cfg.metadata_dir if cfg.save_parquet else None,
            "tensor_dir": cfg.tensor_dir if cfg.save_tensors_separately else None,
            "tensor_manifest_path": cfg.tensor_manifest_path if cfg.write_tensor_manifest else None,
        },
    }
    atomic_json_save(run_manifest, cfg.run_manifest_path)
    print(f"Wrote run manifest: {cfg.run_manifest_path}")


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser("M2M-Bench Task0 data curation")
    parser.add_argument("--processed-dir", type=str, default=Task0Config.processed_dir)
    parser.add_argument("--output-dir", type=str, default=Task0Config.output_dir)
    parser.add_argument("--lincs-meta-name", type=str, default=Task0Config.lincs_meta_name)
    parser.add_argument("--lincs-tensor-name", type=str, default=Task0Config.lincs_tensor_name)
    parser.add_argument("--sc-meta-glob", type=str, default=Task0Config.sc_meta_glob)
    parser.add_argument("--max-pairs-per-group-level1", type=int, default=None)
    parser.add_argument("--reload-check-n", type=int, default=Task0Config.reload_check_n)
    parser.add_argument("--qc-seed", type=int, default=Task0Config.qc_seed)
    parser.add_argument("--qc-zero-norm-thresh", type=float, default=Task0Config.qc_zero_norm_thresh)
    parser.add_argument("--qc-max-zero-rate-lincs", type=float, default=Task0Config.qc_max_zero_rate_lincs)
    parser.add_argument("--qc-spotcheck-n-lincs", type=int, default=Task0Config.qc_spotcheck_n_lincs)
    parser.add_argument("--qc-spotcheck-max-mismatch-rate", type=float, default=Task0Config.qc_spotcheck_max_mismatch_rate)
    parser.add_argument("--save-bundle", action=argparse.BooleanOptionalAction, default=Task0Config.save_bundle)
    parser.add_argument("--save-parquet", action=argparse.BooleanOptionalAction, default=Task0Config.save_parquet)
    parser.add_argument("--save-unified-csv", action=argparse.BooleanOptionalAction, default=Task0Config.save_unified_csv)
    parser.add_argument(
        "--save-tensors-separately",
        action=argparse.BooleanOptionalAction,
        default=Task0Config.save_tensors_separately,
    )
    parser.add_argument(
        "--keep-tensors-in-bundle",
        action=argparse.BooleanOptionalAction,
        default=Task0Config.keep_tensors_in_bundle,
    )
    parser.add_argument(
        "--write-tensor-manifest",
        action=argparse.BooleanOptionalAction,
        default=Task0Config.write_tensor_manifest,
    )
    parser.add_argument(
        "--ignore-genetic-dose-time-in-pair-meta",
        action=argparse.BooleanOptionalAction,
        default=Task0Config.ignore_genetic_dose_time_in_pair_meta,
    )
    parser.add_argument("--strict-load", action=argparse.BooleanOptionalAction, default=Task0Config.strict_load)
    parser.add_argument(
        "--post-save-reload-check",
        action=argparse.BooleanOptionalAction,
        default=Task0Config.post_save_reload_check,
    )
    parser.add_argument("--enable-qc-gates", action=argparse.BooleanOptionalAction, default=Task0Config.enable_qc_gates)
    parser.add_argument("--compute-sha256", action=argparse.BooleanOptionalAction, default=Task0Config.compute_sha256)
    return parser


def main() -> None:
    args = _build_parser().parse_args()
    cfg = Task0Config(
        processed_dir=args.processed_dir,
        output_dir=args.output_dir,
        lincs_meta_name=args.lincs_meta_name,
        lincs_tensor_name=args.lincs_tensor_name,
        sc_meta_glob=args.sc_meta_glob,
        save_bundle=args.save_bundle,
        save_parquet=args.save_parquet,
        save_unified_csv=args.save_unified_csv,
        save_tensors_separately=args.save_tensors_separately,
        keep_tensors_in_bundle=args.keep_tensors_in_bundle,
        write_tensor_manifest=args.write_tensor_manifest,
        ignore_genetic_dose_time_in_pair_meta=args.ignore_genetic_dose_time_in_pair_meta,
        strict_load=args.strict_load,
        max_pairs_per_group_level1=args.max_pairs_per_group_level1,
        post_save_reload_check=args.post_save_reload_check,
        reload_check_n=args.reload_check_n,
        enable_qc_gates=args.enable_qc_gates,
        qc_seed=args.qc_seed,
        qc_zero_norm_thresh=args.qc_zero_norm_thresh,
        qc_max_zero_rate_lincs=args.qc_max_zero_rate_lincs,
        qc_spotcheck_n_lincs=args.qc_spotcheck_n_lincs,
        qc_spotcheck_max_mismatch_rate=args.qc_spotcheck_max_mismatch_rate,
        compute_sha256=args.compute_sha256,
    )
    run_task0(cfg)


if __name__ == "__main__":
    main()
