import os
import glob
import json
import time
import hashlib
import numpy as np
import pandas as pd
import torch
from tqdm import tqdm

# ================= Configuration =================
BASE_DIR = "/mnt/NAS_21T/ProjectData/OSMOSIS/processed"
LINCS_DIR = os.path.join(BASE_DIR, "LINCS_Processed")
SC_DIR = os.path.join(BASE_DIR, "scPerturb_Processed")

LINCS_DATA_PATH = os.path.join(LINCS_DIR, "LINCS_Engine1_TrainData.pt")

OUTPUT_DIR = "/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets"
os.makedirs(OUTPUT_DIR, exist_ok=True)

SAVE_PATH = os.path.join(OUTPUT_DIR, "Chem2Gen_Benchmark_MixedSources.pt")
SAVE_META_PARQUET = True  # recommended

IGNORE_GENETIC_DOSE_TIME_IN_PAIR_META = True
STRICT_LOAD = True

# Optional caps to avoid accidental blow-ups (set None to disable)
MAX_PAIRS_PER_GROUP_LEVEL1 = None  # cap on instance-level merged pairs per pair_id

# -------- New: independent tensor saving --------
SAVE_TENSORS_SEPARATELY = True
TENSOR_OUT_DIR = os.path.join(OUTPUT_DIR, "Level1_Tensors_Split")  # directory to store split tensors
os.makedirs(TENSOR_OUT_DIR, exist_ok=True)

# Keep tensors in the big bundle as well (double storage, but safest)
KEEP_TENSORS_IN_BUNDLE = True

# Save manifest for easy auditing / rollback
WRITE_TENSOR_MANIFEST = True
MANIFEST_PATH = os.path.join(TENSOR_OUT_DIR, "manifest_level1_tensors.json")

# Post-save validation
POST_SAVE_RELOAD_CHECK = True
RELOAD_CHECK_N = 2048  # random rows to verify per tensor (set 0 to disable)

# Strong QC gates (recommended ON)
ENABLE_QC_GATES = True
QC_SEED = 42
QC_ZERO_NORM_THRESH = 1e-12
QC_MAX_ZERO_RATE_LINCS = 1e-3          # should be tiny; your failure case was ~0.39
QC_SPOTCHECK_N_LINCS = 2000            # rows to spot-check against raw for LINCS
QC_SPOTCHECK_MAX_MISMATCH_RATE = 1e-4  # tolerate extremely small numerical differences

# Optional: compute sha256 of each saved tensor file (expensive for huge files)
COMPUTE_SHA256 = False

# ================= Helper Functions =================

def standardize_str(series: pd.Series) -> pd.Series:
    return series.astype(str).str.upper().str.strip()

def determine_lincs_modality(pert_type) -> str:
    pt = str(pert_type).lower()
    if pt in ["trt_cp", "drug"]:
        return "Chemical"
    if pt in ["trt_sh", "sh", "trt_oe", "oe", "trt_xpr", "crispr"]:
        return "Genetic"
    return "Other"

def determine_sc_modality(pert_type) -> str:
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

def make_condition_id(df: pd.DataFrame, dose_col="dose_val", time_col="time_val") -> pd.Series:
    dose = df[dose_col].astype(str).replace("nan", "NA")
    time = df[time_col].astype(str).replace("nan", "NA")
    return "DOSE=" + dose + "|TIME=" + time

def add_pair_id(df: pd.DataFrame) -> pd.DataFrame:
    df = df.copy()
    df["pair_id"] = df["cell_std"].astype(str) + ":" + df["target_std"].astype(str)
    return df

def atomic_torch_save(obj, path: str):
    """
    Atomic save: write to temp file -> flush+fsync -> os.replace
    Prevents partially-written artifacts when IO is unstable.
    """
    tmp_path = path + ".tmp"
    os.makedirs(os.path.dirname(path), exist_ok=True)
    with open(tmp_path, "wb") as f:
        torch.save(obj, f)
        f.flush()
        os.fsync(f.fileno())
    os.replace(tmp_path, path)

def file_sha256(path: str, chunk_size: int = 1024 * 1024) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as f:
        while True:
            b = f.read(chunk_size)
            if not b:
                break
            h.update(b)
    return h.hexdigest()

def tensor_basic_stats(t: torch.Tensor) -> dict:
    t_cpu = t.detach().cpu()
    return {
        "shape": list(t_cpu.shape),
        "dtype": str(t_cpu.dtype).replace("torch.", ""),
        "numel": int(t_cpu.numel()),
        "bytes": int(t_cpu.numel() * t_cpu.element_size()),
    }

def save_tensor_split(tensor: torch.Tensor, out_path: str, reload_check_n: int = RELOAD_CHECK_N, rng: np.random.Generator | None = None):
    atomic_torch_save(tensor, out_path)
    size = os.path.getsize(out_path)

    sha = None
    if COMPUTE_SHA256:
        sha = file_sha256(out_path)

    # reload check (random rows)
    if POST_SAVE_RELOAD_CHECK and reload_check_n and tensor.ndim >= 1:
        tt = torch.load(out_path, map_location="cpu")
        if rng is None:
            rng = np.random.default_rng(QC_SEED)

        n = int(tensor.shape[0])
        k = int(min(reload_check_n, n))
        if k > 0:
            idx = rng.choice(n, size=k, replace=False)
            idx_t = torch.as_tensor(idx, dtype=torch.long)
            a = tensor.detach().cpu().index_select(0, idx_t)
            b = tt.detach().cpu().index_select(0, idx_t)
            if not torch.allclose(a, b):
                # tighten diagnosis
                max_abs = float((a - b).abs().max().item())
                raise RuntimeError(f"[PostSaveCheck] Reload mismatch for {out_path}. max_abs_diff={max_abs}")

    return {"path": out_path, "file_size": int(size), "sha256": sha}

def save_tensor_dict_split(tensors: dict, tag: str, out_dir: str, rng: np.random.Generator) -> dict:
    """
    tensors: dict like {"x": Tensor, "y_path": Tensor, "y_gene": Tensor}
    tag: "L1_chem" or "L1_gene"
    """
    info = {}
    for k, v in tensors.items():
        if not torch.is_tensor(v):
            continue
        fname = f"{tag}_{k}.pt"
        out_path = os.path.join(out_dir, fname)
        meta = save_tensor_split(v, out_path, rng=rng)
        meta.update(tensor_basic_stats(v))
        info[f"{tag}:{k}"] = meta
        print(f"    ✅ Saved split tensor: {out_path} | {meta['shape']} | {meta['dtype']} | {meta['file_size']/1e9:.2f} GB")
    return info

def load_lincs_meta() -> pd.DataFrame:
    print(">>> [1/Loading] LINCS Meta...")
    meta_path = os.path.join(LINCS_DIR, "LINCS_Engine1_Meta.csv")
    df = pd.read_csv(meta_path)

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
    df = df.loc[mask].copy()
    df["source_db"] = "LINCS"
    df["chunk_file"] = "LINCS_Engine1_TrainData.pt"
    df["chunk_idx"] = -1
    return df

def load_scperturb_meta() -> pd.DataFrame:
    print(">>> [2/Loading] scPerturb Meta...")
    meta_files = sorted(glob.glob(os.path.join(SC_DIR, "*Meta*.csv")))
    if len(meta_files) == 0:
        raise FileNotFoundError(f"No scPerturb meta files found under: {SC_DIR}")

    df_list = []
    for f in tqdm(meta_files):
        temp = pd.read_csv(f)
        fname = os.path.basename(f)

        if "part" in fname:
            chunk_suffix = fname.split("Meta_")[1].replace(".csv", "")
            data_filename = f"scPerturb_Engine1_TrainData_{chunk_suffix}.pt"
        else:
            data_filename = fname.replace("_Meta.csv", ".pt")

        temp["chunk_file"] = data_filename
        temp["chunk_idx"] = np.arange(len(temp))
        df_list.append(temp)

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
    df = df.loc[mask].copy()
    df["source_db"] = "scPerturb"
    return df

_LINCS_CACHE = None

def fetch_tensors_mixed(meta_df: pd.DataFrame):
    required = {"source_db", "global_idx", "chunk_file", "chunk_idx"}
    missing_cols = required - set(meta_df.columns)
    if missing_cols:
        raise ValueError(f"meta_df missing required columns: {missing_cols}")

    retrieved = [None] * len(meta_df)

    # LINCS
    lincs_mask = meta_df["source_db"].eq("LINCS").values
    if lincs_mask.any():
        global _LINCS_CACHE
        if _LINCS_CACHE is None:
            # map_location=cpu to avoid any device-related surprises
            _LINCS_CACHE = torch.load(LINCS_DATA_PATH, map_location="cpu")

        lincs_data = _LINCS_CACHE
        global_idxs = meta_df.loc[lincs_mask, "global_idx"].to_numpy()
        global_idxs_t = torch.as_tensor(global_idxs, dtype=torch.long)

        x_sub = lincs_data["x_baseline"][global_idxs_t]
        yp_sub = lincs_data["y_delta_pathway"][global_idxs_t]
        yg_sub = lincs_data["y_delta_gene"][global_idxs_t]

        row_positions = np.flatnonzero(lincs_mask)
        for i, pos in enumerate(row_positions):
            retrieved[pos] = (x_sub[i], yp_sub[i], yg_sub[i])

    # scPerturb
    sc_mask = meta_df["source_db"].eq("scPerturb").values
    if sc_mask.any():
        sc_req = meta_df.loc[sc_mask].copy()
        for fname, group in tqdm(sc_req.groupby("chunk_file"), desc="       Loading scPerturb chunks"):
            path = os.path.join(SC_DIR, fname)
            if not os.path.exists(path):
                msg = f"Missing scPerturb chunk: {path}"
                if STRICT_LOAD:
                    raise FileNotFoundError(msg)
                print("[WARN]", msg)
                continue

            chunk = torch.load(path, map_location="cpu")
            local_idxs = group["chunk_idx"].to_numpy()
            local_idxs_t = torch.as_tensor(local_idxs, dtype=torch.long)

            x_sub = chunk["x_baseline"][local_idxs_t]
            yp_sub = chunk["y_delta_pathway"][local_idxs_t]
            yg_sub = chunk["y_delta_gene"][local_idxs_t]

            # group.index are original indices in meta_df; map to positional indices
            pos_in_meta = meta_df.index.get_indexer(group.index)
            for i, pos in enumerate(pos_in_meta):
                retrieved[pos] = (x_sub[i], yp_sub[i], yg_sub[i])

            del chunk

    missing = [i for i, v in enumerate(retrieved) if v is None]
    if missing:
        msg = f"Tensor retrieval missing {len(missing)} rows (first 5 positions: {missing[:5]})"
        if STRICT_LOAD:
            raise RuntimeError(msg)
        print("[WARN]", msg)

    x = torch.stack([v[0] for v in retrieved])
    y_path = torch.stack([v[1] for v in retrieved])
    y_gene = torch.stack([v[2] for v in retrieved])

    return {"x": x, "y_path": y_path, "y_gene": y_gene}

def build_unified_pool(lincs_df: pd.DataFrame, sc_df: pd.DataFrame) -> pd.DataFrame:
    cols_to_keep = [
        "cell_std", "target_std", "modality", "source_db",
        "global_idx", "chunk_file", "chunk_idx", "pert_type",
        "dose_val", "time_val", "tissue", "disease"
    ]
    for col in cols_to_keep:
        if col not in lincs_df.columns:
            lincs_df[col] = np.nan
        if col not in sc_df.columns:
            sc_df[col] = np.nan

    unified = pd.concat([lincs_df[cols_to_keep], sc_df[cols_to_keep]], ignore_index=True)
    unified = unified.loc[unified["cell_std"] != "UNKNOWN"].copy()

    unified["uid"] = np.arange(len(unified), dtype=np.int64)
    unified["is_bulk"] = _bulk_flag(unified["source_db"])
    unified["cond_id"] = make_condition_id(unified)

    unified = add_pair_id(unified)
    return unified

def meta_from_uids(unified_df: pd.DataFrame, uids: np.ndarray) -> pd.DataFrame:
    sub = unified_df.set_index("uid").loc[uids, ["source_db", "global_idx", "chunk_file", "chunk_idx"]].copy()
    sub = sub.reset_index(drop=True)
    return sub

# ---------- Level1 coverage first, then build instance pairs ----------
def build_level1_pairs_unique(chem_pool: pd.DataFrame, gene_pool: pd.DataFrame) -> pd.DataFrame:
    """
    Pair_ID-level L1 coverage: unique (cell_std, target_std) that exist in BOTH modalities.
    """
    c_pairs = chem_pool[["cell_std", "target_std"]].drop_duplicates()
    g_pairs = gene_pool[["cell_std", "target_std"]].drop_duplicates()
    l1_pairs = c_pairs.merge(g_pairs, on=["cell_std", "target_std"], how="inner")
    l1_pairs = add_pair_id(l1_pairs)
    return l1_pairs

def build_level1_instance_pairs(chem_pool: pd.DataFrame, gene_pool: pd.DataFrame, l1_pairs_unique: pd.DataFrame) -> pd.DataFrame:
    """
    Instance-level many-to-many pairs within L1 coverage.
    Keeps your previous pointer-rich schema.
    """
    chem_l1 = chem_pool.merge(l1_pairs_unique[["cell_std", "target_std"]], on=["cell_std", "target_std"], how="inner")
    gene_l1 = gene_pool.merge(l1_pairs_unique[["cell_std", "target_std"]], on=["cell_std", "target_std"], how="inner")

    pairs = pd.merge(
        chem_l1,
        gene_l1,
        on=["cell_std", "target_std"],
        suffixes=("_chem", "_gene"),
        how="inner",
    )

    if MAX_PAIRS_PER_GROUP_LEVEL1 is not None:
        pairs["_grp"] = pairs["cell_std"].astype(str) + "||" + pairs["target_std"].astype(str)
        pairs = (
            pairs.groupby("_grp", group_keys=False)
                 .head(MAX_PAIRS_PER_GROUP_LEVEL1)
                 .drop(columns=["_grp"])
        )

    pairs["domain_combo"] = pairs["source_db_chem"].astype(str) + "->" + pairs["source_db_gene"].astype(str)
    pairs["bulk_combo"] = pairs["is_bulk_chem"].astype(str) + "->" + pairs["is_bulk_gene"].astype(str)

    pairs["chem_cond_id"] = pairs["cond_id_chem"]
    if IGNORE_GENETIC_DOSE_TIME_IN_PAIR_META:
        pairs["gene_cond_id"] = "IGNORED"
        pairs["dose_val_gene"] = np.nan
        pairs["time_val_gene"] = np.nan
    else:
        pairs["gene_cond_id"] = pairs["cond_id_gene"]

    keep = [
        "uid_chem", "uid_gene",
        "pair_id_chem", "pair_id_gene",
        "cell_std", "target_std",
        "modality_chem", "modality_gene",
        "source_db_chem", "source_db_gene",
        "is_bulk_chem", "is_bulk_gene",
        "domain_combo", "bulk_combo",
        "pert_type_chem", "pert_type_gene",
        "dose_val_chem", "time_val_chem", "dose_val_gene", "time_val_gene",
        "chem_cond_id", "gene_cond_id",
        "global_idx_chem", "chunk_file_chem", "chunk_idx_chem",
        "global_idx_gene", "chunk_file_gene", "chunk_idx_gene",
        "tissue_chem", "disease_chem", "tissue_gene", "disease_gene",
    ]
    keep = [c for c in keep if c in pairs.columns]
    out = pairs[keep].copy()
    out = add_pair_id(out)
    return out

def pair_stats_from_pools(chem_pool: pd.DataFrame, gene_pool: pd.DataFrame) -> pd.DataFrame:
    """
    Pair_ID-level summary for a level's pools.
    Counts instances per modality and records source composition (lightweight).
    """
    c = add_pair_id(chem_pool)
    g = add_pair_id(gene_pool)

    c_stats = (c.groupby(["pair_id", "cell_std", "target_std"])
                 .agg(n_chem=("uid", "size"),
                      chem_sources=("source_db", lambda s: "|".join(sorted(set(map(str, s))))),
                      chem_bulk=("is_bulk", lambda s: int(np.any(s))))
                 .reset_index())
    g_stats = (g.groupby(["pair_id", "cell_std", "target_std"])
                 .agg(n_gene=("uid", "size"),
                      gene_sources=("source_db", lambda s: "|".join(sorted(set(map(str, s))))),
                      gene_bulk=("is_bulk", lambda s: int(np.any(s))))
                 .reset_index())

    out = pd.merge(c_stats, g_stats, on=["pair_id", "cell_std", "target_std"], how="outer")
    for col in ["n_chem", "n_gene"]:
        if col in out.columns:
            out[col] = out[col].fillna(0).astype(int)
    out["has_chem"] = out["n_chem"] > 0
    out["has_gene"] = out["n_gene"] > 0
    return out

# ================= QC Gates =================

def qc_zero_norm_rate_level1_side(level1_pairs: pd.DataFrame, tensors: dict, side: str, label: str):
    """
    side: "chem" or "gene"
    checks LINCS-only subset zero-norm rate for y_gene
    """
    if not ENABLE_QC_GATES:
        return

    assert side in ("chem", "gene")
    src_col = f"source_db_{side}"
    mask = (level1_pairs[src_col].values == "LINCS")
    if mask.sum() == 0:
        print(f"[QC] {label}: no LINCS rows found; skip.")
        return

    y = tensors["y_gene"].detach().cpu()
    n = torch.linalg.vector_norm(y[torch.as_tensor(np.where(mask)[0], dtype=torch.long)], dim=1)
    zero_rate = float((n < QC_ZERO_NORM_THRESH).float().mean().item())
    print(f"[QC] {label}: LINCS zero-norm rate (y_gene) = {zero_rate:.6f} (N={int(mask.sum())})")
    if zero_rate > QC_MAX_ZERO_RATE_LINCS:
        raise RuntimeError(
            f"[QC FAILED] {label}: zero-norm rate too high for LINCS y_gene: {zero_rate:.6f} "
            f"(threshold={QC_MAX_ZERO_RATE_LINCS})."
        )

def qc_spotcheck_lincs_vs_raw(level1_pairs: pd.DataFrame, tensors: dict, side: str, label: str, rng: np.random.Generator):
    """
    spot-check some LINCS rows: fetched y_gene should match raw y_delta_gene at global_idx.
    This catches the class of failures you observed (many zeros in saved tensor but raw is nonzero).
    """
    if not ENABLE_QC_GATES:
        return

    assert side in ("chem", "gene")
    src_col = f"source_db_{side}"
    gidx_col = f"global_idx_{side}"

    mask = (level1_pairs[src_col].values == "LINCS")
    idx_all = np.where(mask)[0]
    if idx_all.size == 0:
        print(f"[QC] {label}: no LINCS rows found; skip spotcheck.")
        return

    k = int(min(QC_SPOTCHECK_N_LINCS, idx_all.size))
    pick = rng.choice(idx_all, size=k, replace=False)

    # IMPORTANT: load raw freshly (not via _LINCS_CACHE) for QC robustness
    lincs_raw = torch.load(LINCS_DATA_PATH, map_location="cpu")
    gidx = level1_pairs.loc[pick, gidx_col].astype(int).to_numpy()
    Yraw = lincs_raw["y_delta_gene"].index_select(0, torch.as_tensor(gidx, dtype=torch.long)).detach().cpu()
    Ynew = tensors["y_gene"].index_select(0, torch.as_tensor(pick, dtype=torch.long)).detach().cpu()

    mx = (Yraw - Ynew).abs().amax(dim=1).numpy()
    mismatch_rate = float((mx > 1e-6).mean())
    print(f"[QC] {label}: LINCS spotcheck mismatch_rate(mx>1e-6) = {mismatch_rate:.6f} (k={k})")

    if mismatch_rate > QC_SPOTCHECK_MAX_MISMATCH_RATE:
        worst = int(np.argmax(mx))
        raise RuntimeError(
            f"[QC FAILED] {label}: spotcheck mismatch too high: {mismatch_rate:.6f} "
            f"(threshold={QC_SPOTCHECK_MAX_MISMATCH_RATE}). "
            f"Worst max_abs_diff={float(mx[worst])}, row={int(pick[worst])}, global_idx={int(gidx[worst])}."
        )

# ================= Main =================

if __name__ == "__main__":
    t0 = time.time()
    rng = np.random.default_rng(QC_SEED)

    # 1) Load metadata
    lincs_df = load_lincs_meta()
    sc_df = load_scperturb_meta()

    # 2) Unified pool
    unified_df = build_unified_pool(lincs_df, sc_df)
    print(f"Total Unified Samples: {len(unified_df)}")
    print(unified_df.groupby(["source_db", "modality"]).size())

    # 3) Split pools by modality (instance-level pools)
    chem_pool = unified_df[unified_df["modality"] == "Chemical"].copy()
    gene_pool = unified_df[unified_df["modality"] == "Genetic"].copy()

    # 4) Level1: coverage first (Pair_ID), then instance pairs
    print(">>> Building Level1 coverage (unique cell:target pairs)...")
    level1_pairs_unique = build_level1_pairs_unique(chem_pool, gene_pool)
    print(f"✅ Level1 unique pairs: {len(level1_pairs_unique)}")

    print(">>> Building Level1 instance-level pairs (many-to-many within each pair)...")
    level1_pairs = build_level1_instance_pairs(chem_pool, gene_pool, level1_pairs_unique)
    print(f"✅ Level1 instance pairs (rows): {len(level1_pairs)}")
    if "domain_combo" in level1_pairs.columns:
        print(level1_pairs.groupby(["domain_combo"]).size().sort_values(ascending=False).head(20))

    # 5) Extract tensors for Level1 (pointer-based, instance-level order preserved)
    print(">>> Extracting Chemical tensors for Level1...")
    chem_fetch = meta_from_uids(unified_df, level1_pairs["uid_chem"].to_numpy())
    chem_tensors = fetch_tensors_mixed(chem_fetch)

    print(">>> Extracting Genetic tensors for Level1...")
    gene_fetch = meta_from_uids(unified_df, level1_pairs["uid_gene"].to_numpy())
    gene_tensors = fetch_tensors_mixed(gene_fetch)

    # --- QC gates (fail fast) ---
    qc_zero_norm_rate_level1_side(level1_pairs, chem_tensors, side="chem", label="Level1-Chem")
    qc_zero_norm_rate_level1_side(level1_pairs, gene_tensors, side="gene", label="Level1-Gene")
    qc_spotcheck_lincs_vs_raw(level1_pairs, chem_tensors, side="chem", label="Level1-Chem", rng=rng)
    qc_spotcheck_lincs_vs_raw(level1_pairs, gene_tensors, side="gene", label="Level1-Gene", rng=rng)

    # 5.5) Save split tensors (independent files) + post-save reload checks
    split_manifest = {}
    if SAVE_TENSORS_SEPARATELY:
        print(f">>> Saving split tensors to: {TENSOR_OUT_DIR}")
        split_manifest.update(save_tensor_dict_split(chem_tensors, tag="L1_chem", out_dir=TENSOR_OUT_DIR, rng=rng))
        split_manifest.update(save_tensor_dict_split(gene_tensors, tag="L1_gene", out_dir=TENSOR_OUT_DIR, rng=rng))

        if WRITE_TENSOR_MANIFEST:
            manifest_obj = {
                "created_at": time.strftime("%Y-%m-%d %H:%M:%S", time.localtime()),
                "lincs_data_path": os.path.realpath(LINCS_DATA_PATH),
                "sc_dir": os.path.realpath(SC_DIR),
                "tensor_out_dir": os.path.realpath(TENSOR_OUT_DIR),
                "reload_check_n": int(RELOAD_CHECK_N) if POST_SAVE_RELOAD_CHECK else 0,
                "compute_sha256": bool(COMPUTE_SHA256),
                "tensors": split_manifest,
            }
            tmp = MANIFEST_PATH + ".tmp"
            with open(tmp, "w", encoding="utf-8") as f:
                json.dump(manifest_obj, f, indent=2)
                f.flush()
                os.fsync(f.fileno())
            os.replace(tmp, MANIFEST_PATH)
            print(f"    ✅ Wrote manifest: {MANIFEST_PATH}")

    # 6) Level2/Level3 derived STRICTLY from Level1 expansion
    l1_targets = sorted(level1_pairs_unique["target_std"].astype(str).unique().tolist())
    l1_cells = sorted(level1_pairs_unique["cell_std"].astype(str).unique().tolist())

    print(f">>> Level1-derived sets: targets={len(l1_targets)} | cells={len(l1_cells)}")

    level2_chem_pool = chem_pool[chem_pool["target_std"].isin(l1_targets)].copy()
    level2_gene_pool = gene_pool[gene_pool["target_std"].isin(l1_targets)].copy()

    level3_chem_pool = chem_pool[chem_pool["cell_std"].isin(l1_cells)].copy()
    level3_gene_pool = gene_pool[gene_pool["cell_std"].isin(l1_cells)].copy()

    # 7) Pair-level summaries for visualization / diagnostics
    level1_pair_stats = pair_stats_from_pools(
        chem_pool[chem_pool["pair_id"].isin(level1_pairs_unique["pair_id"])],
        gene_pool[gene_pool["pair_id"].isin(level1_pairs_unique["pair_id"])]
    )
    level2_pair_stats = pair_stats_from_pools(level2_chem_pool, level2_gene_pool)
    level3_pair_stats = pair_stats_from_pools(level3_chem_pool, level3_gene_pool)

    # 8) Save bundle (atomic)
    if not KEEP_TENSORS_IN_BUNDLE:
        chem_tensors_bundle = None
        gene_tensors_bundle = None
    else:
        chem_tensors_bundle = chem_tensors
        gene_tensors_bundle = gene_tensors

    output = {
        "unified_meta": unified_df,  # large; parquet recommended
        "Level1": {
            "pairs_meta": level1_pairs,            # instance-level many-to-many pairs
            "pairs_unique": level1_pairs_unique,   # coverage-level (unique pair_id)
            "pairs_stats": level1_pair_stats,      # coverage diagnostics
            "chem_tensors": chem_tensors_bundle,
            "gene_tensors": gene_tensors_bundle,
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

    print(f">>> Saving torch bundle (atomic) to: {SAVE_PATH}")
    atomic_torch_save(output, SAVE_PATH)
    print(f"✅ Saved torch bundle to: {SAVE_PATH}")

    if SAVE_META_PARQUET:
        unified_df.to_parquet(os.path.join(OUTPUT_DIR, "unified_meta.parquet"), index=False)
        level1_pairs.to_parquet(os.path.join(OUTPUT_DIR, "level1_pairs_meta.parquet"), index=False)
        level1_pairs_unique.to_parquet(os.path.join(OUTPUT_DIR, "level1_pairs_unique.parquet"), index=False)
        level1_pair_stats.to_parquet(os.path.join(OUTPUT_DIR, "level1_pair_stats.parquet"), index=False)

        level2_chem_pool.to_parquet(os.path.join(OUTPUT_DIR, "level2_chem_pool.parquet"), index=False)
        level2_gene_pool.to_parquet(os.path.join(OUTPUT_DIR, "level2_gene_pool.parquet"), index=False)
        level2_pair_stats.to_parquet(os.path.join(OUTPUT_DIR, "level2_pair_stats.parquet"), index=False)

        level3_chem_pool.to_parquet(os.path.join(OUTPUT_DIR, "level3_chem_pool.parquet"), index=False)
        level3_gene_pool.to_parquet(os.path.join(OUTPUT_DIR, "level3_gene_pool.parquet"), index=False)
        level3_pair_stats.to_parquet(os.path.join(OUTPUT_DIR, "level3_pair_stats.parquet"), index=False)

        print("✅ Saved parquet metadata snapshots (recommended).")

    dt = time.time() - t0
    print(f"\n✅ Done. Elapsed: {dt/60:.1f} min")
    if SAVE_TENSORS_SEPARATELY:
        print(f"   Split tensors: {TENSOR_OUT_DIR}")
        if WRITE_TENSOR_MANIFEST:
            print(f"   Manifest: {MANIFEST_PATH}")
