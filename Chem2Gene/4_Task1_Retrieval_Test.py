import os
import gc
import torch
import numpy as np
import pandas as pd
from tqdm import tqdm
from collections import OrderedDict

from functions import load_chem2gen_bundle_from_split

"""
Chem2Gen-Bench — Task 1 Multi-Scenario Retrieval (v3, exact instance-level, OOM-safe)

Scenarios:
  A_L1Base  (Baseline): fixed (cell, src_c, src_g), retrieve Target among L1 targets (L1 pairs_meta)
  B_CellDep (Cell dep): fixed (target, src_c, src_g), retrieve Cell among cells (Level2 pools)
  C_TgtSpec (Tgt spec): fixed (cell, src_c, src_g), retrieve Target among targets (Level3 pools)

Hard requirements satisfied:
  - BG definition identical to Pairwise:
      build from Level3 raw feature (y_delta_*) with fallback
      (src,mod,cell)->(src,mod)->global; mean via float64 sum/count
  - Retrieval is instance-level: every chem instance becomes a query row;
    output keeps Dose/Time/CondID for downstream covariate analysis
  - No sampling/capping is applied (accuracy preserved)
  - OOM-safe: all heavy operations are streaming/batched
"""

# ================= Configuration =================
BASE_DIR = "/mnt/NAS_21T/ProjectData/Chem2Gen"
BENCHMARK_DIR = "/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets"
OUTPUT_DIR = os.path.join(BASE_DIR, "R_Vis_Ready/Task1_Retrieval_MultiScenario")
os.makedirs(OUTPUT_DIR, exist_ok=True)

LINCS_RAW_PATH = "/mnt/NAS_21T/ProjectData/OSMOSIS/processed/LINCS_Processed/LINCS_Engine1_TrainData.pt"
SC_RAW_DIR = "/mnt/NAS_21T/ProjectData/OSMOSIS/processed/scPerturb_Processed"

SEED = 42
MIN_GALLERY_SIZE = 3
MIN_QUERY_SAMPLES = 3
MIN_BG_SAMPLES = 10

RANK_MODE = "mid"  # optimistic | pessimistic | mid

SC_CHUNK_CACHE_SIZE = 2      # small LRU cache to reduce repeated disk I/O
LINCS_BATCH = 8192           # batch size for LINCS index_select
QUERY_BATCH_SOFT = 8192      # batch size for processing query vectors (memory peak control)

KEY_ALIAS = {
    "x": "x_baseline",
    "x_baseline": "x_baseline",
    "y_gene": "y_delta_gene",
    "y_delta_gene": "y_delta_gene",
    "y_path": "y_delta_pathway",
    "y_delta_pathway": "y_delta_pathway",
}

TRACKS = [
    ("Gene",    "y_gene", "y_delta_gene"),
    ("Pathway", "y_path", "y_delta_pathway"),
]

RNG = np.random.default_rng(SEED)

# ================= Null baseline (label permutation) =================
# This is a *chance-level* baseline for retrieval: keep similarity scores fixed,
# but permute the label->index mapping within each gallery so the "correct" index is random.
# IMPORTANT: use a dedicated RNG so the main analysis results are unaffected.
DO_NULL_BASELINE = True     # set False to disable
NULL_SEED = 0
NULL_SAMPLE_MAX = 5000      # reservoir sample size per (Track,Scenario,View,LabelType)
RNG_NULL = np.random.default_rng(NULL_SEED)


# ================= Raw tensor cache =================
_LINCS_CACHE = None

def get_lincs_cache():
    global _LINCS_CACHE
    if _LINCS_CACHE is None:
        if not os.path.exists(LINCS_RAW_PATH):
            raise FileNotFoundError(f"LINCS raw file not found: {LINCS_RAW_PATH}")
        print("    [Cache] Loading LINCS raw dict into RAM (once)...")
        _LINCS_CACHE = torch.load(LINCS_RAW_PATH, map_location="cpu")
    return _LINCS_CACHE

class LRUChunkCache:
    """Tiny LRU cache for scPerturb chunk dicts to reduce repeated disk I/O."""
    def __init__(self, max_items: int = 2):
        self.max_items = int(max_items)
        self._cache = OrderedDict()

    def get(self, chunk_path: str):
        if chunk_path in self._cache:
            self._cache.move_to_end(chunk_path)
            return self._cache[chunk_path]
        obj = torch.load(chunk_path, map_location="cpu")
        self._cache[chunk_path] = obj
        self._cache.move_to_end(chunk_path)
        while len(self._cache) > self.max_items:
            _, old = self._cache.popitem(last=False)
            del old
        return obj

# ================= Basic utilities =================
def robust_cosine_matrix(Q: np.ndarray, G: np.ndarray) -> np.ndarray:
    if Q.ndim != 2 or G.ndim != 2:
        raise ValueError("Q and G must be 2D matrices.")
    if Q.shape[0] == 0 or G.shape[0] == 0:
        return np.zeros((Q.shape[0], G.shape[0]), dtype=np.float32)
    Qn = Q / (np.linalg.norm(Q, axis=1, keepdims=True) + 1e-9)
    Gn = G / (np.linalg.norm(G, axis=1, keepdims=True) + 1e-9)
    return (Qn @ Gn.T).astype(np.float32)

def tie_aware_rank(scores: np.ndarray, true_idx: int, mode: str = "mid") -> float:
    true_score = scores[true_idx]
    gt = int((scores > true_score).sum())
    eq = int(np.isclose(scores, true_score, rtol=1e-5).sum()) - 1
    if eq < 0:
        eq = 0
    if mode == "optimistic":
        return float(gt + 1)
    if mode == "pessimistic":
        return float(gt + eq + 1)
    return float(gt + 1 + (eq / 2.0))

def norm_success(rank: float, n_gallery: int) -> float:
    if n_gallery <= 1:
        return 1.0
    return float(1.0 - (rank - 1.0) / (n_gallery - 1.0))

# ================= Null baseline helpers =================
def _non_identity_perm(n: int, rng: np.random.Generator) -> np.ndarray:
    """Return a random permutation of [0..n-1] that is not the identity (when possible)."""
    if n <= 1:
        return np.arange(n, dtype=np.int64)
    base = np.arange(n, dtype=np.int64)
    perm = rng.permutation(n).astype(np.int64, copy=False)
    # avoid identity (extremely unlikely for n>=3 but make deterministic)
    for _ in range(10):
        if not np.all(perm == base):
            break
        perm = rng.permutation(n).astype(np.int64, copy=False)
    return perm

def _null_state():
    return {"n": 0, "mean": 0.0, "M2": 0.0, "sample": []}

def _null_update_state(state: dict, x: float, rng: np.random.Generator):
    """Online mean/variance + reservoir sample."""
    n = state["n"] + 1
    state["n"] = n
    dx = x - state["mean"]
    state["mean"] += dx / n
    state["M2"] += dx * (x - state["mean"])
    samp = state["sample"]
    if len(samp) < NULL_SAMPLE_MAX:
        samp.append(float(x))
    else:
        j = int(rng.integers(0, n))
        if j < NULL_SAMPLE_MAX:
            samp[j] = float(x)

def _null_finalize(state: dict):
    n = int(state["n"])
    mu = float(state["mean"]) if n > 0 else float("nan")
    sd = float(np.sqrt(state["M2"] / (n - 1))) if n > 1 else 0.0
    samp = np.asarray(state["sample"], dtype=np.float32)
    if samp.size == 0:
        return n, mu, sd, float("nan"), float("nan"), float("nan")
    p05, p50, p95 = np.quantile(samp, [0.05, 0.5, 0.95]).tolist()
    return n, mu, sd, float(p05), float(p50), float(p95)

# (Track, Scenario, View, Direction, LabelType) -> stats dict
NULL_ACC = {}


# ================= Meta extraction =================
def extract_meta_unified(df: pd.DataFrame):
    dose = df.get("dose_val", pd.Series([np.nan]*len(df))).values
    time = df.get("time_val", pd.Series([np.nan]*len(df))).values
    cond = df.get("cond_id", pd.Series(["NA"]*len(df))).values
    return dose, time, cond

def extract_meta_l1_pairs(df: pd.DataFrame, side: str):
    if side == "chem":
        dose = df.get("dose_val_chem", pd.Series([np.nan]*len(df))).values
        time = df.get("time_val_chem", pd.Series([np.nan]*len(df))).values
        cond = df.get("chem_cond_id", pd.Series(["NA"]*len(df))).values
    else:
        dose = df.get("dose_val_gene", pd.Series([np.nan]*len(df))).values
        time = df.get("time_val_gene", pd.Series([np.nan]*len(df))).values
        cond = df.get("gene_cond_id", pd.Series(["NA"]*len(df))).values
    return dose, time, cond

# ================= Streaming aligned fetch (batch) =================
def iter_fetch_aligned_batches(meta_df: pd.DataFrame,
                              raw_feature_key: str,
                              sc_cache: LRUChunkCache | None = None,
                              lincs_batch: int = LINCS_BATCH):
    """
    Yield (pos, tensor_batch) where:
      - pos: np.ndarray of row positions (0..len(meta_df)-1) within meta_df (reset_index assumed)
      - tensor_batch: torch.Tensor (len(pos), D)
    This is OOM-safe and preserves exact row alignment via 'pos'.
    """
    if meta_df is None or meta_df.empty:
        return

    feature_key = KEY_ALIAS.get(raw_feature_key, raw_feature_key)
    meta_df = meta_df.reset_index(drop=True)

    # LINCS
    lincs_mask = meta_df["source_db"].eq("LINCS").values
    if lincs_mask.any():
        lincs = get_lincs_cache()
        if feature_key not in lincs:
            raise KeyError(f"LINCS missing feature key: {feature_key}")
        pos_all = np.flatnonzero(lincs_mask)
        gidx_all = meta_df.loc[lincs_mask, "global_idx"].astype(int).to_numpy()

        for s in range(0, len(pos_all), int(lincs_batch)):
            pos = pos_all[s:s+int(lincs_batch)]
            gidx = gidx_all[s:s+int(lincs_batch)]
            t = lincs[feature_key].index_select(0, torch.as_tensor(gidx, dtype=torch.long))
            yield pos, t

    # scPerturb
    sc_mask = meta_df["source_db"].eq("scPerturb").values
    if sc_mask.any():
        sc_req = meta_df.loc[sc_mask].copy()
        # group by chunk_file to avoid repeated load
        for chunk_file, grp in sc_req.groupby("chunk_file", sort=False):
            path = os.path.join(SC_RAW_DIR, chunk_file)
            if not os.path.exists(path):
                raise FileNotFoundError(f"Missing scPerturb chunk: {path}")
            chunk = torch.load(path, map_location="cpu") if sc_cache is None else sc_cache.get(path)
            if feature_key not in chunk:
                raise KeyError(f"scPerturb chunk missing feature key: {feature_key} ({chunk_file})")
            pos = grp.index.to_numpy()          # already positions in meta_df because meta_df reset_index(drop=True)
            cidx = grp["chunk_idx"].astype(int).to_numpy()
            t = chunk[feature_key].index_select(0, torch.as_tensor(cidx, dtype=torch.long))
            yield pos, t

# ================= Systema BG (same definition as Pairwise) =================
class SystemaBackground:
    """
    Stores hierarchical means:
      - cell mean:   (source_db, modality, cell_std)
      - domain mean: (source_db, modality)
      - global mean

    Uses:
      1) cell mean if N>=MIN_BG_SAMPLES
      2) else domain mean if N>=MIN_BG_SAMPLES
      3) else global mean
    """
    def __init__(self, dim: int):
        self.dim = int(dim)

        self.cell_sum, self.cell_n = {}, {}
        self.dom_sum, self.dom_n = {}, {}

        self.global_sum = np.zeros(self.dim, dtype=np.float64)
        self.global_n = 0

        self.cell_mean, self.dom_mean = {}, {}
        self.global_mean = np.zeros(self.dim, dtype=np.float32)

    def update_from_torch(self, key_cell, key_dom, t: torch.Tensor):
        if t is None or t.numel() == 0:
            return 0
        # sanitize NaN/Inf defensively
        if not torch.isfinite(t).all():
            t = torch.nan_to_num(t, nan=0.0, posinf=0.0, neginf=0.0)

        n = int(t.shape[0])
        s = t.sum(dim=0, dtype=torch.float64).cpu().numpy()

        self.global_sum += s
        self.global_n += n

        if key_dom not in self.dom_sum:
            self.dom_sum[key_dom] = np.zeros(self.dim, dtype=np.float64)
            self.dom_n[key_dom] = 0
        self.dom_sum[key_dom] += s
        self.dom_n[key_dom] += n

        if key_cell not in self.cell_sum:
            self.cell_sum[key_cell] = np.zeros(self.dim, dtype=np.float64)
            self.cell_n[key_cell] = 0
        self.cell_sum[key_cell] += s
        self.cell_n[key_cell] += n

        return n

    def finalize(self):
        self.global_mean = (self.global_sum / max(1, self.global_n)).astype(np.float32)
        for k, s in self.dom_sum.items():
            self.dom_mean[k] = (s / max(1, self.dom_n[k])).astype(np.float32)
        for k, s in self.cell_sum.items():
            self.cell_mean[k] = (s / max(1, self.cell_n[k])).astype(np.float32)

    def get(self, source_db: str, modality: str, cell_std: str) -> np.ndarray:
        k_cell = (source_db, modality, cell_std)
        k_dom = (source_db, modality)

        if k_cell in self.cell_mean and self.cell_n.get(k_cell, 0) >= MIN_BG_SAMPLES:
            return self.cell_mean[k_cell]
        if k_dom in self.dom_mean and self.dom_n.get(k_dom, 0) >= MIN_BG_SAMPLES:
            return self.dom_mean[k_dom]
        return self.global_mean

def build_bg_from_level3_streaming(l3_chem: pd.DataFrame,
                                  l3_gene: pd.DataFrame,
                                  l3_key: str,
                                  dim: int,
                                  sc_cache: LRUChunkCache) -> SystemaBackground:
    """
    Exact BG build, consistent with pairwise:
      - data source: Level3 pools
      - feature: raw l3_key (y_delta_*)
      - hierarchical mean with fallback
    Implementation uses streaming over aligned batches; no sampling/capping.
    """
    print("    Building Systema BG from Level3 (streaming, exact)...")
    bg = SystemaBackground(dim=dim)

    # concat and normalize schema
    c = l3_chem.copy()
    g = l3_gene.copy()
    # ensure modality is present and standardized
    if "modality" not in c.columns:
        c["modality"] = "Chemical"
    if "modality" not in g.columns:
        g["modality"] = "Genetic"

    l3_all = pd.concat([c, g], ignore_index=True)
    # hard guard: modality values
    mod_ok = set(l3_all["modality"].astype(str).unique()).issubset({"Chemical", "Genetic"})
    if not mod_ok:
        raise ValueError(f"Level3 modality contains unexpected values: {l3_all['modality'].unique()}")

    # fast grouping indices
    l3_idx = l3_all.groupby(["source_db", "modality", "cell_std"], sort=False).indices

    # iterate contexts
    for (src, mod, cell), idx in tqdm(l3_idx.items(), total=len(l3_idx), desc="BG contexts"):
        subset = l3_all.iloc[idx][["source_db", "global_idx", "chunk_file", "chunk_idx"]].copy()
        # add back needed columns for fetch
        subset["source_db"] = src
        # NOTE: subset already has correct pointers; source_db column overwritten ok
        n_used = 0
        for pos, t in iter_fetch_aligned_batches(subset, raw_feature_key=l3_key, sc_cache=sc_cache):
            # t is already aligned to 'subset' rows in this batch; for BG we only need sums
            n_used += bg.update_from_torch(
                key_cell=(src, mod, cell),
                key_dom=(src, mod),
                t=t
            )
            del t
            gc.collect()
        # if n_used==0, bg.get will fallback later; keep silent

    bg.finalize()
    print(f"    BG finalized: global_n={bg.global_n}, domains={len(bg.dom_mean)}, cells={len(bg.cell_mean)}")
    return bg

# ================= Streaming centroids (exact) =================
def compute_centroids_streaming(meta_df: pd.DataFrame,
                                label_col: str,
                                raw_feature_key: str,
                                dim: int,
                                sc_cache: LRUChunkCache) -> tuple[list, np.ndarray, np.ndarray]:
    """
    Exact centroids: sums/counts accumulated by label via streaming.
    Returns:
      labels (in first-occurrence order),
      centroids (K,D) float32,
      counts (K,) float64
    """
    if meta_df is None or meta_df.empty:
        return [], np.zeros((0, dim), dtype=np.float32), np.zeros((0,), dtype=np.float64)

    meta_df = meta_df.reset_index(drop=True).copy()
    labels = pd.Index(meta_df[label_col].astype(str)).unique().tolist()
    if len(labels) == 0:
        return [], np.zeros((0, dim), dtype=np.float32), np.zeros((0,), dtype=np.float64)

    # categorical codes
    cat = pd.Categorical(meta_df[label_col].astype(str), categories=labels)
    codes = cat.codes.astype(np.int64)  # 0..K-1
    K = len(labels)

    sums = torch.zeros((K, dim), dtype=torch.float64)
    cnts = torch.zeros((K,), dtype=torch.float64)

    for pos, t in iter_fetch_aligned_batches(meta_df, raw_feature_key=raw_feature_key, sc_cache=sc_cache):
        pos_t = torch.as_tensor(pos, dtype=torch.long)
        code_t = torch.as_tensor(codes[pos], dtype=torch.long)
        tf = t.to(torch.float64)
        sums.index_add_(0, code_t, tf)
        cnts.index_add_(0, code_t, torch.ones((tf.shape[0],), dtype=torch.float64))
        del t, tf, pos_t, code_t
        gc.collect()

    cnts_np = cnts.cpu().numpy()
    sums_np = sums.cpu().numpy()
    cnts_safe = np.maximum(cnts_np, 1.0)
    cents = (sums_np / cnts_safe[:, None]).astype(np.float32)

    return labels, cents, cnts_np

# ================= Scenario runner =================
def run_track(track_name: str, tensor_key: str, l3_key: str, data: dict) -> pd.DataFrame:
    print(f"\n>>> Track={track_name} | tensor_key={tensor_key} | l3_key={l3_key}")

    sc_cache = LRUChunkCache(max_items=SC_CHUNK_CACHE_SIZE)

    # -------- L1 --------
    l1_meta = data["Level1"]["pairs_meta"].copy().reset_index(drop=True)
    l1_chem = data["Level1"]["chem_tensors"][tensor_key]
    l1_gene = data["Level1"]["gene_tensors"][tensor_key]
    dim = int(l1_chem.shape[1])

    # -------- L2/L3 pools --------
    l2_c = data["Level2"]["chem_pool"].copy()
    l2_g = data["Level2"]["gene_pool"].copy()
    l3_c = data["Level3"]["chem_pool"].copy()
    l3_g = data["Level3"]["gene_pool"].copy()

    # BG (exact, consistent with pairwise)
    bg = build_bg_from_level3_streaming(l3_c, l3_g, l3_key=l3_key, dim=dim, sc_cache=sc_cache)
    gc.collect()

    rows = []

    # --------------------------
    # Scenario A: L1 baseline (Target retrieval within cell)
    # --------------------------
    print("    Scenario A: A_L1Base")
    strata = l1_meta.groupby(["source_db_chem", "source_db_gene", "cell_std"], sort=False)

    for (src_c, src_g, cell), sub in tqdm(strata, desc="A_L1Base"):
        # Dedup to avoid many-to-many overweighting
        sub_chem = sub.drop_duplicates("uid_chem", keep="first")
        sub_gene = sub.drop_duplicates("uid_gene", keep="first")

        if len(sub_chem) < MIN_QUERY_SAMPLES:
            continue

        # gallery targets derived from gene side
        gal_targets = sub_gene["target_std"].unique()
        if len(gal_targets) < MIN_GALLERY_SIZE:
            continue

        # Build gene centroids by target (L1 tensor already aligned by L1 row index)
        gene_idx = sub_gene.index.values
        G_vecs = l1_gene[torch.as_tensor(gene_idx, dtype=torch.long)].numpy()
        g_tgts = sub_gene["target_std"].astype(str).values

        tgt_labels = pd.Index(g_tgts).unique().tolist()
        if len(tgt_labels) < MIN_GALLERY_SIZE:
            continue
        tgt_cat = pd.Categorical(g_tgts, categories=tgt_labels)
        codes = tgt_cat.codes.astype(np.int64)
        K = len(tgt_labels)

        sums = np.zeros((K, dim), dtype=np.float64)
        cnts = np.zeros((K,), dtype=np.float64)
        for i in range(G_vecs.shape[0]):
            j = codes[i]
            sums[j] += G_vecs[i]
            cnts[j] += 1.0
        G_cent = (sums / np.maximum(cnts, 1.0)[:, None]).astype(np.float32)
        tgt2i = {t:i for i,t in enumerate(tgt_labels)}
        perm_tgt = _non_identity_perm(len(tgt_labels), RNG_NULL) if DO_NULL_BASELINE else None
        # Queries: chem instances (dedup)
        chem_idx = sub_chem.index.values
        Q_vecs = l1_chem[torch.as_tensor(chem_idx, dtype=torch.long)].numpy().astype(np.float32)
        q_tgts = sub_chem["target_std"].astype(str).values
        dose, time, cond = extract_meta_l1_pairs(sub_chem, side="chem")

        # BG
        bg_q = bg.get(src_c, "Chemical", cell)
        bg_g = bg.get(src_g, "Genetic", cell)

        for view in ["Standard", "Systema"]:
            Qm = (Q_vecs - bg_q).astype(np.float32) if view == "Systema" else Q_vecs
            Gm = (G_cent - bg_g).astype(np.float32) if view == "Systema" else G_cent
            sims = robust_cosine_matrix(Qm, Gm)

            for i in range(len(sub_chem)):
                true_t = q_tgts[i]
                if true_t not in tgt2i:
                    continue
                true_idx = tgt2i[true_t]
                r = tie_aware_rank(sims[i], true_idx, mode=RANK_MODE)
                # Null baseline (single label permutation within gallery)
                if DO_NULL_BASELINE and perm_tgt is not None:
                    true_idx_perm = int(perm_tgt[true_idx])
                    r0 = tie_aware_rank(sims[i], true_idx_perm, mode=RANK_MODE)
                    key0 = (track_name, "A_L1Base", view, "Chem->Gene", "Target")
                    acc0 = NULL_ACC.get(key0)
                    if acc0 is None:
                        acc0 = {"rank": _null_state(), "mrr": _null_state(), "success": _null_state()}
                        NULL_ACC[key0] = acc0
                    _null_update_state(acc0["rank"], r0, RNG_NULL)
                    _null_update_state(acc0["mrr"], 1.0 / r0, RNG_NULL)
                    _null_update_state(acc0["success"], norm_success(r0, len(tgt_labels)), RNG_NULL)

                rows.append({
                    "Track": track_name,
                    "Scenario": "A_L1Base",
                    "View": view,
                    "Direction": "Chem->Gene",
                    "LabelType": "Target",
                    "Cell": cell,
                    "Target": true_t,
                    "Source_Chem": src_c,
                    "Source_Gene": src_g,
                    "Dose": dose[i], "Time": time[i], "CondID": cond[i],
                    "True_Rank": r,
                    "MRR": 1.0 / r,
                    "Success_Score": norm_success(r, len(tgt_labels)),
                    "N_Gallery": len(tgt_labels),
                })

        del G_vecs, Q_vecs, sims
        gc.collect()

    # --------------------------
    # Scenario B: Cell dependence (Cell retrieval across cells for fixed target)
    # --------------------------
    print("    Scenario B: B_CellDep")
    idx_c_by_tgt = l2_c.groupby("target_std", sort=False).indices
    idx_g_by_tgt = l2_g.groupby("target_std", sort=False).indices

    for tgt in tqdm(list(idx_c_by_tgt.keys()), desc="B_CellDep Targets"):
        if tgt not in idx_g_by_tgt:
            continue
        c_all = l2_c.iloc[idx_c_by_tgt[tgt]]
        g_all = l2_g.iloc[idx_g_by_tgt[tgt]]
        if c_all.empty or g_all.empty:
            continue

        for sc in c_all["source_db"].unique():
            for sg in g_all["source_db"].unique():
                curr_c = c_all[c_all["source_db"] == sc].copy()
                curr_g = g_all[g_all["source_db"] == sg].copy()
                if curr_c.empty or curr_g.empty:
                    continue

                gal_cells = curr_g["cell_std"].astype(str).unique().tolist()
                if len(gal_cells) < MIN_GALLERY_SIZE:
                    continue

                # Gene centroids per cell (exact streaming)
                cell_labels, G_cent, g_cnts = compute_centroids_streaming(
                    curr_g[["source_db", "global_idx", "chunk_file", "chunk_idx", "cell_std"]].copy(),
                    label_col="cell_std",
                    raw_feature_key=tensor_key,
                    dim=dim,
                    sc_cache=sc_cache
                )
                if len(cell_labels) < MIN_GALLERY_SIZE:
                    continue
                cell2i = {c:i for i,c in enumerate(cell_labels)}
                perm_cell = _non_identity_perm(len(cell_labels), RNG_NULL) if DO_NULL_BASELINE else None
                # Queries: chem instances, only those whose cell is in gallery
                curr_c = curr_c[curr_c["cell_std"].astype(str).isin(cell_labels)].copy()
                if len(curr_c) < MIN_QUERY_SAMPLES:
                    continue

                # Precompute BG for gallery labels (gene side) once
                bg_g_labels = np.stack([bg.get(sg, "Genetic", c) for c in cell_labels], axis=0).astype(np.float32)

                # For queries, BG depends on each query cell; precompute map
                uniq_q_cells = curr_c["cell_std"].astype(str).unique().tolist()
                bg_q_map = {c: bg.get(sc, "Chemical", c).astype(np.float32) for c in uniq_q_cells}

                # Iterate query vectors in streaming batches (exact)
                curr_c_min = curr_c[["source_db", "global_idx", "chunk_file", "chunk_idx", "cell_std", "dose_val", "time_val", "cond_id"]].copy()
                curr_c_min = curr_c_min.reset_index(drop=True)

                for pos, tQ in iter_fetch_aligned_batches(curr_c_min, raw_feature_key=tensor_key, sc_cache=sc_cache, lincs_batch=LINCS_BATCH):
                    # batch meta
                    subm = curr_c_min.iloc[pos]
                    q_cells = subm["cell_std"].astype(str).values
                    dose, time, cond = extract_meta_unified(subm)

                    Q = tQ.numpy().astype(np.float32, copy=False)
                    del tQ

                    # BG for queries (per-row)
                    bg_q = np.stack([bg_q_map[c] for c in q_cells], axis=0).astype(np.float32)

                    for view in ["Standard", "Systema"]:
                        Qm = (Q - bg_q).astype(np.float32) if view == "Systema" else Q
                        Gm = (G_cent - bg_g_labels).astype(np.float32) if view == "Systema" else G_cent
                        sims = robust_cosine_matrix(Qm, Gm)

                        for i in range(len(q_cells)):
                            true_idx = cell2i[q_cells[i]]
                            r = tie_aware_rank(sims[i], true_idx, mode=RANK_MODE)
                            # Null baseline (single label permutation within gallery)
                            if DO_NULL_BASELINE and perm_cell is not None:
                                true_idx_perm = int(perm_cell[true_idx])
                                r0 = tie_aware_rank(sims[i], true_idx_perm, mode=RANK_MODE)
                                key0 = (track_name, "B_CellDep", view, "Chem->Gene", "Cell")
                                acc0 = NULL_ACC.get(key0)
                                if acc0 is None:
                                    acc0 = {"rank": _null_state(), "mrr": _null_state(), "success": _null_state()}
                                    NULL_ACC[key0] = acc0
                                _null_update_state(acc0["rank"], r0, RNG_NULL)
                                _null_update_state(acc0["mrr"], 1.0 / r0, RNG_NULL)
                                _null_update_state(acc0["success"], norm_success(r0, len(cell_labels)), RNG_NULL)

                            rows.append({
                                "Track": track_name,
                                "Scenario": "B_CellDep",
                                "View": view,
                                "Direction": "Chem->Gene",
                                "LabelType": "Cell",
                                "Cell": q_cells[i],
                                "Target": str(tgt),
                                "Source_Chem": sc,
                                "Source_Gene": sg,
                                "Dose": dose[i], "Time": time[i], "CondID": cond[i],
                                "True_Rank": r,
                                "MRR": 1.0 / r,
                                "Success_Score": norm_success(r, len(cell_labels)),
                                "N_Gallery": len(cell_labels),
                            })

                        del sims
                        gc.collect()

                    del Q, bg_q
                    gc.collect()

                del curr_c, curr_g, G_cent
                gc.collect()

    # --------------------------
    # Scenario C: Target specificity (Target retrieval within cell)
    # --------------------------
    print("    Scenario C: C_TgtSpec")
    idx_c_by_cell = l3_c.groupby("cell_std", sort=False).indices
    idx_g_by_cell = l3_g.groupby("cell_std", sort=False).indices

    for cell in tqdm(list(idx_c_by_cell.keys()), desc="C_TgtSpec Cells"):
        if cell not in idx_g_by_cell:
            continue
        c_all = l3_c.iloc[idx_c_by_cell[cell]]
        g_all = l3_g.iloc[idx_g_by_cell[cell]]
        if c_all.empty or g_all.empty:
            continue

        for sc in c_all["source_db"].unique():
            for sg in g_all["source_db"].unique():
                curr_c = c_all[c_all["source_db"] == sc].copy()
                curr_g = g_all[g_all["source_db"] == sg].copy()
                if curr_c.empty or curr_g.empty:
                    continue

                gal_tgts = curr_g["target_std"].astype(str).unique().tolist()
                if len(gal_tgts) < 10:
                    continue

                # Gene centroids per target (exact streaming)
                tgt_labels, G_cent, _ = compute_centroids_streaming(
                    curr_g[["source_db", "global_idx", "chunk_file", "chunk_idx", "target_std"]].copy(),
                    label_col="target_std",
                    raw_feature_key=tensor_key,
                    dim=dim,
                    sc_cache=sc_cache
                )
                if len(tgt_labels) < 10:
                    continue
                tgt2i = {t:i for i,t in enumerate(tgt_labels)}
                perm_tgt = _non_identity_perm(len(tgt_labels), RNG_NULL) if DO_NULL_BASELINE else None
                # Queries: chem instances whose target in gallery
                curr_c = curr_c[curr_c["target_std"].astype(str).isin(tgt_labels)].copy()
                if len(curr_c) < MIN_QUERY_SAMPLES:
                    continue

                # BG constants (fixed cell)
                bg_q_const = bg.get(sc, "Chemical", str(cell)).astype(np.float32)
                bg_g_const = bg.get(sg, "Genetic", str(cell)).astype(np.float32)

                curr_c_min = curr_c[["source_db", "global_idx", "chunk_file", "chunk_idx", "target_std", "dose_val", "time_val", "cond_id"]].copy()
                curr_c_min = curr_c_min.reset_index(drop=True)

                for pos, tQ in iter_fetch_aligned_batches(curr_c_min, raw_feature_key=tensor_key, sc_cache=sc_cache, lincs_batch=LINCS_BATCH):
                    subm = curr_c_min.iloc[pos]
                    q_tgts = subm["target_std"].astype(str).values
                    dose, time, cond = extract_meta_unified(subm)

                    Q = tQ.numpy().astype(np.float32, copy=False)
                    del tQ

                    for view in ["Standard", "Systema"]:
                        Qm = (Q - bg_q_const).astype(np.float32) if view == "Systema" else Q
                        Gm = (G_cent - bg_g_const).astype(np.float32) if view == "Systema" else G_cent
                        sims = robust_cosine_matrix(Qm, Gm)

                        for i in range(len(q_tgts)):
                            true_t = q_tgts[i]
                            true_idx = tgt2i[true_t]
                            r = tie_aware_rank(sims[i], true_idx, mode=RANK_MODE)
                            # Null baseline (single label permutation within gallery)
                            if DO_NULL_BASELINE and perm_tgt is not None:
                                true_idx_perm = int(perm_tgt[true_idx])
                                r0 = tie_aware_rank(sims[i], true_idx_perm, mode=RANK_MODE)
                                key0 = (track_name, "C_TgtSpec", view, "Chem->Gene", "Target")
                                acc0 = NULL_ACC.get(key0)
                                if acc0 is None:
                                    acc0 = {"rank": _null_state(), "mrr": _null_state(), "success": _null_state()}
                                    NULL_ACC[key0] = acc0
                                _null_update_state(acc0["rank"], r0, RNG_NULL)
                                _null_update_state(acc0["mrr"], 1.0 / r0, RNG_NULL)
                                _null_update_state(acc0["success"], norm_success(r0, len(tgt_labels)), RNG_NULL)

                            rows.append({
                                "Track": track_name,
                                "Scenario": "C_TgtSpec",
                                "View": view,
                                "Direction": "Chem->Gene",
                                "LabelType": "Target",
                                "Cell": str(cell),
                                "Target": true_t,
                                "Source_Chem": sc,
                                "Source_Gene": sg,
                                "Dose": dose[i], "Time": time[i], "CondID": cond[i],
                                "True_Rank": r,
                                "MRR": 1.0 / r,
                                "Success_Score": norm_success(r, len(tgt_labels)),
                                "N_Gallery": len(tgt_labels),
                            })

                        del sims
                        gc.collect()

                    del Q
                    gc.collect()

                del curr_c, curr_g, G_cent
                gc.collect()

    return pd.DataFrame(rows)

# ================= Main =================
if __name__ == "__main__":
    print(">>> Loading data bundle (split)...")
    data = load_chem2gen_bundle_from_split(
        BENCHMARK_DIR,
        load_unified_meta=True,
        tensor_cache_max_items=2,
        strict=True,
    )

    all_df = []
    for track_name, tensor_key, l3_key in TRACKS:
        df = run_track(track_name, tensor_key, l3_key, data)
        if not df.empty:
            all_df.append(df)
        gc.collect()

    if not all_df:
        print("[WARN] No results generated.")
        raise SystemExit(0)

    master = pd.concat(all_df, ignore_index=True)
    out_path = os.path.join(OUTPUT_DIR, "Task1_Retrieval_MultiScenario_PerQuery.csv")
    master.to_csv(out_path, index=False)
    print(f"\n✅ Saved {len(master)} rows to:\n  {out_path}")

    print("\n=== Mean Success_Score ===")
    print(master.groupby(["Track","Scenario","View"])["Success_Score"].mean().sort_values(ascending=False))


    # ==========================
    # Null baseline summary output (does NOT affect main outputs)
    # ==========================
    if DO_NULL_BASELINE and NULL_ACC:
        null_rows = []
        for (trk, scn, view, direc, lab), acc in NULL_ACC.items():
            n_r, mu_r, sd_r, p05_r, p50_r, p95_r = _null_finalize(acc["rank"])
            n_m, mu_m, sd_m, p05_m, p50_m, p95_m = _null_finalize(acc["mrr"])
            n_s, mu_s, sd_s, p05_s, p50_s, p95_s = _null_finalize(acc["success"])
            # counts should match across metrics; keep rank count as authoritative
            null_rows.append({
                "Track": trk,
                "Scenario": scn,
                "View": view,
                "Direction": direc,
                "LabelType": lab,
                "null_n": n_r,
                "rank_mu": mu_r, "rank_sd": sd_r, "rank_p05": p05_r, "rank_p50": p50_r, "rank_p95": p95_r,
                "mrr_mu": mu_m,  "mrr_sd": sd_m,  "mrr_p05": p05_m,  "mrr_p50": p50_m,  "mrr_p95": p95_m,
                "success_mu": mu_s, "success_sd": sd_s, "success_p05": p05_s, "success_p50": p50_s, "success_p95": p95_s,
            })

        null_df = pd.DataFrame(null_rows).sort_values(["Track","Scenario","View","LabelType"])
        null_out = os.path.join(OUTPUT_DIR, "Task1_Retrieval_MultiScenario_NullSummary.csv")
        null_df.to_csv(null_out, index=False)
        print(f"\n✅ Saved null baseline summary to:\n  {null_out}")
        print("\n[Null] Mean Success_Score (null) by Track/Scenario/View:")
        print(null_df.groupby(["Track","Scenario","View"])["success_mu"].mean().sort_values(ascending=False))

