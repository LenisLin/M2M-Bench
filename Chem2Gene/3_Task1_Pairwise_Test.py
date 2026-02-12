import os
import gc
import torch
import pandas as pd
import numpy as np
from tqdm import tqdm
from collections import OrderedDict
from scipy.spatial.distance import pdist, cdist
from sklearn.metrics.pairwise import cosine_similarity
from scipy.stats import pearsonr

from functions import load_chem2gen_bundle_from_split

"""
Chem2Gen-Bench — Task 1 Pairwise Test (Systema-aligned) [Hardened v1]

Changes vs your original:
  [B] Numerical robustness
    - norm_ratio: denominator ~0 -> NaN (not 0.0)
    - finite checks: BG + clouds are nan_to_num’d to avoid NaN propagation
    - null diagnostics: null_n, null_valid, etc.

  [C] Performance & scalability
    - Level3 pool indexing via groupby().indices to avoid repeated full-table filters
    - scPerturb chunk LRU cache (tiny) to reduce repeated disk I/O
    - BG building streams tensors per chunk (no torch.cat for BG mean)

Outputs:
  - Task1_Pairwise_Metrics_Wide.csv
  - Task1_Pairwise_Metrics_Long.csv
  - Task1_Systema_Background_Stats.csv
"""

# ================= Configuration =================
BASE_DIR = "/mnt/NAS_21T/ProjectData/Chem2Gen"
BENCHMARK_DIR = "/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets"
OUTPUT_DIR = os.path.join(BASE_DIR, "R_Vis_Ready/Task1_Metrics")
os.makedirs(OUTPUT_DIR, exist_ok=True)

LINCS_RAW_PATH = "/mnt/NAS_21T/ProjectData/OSMOSIS/processed/LINCS_Processed/LINCS_Engine1_TrainData.pt"
SC_RAW_DIR = "/mnt/NAS_21T/ProjectData/OSMOSIS/processed/scPerturb_Processed"

SEED = 42
N_BOOTSTRAP = 20
MAX_SAMPLE_N = 100
MIN_UNIQUE_PER_SIDE = 3
MIN_BG_SAMPLES = 10
N_NULL = 200
MAX_BG_PER_CONTEXT = None   # e.g. 5000 / 10000 as safety; None = use all rows
SC_CHUNK_CACHE_SIZE = 2     # small LRU cache for scPerturb chunks
EPS_NORM = 1e-9             # for norm ratio stability

RNG = np.random.default_rng(SEED)

TRACKS = [
    ("Gene", "y_gene", "y_delta_gene", 50),
    ("Pathway", "y_path", "y_delta_pathway", 10),
]

# ================= Lazy Tensor Cache =================
_LINCS_CACHE = None

def get_lincs_cache():
    global _LINCS_CACHE
    if _LINCS_CACHE is None:
        print("    [Cache] Loading LINCS raw tensors into RAM (once)...")
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
        # load
        obj = torch.load(chunk_path, map_location="cpu")
        self._cache[chunk_path] = obj
        self._cache.move_to_end(chunk_path)
        # evict
        while len(self._cache) > self.max_items:
            _, old = self._cache.popitem(last=False)
            # allow GC
            del old
        return obj

def _torch_sanitize_inplace(x: torch.Tensor) -> torch.Tensor:
    """Replace NaN/Inf with 0 to prevent propagation. Returns tensor (may share storage)."""
    if not torch.isfinite(x).all():
        x = torch.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    return x

def _np_sanitize(x: np.ndarray) -> np.ndarray:
    """Replace NaN/Inf with 0 to prevent propagation."""
    if not np.isfinite(x).all():
        x = np.nan_to_num(x, nan=0.0, posinf=0.0, neginf=0.0)
    return x

# ================= Streaming Tensor Fetch (for BG) =================
def iter_fetch_tensors_dynamic(meta_df: pd.DataFrame,
                              feature_key: str,
                              sc_cache: LRUChunkCache | None = None):
    """
    Yield tensors in chunks (LINCS once, scPerturb per chunk_file).
    NOTE: This is for streaming aggregation; ordering is not guaranteed vs meta_df.
    """
    if meta_df is None or meta_df.empty:
        return

    lincs_req = meta_df[meta_df["source_db"] == "LINCS"]
    sc_req = meta_df[meta_df["source_db"] == "scPerturb"]

    # LINCS
    if not lincs_req.empty:
        lincs_data = get_lincs_cache()
        if feature_key not in lincs_data:
            raise KeyError(f"LINCS cache missing key: {feature_key}")
        idx = torch.as_tensor(lincs_req["global_idx"].values, dtype=torch.long)
        t = lincs_data[feature_key][idx]
        yield t

    # scPerturb
    if not sc_req.empty:
        # group by chunk for minimal loads
        for chunk_file, group in sc_req.sort_values("chunk_file").groupby("chunk_file"):
            path = os.path.join(SC_RAW_DIR, chunk_file)
            if not os.path.exists(path):
                continue
            if sc_cache is None:
                chunk = torch.load(path, map_location="cpu")
            else:
                chunk = sc_cache.get(path)
            if feature_key not in chunk:
                raise KeyError(f"scPerturb chunk missing key: {feature_key} ({chunk_file})")
            idx = torch.as_tensor(group["chunk_idx"].values, dtype=torch.long)
            t = chunk[feature_key][idx]
            yield t

def fetch_tensors_dynamic(meta_df: pd.DataFrame,
                          feature_key: str,
                          sc_cache: LRUChunkCache | None = None) -> torch.Tensor | None:
    """Non-streaming convenience wrapper (kept for compatibility)."""
    parts = []
    for t in iter_fetch_tensors_dynamic(meta_df, feature_key, sc_cache=sc_cache):
        parts.append(t)
    if not parts:
        return None
    return torch.cat(parts, dim=0)

# ================= Metrics =================
def robust_cosine(u: np.ndarray, v: np.ndarray) -> float:
    u = np.asarray(u).reshape(1, -1)
    v = np.asarray(v).reshape(1, -1)
    nu = float(np.linalg.norm(u))
    nv = float(np.linalg.norm(v))
    if nu == 0.0 or nv == 0.0:
        return 0.0
    return float(cosine_similarity(u, v)[0, 0])

def compute_energy_distance_boot(X: np.ndarray, Y: np.ndarray,
                                rng: np.random.Generator,
                                n_boot: int = N_BOOTSTRAP,
                                max_n: int = MAX_SAMPLE_N):
    """Bootstrapped Energy Distance between two point clouds X and Y. Returns (mean, std)."""
    X = np.asarray(X)
    Y = np.asarray(Y)
    nx, ny = X.shape[0], Y.shape[0]

    n_sample = int(min(nx, ny, max_n))
    if n_sample < 2:
        return np.nan, np.nan

    dists = []
    for _ in range(int(n_boot)):
        ix = rng.choice(nx, size=n_sample, replace=(nx < n_sample))
        iy = rng.choice(ny, size=n_sample, replace=(ny < n_sample))

        sx = X[ix]
        sy = Y[iy]

        d_xy = float(np.mean(cdist(sx, sy, metric="euclidean")))
        d_xx = float(np.mean(pdist(sx, metric="euclidean")))
        d_yy = float(np.mean(pdist(sy, metric="euclidean")))

        e = 2.0 * d_xy - d_xx - d_yy
        dists.append(max(0.0, float(e)))  # numerical clamp

    return float(np.mean(dists)), float(np.std(dists))

def _jaccard(a_idx: np.ndarray, b_idx: np.ndarray) -> float:
    inter = np.intersect1d(a_idx, b_idx).size
    uni = np.union1d(a_idx, b_idx).size
    return float(inter / uni) if uni > 0 else np.nan

def signal_metrics_union_topk(vec_a: np.ndarray, vec_b: np.ndarray, top_k: int):
    """Signal metrics on union of TopK |z| features + signed overlap."""
    a = np.asarray(vec_a).flatten()
    b = np.asarray(vec_b).flatten()

    k = int(min(top_k, a.size, b.size))
    if k < 2:
        return np.nan, np.nan, np.nan, np.nan, np.nan

    idx_a = np.argsort(np.abs(a))[-k:]
    idx_b = np.argsort(np.abs(b))[-k:]
    union = np.union1d(idx_a, idx_b)

    if union.size >= 2:
        try:
            pcc = float(pearsonr(a[union], b[union])[0])
        except Exception:
            pcc = np.nan
    else:
        pcc = np.nan

    sa = np.sign(a[union])
    sb = np.sign(b[union])
    valid = (sa != 0) & (sb != 0)
    acc = float((sa[valid] == sb[valid]).mean()) if valid.sum() > 0 else np.nan

    jac_abs = _jaccard(idx_a, idx_b)

    idx_up_a = np.argsort(a)[-k:]
    idx_up_b = np.argsort(b)[-k:]
    idx_dn_a = np.argsort(a)[:k]
    idx_dn_b = np.argsort(b)[:k]

    jac_up = _jaccard(idx_up_a, idx_up_b)
    jac_dn = _jaccard(idx_dn_a, idx_dn_b)

    return pcc, acc, jac_abs, jac_up, jac_dn

# ================= Systema Background =================
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
        self.global_mean = np.zeros(self.dim, dtype=np.float64)

    def update_from_torch(self, key_cell, key_dom, t: torch.Tensor):
        """Streaming update using torch sums in float64, no full mat float64 conversion."""
        if t is None or t.numel() == 0:
            return 0
        t = _torch_sanitize_inplace(t)

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
        self.global_mean = self.global_sum / max(1, self.global_n)
        for k, s in self.dom_sum.items():
            self.dom_mean[k] = s / max(1, self.dom_n[k])
        for k, s in self.cell_sum.items():
            self.cell_mean[k] = s / max(1, self.cell_n[k])

    def get(self, source_db: str, modality: str, cell_std: str):
        k_cell = (source_db, modality, cell_std)
        k_dom = (source_db, modality)

        if k_cell in self.cell_mean and self.cell_n.get(k_cell, 0) >= MIN_BG_SAMPLES:
            return self.cell_mean[k_cell], "cell", int(self.cell_n[k_cell])

        if k_dom in self.dom_mean and self.dom_n.get(k_dom, 0) >= MIN_BG_SAMPLES:
            return self.dom_mean[k_dom], "domain", int(self.dom_n[k_dom])

        return self.global_mean, "global", int(self.global_n)

# ================= Utilities =================
def ensure_uid(df: pd.DataFrame, side: str) -> pd.DataFrame:
    """Ensure df has uid_{side}. If absent, build from pointers."""
    col = f"uid_{side}"
    if col in df.columns:
        return df

    src = df[f"source_db_{side}"].astype(str)

    raw_gidx = df.get(f"global_idx_{side}", pd.Series([np.nan] * len(df)))
    gidx = raw_gidx.fillna(-1).astype(int).astype(str)

    cfile = df.get(f"chunk_file_{side}", pd.Series(["NA"] * len(df))).astype(str)
    raw_cidx = df.get(f"chunk_idx_{side}", pd.Series([np.nan] * len(df)))
    cidx = raw_cidx.fillna(-1).astype(int).astype(str)

    cond_lincs = (src == "LINCS")
    uid_arr = np.empty(len(df), dtype=object)
    uid_arr[cond_lincs] = (src[cond_lincs] + ":" + gidx[cond_lincs]).to_numpy()
    uid_arr[~cond_lincs] = (src[~cond_lincs] + ":" + cfile[~cond_lincs] + ":" + cidx[~cond_lincs]).to_numpy()

    df[col] = uid_arr
    return df

def _centroids_for_stratum(l1_meta: pd.DataFrame,
                           cell: str,
                           src_c: str,
                           src_g: str,
                           bulk_t: str,
                           l1_chem: torch.Tensor,
                           l1_gene: torch.Tensor):
    """Build per-target centroids for a (cell,src_c,src_g,bulk) stratum."""
    sub = l1_meta[
        (l1_meta["cell_std"] == cell)
        & (l1_meta["source_db_chem"] == src_c)
        & (l1_meta["source_db_gene"] == src_g)
        & (l1_meta["bulk_type"] == bulk_t)
    ]
    if sub.empty:
        return np.array([]), np.zeros((0, l1_chem.shape[1]), dtype=np.float32), np.zeros((0, l1_gene.shape[1]), dtype=np.float32)

    c_dedup = sub.drop_duplicates("uid_chem")
    g_dedup = sub.drop_duplicates("uid_gene")

    c_idx = torch.as_tensor(c_dedup.index.to_numpy(), dtype=torch.long)
    g_idx = torch.as_tensor(g_dedup.index.to_numpy(), dtype=torch.long)

    c_vecs = l1_chem[c_idx].cpu().numpy()
    g_vecs = l1_gene[g_idx].cpu().numpy()
    c_vecs = _np_sanitize(c_vecs)
    g_vecs = _np_sanitize(g_vecs)

    def build(vecs: np.ndarray, targets: np.ndarray):
        uniq, codes = np.unique(targets, return_inverse=True)
        sums = np.zeros((uniq.size, vecs.shape[1]), dtype=np.float64)
        cnts = np.zeros((uniq.size,), dtype=np.float64)
        np.add.at(sums, codes, vecs)
        np.add.at(cnts, codes, 1.0)
        cents = (sums / cnts[:, None]).astype(np.float32)
        return uniq, cents

    t_c, c_cent = build(c_vecs, c_dedup["target_std"].values)
    t_g, g_cent = build(g_vecs, g_dedup["target_std"].values)

    common = np.intersect1d(t_c, t_g)
    if common.size == 0:
        return np.array([]), np.zeros((0, l1_chem.shape[1]), dtype=np.float32), np.zeros((0, l1_gene.shape[1]), dtype=np.float32)

    pos_c = {t: i for i, t in enumerate(t_c)}
    pos_g = {t: i for i, t in enumerate(t_g)}
    idx_c = np.asarray([pos_c[t] for t in common], dtype=int)
    idx_g = np.asarray([pos_g[t] for t in common], dtype=int)

    return common, c_cent[idx_c], g_cent[idx_g]

def _null_cosine_dist(cent_chem: np.ndarray,
                      cent_gene: np.ndarray,
                      rng: np.random.Generator,
                      n_null: int = N_NULL) -> np.ndarray:
    """Sample a null cosine distribution from mismatched target centroids."""
    n = int(min(cent_chem.shape[0], cent_gene.shape[0]))
    if n < 2:
        return np.array([], dtype=np.float32)

    out = []
    for _ in range(int(n_null)):
        i = int(rng.integers(0, n))
        j = int(rng.integers(0, n - 1))
        if j >= i:
            j += 1
        out.append(robust_cosine(cent_chem[i], cent_gene[j]))
    return np.asarray(out, dtype=np.float32)

def _zscore(x: float, arr: np.ndarray):
    """Return (z, mu, sd, n, valid)."""
    if arr is None or arr.size < 5:
        return np.nan, np.nan, np.nan, int(arr.size) if arr is not None else 0, False
    mu = float(np.mean(arr))
    sd = float(np.std(arr))
    valid = bool(sd > 0)
    z = float((x - mu) / sd) if valid else np.nan
    return z, mu, sd, int(arr.size), valid



def _summarize_null(arr: np.ndarray) -> dict:
    """Return summary stats for a null sample array (mean, sd, p05/p50/p95)."""
    if arr is None or arr.size == 0:
        return {
            "n": 0,
            "mu": np.nan,
            "sd": np.nan,
            "p05": np.nan,
            "p50": np.nan,
            "p95": np.nan,
        }
    arr = np.asarray(arr, dtype=np.float64)
    return {
        "n": int(arr.size),
        "mu": float(np.mean(arr)),
        "sd": float(np.std(arr)),
        "p05": float(np.quantile(arr, 0.05)),
        "p50": float(np.quantile(arr, 0.50)),
        "p95": float(np.quantile(arr, 0.95)),
    }


# ================= Main =================
if __name__ == "__main__":
    print(">>> [Task 1] Pairwise gap quantification (hardened)")
    data = load_chem2gen_bundle_from_split(
    BENCHMARK_DIR,
    load_unified_meta=True,
    tensor_cache_max_items=2,
    strict=True,
)

    # ---- Level1 meta ----
    l1_meta = data["Level1"]["pairs_meta"].copy().reset_index(drop=True)
    l1_meta = ensure_uid(l1_meta, "chem")
    l1_meta = ensure_uid(l1_meta, "gene")

    l1_meta["domain_type"] = np.where(
        l1_meta["source_db_chem"] == l1_meta["source_db_gene"],
        "SameDomain",
        "CrossDomain",
    )
    l1_meta["bulk_type"] = l1_meta["is_bulk_chem"].astype(str) + "->" + l1_meta["is_bulk_gene"].astype(str)

    # ---- Level3 pools ----
    l3_chem = data["Level3"]["chem_pool"].copy()
    l3_gene = data["Level3"]["gene_pool"].copy()
    l3_chem["Modality"] = "Chemical"
    l3_gene["Modality"] = "Genetic"
    l3_all = pd.concat([l3_chem, l3_gene], ignore_index=True)

    # schema checks
    assert {"source_db", "cell_std", "Modality"}.issubset(l3_all.columns)
    if (l3_all["source_db"] == "LINCS").any():
        assert "global_idx" in l3_all.columns
    if (l3_all["source_db"] == "scPerturb").any():
        assert {"chunk_file", "chunk_idx"}.issubset(l3_all.columns)

    # C1: build fast index for (source_db, Modality, cell_std) -> row indices
    print(">>> Building Level3 group indices...")
    l3_idx = l3_all.groupby(["source_db", "Modality", "cell_std"], sort=False).indices

    # contexts needed for BG: only those used in Level1
    l1_needed = pd.concat(
        [
            l1_meta[["source_db_chem", "cell_std"]].rename(columns={"source_db_chem": "source_db"}).assign(Modality="Chemical"),
            l1_meta[["source_db_gene", "cell_std"]].rename(columns={"source_db_gene": "source_db"}).assign(Modality="Genetic"),
        ],
        ignore_index=True,
    ).drop_duplicates()

    all_results = []
    all_bg_stats = []

    # ---- Null baseline exports (one-run permutation-style baseline via mismatched target sampling) ----
    # These lists collect ONLY per-stratum null samples, and do not affect the main Task1 outputs.
    all_null_draws = []           # long format: one row per null draw
    all_null_stratum_stats = []   # one row per stratum summary

    sc_cache = LRUChunkCache(max_items=SC_CHUNK_CACHE_SIZE)

    for track_name, l1_key, l3_key, topk_deg in TRACKS:
        print(f"\n>>> [Track] {track_name} | L1={l1_key} | L3={l3_key} | TopK={topk_deg}")

        # ---- Level1 tensors for this track ----
        l1_chem = data["Level1"]["chem_tensors"].get(l1_key, None)
        l1_gene = data["Level1"]["gene_tensors"].get(l1_key, None)
        if l1_chem is None or l1_gene is None:
            print(f"  [Skip] Missing Level1 tensors for {track_name} ({l1_key}).")
            continue

        assert len(l1_meta) == l1_chem.shape[0] == l1_gene.shape[0], "L1 meta/tensor mismatch"
        dim = int(l1_chem.shape[1])

        # ---- Systema background (streaming) ----
        print(">>> Building Systema backgrounds from Level3 (streaming)...")
        bg = SystemaBackground(dim=dim)

        for _, row in tqdm(l1_needed.iterrows(), total=len(l1_needed), desc=f"BG contexts ({track_name})"):
            s_db, mod, cell = row["source_db"], row["Modality"], row["cell_std"]
            key = (s_db, mod, cell)
            if key not in l3_idx:
                continue

            idx = l3_idx[key]
            subset = l3_all.iloc[idx]

            # C3: subsample huge contexts at meta level to cap fetch
            if MAX_BG_PER_CONTEXT is not None and len(subset) > int(MAX_BG_PER_CONTEXT):
                subset = subset.sample(n=int(MAX_BG_PER_CONTEXT), random_state=SEED)

            n_used_total = 0
            # C2+C3: stream per chunk, update bg by sum/n
            for t in iter_fetch_tensors_dynamic(subset, feature_key=l3_key, sc_cache=sc_cache):
                if t is None or t.numel() == 0:
                    continue
                n_used = bg.update_from_torch(
                    key_cell=(s_db, mod, cell),
                    key_dom=(s_db, mod),
                    t=t
                )
                n_used_total += int(n_used)
                del t
                gc.collect()

            if n_used_total > 0:
                all_bg_stats.append(
                    {
                        "Track": track_name,
                        "source_db": s_db,
                        "Modality": mod,
                        "cell_std": cell,
                        "n_used": int(n_used_total),
                    }
                )

        bg.finalize()
        print(f"  ✅ BG ready: global_n={bg.global_n}, domains={len(bg.dom_mean)}, cells={len(bg.cell_mean)}")

        # ---- Pairwise metrics ----
        group_cols = ["cell_std", "target_std", "domain_type", "bulk_type", "source_db_chem", "source_db_gene"]
        null_cache = {}

        grouper = l1_meta.groupby(group_cols, sort=False)
        for name, gdf in tqdm(grouper, desc=f"Pairwise groups ({track_name})"):
            cell, target, domain, bulk_t, src_c, src_g = name

            # Dedup per side (prevents many-to-many inflation)
            g_chem = gdf.drop_duplicates("uid_chem")
            g_gene = gdf.drop_duplicates("uid_gene")

            idx_chem = torch.as_tensor(g_chem.index.to_numpy(), dtype=torch.long)
            idx_gene = torch.as_tensor(g_gene.index.to_numpy(), dtype=torch.long)

            n_chem = int(idx_chem.numel())
            n_gene = int(idx_gene.numel())
            if n_chem < MIN_UNIQUE_PER_SIDE or n_gene < MIN_UNIQUE_PER_SIDE:
                continue

            # clouds (sanitize)
            cloud_c = l1_chem[idx_chem].cpu().numpy().astype(np.float64)
            cloud_g = l1_gene[idx_gene].cpu().numpy().astype(np.float64)
            cloud_c = _np_sanitize(cloud_c)
            cloud_g = _np_sanitize(cloud_g)

            cent_c = cloud_c.mean(axis=0)
            cent_g = cloud_g.mean(axis=0)

            bg_c, bgc_level, bgc_n = bg.get(src_c, "Chemical", cell)
            bg_g, bgg_level, bgg_n = bg.get(src_g, "Genetic", cell)

            cent_c_sys = cent_c - bg_c
            cent_g_sys = cent_g - bg_g

            cloud_c_sys = cloud_c - bg_c[None, :]
            cloud_g_sys = cloud_g - bg_g[None, :]

            # Standard metrics
            cos_std = robust_cosine(cent_c, cent_g)
            deg_pcc_std, deg_acc_std, jac_abs_std, jac_up_std, jac_dn_std = signal_metrics_union_topk(
                cent_c, cent_g, top_k=topk_deg
            )
            ed_mean_std, ed_sd_std = compute_energy_distance_boot(cloud_c, cloud_g, rng=RNG)

            # Systema metrics
            cos_sys = robust_cosine(cent_c_sys, cent_g_sys)
            deg_pcc_sys, deg_acc_sys, jac_abs_sys, jac_up_sys, jac_dn_sys = signal_metrics_union_topk(
                cent_c_sys, cent_g_sys, top_k=topk_deg
            )
            ed_mean_sys, ed_sd_sys = compute_energy_distance_boot(cloud_c_sys, cloud_g_sys, rng=RNG)

            # Magnitude metrics (B1: denom ~0 => NaN)
            norm_c = float(np.linalg.norm(cent_c))
            norm_g = float(np.linalg.norm(cent_g))
            norm_ratio_std = (norm_c / norm_g) if norm_g > EPS_NORM else np.nan

            norm_c_sys = float(np.linalg.norm(cent_c_sys))
            norm_g_sys = float(np.linalg.norm(cent_g_sys))
            norm_ratio_sys = (norm_c_sys / norm_g_sys) if norm_g_sys > EPS_NORM else np.nan

            # Null baseline (cosine) per stratum
            null_key = (track_name, cell, src_c, src_g, bulk_t)
            if null_key not in null_cache:
                t_common, c_cent, g_cent = _centroids_for_stratum(
                    l1_meta, cell, src_c, src_g, bulk_t, l1_chem, l1_gene
                )
                if t_common.size >= 2:
                    null_std = _null_cosine_dist(c_cent, g_cent, rng=RNG, n_null=N_NULL)

                    c_cent_sys = (c_cent.astype(np.float64) - bg_c[None, :]).astype(np.float32)
                    g_cent_sys = (g_cent.astype(np.float64) - bg_g[None, :]).astype(np.float32)
                    null_sys = _null_cosine_dist(c_cent_sys, g_cent_sys, rng=RNG, n_null=N_NULL)
                else:
                    null_std = np.array([], dtype=np.float32)
                    null_sys = np.array([], dtype=np.float32)

                null_cache[null_key] = {
                    "n_targets_stratum": int(t_common.size),
                    "null_std": null_std,
                    "null_sys": null_sys,
                }

                # ---- Exportable null baseline (one-run) ----
                # This uses mismatched target centroid sampling within the same (cell, src_c, src_g, bulk) stratum.
                # It is lightweight and does NOT change the main Task1 outputs; it only writes extra CSV(s).
                stratum_id = f"{cell}||{src_c}->{src_g}||{bulk_t}"
                s_std = _summarize_null(null_std)
                s_sys = _summarize_null(null_sys)

                all_null_stratum_stats.append(
                    {
                        "Track": track_name,
                        "cell_std": cell,
                        "source_db_chem": src_c,
                        "source_db_gene": src_g,
                        "bulk_type": bulk_t,
                        "stratum_id": stratum_id,
                        "n_targets_stratum": int(t_common.size),

                        "null_n_std": s_std["n"],
                        "null_mu_std": s_std["mu"],
                        "null_sd_std": s_std["sd"],
                        "null_p05_std": s_std["p05"],
                        "null_p50_std": s_std["p50"],
                        "null_p95_std": s_std["p95"],

                        "null_n_sys": s_sys["n"],
                        "null_mu_sys": s_sys["mu"],
                        "null_sd_sys": s_sys["sd"],
                        "null_p05_sys": s_sys["p05"],
                        "null_p50_sys": s_sys["p50"],
                        "null_p95_sys": s_sys["p95"],
                    }
                )

                # Long-form null draws (small: N_NULL per stratum) for quick dashed-line plotting.
                if s_std["n"] > 0:
                    for vv in null_std:
                        all_null_draws.append(
                            {
                                "Track": track_name,
                                "View": "Standard",
                                "cell_std": cell,
                                "source_db_chem": src_c,
                                "source_db_gene": src_g,
                                "bulk_type": bulk_t,
                                "stratum_id": stratum_id,
                                "Value": float(vv),
                            }
                        )
                if s_sys["n"] > 0:
                    for vv in null_sys:
                        all_null_draws.append(
                            {
                                "Track": track_name,
                                "View": "Systema",
                                "cell_std": cell,
                                "source_db_chem": src_c,
                                "source_db_gene": src_g,
                                "bulk_type": bulk_t,
                                "stratum_id": stratum_id,
                                "Value": float(vv),
                            }
                        )

            cached = null_cache[null_key]
            n_targets_stratum = int(cached["n_targets_stratum"])

            z_std, null_mu_std, null_sd_std, null_n_std, null_valid_std = _zscore(cos_std, cached["null_std"])
            z_sys, null_mu_sys, null_sd_sys, null_n_sys, null_valid_sys = _zscore(cos_sys, cached["null_sys"])

            log_ed_std = float(np.log10(ed_mean_std + 1e-8)) if np.isfinite(ed_mean_std) else np.nan
            log_ed_sys = float(np.log10(ed_mean_sys + 1e-8)) if np.isfinite(ed_mean_sys) else np.nan

            all_results.append(
                {
                    "Track": track_name,
                    "cell_std": cell,
                    "target_std": target,
                    "domain_type": domain,
                    "bulk_type": bulk_t,
                    "source_db_chem": src_c,
                    "source_db_gene": src_g,
                    "pair_key": f"{track_name}||{cell}||{target}||{src_c}->{src_g}||{bulk_t}",

                    "n_chem_unique": n_chem,
                    "n_gene_unique": n_gene,
                    "n_targets_stratum": n_targets_stratum,

                    "bg_chem_level": bgc_level,
                    "bg_gene_level": bgg_level,
                    "bg_chem_n": int(bgc_n),
                    "bg_gene_n": int(bgg_n),

                    # cosine
                    "cosine_std": cos_std,
                    "cosine_sys": cos_sys,

                    "cosine_null_mu_std": null_mu_std,
                    "cosine_null_sd_std": null_sd_std,
                    "cosine_null_n_std": int(null_n_std),
                    "cosine_null_valid_std": bool(null_valid_std),

                    "cosine_null_mu_sys": null_mu_sys,
                    "cosine_null_sd_sys": null_sd_sys,
                    "cosine_null_n_sys": int(null_n_sys),
                    "cosine_null_valid_sys": bool(null_valid_sys),

                    "cosine_z_std": z_std,
                    "cosine_z_sys": z_sys,

                    "norm_ratio_std": norm_ratio_std,
                    "norm_ratio_sys": norm_ratio_sys,

                    # signal metrics
                    "deg_pcc_std": deg_pcc_std,
                    "deg_pcc_sys": deg_pcc_sys,
                    "deg_acc_std": deg_acc_std,
                    "deg_acc_sys": deg_acc_sys,
                    "jaccard_abs_std": jac_abs_std,
                    "jaccard_abs_sys": jac_abs_sys,
                    "jaccard_up_std": jac_up_std,
                    "jaccard_up_sys": jac_up_sys,
                    "jaccard_down_std": jac_dn_std,
                    "jaccard_down_sys": jac_dn_sys,

                    # distribution metrics
                    "edist_mean_std": ed_mean_std,
                    "edist_sd_std": ed_sd_std,
                    "edist_mean_sys": ed_mean_sys,
                    "edist_sd_sys": ed_sd_sys,
                    "log10_edist_std": log_ed_std,
                    "log10_edist_sys": log_ed_sys,

                    # gains (sys - std)  (edist gain defined as std - sys so positive = improved/closer)
                    "systema_gain_cosine": cos_sys - cos_std,
                    "systema_gain_cosine_z": (z_sys - z_std) if (np.isfinite(z_sys) and np.isfinite(z_std)) else np.nan,
                    "systema_gain_deg_pcc": deg_pcc_sys - deg_pcc_std,
                    "systema_gain_deg_acc": deg_acc_sys - deg_acc_std,
                    "systema_gain_jaccard_abs": jac_abs_sys - jac_abs_std,
                    "systema_gain_edist_mean": (ed_mean_std - ed_mean_sys) if (np.isfinite(ed_mean_std) and np.isfinite(ed_mean_sys)) else np.nan,
                }
            )

            # cleanup per group
            del cloud_c, cloud_g, cloud_c_sys, cloud_g_sys
            gc.collect()

    # ---- Export ----
    res_df = pd.DataFrame(all_results)

    bg_df = pd.DataFrame(all_bg_stats)
    if not bg_df.empty:
        bg_df.to_csv(os.path.join(OUTPUT_DIR, "Task1_Systema_Background_Stats.csv"), index=False)

    wide_path = os.path.join(OUTPUT_DIR, "Task1_Pairwise_Metrics_Wide.csv")
    res_df.to_csv(wide_path, index=False)

    # Long format for ggplot
    long_rows = []
    if not res_df.empty:
        id_cols = [
            "Track",
            "cell_std",
            "target_std",
            "domain_type",
            "bulk_type",
            "source_db_chem",
            "source_db_gene",
            "pair_key",
            "n_chem_unique",
            "n_gene_unique",
            "n_targets_stratum",
            "bg_chem_level",
            "bg_gene_level",
            "bg_chem_n",
            "bg_gene_n",
            "cosine_null_n_std",
            "cosine_null_valid_std",
            "cosine_null_n_sys",
            "cosine_null_valid_sys",
        ]

        metric_map = {
            "cosine": ("cosine_std", "cosine_sys"),
            "cosine_z": ("cosine_z_std", "cosine_z_sys"),
            "deg_pcc": ("deg_pcc_std", "deg_pcc_sys"),
            "deg_acc": ("deg_acc_std", "deg_acc_sys"),
            "jaccard_abs": ("jaccard_abs_std", "jaccard_abs_sys"),
            "jaccard_up": ("jaccard_up_std", "jaccard_up_sys"),
            "jaccard_down": ("jaccard_down_std", "jaccard_down_sys"),
            "edist_mean": ("edist_mean_std", "edist_mean_sys"),
            "log10_edist": ("log10_edist_std", "log10_edist_sys"),
            "norm_ratio": ("norm_ratio_std", "norm_ratio_sys"),
        }

        for metric, (c_std, c_sys) in metric_map.items():
            if (c_std not in res_df.columns) or (c_sys not in res_df.columns):
                continue
            tmp = res_df[id_cols + [c_std, c_sys]].copy()
            tmp = tmp.rename(columns={c_std: "Standard", c_sys: "Systema"})
            tmp = tmp.melt(id_vars=id_cols, value_vars=["Standard", "Systema"], var_name="View", value_name="Value")
            tmp["Metric"] = metric
            long_rows.append(tmp)

    long_df = pd.concat(long_rows, ignore_index=True) if long_rows else pd.DataFrame()
    long_path = os.path.join(OUTPUT_DIR, "Task1_Pairwise_Metrics_Long.csv")
    long_df.to_csv(long_path, index=False)

    # ---- Extra export: one-run null baseline for dashed-line visualization ----
    # These files are OPTIONAL for downstream analysis and do not change existing outputs.
    if all_null_stratum_stats:
        null_stratum_df = pd.DataFrame(all_null_stratum_stats)
        null_stratum_path = os.path.join(OUTPUT_DIR, "Task1_Pairwise_NullCosine_StratumSummary.csv")
        null_stratum_df.to_csv(null_stratum_path, index=False)

    if all_null_draws:
        null_draws_df = pd.DataFrame(all_null_draws)
        null_draws_path = os.path.join(OUTPUT_DIR, "Task1_Pairwise_NullCosine_Draws_Long.csv")
        null_draws_df.to_csv(null_draws_path, index=False)

        # Global summary per (Track, View) to support a single dashed baseline line.
        def _q(x, qv):
            return float(np.quantile(x, qv)) if len(x) else np.nan

        gsum = (
            null_draws_df
            .groupby(["Track", "View"], as_index=False)["Value"]
            .agg(
                null_n="count",
                null_mu="mean",
                null_sd="std",
                null_p05=lambda x: _q(x.values, 0.05),
                null_p50=lambda x: _q(x.values, 0.50),
                null_p95=lambda x: _q(x.values, 0.95),
            )
        )
        gsum_path = os.path.join(OUTPUT_DIR, "Task1_Pairwise_NullCosine_GlobalSummary.csv")
        gsum.to_csv(gsum_path, index=False)

    print("\n✅ Task 1 Pairwise Complete (hardened)")
    print(f"   BG  : {os.path.join(OUTPUT_DIR, 'Task1_Systema_Background_Stats.csv')}")
    print(f"   Wide: {wide_path}")
    print(f"   Long: {long_path}")

    if not res_df.empty:
        print("\n[Summary] Mean cosine by domain_type (std vs sys):")
        print(res_df.groupby(["Track", "domain_type"])[["cosine_std", "cosine_sys", "systema_gain_cosine"]].mean())
        print("\n[Diagnostics] Null validity rate (std/sys):")
        print(res_df.groupby(["Track"])[["cosine_null_valid_std", "cosine_null_valid_sys"]].mean())
