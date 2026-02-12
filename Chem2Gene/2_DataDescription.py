import os
import torch
import pandas as pd
import numpy as np
import gc
from sklearn.decomposition import PCA

from functions import load_chem2gen_bundle_from_split

# Try UMAP
try:
    import umap
    HAS_UMAP = True
    print(">>> UMAP library detected.")
except ImportError:
    HAS_UMAP = False
    print(">>> UMAP not found. Will use PCA for final 2D projection.")

# ================= Configuration =================
BENCHMARK_DIR = "/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets"
GENE_META_PATH = "/mnt/NAS_21T/ProjectData/OSMOSIS/processed/LINCS_Processed/LINCS_gene_alignment.csv"
OUTPUT_DIR = "/mnt/NAS_21T/ProjectData/Chem2Gen/R_Vis_Ready/Data_Description"
os.makedirs(OUTPUT_DIR, exist_ok=True)

# Total cap for embedding points (across BOTH modalities)
MAX_VIS_SAMPLES_TOTAL = 100000
# Cap for centroid embedding (Pair_ID centroids)
MAX_VIS_CENTROIDS_TOTAL = 60000

RANDOM_SEED = 42


# ================= Helper Functions =================
def safe_col(df: pd.DataFrame, col: str, default="Unknown") -> pd.Series:
    if col not in df.columns:
        return pd.Series([default] * len(df), index=df.index)
    return df[col].fillna(default)

def ensure_pair_id(df: pd.DataFrame, cell_col="cell_std", target_col="target_std", out_col="pair_id") -> pd.DataFrame:
    df = df.copy()
    if out_col not in df.columns:
        df[out_col] = df[cell_col].astype(str) + ":" + df[target_col].astype(str)
    return df

def unique_pairs(df: pd.DataFrame) -> pd.DataFrame:
    df = ensure_pair_id(df)
    return df[["pair_id", "cell_std", "target_std"]].drop_duplicates()

def get_centroids_from_matrix(X: np.ndarray, meta: pd.DataFrame, pair_col="pair_id") -> pd.DataFrame:
    """
    Mean profile per pair_id (centroid). X rows align to meta rows.
    """
    df = pd.DataFrame(X)
    df[pair_col] = meta[pair_col].astype(str).values
    return df.groupby(pair_col, sort=False).mean()

def compute_reduction(tensor_np: np.ndarray, final_dim=2):
    """
    Standard workflow:
    1) If dim > 50: PCA -> 50 PCs
    2) UMAP -> 2D (or PCA -> 2D fallback)
    """
    n_samples, n_features = tensor_np.shape

    if n_features > 50:
        print(f"       ...Pre-reducing {n_features} features to 50 PCs via PCA...")
        pca = PCA(n_components=50, random_state=RANDOM_SEED)
        X_reduced = pca.fit_transform(tensor_np)
        input_for_manifold = X_reduced
        method_prefix = "PCA50+"
    else:
        print(f"       ...Input dimension {n_features} <= 50. Using directly.")
        input_for_manifold = tensor_np
        method_prefix = ""

    if HAS_UMAP:
        reducer = umap.UMAP(
            n_components=final_dim,
            random_state=RANDOM_SEED,
            n_neighbors=30,
            min_dist=0.3,
            low_memory=True
        )
        emb = reducer.fit_transform(input_for_manifold)
        method = f"{method_prefix}UMAP"
    else:
        reducer = PCA(n_components=final_dim, random_state=RANDOM_SEED)
        emb = reducer.fit_transform(input_for_manifold)
        method = f"{method_prefix}PCA"

    return emb, method

def write_csv(df: pd.DataFrame, fname: str):
    out = os.path.join(OUTPUT_DIR, fname)
    df.to_csv(out, index=False)
    print(f"Saved: {out}")

def pick_first_index_per_key(keys: pd.Series) -> np.ndarray:
    """
    Return indices of first occurrence for each unique key (stable).
    """
    # pandas factorize preserves order of appearance
    _, idx = np.unique(keys.to_numpy(), return_index=True)
    return np.sort(idx)

def downsample_indices(idx: np.ndarray, cap: int, seed=RANDOM_SEED) -> np.ndarray:
    if cap is None or len(idx) <= cap:
        return idx
    rng = np.random.default_rng(seed)
    return rng.choice(idx, size=cap, replace=False)

def infer_l1_pairs_unique(data: dict) -> pd.DataFrame:
    """
    Prefer Level1["pairs_unique"] from new bundle;
    else derive from pairs_meta.
    """
    if "Level1" in data and isinstance(data["Level1"], dict):
        if "pairs_unique" in data["Level1"] and isinstance(data["Level1"]["pairs_unique"], pd.DataFrame):
            df = data["Level1"]["pairs_unique"].copy()
            df = df.rename(columns={"Pair_ID": "pair_id"}) if "Pair_ID" in df.columns and "pair_id" not in df.columns else df
            df = ensure_pair_id(df)
            return df[["pair_id", "cell_std", "target_std"]].drop_duplicates()

    # fallback
    l1_pairs = data["Level1"]["pairs_meta"]
    if not isinstance(l1_pairs, pd.DataFrame):
        raise ValueError("Cannot infer Level1 unique pairs: pairs_meta is not a DataFrame.")
    df = l1_pairs[["cell_std", "target_std"]].drop_duplicates().copy()
    df = ensure_pair_id(df)
    return df[["pair_id", "cell_std", "target_std"]].drop_duplicates()

def pair_stats_from_pools(chem_pool: pd.DataFrame, gene_pool: pd.DataFrame) -> pd.DataFrame:
    """
    Pair_ID-level summary for pools.
    """
    c = ensure_pair_id(chem_pool)
    g = ensure_pair_id(gene_pool)

    c_stats = (c.groupby(["pair_id", "cell_std", "target_std"])
                 .agg(n_chem=("pair_id", "size"),
                      chem_sources=("source_db", lambda s: "|".join(sorted(set(map(str, s))))),
                      chem_tissues=("tissue", lambda s: "|".join(sorted(set(map(str, s.dropna())))) if "tissue" in c.columns else "Unknown"))
                 .reset_index())
    g_stats = (g.groupby(["pair_id", "cell_std", "target_std"])
                 .agg(n_gene=("pair_id", "size"),
                      gene_sources=("source_db", lambda s: "|".join(sorted(set(map(str, s))))),
                      gene_tissues=("tissue", lambda s: "|".join(sorted(set(map(str, s.dropna())))) if "tissue" in g.columns else "Unknown"))
                 .reset_index())

    out = pd.merge(c_stats, g_stats, on=["pair_id", "cell_std", "target_std"], how="outer")
    for col in ["n_chem", "n_gene"]:
        if col in out.columns:
            out[col] = out[col].fillna(0).astype(int)
    out["has_chem"] = out["n_chem"] > 0
    out["has_gene"] = out["n_gene"] > 0
    return out

def extract_time_from_pairs(df_pairs, cols, mod):
    found = next((c for c in cols if c in df_pairs.columns), None)
    if found is None:
        return pd.DataFrame()
    vc = df_pairs[found].value_counts(dropna=False).reset_index()
    vc.columns = ["Time", "Count"]
    vc["Modality"] = mod
    return vc

def build_meta_from_l1pairs(sub_df: pd.DataFrame, modality: str):
    sub_df = sub_df.copy()
    sub_df = ensure_pair_id(sub_df, out_col="pair_id")
    if modality == "Chemical":
        tissue = safe_col(sub_df, "tissue_chem", "Unknown")
        src = safe_col(sub_df, "source_db_chem", "Unknown")
    else:
        tissue = safe_col(sub_df, "tissue_gene", "Unknown")
        src = safe_col(sub_df, "source_db_gene", "Unknown")
    return pd.DataFrame({
        "cell_std": sub_df["cell_std"].astype(str).values,
        "target_std": sub_df["target_std"].astype(str).values,
        "pair_id": sub_df["pair_id"].astype(str).values,
        "tissue": tissue.astype(str).values,
        "SourceDB": src.astype(str).values,
        "Modality": modality,
    })

def process_embedding_instance_dedup(tkey: str, tensor_name: str):
    print(f"    -> Embedding (instance-dedup): {tensor_name} [{tkey}]")

    t_c = chem_tensors[tkey][idx_chem_use].float().numpy()
    t_g = gene_tensors[tkey][idx_gene_use].float().numpy()

    X = np.concatenate([t_c, t_g], axis=0).astype(np.float32)
    del t_c, t_g
    gc.collect()

    meta_c = build_meta_from_l1pairs(l1_pairs.iloc[idx_chem_use], "Chemical")
    meta_g = build_meta_from_l1pairs(l1_pairs.iloc[idx_gene_use], "Genetic")
    meta = pd.concat([meta_c, meta_g], ignore_index=True)

    emb, method = compute_reduction(X)
    df_emb = pd.DataFrame(emb, columns=["Dim1", "Dim2"])
    df_out = pd.concat([df_emb, meta.reset_index(drop=True)], axis=1)
    df_out["Type"] = tensor_name
    df_out["Method"] = method
    df_out["Granularity"] = "InstanceDedup"
    return df_out

def process_embedding_pair_centroids(tkey: str, tensor_name: str):
    print(f"    -> Embedding (pair-centroids): {tensor_name} [{tkey}]")

    # Use dedup instances for centroid calculation
    # (dedup by uid removes L1 cartesian bias; then centroid by pair_id)
    t_c = chem_tensors[tkey][idx_chem_unique].float().numpy().astype(np.float32)
    t_g = gene_tensors[tkey][idx_gene_unique].float().numpy().astype(np.float32)

    meta_c = build_meta_from_l1pairs(l1_pairs.iloc[idx_chem_unique], "Chemical")
    meta_g = build_meta_from_l1pairs(l1_pairs.iloc[idx_gene_unique], "Genetic")

    # centroid in each modality by pair_id
    c_cent = get_centroids_from_matrix(t_c, meta_c, pair_col="pair_id")
    g_cent = get_centroids_from_matrix(t_g, meta_g, pair_col="pair_id")

    # keep only pairs that exist in BOTH modalities (true L1 coverage)
    common_pairs = sorted(set(c_cent.index.astype(str)).intersection(set(g_cent.index.astype(str))))
    if len(common_pairs) == 0:
        print("       ! No common pair_id between modalities for centroid embedding.")
        return pd.DataFrame()

    c_cent = c_cent.loc[common_pairs]
    g_cent = g_cent.loc[common_pairs]

    # cap total centroids
    cap_pairs = max(int(MAX_VIS_CENTROIDS_TOTAL // 2), 1)
    if len(common_pairs) > cap_pairs:
        rng = np.random.default_rng(RANDOM_SEED)
        keep_pairs = rng.choice(common_pairs, size=cap_pairs, replace=False)
        keep_pairs = sorted(list(keep_pairs))
        c_cent = c_cent.loc[keep_pairs]
        g_cent = g_cent.loc[keep_pairs]
        common_pairs = keep_pairs

    X = np.concatenate([c_cent.to_numpy(), g_cent.to_numpy()], axis=0).astype(np.float32)

    # metadata for centroids: take (cell,target) from l1_pairs_unique
    l1u = l1_pairs_unique.set_index("pair_id").loc[common_pairs]
    meta_base = pd.DataFrame({
        "pair_id": common_pairs,
        "cell_std": l1u["cell_std"].astype(str).values,
        "target_std": l1u["target_std"].astype(str).values
    })

    meta_c2 = meta_base.copy()
    meta_c2["Modality"] = "Chemical"
    meta_g2 = meta_base.copy()
    meta_g2["Modality"] = "Genetic"
    meta = pd.concat([meta_c2, meta_g2], ignore_index=True)
    meta["tissue"] = "Unknown"
    meta["SourceDB"] = "MixedOrUnknown"

    emb, method = compute_reduction(X)
    df_emb = pd.DataFrame(emb, columns=["Dim1", "Dim2"])
    df_out = pd.concat([df_emb, meta.reset_index(drop=True)], axis=1)
    df_out["Type"] = tensor_name
    df_out["Method"] = method
    df_out["Granularity"] = "PairCentroid"
    return df_out


# ================= Main =================
if __name__ == "__main__":
    print(">>> [Init] Loading Data...")
    data = load_chem2gen_bundle_from_split(
    BENCHMARK_DIR,
    load_unified_meta=True,
    tensor_cache_max_items=2,
    strict=True,
)
    # ----- Align to new schema (and still work with old) -----
    l1_pairs = data["Level1"]["pairs_meta"]  # instance-level many-to-many pairs
    chem_tensors = data["Level1"]["chem_tensors"]
    gene_tensors = data["Level1"]["gene_tensors"]

    # Level3 pools are used for data description coverage
    l3_chem = data["Level3"]["chem_pool"].copy()
    l3_gene = data["Level3"]["gene_pool"].copy()

    # Ensure key columns
    for df in [l3_chem, l3_gene]:
        if "cell_std" not in df.columns or "target_std" not in df.columns:
            raise ValueError("Level3 pools must contain cell_std and target_std.")

    # Standardize columns expected by downstream R
    l3_chem["Modality"] = "Chemical"
    l3_gene["Modality"] = "Genetic"

    # unify source_db/tissue
    l3_chem["source_db"] = safe_col(l3_chem, "source_db", "Unknown")
    l3_gene["source_db"] = safe_col(l3_gene, "source_db", "Unknown")
    l3_chem["tissue"] = safe_col(l3_chem, "tissue", "Unknown")
    l3_gene["tissue"] = safe_col(l3_gene, "tissue", "Unknown")

    l3_chem = ensure_pair_id(l3_chem)
    l3_gene = ensure_pair_id(l3_gene)

    # L1 unique pairs (coverage definition)
    l1_pairs_unique = infer_l1_pairs_unique(data)
    l1_pairs_unique = ensure_pair_id(l1_pairs_unique)
    l1_pair_set = set(l1_pairs_unique["pair_id"].astype(str).tolist())
    l1_targets = set(l1_pairs_unique["target_std"].astype(str).tolist())

    # Build Level3 global table for description
    cols_global = ["Modality", "source_db", "tissue", "cell_std", "target_std", "pair_id"]
    df_global = pd.concat([l3_chem[cols_global], l3_gene[cols_global]], ignore_index=True)

    print(f">>> Level3 records: chem={len(l3_chem)} | gene={len(l3_gene)} | total={len(df_global)}")
    print(f">>> Level1 unique pairs: {len(l1_pairs_unique)} | unique targets: {len(l1_targets)}")

    # -------------------------------------------------------
    # 1) Global Hierarchy (Level3 pools composition)
    # -------------------------------------------------------
    print(">>> [1/9] Generating Hierarchy Data (Level3)...")

    # Record-count version (kept as the default file name to match your R script expectation)
    df_polar_record = (
        df_global.groupby(["Modality", "tissue", "cell_std"])
        .size().reset_index(name="Count")
    )
    write_csv(df_polar_record, "Hierarchy_Polar.csv")  # BACKWARD-COMPAT default (record counts)
    write_csv(df_polar_record, "Hierarchy_Polar_RecordCount.csv")

    # Pair-count version (coverage, recommended for summary figures)
    df_polar_pair = (
        df_global.drop_duplicates(["Modality", "source_db", "pair_id"])
        .groupby(["Modality", "tissue", "cell_std"])
        .size().reset_index(name="Count")
    )
    write_csv(df_polar_pair, "Hierarchy_Polar_PairCount.csv")

    # -------------------------------------------------------
    # 2) Fig1B1 NEW: SourceDB panels + stacked bars (Pair_ID coverage)
    # -------------------------------------------------------
    print(">>> [2/9] Generating Fig1B1 stacked-bar inputs (SourceDB panels)...")

    # Bar1: modality ratio per SourceDB (Pair_ID coverage)
    mod_bar = (
        df_global.drop_duplicates(["source_db", "Modality", "pair_id"])
        .groupby(["source_db", "Modality"])
        .size().reset_index(name="Count")
    )
    mod_bar["Bar"] = "Modality"
    mod_bar = mod_bar.rename(columns={"source_db": "SourceDB", "Modality": "Segment"})

    # Bar2: Level1 vs Level3 per SourceDB (Pair_ID coverage, modality-agnostic)
    rows = []
    for src, sub in df_global.groupby("source_db"):
        pairs_src = set(sub["pair_id"].astype(str).unique().tolist())
        n_l3 = len(pairs_src)
        n_l1 = len(pairs_src.intersection(l1_pair_set))
        rows.append({"SourceDB": src, "Bar": "Level", "Segment": "Level1", "Count": n_l1})
        rows.append({"SourceDB": src, "Bar": "Level", "Segment": "Level3_only", "Count": max(n_l3 - n_l1, 0)})

    lvl_bar = pd.DataFrame(rows)

    fig1b1 = pd.concat([mod_bar[["SourceDB", "Bar", "Segment", "Count"]],
                        lvl_bar[["SourceDB", "Bar", "Segment", "Count"]]], ignore_index=True)
    fig1b1["CountUnit"] = "PairID"
    write_csv(fig1b1, "Fig1B1_SourcePanels_StackedBars_PairCount.csv")

    # Optional: record-count version (may be useful for sanity checks)
    mod_bar_rec = (
        df_global.groupby(["source_db", "Modality"])
        .size().reset_index(name="Count")
        .rename(columns={"source_db": "SourceDB", "Modality": "Segment"})
    )
    mod_bar_rec["Bar"] = "Modality"
    rows = []
    for src, sub in df_global.groupby("source_db"):
        # record-level membership for L1: record belongs to a pair_id in L1 coverage
        n_l3 = len(sub)
        n_l1 = int(np.sum(sub["pair_id"].astype(str).isin(l1_pair_set)))
        rows.append({"SourceDB": src, "Bar": "Level", "Segment": "Level1", "Count": n_l1})
        rows.append({"SourceDB": src, "Bar": "Level", "Segment": "Level3_only", "Count": max(n_l3 - n_l1, 0)})
    lvl_bar_rec = pd.DataFrame(rows)
    fig1b1_rec = pd.concat([mod_bar_rec[["SourceDB", "Bar", "Segment", "Count"]],
                            lvl_bar_rec[["SourceDB", "Bar", "Segment", "Count"]]], ignore_index=True)
    fig1b1_rec["CountUnit"] = "Record"
    write_csv(fig1b1_rec, "Fig1B1_SourcePanels_StackedBars_RecordCount.csv")

    # -------------------------------------------------------
    # 3) Butterfly Data (L1 targets within Level3 pools)
    # -------------------------------------------------------
    print(">>> [3/9] Generating Butterfly Data...")

    l3_chem_matched = l3_chem[l3_chem["target_std"].astype(str).isin(l1_targets)].copy()
    l3_gene_matched = l3_gene[l3_gene["target_std"].astype(str).isin(l1_targets)].copy()

    # Record counts
    chem_counts_rec = l3_chem_matched["target_std"].value_counts()
    gen_counts_rec = l3_gene_matched["target_std"].value_counts()
    df_butter_rec = (
        pd.DataFrame({"Chemical_Count": chem_counts_rec, "Genetic_Count": gen_counts_rec})
        .fillna(0)
    )
    df_butter_rec["Total"] = df_butter_rec["Chemical_Count"] + df_butter_rec["Genetic_Count"]
    df_butter_rec.index.name = "Target"
    df_butter_rec = df_butter_rec.reset_index()
    df_butter_rec.sort_values("Total", ascending=False).head(30).to_csv(
        os.path.join(OUTPUT_DIR, "Butterfly_L1_Targets.csv"), index=False
    )
    write_csv(df_butter_rec.sort_values("Total", ascending=False).head(30), "Butterfly_L1_Targets_RecordCount.csv")

    # Pair counts (coverage)
    chem_counts_pair = (l3_chem_matched.drop_duplicates(["target_std", "pair_id"])
                        .groupby("target_std").size())
    gen_counts_pair = (l3_gene_matched.drop_duplicates(["target_std", "pair_id"])
                       .groupby("target_std").size())
    df_butter_pair = (pd.DataFrame({"Chemical_Count": chem_counts_pair, "Genetic_Count": gen_counts_pair})
                      .fillna(0))
    df_butter_pair["Total"] = df_butter_pair["Chemical_Count"] + df_butter_pair["Genetic_Count"]
    df_butter_pair.index.name = "Target"
    df_butter_pair = df_butter_pair.reset_index()
    write_csv(df_butter_pair.sort_values("Total", ascending=False).head(30), "Butterfly_L1_Targets_PairCount.csv")

    # All targets (record + pair)
    chem_all_rec = l3_chem["target_std"].value_counts()
    gen_all_rec = l3_gene["target_std"].value_counts()
    df_all_targets_rec = pd.DataFrame({"Total": chem_all_rec.add(gen_all_rec, fill_value=0)})
    df_all_targets_rec.index.name = "Target"
    write_csv(df_all_targets_rec.reset_index(), "All_Target_Counts_RecordCount.csv")

    chem_all_pair = (l3_chem.drop_duplicates(["target_std", "pair_id"]).groupby("target_std").size())
    gen_all_pair = (l3_gene.drop_duplicates(["target_std", "pair_id"]).groupby("target_std").size())
    df_all_targets_pair = pd.DataFrame({"Total": chem_all_pair.add(gen_all_pair, fill_value=0)})
    df_all_targets_pair.index.name = "Target"
    write_csv(df_all_targets_pair.reset_index(), "All_Target_Counts_PairCount.csv")

    # -------------------------------------------------------
    # 4) Treatment Meta (from Level1 instance pairs)
    # -------------------------------------------------------
    print(">>> [4/9] Extracting Treatment Meta (Level1 instance pairs)...")

    time_cols_chem = ["time_val_chem", "pert_time_chem", "time_point_chem"]
    time_cols_gene = ["time_val_gene", "pert_time_gene", "time_point_gene"]

    t_chem = extract_time_from_pairs(l1_pairs, time_cols_chem, "Chemical")
    t_gene = extract_time_from_pairs(l1_pairs, time_cols_gene, "Genetic")

    if t_gene.empty:
        t_gene = pd.DataFrame({"Time": ["IGNORED"], "Count": [len(l1_pairs)], "Modality": ["Genetic"]})

    write_csv(pd.concat([t_chem, t_gene], ignore_index=True), "Meta_Time_Dist.csv")

    dose_cols = ["dose_val_chem", "pert_dose_chem", "dose_chem", "concentration_chem"]
    found_dose = next((c for c in dose_cols if c in l1_pairs.columns), None)
    if found_dose is not None:
        write_csv(l1_pairs[[found_dose]].rename(columns={found_dose: "Dose"}), "Meta_Dose_Chem.csv")
    else:
        write_csv(pd.DataFrame({"Dose": [np.nan] * len(l1_pairs)}), "Meta_Dose_Chem.csv")

    # -------------------------------------------------------
    # 5) Embeddings (Instance-dedup + Pair-centroids) from Level1 tensors
    # -------------------------------------------------------
    print(f">>> [5/9] Computing Reductions (cap total={MAX_VIS_SAMPLES_TOTAL})...")

    keys = set(chem_tensors.keys()).intersection(set(gene_tensors.keys()))
    key_base = "x" if "x" in keys else next(iter(keys))
    key_d_gene = "y_gene" if "y_gene" in keys else None
    key_d_path = "y_path" if "y_path" in keys else None

    # --- Dedup by uid to avoid L1 many-to-many pairing bias ---
    # Chemical side dedup indices in l1_pairs table
    if "uid_chem" in l1_pairs.columns:
        idx_chem_unique = pick_first_index_per_key(l1_pairs["uid_chem"].astype(str))
    else:
        idx_chem_unique = np.arange(len(l1_pairs))
    # Genetic side
    if "uid_gene" in l1_pairs.columns:
        idx_gene_unique = pick_first_index_per_key(l1_pairs["uid_gene"].astype(str))
    else:
        idx_gene_unique = np.arange(len(l1_pairs))

    # cap per modality so total <= MAX_VIS_SAMPLES_TOTAL
    cap_per_mod = max(int(MAX_VIS_SAMPLES_TOTAL // 2), 1)
    idx_chem_use = downsample_indices(idx_chem_unique, cap_per_mod, seed=RANDOM_SEED)
    idx_gene_use = downsample_indices(idx_gene_unique, cap_per_mod, seed=RANDOM_SEED + 1)

    outs = []
    outs.append(process_embedding_instance_dedup(key_base, "Baseline"))
    if key_d_gene is not None:
        outs.append(process_embedding_instance_dedup(key_d_gene, "Delta_Gene"))
    if key_d_path is not None:
        outs.append(process_embedding_instance_dedup(key_d_path, "Delta_Pathway"))

    emb_inst = pd.concat([x for x in outs if isinstance(x, pd.DataFrame) and len(x) > 0], ignore_index=True)
    write_csv(emb_inst, "Embedding_Combined_InstanceDedup.csv")

    outs2 = []
    outs2.append(process_embedding_pair_centroids(key_base, "Baseline"))
    if key_d_gene is not None:
        outs2.append(process_embedding_pair_centroids(key_d_gene, "Delta_Gene"))
    if key_d_path is not None:
        outs2.append(process_embedding_pair_centroids(key_d_path, "Delta_Pathway"))

    emb_cent = pd.concat([x for x in outs2 if isinstance(x, pd.DataFrame) and len(x) > 0], ignore_index=True)
    if len(emb_cent) > 0:
        write_csv(emb_cent, "Embedding_Combined_PairCentroids.csv")

    # -------------------------------------------------------
    # 6) Sankey classification (Level3 pool, 3-class decomposition w.r.t L1 targets/pairs)
    # -------------------------------------------------------
    print(">>> [6/9] Generating Sankey Classification (Level3, Pair_ID)...")

    # Define classes for each modality in Level3:
    # L1: pair_id in L1 coverage
    # L2-like: target in L1 targets but pair_id not in L1 (target-match only)
    # L3-only: target not in L1 targets (cell context only / cross-target)
    rows = []
    for mod_name, df_mod in [("Chemical", l3_chem), ("Genetic", l3_gene)]:
        dfm = ensure_pair_id(df_mod)
        # pair-count
        dfm_u = dfm.drop_duplicates(["pair_id"]).copy()
        cls = np.where(dfm_u["pair_id"].astype(str).isin(l1_pair_set), "L1: Exact (cell+target)",
                       np.where(dfm_u["target_std"].astype(str).isin(l1_targets), "L2: Target Match Only",
                                "L3: Cell Context Only"))
        dfm_u["Class"] = cls
        out = dfm_u.groupby("Class").size().reset_index(name="Count")
        for _, r in out.iterrows():
            rows.append({
                "Modality": mod_name,
                "Source": "Level3 pools (cells from L1 expansion)",
                "Class": r["Class"],
                "Count": int(r["Count"])
            })
    sankey_df = pd.DataFrame(rows)
    write_csv(sankey_df, "Sankey_Classification.csv")
    write_csv(sankey_df, "Sankey_Classification_PairCount.csv")

    # Optional: source-specific sankey (for paneling in R if you want)
    rows = []
    for mod_name, df_mod in [("Chemical", l3_chem), ("Genetic", l3_gene)]:
        dfm = ensure_pair_id(df_mod)
        for src, sub in dfm.groupby("source_db"):
            sub_u = sub.drop_duplicates(["pair_id"]).copy()
            cls = np.where(sub_u["pair_id"].astype(str).isin(l1_pair_set), "L1: Exact (cell+target)",
                           np.where(sub_u["target_std"].astype(str).isin(l1_targets), "L2: Target Match Only",
                                    "L3: Cell Context Only"))
            sub_u["Class"] = cls
            out = sub_u.groupby("Class").size().reset_index(name="Count")
            for _, r in out.iterrows():
                rows.append({
                    "Modality": mod_name,
                    "SourceDB": src,
                    "Class": r["Class"],
                    "Count": int(r["Count"])
                })
    sankey_src = pd.DataFrame(rows)
    write_csv(sankey_src, "Sankey_Classification_BySourceDB_PairCount.csv")

    # -------------------------------------------------------
    # 7) UpSet (Level3; record-based presence matrix)
    # -------------------------------------------------------
    print(">>> [7/9] Generating UpSet (Level3)...")
    top_cells = df_global["cell_std"].value_counts().head(20).index.tolist()
    data_for_upset = df_global[df_global["cell_std"].isin(top_cells)]
    upset_df = pd.crosstab(data_for_upset["target_std"], data_for_upset["cell_std"])
    upset_df = (upset_df > 0).astype(int)
    upset_df.to_csv(os.path.join(OUTPUT_DIR, "UpSet_Top20Cells.csv"), index_label="Target")
    print(f"Saved: {os.path.join(OUTPUT_DIR, 'UpSet_Top20Cells.csv')}")

    # -------------------------------------------------------
    # 8) Gene Heatmap (TopVar50; de-biased by uid-dedup before variance)
    # -------------------------------------------------------
    print(">>> [8/9] Extracting Gene Heatmap (TopVar50, uid-dedup)...")
    if key_d_gene is None:
        print("    ! No y_gene found.")
    else:
        # Dedup indices again
        idx_c = idx_chem_unique
        idx_g = idx_gene_unique

        t_chem_all = chem_tensors[key_d_gene][idx_c].float().numpy().astype(np.float32)
        t_gen_all  = gene_tensors[key_d_gene][idx_g].float().numpy().astype(np.float32)

        combined_var = np.var(t_chem_all, axis=0) + np.var(t_gen_all, axis=0)
        top_50_idx = np.argsort(combined_var)[::-1][:50]

        try:
            gene_meta_df = pd.read_csv(GENE_META_PATH)
            if "pr_gene_symbol" in gene_meta_df.columns and len(gene_meta_df) >= t_chem_all.shape[1]:
                gene_col_names = gene_meta_df["pr_gene_symbol"].astype(str).values[top_50_idx]
            else:
                gene_col_names = [f"G_{i}" for i in top_50_idx]
        except Exception:
            gene_col_names = [f"G_{i}" for i in top_50_idx]

        meta_c = l1_pairs.iloc[idx_c][["cell_std", "target_std"]].copy()
        meta_g = l1_pairs.iloc[idx_g][["cell_std", "target_std"]].copy()
        meta_c = ensure_pair_id(meta_c, out_col="pair_id")
        meta_g = ensure_pair_id(meta_g, out_col="pair_id")

        df_chem_g = pd.DataFrame(t_chem_all[:, top_50_idx], columns=gene_col_names)
        df_chem_g["Pair_ID"] = meta_c["pair_id"].values
        df_chem_g = df_chem_g.groupby("Pair_ID").mean()

        df_gen_g = pd.DataFrame(t_gen_all[:, top_50_idx], columns=gene_col_names)
        df_gen_g["Pair_ID"] = meta_g["pair_id"].values
        df_gen_g = df_gen_g.groupby("Pair_ID").mean()

        common_idx = df_chem_g.index.intersection(df_gen_g.index)
        df_chem_g.loc[common_idx].to_csv(os.path.join(OUTPUT_DIR, "Heatmap_Gene_Chem_TopVar50.csv"))
        df_gen_g.loc[common_idx].to_csv(os.path.join(OUTPUT_DIR, "Heatmap_Gene_Gen_TopVar50.csv"))

        meta_out = pd.DataFrame({"Pair_ID": common_idx})
        meta_out[["cell_std", "target_std"]] = meta_out["Pair_ID"].str.split(":", expand=True)
        meta_out.to_csv(os.path.join(OUTPUT_DIR, "Heatmap_Meta.csv"), index=False)

        print(f"Saved: {os.path.join(OUTPUT_DIR, 'Heatmap_Gene_Chem_TopVar50.csv')}")
        print(f"Saved: {os.path.join(OUTPUT_DIR, 'Heatmap_Gene_Gen_TopVar50.csv')}")
        print(f"Saved: {os.path.join(OUTPUT_DIR, 'Heatmap_Meta.csv')}")

        del t_chem_all, t_gen_all
        gc.collect()

    # -------------------------------------------------------
    # 9) Pathway Heatmap (centroids by Pair_ID; de-biased by uid-dedup)
    # -------------------------------------------------------
    print(">>> [9/9] Extracting Pathway Heatmap (uid-dedup)...")
    if key_d_path is not None:
        idx_c = idx_chem_unique
        idx_g = idx_gene_unique

        t_p_c = chem_tensors[key_d_path][idx_c].float().numpy().astype(np.float32)
        t_p_g = gene_tensors[key_d_path][idx_g].float().numpy().astype(np.float32)
        p_cols = [f"HM_{i}" for i in range(t_p_c.shape[1])]

        meta_c = l1_pairs.iloc[idx_c][["cell_std", "target_std"]].copy()
        meta_g = l1_pairs.iloc[idx_g][["cell_std", "target_std"]].copy()
        meta_c = ensure_pair_id(meta_c, out_col="pair_id")
        meta_g = ensure_pair_id(meta_g, out_col="pair_id")

        df_pc = pd.DataFrame(t_p_c, columns=p_cols)
        df_pc["Pair_ID"] = meta_c["pair_id"].values
        df_pc = df_pc.groupby("Pair_ID").mean()

        df_pg = pd.DataFrame(t_p_g, columns=p_cols)
        df_pg["Pair_ID"] = meta_g["pair_id"].values
        df_pg = df_pg.groupby("Pair_ID").mean()

        common_idx2 = df_pc.index.intersection(df_pg.index)
        df_pc.loc[common_idx2].to_csv(os.path.join(OUTPUT_DIR, "Heatmap_Pathway_Chem.csv"))
        df_pg.loc[common_idx2].to_csv(os.path.join(OUTPUT_DIR, "Heatmap_Pathway_Gen.csv"))

        print(f"Saved: {os.path.join(OUTPUT_DIR, 'Heatmap_Pathway_Chem.csv')}")
        print(f"Saved: {os.path.join(OUTPUT_DIR, 'Heatmap_Pathway_Gen.csv')}")

        del t_p_c, t_p_g
        gc.collect()
    else:
        print("    ! No y_path found.")

    print(f"✅ Data processing complete. All files saved to: {OUTPUT_DIR}")