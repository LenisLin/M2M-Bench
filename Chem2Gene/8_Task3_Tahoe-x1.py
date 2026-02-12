import os
import argparse
import subprocess
import pickle
import numpy as np
import pandas as pd
import torch
import scanpy as sc
from anndata import AnnData
from scipy import sparse
import yaml
# python 8_Task3_Tahoe-x1.py --model_size 3b --batch_size 8

def load_chem2gen_eval(eval_dir: str):
    var_names = pd.read_csv(os.path.join(eval_dir, "shared_var_names.csv"))["gene_symbol"].tolist()

    g_counts = torch.load(os.path.join(eval_dir, "CRISPR_counts.pt"), map_location="cpu")
    g_meta = pd.read_csv(os.path.join(eval_dir, "CRISPR_meta.csv"), index_col=0)

    c_counts = torch.load(os.path.join(eval_dir, "Drug_counts.pt"), map_location="cpu")
    c_meta = pd.read_csv(os.path.join(eval_dir, "Drug_meta.csv"), index_col=0)

    assert g_counts.shape[0] == g_meta.shape[0]
    assert c_counts.shape[0] == c_meta.shape[0]
    assert g_counts.shape[1] == len(var_names)
    assert c_counts.shape[1] == len(var_names)

    return var_names, g_counts, g_meta, c_counts, c_meta


def load_symbol_to_ensembl(pkl_path: str) -> dict:
    with open(pkl_path, "rb") as f:
        d = pickle.load(f)
    if not isinstance(d, dict):
        raise TypeError(f"Expected dict in {pkl_path}, got {type(d)}")
    return d


def build_h5ad_for_tahoe(
    counts: torch.Tensor,
    obs_df: pd.DataFrame,
    gene_symbols: list[str],
    symbol2ens: dict,
    out_h5ad: str,
    gene_id_key: str = "ensembl_id",
    cell_type_key: str = "cell_type",
    default_cell_type: str = "unknown",
    make_sparse: bool = True,
):
    os.makedirs(os.path.dirname(out_h5ad), exist_ok=True)

    ens = [symbol2ens.get(gs, None) for gs in gene_symbols]
    keep_mask = np.array([e is not None for e in ens], dtype=bool)

    kept_gene_symbols = [g for g, k in zip(gene_symbols, keep_mask) if k]
    kept_ens = [e for e in ens if e is not None]

    if len(kept_gene_symbols) < 500:
        raise RuntimeError(
            f"Too few genes mapped to Ensembl: {len(kept_gene_symbols)}. "
            f"Check your mapping dict and gene symbol convention."
        )

    if counts.dtype != torch.float32:
        counts = counts.float()

    keep_idx = torch.from_numpy(np.where(keep_mask)[0]).long()
    counts = counts[:, keep_idx]

    X = counts.numpy()
    if make_sparse:
        X = sparse.csr_matrix(X, dtype=np.float32)
    else:
        X = np.asarray(X, dtype=np.float32)

    var = pd.DataFrame(index=pd.Index(kept_gene_symbols, name="gene_symbol"))
    var[gene_id_key] = kept_ens

    obs = obs_df.copy()
    if cell_type_key not in obs.columns:
        obs[cell_type_key] = default_cell_type

    adata = AnnData(X=X, obs=obs, var=var)
    adata.obs_names = obs.index.astype(str)
    adata.var_names_make_unique()
    adata.write_h5ad(out_h5ad)

    return len(gene_symbols), len(kept_gene_symbols)

def write_predict_yaml(
    out_yaml: str,
    model_name: str,
    hf_repo_id: str,
    hf_model_size: str,
    adata_input: str,
    adata_output: str,
    cell_type_key: str,
    gene_id_key: str,
    seq_len_dataset: int,
    return_gene_embeddings: bool,
    model_dir: str,
):
    cfg = {
        "model_name": model_name, 
        "paths": {
            "model_name": model_name,
            "adata_input": adata_input,
            "adata_output": adata_output,
        },
        "data": {"cell_type_key": cell_type_key, "gene_id_key": gene_id_key},
        "predict": {
            "seq_len_dataset": int(seq_len_dataset),
            "return_gene_embeddings": bool(return_gene_embeddings),
        },
    }

    if model_dir is not None:
        cfg["paths"]["model_dir"] = model_dir
    else:
        cfg["paths"]["hf_repo_id"] = hf_repo_id
        cfg["paths"]["hf_model_size"] = hf_model_size

    os.makedirs(os.path.dirname(out_yaml), exist_ok=True)
    with open(out_yaml, "w") as f:
        yaml.safe_dump(cfg, f, sort_keys=False)
    return out_yaml


def run_predict_embeddings(cfg_yaml: str, batch_size: int | None = None):
    cmd = ["python", "benchmark/tahoe-x1/scripts/inference/predict_embeddings.py", cfg_yaml]
    if batch_size is not None:
        # Tahoe-x1 README: 可以用 --batch_size=128 覆盖。:contentReference[oaicite:2]{index=2}
        cmd += [f"--predict.batch_size={batch_size}"]

    print("\n[Tahoe-x1] Running:\n  " + " ".join(cmd) + "\n")
    subprocess.run(cmd, check=True)


def export_embedding(h5ad_path: str, obsm_key: str, out_npy: str):
    adata = sc.read_h5ad(h5ad_path)
    if obsm_key not in adata.obsm:
        raise KeyError(f"Missing obsm['{obsm_key}'] in {h5ad_path}. Available: {list(adata.obsm.keys())}")
    emb = np.asarray(adata.obsm[obsm_key])
    os.makedirs(os.path.dirname(out_npy), exist_ok=True)
    np.save(out_npy, emb)
    print(f"[OK] Saved: {out_npy} | key={obsm_key} | shape={emb.shape}")


def main():
    ap = argparse.ArgumentParser("Chem2Gen → Tahoe-x1 embeddings (run INSIDE container)")

    ap.add_argument("--eval_dir", type=str,
                    default="/data/Benchmark_Datasets/Evaluation_Set_K562")
    ap.add_argument("--work_dir", type=str,
                    default="/data/Model_Evaluation_Results/TahoeX1_workdir_K562")
    ap.add_argument("--out_dir", type=str,
                    default="/data/Model_Evaluation_Results/TahoeX1")

    ap.add_argument("--gene_map_pkl", type=str,
                    default="/workspace/benchmark/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl")

    ap.add_argument("--hf_repo_id", type=str, default="tahoebio/Tahoe-x1")
    ap.add_argument("--model_size", type=str, default="70m", choices=["70m", "1b", "3b"])
    ap.add_argument("--seq_len_dataset", type=int, default=2048)

    ap.add_argument("--cell_type_key", type=str, default="cell_type")
    ap.add_argument("--gene_id_key", type=str, default="ensembl_id")
    ap.add_argument("--make_sparse", action=argparse.BooleanOptionalAction, default=True)

    ap.add_argument("--batch_size", type=int, default=1,
                    help="Optional; passed to predict_embeddings.py as --batch_size=...")

    ap.add_argument("--return_gene_embeddings", action=argparse.BooleanOptionalAction, default=False)

    args = ap.parse_args()

    os.makedirs(args.work_dir, exist_ok=True)
    os.makedirs(args.out_dir, exist_ok=True)

    var_names, g_counts, g_meta, c_counts, c_meta = load_chem2gen_eval(args.eval_dir)
    symbol2ens = load_symbol_to_ensembl(args.gene_map_pkl)

    inputs_dir = os.path.join(args.work_dir, "inputs")
    cfg_dir = os.path.join(args.work_dir, "configs")
    out_h5ad_dir = os.path.join(args.work_dir, "embedded_h5ad")

    os.makedirs(inputs_dir, exist_ok=True)
    os.makedirs(cfg_dir, exist_ok=True)
    os.makedirs(out_h5ad_dir, exist_ok=True)
    
    model_name = f"Tx1-{args.model_size}"
    args.model_dir = "/workspace/benchmark/Tahoe-x1/" + args.model_size + "-model"

    # -------- CRISPR --------
    g_in = os.path.join(inputs_dir, "Chem2Gen_K562_CRISPR_tx1_input.h5ad")
    g_out = os.path.join(out_h5ad_dir, "Chem2Gen_K562_CRISPR_tx1_emb.h5ad")
    g_cfg = os.path.join(cfg_dir, "predict_crispr.yaml")

    g_all, g_kept = build_h5ad_for_tahoe(
        g_counts, g_meta, var_names, symbol2ens, g_in,
        gene_id_key=args.gene_id_key, cell_type_key=args.cell_type_key,
        make_sparse=args.make_sparse
    )
    print(f"[Input] CRISPR genes: {g_all} -> {g_kept} mapped")

    write_predict_yaml(
        g_cfg, model_name=model_name,
        hf_repo_id=args.hf_repo_id, hf_model_size=args.model_size,
        adata_input=g_in, adata_output=g_out,
        cell_type_key=args.cell_type_key, gene_id_key=args.gene_id_key,
        seq_len_dataset=args.seq_len_dataset,
        return_gene_embeddings=args.return_gene_embeddings,
        model_dir=args.model_dir
    )

    run_predict_embeddings(g_cfg, batch_size=args.batch_size)

    export_embedding(
        g_out, model_name,
        os.path.join(args.out_dir, f"TahoeX1_{args.model_size}_CRISPR_embeddings.npy")
    )

    # -------- Drug --------
    c_in = os.path.join(inputs_dir, "Chem2Gen_K562_Drug_tx1_input.h5ad")
    c_out = os.path.join(out_h5ad_dir, "Chem2Gen_K562_Drug_tx1_emb.h5ad")
    c_cfg = os.path.join(cfg_dir, "predict_drug.yaml")

    c_all, c_kept = build_h5ad_for_tahoe(
        c_counts, c_meta, var_names, symbol2ens, c_in,
        gene_id_key=args.gene_id_key, cell_type_key=args.cell_type_key,
        make_sparse=args.make_sparse
    )
    print(f"[Input] Drug genes: {c_all} -> {c_kept} mapped")

    write_predict_yaml(
        c_cfg, model_name=model_name,
        hf_repo_id=args.hf_repo_id, hf_model_size=args.model_size,
        adata_input=c_in, adata_output=c_out,
        cell_type_key=args.cell_type_key, gene_id_key=args.gene_id_key,
        seq_len_dataset=args.seq_len_dataset,
        return_gene_embeddings=args.return_gene_embeddings,
        model_dir=args.model_dir
    )

    run_predict_embeddings(c_cfg, batch_size=args.batch_size)

    export_embedding(
        c_out, model_name,
        os.path.join(args.out_dir, f"TahoeX1_{args.model_size}_Drug_embeddings.npy")
    )


if __name__ == "__main__":
    main()
