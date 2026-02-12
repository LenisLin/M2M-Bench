import os
import argparse
import subprocess
import numpy as np
import pandas as pd
import torch
import shlex
import scanpy as sc
from anndata import AnnData
from scipy import sparse


# -------------------------
# Chem2Gen loader (minimal)
# -------------------------
def load_chem2gen_eval(eval_dir: str):
    var_names = pd.read_csv(os.path.join(eval_dir, "shared_var_names.csv"))["gene_symbol"].tolist()

    g_counts = torch.load(os.path.join(eval_dir, "CRISPR_counts.pt"), map_location="cpu")
    g_meta = pd.read_csv(os.path.join(eval_dir, "CRISPR_meta.csv"), index_col=0)

    c_counts = torch.load(os.path.join(eval_dir, "Drug_counts.pt"), map_location="cpu")
    c_meta = pd.read_csv(os.path.join(eval_dir, "Drug_meta.csv"), index_col=0)

    # basic sanity
    assert g_counts.shape[0] == g_meta.shape[0], f"CRISPR counts/meta mismatch: {g_counts.shape[0]} vs {g_meta.shape[0]}"
    assert c_counts.shape[0] == c_meta.shape[0], f"Drug counts/meta mismatch: {c_counts.shape[0]} vs {c_meta.shape[0]}"
    assert g_counts.shape[1] == len(var_names), f"CRISPR genes mismatch: {g_counts.shape[1]} vs {len(var_names)}"
    assert c_counts.shape[1] == len(var_names), f"Drug genes mismatch: {c_counts.shape[1]} vs {len(var_names)}"

    return var_names, g_counts, g_meta, c_counts, c_meta


def build_h5ad_from_counts(
    counts: torch.Tensor,
    obs_df: pd.DataFrame,
    var_gene_symbols: list[str],
    out_h5ad: str,
    gene_column: str = "gene_name",
    make_sparse: bool = True,
):
    """
    Build an AnnData that matches STATE emb.transform requirement:
      - CSR sparse matrix in adata.X
      - gene name column exists in adata.var (default: gene_name)
    """
    os.makedirs(os.path.dirname(out_h5ad), exist_ok=True)

    # Ensure float32 (counts may be float already in your spec)
    if counts.dtype != torch.float32:
        counts = counts.float()
    X = counts.numpy()  # (n_cells, n_genes)

    if make_sparse:
        # STATE requires CSR (per repo notes)
        X = sparse.csr_matrix(X, dtype=np.float32)
    else:
        # not recommended for STATE; keep for debugging only
        X = np.asarray(X, dtype=np.float32)

    var = pd.DataFrame(index=pd.Index(var_gene_symbols, name="gene_symbol"))
    var[gene_column] = var_gene_symbols  # critical for STATE

    adata = AnnData(X=X, obs=obs_df.copy(), var=var)
    adata.var_names_make_unique()
    adata.obs_names = obs_df.index.astype(str)

    # Write
    adata.write_h5ad(out_h5ad)
    return out_h5ad


def run_state_emb_transform(
    state_exe: str,
    model_folder: str,
    input_h5ad: str,
    output_h5ad: str,
    checkpoint: str | None = None,
    gene_column: str = "gene_name",
    extra_args: list[str] | None = None,
):
    """
    Calls:
      state emb transform --model-folder ... [--checkpoint ...] --input ... --output ... --gene-column ...
    """
    os.makedirs(os.path.dirname(output_h5ad), exist_ok=True)

    prefix = shlex.split(state_exe)

    cmd = prefix + [
        "emb", "transform",
        "--model-folder", model_folder,
        "--input", input_h5ad,
        "--output", output_h5ad
    ]
    if checkpoint:
        cmd += ["--checkpoint", checkpoint]
    if extra_args:
        cmd += extra_args

    print("\n[STATE] Running:\n  " + " ".join(cmd) + "\n")
    subprocess.run(cmd, check=True)


def extract_state_embedding_to_npy(
    output_h5ad: str,
    out_npy: str,
    prefer_key: str = "X_state",
):
    """
    Load output h5ad and extract embedding (default: .obsm['X_state']).
    Falls back to any obsm key containing 'state' if needed.
    """
    adata = sc.read_h5ad(output_h5ad)

    key = None
    if prefer_key in adata.obsm:
        key = prefer_key
    else:
        # fallback heuristic
        candidates = [k for k in adata.obsm.keys() if "state" in k.lower()]
        if len(candidates) > 0:
            key = candidates[0]

    if key is None:
        raise KeyError(
            f"Cannot find STATE embedding in .obsm. Available keys: {list(adata.obsm.keys())}"
        )

    emb = np.asarray(adata.obsm[key])
    os.makedirs(os.path.dirname(out_npy), exist_ok=True)
    np.save(out_npy, emb)
    print(f"[STATE] Saved embeddings: {out_npy} | key={key} | shape={emb.shape}")
    return emb.shape, key


def main():
    parser = argparse.ArgumentParser("Chem2Gen → STATE(SE) embedding extraction")

    parser.add_argument("--eval_dir", type=str,
                        default="/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets/Evaluation_Set_K562")
    parser.add_argument("--out_dir", type=str,
                        default="/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results")
    parser.add_argument("--work_dir", type=str,
                        default="/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/STATE/workdir_K562")

    # STATE CLI + model
    parser.add_argument("--state_exe", type=str, default="state", #  "uvx -q --from git+https://github.com/ArcInstitute/state@main"
                        help="Executable name/path. If you installed via `uv`, you can also use `uv run state` manually.")
    parser.add_argument("--model_folder", type=str, default="/home/lenislin/Experiment/projects/Chem2Gen/benchmark/state/arcinstitute/SE-600M",
                        help="Local path to the SE model folder (e.g., downloaded arcinstitute/SE-600M).")
    parser.add_argument("--checkpoint", type=str, default="/home/lenislin/Experiment/projects/Chem2Gen/benchmark/state/arcinstitute/SE-600M/se600m_epoch16.ckpt",
                        help="Optional checkpoint path (e.g., se600m_epoch16.ckpt). If omitted, STATE may use defaults in the model folder.")
    parser.add_argument("--gene_column", type=str, default="gene_name",
                        help="Column in adata.var storing gene symbols (STATE requires gene_name unless overridden).")

    parser.add_argument("--make_sparse", action=argparse.BooleanOptionalAction, default=True,
                        help="Store AnnData.X as CSR sparse matrix (recommended/required for STATE).")

    args = parser.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.work_dir, exist_ok=True)

    # 1) Load Chem2Gen eval tensors + meta
    var_names, g_counts, g_meta, c_counts, c_meta = load_chem2gen_eval(args.eval_dir)

    # 2) Build input h5ad (CSR + gene_name in var)
    inputs_dir = os.path.join(args.work_dir, "inputs")
    g_in = os.path.join(inputs_dir, "Chem2Gen_K562_CRISPR.h5ad")
    c_in = os.path.join(inputs_dir, "Chem2Gen_K562_Drug.h5ad")

    build_h5ad_from_counts(g_counts, g_meta, var_names, g_in,
                           gene_column=args.gene_column, make_sparse=args.make_sparse)
    build_h5ad_from_counts(c_counts, c_meta, var_names, c_in,
                           gene_column=args.gene_column, make_sparse=args.make_sparse)

    # 3) Run STATE emb transform (writes embedding into output h5ad)
    out_h5ad_dir = os.path.join(args.work_dir, "embedded_h5ad")
    g_out_h5ad = os.path.join(out_h5ad_dir, "Chem2Gen_K562_CRISPR_state.h5ad")
    c_out_h5ad = os.path.join(out_h5ad_dir, "Chem2Gen_K562_Drug_state.h5ad")

    run_state_emb_transform(args.state_exe, args.model_folder, g_in, g_out_h5ad,
                            checkpoint=args.checkpoint, gene_column=args.gene_column)
    run_state_emb_transform(args.state_exe, args.model_folder, c_in, c_out_h5ad,
                            checkpoint=args.checkpoint, gene_column=args.gene_column)

    # 4) Extract to .npy (cell order preserved)
    g_out_npy = os.path.join(args.out_dir, "STATE_CRISPR_embeddings.npy")
    c_out_npy = os.path.join(args.out_dir, "STATE_Drug_embeddings.npy")

    extract_state_embedding_to_npy(g_out_h5ad, g_out_npy, prefer_key="X_state")
    extract_state_embedding_to_npy(c_out_h5ad, c_out_npy, prefer_key="X_state")


if __name__ == "__main__":
    main()
