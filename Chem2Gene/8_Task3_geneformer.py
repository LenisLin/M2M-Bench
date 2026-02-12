#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import os
import argparse
import pickle
import shutil
import tempfile
import numpy as np
import pandas as pd
import torch
from scipy import sparse
import anndata as ad

# Geneformer
from geneformer import TranscriptomeTokenizer, EmbExtractor


def load_chem2gen_eval(eval_dir: str):
    var_names = pd.read_csv(os.path.join(eval_dir, "shared_var_names.csv"))["gene_symbol"].astype(str).tolist()

    g_counts = torch.load(os.path.join(eval_dir, "CRISPR_counts.pt"), map_location="cpu")
    g_meta = pd.read_csv(os.path.join(eval_dir, "CRISPR_meta.csv"), index_col=0)

    c_counts = torch.load(os.path.join(eval_dir, "Drug_counts.pt"), map_location="cpu")
    c_meta = pd.read_csv(os.path.join(eval_dir, "Drug_meta.csv"), index_col=0)

    assert g_counts.shape[0] == len(g_meta), "CRISPR counts rows != CRISPR_meta rows"
    assert c_counts.shape[0] == len(c_meta), "Drug counts rows != Drug_meta rows"
    assert g_counts.shape[1] == len(var_names), "CRISPR gene cols != shared_var_names"
    assert c_counts.shape[1] == len(var_names), "Drug gene cols != shared_var_names"

    return var_names, g_counts, g_meta, c_counts, c_meta


def load_gene_symbol_to_ensembl(pkl_path: str) -> dict:
    with open(pkl_path, "rb") as f:
        d = pickle.load(f)
    if not isinstance(d, dict):
        raise TypeError(f"Expected dict in {pkl_path}, got {type(d)}")
    return d


def build_gene_mapping(shared_gene_symbols: list[str], symbol2ens: dict):
    keep_cols, kept_symbols, kept_ens = [], [], []
    for j, sym in enumerate(shared_gene_symbols):
        ens = symbol2ens.get(sym)
        if ens is not None and isinstance(ens, str) and ens.startswith("ENSG"):
            keep_cols.append(j)
            kept_symbols.append(sym)
            kept_ens.append(ens)

    keep_cols = np.asarray(keep_cols, dtype=np.int64)
    if keep_cols.size == 0:
        raise RuntimeError("No genes can be mapped to Ensembl IDs using the provided dictionary.")
    return keep_cols, kept_symbols, kept_ens


def make_h5ad_from_counts(counts_t: torch.Tensor,
                          meta: pd.DataFrame,
                          keep_cols: np.ndarray,
                          kept_symbols: list[str],
                          kept_ens: list[str],
                          out_h5ad: str):
    """
    Geneformer tokenizer for h5ad requires:
      - adata.var["ensembl_id"] exists
      - adata.obs["n_counts"] exists
    """
    os.makedirs(os.path.dirname(out_h5ad), exist_ok=True)

    X = counts_t[:, keep_cols].cpu().numpy().astype(np.float32)
    X = sparse.csr_matrix(X, dtype=np.float32)

    obs = meta.copy()
    obs["n_counts"] = np.asarray(X.sum(axis=1)).reshape(-1).astype(np.float32)
    if "filter_pass" not in obs.columns:
        obs["filter_pass"] = True

    var = pd.DataFrame(index=pd.Index(kept_symbols, name="gene_symbol"))
    var["ensembl_id"] = kept_ens

    adata = ad.AnnData(X=X, obs=obs, var=var)
    adata.obs_names = meta.index.astype(str)
    adata.write_h5ad(out_h5ad)

    return out_h5ad


def tokenize_h5ad_to_dataset_singlefile(
    h5ad_path: str,
    token_out_dir: str,
    output_prefix: str,
    nproc: int,
    chunk_size: int,
):
    """
    方案A：确保每次 tokenization 只看到一个 h5ad 文件。
    做法：把该 h5ad 复制到一个临时目录，然后把 data_directory 指向该临时目录。

    输出: {token_out_dir}/{output_prefix}.dataset
    """
    os.makedirs(token_out_dir, exist_ok=True)

    tokenizer = TranscriptomeTokenizer(
        model_version="V2",
        nproc=nproc,
        chunk_size=chunk_size,
        use_h5ad_index=True,
        custom_attr_name_dict={
            "target": "target",
            "benchmark_group": "benchmark_group",
            "specificity_tier": "specificity_tier",
            "clean_target_mapped": "clean_target_mapped",
        },
    )

    # 临时目录仅放这一个 h5ad
    tmp_dir = tempfile.mkdtemp(prefix=f"geneformer_tokenize_{output_prefix}_")
    try:
        tmp_h5ad = os.path.join(tmp_dir, os.path.basename(h5ad_path))
        shutil.copy2(h5ad_path, tmp_h5ad)

        tokenizer.tokenize_data(
            data_directory=tmp_dir,
            output_directory=token_out_dir,
            output_prefix=output_prefix,
            use_generator=True,
            file_format="h5ad",
        )
    finally:
        shutil.rmtree(tmp_dir, ignore_errors=True)

    dataset_path = os.path.join(token_out_dir, f"{output_prefix}.dataset")
    if not os.path.exists(dataset_path):
        candidates = [
            p for p in os.listdir(token_out_dir)
            if p.endswith(".dataset") and p.startswith(output_prefix)
        ]
        if not candidates:
            raise FileNotFoundError(f"Tokenization finished but dataset not found under {token_out_dir}")
        dataset_path = os.path.join(token_out_dir, candidates[0])

    return dataset_path


def extract_cell_embeddings(model_dir: str,
                            dataset_path: str,
                            out_dir: str,
                            out_prefix: str,
                            forward_batch_size: int,
                            nproc: int,
                            emb_mode: str = "cls",
                            emb_layer: int = -1):
    os.makedirs(out_dir, exist_ok=True)

    embex = EmbExtractor(
        model_type="Pretrained",
        num_classes=0,
        emb_mode=emb_mode,
        max_ncells=None,           # 全量
        emb_layer=emb_layer,
        forward_batch_size=forward_batch_size,
        nproc=nproc,
        model_version="V2",
    )

    embs_df, embs = embex.extract_embs(
        model_directory=model_dir,
        input_data_file=dataset_path,
        output_directory=out_dir,
        output_prefix=out_prefix,
        output_torch_embs=True,
    )
    return embs_df, embs


def align_and_save_npy(embs_df: pd.DataFrame,
                       embs_tensor,
                       meta_index: pd.Index,
                       out_npy: str):
    meta_ids = meta_index.astype(str)

    if embs_df.index.astype(str).isin(meta_ids).all():
        emb_ids = embs_df.index.astype(str)
    else:
        candidate_cols = ["cell_id", "cell", "barcode", "obs_names", "obs_name", "sample"]
        found = None
        for col in candidate_cols:
            if col in embs_df.columns:
                found = col
                break
        if found is None:
            raise RuntimeError(
                f"Cannot find cell id in embs_df index or common columns. embs_df columns={embs_df.columns.tolist()}"
            )
        emb_ids = embs_df[found].astype(str)

    emb = embs_tensor.detach().cpu().numpy().astype(np.float32)
    emb_df = pd.DataFrame(emb, index=pd.Index(emb_ids, name="cell_id"))

    aligned = emb_df.reindex(meta_ids)
    missing = aligned.isna().any(axis=1).sum()
    if missing > 0:
        raise RuntimeError(
            f"{missing} cells in meta have no embedding after Geneformer tokenization/extraction. "
            f"This usually happens if tokenizer filtered those cells."
        )

    os.makedirs(os.path.dirname(out_npy), exist_ok=True)
    np.save(out_npy, aligned.to_numpy(dtype=np.float32))
    print(f"[Saved] {out_npy} | shape={aligned.shape}")


def main():
    ap = argparse.ArgumentParser("Chem2Gen → Geneformer cell embedding extraction (Strategy A)")
    ap.add_argument("--eval_dir", type=str,
                    default="/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets/Evaluation_Set_K562")
    ap.add_argument("--out_dir", type=str,
                    default="/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results")
    ap.add_argument("--work_dir", type=str,
                    default="/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results/Geneformer_workdir_K562")

    ap.add_argument("--gene_dict_pkl", type=str,
                    default="/home/lenislin/Experiment/projects/Chem2Gen/benchmark/Geneformer/geneformer/gene_name_id_dict_gc104M.pkl")
    ap.add_argument("--model_dir", type=str,
                    default="/home/lenislin/Experiment/projects/Chem2Gen/benchmark/Geneformer/Geneformer-V2-104M")

    ap.add_argument("--nproc", type=int, default=8)
    ap.add_argument("--chunk_size", type=int, default=512)
    ap.add_argument("--forward_batch_size", type=int, default=32)
    ap.add_argument("--emb_mode", type=str, default="cls", choices=["cls", "cell"])
    ap.add_argument("--emb_layer", type=int, default=-1)

    # 可选：强制重新tokenize（避免旧的 .dataset 混入）
    ap.add_argument("--force_retokenize", action="store_true",
                    help="If set, delete existing tokenized datasets before tokenizing.")

    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)
    os.makedirs(args.work_dir, exist_ok=True)

    # 1) load data
    shared_genes, g_counts, g_meta, c_counts, c_meta = load_chem2gen_eval(args.eval_dir)

    # 2) mapping
    symbol2ens = load_gene_symbol_to_ensembl(args.gene_dict_pkl)
    keep_cols, kept_symbols, kept_ens = build_gene_mapping(shared_genes, symbol2ens)
    print(f"[Info] Gene mapping: {len(shared_genes)} input genes -> {len(keep_cols)} mappable genes")

    # 3) build h5ad
    h5ad_root = os.path.join(args.work_dir, "h5ad")
    h5ad_crispr_dir = os.path.join(h5ad_root, "CRISPR")
    h5ad_drug_dir = os.path.join(h5ad_root, "Drug")
    os.makedirs(h5ad_crispr_dir, exist_ok=True)
    os.makedirs(h5ad_drug_dir, exist_ok=True)

    h5ad_g = os.path.join(h5ad_crispr_dir, "Chem2Gen_K562_CRISPR_geneformer_input.h5ad")
    h5ad_c = os.path.join(h5ad_drug_dir, "Chem2Gen_K562_Drug_geneformer_input.h5ad")

    if not os.path.exists(h5ad_g):
        make_h5ad_from_counts(g_counts, g_meta, keep_cols, kept_symbols, kept_ens, h5ad_g)
    if not os.path.exists(h5ad_c):
        make_h5ad_from_counts(c_counts, c_meta, keep_cols, kept_symbols, kept_ens, h5ad_c)

    # 4) tokenize (Strategy A: single-file visibility)
    token_root = os.path.join(args.work_dir, "tokenized")
    token_crispr_dir = os.path.join(token_root, "CRISPR")
    token_drug_dir = os.path.join(token_root, "Drug")
    os.makedirs(token_crispr_dir, exist_ok=True)
    os.makedirs(token_drug_dir, exist_ok=True)

    ds_g = os.path.join(token_crispr_dir, "Chem2Gen_K562_CRISPR.dataset")
    ds_c = os.path.join(token_drug_dir, "Chem2Gen_K562_Drug.dataset")

    if args.force_retokenize:
        if os.path.exists(ds_g):
            shutil.rmtree(ds_g, ignore_errors=True)
        if os.path.exists(ds_c):
            shutil.rmtree(ds_c, ignore_errors=True)

    if not os.path.exists(ds_g):
        ds_g = tokenize_h5ad_to_dataset_singlefile(
            h5ad_g, token_crispr_dir, "Chem2Gen_K562_CRISPR", args.nproc, args.chunk_size
        )
    if not os.path.exists(ds_c):
        ds_c = tokenize_h5ad_to_dataset_singlefile(
            h5ad_c, token_drug_dir, "Chem2Gen_K562_Drug", args.nproc, args.chunk_size
        )

    print(f"[Info] Tokenized datasets:\n  {ds_g}\n  {ds_c}")

    # 5) embedding extraction
    emb_dir = os.path.join(args.work_dir, "embeddings_raw")
    os.makedirs(emb_dir, exist_ok=True)

    embs_df_g, embs_g = extract_cell_embeddings(
        model_dir=args.model_dir,
        dataset_path=ds_g,
        out_dir=emb_dir,
        out_prefix="Geneformer_Chem2Gen_K562_CRISPR",
        forward_batch_size=args.forward_batch_size,
        nproc=args.nproc,
        emb_mode=args.emb_mode,
        emb_layer=args.emb_layer,
    )
    embs_df_c, embs_c = extract_cell_embeddings(
        model_dir=args.model_dir,
        dataset_path=ds_c,
        out_dir=emb_dir,
        out_prefix="Geneformer_Chem2Gen_K562_Drug",
        forward_batch_size=args.forward_batch_size,
        nproc=args.nproc,
        emb_mode=args.emb_mode,
        emb_layer=args.emb_layer,
    )

    # 6) align & save
    out_g = os.path.join(args.out_dir, "Geneformer_CRISPR_emb.npy")
    out_c = os.path.join(args.out_dir, "Geneformer_Drug_emb.npy")
    align_and_save_npy(embs_df_g, embs_g, g_meta.index, out_g)
    align_and_save_npy(embs_df_c, embs_c, c_meta.index, out_c)

    print("[Done] Geneformer embeddings extracted and saved (Strategy A).")


if __name__ == "__main__":
    main()
