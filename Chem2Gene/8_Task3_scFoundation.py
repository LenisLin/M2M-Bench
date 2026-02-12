import os
import re
import glob
import argparse
import subprocess
import shutil

import numpy as np
import pandas as pd
import torch

# anndata / scipy 是做 h5ad & sparse 的最稳选择
import anndata as ad
from scipy import sparse

# python 9_Task3_scFoundation.py \
#   --scfoundation_root /home/lenislin/Experiment/projects/Chem2Gen/benchmark/scFoundation \
#   --eval_data_dir /mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets/Evaluation_Set_K562 \
#   --output_dir /mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results \
#   --version rde \
#   --tgthighres f1 \
#   --pool_type all \
#   --pre_normalized F


def load_counts_and_genes(counts_pt: str, shared_var_names_csv: str):
    counts = torch.load(counts_pt, map_location="cpu")
    if torch.is_tensor(counts):
        counts = counts.cpu().numpy()
    genes = pd.read_csv(shared_var_names_csv)["gene_symbol"].astype(str).tolist()
    assert counts.shape[1] == len(genes), (
        f"Shape mismatch: counts has {counts.shape[1]} genes, "
        f"but shared_var_names has {len(genes)}"
    )
    return counts, genes


def load_scfoundation_gene_list(gene_list_tsv: str, to_upper: bool = False):
    """
    scFoundation 的 19,264 基因列表文件（OS_scRNA_gene_index.19264.tsv）
    文件具体列数可能变化；通常第一列是 gene symbol。
    """
    df = pd.read_csv(gene_list_tsv, sep="\t", header=None)
    genes = df.iloc[:, 0].astype(str).tolist()
    if to_upper:
        genes = [g.upper() for g in genes]
    return genes


def align_to_scfoundation_genes(counts: np.ndarray, genes: list, sf_genes: list, to_upper: bool = False):
    """
    将 (n_cells, n_input_genes) 的 counts 对齐到 (n_cells, 19264) 的 scFoundation 基因顺序。
    用稀疏 COO 重映射，避免构建巨大 dense 中间矩阵。
    """
    if to_upper:
        genes = [g.upper() for g in genes]
        sf_genes = [g.upper() for g in sf_genes]

    sf_index = {g: i for i, g in enumerate(sf_genes)}

    input_cols = []
    target_cols = []
    for j, g in enumerate(genes):
        if g in sf_index:
            input_cols.append(j)
            target_cols.append(sf_index[g])

    input_cols = np.array(input_cols, dtype=np.int32)
    target_cols = np.array(target_cols, dtype=np.int32)

    if input_cols.size == 0:
        raise RuntimeError("No overlap between input genes and scFoundation 19264 gene list.")

    # 构建 input 的稀疏矩阵
    X = sparse.csr_matrix(counts)

    # 先取交集列（保持稀疏）
    X_sub = X[:, input_cols].tocoo()

    # 将子矩阵列索引 remap 到 scFoundation 全量列坐标
    new_cols = target_cols[X_sub.col]

    X_full = sparse.coo_matrix(
        (X_sub.data, (X_sub.row, new_cols)),
        shape=(counts.shape[0], len(sf_genes)),
        dtype=np.float32
    ).tocsr()

    return X_full, input_cols.size


def write_h5ad(X_csr, var_names, out_h5ad: str):
    obs = pd.DataFrame(index=[f"cell_{i}" for i in range(X_csr.shape[0])])

    # scFoundation/社区使用中经常需要 total_count（每 cell 总 reads/UMI）
    obs["total_count"] = np.asarray(X_csr.sum(axis=1)).reshape(-1)

    var = pd.DataFrame(index=pd.Index(var_names, name="gene_symbol"))

    adata = ad.AnnData(X=X_csr, obs=obs, var=var)
    adata.write_h5ad(out_h5ad)
    return out_h5ad


def run_get_embedding(
    scfoundation_root: str,
    task_name: str,
    input_type: str,
    output_type: str,
    pool_type: str,
    tgthighres: str,
    version: str,
    pre_normalized: str,
    data_path: str,
    save_path: str,
    python_bin: str = "python",
):
    get_embedding_py = os.path.join(scfoundation_root, "model", "get_embedding.py")
    if not os.path.isfile(get_embedding_py):
        raise FileNotFoundError(f"Cannot find: {get_embedding_py}")

    os.makedirs(save_path, exist_ok=True)

    cmd = [
        python_bin, get_embedding_py,
        "--task_name", task_name,
        "--input_type", input_type,
        "--output_type", output_type,
        "--pool_type", pool_type,
        "--tgthighres", tgthighres,
        "--data_path", data_path,
        "--save_path", save_path,
        "--pre_normalized", pre_normalized,
        "--version", version,
    ]

    print("\n>>> Running scFoundation get_embedding.py:")
    print(" ".join(cmd))

    # 为了让 get_embedding.py 的相对路径逻辑更稳定，建议 cwd= scFoundation/model
    workdir = os.path.join(scfoundation_root, "model")
    subprocess.run(cmd, check=True, cwd=workdir)


def pick_latest_embedding_file(save_path: str, task_name: str):
    """
    scFoundation 输出文件命名在不同版本/参数下可能略有差异；
    这里用通用规则：在 save_path 里找包含 task_name 和 cell_embedding 的 .npy。
    """
    patterns = [
        os.path.join(save_path, f"{task_name}*cell_embedding*.npy"),
        os.path.join(save_path, f"*{task_name}*cell_embedding*.npy"),
        os.path.join(save_path, f"{task_name}*.npy"),
    ]
    candidates = []
    for p in patterns:
        candidates.extend(glob.glob(p))

    candidates = list(set(candidates))
    if not candidates:
        raise FileNotFoundError(
            f"No embedding .npy found under {save_path} for task_name={task_name}. "
            f"Please check scFoundation stdout for the exact output filename."
        )

    candidates.sort(key=lambda x: os.path.getmtime(x), reverse=True)
    return candidates[0]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--scfoundation_root", type=str, default="/home/lenislin/Experiment/projects/Chem2Gen/benchmark/scFoundation",
                    help="Path to cloned scFoundation repo root.")
    ap.add_argument("--eval_data_dir", type=str, default="/mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets/Evaluation_Set_K562",
                    help="Your Chem2Gen Evaluation_Set_K562 directory.")
    ap.add_argument("--output_dir", type=str, default="/mnt/NAS_21T/ProjectData/Chem2Gen/Model_Evaluation_Results",
                    help="Output dir to save final embeddings.")
    ap.add_argument("--python_bin", type=str, default="python",
                    help="Python executable used to run get_embedding.py")

    ap.add_argument("--version", type=str, default="rde",
                    help="scFoundation version flag (e.g., rde or ce).")
    ap.add_argument("--tgthighres", type=str, default="f1",
                    help="tgthighres flag used by get_embedding.py")
    ap.add_argument("--pool_type", type=str, default="all",
                    help="pool_type flag used by get_embedding.py")
    ap.add_argument("--pre_normalized", type=str, default="F",
                    help="pre_normalized flag used by get_embedding.py (T/F).")

    ap.add_argument("--to_upper", action="store_true",
                    help="Uppercase gene symbols before matching to scFoundation list.")

    args = ap.parse_args()

    os.makedirs(args.output_dir, exist_ok=True)

    # scFoundation 19264 gene list
    sf_gene_list_tsv = os.path.join(args.scfoundation_root, "OS_scRNA_gene_index.19264.tsv")
    if not os.path.isfile(sf_gene_list_tsv):
        raise FileNotFoundError(
            f"Cannot find {sf_gene_list_tsv}. "
            f"Make sure you cloned scFoundation correctly."
        )
    sf_genes = load_scfoundation_gene_list(sf_gene_list_tsv, to_upper=args.to_upper)
    print(f">>> Loaded scFoundation gene list: {len(sf_genes)} genes")

    shared_var = os.path.join(args.eval_data_dir, "shared_var_names.csv")

    sources = [
        ("CRISPR", os.path.join(args.eval_data_dir, "CRISPR_counts.pt")),
        ("Drug",   os.path.join(args.eval_data_dir, "Drug_counts.pt")),
    ]

    for source_name, pt_path in sources:
        print(f"\n==============================")
        print(f">>> Processing source: {source_name}")
        print(f"==============================")

        counts, genes = load_counts_and_genes(pt_path, shared_var)
        print(f"    Raw counts shape: {counts.shape} | genes: {len(genes)}")

        X_full, n_overlap = align_to_scfoundation_genes(
            counts=counts,
            genes=genes,
            sf_genes=sf_genes,
            to_upper=args.to_upper
        )
        print(f"    Gene overlap with scFoundation list: {n_overlap} / {len(genes)}")
        print(f"    Aligned matrix: {X_full.shape} (should be n_cells x 19264)")

        # 写出 aligned h5ad（作为 get_embedding 的输入）
        work_dir = os.path.join(args.output_dir, "_scfoundation_workdir", source_name)
        os.makedirs(work_dir, exist_ok=True)

        aligned_h5ad = os.path.join(work_dir, f"aligned_{source_name}.h5ad")
        write_h5ad(X_full, sf_genes, aligned_h5ad)
        print(f"    ✅ Wrote aligned h5ad: {aligned_h5ad}")

        # 调 get_embedding.py
        task_name = f"Chem2Gen_K562_{source_name}"
        save_path = os.path.join(work_dir, "scfoundation_out")
        run_get_embedding(
            scfoundation_root=args.scfoundation_root,
            task_name=task_name,
            input_type="singlecell",
            output_type="cell",
            pool_type=args.pool_type,
            tgthighres=args.tgthighres,
            version=args.version,
            pre_normalized=args.pre_normalized,
            data_path=aligned_h5ad,
            save_path=save_path,
            python_bin=args.python_bin,
        )

        # 找到输出 embedding，并复制到你的 OUTPUT_DIR，统一命名
        emb_file = pick_latest_embedding_file(save_path, task_name)
        final_name = os.path.join(args.output_dir, f"scFoundation_{source_name}_embeddings.npy")
        shutil.copy2(emb_file, final_name)

        print(f"    ✅ Found embedding file: {emb_file}")
        print(f"    ✅ Copied to: {final_name}")

    print("\nAll done.")


if __name__ == "__main__":
    main()
