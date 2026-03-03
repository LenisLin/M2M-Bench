# SCRIPT_HEADER_CONTRACT
# Script: scripts/s3_build_task2_snapshot.py
# Purpose: 派生 Task2 LINCS pairs 并重建 K562 数据以构建 Task2 Snapshot
# Inputs:
#   - Task1 Snapshot (LINCS): data/task1_snapshot_v1/lincs/
#   - Raw K562 Data: /mnt/NAS_21T/ProjectData/Chem2Gen/Benchmark_Datasets/Evaluation_Set_K562
# Outputs:
#   - task2_lincs_pairs.csv: data/task2_snapshot_v1/lincs/
#   - pair_list.parquet, delta_meta.csv, fm_delta_meta.csv: data/task2_snapshot_v1/k562/derived/
#   - task2_post_build_inventory.csv: runs/<run_id>/s3_build_task2_snapshot/
#   - task2_pairs_coverage.csv: runs/<run_id>/s3_build_task2_snapshot/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Writes to frozen snapshot directory: data/task2_snapshot_v1/ (One-time build)
#   - Creates isolated run directory: runs/<run_id>/s3_build_task2_snapshot/
# Config Dependencies:
#   - config/config.yaml
# Execution:
#   - python scripts/s3_build_task2_snapshot.py --run-id <run_id> --seed 619
# Failure Modes:
#   - FM valid_mask drops leading to zero pairs -> exit non-zero
#   - Missing required input files -> exit non-zero
# Last Updated: 2026-03-02