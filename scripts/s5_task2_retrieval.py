# SCRIPT_HEADER_CONTRACT
# Script: scripts/s5_task2_retrieval.py
# Purpose: 计算 Task2 实例级双向检索指标 (C2G, G2C) 及其 chance-correction
# Inputs:
#   - Task2 Snapshot: data/task2_snapshot_v1/
# Outputs:
#   - task2_retrieval_per_query.parquet: runs/<run_id>/s5_task2_retrieval/
#   - task2_retrieval_summary.csv: runs/<run_id>/s5_task2_retrieval/
#   - task2_chance_identity_check.csv: runs/<run_id>/s5_task2_retrieval/
#   - task2_retrieval_attrition.csv: runs/<run_id>/s5_task2_retrieval/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s5_task2_retrieval/
# Config Dependencies:
#   - config/config.yaml
# Execution:
#   - python scripts/s5_task2_retrieval.py --run-id <run_id> --seed 619
# Failure Modes:
#   - Chance identity check tolerance > 1e-12 -> exit non-zero
# Last Updated: 2026-03-02