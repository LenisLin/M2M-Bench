# SCRIPT_HEADER_CONTRACT
# Script: scripts/s1_task1_internal_metrics.py
# Purpose: 计算 Task1 internal 组级一致性与检索指标 (Modality Concordance)
# Inputs:
#   - Task1 Snapshot: data/task1_snapshot_v1/
# Outputs:
#   - task1_retrieval_per_query.parquet: runs/<run_id>/s1_task1_internal_metrics/
#   - task1_retrieval_summary.csv: runs/<run_id>/s1_task1_internal_metrics/
#   - task1_chance_identity_check.csv: runs/<run_id>/s1_task1_internal_metrics/
#   - task1_leaderboard_long.csv: runs/<run_id>/s1_task1_internal_metrics/
#   - task1_attrition.csv: runs/<run_id>/s1_task1_internal_metrics/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s1_task1_internal_metrics/
# Config Dependencies:
#   - config/config.yaml
# Execution:
#   - python scripts/s1_task1_internal_metrics.py --run-id <run_id> --seed 619
# Failure Modes:
#   - LOO strict_recompute violated -> exit non-zero
#   - Chance identity check tolerance > 1e-12 -> exit non-zero
# Last Updated: 2026-03-02