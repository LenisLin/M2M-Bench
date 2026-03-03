# SCRIPT_HEADER_CONTRACT
# Script: scripts/s2_task1_cross_metrics.py
# Purpose: 计算 Task1 cross matched pairs 一致性（跨数据集）
# Inputs:
#   - Task1 Snapshot: data/task1_snapshot_v1/cross_contract/
# Outputs:
#   - task1_cross_retrieval_per_query.parquet: runs/<run_id>/s2_task1_cross_metrics/
#   - task1_cross_retrieval_summary.csv: runs/<run_id>/s2_task1_cross_metrics/
#   - task1_cross_chance_identity_check.csv: runs/<run_id>/s2_task1_cross_metrics/
#   - task1_cross_leaderboard_long.csv: runs/<run_id>/s2_task1_cross_metrics/
#   - task1_cross_alignment_proof.csv: runs/<run_id>/s2_task1_cross_metrics/
#   - task1_cross_attrition.csv: runs/<run_id>/s2_task1_cross_metrics/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s2_task1_cross_metrics/
# Config Dependencies:
#   - config/config.yaml
# Execution:
#   - python scripts/s2_task1_cross_metrics.py --run-id <run_id> --seed 619
# Failure Modes:
#   - Eligibility gate mismatch (target alignment failed) -> exit non-zero
#   - Chance identity check tolerance > 1e-12 -> exit non-zero
# Last Updated: 2026-03-02