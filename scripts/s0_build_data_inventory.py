# SCRIPT_HEADER_CONTRACT
# Script: scripts/s0_build_data_inventory.py
# Purpose: 盘点原始输入生态（评估数据源与支持度，不依赖 Task2 snapshot）
# Inputs:
#   - Task1 Snapshot: data/task1_snapshot_v1/
# Outputs:
#   - result_data_inventory_long.csv: runs/<run_id>/s0_build_data_inventory/
#   - data_source_manifest.csv: runs/<run_id>/s0_build_data_inventory/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s0_build_data_inventory/
# Config Dependencies:
#   - config/config.yaml
# Execution:
#   - python scripts/s0_build_data_inventory.py --run-id <run_id> --seed 619
# Failure Modes:
#   - Missing input snapshot directories -> exit non-zero
# Last Updated: 2026-03-02