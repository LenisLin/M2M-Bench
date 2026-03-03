# SCRIPT_HEADER_CONTRACT
# Script: scripts/s4_task2_group_concordance.py
# Purpose: 计算 Task2 组级一致性 (Cosine, PCC, E-distance) 比较 Chemical vs Genetic
# Inputs:
#   - Task2 Snapshot: data/task2_snapshot_v1/
# Outputs:
#   - task2_group_concordance.csv: runs/<run_id>/s4_task2_group_concordance/
#   - task2_group_attrition.csv: runs/<run_id>/s4_task2_group_concordance/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s4_task2_group_concordance/
# Config Dependencies:
#   - config/config.yaml
# Execution:
#   - python scripts/s4_task2_group_concordance.py --run-id <run_id> --seed 619
# Failure Modes:
#   - Input missing mech_key pairs -> exit non-zero
# Last Updated: 2026-03-02