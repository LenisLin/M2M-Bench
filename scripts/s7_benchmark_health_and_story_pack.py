# SCRIPT_HEADER_CONTRACT
# Script: scripts/s7_benchmark_health_and_story_pack.py
# Purpose: 基准健康度诊断 (Dominance, Gini) 及汇总 Story Pack
# Inputs:
#   - Runs from S1, S2, S4, S5
# Outputs:
#   - health_dominance_curves.csv: runs/<run_id>/s7_benchmark_health/
#   - health_gini_effective_n.csv: runs/<run_id>/s7_benchmark_health/
#   - task2_denominator_summary.csv: runs/<run_id>/s7_benchmark_health/
#   - task2_stratified_summary.csv: runs/<run_id>/s7_benchmark_health/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s7_benchmark_health/
# Config Dependencies:
#   - config/config.yaml
# Execution:
#   - python scripts/s7_benchmark_health_and_story_pack.py --run-id <run_id> --seed 619
# Failure Modes:
#   - Dominance weight computed from uncorrected/negative metric -> exit non-zero
# Last Updated: 2026-03-02