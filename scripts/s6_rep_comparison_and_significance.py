# SCRIPT_HEADER_CONTRACT
# Script: scripts/s6_rep_comparison_and_significance.py
# Purpose: 以 Gene 为 baseline 进行表示比较的显著性检验 (Paired testing) 及富集分析
# Inputs:
#   - Runs from S1, S2, S4, S5
# Outputs:
#   - sig_comparison_task1.csv: runs/<run_id>/s6_rep_comparison/
#   - sig_comparison_task2.csv: runs/<run_id>/s6_rep_comparison/
#   - task1_enrichment_cell_line.csv: runs/<run_id>/s6_rep_comparison/
#   - task1_enrichment_target.csv: runs/<run_id>/s6_rep_comparison/
#   - task2_enrichment_*.csv: runs/<run_id>/s6_rep_comparison/
#   - AVCP Artifacts: run_manifest.json, audit_assertions.json, manifest.json
# Side Effects:
#   - Creates isolated run directory: runs/<run_id>/s6_rep_comparison/
# Config Dependencies:
#   - config/config.yaml
# Execution:
#   - python scripts/s6_rep_comparison_and_significance.py --run-id <run_id> --seed 619
# Failure Modes:
#   - Missing canonical_query_uid for paired testing -> exit non-zero
#   - Wilcoxon input length mismatch -> exit non-zero
# Last Updated: 2026-03-02