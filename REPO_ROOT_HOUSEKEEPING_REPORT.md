# Repo Root Housekeeping Report

Checked on 2026-03-13 for `/home/lenislin/Experiment/projects/M2M`.

## [GITIGNORE STATUS]

- Repo-level `.gitignore` is now committed with root ignore coverage for `/runs/` and `/reports/`.
- An additional `/reports` entry was added because the current `reports` path is a root symlink, and `git check-ignore -v` resolves it via that pattern while `runs/0310` resolves via `/runs/`.
- `git ls-files runs reports` returned no indexed paths before cleanup, so no `git rm --cached` step was needed.
- Commit created: `f8407c8` `chore(gitignore): ignore runs and reports`.

## [REPO ROOT RE-AUDIT]

- Re-audited the repo root for operational cleanup, migration, reconciliation, verification, and decision artifacts.
- Checked active-reference signals across `README.md`, `project.yaml`, `config/`, `docs/`, `prompts/`, `scripts/`, `src/`, and `tests/`.
- Re-audit result:
  - `MOVE_TO_NAS_OPS`: `5`
  - `DELETE_LOCAL_REDUNDANT`: `0`
  - `REVIEW_REQUIRED`: `8`
  - `KEEP_IN_REPO`: `1`
- `runs/` and `reports` were not restructured; they were only brought under repo ignore control in line with the accepted storage strategy.

## [MISPLACED OR REDUNDANT FILES]

| Path | Classification | Basis |
| --- | --- | --- |
| `CLEANUP_AUDIT_REPORT.md` | `MOVE_TO_NAS_OPS` | Unreferenced repo-root cleanup audit artifact. |
| `cleanup_actions_proposed.csv` | `MOVE_TO_NAS_OPS` | Unreferenced machine-generated proposal table. |
| `cleanup_inventory.csv` | `MOVE_TO_NAS_OPS` | Unreferenced machine-generated inventory. |
| `REPO_HYGIENE_FINAL_REPORT.md` | `MOVE_TO_NAS_OPS` | Operational summary, not active contract/governance documentation. |
| `repo_hygiene_final_actions.csv` | `MOVE_TO_NAS_OPS` | Operational action table paired with the relocated repo hygiene report. |
| `POST_MIGRATION_FINALIZATION_REPORT.md` | `REVIEW_REQUIRED` | Still referenced by active docs. |
| `RUNS_ROOT_RECONCILIATION_REPORT.md` | `REVIEW_REQUIRED` | Still referenced by active governance docs. |
| `STORAGE_MIGRATION_EXECUTION_REPORT.md` | `REVIEW_REQUIRED` | Still referenced by active docs. |
| `authoritative_active_results_verification.csv` | `REVIEW_REQUIRED` | Still referenced by active docs/governance. |
| `path_normalization_decision.md` | `REVIEW_REQUIRED` | Still referenced by active docs. |
| `runs_root_state.csv` | `REVIEW_REQUIRED` | Supporting state table for the still-referenced runs reconciliation report. |
| `storage_migration_actions_executed.csv` | `REVIEW_REQUIRED` | Still referenced by active docs. |
| `storage_migration_verification.csv` | `REVIEW_REQUIRED` | Still referenced by active docs. |
| `docs/governance/local_storage_policy.md` | `KEEP_IN_REPO` | Active governance document stating the accepted storage contract. |

## [FILES MOVED TO NAS OPS]

- Moved to `/mnt/NAS_21T/ProjectData/M2M/manifests/local_repo_housekeeping/2026-03-13_repo_root_housekeeping/`:
  - `CLEANUP_AUDIT_REPORT.md`
  - `cleanup_actions_proposed.csv`
  - `cleanup_inventory.csv`
  - `REPO_HYGIENE_FINAL_REPORT.md`
  - `repo_hygiene_final_actions.csv`
- The move emitted permission-preservation warnings from `mv`, but file relocation completed and the source files are no longer present at repo root.

## [FILES DELETED LOCALLY]

- None.
- No repo-root artifact met the bar for a direct local delete without either losing useful ops history or forcing unrelated doc changes.

## [REVIEW_REQUIRED ITEMS]

- `POST_MIGRATION_FINALIZATION_REPORT.md`
  - Blocker: referenced by `docs/progress_update.md` and `docs/governance/state.md`.
- `RUNS_ROOT_RECONCILIATION_REPORT.md`
  - Blocker: referenced by `docs/governance/state.md`.
- `STORAGE_MIGRATION_EXECUTION_REPORT.md`
  - Blocker: referenced by `docs/progress_update.md`.
- `authoritative_active_results_verification.csv`
  - Blocker: referenced by `docs/progress_update.md` and `docs/governance/state.md`.
- `path_normalization_decision.md`
  - Blocker: referenced by `docs/progress_update.md`.
- `runs_root_state.csv`
  - Blocker: companion evidence for the still-retained runs reconciliation report.
- `storage_migration_actions_executed.csv`
  - Blocker: referenced by `docs/progress_update.md`.
- `storage_migration_verification.csv`
  - Blocker: referenced by `docs/progress_update.md`.

## [COMMITS CREATED]

- `f8407c8` `chore(gitignore): ignore runs and reports`
- `54c8fd7` `chore(cleanup): add repo root housekeeping actions`
- Additional task-output commits are recorded separately in `COMMIT_LEDGER.md`.

## [FINAL GIT STATUS]

- `.gitignore` now ignores repo-root `runs` and `reports` without altering the accepted child-symlink layout under `runs/`.
- No authoritative benchmark scripts listed in the task constraints were edited or removed by this housekeeping pass.
- No active tracked docs or code were deleted as part of repo-root cleanup.
- The only repo-root artifacts removed from the working tree were untracked operational files relocated to NAS ops.
- Unrelated pre-existing worktree changes were not staged together with task files.
