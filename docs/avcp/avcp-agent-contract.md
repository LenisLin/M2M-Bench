# AVCP Agent Contract (CODEX) — v1.0

**Scope**  
This contract defines non-negotiable operational constraints for CODEX when working inside the AVCP-template workspace on **M2M-Bench** (ecological benchmark; report facts only). It exists to prevent cross-project contamination, enforce reproducibility, and ensure third-party-auditable outputs.

Repository reference (informational): https://github.com/LenisLin/AVCP-template.git

---

## 0) Guiding Principles

1) **Zero contamination**: never write outside the confirmed PROJECT_ROOT.  
2) **Run isolation**: every execution writes to a unique run directory.  
3) **Falsifiability-first evidence**: audit/proof fields must be data-derivable or accompanied by verifiable recomputation checks.  
4) **Report facts, not narratives**: no external MOA/pathway databases; no causal claims; no regression “explanations” unless explicitly requested as optional appendix.

---

## 1) Definitions

- **PROJECT_ROOT**: the single directory that contains the AVCP workspace for this run.  
- **Root sentinel**: a file stored at PROJECT_ROOT that proves the directory is the intended workspace.
  - Required filename (choose one and standardize): `.avcp_root` (preferred)  
  - Required content: `project_name`, `created_at`, optional `owner`.
- **OUTPUT_ROOT**: the base output directory for *this run only*.
- **RUN_ID**: unique identifier, e.g. `m2m_YYYYMMDD_HHMMSS`.
- **Shared directory**: any directory not under `OUTPUT_ROOT`, especially `data/interim_viz/` top-level paths.
- **Blocked output**: a required input is missing; outputs must be explicitly non-usable.

---

## 2) Hard Constraints (Must)

### 2.1 Root safety check (fail-fast)
Before any read/write/compute:
- Print:
  - `pwd`
  - `PROJECT_ROOT` (absolute path)
  - if git: `git rev-parse --show-toplevel`, `git rev-parse HEAD` (or explicitly “git unavailable”)
- Verify:
  1) Root sentinel file exists at PROJECT_ROOT
  2) `pwd` is inside PROJECT_ROOT
  3) All file writes resolve under OUTPUT_ROOT
- If any check fails: **exit with non-zero code** and do not write artifacts.

### 2.2 Run-isolated outputs only
All outputs must be written under:

`OUTPUT_ROOT = <PROJECT_ROOT>/data/interim_viz/m2m_runs/<RUN_ID>/`

Prohibited:
- writing to `<PROJECT_ROOT>/data/interim_viz/` top-level directly
- writing outside PROJECT_ROOT
- overwriting previous RUN_ID outputs unless explicitly instructed

### 2.3 Mandatory dry-run phase
Before execution that writes files:
- Provide a **dry-run**:
  - PROJECT_ROOT / OUTPUT_ROOT / RUN_ID
  - Full list of intended output files (absolute + relative)
  - Input file paths (absolute + relative)
- Wait for explicit human confirmation **in the conversation** before running.

### 2.4 Falsifiable audit evidence (no “constant-by-task” proof fields)
For any audit/proof fields (examples):
- `query_in_gallery_bool`
- `query_in_centroid_members_bool`
- `centroid_recompute_proof`
- `loo_policy`

Rules:
- These fields **must not** be set purely by `if task == ...` logic.
- They must be:
  1) computed from data that can be reloaded and rechecked (e.g. centroid member lists / hashes), **or**
  2) validated by an explicit recomputation table (random sample) that compares stored vs recomputed values and must match 100%.
- Proof documents must contain **one** authoritative definition (no contradictory “conceptual” vs “actual” definitions).

### 2.5 Assertions must include rules text
`audit_assertions.json` must include, for every critical assertion:
- `name`, `pass`, `details`
- `details.rules`: human-readable AND-conditions (as strings or list)
- `counterexamples` (first K rows) when failing

### 2.6 Blocked behavior must prevent mis-use
If required inputs are missing:
- Write only a blocked meta file with `blocked_reason`, and either:
  - do **not** write the corresponding CSV/Parquet outputs, **or**
  - write outputs with:
    - `analysis_status="blocked"` column (constant)
    - all derived fields set to NA (not copied from another column)
    - `is_valid_* = false`

---

## 3) Forbidden Behaviors (Must Not)

1) **No cross-project writes**: never write into other repos, sibling folders, or shared global caches unless explicitly requested.
2) **No silent fallback**: do not substitute missing inputs with other inputs while still producing “normal-looking” outputs.
3) **No metric mixing**: never aggregate across `metric_family` unless explicitly intended and clearly labeled.
4) **No external annotation dependencies**: do not use MOA/pathway databases or external enrichment tools unless explicitly requested as optional appendix.
5) **No retroactive editing of evidence**: do not manually edit evidence outputs; only regenerate via code and record the run.

---

## 4) Required Deliverables Per Run (Minimum)

### 4.1 Run manifest
`run_manifest.json` must include:
- PROJECT_ROOT, OUTPUT_ROOT, RUN_ID
- timestamps
- script paths + script sha256 (at minimum for the executed scripts)
- git commit if available
- full list of input file paths
- full list of output files generated

### 4.2 Evidence bundle (third-party review)
- `evidence_report.md`: must contain *only facts* (no “I think”)
- `evidence_manifest.md`: file sha256, sizes, and 3 one-click verification commands  
  - Prefer `duckdb`; if unavailable, provide python equivalents + actual outputs.

### 4.3 Sidecar meta for key CSV tables
Every released CSV must have a `_meta.json` with:
- `script`, `script_sha256`, `run_id`
- `is_valid_*` flags
- `blocked_reason` if applicable
- `group_key_schema` and any policy bindings (e.g., headline policy)

---

## 5) Communication Protocol (AVCP-style outputs)

Every significant response should follow this structure:

- `[STATE SNAPSHOT]`: where you are, what files touched, what paths used  
- `[PLAN]`: numbered steps with clear outputs  
- `[PATCH SET]`: files changed/added, with exact paths  
- `[TEST]`: commands run + results  
- `[EVIDENCE]`: paths to artifacts + minimal excerpts (tables/JSON snippets)  
- `[RISKS & ASSUMPTIONS]`: only if needed; must be factual and scoped  
- `[NEXT]`: optional

---

## 6) Acceptance Criteria (Run Pass/Fail)

A run is acceptable **only if**:
1) Root sentinel checks pass and all writes are under OUTPUT_ROOT  
2) No silent fallback (blocked outputs are explicitly blocked/non-usable)  
3) Falsifiable evidence rules are satisfied (no constant-by-task proof fields)  
4) `audit_assertions.json` critical assertions pass (or policy-gated warnings are correctly bound to policy outputs)  
5) Report/manifest contain sufficient hashes and one-click checks

---

## 7) Escalation

If any constraint is impossible due to environment limitations (e.g., missing duckdb):
- State the limitation explicitly
- Provide an equivalent, reproducible alternative (e.g., python query + captured output)
- Do not proceed with a degraded approach if it risks producing misleading “valid-looking” outputs.

---

**End of contract.**
