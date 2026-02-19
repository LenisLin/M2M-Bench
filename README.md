# avcp-template

Template repository for AVCP-based agentic development.

<!-- AVCP:README:START -->
## AVCP Template

Research engineering starter with AVCP prompts, repo-memory docs, and enforced bridge contracts.

> [!TIP]
> AVCP = **Repo-as-Memory** for AI engineering. Keep truth in versioned files, not chat context.

AVCP (Agentic Version Control Protocol) is a practical operating model for building software with coding agents (Codex / Claude Code / Cursor) under explicit, testable contracts.

### ✨ Why AVCP

| Risk | Typical symptom | AVCP control |
|---|---|---|
| 🧠 Context amnesia | decisions disappear across sessions | `docs/state.md` + `docs/decisions.md` |
| 🎭 Hallucinated implementation | invented contracts, guessed logic | tier gates + specs-first workflow |
| 🤫 Silent failures | pipelines "succeed" with wrong outputs | fail-fast validation + manifests |
| 📄 Docs drift | README/changelog become stale | generated README + safe changelog updates |

### 🗂 Repository Layout

```text
prompts/                 # pinned system prompt
config/                  # runtime configuration source-of-truth
docs/                    # memory, constraints, decisions, contracts
src/avcp_template/       # installable package
scripts/dev/             # README/changelog maintenance tooling
tests/                   # verification suite
```

## 🧠 AVCP Operating Contracts

### 1) Agent cognition contract
- Use `prompts/AVCP_SYSTEM_PROMPT_MIN.md` as the pinned operating protocol.
- Start each coding session by loading:
  - `docs/state.md`
  - `docs/constraints.md`
  - `docs/decisions.md`
  - `docs/api_specs.md`
  - `docs/data_contracts.md`
  - `docs/avcp_guidelines.md`

### 2) Configuration contract
- Runtime paths/parameters must come from `config/config.yaml`.
- Hardcoded paths are contract violations.

### 3) Data handoff contract
- Use `src/avcp_template/io/bridge.py::save_for_r()` for Python->R handoffs.
- Require explicit primary key + `<stem>_meta.json` sidecar.

### 4) Documentation contract
- `README.md` is derived from `project.yaml` + `docs/readme.template.md`.
- Changelog updates go through `scripts/dev/update_changelog.py` only.

### 5) Role + evidence contract
- AI must remain objective: no flattery, no fabricated conclusions.
- Non-trivial conclusions must include explicit evidence and uncertainty.
- Reference: `docs/avcp_guidelines.md#4.2 AI Role Positioning: Objectivity and Evidence`.

## 🚀 Scenario A: Start a New Project

### A1. Local bootstrap (human)

```bash
python -m pip install -e ".[dev]"
ruff check .
ruff format --check .
mypy .
pytest -q
python scripts/dev/generate_readme.py --check
```

### A2. First agent prompt (copy-paste)

```text
Read and enforce prompts/AVCP_SYSTEM_PROMPT_MIN.md.
Before coding, read docs/state.md, docs/constraints.md, docs/decisions.md,
docs/api_specs.md, docs/data_contracts.md, docs/avcp_guidelines.md.
Reply with:
1) [STATE SNAPSHOT]
2) [PLAN]
3) [PATCH SET]
4) [TEST]
5) [EVIDENCE]
No invented claims. If uncertain, escalate via gates.
```

### A3. Initialize project metadata

Ask the agent to:
1. update `project.yaml` (`name`, `title`, `domain`, `stage`, `owner`, `license`, `entrypoints`),
2. update sprint context in `docs/state.md`,
3. regenerate README via `python scripts/dev/generate_readme.py`.

### A4. Daily delivery loop

1. **Lock intent first**
- Tier-1/Tier-2 changes update docs/specs before code.

2. **Implement in small patches**
- Scripts under `scripts/` follow `docs/avcp_guidelines.md#4.1 Script Header Contract`.

3. **Verify before commit**
- Run lint/type/test/README checks.

4. **Record change safely**

```bash
python scripts/dev/update_changelog.py --entry "feat(scope): concise description"
```

5. **Commit with Conventional Commits**
- Example: `fix(bridge): enforce primary-key uniqueness on export`.

## 🛠 Scenario B: Migrate an Existing Project

### B0. Pick migration strategy

Choose one:
1. docs-first (recommended),
2. config-hardening,
3. full contract migration.

### B1. Inject AVCP skeleton

Copy into the existing repository root:
1. `prompts/`
2. `docs/`
3. `config/`
4. `scripts/dev/`

### B2. Migration prompt (copy-paste)

```text
This repository is migrating to AVCP.
Read prompts/AVCP_SYSTEM_PROMPT_MIN.md and docs/avcp_guidelines.md first.
Produce a docs-first migration snapshot:
- current state -> docs/state.md
- hard constraints -> docs/constraints.md
- open design decisions -> docs/decisions.md
No large refactor in this step.
```

### B3. Extract hardcoded runtime config (Tier-1)

Ask the agent to:
1. scan `src/` for hardcoded paths/parameters,
2. move them to `config/config.yaml`,
3. add fail-fast validation,
4. update `docs/api_specs.md` / `docs/data_contracts.md` if behavior/schema changes.

### B4. Standardize data output contracts

For cross-stage outputs:
1. adopt `save_for_r()` pattern or equivalent wrapper,
2. enforce primary key + schema/provenance sidecar,
3. document schemas in `docs/data_contracts.md`.

### B5. Turn on derived-doc workflow

```bash
python scripts/dev/generate_readme.py
python scripts/dev/generate_readme.py --check
```

Then enforce this in CI/pre-commit.

### B6. Migration done criteria

Migration is complete when:
1. startup prompt is repeatable,
2. decisions/constraints/specs are current in `docs/`,
3. runtime config is centralized in `config/config.yaml`,
4. output contracts are explicit and testable,
5. CI passes lint/type/test/README checks.

## 📋 Prompt Pack (Fast Reuse)

### Tier-0/1 implementation

```text
Consult docs/avcp_guidelines.md first.
Implement <feature> as Tier-1:
- update docs/specs for interface/schema changes,
- provide unified diff,
- run ruff + mypy + pytest,
- update changelog via scripts/dev/update_changelog.py.
```

### Tier-2 proposal (before coding)

```text
This is Tier-2.
Do not code yet.
Write assumptions, formal logic, pseudo-code, and validation plan in docs/decisions.md.
After human lock, implement with tests.
```

### Evidence-first reporting

```text
For each conclusion provide:
1) conclusion,
2) evidence list (files/outputs/metrics),
3) confidence,
4) unresolved uncertainty + next verification action.
Do not output unsupported conclusions.
```

## ⚠️ Common Pitfalls and Controls

| Pitfall | Control |
|---|---|
| Chat-only decisions | persist accepted decisions to `docs/decisions.md` |
| Silent fallback behavior | fail-fast + explicit warnings + manifest fields |
| Schema drift | update `docs/data_contracts.md` in same patch |
| README drift | edit `project.yaml` / `docs/readme.template.md` then regenerate |
| Script ambiguity | enforce `4.1 Script Header Contract` |
| Unsupported claims | enforce evidence-linked conclusions (`4.2`) |

## ✅ Recommended Commands

```bash
python -m pip install -e ".[dev]"
ruff check .
ruff format --check .
mypy .
pytest -q
python scripts/dev/generate_readme.py --check
python scripts/dev/update_changelog.py --entry "chore(docs): update guide"
```

## 📌 Project At A Glance

- **Name:** `avcp-template`
- **Domain:** Agentic research engineering
- **Stage:** Data / Pipeline / Analysis / Manuscript
- **Owner:** TODO-owner
- **License:** TODO-license

## ⚙️ Entrypoints

- `python -m pip install -e ".[dev]"`
- `pytest -q`
- `python scripts/dev/generate_readme.py --check`
## 📦 Outputs

- `data/interim_viz/*.parquet`
- `data/interim_viz/*.csv`
- `data/interim_viz/*_meta.json`
## 🔗 AVCP References

- Pinned system prompt: `prompts/AVCP_SYSTEM_PROMPT_MIN.md`
- Guidelines and gates: `docs/avcp_guidelines.md`
<!-- AVCP:README:END -->

## Notes

This repository follows AVCP repo-as-memory conventions with docs and config as source-of-truth artifacts.
