# avcp-template

Template repository for AVCP-based agentic development.

<!-- AVCP:README:START -->
## AVCP Template

Research engineering starter with AVCP prompts, repo-memory docs, and enforced bridge contracts.

This repository is an AVCP (Agentic Version Control Protocol) template for building software with coding agents (Codex / Claude Code / Cursor) under explicit engineering contracts.

AVCP solves three recurring failure modes in AI-assisted development:
1. Context amnesia: decisions disappear across sessions.
2. Hallucinated implementation: AI invents contracts, schemas, or logic.
3. Silent failures: pipelines continue with hidden errors.

The core mechanism is Repo-as-Memory: project truth is persisted in versioned files, not in chat history.

### AVCP Core Mechanisms

1. Agent cognition entrypoint (`prompts/` + `docs/`)
- `prompts/AVCP_SYSTEM_PROMPT_MIN.md` is the pinned operating protocol.
- `docs/state.md`, `docs/constraints.md`, `docs/decisions.md`, `docs/api_specs.md`, `docs/data_contracts.md` are long-term project memory.

2. Single Source of Truth (`config/`)
- Paths and runtime parameters must come from `config/config.yaml`.
- Hardcoded paths in scripts are treated as contract violations.

3. Contract-enforced data handoff (`src/avcp_template/io/bridge.py`)
- Use `save_for_r()` for Python->R handoff artifacts.
- Enforces output format + `<stem>_meta.json` sidecar + primary key invariants.

4. Documentation as controlled artifacts (`scripts/dev/`)
- `README.md` is derived from `project.yaml` + `docs/readme.template.md`.
- Changelog updates use `scripts/dev/update_changelog.py` (no blind append).

5. Risk-tiered execution gates
- Tier-0: low-risk refactor/tests/logging/perf.
- Tier-1: interface/schema/dependency changes require docs/spec updates first.
- Tier-2: new scientific/statistical logic must be specified and reviewed before coding.

6. AI role discipline: objective + evidence-backed
- The AI is an objective engineering collaborator, not a persuasive assistant.
- No flattery, no fabricated conclusions, and no invented certainty.
- Non-trivial conclusions must include explicit evidence list and confidence.

### Repository Layout

```text
prompts/                 # pinned system prompt
config/                  # runtime configuration source-of-truth
docs/                    # project memory, contracts, decisions
src/avcp_template/       # installable package
scripts/dev/             # changelog + README generators
tests/                   # unit/integration tests
```

### Scenario 1: Start a New Project from This Template

#### Step 1: Human bootstrap

```bash
python -m pip install -e ".[dev]"
ruff check .
ruff format --check .
mypy .
pytest -q
python scripts/dev/generate_readme.py --check
```

#### Step 2: Initialize agent cognition (first chat message)

Send this as your first prompt in Codex or Claude Code:

```text
Please read and strictly follow prompts/AVCP_SYSTEM_PROMPT_MIN.md.
Before coding, read docs/state.md, docs/constraints.md, docs/decisions.md,
docs/api_specs.md, docs/data_contracts.md, docs/avcp_guidelines.md.
Then reply with:
1) a [STATE SNAPSHOT],
2) a short [PLAN],
3) proposed [PATCH SET],
4) test plan.
Do not invent algorithms or contracts; escalate uncertain items via gates.
```

Append this role requirement in the same first prompt:

```text
Remain strictly objective and evidence-based:
- do not flatter the user,
- do not fabricate findings,
- for each conclusion, provide a numbered evidence list with concrete references.
If evidence is insufficient, say so explicitly and propose next verification steps.
```

#### Step 3: Project initialization tasks

Use the agent to fill metadata and regenerate README:

```text
Update project.yaml fields (name/title/domain/stage/owner/license/entrypoints),
refresh docs/state.md current sprint goal,
and run python scripts/dev/generate_readme.py.
```

#### Step 4: Build with daily AVCP loop

1. Lock intent first
- For Tier-1/Tier-2 tasks, update `docs/decisions.md` and contracts before code.

2. Implement with script contract
- Any script under `scripts/` must follow `docs/avcp_guidelines.md#4.1 Script Header Contract`.

3. Validate before commit
- Run lint/type/tests/README check.

4. Record change safely
- Update changelog via:

```bash
python scripts/dev/update_changelog.py --entry "feat(scope): concise description"
```

5. Commit with Conventional Commit message
- Example: `feat(bridge): enforce primary-key validation for export`.

### Scenario 2: Migrate an Existing Project into AVCP

#### Step 0: Define migration scope

Decide if migration is:
1. Docs-only first (recommended),
2. config hardening,
3. full code contract migration.

#### Step 1: Inject AVCP skeleton

Copy the following into your existing repository root:
1. `prompts/`
2. `docs/`
3. `config/`
4. `scripts/dev/`

#### Step 2: Load migration directive into agent

```text
This repository is migrating to AVCP.
Read prompts/AVCP_SYSTEM_PROMPT_MIN.md and docs/avcp_guidelines.md first.
Scan current structure and write a migration snapshot into docs/state.md,
constraints into docs/constraints.md,
and open decisions into docs/decisions.md.
No large refactor yet.
```

#### Step 3: Remove hardcoded runtime values (Tier-1)

```text
Scan src/ for hardcoded paths/parameters.
Move them to config/config.yaml.
Refactor code to load config values with fail-fast checks.
Update docs/api_specs.md or docs/data_contracts.md when behavior changes.
```

#### Step 4: Standardize data output contracts

For cross-stage outputs:
1. Adopt `save_for_r()` pattern or equivalent contract wrapper.
2. Enforce primary key, schema capture, provenance sidecar.
3. Document schema in `docs/data_contracts.md`.

#### Step 5: Turn on derived documentation workflow

```bash
python scripts/dev/generate_readme.py
python scripts/dev/generate_readme.py --check
```

Then enforce this in CI and pre-commit to prevent README drift.

#### Step 6: Cutover acceptance checklist

Migration is complete when:
1. Agent startup prompt is fixed and repeatable.
2. Decisions/constraints/specs live in `docs/` and are current.
3. Runtime config is centralized in `config/config.yaml`.
4. Output contracts are explicit and machine-checkable.
5. CI passes lint + type + tests + README check.

### Prompt Recipes (Copy-Paste)

#### A) Feature implementation (Tier-0/1)

```text
Consult docs/avcp_guidelines.md first.
Implement <feature> as Tier-1:
- update docs/specs if interface/schema changes,
- provide unified diff,
- run ruff, mypy, pytest,
- update changelog using scripts/dev/update_changelog.py.
```

#### B) Tier-2 research/statistics task

```text
This is Tier-2.
Do not code yet.
Write formal logic, assumptions, pseudo-code, and validation plan in docs/decisions.md.
After human lock, implement with tests.
```

#### C) Existing repo cleanup

```text
Perform a docs-first migration plan.
List hardcoded configs, implicit contracts, and silent-failure risks.
Propose minimal-risk patch batches with verification commands.
```

#### D) Evidence-first reporting

```text
For each conclusion, provide:
1) conclusion statement,
2) evidence list (file paths / outputs / metrics),
3) confidence level,
4) unresolved uncertainty and next verification action.
Do not state unsupported conclusions.
```

### Common Failure Modes and Controls

1. Chat-only decisions
- Control: persist every accepted decision in `docs/decisions.md`.

2. Silent fallback behavior
- Control: enforce fail-fast validation + explicit warnings + manifest fields.

3. Untracked schema drift
- Control: update `docs/data_contracts.md` in same patch as schema change.

4. Manual README edits causing drift
- Control: edit `project.yaml` or `docs/readme.template.md`, then regenerate.

5. Script behavior undocumented
- Control: enforce `4.1 Script Header Contract` for every executable script.

6. Unsupported claims in summaries
- Control: require evidence-linked conclusions using `docs/avcp_guidelines.md#4.2 AI Role Positioning: Objectivity and Evidence`.

### Recommended Dev Commands

```bash
python -m pip install -e ".[dev]"
ruff check .
ruff format --check .
mypy .
pytest -q
python scripts/dev/generate_readme.py --check
python scripts/dev/update_changelog.py --entry "chore(docs): update operational guide"
```

### Project At A Glance

- **Name:** `avcp-template`
- **Domain:** Agentic research engineering
- **Stage:** Data / Pipeline / Analysis / Manuscript
- **Owner:** TODO-owner
- **License:** TODO-license

### Entrypoints

- `python -m pip install -e ".[dev]"`
- `pytest -q`
- `python scripts/dev/generate_readme.py --check`
### Outputs

- `data/interim_viz/*.parquet`
- `data/interim_viz/*.csv`
- `data/interim_viz/*_meta.json`
### AVCP References

- Pinned system prompt: `prompts/AVCP_SYSTEM_PROMPT_MIN.md`
- Guidelines and gates: `docs/avcp_guidelines.md`
<!-- AVCP:README:END -->

## Notes

This repository follows AVCP repo-as-memory conventions with docs and config as source-of-truth artifacts.
