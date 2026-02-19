# avcp-template

Template repository for AVCP-based agentic development.

<!-- AVCP:README:START -->
## AVCP Template

Research engineering starter with AVCP prompts, repo-memory docs, and enforced bridge contracts.

### Project At A Glance

- **Name:** `avcp-template`
- **Domain:** Agentic research engineering
- **Stage:** Data / Pipeline / Analysis / Manuscript
- **Owner:** TODO-owner
- **License:** TODO-license

### Repository Layout

```text
prompts/       # AVCP pinned system prompt
config/        # runtime/project configuration (source of truth)
docs/          # repo-as-memory docs and contracts
src/avcp_template/  # installable Python package
scripts/dev/   # maintenance scripts (changelog/README generators)
tests/         # pytest suite
```

### Quickstart

```bash
python -m pip install -e ".[dev]"
ruff check .
mypy .
pytest -q
```

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
