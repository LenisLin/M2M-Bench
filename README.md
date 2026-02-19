# avcp-template

Template repository for AVCP-based agentic development:
- Minimal pinned system prompt: `prompts/AVCP_SYSTEM_PROMPT_MIN.md`
- Modular guidelines: `docs/avcp_guidelines.md`
- Repo-as-memory docs under `docs/`
- Code-enforced Python→R bridge: `src/io/bridge.py::save_for_r`
- Safe changelog update script: `scripts/dev/update_changelog.py`

## Quickstart
```bash
python -m pip install -e .[dev]
pytest -q
```
