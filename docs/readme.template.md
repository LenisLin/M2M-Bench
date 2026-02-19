## {{ title }}

{{ one_liner }}

### Project At A Glance

- **Name:** `{{ name }}`
- **Domain:** {{ domain }}
- **Stage:** {{ stage }}
- **Owner:** {{ owner }}
- **License:** {{ license }}

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

{% for command in entrypoints %}- `{{ command }}`
{% endfor %}
{% if datasets %}
### Datasets

{% for dataset in datasets %}- `{{ dataset }}`
{% endfor %}
{% endif %}
{% if outputs %}
### Outputs

{% for output in outputs %}- `{{ output }}`
{% endfor %}
{% endif %}
### AVCP References

- Pinned system prompt: `prompts/AVCP_SYSTEM_PROMPT_MIN.md`
- Guidelines and gates: `docs/avcp_guidelines.md`
