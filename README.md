# Little Chemik – backend

FastAPI service that analyses PDB structures and prepares simulation-ready
AMBER systems with the vendored FORGE builder.

```bash
python -m venv .venv
.venv/Scripts/pip install -r requirements-dev.txt     # Linux: .venv/bin/pip
cp .env.example .env
.venv/Scripts/python -m uvicorn app.main:app --reload --port 8000
```

API docs while running: http://localhost:8000/docs

```bash
python -m ruff check .     # lint
pytest                     # unit + integration
pytest -m golden           # whole pipeline vs. recorded baseline (~3 min)
bash scripts/check_all.sh  # everything, backend and frontend
```

Layout: `app/api` (HTTP), `app/services` (logic), `app/builder` (vendored
FORGE builder, see `app/builder/LOCAL_PATCHES.md`), `tests/` (see
`tests/README.md`), `scripts/` (`check_all.sh`, `incidents.py`).

Project documentation (architecture, workflow, API, deployment, operations)
is in `docs/` next to this repository.
