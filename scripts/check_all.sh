#!/usr/bin/env bash
# Runs every check that must pass before a change is committed or deployed:
#   1. backend unit + integration tests (incl. the OpenAPI contract snapshot)
#   2. golden pipeline tests (whole preparation pipeline vs. recorded baseline)
#   3. frontend type check and production build
#
# Usage (from anywhere):  bash little-chemik-back/scripts/check_all.sh [--skip-golden]
set -euo pipefail

BACK="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
FRONT="$BACK/../little-chemik-front"

if [ -x "$BACK/.venv/Scripts/python.exe" ]; then
    PY="$BACK/.venv/Scripts/python.exe"
elif [ -x "$BACK/.venv/bin/python" ]; then
    PY="$BACK/.venv/bin/python"
else
    PY=python
fi

cd "$BACK"
echo "== backend: unit + integration"
"$PY" -m pytest -q -p no:logging

if [ "${1:-}" != "--skip-golden" ]; then
    echo "== backend: golden pipeline"
    "$PY" -m pytest -q -p no:logging -m golden
fi

cd "$FRONT"
echo "== frontend: type check"
node_modules/.bin/tsc --noEmit
echo "== frontend: build"
node_modules/.bin/vite build --outDir "$(mktemp -d)" --emptyOutDir > /dev/null

echo "== all checks passed"
