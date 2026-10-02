"""
Pins the public HTTP contract (paths, methods, parameters, request and
response schemas) to a committed OpenAPI snapshot, so a refactor cannot
rename, drop or reshape an endpoint the frontend relies on without noticing.

After an intended API change, re-record the snapshot:

    GOLDEN_UPDATE=1 pytest tests/integration/test_api_contract.py
"""

import json
import os
from pathlib import Path

import pytest

pytestmark = pytest.mark.integration

SNAPSHOT = Path(__file__).resolve().parents[1] / "fixtures" / "openapi.json"

# Free-text fields that do not affect the contract.
_DOC_KEYS = {"summary", "description", "title", "example", "examples"}


def _contract(node):
    if isinstance(node, dict):
        return {k: _contract(v) for k, v in sorted(node.items()) if k not in _DOC_KEYS}
    if isinstance(node, list):
        return [_contract(v) for v in node]
    return node


def test_openapi_matches_snapshot(client):
    schema = client.get("/openapi.json").json()
    actual = _contract({"paths": schema["paths"], "components": schema.get("components", {})})

    if os.environ.get("GOLDEN_UPDATE") == "1":
        SNAPSHOT.write_text(json.dumps(actual, indent=1) + "\n", encoding="utf-8")
        pytest.skip("OpenAPI snapshot recorded")

    expected = json.loads(SNAPSHOT.read_text(encoding="utf-8"))

    missing = sorted(set(expected["paths"]) - set(actual["paths"]))
    added = sorted(set(actual["paths"]) - set(expected["paths"]))
    assert not missing, f"endpoints removed: {missing}"
    assert not added, f"endpoints added (re-record the snapshot if intended): {added}"
    for path, item in expected["paths"].items():
        assert actual["paths"][path] == item, f"contract of {path} changed"
    assert actual["components"] == expected["components"], "request/response schemas changed"
