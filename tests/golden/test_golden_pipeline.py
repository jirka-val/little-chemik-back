"""
Golden (characterization) tests of the whole preparation pipeline.

They pin what the app produces today for a set of reference structures, so a
refactor that silently changes any output fails here. They are not a spec: a
difference means "look at it", and an intended change is accepted by
re-recording the baseline.

    pytest -m golden                      # compare against the baseline
    GOLDEN_UPDATE=1 pytest -m golden      # re-record after an intended change

What is stored:
- expected/<case>/responses.json  sha256 of every JSON response along the way (committed)
- expected/<case>/files.json      sha256 + line/atom count of every output file (committed)
- .baseline/<case>/               full responses and output files, gzipped (git-ignored,
                                  local only - tens of MB, used to show what changed)

When a hash differs and the local full baseline exists, the files are
compared line by line with a small numeric tolerance, so floating-point noise
from another platform/library build is reported but does not fail, while a
real change shows the first differing lines.
"""

from __future__ import annotations

import gzip
import hashlib
import json
import math
import os
import re
import shutil
from pathlib import Path
from typing import Any, Dict, List

import pytest

from app.services.ff_catalog_service import catalog_service
from app.services.forcefield_service import ForceFieldService
from app.workspaces.manager import workspace_manager

from .pipeline import PipelineCase, run_pipeline

pytestmark = pytest.mark.golden

HERE = Path(__file__).resolve().parent
INPUTS = HERE.parent / "fixtures" / "pdb" / "golden"
EXPECTED = HERE / "expected"
BASELINE = HERE / ".baseline"
FROZEN_CATALOG = HERE / "ff_catalog_frozen.json.gz"
UPDATE = os.environ.get("GOLDEN_UPDATE") == "1"

# Relative tolerance for numeric tokens when the exact hash differs.
REL_TOL = 1e-6
ABS_TOL = 1e-4


def _case(name: str, pdb: str, **kwargs) -> PipelineCase:
    return PipelineCase(name=name, pdb_path=INPUTS / pdb, **kwargs)


CASES: List[PipelineCase] = [
    _case("1crn_protein", "1CRN.pdb"),
    _case("1crn_protein_vacuum", "1CRN.pdb", add_solvent=False),
    _case("1rna_rna", "1RNA.pdb"),
    _case("2tra_trna_altlocs_mg", "2TRA.pdb"),
    _case("3skr_chain_copies", "3SKR.pdb"),
    _case("4lzt_protein_altlocs", "4LZT.pdb"),
    _case("4lzt_keep_water_hmr", "4LZT.pdb", crystal_water_mode="keep_water", hmr=True),
    _case("2oue", "2OUE.pdb"),
    _case("3dvz", "3DVZ.pdb"),
    _case("1ca2_zinc", "HIS_1CA2_carbonic_anhydrase_zinc.pdb"),
    _case("1ubq_his", "HIS_1UBQ_ubiquitin_his68.pdb"),
    _case("1jj2_k_missing_sidechain", "1JJ2_chainK_78-92_gap83-89.pdb"),
    _case("1jj2_d_gap", "1JJ2_chainD_20-45_gap30-34.pdb"),
]


@pytest.fixture
def frozen_force_fields(monkeypatch, tmp_path):
    """Pin the FF catalog to the committed snapshot and keep FF file writes out of data/ff_cache*."""
    with gzip.open(FROZEN_CATALOG, "rt", encoding="utf-8") as fh:
        snapshot = json.load(fh)
    monkeypatch.setattr(catalog_service, "_snapshot_cache", snapshot)
    monkeypatch.setattr(catalog_service, "_buildable_ions_cache", None)
    monkeypatch.setattr(ForceFieldService, "CACHE_DIR", tmp_path / "ff_cache")
    monkeypatch.setattr(ForceFieldService, "FORGE_CACHE_DIR", tmp_path / "ff_cache_forge")


# ---------------------------------------------------------------------------
# Normalization
# ---------------------------------------------------------------------------

# The prmtop header carries the generation time.
_PRMTOP_DATE = re.compile(r"(%VERSION .*DATE = )\S+\s+\S+")


def _normalize_text(text: str, ws: str) -> str:
    text = text.replace("\r\n", "\n").replace(ws, "<ws>")
    return _PRMTOP_DATE.sub(r"\1<date>", text, count=1)


def _normalize_json(value: Any, ws: str) -> Any:
    return json.loads(_normalize_text(json.dumps(value, sort_keys=True), ws))


def _json_hash(value: Any) -> str:
    return hashlib.sha256(json.dumps(value, sort_keys=True).encode("utf-8")).hexdigest()


def _fingerprint(text: str) -> Dict[str, Any]:
    return {
        "sha256": hashlib.sha256(text.encode("utf-8")).hexdigest(),
        "lines": text.count("\n"),
        "atoms": sum(1 for line in text.splitlines() if line.startswith(("ATOM", "HETATM"))),
    }


# ---------------------------------------------------------------------------
# Tolerant comparison
# ---------------------------------------------------------------------------

_NUMBER = re.compile(r"[-+]?(?:\d+\.\d*|\.\d+|\d+)(?:[eE][-+]?\d+)?")


def _lines_match(a: str, b: str) -> bool:
    if a == b:
        return True
    na, nb = _NUMBER.findall(a), _NUMBER.findall(b)
    if len(na) != len(nb) or _NUMBER.sub("#", a) != _NUMBER.sub("#", b):
        return False
    return all(math.isclose(float(x), float(y), rel_tol=REL_TOL, abs_tol=ABS_TOL) for x, y in zip(na, nb))


def _diff_report(expected: str, actual: str, limit: int = 8) -> List[str]:
    exp_lines, act_lines = expected.splitlines(), actual.splitlines()
    problems = []
    if len(exp_lines) != len(act_lines):
        problems.append(f"line count {len(exp_lines)} -> {len(act_lines)}")
    for i, (a, b) in enumerate(zip(exp_lines, act_lines), start=1):
        if not _lines_match(a, b):
            problems.append(f"line {i}:\n  expected: {a}\n  actual:   {b}")
            if len(problems) >= limit:
                break
    return problems


def _first_difference(expected: Any, actual: Any, path: str) -> str | None:
    if type(expected) is not type(actual):
        return f"{path}: {expected!r:.200} -> {actual!r:.200}"
    if isinstance(expected, dict):
        for key in sorted(set(expected) | set(actual)):
            if key not in expected or key not in actual:
                return f"{path}.{key}: {'added' if key in actual else 'removed'}"
            found = _first_difference(expected[key], actual[key], f"{path}.{key}")
            if found:
                return found
        return None
    if isinstance(expected, list):
        if len(expected) != len(actual):
            return f"{path}: length {len(expected)} -> {len(actual)}"
        for i, (a, b) in enumerate(zip(expected, actual)):
            found = _first_difference(a, b, f"{path}[{i}]")
            if found:
                return found
        return None
    return None if expected == actual else f"{path}: {expected!r:.200} -> {actual!r:.200}"


# ---------------------------------------------------------------------------
# Test
# ---------------------------------------------------------------------------

@pytest.mark.parametrize("case", CASES, ids=[c.name for c in CASES])
def test_pipeline_matches_baseline(client, frozen_force_fields, case: PipelineCase):
    result = run_pipeline(client, case)
    ws = result.workspace_id  # type: ignore[attr-defined]
    try:
        responses = _normalize_json(result.responses, ws)
        files = {name: _normalize_text(text, ws) for name, text in result.files.items()}
        fingerprints = {name: _fingerprint(text) for name, text in sorted(files.items())}

        case_dir = EXPECTED / case.name
        baseline_dir = BASELINE / case.name

        if UPDATE:
            case_dir.mkdir(parents=True, exist_ok=True)
            response_hashes = {key: _json_hash(value) for key, value in sorted(responses.items())}
            (case_dir / "responses.json").write_text(json.dumps(response_hashes, indent=1) + "\n", encoding="utf-8")
            baseline_dir.mkdir(parents=True, exist_ok=True)
            with gzip.open(baseline_dir / "responses.json.gz", "wt", encoding="utf-8") as fh:
                json.dump(responses, fh, sort_keys=True)
            (case_dir / "files.json").write_text(json.dumps(fingerprints, indent=1, sort_keys=True) + "\n", encoding="utf-8")
            for name, text in files.items():
                target = baseline_dir / f"{name}.gz"
                target.parent.mkdir(parents=True, exist_ok=True)
                with gzip.open(target, "wt", encoding="utf-8", newline="\n") as fh:
                    fh.write(text)
            pytest.skip(f"baseline recorded for {case.name}")

        if not case_dir.exists():
            pytest.fail(f"No baseline for {case.name}; run with GOLDEN_UPDATE=1 first.")

        expected_response_hashes = json.loads((case_dir / "responses.json").read_text(encoding="utf-8"))
        full_responses_file = baseline_dir / "responses.json.gz"
        expected_responses = None
        if full_responses_file.exists():
            with gzip.open(full_responses_file, "rt", encoding="utf-8") as fh:
                expected_responses = json.load(fh)
        expected_files = json.loads((case_dir / "files.json").read_text(encoding="utf-8"))

        failures: List[str] = []

        for key in sorted(set(expected_response_hashes) | set(responses)):
            if key in responses and expected_response_hashes.get(key) == _json_hash(responses[key]):
                continue
            if expected_responses is None:
                failures.append(f"response '{key}' differs (no local full baseline to show where)")
            else:
                failures.append("response " + (_first_difference(expected_responses.get(key), responses.get(key), key)
                                               or f"{key}: hash differs"))

        if set(expected_files) != set(fingerprints):
            failures.append(f"output files differ: expected {sorted(expected_files)}, got {sorted(fingerprints)}")

        for name, fp in fingerprints.items():
            exp = expected_files.get(name)
            if exp is None or exp["sha256"] == fp["sha256"]:
                continue
            baseline_file = baseline_dir / f"{name}.gz"
            if not baseline_file.exists():
                failures.append(f"{name}: hash differs ({exp['lines']} -> {fp['lines']} lines, "
                                f"{exp['atoms']} -> {fp['atoms']} atoms); no local full baseline to diff")
                continue
            with gzip.open(baseline_file, "rt", encoding="utf-8") as fh:
                problems = _diff_report(fh.read(), files[name])
            if problems:
                failures.append(f"{name}:\n" + "\n".join(problems))
            else:
                print(f"[golden] {case.name}/{name}: differs only within numeric tolerance")

        assert not failures, f"{case.name} changed:\n" + "\n".join(failures)
    finally:
        shutil.rmtree(workspace_manager.get_workspace_dir(ws), ignore_errors=True)
