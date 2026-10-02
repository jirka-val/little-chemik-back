# Tests – Little Chemik backend

## Layout

```
tests/
  conftest.py          shared fixtures (TestClient, workspace factory, PDB fixtures,
                       offline_forge_ff for testing the FORGE builder without network)
  fixtures/pdb/        small deterministic PDB files used by the tests
  fixtures/pdb/golden/ input structures of the golden tests
  fixtures/openapi.json  snapshot of the HTTP contract (test_api_contract.py)
  unit/                pure logic: no I/O, no network, no TestClient
  integration/         through the FastAPI TestClient, local workspace lifecycle
  golden/              whole pipeline vs. recorded baseline   (marker: golden)
  network/             needs the real RCSB/IDA APIs            (marker: network)
  performance/         long-running load tests                 (marker: slow)
```

**Where a test goes:** if it imports and calls a function or class directly
(no `client.post(...)`), it belongs in `unit/`. If it goes through the
`client` fixture (FastAPI `TestClient`), it belongs in `integration/`, unless
it needs the real network. Anything that downloads from RCSB or IDA belongs
in `network/`.

## Running

```bash
# default run: unit + integration, fast, no network (~1 min)
pytest

# one layer
pytest tests/unit -v
pytest tests/integration -v

# golden pipeline tests (~3 min)
pytest -m golden

# network tests (need internet)
pytest -m network

# performance tests
pytest -m slow

# everything
pytest -m "unit or integration or golden or network or slow"
```

Markers are registered in `pytest.ini` (`--strict-markers`, so a typo in a
marker fails immediately).

`scripts/check_all.sh` runs lint, the default run, the golden tests and the
frontend checks in one go.

## Golden tests

`golden/test_golden_pipeline.py` runs 13 reference structures through the
same requests the frontend makes in Guided mode and compares all responses
and output files (PDB, prmtop, crd, mdin, export ZIP, viewer download) with
a recorded baseline. The force-field catalog is frozen in
`golden/ff_catalog_frozen.json.gz`.

- `golden/expected/` (committed): sha256 of every response and file.
- `golden/.baseline/` (git-ignored): the full outputs, used to print the
  first differing lines; without it a failure only names the changed file.
- After an intended change: `GOLDEN_UPDATE=1 pytest -m golden`, then commit
  the new hashes with the change.

The same `GOLDEN_UPDATE=1` re-records the OpenAPI snapshot when run on
`integration/test_api_contract.py`.

## Key safety rule

**No test may write to `data/ff_cache/` or `data/ff_cache_forge/`.** These
hold real force-field files downloaded from the IDA API. Tests that need a
force field for the FORGE builder (`app/builder`) use the `offline_forge_ff`
fixture from `conftest.py`: it builds an isolated copy in `tmp_path` from
what is cached locally, and skips the test (`pytest.skip`) when that force
field is not cached, instead of failing or (worse) overwriting the real
cache with empty or wrong content. The golden tests redirect both cache
directories to `tmp_path` and take force-field files from the frozen catalog.

For the same reason `unit/test_forcefield_service.py` always monkeypatches
`ForceFieldService.CACHE_DIR`/`FORGE_CACHE_DIR` to `tmp_path` through the
`service` fixture.

## PDB fixtures (`fixtures/pdb/`)

| file | what it tests |
|---|---|
| `alanine_single.pdb` | smallest valid protein fragment |
| `protein_gap.pdb` | GLU83 → [84-88 missing] → PHE89, synthetic version of 1JJ2 chain K; terminals at a sequence gap |
| `rna_gap.pdb` | **real excerpt from 1RNA** (chain A, residues 1-3 and 6-8, 4-5 left out); needs real geometry, because the FORGE builder fails on degenerate/collinear atoms with `Cannot normalize near-zero central dihedral bond` |
| `altloc_sample.pdb` | SER42 with two alternative OG conformations (occupancy 0.6/0.4) |
| `two_model_nmr.pdb` | the same residue in 2 MODEL blocks (NMR ensemble) |

A new fixture PDB that goes through `ForgeStructureService` (anything in
`integration/test_prepare_*` or `network/`) needs real, non-degenerate
geometry: computed by hand or, more reliably, cut out of a real downloaded
structure like `rna_gap.pdb`.

## Solvation and the periodic box

`network/test_solvation_box.py` is the counterpart of the old
`test_performance.py::test_03_solvation_performance` (PDBFixer `addSolvent`),
now on top of `ForgeStructureService`. Learned while writing it: the builder
looks up the water force field strictly under `mol_type="W3"` (not the
generic `"W"` used by the older `TopologyService`/`pdb_service` path for the
AMBER topology), and even an empty salt list triggers neutralisation with the
default K+/Cl- ions (`mol_type="I1"`). For solvation `ff_selections` therefore
always needs an `I1` entry too, otherwise it fails with
`KeyError: Ion parameters missing for I1:K+`.

`TestSolvationCreatesBox` on 1RNA runs for real and finishes (OL3 + TIP3P +
JC-TIP3P-I1 are cached locally). `TestLargeStructureSolvationPerformance` on
1JJ2 mirrors the old test, but **1JJ2 never finishes solvation**: after ~93 s
it hits a legitimate `missing_dof` at `K:83 (CGLU:CG)` (GLU83), before
solvation starts. That is the expected, correct behaviour (1JJ2 really is an
incomplete structure), not a bug, so the test measures the time to that
point. `performance/test_forge_performance.py::TestForgeBuildPerformance`
confirms the same finding and timing.

## Known, documented gaps (not bugs in this code)

- `network/test_reference_structures.py::TestDnaReference::test_1bna_builds_successfully`
  is `xfail`: the locally cached `FF99BSC0` has no parameter for `DT:H72`.
  A gap in that force field's data, not in `ForgeStructureService`. Remove
  the `xfail` once the force field is fixed.
- `unit/test_pdb_topology_dict.py::TestWithoutSidecar::test_protein_defaults_to_mol_type_r`
  documents a pre-existing shortcoming outside the FORGE integration:
  without the `forge_meta` sidecar, `parse_pdb_to_topology_dict` gives
  proteins `mol_type="R"` instead of `"P"`.

## Adding tests

1. Pick the layer by the rule above.
2. If the test needs the FORGE builder (`ForgeStructureService`/`run_forge_workflow`),
   always use the `offline_forge_ff` fixture; never call `ForceFieldService`
   without a monkeypatch in a test that runs against the real `data/`.
3. If the test needs a PDB with real geometry, cut it out of a real
   downloaded structure (see `rna_gap.pdb`), or make sure your synthetic
   input does not go through the builder (pure token/text logic in
   `services/analysis` copes with synthetic coordinates).
4. Register any new marker in `pytest.ini` (`--strict-markers` rejects it otherwise).
