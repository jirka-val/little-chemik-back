# Local changes to the vendored FORGE builder

`app/builder/` is the FORGE builder v1.0.0 as delivered (see `VERSION`,
`INTEGRATION_CONTRACT.md` and `MANIFEST.sha256`). It is kept as close to the
delivered code as possible. Changes we had to make are listed here, so they
can be sent upstream or re-applied after an update.

## Layout differences from the delivered package

| Delivered | Here | Why |
|---|---|---|
| `src/*.py` | `*.py` (flat) | imported by bare module name; `app/services/structure/__init__.py` puts this directory on `sys.path` |
| `data/*.json` | `little-chemik-back/data/*.json` | shared with the rest of the backend; contents unchanged |
| `requirements.txt` | merged into `little-chemik-back/requirements.txt` | one dependency list for the image |
| `README.md`, `examples/` | not copied | not needed at runtime |

`MANIFEST.sha256` still lists the delivered paths and hashes. Files that
are unchanged match it (ignoring CRLF line endings):
`forge_molecule_builder.py`, `forge_molecule_mm.py`, `INTEGRATION_CONTRACT.md`,
`VERSION` and the four data files.

## Code changes

Each row: the file, the commit that changed it, and why. `git diff <base> HEAD -- app/builder/<file>`
shows the exact change (base = the commit where the file still matched the manifest).

| File | Base | Change | Commit |
|---|---|---|---|
| `forge_molecule_parser.py` | 5a1d78c | `Residue.terminus_reason`: why a residue became a terminus (gap, TER, geometry, chain end), so the PDB writer can put a TER record at a gap boundary | 73a5a19 |
| `forge_molecule_ions.py` | 56c5076 | renumber only solvent/ion chains after waters are replaced by ions; polymer chains keep their input numbering | d4125d5 |
| `forge_molecule_state_assignment.py` | 5a1d78c | user-forced protonation states (`forced`, `is_forced`): a forced residue keeps its state, its contacts still constrain the others | 593ba9e |
| `forge_workflow.py` | 56c5076 | `WorkflowSettings.protonation_overrides`, passed to state assignment as forced states | 593ba9e |
| `forge_molecule_solvation.py` | 56c5076 | drop template waters that clash across the periodic boundary | b678256 |

## Checking for unintended changes

```bash
cd app/builder
for f in forge_molecule_builder.py forge_molecule_mm.py; do
  grep "src/$f" MANIFEST.sha256 | tr -d '\r' | cut -d' ' -f1
  tr -d '\r' < "$f" | sha256sum
done
```

Any builder file that differs from the manifest and is not in the table
above is an unintended change.
