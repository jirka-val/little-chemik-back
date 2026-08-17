# FORGE builder v0 integration contract

## Scope

This package is the production-facing v0 structure-treatment core. It parses
an upstream-cleaned FORGE structure, assigns covalent/protonation states,
builds all atoms reachable without an external degree of freedom, and—only
after the build plan is complete—solvates and adds ions.

The package intentionally excludes development tests, benchmarks, PREP
template generators, historical snapshots, example structures and temporary
input-normalization helpers.

## Public Python API

The preferred entry point is:

```python
from forge_workflow import (
    WorkflowResources,
    WorkflowSettings,
    WorkflowResult,
    run_forge_workflow,
)
```

Load immutable resources once per worker process:

```python
resources = WorkflowResources.from_paths(
    converting_dictionary="data/converting_dictionary.json",
    building_template="data/building_template_v1.json",
    state_definitions="data/protonation_states_v1.json",
    solvation_template="data/solvation_template_v1.json",
    force_field_root="/path/to/UMFFF/forcefields",
)
```

Run one structure:

```python
result = run_forge_workflow(
    cleaned_structure_json,
    resources,
    salts=ordered_salt_specifications,
    settings=WorkflowSettings(),
)
```

`WorkflowResult.status` is one of:

- `complete`: all requested layers completed;
- `missing_dof`: building stopped safely before the first unresolved degree
  of freedom. No ion cleanup, coordinate recentering, solvation or ion
  placement has occurred.

When `status == "missing_dof"`, inspect `remaining_plan`. The first step is a
`PlannedMissingDOFStep` and explicitly identifies the requested dihedral. A
missing DOF is a normal modeled boundary, not a parser failure.

The benchmark-only operation that substitutes a template dihedral for missing
DOF is not part of this package and must not be reproduced in production.

## Core modules

- `forge_molecule_parser.py`: FORGE JSON/PDB records to `Molecule`.
- `forge_molecule_state_assignment.py`: covalent, protonation and tautomer
  state assignment.
- `forge_molecule_builder.py`: planner, dependency index and coordinate
  executor, including MM free-rotor H search.
- `forge_molecule_solvation.py`: periodic box construction and W3 water
  placement.
- `forge_molecule_ions.py`: crystal-ion cleanup, neutralization and salt
  excess.
- `forge_workflow.py`: supported orchestration boundary for integration.

Lower-level functions remain importable, but integration code should prefer
`run_forge_workflow()` so layer ordering and missing-DOF behavior stay
consistent.

## Required structure input

The input is a mapping with the same two top-level fields currently produced
by upstream FORGE:

```text
pdb_text
missing_atoms
```

`pdb_text` contains the selected coordinate model. `missing_atoms` contains
the chain/token analysis used to assign `mol_type`, `ff_resname`, expected
atoms, missing atoms, connectivity parts and gap/broken status.

The input object is treated as library data supplied by upstream. The builder
does not silently repair residue identity, aliases, AltLocs, chain gaps or
unknown chemistry.

## Mandatory upstream invariants

Before calling the builder, upstream must guarantee:

1. Exactly one coordinate model has been selected.
2. Alternate locations have been resolved to one unambiguous atom record.
3. Requested biological-assembly operations from REMARK 350 are complete.
4. Every supported residue and atom name has been normalized through the
   agreed alias library.
5. Generic `HIS` has been converted to the provisional/default supported
   identity expected by state assignment (currently HIE) while preserving its
   original chain position.
6. Waters and ions have assigned `mol_type` and `ff_resname` compatible with
   the converting dictionary.
7. Every non-passthrough residue is present in the supplied converting
   dictionary. Unsupported residues are removed upstream and returned to the
   user as structured warnings.
8. No residue contains duplicate coordinate records for the same canonical
   atom.
9. Chain and residue order represents polymer order, not merely arbitrary PDB
   record order.
10. Polymer discontinuities are explicit and terminality is internally
    consistent as described below.

Violation of these conditions is an upstream error. The production workflow
should report it rather than applying the development-only normalization shim.

## Required gap and terminality policy

This is the remaining upstream integration task requiring agreement.

If an unresolved polymer segment is not reconstructed, coordinate fragments
on its two sides must not be treated as directly bonded neighbors. Upstream
must detect such spans from available evidence, including:

- REMARK 465 missing residues;
- sequence/coordinate discontinuity;
- explicit connectivity information;
- chemically impossible inter-residue distance;
- explicit TER or equivalent structure metadata.

For the current conservative policy, both sides become artificial fragment
termini:

```text
RNA/DNA: previous observed residue -> 3' terminal variant
         next observed residue     -> 5' terminal variant

protein: previous observed residue -> C-terminal variant
         next observed residue     -> N-terminal variant
```

The corresponding tokens must carry an explicit gap/broken boundary so no
inter-residue build rule, torsion member, topology edge or MM exclusion can
cross the absent span.

Example from 1JJ2:

```text
U125, [C126-U127 missing], A128 -> RU3 ... RA5
GLU83, [VEDGG84-88 missing], PHE89 -> CGLU ... NPHE
```

If a future policy instead reconstructs a missing span or applies chemical
caps, that operation must run before this builder and must emit residue
identities/topology consistent with the converting dictionary.

## Runtime data contract

The package contains the exact validated versions of:

- `converting_dictionary.json`;
- `building_template_v1.json` covering R, D and P;
- `protonation_states_v1.json`;
- `solvation_template_v1.json` containing the complete inline W3/TIP3P box
  geometry used for solvent placement.

Force fields are not bundled. `force_field_root` must contain one compatible
UMFFF directory for every `mol_type` present or generated in the workflow.
Directory names end in the corresponding `mol_type`, and each contains the
residue library, nonbonded, bonded and force-field/default files expected by
`SolvationVdwParameters.from_force_field_root()`.

All force-field residue names, atom names and atom types must match the
converting dictionary and templates. Missing MM parameters are a hard error.

## Salt input

Ordered salt specifications have the form:

```json
{
  "salts": [
    {
      "cation": {"mol_type": "I1", "resname": "K+"},
      "anion": {"mol_type": "I1", "resname": "Cl-"},
      "concentration": 0.15
    }
  ]
}
```

`examples/ion_composition_example.json` is an example request, not a runtime
library. Production callers normally construct this object from user input.

The first salt supplies neutralizing counterions. Salt excess is then added in
input order. An empty list performs net-neutralization only using the default
K+/Cl- identities.

## Output contract

`WorkflowResult` returns Python objects, not formatted text:

- the current `Molecule`;
- the remaining `BuildPlan`;
- structured state-assignment report;
- when complete, structured crystal-ion, solvation and ion-addition reports.

Warnings and protonation conflicts remain structured in the reports and in
`Molecule.warnings`. The web/backend layer decides how to expose them to the
user.

The core API does not write PDB. PDB is a lossy three-character-residue
visualization/export format; the bundled CLI contains an example writer only.

## Known v0 limits

- Missing polymer spans are not reconstructed.
- A flexible missing side chain or base can legitimately require external
  DOF values.
- Interactive, rotamer-library, ML and environment-scored missing-DOF solvers
  are future extensions.
- Noncanonical protonation families absent from the runtime templates do not
  exist for this workflow.
- The current solvation layer places W3 waters; W4/W5 extra sites remain a
  downstream force-field operation.
- Large-system neutralization is the principal remaining performance hotspot.

## Integration smoke tests

Before accepting a new release, integration should at minimum exercise:

- one RNA structure (1RNA or 3DVZ);
- one DNA structure (1BNA);
- one protein with state assignment (1IZ7 or 8TCA);
- one structural-ion example (2HO7 or 6N65);
- one expected missing-DOF example.

The development workspace retains the complete 80-test regression suite and
large-system profiling fixtures; they are intentionally not distributed in
this integration package.
