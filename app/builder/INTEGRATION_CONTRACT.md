# FORGE builder v1 integration contract

## Scope

This package is the production-facing v1 structure-treatment core. It parses
an upstream-cleaned FORGE structure, assigns covalent/protonation states,
builds all atoms reachable without an external degree of freedom, supports
interactive completion of safe residue-local open branches, and—only after the
accepted build plan is complete—solvates and adds ions.

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

Safe `residue_local_open_branch` components in this remaining plan can be
completed through the side-chain API below. Other missing DOFs remain explicit
boundaries and must not be silently bridged.

The benchmark-only operation that substitutes a template dihedral for missing
DOF is not part of this package and must not be reproduced in production.

## Interactive missing-side-chain API

The backend should prepare one immutable execution index and one set of FF
defaults after the ordinary executor stops at its first missing DOF:

```python
from forge_molecule_builder import (
    build_missing_sidechain_gui_payload,
    complete_missing_sidechains,
    optimize_missing_sidechain_dofs,
    optimize_missing_sidechain_preview,
    prepare_missing_sidechain_local_optimization,
    prepare_sidechain_execution_index,
    update_missing_sidechains,
)

index = prepare_sidechain_execution_index(result.remaining_plan)
mm_parameters = resources.force_field_parameters.merged_builder_provider()
default_dofs = optimize_missing_sidechain_dofs(
    result.molecule,
    resources.building_template,
    result.remaining_plan,
    mm_parameters,
    execution_index=index,
)
gui_payload = build_missing_sidechain_gui_payload(
    result.molecule,
    result.remaining_plan,
    default_dofs,
    execution_index=index,
)
working_molecule = complete_missing_sidechains(
    result.molecule,
    resources.building_template,
    result.remaining_plan,
    gui_payload,
    execution_index=index,
    mm_params=mm_parameters,
)
```

`gui_payload["sidechains"]` contains independent residue panels. Every DOF has
both immutable `default_degrees` and mutable `value_degrees`; therefore Reset
does not rerun MM optimization. A slider change calls
`update_missing_sidechains()` and returns only the transitively affected atom
coordinates.

The optional per-residue `Opt` button performs local pattern refinement from
the current slider values without repeating the 15-degree global scan:

```python
context = prepare_missing_sidechain_local_optimization(
    working_molecule,
    resources.building_template,
    result.remaining_plan,
    selected_residue,
    mm_parameters,
    execution_index=index,
)
working_molecule, response = optimize_missing_sidechain_preview(
    working_molecule,
    resources.building_template,
    result.remaining_plan,
    gui_payload,
    selected_residue,
    mm_parameters,
    execution_index=index,
    optimization_context=context,
)
```

Only the selected residue is optimized. Every other side-chain is fixed at its
current GUI geometry; no cluster identity is exposed to the client. `response`
is JSON-compatible and contains new DOF values, energy diagnostics and an
`updated_atoms` coordinate patch. The input payload and FF reset defaults are
not changed.

A prepared local context may be reused while only that selected residue
changes. It must be discarded after any coordinate outside the selected
residue changes. Omitting `optimization_context` is always correct and prepares
a fresh environment for the call.

The canonical build plan is intentionally not mutated by preview operations.
`SidechainExecutionIndex.sidechain_steps` and `.residual_steps` provide the
explicit partition needed by backend orchestration after the user accepts the
local completions.

## Core modules

- `forge_molecule_parser.py`: FORGE JSON/PDB records to `Molecule`.
- `forge_molecule_state_assignment.py`: covalent, protonation and tautomer
  state assignment.
- `forge_molecule_builder.py`: planner, dependency index and coordinate
  executor plus molecule/build-plan adapters for missing-side-chain work.
- `forge_molecule_mm.py`: shared UMFFF parameter loading, local switched
  LJ/Coulomb and torsion scoring, spatial indexing, free-rotor search and
  periodic-DOF optimizers.
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

The interactive side-chain API additionally returns JSON-compatible DOF
records and coordinate patches; it does not require retransmission of the
whole molecule for each slider or `Opt` action.

The core API does not write PDB. PDB is a lossy three-character-residue
visualization/export format; the bundled CLI contains an example writer only.

## Known v1 limits

- Missing polymer spans are not reconstructed.
- Safe single-anchor, acyclic residue-local open branches can be completed and
  optimized. Multi-anchor branches, cyclic missing fragments and missing
  polymer spans remain explicit unresolved DOFs.
- The local `Opt` action deliberately optimizes only one GUI residue against
  the current fixed environment; it does not jointly optimize neighboring
  side-chains.
- Rotamer-library and ML missing-DOF solvers remain future extensions.
- Noncanonical protonation families absent from the runtime templates do not
  exist for this workflow.
- The current solvation layer places W3 waters; W4/W5 extra sites remain a
  downstream force-field operation.
- Large systems remain more expensive, but builder, MM, solvation and ion
  placement use local topology/spatial indexes and avoid the earlier global
  quadratic update paths.

## Integration smoke tests

Before accepting a new release, integration should at minimum exercise:

- one RNA structure (1RNA or 3DVZ);
- one DNA structure (1BNA);
- one protein with state assignment (1IZ7 or 8TCA);
- one structural-ion example (2HO7 or 6N65);
- one expected missing-DOF example.

The development workspace retains the complete 105-test regression suite and
large-system profiling fixtures; they are intentionally not distributed in
this integration package.
