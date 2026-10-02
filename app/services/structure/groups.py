"""Molecule-type groups exactly as the FORGE builder distinguishes them
(top-level keys of data/converting_dictionary.json)."""

ION_GROUPS = frozenset({"I1", "I1+", "Im", "Im+"})
WATER_GROUPS = frozenset({"W3", "W4", "W5"})
POLYMER_GROUPS = frozenset({"R", "D", "P"})

# ff_selections from the frontend are keyed by the detected type ("W", "I")
# rather than by these - see forge_service._resolve_mol_type.
BUILDER_GROUPS = POLYMER_GROUPS | WATER_GROUPS | ION_GROUPS
