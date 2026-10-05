"""Molecule-type groups exactly as the FORGE builder distinguishes them
(top-level keys of data/converting_dictionary.json)."""

ION_GROUPS = frozenset({"I1", "I1+", "Im", "Im+"})
WATER_GROUPS = frozenset({"W3", "W4", "W5"})
POLYMER_GROUPS = frozenset({"R", "D", "P"})

# ff_selections from the frontend are keyed by the detected type ("W", "I")
# rather than by these - see forge_service._resolve_mol_type.
BUILDER_GROUPS = POLYMER_GROUPS | WATER_GROUPS | ION_GROUPS


def drop_unused_ion_groups(ff_selections: dict, used_groups) -> dict:
    """
    The force-field panel offers every ion group (I1/I1+/Im/Im+), so the
    selections can carry ion force fields the structure never uses. The
    builder and the topology get only the ion groups in `used_groups`;
    non-ion entries are kept as they are.
    """
    used = set(used_groups)
    return {key: ff for key, ff in ff_selections.items() if key not in ION_GROUPS or key in used}
