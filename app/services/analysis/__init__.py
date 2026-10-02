"""
Structure analysis on raw PDB text (no workspace, no I/O):

- residues         force-field residue definitions (converting_dictionary.json)
- sequence         sequence tokens with gaps, breaks and missing atoms
- altlocs          models, symmetry and AltLoc recommendations
- chain_copies     identical copies of one molecule
- structure_prep   applying the user's choices from the Conformations step
- ff_requirements  which force-field groups and ions a structure needs
"""

from .altlocs import analyze_pdb_altlocs, find_altloc_selection_breaks
from .chain_copies import find_identical_chain_copies
from .ff_requirements import list_ion_options, required_ff_groups, resolve_ion_mol_type
from .residues import load_converting_dictionary
from .sequence import build_sequence_tokens
from .structure_prep import clean_pdb_altlocs, process_structure

__all__ = [
    "analyze_pdb_altlocs",
    "build_sequence_tokens",
    "clean_pdb_altlocs",
    "find_altloc_selection_breaks",
    "find_identical_chain_copies",
    "list_ion_options",
    "load_converting_dictionary",
    "process_structure",
    "required_ff_groups",
    "resolve_ion_mol_type",
]
