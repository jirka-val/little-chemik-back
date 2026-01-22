# import io
# import json
# import logging
# from typing import Dict, Any, List, Set
# from pdbfixer import PDBFixer
#
# logger = logging.getLogger(__name__)
#
#
# class ValidationService:
#     def __init__(self, dict_path: str = "converting_dictionary.json"):
#         try:
#             with open(dict_path, 'r') as f:
#                 self.ff_dict = json.load(f)  #
#
#             self.supported_residues = set()
#             for category in self.ff_dict.values():
#                 self.supported_residues.update(category.keys())  #
#             logger.info(f"Načteno {len(self.supported_residues)} podporovaných reziduí.")
#         except Exception as e:
#             logger.error(f"Nepodařilo se načíst konverzní slovník: {e}")
#             self.ff_dict = {}
#             self.supported_residues = set()
#
#     def validate_pdb_content(self, pdb_content: str, label: str = "current_state") -> Dict[str, Any]:
#         try:
#             # OPRAVA: Voláme správný název metody
#             alt_loc_info = self._detect_alt_loc_structured(pdb_content)
#
#             f = io.StringIO(pdb_content)
#             fixer = PDBFixer(pdbfile=f)
#
#             fixer.findMissingResidues()
#             fixer.findNonstandardResidues()
#             fixer.findMissingAtoms()
#
#             errors = []
#             warnings = []
#
#             # --- KONTROLA ALTERNATIVNÍCH KONFORMACÍ ---
#             if alt_loc_info:
#                 errors.append({
#                     "issue": "alt_locs_detected",
#                     "message": "Detekovány alternativní konformace.",
#                     "details": alt_loc_info,
#                     "action_required": "Vyberte jednu konformaci pro každé reziduum."
#                 })
#
#             # --- KONTROLA PODLE CONVERTING DICTIONARY (S MAPOVÁNÍM KONCŮ) ---
#             unsupported_res = set()
#             for chain in fixer.topology.chains():
#                 residues = list(chain.residues())
#                 for i, residue in enumerate(residues):
#                     res_name = residue.name
#
#                     # Logika pro mapování terminálních reziduí (DNA/RNA/Proteiny)
#                     possible_names = [res_name]
#
#                     # Pokud je to DNA (začíná na D) nebo RNA a je na kraji
#                     if len(residues) > 1:
#                         if i == 0:  # Začátek řetězce
#                             possible_names.append(f"{res_name}5")  # DNA/RNA konec
#                             possible_names.append(f"N{res_name}")  # Proteinový konec
#                         elif i == len(residues) - 1:  # Konec řetězce
#                             possible_names.append(f"{res_name}3")  # DNA/RNA konec
#                             possible_names.append(f"C{res_name}")  # Proteinový konec
#                     else:
#                         # Samostatný nukleotid
#                         possible_names.append(f"{res_name}N")
#
#                     # Zkontrolujeme, jestli alespoň jedno z možných jmen je ve slovníku
#                     if not any(name in self.supported_residues for name in possible_names):
#                         unsupported_res.add(res_name)
#
#             if unsupported_res:
#                 errors.append({
#                     "issue": "Nepodporovaná rezidua pro Force Field.",
#                     "details": sorted(list(unsupported_res)),
#                     "action_required": "Tato rezidua nejsou definována v converting_dictionary.json."
#                 })
#
#             # --- KONTROLA TĚŽKÝCH ATOMŮ ---
#             for residue, atoms in fixer.missingAtoms.items():
#                 atom_names = [atom.name for atom in atoms]
#                 errors.append({
#                     "resn": residue.name,
#                     "id": residue.index,
#                     "issue": f"Chybějící těžké atomy: {', '.join(atom_names)}"
#                 })
#
#             # --- KONTROLA DISKONTINUITY ---
#             for chain_res, res_names in fixer.missingResidues.items():
#                 chain_id = chain_res[0].id if hasattr(chain_res[0], 'id') else str(chain_res[0])
#                 errors.append({
#                     "chain": chain_id,
#                     "issue": f"V řetězci chybí úsek o délce {len(res_names)} reziduí."
#                 })
#
#             # --- KONTROLA VODÍKŮ ---
#             has_hydrogens = any(atom.element.symbol == 'H' for atom in fixer.topology.atoms())
#             if not has_hydrogens:
#                 errors.append({"issue": "Chybí vodíky. Je nutná protonace (Add Hydrogens)."})
#
#             return {
#                 "label": label,
#                 "is_ready_for_hpc": len(errors) == 0,
#                 "stats": {
#                     "total_errors": len(errors),
#                     "total_warnings": len(warnings),
#                     "atom_count": sum(1 for _ in fixer.topology.atoms())
#                 },
#                 "errors": errors,
#                 "warnings": warnings
#             }
#
#         except Exception as e:
#             logger.error(f"Chyba diagnostiky: {str(e)}", exc_info=True)
#             return {"error": f"Diagnostika selhala: {str(e)}"}
#
#     def _detect_alt_loc_structured(self, pdb_content: str) -> List[Dict[str, Any]]:
#         """Vytvoří mapu reziduí s jejich variantami pro frontend."""
#         res_map = {}
#         for line in pdb_content.splitlines():
#             if line.startswith(("ATOM", "HETATM")) and len(line) > 16:
#                 alt_id = line[16].strip()
#                 if alt_id:
#                     res_name = line[17:20].strip()
#                     res_id = line[22:26].strip()
#                     chain = line[21].strip()
#                     key = f"{chain}-{res_id}-{res_name}"
#
#                     if key not in res_map:
#                         res_map[key] = {"chain": chain, "res_id": res_id, "res_name": res_name, "variants": set()}
#                     res_map[key]["variants"].add(alt_id)
#
#         # Vrátíme jen ty, které mají skutečně více než jednu variantu (např. A i B)
#         return [v for v in res_map.values() if len(v["variants"]) > 1]
#
#     def apply_alt_loc_selection(self, pdb_content: str, selections: Dict[str, str]) -> str:
#         output = []
#         for line in pdb_content.splitlines():
#             if line.startswith(("ATOM", "HETATM")):
#                 alt_id = line[16].strip()
#                 if not alt_id:
#                     output.append(line)
#                     continue
#
#                 res_name = line[17:20].strip()
#                 res_id = line[22:26].strip()
#                 chain = line[21].strip()
#                 key = f"{chain}-{res_id}-{res_name}"
#
#                 # Pokud máme pro reziduum výběr, aplikujeme ho
#                 if key in selections:
#                     if alt_id == selections[key]:
#                         new_line = line[:16] + " " + line[17:]
#                         output.append(new_line)
#                 else:
#                     # Defaultně ponecháme variantu A a vymažeme příznak
#                     if alt_id == 'A':
#                         output.append(line[:16] + " " + line[17:])
#             else:
#                 output.append(line)
#         return "\n".join(output)