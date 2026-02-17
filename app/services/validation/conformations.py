from typing import List, Dict, Any, Set
import numpy as np
import logging

logger = logging.getLogger(__name__)


class ConformationManager:
    @staticmethod
    def detect_alt_locs(pdb_content: str) -> List[Dict[str, Any]]:
        """
        Profesionální detekce AltLocs.
        Mapuje všechny varianty (A, B, C...) a ukládá souřadnice pro kontrolu kontinuity.
        """
        res_map = {}
        lines = pdb_content.splitlines()

        for line in lines:
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                alt_id = line[16].strip()
                if alt_id:
                    # Extrakce dat dle PDB standardu
                    res_name = line[17:20].strip()
                    chain = line[21].strip()
                    res_id = line[22:26].strip()
                    ins_code = line[26].strip()
                    atom_name = line[12:16].strip()

                    full_res_id = f"{res_id}{ins_code}"
                    key = f"{chain}-{full_res_id}"

                    if key not in res_map:
                        res_map[key] = {
                            "chain": chain,
                            "res_id": full_res_id,
                            "res_name": res_name,
                            "variants": set(),
                            "atom_coords": {}  # {variant: {atom_name: coords}}
                        }

                    res_map[key]["variants"].add(alt_id)

                    # Ukládáme souřadnice pro kontrolu kontinuity (páteř molekuly)
                    if atom_name in ["N", "C", "CA"]:
                        try:
                            coords = np.array([
                                float(line[30:38]),
                                float(line[38:46]),
                                float(line[46:54])
                            ])
                            if alt_id not in res_map[key]["atom_coords"]:
                                res_map[key]["atom_coords"][alt_id] = {}
                            res_map[key]["atom_coords"][alt_id][atom_name] = coords
                        except ValueError:
                            continue

        # Filtrace a příprava pro frontend
        detected = []
        sorted_keys = sorted(res_map.keys(), key=lambda k: (res_map[k]["chain"], res_map[k]["res_id"]))

        for key in sorted_keys:
            data = res_map[key]
            if len(data["variants"]) > 1:
                # Převedeme set na seřazený list, ale souřadnice necháme skryté v objektu
                data["variants"] = sorted(list(data["variants"]))
                # Pro frontend vyčistíme souřadnice, abychom neposílali MB dat,
                # ale v rámci třídy je můžeme použít pro validaci.
                frontend_item = {k: v for k, v in data.items() if k != "atom_coords"}
                detected.append(frontend_item)

        return detected

    @staticmethod
    def filter_pdb_by_selection(pdb_content: str, selections: Dict[str, str]) -> str:
        """
        Aplikuje výběr varianty a vyčistí AltLoc sloupec (pozice 17).
        """
        output = []
        for line in pdb_content.splitlines():
            if line.startswith(("ATOM", "HETATM")):
                alt_id = line[16].strip()
                if not alt_id:
                    output.append(line)
                    continue

                chain = line[21].strip()
                res_id = line[22:27].strip()
                key = f"{chain}-{res_id}"

                should_keep = False
                if key in selections:
                    if alt_id == selections[key]:
                        should_keep = True
                elif alt_id == 'A':  # Fallback pro nespecifikované
                    should_keep = True

                if should_keep:
                    # Vymazání AltLoc ID pro kompatibilitu se simulátory
                    output.append(line[:16] + " " + line[17:])
            else:
                output.append(line)

        return "\n".join(output)

    @staticmethod
    def validate_continuity(pdb_content: str, selections: Dict[str, str]) -> List[Dict[str, Any]]:
        struct_data = {}
        # Normalizace klíčů (odstranění mezer, aby "A-42-SER" sedělo)
        normalized_selections = {k.replace(" ", ""): v for k, v in selections.items()}

        for line in pdb_content.splitlines():
            if line.startswith(("ATOM", "HETATM")) and len(line) >= 54:
                atom_name = line[12:16].strip()
                if atom_name not in ["N", "C"]:
                    continue

                alt_id = line[16].strip()
                chain = line[21].strip()
                res_id_raw = line[22:27].strip()
                res_name = line[17:20].strip()  # <--- PŘIDÁNO: Načtení jména rezidua

                # KLÍČ MUSÍ ODPOVÍDAT FRONTENDU: "Chain-ID-Name"
                lookup_key = f"{chain}-{res_id_raw}-{res_name}".replace(" ", "")
                target_variant = normalized_selections.get(lookup_key)

                # DEBUG: Odkomentuj tohle, pokud chceš v logu vidět, co se děje
                # if alt_id:
                #     print(f"DEBUG: PDB AltLoc: {alt_id}, Key: {lookup_key}, Match: {target_variant}")

                keep = False
                if not alt_id:
                    keep = True
                elif target_variant and alt_id == target_variant:
                    keep = True
                elif not target_variant and alt_id == 'A':
                    keep = True

                if keep:
                    if chain not in struct_data:
                        struct_data[chain] = []
                    if not struct_data[chain] or struct_data[chain][-1]["res_id"] != res_id_raw:
                        struct_data[chain].append({"res_id": res_id_raw, "N": None, "C": None})

                    try:
                        coords = np.array([float(line[30:38]), float(line[38:46]), float(line[46:54])])
                        struct_data[chain][-1][atom_name] = coords
                    except:
                        continue

        warnings = []
        # 2. Výpočet vzdáleností mezi sousedními rezidui
        for chain, residues in struct_data.items():
            for i in range(len(residues) - 1):
                res_curr = residues[i]
                res_next = residues[i + 1]

                # Kontrolujeme vzdálenost mezi C atomem aktuálního a N atomem následujícího rezidua
                if res_curr["C"] is not None and res_next["N"] is not None:
                    dist = np.linalg.norm(res_curr["C"] - res_next["N"])

                    # Standardní peptidová vazba je ~1.33 A. Hranice 1.6 A spolehlivě detekuje chyby.
                    if dist > 1.6:
                        warnings.append({
                            "type": "CONTINUITY_GAP",
                            "message": f"Mezera {dist:.2f}Å detekována mezi rezidui {res_curr['res_id']} a {res_next['res_id']} v řetězci {chain}.",
                            "details": {
                                "chain": chain,
                                "res_i": res_curr['res_id'],
                                "res_j": res_next['res_id'],
                                "distance": round(float(dist), 3)
                            }
                        })
        return warnings