"""Residue definitions from data/converting_dictionary.json: which force-field
residue (and terminal variant) a PDB residue maps to."""

from __future__ import annotations

import json
from functools import lru_cache
from typing import Dict, List, Optional, Tuple

from app.core.config import settings
from app.utils.adams4sims_processing_library.utils.alias import name_alias, resn_alias

@lru_cache(maxsize=1)
def load_converting_dictionary() -> dict:
    """
    NAČTE KONVERZNÍ SLOVNÍK ZE SOUBORU JSON V KOŘENOVÉM ADRESÁŘI A ZAJIŠŤUJE JEHO CACHOVÁNÍ PRO RYCHLÝ PŘÍSTUP.
    """
    dict_path = settings.BASE_DIR / "data" / "converting_dictionary.json"

    try:
        with open(dict_path, "r", encoding="utf-8") as f:
            return json.load(f)
    except Exception:
        return {}


def _parse_residues_from_pdb(pdb_text: str, chain: Optional[str]) -> List[Tuple[str, int, str, str, List[str]]]:
    """
    EXTRAHUJE DATA O REZIDUÍCH A JEJICH ATOMECH Z PDB FORMÁTU, PŘIČEMŽ PROVÁDÍ ALIASING NÁZVŮ ATOMŮ PODLE SLOVNÍKU.
    """
    residue_data = {}
    ordered_keys = []

    for line in pdb_text.splitlines():
        if not (line.startswith("ATOM") or line.startswith("HETATM")):
            continue

        resname = line[17:20].strip()
        ch = (line[21] or "").strip() or "?"
        resseq_raw = line[22:26].strip()
        icode = (line[26] or " ").strip()
        atom_name_raw = line[12:16].strip()

        atom_name = name_alias(resname, atom_name_raw)

        if not resseq_raw:
            continue

        try:
            resseq = int(resseq_raw)
        except ValueError:
            continue

        if chain and ch != chain:
            continue

        key = (ch, resseq, icode)
        if key not in residue_data:
            residue_data[key] = {"resname": resname, "atoms": []}
            ordered_keys.append(key)

        if atom_name not in residue_data[key]["atoms"]:
            residue_data[key]["atoms"].append(atom_name)

    # Řetězce abecedně, uvnitř řetězce ale pořadí ze souboru - to je pořadí
    # vazeb. Řazení podle (resseq, icode) by rozbilo insertion kódy, které
    # v PDB stojí PŘED reziduem se stejným číslem (chymotrypsinové číslování
    # u trypsinu: 183, 184A, 184, 185) - builder by pak vázal 183 -> 184 -> 184A.
    ordered_keys.sort(key=lambda x: x[0])

    return [
        (k[0], k[1], k[2], residue_data[k]["resname"], residue_data[k]["atoms"])
        for k in ordered_keys
    ]


def _infer_group(resname: str, conv: Dict) -> Optional[str]:
    """
    IDENTIFIKUJE CHEMICKOU KATEGORII REZIDUA PROHLEDÁVÁNÍM KLÍČŮ V KONVERZNÍM SLOVNÍKU NEBO POMOCÍ PREFIXŮ.
    """
    aliased = resn_alias(resname)
    for category in conv.keys():
        if isinstance(conv[category], dict):
            if resname in conv[category] or aliased in conv[category]:
                return category

    # Generické "HIS" je jediné standardní reziduum, jehož nevyřešený PDB
    # název NENÍ sám o sobě klíčem v converting_dictionary.json (jen jeho
    # HID/HIE/HIP tautomerní varianty jsou) - viz INTEGRATION_CONTRACT.md
    # invarianta #5. Bez týhle výjimky by výše uvedená smyčka pro "HIS"
    # nikdy nenašla kategorii, reziduum by nespadlo do main_chain a builder
    # by ho dostal jako nesouvisející HETATM (potvrzeno pádem na reálném
    # 1JJ2: "Passthrough atom HIS:N ... lacks converting identity") - i když
    # _pick_variant/_get_res_def níže samo o sobě HID/HIE/HIP podle vodíků
    # správně dohledá, tenhle chybějící "group" je to, co reziduum vyřazuje
    # z hlavního řetězce.
    if resname == "HIS" or aliased == "HIS":
        return "P"

    if resname.startswith("D"):
        return "D"
    if resname in {"A", "C", "G", "U"} or resname.startswith("R"):
        return "R"

    return None


def _get_res_def(group: Optional[str], ff_name: str, conv: Dict) -> Optional[Dict]:
    """
    VYHLEDÁ DEFINICI REZIDUA (ATOMY A KONEKTIVITU) V KONKRÉTNÍ KATEGORII NEBO PROHLEDÁNÍM CELÉHO SLOVNÍKU.
    """
    if group and group in conv and ff_name in conv[group]:
        return conv[group][ff_name]

    for category in conv.values():
        if isinstance(category, dict) and ff_name in category:
            return category[ff_name]
    return None


def _pick_variant(group: Optional[str], pdb_resname: str, atoms: List[str], conv: Dict, terminal: str) -> Tuple[
    Optional[str], bool, Optional[str]]:
    """
    URČUJE VHODNOU FF VARIANTU NA ZÁKLADĚ TERMINÁLNÍ POZICE A PŘÍTOMNÝCH ATOMŮ (NAPŘ. PRO HISTIDIN).
    """
    search_group = resn_alias(pdb_resname)

    # 1. INTELIGENTNÍ DETEKCE HISTIDINU
    # Pokud máme v PDB 'HIS' (nebo aliasovaný HID, HIE, HIP), zkusíme variantu potvrdit podle reálných vodíků
    if pdb_resname == "HIS" or search_group in ["HID", "HIE", "HIP"]:
        # HIP má oba vodíky (HD1 na delta-dusíku a HE2 na epsilon-dusíku)
        if "HD1" in atoms and "HE2" in atoms:
            search_group = "HIP"
        # HIE má vodík jen na epsilon dusíku
        elif "HE2" in atoms:
            search_group = "HIE"
        # HID má vodík na delta dusíku (nebo je to výchozí stav, pokud vodíky chybí)
        else:
            search_group = "HID"

    candidates: List[str] = []

    # 2. SESTAVENÍ KANDIDÁTŮ PRO KONCE ŘETĚZCŮ
    if group == "P":
        # Proteiny (skupina P) mají N-konec a C-konec (např. NALA, CALA, NHID, CHID)
        if terminal in ("5", "53"):
            candidates.append(f"N{search_group}")
        elif terminal == "3":
            candidates.append(f"C{search_group}")
    else:
        # Nukleové kyseliny mají koncovky 5 a 3 (např. RU5, RU3) a variantu N
        # pro reziduum, které je 5'- i 3'-koncem zároveň (osamocený nukleotid,
        # např. poslední reziduum za geometrickým zlomem). Samotná 5 nebo 3
        # varianta by tam nechala neceločíselný náboj (OL3: 5' = -0.3081,
        # 3' = -0.6919) a solvatace pak padá na "Fixed system charge ...
        # is not sufficiently close to an integer".
        if terminal == "53":
            candidates.append(f"{search_group}N")
        elif terminal == "5":
            candidates.append(f"{search_group}5")
        elif terminal == "3":
            candidates.append(f"{search_group}3")

    # Vždy přidáme jako fallback základní variantu uprostřed řetězce
    candidates.append(search_group)

    # 3. OVĚŘENÍ PROTI SLOVNÍKU
    for k in candidates:
        res_def = _get_res_def(group, k, conv)
        if res_def:
            return k, True, search_group

    return search_group, False, search_group
