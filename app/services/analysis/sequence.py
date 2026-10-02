"""Sequence tokens for the sequence panel: residues, gaps, chain breaks,
terminals and missing/extra atoms per residue."""

from __future__ import annotations

from typing import Any, Dict, List, Optional, Set, Tuple

from app.utils.adams4sims_processing_library.utils.alias import name_alias

from .residues import _get_res_def, _infer_group, _parse_residues_from_pdb, _pick_variant, load_converting_dictionary

_BREAK_REASON_LABEL = {
    "gap": "sequence gap",
    "ter": "explicit chain break (TER record)",
    "geometry": "chemically implausible inter-residue distance",
}

# Boundary atoms that must be within bonding distance of the corresponding
# atom on the neighbour for the two residues to plausibly be covalently
# linked. Protein: previous C -> next N. Nucleic: previous O3' -> next P.
_BOUNDARY_ATOM_NAMES = {"P": ("C", "N"), "R": ("O3'", "P"), "D": ("O3'", "P")}
_BOND_DISTANCE_LIMIT_ANGSTROM = {"P": 1.9, "R": 2.1, "D": 2.1}

# Heavy (non-hydrogen) atoms a terminal variant is *expected* to still be
# missing purely because of the capping itself (e.g. the extra carboxylate
# oxygen on a protein C-terminus) - not a sign of an unbuildable gap.
_TERMINAL_EXTRA_HEAVY_ATOMS = {"P": {"OXT"}, "R": set(), "D": set()}

# Symmetric case, but on the *extra* side: heavy atoms a 5'-terminal nucleic
# variant is allowed to still HAVE even though its own FF template doesn't
# require them. A residue right after a sequence gap (or one that is simply
# a genuine 5'-phosphorylated terminus in the source structure) keeps the
# phosphate group that would normally link back to whatever precedes it -
# that's expected, not a sign of an unrecognized/mismatched residue.
_TERMINAL_ALLOWED_EXTRA_HEAVY_ATOMS = {"P": set(), "R": {"P", "OP1", "OP2"}, "D": {"P", "OP1", "OP2"}}

# The single backbone atom a residue right before a gap must still have for
# the builder's interactive side-chain completion (app/builder - see
# INTEGRATION_CONTRACT.md "residue_local_open_branch") to have any anchor to
# build from at all - same reference atoms as _BOUNDARY_ATOM_NAMES' "prev"
# side. Missing *this* atom means the residue genuinely has no usable
# connection point and must still be excluded (see _reterminate_as_gap_end).
# Missing anything else (e.g. a protein side chain past CB, or a nucleic base)
# is exactly what the interactive builder can now resolve, so it must no
# longer trigger exclusion by itself.
_GAP_BOUNDARY_ANCHOR_ATOM = {"P": "C", "R": "O3'", "D": "O3'"}


def _parse_remark465(pdb_text: str) -> Dict[Tuple[str, int, str], str]:
    """
    Autoritativní seznam reziduí, která nebyla v experimentu lokalizována,
    přímo z hlavičky PDB (REMARK 465) - viz INTEGRATION_CONTRACT.md, kde je
    tohle první z vyjmenovaných důkazů pro detekci polymerní mezery. Používá
    se k obohacení gap warningů o skutečné identity chybějících reziduí,
    místo pouhého odvození z díry v číslování.
    """
    missing: Dict[Tuple[str, int, str], str] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("REMARK 465"):
            continue
        tokens = line[10:].split()
        if len(tokens) < 3:
            continue
        resname, chain, seq_tok = tokens[-3], tokens[-2], tokens[-1]
        if len(chain) != 1:
            continue
        icode = ""
        if seq_tok and seq_tok[-1].isalpha():
            icode = seq_tok[-1]
            seq_tok = seq_tok[:-1]
        if not seq_tok.lstrip("-").isdigit():
            continue
        missing[(chain, int(seq_tok), icode)] = resname
    return missing


def _parse_ter_chain_breaks(pdb_text: str) -> Set[Tuple[str, int, str]]:
    """
    Vrátí (chain, resseq, icode) posledního rezidua PŘED každým TER záznamem,
    který není posledním výskytem daného řetězce v souboru - tedy řetězec
    pokračuje dalšími ATOM/HETATM záznamy i po tomto TER, což signalizuje
    fyzický zlom polymeru uprostřed jednoho PDB chain ID (viz
    INTEGRATION_CONTRACT.md - "explicit TER or equivalent structure
    metadata"). Neparsujeme vlastní sloupce TER záznamu (bývají nespolehlivé
    u hetero-ukončených řetězců) - řetězec a reziduum, které TER uzavírá,
    odvozujeme z posledního předchozího ATOM/HETATM záznamu.
    """
    lines = pdb_text.splitlines()
    last_atom: Optional[Tuple[str, int, str]] = None
    ter_events: List[Tuple[str, int, str, int]] = []

    for idx, line in enumerate(lines):
        if line.startswith("ATOM") or line.startswith("HETATM"):
            ch = (line[21] or "").strip() or "?"
            resseq_raw = line[22:26].strip()
            icode = (line[26] or " ").strip()
            if not resseq_raw:
                continue
            try:
                resseq = int(resseq_raw)
            except ValueError:
                continue
            last_atom = (ch, resseq, icode)
        elif line.startswith("TER") and last_atom is not None:
            ter_events.append((*last_atom, idx))

    breaks: Set[Tuple[str, int, str]] = set()
    for ch, resseq, icode, ter_idx in ter_events:
        for line in lines[ter_idx + 1:]:
            if (line.startswith("ATOM") or line.startswith("HETATM")) and (line[21] or "").strip() == ch:
                breaks.add((ch, resseq, icode))
                break
    return breaks


def _parse_boundary_atom_coords(pdb_text: str) -> Dict[Tuple[str, int, str, str], Tuple[float, float, float]]:
    """
    Souřadnice jen pro atomy, které tvoří kostru meziresiduové vazby (protein
    C/N, nukleové kyseliny O3'/P) - použito výhradně pro kontrolu chemicky
    nemožné meziresiduové vzdálenosti (INTEGRATION_CONTRACT.md). Netáhneme si
    sem souřadnice všech atomů, ať je to levné i na velkých strukturách.
    """
    wanted_names = {"C", "N", "O3'", "P"}
    coords: Dict[Tuple[str, int, str, str], Tuple[float, float, float]] = {}
    for line in pdb_text.splitlines():
        if not line.startswith("ATOM"):
            continue
        resname = line[17:20].strip()
        atom_name = name_alias(resname, line[12:16].strip())
        if atom_name not in wanted_names:
            continue
        ch = (line[21] or "").strip() or "?"
        resseq_raw = line[22:26].strip()
        icode = (line[26] or " ").strip()
        if not resseq_raw:
            continue
        try:
            resseq = int(resseq_raw)
            x, y, z = float(line[30:38]), float(line[38:46]), float(line[46:54])
        except ValueError:
            continue
        coords[(ch, resseq, icode, atom_name)] = (x, y, z)
    return coords


def _is_chemically_impossible_bond(
    group: Optional[str],
    prev_key: Tuple[str, int, str],
    curr_key: Tuple[str, int, str],
    coords: Dict[Tuple[str, int, str, str], Tuple[float, float, float]],
) -> bool:
    if group not in _BOUNDARY_ATOM_NAMES:
        return False
    prev_atom, curr_atom = _BOUNDARY_ATOM_NAMES[group]
    p = coords.get(prev_key + (prev_atom,))
    c = coords.get(curr_key + (curr_atom,))
    if p is None or c is None:
        return False
    dist = ((p[0] - c[0]) ** 2 + (p[1] - c[1]) ** 2 + (p[2] - c[2]) ** 2) ** 0.5
    return dist > _BOND_DISTANCE_LIMIT_ANGSTROM[group]


def _is_bonded(
    group: Optional[str],
    prev_key: Tuple[str, int, str],
    curr_key: Tuple[str, int, str],
    coords: Dict[Tuple[str, int, str, str], Tuple[float, float, float]],
) -> bool:
    """Opak _is_chemically_impossible_bond: obě kotvy známé a ve vazebné vzdálenosti."""
    if group not in _BOUNDARY_ATOM_NAMES:
        return False
    prev_atom, curr_atom = _BOUNDARY_ATOM_NAMES[group]
    p = coords.get(prev_key + (prev_atom,))
    c = coords.get(curr_key + (curr_atom,))
    if p is None or c is None:
        return False
    dist = ((p[0] - c[0]) ** 2 + (p[1] - c[1]) ** 2 + (p[2] - c[2]) ** 2) ** 0.5
    return dist <= _BOND_DISTANCE_LIMIT_ANGSTROM[group]


def _reterminate_as_gap_end(
    token: Dict[str, Any],
    conv: Dict,
    warnings: List[str],
    next_chain: str,
    next_resseq: int,
    break_reason: str = "gap",
    missing_residue_labels: Optional[List[str]] = None,
) -> None:
    """
    Přepíše už zapsaný token (poslední residuum PŘED přerušením řetězce v
    hlavním řetězci) na jeho umělou terminální variantu (C-konec pro protein /
    3'-konec pro RNA-DNA).

    Bez tohoto kroku by ff_resname zůstal ve "středové" variantě, která v konverzním
    slovníku očekává navazující sousední residuum - a downstream builder (app/builder)
    by pak mohl přes chybějící úsek vytvořit nesmyslnou vazbu. Builder sám gap
    nedostavuje (viz INTEGRATION_CONTRACT.md), takže tohle rozhodnutí musí padnout tady.
    """
    group = token["group"]
    # Reziduum, které už je 5'-koncem (první v řetězci nebo hned za jiným
    # zlomem), je teď osamocené - u nukleových kyselin potřebuje variantu N,
    # jinak by přepis na 3' nechal neceločíselný náboj (viz _pick_variant).
    already_5prime = token.get("terminal") in ("5", "53")
    terminal = "53" if already_5prime and group != "P" else "3"
    new_ff_resname, new_known, _ = _pick_variant(group, token["pdb_resname"], token["atoms"], conv, terminal)
    token["ff_resname"] = new_ff_resname
    token["known"] = new_known
    token["terminal"] = terminal
    token["missing_atoms"] = _check_missing_atoms(group, new_ff_resname, token["atoms"], conv)
    extra_atoms = _check_extra_atoms(group, new_ff_resname, token["atoms"], conv)
    if terminal == "53":
        allowed_extra_heavy = _TERMINAL_ALLOWED_EXTRA_HEAVY_ATOMS.get(group, set())
        extra_atoms = [a for a in extra_atoms if a not in allowed_extra_heavy]
    token["extra_atoms"] = extra_atoms
    conn_info = _check_connectivity_integrity(group, new_ff_resname, token["atoms"], conv)
    token["is_broken"] = conn_info["is_broken"]
    token["connectivity_parts"] = conn_info["components"]
    token["terminus_reason"] = break_reason

    label = "C-terminus" if group == "P" else "3'-terminus"
    reason_text = _BREAK_REASON_LABEL.get(break_reason, break_reason)
    detail = f" ({', '.join(missing_residue_labels)} missing per REMARK 465)" if missing_residue_labels else ""
    warnings.append(
        f"{token['chain']}:{token['resseq']} ({token['pdb_resname']}) treated as artificial {label} "
        f"— {reason_text} before {next_chain}:{next_resseq}{detail}."
    )

    # Terminal capping only ever ADDS hydrogens (extra NH3+/OH/carboxylate H)
    # or, for a protein C-terminus, the single OXT heavy atom - it never
    # requires rebuilding a side chain or base from scratch. A residue that's
    # still missing OTHER heavy atoms here was already incomplete in the
    # source structure - but that is no longer automatically unbuildable: the
    # builder's interactive side-chain completion (see
    # app/builder/INTEGRATION_CONTRACT.md, "residue_local_open_branch") can
    # safely resolve a single missing side-chain/base branch through the GUI,
    # as long as the residue still has its own backbone connection point
    # (_GAP_BOUNDARY_ANCHOR_ATOM - same C/O3' reference atom used above for
    # the boundary bond-distance check). Only exclude the residue when even
    # that anchor is gone - there is then genuinely nothing for the builder to
    # attach anything to, regardless of GUI support.
    allowed_extra = _TERMINAL_EXTRA_HEAVY_ATOMS.get(group, set())
    heavy_missing = [a for a in token["missing_atoms"] if a[:1] != "H" and a not in allowed_extra]
    anchor_atom = _GAP_BOUNDARY_ANCHOR_ATOM.get(group)
    token["gap_boundary_incomplete"] = bool(anchor_atom) and anchor_atom in heavy_missing
    if token["gap_boundary_incomplete"]:
        warnings.append(
            f"{token['chain']}:{token['resseq']} ({token['pdb_resname']}) is missing its own backbone "
            f"anchor atom ({anchor_atom}) even as an artificial {label} — the builder has nothing to "
            f"attach to here. If that happens, consider excluding this residue from the model and "
            f"shifting the terminus to the previous main-chain residue instead."
        )
    elif heavy_missing:
        warnings.append(
            f"{token['chain']}:{token['resseq']} ({token['pdb_resname']}) is still missing heavy atoms "
            f"{heavy_missing} even as an artificial {label} — the backbone anchor is present, so the "
            f"builder will offer this as an interactive side-chain completion instead of failing outright."
        )


def build_sequence_tokens(pdb_text: str, chain: Optional[str] = None, fill_gaps: bool = True):
    """
    SESTAVUJE KOMPLETNÍ SEZNAM TOKENS PRO DANÝ ŘETĚZEC, PROVÁDÍ ANALÝZU VARIANT A DETEKCI CHYBĚJÍCÍCH ČÁSTÍ STRUKTURY.
    """
    conv = load_converting_dictionary()
    all_residues = _parse_residues_from_pdb(pdb_text, chain)

    if not all_residues:
        return {"chain": chain, "tokens": [], "warnings": ["No residues found in PDB."]}

    residues_to_process = [r for r in all_residues if r[0] == chain] if chain else all_residues

    chains_dict = {}
    for r in residues_to_process:
        ch_id = r[0]
        if ch_id not in chains_dict:
            chains_dict[ch_id] = []
        chains_dict[ch_id].append(r)

    # Doplňkové důkazy o přerušení polymeru, které nejsou vidět jen z díry v
    # číslování reziduí - viz INTEGRATION_CONTRACT.md "Required gap and
    # terminality policy". Parsují se jednou za celý soubor, ne per-chain.
    remark465 = _parse_remark465(pdb_text) if fill_gaps else {}
    ter_breaks = _parse_ter_chain_breaks(pdb_text) if fill_gaps else set()
    boundary_coords = _parse_boundary_atom_coords(pdb_text) if fill_gaps else {}

    tokens = []
    warnings = []
    global_pos = 0

    for ch_id, residues in chains_dict.items():
        main_chain = []
        ligands = []

        for r in residues:
            resname = r[3]
            group = _infer_group(resname, conv)
            if group in ["R", "D", "P"]:
                main_chain.append(r)
            else:
                ligands.append(r)

        # Pořadí ze souboru (viz _parse_residues_from_pdb), žádné řazení podle čísel.
        first_main_key = (main_chain[0][1], main_chain[0][2]) if main_chain else None
        last_main_key = (main_chain[-1][1], main_chain[-1][2]) if main_chain else None
        processed_ordered = main_chain + ligands
        prev_resseq = None
        prev_icode = ""
        prev_main_token = None
        # Historie doposud přidaných hlavních (main-chain) tokenů tohoto
        # řetězce - umožňuje při kaskádovém vylučování neúplných okrajových
        # reziduí (viz níže) sáhnout i za bezprostředně předchozí reziduum.
        main_chain_history: List[Dict[str, Any]] = []

        for ch, resseq, icode, resname, atoms in processed_ordered:
            is_main = any(r[1] == resseq and r[3] == resname for r in main_chain)
            after_gap = False
            break_reason = None
            group = _infer_group(resname, conv)

            if fill_gaps and is_main and prev_resseq is not None:
                gap_found = False
                gap_labels: List[str] = []

                prev_key = (ch, prev_resseq, prev_icode)
                curr_key = (ch, resseq, icode or "")
                prev_group = prev_main_token["group"] if prev_main_token else group
                # Díra v číslování není důkaz chybějícího úseku, když jsou
                # sousedé prokazatelně vázaní (C-N / O3'-P v délce vazby).
                # Chymotrypsinové číslování trypsinu přeskakuje čísla
                # (34 -> 37, 217 -> 219) u souvislého řetězce - dřív se tam
                # řetězec uměle rozřízl na nabité konce.
                numbering_only_skip = (
                    resseq > prev_resseq + 1
                    and _is_bonded(prev_group, prev_key, curr_key, boundary_coords)
                )

                if resseq > prev_resseq + 1 and not numbering_only_skip:
                    # OPRAVA: Zkontroluj, zda GAP je OPRAVDU prázdný nebo tam jen je neznámé reziduum
                    for missing_seq in range(prev_resseq + 1, resseq):
                        # PDB číslování reziduí nikdy nepoužívá sekvenční číslo 0 (konvence
                        # při přechodu ze záporného na kladné číslování, např. -1 -> 1 u
                        # konstruktů s uměle přidaným 5'/N-koncovým leaderem - viz 2OUE
                        # chain A). "Chybějící" reziduum 0 proto není skutečná mezera a
                        # nesmí vytvořit GAP placeholder token.
                        if missing_seq == 0:
                            continue

                        # Hledej, zda existuje JAKÉKOLI reziduum se sekvencí missing_seq v PDB
                        residue_exists_in_pdb = any(r[1] == missing_seq and r[0] == ch for r in residues)

                        # Jen pokud OPRAVDU chybí v PDB -> vytvoř GAP token
                        if not residue_exists_in_pdb:
                            gap_found = True
                            missing_resname = remark465.get((ch, missing_seq, ""), "?")
                            gap_labels.append(f"{missing_resname}{missing_seq}")
                            global_pos += 1
                            # resseq/pdb_resname nesou SKUTEČNOU identitu chybějícího rezidua
                            # (missing_seq je vždy known - je to přímo číslo, na kterém
                            # smyčka zrovna je; pdb_resname z REMARK 465, pokud ho PDB
                            # hlavička uvádí, jinak zůstává "?"). Dřív se sem tvrdě
                            # zapisovalo resseq=None/pdb_resname="0" pro KAŽDÉ chybějící
                            # reziduum bez rozdílu - frontend to pak nedokázal ukázat jinak
                            # než jako nerozlišitelnou řadu "UNK ??" řádků u delších mezer.
                            tokens.append({
                                "position": global_pos, "chain": ch, "resseq": missing_seq, "icode": None,
                                "pdb_resname": missing_resname, "is_gap": True, "group": None, "ff_resname": None,
                                "known": False, "atoms": [], "missing_atoms": [], "extra_atoms": []
                            })
                    break_reason = "gap"
                else:
                    # Číslování je souvislé, ale řetězec může být přesto fyzicky
                    # přerušený - explicitní TER uprostřed řetězce nebo chemicky
                    # nemožná meziresiduová vzdálenost (další dva důkazy jmenované
                    # v INTEGRATION_CONTRACT.md vedle díry v číslování).
                    if prev_key in ter_breaks:
                        gap_found = True
                        break_reason = "ter"
                    elif _is_chemically_impossible_bond(prev_group, prev_key, curr_key, boundary_coords):
                        gap_found = True
                        break_reason = "geometry"

                if gap_found:
                    # Residuum před přerušením i residuum za ním se stávají uměle
                    # terminálními, ať se přes chybějící/přerušený úsek nepočítá
                    # žádná vazba.
                    after_gap = True
                    excluded_here: List[Dict[str, Any]] = []
                    # Pokud reziduum bezprostředně před přerušením zůstane
                    # neúplné (chybí těžké atomy) i po přeznačení na
                    # terminální variantu, builder ho stejně nikdy nedostaví
                    # (nemá kotvu pro vnitřní dihedral) - ponechat ho v
                    # modelu by z čistého přerušení udělalo neřešitelný
                    # missing_dof pád. Takové reziduum se z modelu vyřadí a
                    # terminalita se zkusí o krok blíž k začátku řetězce -
                    # opakovaně, dokud nenarazíme na použitelné reziduum.
                    while main_chain_history:
                        candidate = main_chain_history[-1]
                        _reterminate_as_gap_end(
                            candidate, conv, warnings, ch, resseq,
                            break_reason=break_reason, missing_residue_labels=gap_labels or None,
                        )
                        if not candidate["gap_boundary_incomplete"]:
                            break
                        excluded_here.append(candidate)
                        for idx, existing in enumerate(tokens):
                            if existing is candidate:
                                del tokens[idx]
                                break
                        main_chain_history.pop()

                    if excluded_here:
                        names = ", ".join(
                            f"{c['chain']}:{c['resseq']} ({c['pdb_resname']})" for c in excluded_here
                        )
                        warnings.append(
                            f"Excluded {names} from the model — still missing heavy atoms even as an "
                            f"artificial terminus, so the terminus was shifted further back."
                        )
                        if not main_chain_history:
                            warnings.append(
                                f"{ch}: entire leading segment before {ch}:{resseq} was excluded — no "
                                f"heavy-atom-complete residue remained to anchor a terminus."
                            )

            global_pos += 1

            terminal = ""
            terminus_reason = None
            if is_main:
                if after_gap:
                    terminal = "5"
                    terminus_reason = break_reason
                elif (resseq, icode) == first_main_key:
                    terminal = "5"
                    terminus_reason = "chain_end"
                elif (resseq, icode) == last_main_key:
                    terminal = "3"
                    terminus_reason = "chain_end"
                # Poslední reziduum řetězce hned za zlomem (nebo řetězec o
                # jediném reziduu) je 5'- i 3'-koncem zároveň.
                if terminal == "5" and (resseq, icode) == last_main_key:
                    terminal = "53"

            ff_resname, known, search_group = _pick_variant(group, resname, atoms, conv, terminal)
            missing_atoms = _check_missing_atoms(group, ff_resname, atoms, conv)
            extra_atoms = _check_extra_atoms(group, ff_resname, atoms, conv)
            if terminal in ("5", "53"):
                allowed_extra = _TERMINAL_ALLOWED_EXTRA_HEAVY_ATOMS.get(group, set())
                extra_atoms = [a for a in extra_atoms if a not in allowed_extra]

            if not known:
                warnings.append(f"Unknown residue '{resname}' at {ch}:{resseq}{icode or ''}")
            else:
                if missing_atoms:
                    warnings.append(f"Incomplete residue '{resname}' at {ch}:{resseq}: Missing {missing_atoms}")
                if extra_atoms:
                    warnings.append(f"Unexpected residue '{resname}' at {ch}:{resseq}: Extra atoms not in template {extra_atoms}")

            if after_gap:
                label = "N-terminus" if group == "P" else "5'-terminus"
                reason_text = _BREAK_REASON_LABEL.get(break_reason, break_reason)
                warnings.append(
                    f"{ch}:{resseq} ({resname}) treated as artificial {label} "
                    f"— {reason_text} before this residue."
                )

            conn_info = _check_connectivity_integrity(group, ff_resname, atoms, conv)

            token = {
                "position": global_pos,
                "chain": ch,
                "resseq": resseq,
                "icode": icode or "",
                "pdb_resname": resname,
                "is_gap": False,
                "group": group,
                "ff_resname": ff_resname,
                "known": known,
                "atoms": atoms,
                "missing_atoms": missing_atoms,
                "extra_atoms": extra_atoms,
                "is_broken": conn_info["is_broken"],
                "connectivity_parts": conn_info["components"],
                "terminus_reason": terminus_reason,
                "terminal": terminal,
            }
            tokens.append(token)

            if is_main:
                prev_resseq = resseq
                prev_icode = icode or ""
                prev_main_token = token
                main_chain_history.append(token)

    return {
        "chains": {
            ch_id: {
                "chain": ch_id,
                "tokens": [t for t in tokens if t["chain"] == ch_id],
                "warnings": [w for w in warnings if f"chain {ch_id}" in w or f"{ch_id}:" in w]
            }
            for ch_id in chains_dict.keys()
        }
    }


def _check_missing_atoms(group: Optional[str], ff_name: str, atoms: List[str], conv: Dict) -> List[str]:
    """
    POROVNÁVÁ SEZNAM ATOMŮ Z PDB SE ŠABLONOU V KONVERZNÍM SLOVNÍKU A VRACÍ SEZNAM VŠECH CHYBĚJÍCÍCH ELEMENTŮ.
    """
    res_def = _get_res_def(group, ff_name, conv)
    if not res_def:
        return []

    required_atoms = set(res_def.get("atom", {}).keys())
    actual_atoms = set(atoms)

    return sorted([a for a in required_atoms if a not in actual_atoms])


def _check_extra_atoms(group: Optional[str], ff_name: str, atoms: List[str], conv: Dict) -> List[str]:
    """
    ZRCADLOVÁ FUNKCE K _check_missing_atoms - POROVNÁ SEZNAM ATOMŮ Z PDB SE
    ŠABLONOU A VRÁTÍ ATOMY, KTERÉ V PDB JSOU NAVÍC OPROTI ŠABLONĚ (tj. atomy,
    které builder interně eviduje jako observed_extra_atoms - viz
    forge_molecule_parser.py). Bez šablony (neznámé reziduum) nelze nic
    porovnat, vrací se prázdný seznam stejně jako u _check_missing_atoms.
    """
    res_def = _get_res_def(group, ff_name, conv)
    if not res_def:
        return []

    required_atoms = set(res_def.get("atom", {}).keys())
    actual_atoms = set(atoms)

    return sorted([a for a in actual_atoms if a not in required_atoms])


def _check_connectivity_integrity(group: Optional[str], ff_name: str, atoms: List[str], conv: Dict) -> Dict[str, Any]:
    """
    ANALYZUJE GEOMETRICKOU INTEGRITU REZIDUA POMOCÍ GRAFU KONEKTIVITY A IDENTIFIKUJE IZOLOVANÉ SKUPINY ATOMŮ.
    """
    res_def = _get_res_def(group, ff_name, conv)
    if not res_def:
        return {"is_broken": False, "components": [atoms] if atoms else []}

    conn_map = res_def.get("connectivity", {})
    if not conn_map:
        return {"is_broken": False, "components": [atoms] if atoms else []}

    present_atoms = set(atoms)
    graph = {atom: set() for atom in present_atoms}
    for u, neighbors in conn_map.items():
        if u in present_atoms:
            for v in neighbors:
                if v in present_atoms:
                    graph[u].add(v)
                    graph[v].add(u)

    visited = set()
    components = []
    for start_node in sorted(list(present_atoms)):
        if start_node not in visited:
            component = []
            queue = [start_node]
            visited.add(start_node)
            while queue:
                u = queue.pop(0)
                component.append(u)
                for v in graph.get(u, []):
                    if v not in visited:
                        visited.add(v)
                        queue.append(v)
            components.append(sorted(component))

    return {
        "is_broken": len(components) > 1,
        "components": components
    }
