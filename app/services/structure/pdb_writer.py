"""Writes a builder Molecule as fixed-column PDB text plus the forge_meta sidecar."""

from __future__ import annotations

from typing import Any, Dict, List, Optional




from forge_molecule_parser import (  # noqa: E402
    Molecule,
    Residue,
    format_pdb_atom_line,
    infer_element,
)

from app.services.structure.groups import ION_GROUPS, WATER_GROUPS


class ForgeWriterError(RuntimeError):
    """Vyhozeno, když by výsledná molekula nešla bezpečně zapsat do fixed-column PDB."""


def _pdb_safe_resname(residue: Residue) -> str:
    """
    Vrátí resname zapsatelné do 3sloupcového PDB pole. RNA/DNA varianty
    (RU3/RA5/...) do 3 znaků vždy vejdou. Proteinové terminální varianty
    (CGLU/NPHE/...) jsou 4 znaky - format_pdb_atom_line je NEOŘEZÁVÁ, jen by
    tiše posunul všechny další sloupce na řádku, takže se vždy vrací
    original_resname (builder ho u téhle mutace nepřepisuje).
    """
    if len(residue.ff_resname) <= 3:
        return residue.ff_resname
    original = residue.original_resname
    if original and len(original) <= 3:
        return original
    raise ForgeWriterError(
        f"Residue {residue.chain_id}:{residue.resseq}{residue.icode} "
        f"({residue.ff_resname!r}) has no PDB-safe (<=3 char) representation."
    )


def _format_ter_line(serial: int, resname: str, chain_id: str, resseq: int, icode: str) -> str:
    """
    Standardní PDB TER záznam - signalizuje downstream nástrojům (Mol* mimo
    jiné), že tady polymerní řetězec končí. Bez něj hrozí, že se poslední
    reziduum jednoho chainu a první reziduum dalšího vyhodnotí jako přerušený
    (gap) polymer téhož řetězce, což se ve vieweru projeví tečkovanou
    "vazbou" mezi dvěma chainy i po nastavení jejich terminality - viz
    molecule_to_pdb().
    """
    return f"TER   {serial:5d}      {resname:>3s} {chain_id[:1]:1s}{resseq:4d}{icode[:1]:1s}"


def _cryst1_line(molecule: Molecule) -> Optional[str]:
    import math

    if molecule.periodic_box is None:
        return None
    vectors = molecule.periodic_box.vectors
    lengths = [math.sqrt(sum(x * x for x in vector)) for vector in vectors]

    def angle(left: int, right: int) -> float:
        dot = sum(vectors[left][i] * vectors[right][i] for i in range(3))
        cosine = max(-1.0, min(1.0, dot / (lengths[left] * lengths[right])))
        return math.degrees(math.acos(cosine))

    alpha, beta, gamma = angle(1, 2), angle(0, 2), angle(0, 1)
    return (
        f"CRYST1{lengths[0]:9.3f}{lengths[1]:9.3f}{lengths[2]:9.3f}"
        f"{alpha:7.2f}{beta:7.2f}{gamma:7.2f} P 1           1"
    )


def _is_artificial_break_terminus(residue: "Residue") -> bool:
    """
    True pro reziduum, které se stalo N- nebo C-terminálním kvůli přerušení
    uprostřed řetězce (mezera v číslování, explicitní TER, chemicky nemožná
    vzdálenost - viz analysis/sequence.py), NE proto, že by šlo o
    skutečný začátek/konec celého řetězce ("chain_end" - ten už řeší
    end-of-chain TER na konci molecule_to_pdb).
    """
    return residue.terminus_reason is not None and residue.terminus_reason != "chain_end"


def _chain_sort_key(molecule: Molecule, chain_id: str):
    residues = molecule.chains[chain_id].residues
    is_water = bool(residues) and all(r.group in WATER_GROUPS for r in residues)
    is_ion = bool(residues) and all(r.group in ION_GROUPS for r in residues)
    return (2 if is_water else 1 if is_ion else 0, chain_id)


def _pdb_serial(serial: int) -> int:
    """
    Sériové číslo atomu ve fixed-column PDB smí mít nejvýš 5 číslic
    (sloupce 7-11). format_pdb_atom_line() to samo nehlídá - `f"{serial:5d}"`
    u čísla >= 100000 tiše přeteče na 6 znaků a posune všechny další sloupce
    na řádku o jeden doprava, takže resname/chain/souřadnice skončí na
    špatné pozici. Potvrzeno pádem na solvatovaném 1JJ2 (925 179 atomů):
    OpenMM (přes PDBFixer ve StructureChecker) na takhle posunutém řádku
    spadne na "Misaligned residue name". Sériové číslo je čistě kosmetický
    popisek (nic downstream ho nepoužívá jako identitu - všude se pracuje
    přes chain/resseq/atom name), takže cyklické zabalení zpátky do rozsahu
    1-99999 je bezpečné a zachová platný fixed-column formát i nad hranicí
    legacy PDB limitu.
    """
    return ((serial - 1) % 99999) + 1


def molecule_to_pdb(molecule: Molecule) -> str:
    """
    Zapíše výsledek FORGE builderu do PDB textu pro Molstar/downstream nástroje.
    Reimplementace vzoru z (nevendorovaného) forge_builder_v0/examples/forge_workflow_cli.py,
    doplněná o bezpečnou volbu resname (viz _pdb_safe_resname) - ta v příkladovém
    CLI writeru chybí a u proteinových terminálních variant/nahrazených iontů by
    tiše poškodila fixed-column formát.
    """
    lines: List[str] = []
    cryst1 = _cryst1_line(molecule)
    if cryst1:
        lines.append(cryst1)

    serial = 1
    for chain_id in sorted(molecule.chains, key=lambda cid: _chain_sort_key(molecule, cid)):
        # _chain_sort_key()[0] je 0 jen pro "normální" (ne čistě voda/ionty)
        # chainy - TER dává smysl jen pro tyhle polymerní řetězce, HETATM
        # voda/ionty žádnou spojitost/gap logiku ve vieweru nespouští.
        chain_is_polymer = _chain_sort_key(molecule, chain_id)[0] == 0
        chain_residues = molecule.chains[chain_id].residues
        last_written: Optional[tuple] = None
        for res_idx, residue in enumerate(chain_residues):
            is_ion = residue.group in ION_GROUPS
            hetero = is_ion or residue.group in WATER_GROUPS
            for atom in residue.atoms.values():
                if atom.coord is None:
                    continue
                if is_ion:
                    # Monatomární ionty: FF resname (Mg2+, Cs+, ...) do 3 sloupců
                    # nevejde. Builder sám tenhle případ řeší přes element symbol
                    # (viz forge_molecule_ions._ion_element) - držíme se stejné
                    # konvence, ať je psaní iontů konzistentní s tím, co builder
                    # sám interně považuje za jejich identitu.
                    element = atom.element or infer_element(atom.name)
                    if not element:
                        raise ForgeWriterError(
                            f"Ion {residue.chain_id}:{residue.resseq} ({residue.ff_resname}) "
                            "has no resolvable element symbol."
                        )
                    atom_name_out = element.upper()
                    resname_out = element.upper()
                else:
                    atom_name_out = atom.name
                    resname_out = _pdb_safe_resname(residue)

                lines.append(
                    format_pdb_atom_line(
                        serial=_pdb_serial(serial),
                        record_name="HETATM" if hetero else "ATOM",
                        atom_name=atom_name_out,
                        resname=resname_out,
                        chain_id=residue.chain_id,
                        resseq=residue.resseq,
                        icode=residue.icode,
                        coord=atom.coord,
                        occupancy=atom.occupancy if atom.occupancy is not None else 1.0,
                        bfactor=atom.bfactor if atom.bfactor is not None else 0.0,
                        element=atom.element,
                        altloc="",
                    )
                )
                serial += 1
                if chain_is_polymer:
                    last_written = (resname_out, residue.resseq, residue.icode)

            # Hranice mezery uprostřed řetězce (GLU83/PHE89 na 1JJ2 apod.):
            # tohle reziduum i to bezprostředně následující jsou OBĚ umělé
            # terminusy (viz _is_artificial_break_terminus) přesně tehdy, když
            # mezi nimi je gap/TER/geometrický zlom - "chain_end" (skutečný
            # začátek/konec řetězce) tuhle podmínku nikdy nesplní, protože je
            # vždy jen na jednom z dvojice sousedů. Bez explicitního TER by
            # Molstar (a další downstream nástroje) mohly tenhle úsek
            # vyhodnotit jako spojitý polymer navzdory nastavené terminalitě.
            if (
                chain_is_polymer
                and last_written is not None
                and _is_artificial_break_terminus(residue)
                and res_idx + 1 < len(chain_residues)
                and _is_artificial_break_terminus(chain_residues[res_idx + 1])
            ):
                resname_out, resseq, icode = last_written
                lines.append(_format_ter_line(_pdb_serial(serial), resname_out, chain_id, resseq, icode))
                serial += 1

        if chain_is_polymer and last_written is not None:
            resname_out, resseq, icode = last_written
            lines.append(_format_ter_line(_pdb_serial(serial), resname_out, chain_id, resseq, icode))
            serial += 1

    for record in molecule.passthrough_atoms:
        is_ion = record.group in ION_GROUPS
        if is_ion:
            element = record.element or infer_element(record.atom_name)
            atom_name_out = (element or record.atom_name).upper()
            resname_out = (element or record.resname[:3]).upper()
        else:
            atom_name_out = record.atom_name
            resname_out = record.resname[-3:] if len(record.resname) > 3 else record.resname

        lines.append(
            format_pdb_atom_line(
                serial=_pdb_serial(serial),
                record_name="HETATM",
                atom_name=atom_name_out,
                resname=resname_out,
                chain_id=record.chain_id,
                resseq=record.resseq,
                icode=record.icode,
                coord=record.coord,
                occupancy=record.occupancy if record.occupancy is not None else 1.0,
                bfactor=record.bfactor if record.bfactor is not None else 0.0,
                element=record.element,
                altloc="",
            )
        )
        serial += 1

    lines.extend(("END", ""))
    return "\n".join(lines)


def build_forge_meta(molecule: Molecule) -> Dict[str, Dict[str, Any]]:
    """
    Sidecar metadata (chain:resseq:icode -> autoritativní ff_resname/group) pro
    TopologyService. Builder do PDB textu zapisuje jen 3znakovou reprezentaci
    (viz _pdb_safe_resname), takže proteinové terminální varianty jako CGLU/NPHE
    by se po zpětném parsování PDB ztratily. Tohle je jediné místo, kde je
    plná (i 4znaková) identita reziduí po doběhnutí state-assignmentu dostupná.
    """
    meta: Dict[str, Dict[str, Any]] = {}
    for chain in molecule.chains.values():
        for residue in chain.residues:
            key = f"{residue.chain_id}:{residue.resseq}:{residue.icode}"
            meta[key] = {"ff_resname": residue.ff_resname, "group": residue.group}
    return meta
