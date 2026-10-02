"""Applies the chosen model, symmetry, AltLocs and removed chain copies to
the PDB text (the Conformations step)."""

from __future__ import annotations

from typing import List, Optional


def clean_pdb_altlocs(pdb_text: str, user_selection: dict) -> str:
    """
    PROJDE PDB SOUBOR A SMAŽE VŠECHNY ALTERNATIVNÍ ATOMY I JEJICH ANISOU ZÁZNAMY,
    KTERÉ UŽIVATEL NEVYBRAL. VYBRANÝM ATOMŮM ODSTRANÍ ALTLOC INDIKÁTOR
    A NASTAVÍ OBSAZENOST (OCCUPANCY) ZPĚT NA 1.00.
    """
    cleaned_lines = []

    for line in pdb_text.splitlines():
        # OPRAVA 1: Přidali jsme kontrolu i pro ANISOU řádky
        if line.startswith("ATOM") or line.startswith("HETATM") or line.startswith("ANISOU"):
            alt_loc = line[16]

            # Pokud má atom/anisou alternativní pozici
            if alt_loc != ' ':
                chain = line[21].strip() or "?"
                resseq_raw = line[22:26].strip()
                resname = line[17:20].strip()

                # Vygenerujeme stejný unikátní klíč, jaký nám posílá frontend (např. "A_45_TYR")
                unique_key = f"{chain}_{resseq_raw}_{resname}"

                # Pokud uživatel pro tento zbytek provedl volbu
                if unique_key in user_selection:
                    chosen_altloc = user_selection[unique_key]

                    # Pokud se AltLoc atomu neshoduje s volbou uživatele -> Smazat řádek (nepřidá se do cleaned_lines)
                    if alt_loc != chosen_altloc:
                        continue
                    else:
                        # Pokud je to ten správný atom/anisou, nahradíme jeho AltLoc znak ('A'/'B') za mezeru
                        line = line[:16] + ' ' + line[17:]

                        # OPRAVA 2: Pokud jde o ATOM/HETATM, musíme vrátit obsazenost (Occupancy) na 1.00
                        # (ANISOU řádky pole Occupancy nemají, proto ta podmínka)
                        if line.startswith("ATOM") or line.startswith("HETATM"):
                            if len(line) >= 60:
                                # Obsazenost je v PDB formátu přesně na pozicích 55-60 (index 54:60)
                                line = line[:54] + "  1.00" + line[60:]

        # Přidáme upravený (nebo nedotčený) řádek do nového čistého souboru
        cleaned_lines.append(line)

    return '\n'.join(cleaned_lines)


def process_structure(
    pdb_text: str,
    target_model: int,
    apply_symmetry: bool,
    selection: dict,
    remove_chains: Optional[List[str]] = None,
) -> str:
    """
    Kombinovaná funkce pro kompletní fyzickou přípravu PDB souboru:
    1. Ponechá pouze vybraný MODEL (např. NMR ensemble).
    2. Vymaže nevybrané alternativní pozice a upraví obsazenost na 1.0.
    3. Odstraní řetězce z remove_chains (nechtěné identické kopie molekuly,
       viz find_identical_chain_copies) včetně jejich iontů/ligandů/vod.
    4. Pokud apply_symmetry=True, vybuduje plnou biologickou jednotku (Biological Assembly)
       pomocí BIOMT matic a matematicky dopočítá atomy.
    """
    lines = pdb_text.splitlines()
    removed_chains = {c.strip() or "?" for c in (remove_chains or [])}

    # --- KROK 1: Přečtení BIOMT matic z REMARK 350 ---
    matrices = {}
    for line in lines:
        if line.startswith("REMARK 350   BIOMT"):
            try:
                row_idx = int(line[18:19])  # Řádek matice (1, 2 nebo 3)
                mat_num = int(line[20:24].strip())
                parts = line[24:].split()
                if len(parts) >= 4:
                    if mat_num not in matrices:
                        matrices[mat_num] = [[0.0] * 4 for _ in range(3)]
                    matrices[mat_num][row_idx - 1] = [float(x) for x in parts[:4]]
            except Exception:
                continue

    biomt_list = list(matrices.values())
    # Pokud v PDB není BIOMT matice, dodáme základní (Identitu - x*1 = x)
    if not biomt_list:
        biomt_list = [[[1.0, 0.0, 0.0, 0.0], [0.0, 1.0, 0.0, 0.0], [0.0, 0.0, 1.0, 0.0]]]

    # --- KROK 2: Extrakce konkrétního modelu a AltLocs filtrace ---
    in_model = False
    has_models = any(l.startswith("MODEL ") for l in lines)

    # Pokud soubor nemá vůbec tagy MODEL (typická krystalografie), čteme atomy rovnou
    if not has_models:
        in_model = True

    base_atoms = []
    for line in lines:
        if line.startswith("MODEL "):
            try:
                current_model = int(line[6:].strip())
                in_model = (current_model == target_model)
            except ValueError:
                pass
            continue
        elif line.startswith("ENDMDL"):
            in_model = False
            continue

        if in_model and (line.startswith("ATOM  ") or line.startswith("HETATM")):
            alt_loc = line[16]
            chain = line[21]
            if (chain.strip() or "?") in removed_chains:
                continue
            resseq = line[22:26].strip()
            resname = line[17:20].strip()

            # Vytvoření unikátního klíče z frontendu
            key = f"{chain.strip() or '?'}_{resseq}_{resname}"

            if alt_loc != ' ':
                chosen_alt = selection.get(key)
                # Pokud na této pozici sedí AltLoc a uživatel vybral něco jiného, přeskočíme atom
                if chosen_alt and alt_loc != chosen_alt:
                    continue
                # Úprava řádku: Vymazání AltLoc písmene (mezera) a fixace obsazenosti na 1.00
                line = line[:16] + ' ' + line[17:54] + "  1.00" + line[60:]

            base_atoms.append(line)

    # --- KROK 3: Fyzické budování biologické jednotky a přečíslování ---
    # Serial má v PDB jen 5 sloupců; nad 99999 se přetáčí (jako AMBER/FORGE),
    # jinak by se posunuly všechny další sloupce (incident 20261001-061334-ce249e).
    final_lines = []

    if apply_symmetry and len(biomt_list) > 1:
        # Pokud násobíme řetězce (A), musíme zkopírované kusy přejmenovat na nová písmena (B, C atd.)
        unique_chains = []
        for line in base_atoms:
            ch = line[21]
            if ch not in unique_chains:
                unique_chains.append(ch)

        chain_alphabet = "ABCDEFGHIJKLMNOPQRSTUVWXYZabcdefghijklmnopqrstuvwxyz0123456789"
        used_chains = set(unique_chains)
        available_chains = [c for c in chain_alphabet if c not in used_chains]

        # Mapa pro dynamické přidělování názvů řetězců
        matrix_chain_maps = {0: {ch: ch for ch in unique_chains}}  # 1. matice nemění název
        for i in range(1, len(biomt_list)):
            matrix_chain_maps[i] = {}
            for ch in unique_chains:
                if available_chains:
                    new_c = available_chains.pop(0)
                    matrix_chain_maps[i][ch] = new_c
                else:
                    matrix_chain_maps[i][ch] = ch  # Došly znaky, použijeme starý

        # Násobení: Vezmeme každou transformační matici a převalíme přes ni všechny vyčištěné atomy
        atom_serial = 1
        for i, matrix in enumerate(biomt_list):
            for line in base_atoms:
                try:
                    # Rozparsování původních souřadnic
                    x = float(line[30:38])
                    y = float(line[38:46])
                    z = float(line[46:54])

                    # Maticové násobení + translace
                    nx = matrix[0][0] * x + matrix[0][1] * y + matrix[0][2] * z + matrix[0][3]
                    ny = matrix[1][0] * x + matrix[1][1] * y + matrix[1][2] * z + matrix[1][3]
                    nz = matrix[2][0] * x + matrix[2][1] * y + matrix[2][2] * z + matrix[2][3]

                    # Nové ID řetězce
                    orig_ch = line[21]
                    new_ch = matrix_chain_maps[i].get(orig_ch, orig_ch)

                    # Složení nového řádku zpět podle striktního PDB formátu
                    new_line = (
                            line[:6] +
                            f"{atom_serial % 100000:5d}" +
                            line[11:21] +
                            new_ch +
                            line[22:30] +
                            f"{nx:8.3f}{ny:8.3f}{nz:8.3f}" +
                            line[54:]
                    )
                    final_lines.append(new_line)
                    atom_serial += 1
                except Exception:
                    final_lines.append(line)
    else:
        # Přečíslování atomů
        atom_serial = 1
        for line in base_atoms:
            try:
                new_line = line[:6] + f"{atom_serial % 100000:5d}" + line[11:]
                final_lines.append(new_line)
                atom_serial += 1
            except Exception:
                final_lines.append(line)

    final_lines.append("END")
    return "\n".join(final_lines)
