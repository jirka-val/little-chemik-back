import numpy as np
import math
from typing import List, Any

#  Konstanty
NM_TO_ANGSTROM = 10.0
TIP5P_ANGLE_RAD = math.radians(54.735)


def water_extra_points(ff: Any, water_crd: List[List[float]]) -> List[List[float]]:
    """
    Vypočítá pozice Extra Points (EP) pro jeden model vody (4-site nebo 5-site).

    Očekává:
    - ff: Načtený objekt silového pole (instance z knihovny FF_IDA)
    - water_crd: Seznam 3 souřadnic v pořadí [Kyslík, Vodík1, Vodík2]

    Vrací:
    - Seznam souřadnic pro vypočítané Extra Points (1 nebo 2 EP).
    - Pokud model EP nevyžaduje (např. běžný 3-site), vrací prázdný seznam [].
    """
    # 1. Kontrola topologie vody
    if 'WAT' not in ff.units:
        raise ValueError("Error: 'WAT' topology missing in the provided force field.")

    nEP = ff.units['WAT']['atoms']['type'].count('WEP')

    # 3-site model nemá žádné Extra Points, vracíme prázdné pole
    if nEP == 0:
        return []

    if nEP > 2:
        raise NotImplementedError("Error: Models with more than 2 EPs are not supported.")

    # 2. Bezpečné vyhledání délky vazby (nezáleží na pořadí atomů WEP-WO vs WO-WEP)
    bond_types = ff.b.get('bondtypes', {})
    bond = bond_types.get(('WEP', 'WO')) or bond_types.get(('WO', 'WEP'))

    if not bond:
        raise ValueError("Error: No WEP-WO or WO-WEP bonded term found in force field.")

    # Převod vzdálenosti EP z nanometrů na Ångströmy
    r = bond[1] * NM_TO_ANGSTROM

    # 3. Načtení a příprava vektorů
    O = np.asarray(water_crd[0], dtype=float)
    H1 = np.asarray(water_crd[1], dtype=float)
    H2 = np.asarray(water_crd[2], dtype=float)

    # Výpočet směrových vektorů z kyslíku na vodíky
    OH1 = H1 - O
    OH2 = H2 - O

    norm_OH1 = np.linalg.norm(OH1)
    norm_OH2 = np.linalg.norm(OH2)

    if norm_OH1 == 0 or norm_OH2 == 0:
        raise ValueError("Error: Invalid water coordinates (overlapping atoms).")

    OH1 /= norm_OH1
    OH2 /= norm_OH2

    # Vektor půlící úhel (osa symetrie molekuly vody)
    v1 = OH1 + OH2
    v1 /= np.linalg.norm(v1)

    # Vektor kolmý na rovinu molekuly
    v2 = np.cross(OH1, OH2)
    v2 /= np.linalg.norm(v2)

    # 4. Výpočet finálních pozic EP
    if nEP == 1:
        # 4-site water model (např. TIP4P) - 1 EP na ose půlící úhel H-O-H
        EP = O + r * v1
        return [EP.tolist()]

    elif nEP == 2:
        # 5-site water model (např. TIP5P) - 2 EP v rovině kolmé na molekulu
        c = math.cos(TIP5P_ANGLE_RAD)
        s = math.sin(TIP5P_ANGLE_RAD)

        EP1 = O - r * v1 * c + r * v2 * s
        EP2 = O - r * v1 * c - r * v2 * s
        return [EP1.tolist(), EP2.tolist()]