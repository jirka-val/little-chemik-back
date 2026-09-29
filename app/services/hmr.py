"""
Hydrogen Mass Repartitioning (HMR) nad AMBER topologickým slovníkem
(výstup AMBER_topology.create_AMBER_topology, před zápisem .prmtop).

Stejné chování jako ParmEd `HMassRepartition` (výchozí volby): každý vodík
dostane hmotnost `h_mass` (3.024 Da) a o přidanou hmotnost se sníží těžký
atom, na který je vodík vázaný - celková hmotnost systému se nemění. Voda se
nechává být (je rigidní přes SETTLE, HMR jí nic nepřinese a rozbila by
parametry vodního modelu).

HMR NENÍ obyčejný přepínač v mdin (viz FORGE design, W4.5): mdin s dt=4 fs
bez přerozdělených hmot v topologii by simulaci rozstřelil, proto se to
aplikuje už tady při zápisu topologie.
"""

from __future__ import annotations

import bisect
from typing import Any, Dict

DEFAULT_HMR_H_MASS = 3.024
WATER_RESIDUE_NAMES = {"WAT", "HOH", "SOL"}


def apply_hmr(topology: Dict[str, Any], h_mass: float = DEFAULT_HMR_H_MASS) -> int:
    """
    Přerozdělí hmoty v `topology["MASS"]` na místě. Vrací počet
    přerozdělených vodíků.

    Vodíky se hledají přes BONDS_INC_HYDROGEN (indexy jsou v AMBER formátu
    3 * atom) a ATOMIC_NUMBER == 1. Vyhodí ValueError, pokud by těžkému
    atomu vyšla nekladná hmotnost (nesmyslná topologie).
    """
    masses = topology["MASS"]
    atomic_numbers = topology["ATOMIC_NUMBER"]
    residue_pointer = topology["RESIDUE_POINTER"]  # 1-based index prvního atomu rezidua
    residue_label = topology["RESIDUE_LABEL"]

    def is_water(atom: int) -> bool:
        res = bisect.bisect_right(residue_pointer, atom + 1) - 1
        return residue_label[res].strip() in WATER_RESIDUE_NAMES

    bonds = topology["BONDS_INC_HYDROGEN"]
    repartitioned = 0
    for i in range(0, len(bonds), 3):
        a, b = bonds[i] // 3, bonds[i + 1] // 3
        a_is_h = atomic_numbers[a] == 1
        b_is_h = atomic_numbers[b] == 1
        if a_is_h == b_is_h:
            # H-H (např. rigidní voda) nebo vazba bez vodíku - nic k přerozdělení.
            continue
        hydrogen, heavy = (a, b) if a_is_h else (b, a)
        if is_water(hydrogen):
            continue

        transfer = h_mass - masses[hydrogen]
        if masses[heavy] - transfer <= 0:
            raise ValueError(
                f"HMR would leave atom {heavy + 1} ({topology['ATOM_NAME'][heavy]}) with non-positive mass."
            )
        masses[hydrogen] = h_mass
        masses[heavy] -= transfer
        repartitioned += 1

    return repartitioned
