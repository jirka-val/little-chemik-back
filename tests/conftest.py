"""
Sdílené fixtures pro celou testovací sadu.

Rozvržení testů (viz tests/README.md pro detaily):
- tests/unit/         čistá logika, žádné I/O, žádná síť
- tests/integration/  přes FastAPI TestClient, lokální workspace lifecycle
- tests/network/      potřebuje skutečné RCSB/IDA API (mark: network)
- tests/performance/  dlouho běžící zátěžové testy (mark: slow)

Klíčová bezpečnostní zásada těchto testů: NIKDY nezapisovat do reálné
`data/ff_cache/` ani `data/ff_cache_forge/` (jsou to skutečná, draze stažená
data z IDA API). Kde testy potřebují silové pole pro FORGE builder, použij
fixturu `offline_forge_ff`, která si postaví izolovanou kopii v tmp_path a
test bez ní elegantně přeskočí (`pytest.skip`), místo aby cokoliv riskovala.
"""

from __future__ import annotations

import shutil
import sys
from pathlib import Path
from typing import Any, Callable, Dict

import pytest
from fastapi.testclient import TestClient

BACKEND_ROOT = Path(__file__).resolve().parents[1]
FIXTURES_PDB_DIR = Path(__file__).resolve().parent / "fixtures" / "pdb"
REAL_FF_CACHE_DIR = BACKEND_ROOT / "data" / "ff_cache"

if str(BACKEND_ROOT) not in sys.path:
    sys.path.insert(0, str(BACKEND_ROOT))

from app.main import app  # noqa: E402
from app.workspaces.manager import workspace_manager  # noqa: E402
from app.services.forcefield_service import ForceFieldService  # noqa: E402


# ---------------------------------------------------------------------------
# HTTP klient
# ---------------------------------------------------------------------------

@pytest.fixture(scope="session")
def client() -> TestClient:
    """Sdílený synchronní TestClient - FastAPI app se inicializuje jen jednou za session."""
    return TestClient(app)


# ---------------------------------------------------------------------------
# PDB fixture data (tests/fixtures/pdb/*.pdb)
# ---------------------------------------------------------------------------

def _load_pdb(name: str) -> str:
    path = FIXTURES_PDB_DIR / name
    if not path.exists():
        raise FileNotFoundError(f"Missing PDB fixture: {path}")
    return path.read_text(encoding="utf-8")


@pytest.fixture
def pdb_alanine_single() -> str:
    """Jediné ALA reziduum - nejmenší platný protein fragment."""
    return _load_pdb("alanine_single.pdb")


@pytest.fixture
def pdb_1jj2_glu_gap_boundary() -> str:
    """
    Reálný výřez ze staženého 1JJ2 (chain K, residua 78-92), přesně obsahující
    ALA78..ALA82-GLU83-[reálná mezera 84-88 z deponované struktury]-PHE89..
    ASP92 - GLU83 v ATOM záznamech skutečně má jen N/CA/C/O/CB (CG/CD/OE1/OE2
    chybí), přesně scénář z app/builder/INTEGRATION_CONTRACT.md ("GLU83, ...
    -> CGLU"). Na rozdíl od pdb_glu_open_branch (syntetická geometrie) je
    tohle 1:1 z RCSB (https://files.rcsb.org/download/1JJ2.pdb) - použij tohle
    kdykoliv je potřeba ověřit chování na reálné, ne izolované, hraniční
    situaci (gap + neúplný postranní řetězec zároveň).
    """
    return _load_pdb("protein_1jj2_glu_gap_boundary.pdb")


@pytest.fixture
def pdb_glu_open_branch() -> str:
    """
    Jediné GLU reziduum se stejnou (reálnou, ne syntetickou) N/CA/C/O/CB
    geometrií jako alanine_single.pdb, jen přejmenované z ALA na GLU - CG/CD/
    OE1/OE2 tak chybí a builder na nich narazí na přesně jeden bezpečný
    "residue_local_open_branch" (3 DOF, viz analýza v sidechain_service
    diagnostice) - stejný typ chybějícího postranního řetězce jako u GLU83 ve
    skutečném 1JJ2 chain K, jen izolovaný do jednoho rezidua pro rychlý test.
    """
    return _load_pdb("glu_open_branch.pdb")


@pytest.fixture
def pdb_protein_gap() -> str:
    """GLU83 ... [84-88 chybí] ... PHE89 - syntetická obdoba 1JJ2 chain K."""
    return _load_pdb("protein_gap.pdb")


@pytest.fixture
def pdb_protein_gap_incomplete_boundary() -> str:
    """
    ALA82 (kompletní) - GLU83 (jen N/CA/C/O/CB, bez CG/CD/OE1/OE2) - [84-88
    chybí] - PHE89 (kompletní) - věrná reprodukce skutečného 1JJ2 chain K
    (ověřeno přímo na staženém 1JJ2.pdb z RCSB, ne jen odhadem). Na rozdíl od
    protein_gap.pdb má GLU83 tady záměrně neúplný postranní řetězec, aby šlo
    otestovat kaskádové vylučování neúplných okrajových reziduí.
    """
    return _load_pdb("protein_gap_incomplete_boundary.pdb")


@pytest.fixture
def pdb_rna_gap() -> str:
    """
    Výřez ze skutečné struktury 1RNA (chain A, residua 1-3 a 6-8, 4-5 vynechána)
    - reálná geometrie, ne syntetické souřadnice. FORGE builder na degenerovaných/
    kolineárních souřadnicích spadne ("Cannot normalize near-zero central
    dihedral bond"), takže tahle fixtura musí mít skutečné, ne vymyšlené atomy.
    """
    return _load_pdb("rna_gap.pdb")


@pytest.fixture
def pdb_altloc_sample() -> str:
    """SER42 se dvěma alternativními konformacemi OG (occupancy 0.6/0.4)."""
    return _load_pdb("altloc_sample.pdb")


@pytest.fixture
def pdb_two_model_nmr() -> str:
    """Stejné reziduum ve dvou MODEL blocích (NMR ensemble)."""
    return _load_pdb("two_model_nmr.pdb")


# ---------------------------------------------------------------------------
# Workspace lifecycle
# ---------------------------------------------------------------------------

@pytest.fixture
def make_workspace() -> Callable[[str], str]:
    """
    Vrátí funkci, která z PDB textu vytvoří skutečný workspace (přes tu samou
    cestu, kterou používá produkční kód) a vrátí jeho workspace_id.
    """
    def _make(pdb_text: str, filename: str = "structure.pdb") -> str:
        return workspace_manager.create_from_string(pdb_text, filename)
    return _make


# ---------------------------------------------------------------------------
# FORGE force-field bridging - offline, bez sítě, bez zásahu do reálné cache
# ---------------------------------------------------------------------------

@pytest.fixture
def offline_forge_ff(monkeypatch, tmp_path_factory):
    """
    Monkeypatchuje ForceFieldService.prepare_forge_force_field_directory tak,
    aby žádný test nesahal na síť ani nepřepisoval reálnou data/ff_cache/ (to
    jsou skutečná, draze stažená data z IDA API - testovací ff_data se skoro
    jistě nezakóduje/nedekóduje 1:1 a přepsala by je prázdným/špatným obsahem).

    Místo toho si pro daný ff_name/mol_type postaví izolovanou kopii přímo
    z toho, co je už lokálně nacachované v data/ff_cache/{ff_name}/ - přesně
    tou samou transformací (přejmenování na residue_lib_*/forcefield_*), jakou
    dělá ForceFieldService.prepare_forge_force_field_directory ve skutečnosti.

    Pokud lokální cache pro daný FF neexistuje (čerstvý checkout bez
    doposud stažených silových polí), test se elegantně přeskočí místo pádu.
    """
    def _fake_prepare(self, ff_data: Dict[str, Any], mol_type: str) -> Path:
        ff_name = (ff_data.get("display_name") or ff_data.get("ff_name") or "unknown_ff").replace(" ", "_")
        source_dir = REAL_FF_CACHE_DIR / ff_name
        if not source_dir.exists():
            pytest.skip(
                f"data/ff_cache/{ff_name} not present locally - this test needs a "
                f"force field that was already downloaded through the IDA API at least once."
            )

        target_dir = tmp_path_factory.mktemp("forge_ff") / f"{ff_name}_{mol_type}"
        target_dir.mkdir(parents=True, exist_ok=True)
        shutil.copyfile(source_dir / f"{ff_name}.rtp", target_dir / f"residue_lib_{ff_name}.rtp")
        shutil.copyfile(source_dir / f"{ff_name}.rtp", target_dir / f"forcefield_{ff_name}.itp")
        shutil.copyfile(source_dir / f"nonbonded_{ff_name}.itp", target_dir / f"nonbonded_{ff_name}.itp")
        shutil.copyfile(source_dir / f"bonded_{ff_name}.itp", target_dir / f"bonded_{ff_name}.itp")
        return target_dir

    monkeypatch.setattr(ForceFieldService, "prepare_forge_force_field_directory", _fake_prepare)


def has_local_ff_cache(ff_name: str) -> bool:
    """Pomocná funkce pro `pytest.mark.skipif` na testy vázané na konkrétní FF."""
    return (REAL_FF_CACHE_DIR / ff_name).exists()
