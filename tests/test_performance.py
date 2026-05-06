import pytest
import time
import logging
from pathlib import Path
from app.services.pdb_service import PDBService
from app.workspaces.manager import workspace_manager
from app.services.structure.hydrogenation import HydrogenationService  # Přidán import

logger = logging.getLogger(__name__)

# Globální proměnné pro předávání dat mezi testy
PDB_CODE = "1jj2"
SHARED_WORKSPACE_ID = None
PDB_CONTENT_CACHE = None


@pytest.mark.asyncio
async def test_01_download_and_workspace_creation():
    """
    KROK 1: Měří čas stažení velkého PDB souboru a jeho uložení na disk.
    """
    global SHARED_WORKSPACE_ID, PDB_CONTENT_CACHE
    pdb_service = PDBService()

    logger.info(f"--- Spouštím test stahování pro {PDB_CODE.upper()} ---")
    start_time = time.time()

    try:
        pdb_content = await pdb_service.get_remote_pdb_content(PDB_CODE)
        download_time = time.time() - start_time
        PDB_CONTENT_CACHE = pdb_content

        size_mb = len(pdb_content) / (1024 * 1024)
        logger.info(f"Staženo za: {download_time:.3f} s. Velikost: {size_mb:.2f} MB")

        assert download_time < 10.0, "Stahování trvalo příliš dlouho."
        assert len(pdb_content) > 0, "PDB obsah je prázdný."

        ws_start = time.time()
        SHARED_WORKSPACE_ID = workspace_manager.create_from_string(pdb_content, "structure.pdb")
        ws_time = time.time() - ws_start

        logger.info(f"Workspace {SHARED_WORKSPACE_ID} vytvořen za: {ws_time:.3f} s")
        assert workspace_manager.workspace_exists(SHARED_WORKSPACE_ID)

    except Exception as e:
        pytest.fail(f"Test 1 selhal: {e}")


def test_02_molecule_type_detection():
    """
    KROK 2: Měří čas detekce typů molekul.
    """
    global PDB_CONTENT_CACHE
    if not PDB_CONTENT_CACHE:
        pytest.skip("PDB obsah není k dispozici z předchozího testu.")

    pdb_service = PDBService()

    logger.info(f"--- Spouštím test detekce typů pro {PDB_CODE.upper()} ---")
    start_time = time.time()

    try:
        detected_types = pdb_service.get_molecule_types(PDB_CONTENT_CACHE)
        detection_time = time.time() - start_time

        logger.info(f"Typy detekovány za: {detection_time:.3f} s. Nalezeno: {detected_types}")
        assert detection_time < 1.0, "Detekce typů je příliš pomalá."

    except Exception as e:
        pytest.fail(f"Test 2 selhal: {e}")


def test_03_solvation_performance():
    """
    KROK 3: Zátěžový test solvatace.
    Testuje schopnost PDBFixeru nabalit vodu na velkou strukturu (1JJ2).
    Měří výpočetní čas a kontroluje stabilitu paměti OpenMM/PDBFixeru.
    """
    global PDB_CONTENT_CACHE, SHARED_WORKSPACE_ID
    if not PDB_CONTENT_CACHE:
        pytest.skip("PDB obsah není k dispozici.")

    hydro_service = HydrogenationService()

    logger.info(f"--- Spouštím test SOLVATACE pro {PDB_CODE.upper()} ---")
    start_time = time.time()

    try:
        # Využijeme parametry, které způsobí maximální zátěž (padding 1.0 nm přidá obrovské množství vody u 1JJ2)
        solvated_pdb_content = hydro_service.prepare_structure(
            pdb_content=PDB_CONTENT_CACHE,
            ph=7.0,
            crystal_water_mode="remove_all",
            add_solvent=True,
            box_padding_nm=1.0,  # Standardní, ale u 1JJ2 to vytvoří gigantický box
            ionic_strength=0.15,
            positive_ion="Na+",
            negative_ion="Cl-"
        )

        solvation_time = time.time() - start_time
        size_mb = len(solvated_pdb_content) / (1024 * 1024)

        logger.info(f"Solvatace úspěšně dokončena za: {solvation_time:.2f} s.")
        logger.info(f"Velikost výsledného solvatovaného PDB: {size_mb:.2f} MB")

        assert len(solvated_pdb_content) > len(PDB_CONTENT_CACHE), "Solvatované PDB by mělo být větší."

        # Volitelně můžeme to solvatované PDB zapsat do workspace pro případný další test topologie
        workspace_manager.create_from_string(solvated_pdb_content, "structure_solvated.pdb")

    except Exception as e:
        solvation_time = time.time() - start_time
        logger.error(f"Solvatace (OpenMM PDBFixer) ZKOLABOVALA po {solvation_time:.2f} s!")
        pytest.fail(f"Test 3 selhal - pád solvatace: {e}")