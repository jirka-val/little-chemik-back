"""
Výkonnostní/zátěžové testy pro ForgeStructureService na velké struktuře (1JJ2).

Nahrazuje starý tests/test_performance.py, který měřil HydrogenationService
(PDBFixer/OpenMM) - ten už žádný live endpoint nepoužívá. Marker `slow` +
`network` - mimo výchozí běh, zapneš přes `pytest -m "slow and network"`.
"""

import logging
import time

import pytest

from app.services.pdb_service import PDBService
from app.services.structure.forge_service import ForgeStructureService
from app.workspaces.manager import workspace_manager

logger = logging.getLogger(__name__)

pytestmark = [pytest.mark.slow, pytest.mark.network]

PDB_CODE = "1jj2"


@pytest.fixture(scope="module")
def large_structure_pdb():
    import asyncio
    pdb_service = PDBService()
    return asyncio.run(pdb_service.get_remote_pdb_content(PDB_CODE))


class TestDownloadPerformance:
    def test_download_and_workspace_creation_is_fast(self, large_structure_pdb):
        start = time.time()
        ws_id = workspace_manager.create_from_string(large_structure_pdb, "structure.pdb")
        duration = time.time() - start

        logger.info(f"Workspace {ws_id} created in {duration:.3f}s")
        assert workspace_manager.workspace_exists(ws_id)
        assert duration < 2.0


class TestMoleculeTypeDetectionPerformance:
    def test_type_detection_is_fast(self, large_structure_pdb):
        pdb_service = PDBService()
        start = time.time()
        detected = pdb_service.get_molecule_types(large_structure_pdb)
        duration = time.time() - start

        logger.info(f"Detected types {detected} in {duration:.3f}s")
        assert duration < 1.0


class TestForgeBuildPerformance:
    def test_forge_build_completes_within_budget(self, large_structure_pdb, offline_forge_ff):
        """
        1JJ2 je velký ribozomální komplex - RNA + protein dohromady. Tenhle
        test tedy potřebuje lokálně nacachované FF pro OBOJÍ ("R" i "P");
        pokud lokálně chybí proteinové silové pole, offline_forge_ff fixtura
        test korektně přeskočí místo matoucího pádu hluboko v builderu.

        REÁLNÝ NÁLEZ (proměřeno se skutečným FF14SB staženým přes IDA API):
        1JJ2 narazí na legitimní missing_dof přesně na K:83 (CGLU:CG) - GLU83,
        stejné reziduum probírané na začátku konverzace o 1JJ2 gapu ("dokonce
        chybí i kus GLU83"). Boční řetězec za CB v reálné struktuře chybí
        natolik, že ho nejde jednoznačně dostavět - builder se podle očekávání
        zastaví (viz tests/network/test_solvation_box.py pro plný detail a
        stejný nález). Test proto počítá i tenhle výsledek jako úspěch a
        měří čas do něj, ne čas do "complete" stavu, který 1JJ2 bez ručního
        zásahu strukturně nemůže dosáhnout.
        """
        from app.services.structure.forge_service import ForgeMissingDOFError

        service = ForgeStructureService()
        start = time.time()

        try:
            result = service.prepare_structure(
                pdb_text=large_structure_pdb,
                ff_selections={"R": {"display_name": "OL3"}, "P": {"display_name": "FF14SB"}},
                ph=7.0,
                add_solvent_and_ions=False,
            )
        except ForgeMissingDOFError as exc:
            duration = time.time() - start
            logger.info(
                f"FORGE build for {PDB_CODE} hit missing_dof at "
                f"{exc.payload['chain']}:{exc.payload['resseq']} after {duration:.2f}s"
            )
            assert duration < 300.0
            return

        duration = time.time() - start
        logger.info(f"FORGE build for {PDB_CODE} completed in {duration:.2f}s")
        assert len(result.pdb_text) > len(large_structure_pdb)
        assert duration < 300.0  # 5 minut - hrubá pojistka, ne přesný cíl
