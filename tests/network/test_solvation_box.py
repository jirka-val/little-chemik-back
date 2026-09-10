"""
Solvatace a tvorba periodického boxu přes ForgeStructureService - obdoba
starého tests/test_performance.py::test_03_solvation_performance (který dělal
fixer.addSolvent() přes PDBFixer/OpenMM), teď nad FORGE builderem.

Důležitý detail objevený při psaní tohoto testu: builder interně hledá
vodní silové pole striktně pod mol_type "W3" (viz
forge_molecule_solvation.py::solvate_molecule -> parameters.sigma("W3", ...)),
ne pod obecným "W", které používá starší TopologyService/pdb_service pipeline
pro AMBER topologii. ff_selections pro ForgeStructureService proto musí
používat builderovu vlastní vokabulář mol_type ("R"/"D"/"P"/"W3"/"W4"/"W5"/
"I1"/"I1+"/"Im"/"Im+"), ne zjednodušené kódy ("W"/"I") z pdb_service.get_molecule_types.
"""

import pytest

from app.workspaces.manager import workspace_manager

pytestmark = pytest.mark.network


class TestSolvationCreatesBox:
    """
    1RNA - máme lokálně RNA (OL3) i vodní (TIP3P) silové pole, takže tohle
    běží dnes reálně, ne jen jako skip.
    """

    def test_solvation_produces_periodic_box_and_waters(self, client, offline_forge_ff):
        ws_id = client.get("/api/molecules/fetch-pdb/1RNA").json()["workspace_id"]

        response = client.post("/api/validation/prepare", json={
            "workspace_id": ws_id,
            "ff_selections": {
                "R": {"display_name": "OL3"},
                "W3": {"display_name": "TIP3P"},
                "I1": {"display_name": "JC-TIP3P-I1"},
            },
            "ph": 7.0,
            "add_solvent": True,
            "box_padding_nm": 1.0,
            "box_shape": "cube",
        })
        assert response.status_code == 200, response.text

        pdb_text = workspace_manager.get_file_path(ws_id, "structure.pdb").read_text(encoding="utf-8")
        assert pdb_text.startswith("CRYST1"), "Chybí periodický box (CRYST1 řádek)"

        water_oxygens = sum(
            1 for l in pdb_text.splitlines()
            if l.startswith("HETATM") and l[17:20].strip() == "WAT" and l[12:16].strip() == "O"
        )
        assert water_oxygens > 0, "Solvatace nepřidala žádné molekuly vody"

    def test_solvation_report_reflects_requested_box_shape(self, client, offline_forge_ff):
        from app.services.structure.forge_service import ForgeStructureService
        from app.services.pdb_service import PDBService
        import asyncio

        pdb_text = asyncio.run(PDBService().get_remote_pdb_content("1rna"))

        result = ForgeStructureService().prepare_structure(
            pdb_text=pdb_text,
            ff_selections={
                "R": {"display_name": "OL3"},
                "W3": {"display_name": "TIP3P"},
                "I1": {"display_name": "JC-TIP3P-I1"},
            },
            ph=7.0,
            add_solvent_and_ions=True,
            box_shape="cubic",
            box_padding_angstrom=10.0,
        )

        assert result.solvation is not None
        assert result.solvation.box_shape == "cubic"
        assert result.solvation.total_waters > 0
        assert result.solvation.padding_angstrom == pytest.approx(10.0)


@pytest.mark.slow
class TestLargeStructureSolvationPerformance:
    """
    1JJ2 - parita se starým výkonnostním testem (velký box, hodně vody).

    DŮLEŽITÉ ZJIŠTĚNÍ Z REÁLNÉHO BĚHU (staženo skutečné FF14SB přes IDA API,
    proměřeno na tomhle stroji): 1JJ2 je přesně ten typ struktury, který
    INTEGRATION_CONTRACT.md popisuje jako "one expected missing-DOF example" -
    builder legitimně narazí na chybějící stupeň volnosti a ZASTAVÍ SE dřív,
    než vůbec začne solvatace/box - přesně podle kontraktu ("no ion cleanup,
    coordinate recentering, solvation or ion placement has occurred").

    Konkrétně narazí přesně na K:83 (CGLU:CG) po ~92 s - tj. na GLU83, to samé
    reziduum z úplně první diskuze o 1JJ2 v týhle konverzaci ("dokonce chybí
    i kus GLU83"). Fáze 3 správně rozpoznala terminální variantu (CGLU), ale
    boční řetězec od CG dál v reálné struktuře chybí natolik, že ho nejde
    jednoznačně dostavět - takže se builder podle očekávání zastaví, místo
    aby něco vymyslel. Tohle je živé, end-to-end potvrzení přesně toho
    scénáře, který na začátku popsal šéf.

    Solvatace samotná se tedy na 1JJ2 beze změny vstupu NEDÁ změřit - je to
    očekávaný, správný výsledek, ne bug. Test proto ověřuje TOHLE chování
    (missing_dof v rozumném čase), ne úspěšnou solvataci. Pro reálné číslo
    "solvatace + box na velké struktuře" viz TestSolvationCreatesBox výše
    (1RNA, dokončí se v pořádku).
    """

    def test_1jj2_hits_missing_dof_before_solvation_within_time_budget(self, client, offline_forge_ff):
        import time
        import logging
        from app.services.structure.forge_service import ForgeStructureService, ForgeMissingDOFError
        from app.services.pdb_service import PDBService
        import asyncio

        logger = logging.getLogger(__name__)
        pdb_text = asyncio.run(PDBService().get_remote_pdb_content("1jj2"))

        start = time.time()
        try:
            result = ForgeStructureService().prepare_structure(
                pdb_text=pdb_text,
                ff_selections={
                    "R": {"display_name": "OL3"},
                    "P": {"display_name": "FF14SB"},
                    "W3": {"display_name": "TIP3P"},
                    "I1": {"display_name": "JC-TIP3P-I1"},
                },
                ph=7.0,
                add_solvent_and_ions=True,
                box_padding_angstrom=10.0,
            )
        except ForgeMissingDOFError as exc:
            duration = time.time() - start
            logger.info(
                f"1JJ2 hit missing_dof at {exc.payload['chain']}:{exc.payload['resseq']} "
                f"({exc.payload['ff_resname']}:{exc.payload['atom_name']}) after {duration:.2f}s"
            )
            assert duration < 300.0  # hrubá pojistka, ne přesný cíl
            return

        # Kdyby jednou 1JJ2 nebo tahle verze builderu missing_dof neměly,
        # ověříme aspoň, že solvatace/box skutečně proběhly.
        duration = time.time() - start
        logger.info(f"1JJ2 solvation completed in {duration:.2f}s, {result.solvation.total_waters} waters added")
        assert result.solvation.total_waters > 0
        assert len(result.pdb_text) > len(pdb_text)
