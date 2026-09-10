"""
Integrační testy pro POST /api/sidechains/* - interaktivní dostavování
bezpečných "residue_local_open_branch" side-chainů (viz
app/services/structure/sidechain_service.py a
app/builder/INTEGRATION_CONTRACT.md).

Používá `offline_forge_ff` (žádná síť, žádný zápis do reálné data/ff_cache/) a
fixturu `pdb_glu_open_branch` - jediné GLU reziduum se stejnou reálnou N/CA/C/
O/CB geometrií jako alanine_single.pdb, jen přejmenované z ALA (viz
conftest.py) - CG/CD/OE1/OE2 tak chybí a builder na nich narazí na přesně
jeden bezpečný open-branch (ověřeno přímo proti sidechain_service, ne jen
odhadem: 3 DOF - CG/CD/OE1 dihedraly, 8 chybějících atomů).
"""

import json

import pytest

from app.workspaces.manager import workspace_manager

pytestmark = pytest.mark.integration

_FF_SELECTIONS = {"P": {"display_name": "FF14SB"}}


def _start_payload(ws_id: str, **overrides):
    payload = {
        "workspace_id": ws_id,
        "ff_selections": _FF_SELECTIONS,
        "ph": 7.0,
        "add_solvent": False,
    }
    payload.update(overrides)
    return payload


class TestSidechainStart:
    def test_start_opens_gui_session_for_open_branch(
        self, client, offline_forge_ff, make_workspace, pdb_glu_open_branch
    ):
        ws_id = make_workspace(pdb_glu_open_branch)

        response = client.post(f"/api/sidechains/start/{ws_id}", json=_start_payload(ws_id))

        assert response.status_code == 200, response.text
        body = response.json()
        assert body["status"] == "missing_dof"
        assert body["preview_filename"] == "structure_preview.pdb"

        sidechains = body["gui_payload"]["sidechains"]
        assert len(sidechains) == 1
        residue = sidechains[0]
        assert residue["residue"]["chain_id"] == "A"
        dof_names = sorted(dof["dof_key"]["atom"] for dof in residue["dofs"])
        assert dof_names == ["CD", "CG", "OE1"]
        for dof in residue["dofs"]:
            assert dof["value_degrees"] == dof["default_degrees"]

        preview_path = workspace_manager.get_file_path(ws_id, "structure_preview.pdb")
        assert preview_path.exists()
        preview_atom_count = sum(1 for l in preview_path.read_text().splitlines() if l.startswith("ATOM"))
        assert preview_atom_count > 5  # víc než jen původní N/CA/C/O/CB - preview má FF-optimální dostavěné side-chainy

        # structure.pdb (kanonický soubor workspace) nesmí být commitem nedotčen,
        # dokud uživatel GUI relaci nepřijme.
        canonical = workspace_manager.get_file_path(ws_id, "structure.pdb").read_text()
        assert sum(1 for l in canonical.splitlines() if l.startswith("ATOM")) == 5

    def test_start_behaves_like_prepare_when_nothing_is_missing(
        self, client, offline_forge_ff, make_workspace, pdb_alanine_single
    ):
        ws_id = make_workspace(pdb_alanine_single)

        response = client.post(f"/api/sidechains/start/{ws_id}", json=_start_payload(ws_id))

        assert response.status_code == 200, response.text
        body = response.json()
        assert body["status"] == "complete"
        assert body["message"] == "Structure successfully prepared."
        assert not workspace_manager.get_file_path(ws_id, "structure_preview.pdb").exists()


class TestSidechainUpdateAndOptimize:
    def _start(self, client, ws_id):
        response = client.post(f"/api/sidechains/start/{ws_id}", json=_start_payload(ws_id))
        assert response.status_code == 200, response.text
        return response.json()["gui_payload"]

    def test_update_returns_small_coordinate_patch(
        self, client, offline_forge_ff, make_workspace, pdb_glu_open_branch
    ):
        ws_id = make_workspace(pdb_glu_open_branch)
        gui_payload = self._start(client, ws_id)
        preview_path = workspace_manager.get_file_path(ws_id, "structure_preview.pdb")
        preview_before = preview_path.read_text()

        dof = gui_payload["sidechains"][0]["dofs"][0]
        dof["value_degrees"] = dof["default_degrees"] + 30.0

        response = client.post(f"/api/sidechains/update/{ws_id}", json={
            "sidechain_data": gui_payload,
            "changed_dof_keys": [dof["dof_key"]],
        })

        assert response.status_code == 200, response.text
        # structure_preview.pdb musí odrážet novou geometrii, ne zůstat na
        # počátečním stavu z /start - jinak by frontendový plný-reload fallback
        # (viewer.applyCoordinatePatch) ukazoval zastaralou strukturu.
        assert preview_path.read_text() != preview_before
        patch = response.json()
        assert len(patch["updated_atoms"]) > 0
        # Patch je jen atomy postiženého rezidua/DOF, ne celá struktura -
        # přesně tenhle rozdíl (řádově jednotky/desítky atomů vs. celá
        # struktura) byl smyslem celé side-chain GUI funkce.
        assert len(patch["updated_atoms"]) < 20
        for entry in patch["updated_atoms"]:
            assert entry["atom"]["chain_id"] == "A"
            assert len(entry["coord"]) == 3

    def test_optimize_returns_dof_values_and_patch(
        self, client, offline_forge_ff, make_workspace, pdb_glu_open_branch
    ):
        ws_id = make_workspace(pdb_glu_open_branch)
        gui_payload = self._start(client, ws_id)
        residue = gui_payload["sidechains"][0]["residue"]

        # /start už dostaví na FF optimum, takže Opt hned po startu je legitimně
        # no-op (nic se nezmění) - nejdřív reziduum ručně "rozhýbeme" (jako
        # slider tah), ať má Opt co reálně re-optimalizovat.
        dof = gui_payload["sidechains"][0]["dofs"][0]
        dof["value_degrees"] = dof["default_degrees"] + 45.0
        update_response = client.post(f"/api/sidechains/update/{ws_id}", json={
            "sidechain_data": gui_payload,
            "changed_dof_keys": [dof["dof_key"]],
        })
        assert update_response.status_code == 200, update_response.text

        response = client.post(f"/api/sidechains/optimize/{ws_id}", json={
            "sidechain_data": gui_payload,
            "residue": residue,
        })

        assert response.status_code == 200, response.text
        body = response.json()
        assert len(body["dofs"]) == 3
        assert len(body["updated_atoms"]) > 0

    def test_update_without_active_session_returns_404(self, client, offline_forge_ff, make_workspace, pdb_glu_open_branch):
        ws_id = make_workspace(pdb_glu_open_branch)

        response = client.post(f"/api/sidechains/update/{ws_id}", json={
            "sidechain_data": {"sidechains": []},
            "changed_dof_keys": [],
        })

        assert response.status_code == 404
        assert response.json()["code"] == "not_found"


class TestSidechainCommit:
    def test_commit_finishes_structure_and_clears_session(
        self, client, offline_forge_ff, make_workspace, pdb_glu_open_branch
    ):
        ws_id = make_workspace(pdb_glu_open_branch)
        start_response = client.post(f"/api/sidechains/start/{ws_id}", json=_start_payload(ws_id))
        assert start_response.status_code == 200, start_response.text

        response = client.post(f"/api/sidechains/commit/{ws_id}")

        assert response.status_code == 200, response.text
        body = response.json()
        assert body["message"] == "Structure successfully prepared."

        updated = workspace_manager.get_file_path(ws_id, "structure.pdb").read_text()
        updated_atom_count = sum(1 for l in updated.splitlines() if l.startswith("ATOM"))
        assert updated_atom_count > 5  # side-chain + vodíky doplněny

        assert not workspace_manager.get_file_path(ws_id, "structure_preview.pdb").exists()

        # Relace je po commitu pryč - další commit musí vrátit 404, ne omylem zopakovat práci.
        second_commit = client.post(f"/api/sidechains/commit/{ws_id}")
        assert second_commit.status_code == 404


class TestSidechainCancel:
    def test_cancel_leaves_canonical_structure_untouched(
        self, client, offline_forge_ff, make_workspace, pdb_glu_open_branch
    ):
        ws_id = make_workspace(pdb_glu_open_branch)
        original = workspace_manager.get_file_path(ws_id, "structure.pdb").read_text()
        start_response = client.post(f"/api/sidechains/start/{ws_id}", json=_start_payload(ws_id))
        assert start_response.status_code == 200, start_response.text

        response = client.post(f"/api/sidechains/cancel/{ws_id}")
        assert response.status_code == 200, response.text

        unchanged = workspace_manager.get_file_path(ws_id, "structure.pdb").read_text()
        assert unchanged == original

        # relace je pryč - navazující update musí selhat na 404, ne tiše dojet na starém stavu
        after_cancel = client.post(f"/api/sidechains/commit/{ws_id}")
        assert after_cancel.status_code == 404


class TestSidechainRealGapBoundaryRegression:
    """
    Regresní test pro skutečný bug: analysis_service dřív kaskádovitě mazala
    z modelu KAŽDÉ reziduum na okraji gapu, které mělo jakýkoliv chybějící
    těžký atom - i když mělo svou vlastní backbone kotvu (C) a bylo tedy
    přesně tím "residue_local_open_branch" případem, který nový builder umí
    interaktivně dostavit. Na reálném 1JJ2 to znamenalo, že GLU83 (chain K)
    zmizelo z modelu ještě před buildem a /sidechains/start rovnou vrátilo
    "complete" bez jakékoliv GUI relace - přesně to, co uživatel nahlásil.
    """

    def test_glu83_survives_gap_boundary_exclusion_and_opens_gui(
        self, client, offline_forge_ff, make_workspace, pdb_1jj2_glu_gap_boundary
    ):
        ws_id = make_workspace(pdb_1jj2_glu_gap_boundary)

        response = client.post(f"/api/sidechains/start/{ws_id}", json=_start_payload(ws_id))

        assert response.status_code == 200, response.text
        body = response.json()
        assert body["status"] == "missing_dof"

        sidechains = body["gui_payload"]["sidechains"]
        assert len(sidechains) == 1
        residue = sidechains[0]["residue"]
        assert residue["chain_id"] == "K"
        assert residue["resseq"] == 83
        assert residue["ff_resname"] == "CGLU"
        dof_names = sorted(dof["dof_key"]["atom"] for dof in sidechains[0]["dofs"])
        assert dof_names == ["CD", "CG", "OE1"]
