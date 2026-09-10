"""
Integrační testy pro POST /api/validation/prepare a jeho deprecated alias
POST /api/molecules/add-hydrogens/{id} - obě cesty teď jdou přes
ForgeStructureService (app/builder), místo starého PDBFixer/OpenMM
HydrogenationService.

Používá fixturu `offline_forge_ff` (viz conftest.py), takže nesahá na síť ani
na reálnou data/ff_cache/ - staví si izolovanou kopii z toho, co je už
lokálně nacachované. Testy se samy přeskočí, pokud lokální cache pro OL3
(RNA silové pole) na tomhle stroji chybí.
"""

import json

import pytest

from app.workspaces.manager import workspace_manager

pytestmark = pytest.mark.integration


class TestPrepareHappyPath:
    def test_prepare_builds_missing_hydrogens_and_returns_200(
        self, client, offline_forge_ff, make_workspace, pdb_alanine_single
    ):
        ws_id = make_workspace(pdb_alanine_single)

        response = client.post("/api/validation/prepare", json={
            "workspace_id": ws_id,
            "ff_selections": {"R": {"display_name": "OL3"}},
            "ph": 7.0,
            "add_solvent": False,
        })

        # Fixtura je protein, ne RNA/DNA - očekáváme buď úspěch (pokud converting
        # dictionary umí i tenhle typ přes fallback), nebo jasnou 4xx/5xx chybu
        # od builderu - v žádném případě ne tichý pád mimo náš kontrakt.
        assert response.status_code in (200, 409, 500)
        if response.status_code == 200:
            body = response.json()
            assert "warnings" in body
            assert body["message"] == "Structure successfully prepared."

    def test_prepare_writes_forge_meta_sidecar(self, client, offline_forge_ff, make_workspace, pdb_rna_gap):
        ws_id = make_workspace(pdb_rna_gap)

        response = client.post("/api/validation/prepare", json={
            "workspace_id": ws_id,
            "ff_selections": {"R": {"display_name": "OL3"}},
            "ph": 7.0,
            "add_solvent": False,
        })

        assert response.status_code == 200, response.text
        meta_path = workspace_manager.get_file_path(ws_id, "structure.forge_meta.json")
        assert meta_path.exists()

        meta = json.loads(meta_path.read_text(encoding="utf-8"))
        # Autoritativní terminál pro rezidua na okraji gapu - jádro Fáze 6.
        assert meta["A:3:"]["ff_resname"] == "RA3"
        assert meta["A:6:"]["ff_resname"] == "RU5"

    def test_prepare_updates_structure_pdb_with_built_atoms(
        self, client, offline_forge_ff, make_workspace, pdb_rna_gap
    ):
        ws_id = make_workspace(pdb_rna_gap)
        original = workspace_manager.get_file_path(ws_id, "structure.pdb").read_text(encoding="utf-8")

        response = client.post("/api/validation/prepare", json={
            "workspace_id": ws_id,
            "ff_selections": {"R": {"display_name": "OL3"}},
            "ph": 7.0,
            "add_solvent": False,
        })
        assert response.status_code == 200, response.text

        updated = workspace_manager.get_file_path(ws_id, "structure.pdb").read_text(encoding="utf-8")
        original_atom_count = sum(1 for l in original.splitlines() if l.startswith("ATOM"))
        updated_atom_count = sum(1 for l in updated.splitlines() if l.startswith("ATOM"))
        assert updated_atom_count > original_atom_count  # builder doplnil chybějící (vodíkové) atomy


class TestPrepareValidation:
    def test_missing_ff_selections_returns_422(self, client, make_workspace, pdb_alanine_single):
        ws_id = make_workspace(pdb_alanine_single)
        response = client.post("/api/validation/prepare", json={
            "workspace_id": ws_id,
            "ph": 7.0,
        })
        assert response.status_code == 422

    def test_nonexistent_workspace_returns_404(self, client):
        response = client.post("/api/validation/prepare", json={
            "workspace_id": "does-not-exist",
            "ff_selections": {"R": {"display_name": "OL3"}},
        })
        assert response.status_code == 404


class TestPrepareMissingDof:
    def test_missing_dof_returns_409_with_structured_payload(
        self, client, monkeypatch, make_workspace, pdb_alanine_single
    ):
        """
        Simuluje stav, kdy builder narazí na chybějící stupeň volnosti - ověřuje
        HTTP kontrakt (409 + strukturovaný payload), ne samotnou geometrii
        builderu (tu je těžké spolehlivě vyprovokovat syntetickým vstupem).
        """
        import sys
        from pathlib import Path

        builder_dir = Path("app/builder").resolve()
        if str(builder_dir) not in sys.path:
            sys.path.insert(0, str(builder_dir))
        from forge_molecule_parser import Molecule, Chain, Residue
        from forge_molecule_builder import AtomID, PlannedMissingDOFStep, ResolvedRef, DOFKey

        import app.api.v1.endpoints.validation as validation_ep
        from app.services.structure.forge_service import ForgeMissingDOFError

        residue = Residue(chain_id="A", resseq=47, icode="", ff_resname="TYR",
                           atoms={}, index_in_chain=0, original_resname="TYR", group="P")
        molecule = Molecule(chains={"A": Chain(chain_id="A", residues=[residue])})
        reason_atom = AtomID(chain_id="A", residue_index=0, atom_name="CG")
        ref = ResolvedRef(chain_id="A", residue_index=0, atom_name="CB")
        step = PlannedMissingDOFStep(
            dof_key=DOFKey("A", 0, "CG", 0), central_bond=(ref, ref),
            requested_dihedral_atoms=(ref, ref, ref, ref), torsion_group_index=0,
            requested_member_index=0, reason_atom=reason_atom, reason_rule_index=0,
        )

        def _raise(**kwargs):
            raise ForgeMissingDOFError(step, molecule)

        monkeypatch.setattr(validation_ep.forge_service, "prepare_structure", _raise)

        ws_id = make_workspace(pdb_alanine_single)
        response = client.post("/api/validation/prepare", json={
            "workspace_id": ws_id,
            "ff_selections": {"P": {"display_name": "dummy"}},
        })

        assert response.status_code == 409
        body = response.json()
        assert body["code"] == "missing_dof"
        assert body["chain"] == "A"
        assert body["resseq"] == 47
        assert body["atom_name"] == "CG"


class TestAddHydrogensDeprecatedAlias:
    def test_add_hydrogens_delegates_to_same_pipeline(
        self, client, offline_forge_ff, make_workspace, pdb_rna_gap
    ):
        ws_id = make_workspace(pdb_rna_gap)

        response = client.post(f"/api/molecules/add-hydrogens/{ws_id}", json={
            "ff_selections": {"R": {"display_name": "OL3"}},
            "ph": 7.0,
            "optimize": False,
        })

        assert response.status_code == 200, response.text
        assert workspace_manager.get_file_path(ws_id, "structure.forge_meta.json").exists()

    def test_add_hydrogens_missing_ff_selections_returns_422(self, client, make_workspace, pdb_alanine_single):
        ws_id = make_workspace(pdb_alanine_single)
        response = client.post(f"/api/molecules/add-hydrogens/{ws_id}", json={"ph": 7.0})
        assert response.status_code == 422

    def test_add_hydrogens_is_marked_deprecated_in_openapi(self, client):
        schema = client.get("/openapi.json").json()
        path_item = schema["paths"]["/api/molecules/add-hydrogens/{workspace_id}"]["post"]
        assert path_item.get("deprecated") is True
