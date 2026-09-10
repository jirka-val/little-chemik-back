"""
Integrační test pro propojení /api/validation/prepare -> structure.forge_meta.json
-> TopologyService._load_forge_meta (Fáze 6).

Záměrně NEspouští celou generate_topology()/AMBER_topology pipeline - ta by si
žádala reálná silová pole ve starém ff_cache layoutu a je to samostatná,
mnohem těžší závislost. Tenhle test ověřuje přesně tu spojku, kterou jsme
přidali: že TopologyService umí přečíst sidecar, který ForgeStructureService
skutečně zapsal, a že bez sidecaru se nic nerozbije.
"""

import pytest

from app.services.topology_service import TopologyService
from app.workspaces.manager import workspace_manager

pytestmark = pytest.mark.integration


class TestSidecarRoundTrip:
    def test_topology_service_reads_sidecar_written_by_prepare(
        self, client, offline_forge_ff, make_workspace, pdb_rna_gap
    ):
        ws_id = make_workspace(pdb_rna_gap)
        response = client.post("/api/validation/prepare", json={
            "workspace_id": ws_id,
            "ff_selections": {"R": {"display_name": "OL3"}},
            "ph": 7.0,
            "add_solvent": False,
        })
        assert response.status_code == 200, response.text

        workspace_dir = workspace_manager.get_workspace_dir(ws_id)
        meta = TopologyService()._load_forge_meta(workspace_dir, "structure.pdb")

        assert meta is not None
        assert meta["A:3:"]["ff_resname"] == "RA3"
        assert meta["A:6:"]["ff_resname"] == "RU5"

    def test_load_forge_meta_returns_none_when_sidecar_missing(self, make_workspace, pdb_alanine_single):
        """Workspace, ktery nikdy neprosel /prepare - musi bezpecne spadnout na None, ne vyjimku."""
        ws_id = make_workspace(pdb_alanine_single)
        workspace_dir = workspace_manager.get_workspace_dir(ws_id)

        meta = TopologyService()._load_forge_meta(workspace_dir, "structure.pdb")
        assert meta is None

    def test_load_forge_meta_survives_corrupted_sidecar(self, make_workspace, pdb_alanine_single):
        ws_id = make_workspace(pdb_alanine_single)
        workspace_dir = workspace_manager.get_workspace_dir(ws_id)
        (workspace_dir / "structure.forge_meta.json").write_text("{not valid json", encoding="utf-8")

        meta = TopologyService()._load_forge_meta(workspace_dir, "structure.pdb")
        assert meta is None
