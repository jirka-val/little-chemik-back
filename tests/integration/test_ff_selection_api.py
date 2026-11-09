"""
FF endpoint posílá FF bez obsahů souborů; příprava a topologie přijmou jen FF
z katalogu (app/services/ff_selection.py).
"""

import json

import pytest

from app.services.ff_selection import FF_FILE_KEYS

pytestmark = pytest.mark.integration


def _all_ff_objects(response_json):
    yield from response_json["forcefields"]
    for group in response_json["forcefields_by_group"].values():
        yield from group
    for profile in response_json["water_profiles"]:
        yield profile["water_ff"]
        yield from profile["ion_ffs"].values()


def test_forcefield_list_has_no_file_contents(client, make_workspace, pdb_alanine_single):
    ws_id = make_workspace(pdb_alanine_single)
    response = client.get(f"/api/forcefields/{ws_id}", params={"mode": "expert"})
    assert response.status_code == 200, response.text
    for ff in _all_ff_objects(response.json()):
        assert not set(FF_FILE_KEYS) & set(ff), ff.get("display_name")
    assert len(response.content) < 1_000_000


@pytest.mark.parametrize("path", ["/api/sidechains/start/{ws}", "/api/topology/{ws}/generate"])
def test_force_field_outside_the_catalog_is_rejected(client, make_workspace, pdb_alanine_single, path):
    ws_id = make_workspace(pdb_alanine_single)
    body = {"workspace_id": ws_id, "ph": 7.0, "add_solvent": False,
            "ff_selections": {"P": {"id": -1, "display_name": "../../evil", "bonded_ff_file": "INJECTED"}}}
    response = client.post(path.format(ws=ws_id), json=body)
    assert response.status_code == 400, response.text
    assert response.json()["code"] == "unknown_force_field"
    assert "INJECTED" not in json.dumps(response.json())


def test_main_file_download_comes_from_the_catalog(client):
    from app.services.ff_catalog_service import catalog_service
    catalog = catalog_service.get_forcefields()
    if not catalog:
        pytest.skip("FF catalog (data/ff_catalog.json) not downloaded")
    ff = next(f for f in catalog if f.get("force_field_file"))
    response = client.get(f"/api/forcefields/catalog/{ff['id']}/main-file")
    assert response.status_code == 200
    assert "attachment" in response.headers["content-disposition"]
    assert response.text.strip()
    assert client.get("/api/forcefields/catalog/-5/main-file").status_code == 404
