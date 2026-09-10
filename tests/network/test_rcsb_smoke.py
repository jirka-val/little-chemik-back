"""
Síťové testy - skutečné stahování z RCSB. Ve výchozím běhu `pytest` se
přeskakují (viz pytest.ini `addopts`); zapneš je přes `pytest -m network`.

Migrace síťových částí ze starého tests/test_api.py.
"""

import asyncio
import time

import pytest
from httpx import AsyncClient, ASGITransport

from app.main import app

pytestmark = pytest.mark.network


class TestFetchFromRcsb:
    def test_fetch_valid_pdb_code(self, client):
        response = client.get("/api/molecules/fetch-pdb/1crn")
        assert response.status_code == 200
        assert "workspace_id" in response.json()

    def test_fetch_invalid_pdb_code_returns_404(self, client):
        response = client.get("/api/molecules/fetch-pdb/NON_EXISTENT_9999")
        assert response.status_code == 404


class TestRemoteAnalysis:
    def test_3dvz_missing_atoms_detected(self, client):
        ws_id = client.get("/api/molecules/fetch-pdb/3dvz").json()["workspace_id"]
        response = client.get(f"/api/analysis/sequence/{ws_id}")
        assert response.status_code == 200

        tokens = response.json()["sequence"]["chains"]["A"]["tokens"]
        u2647 = next((t for t in tokens if str(t.get("resseq")) == "2647"), None)
        assert u2647 is not None
        assert len(u2647["missing_atoms"]) > 0

    def test_3dvz_rna_variants_identified(self, client):
        ws_id = client.get("/api/molecules/fetch-pdb/3dvz").json()["workspace_id"]
        response = client.get(f"/api/analysis/sequence/{ws_id}")
        tokens = response.json()["sequence"]["chains"]["A"]["tokens"]

        first = next(t for t in tokens if str(t.get("resseq")) == "2647")
        assert first["ff_resname"] == "RU5"
        assert first["known"] is True

        internal = next(t for t in tokens if str(t.get("resseq")) == "2648")
        assert internal["ff_resname"] == "RG"
        assert internal["known"] is True

    def test_analyze_remote_pdb_endpoint_auto_cleans_structure(self, client):
        """GET /api/analysis/analyze-pdb/{code} - stáhne, vyčistí a vrátí i pdb_text."""
        response = client.get("/api/analysis/analyze-pdb/1RNA")
        assert response.status_code == 200
        body = response.json()
        assert "pdb_text" in body and len(body["pdb_text"]) > 0
        assert "missing_atoms" in body


@pytest.mark.slow
class TestConcurrency:
    def test_heavy_cpu_request_does_not_block_light_request(self, client):
        """
        Ověřuje, že těžký CPU-bound request (sekvenční analýza) neblokuje event
        loop natolik, aby lehký request (download) musel čekat na jeho dokončení.
        """
        ws_id = client.get("/api/molecules/fetch-pdb/1crn").json()["workspace_id"]

        async def _run():
            async with AsyncClient(transport=ASGITransport(app=app), base_url="http://test") as ac:
                async def fetch_heavy():
                    start = time.time()
                    res = await ac.get(f"/api/analysis/sequence/{ws_id}")
                    return time.time() - start, res.status_code

                async def fetch_light():
                    await asyncio.sleep(0.05)
                    start = time.time()
                    res = await ac.get(f"/api/download/{ws_id}")
                    return time.time() - start, res.status_code

                return await asyncio.gather(fetch_heavy(), fetch_light())

        (heavy_duration, heavy_status), (light_duration, light_status) = asyncio.run(_run())

        assert heavy_status == 200
        assert light_status == 200
        assert light_duration < heavy_duration or light_duration < 0.1
