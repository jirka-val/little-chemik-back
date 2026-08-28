# app/core/http_client.py
#
# Sdílený httpx.AsyncClient pro odchozí volání na externí služby (RCSB).
# Bez tohohle by si každý request vytvářel vlastní klienta a s ním nové
# TCP/TLS spojení místo využití keep-alive connection poolu.
import httpx

external_http_client = httpx.AsyncClient()


async def close_external_http_client():
    await external_http_client.aclose()
