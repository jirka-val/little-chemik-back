import asyncio
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.api.api import api_router
from app.core.logging import setup_logging
from app.core.exceptions import AppBaseException, app_exception_handler
from app.core.config import settings
from app.core.http_client import close_external_http_client
from app.core.incident_middleware import IncidentJournalMiddleware
from app.core.console_middleware import ConsoleWorkspaceMiddleware

from app.workspaces.tasks.garbage_collector import cleanup_old_workspaces
from app.workspaces.tasks.ff_catalog_refresher import refresh_ff_catalog_periodically
from app.services.ff_catalog_service import catalog_service

setup_logging()


@asynccontextmanager
async def lifespan(app: FastAPI):
    # Úlohy na pozadí: úklid starých workspace a noční refresh FF katalogu.
    cleanup_task = asyncio.create_task(cleanup_old_workspaces())

    # Bootstrap FF katalogu - pokud po čerstvém deployi ještě neexistuje
    # žádný lokální snapshot, uděláme jeden synchronní refresh, ať FF panel
    # hned po startu nevrátí prázdný seznam. Dál se stará noční background job.
    await asyncio.to_thread(catalog_service.ensure_catalog)
    ff_catalog_task = asyncio.create_task(refresh_ff_catalog_periodically())

    yield

    # Vypnutí serveru (Ctrl+C, docker stop): smyčky na pozadí bezpečně ukončit.
    cleanup_task.cancel()
    ff_catalog_task.cancel()
    await close_external_http_client()


app = FastAPI(
    title=settings.PROJECT_NAME,
    version=settings.VERSION,
    lifespan=lifespan,
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.BACKEND_CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
    # Frontend čte jméno staženého souboru z Content-Disposition (mdin, export).
    # Cross-origin (dev: :5173 -> :8000) ho prohlížeč bez expose neukáže a
    # frontend pak spadl na výchozí "production.mdin".
    expose_headers=["Content-Disposition"],
)

# Deník akcí a snímky struktury pro hlášení chyb (INCIDENT_REPORTS_ENABLED).
app.add_middleware(IncidentJournalMiddleware)

# Zprávy Console panelu vidí jen workspace, který je vyvolal.
app.add_middleware(ConsoleWorkspaceMiddleware)

app.add_exception_handler(AppBaseException, app_exception_handler)

app.include_router(api_router, prefix="/api")

@app.get("/", tags=["Health Check"])
async def root():
    return {"status": "online", "version": settings.VERSION}


# Stejná odpověď pod /api - jen /api/* jde přes nginx proxy na backend, takže
# tohle je adresa pro healthcheck zvenku (docker, monitoring, proxy).
@app.get("/api/health", tags=["Health Check"])
async def api_health():
    return await root()