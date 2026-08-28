import asyncio  # <-- NOVÉ: Knihovna pro asynchronní úlohy na pozadí
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.api.api import api_router
from app.core.logging import setup_logging
from app.core.exceptions import AppBaseException, app_exception_handler
from app.core.config import settings
from app.core.http_client import close_external_http_client

# <-- NOVÉ: Import naší vytvořené uklízečky
from app.workspaces.tasks.garbage_collector import cleanup_old_workspaces
from app.workspaces.tasks.ff_catalog_refresher import refresh_ff_catalog_periodically
from app.services.ff_catalog_service import catalog_service

setup_logging()

app = FastAPI(
    title=settings.PROJECT_NAME,
    version=settings.VERSION
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=settings.BACKEND_CORS_ORIGINS,
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ... horní část kódu zůstává stejná (importy, middleware, atd.) ...

app.add_exception_handler(AppBaseException, app_exception_handler)

app.include_router(api_router, prefix="/api")

# Vytvoříme si globální proměnné, abychom si pamatovali naše background procesy
cleanup_task = None
ff_catalog_task = None

@app.on_event("startup")
async def startup_event():
    global cleanup_task, ff_catalog_task
    # Spustíme naši uklízečku a uložíme si ji do proměnné
    cleanup_task = asyncio.create_task(cleanup_old_workspaces())

    # Bootstrap FF katalogu - pokud po čerstvém deployi ještě neexistuje
    # žádný lokální snapshot, uděláme jeden synchronní refresh, ať FF panel
    # hned po startu nevrátí prázdný seznam. Dál se stará noční background job.
    await asyncio.to_thread(catalog_service.ensure_catalog)
    ff_catalog_task = asyncio.create_task(refresh_ff_catalog_periodically())

# <-- NOVÉ: Přidáme událost vypnutí serveru
@app.on_event("shutdown")
async def shutdown_event():
    global cleanup_task, ff_catalog_task
    # Když zmáčknete Ctrl+C, server tyto smyčky bezpečně odstřelí
    if cleanup_task:
        cleanup_task.cancel()
    if ff_catalog_task:
        ff_catalog_task.cancel()
    await close_external_http_client()

@app.get("/", tags=["Health Check"])
async def root():
    return {"status": "online", "version": settings.VERSION}