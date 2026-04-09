import asyncio  # <-- NOVÉ: Knihovna pro asynchronní úlohy na pozadí
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.api.api import api_router
from app.core.logging import setup_logging
from app.core.exceptions import AppBaseException, app_exception_handler
from app.core.config import settings

# <-- NOVÉ: Import naší vytvořené uklízečky
from app.workspaces.tasks.garbage_collector import cleanup_old_workspaces

setup_logging()

app = FastAPI(
    title=settings.PROJECT_NAME,
    version=settings.VERSION
)

app.add_middleware(
    CORSMiddleware,
    allow_origins=[
        "http://147.251.115.223",  # Produkce
        "http://localhost:5173",   # Standardní Vite port
        "http://localhost:5174",   # Tvůj aktuální Vite port
    ],
    allow_credentials=True,
    allow_methods=["*"],
    allow_headers=["*"],
)

# ... horní část kódu zůstává stejná (importy, middleware, atd.) ...

app.add_exception_handler(AppBaseException, app_exception_handler)

app.include_router(api_router, prefix="/api")

# Vytvoříme si globální proměnnou, abychom si pamatovali náš proces
cleanup_task = None

@app.on_event("startup")
async def startup_event():
    global cleanup_task
    # Spustíme naši uklízečku a uložíme si ji do proměnné
    cleanup_task = asyncio.create_task(cleanup_old_workspaces())

# <-- NOVÉ: Přidáme událost vypnutí serveru
@app.on_event("shutdown")
async def shutdown_event():
    global cleanup_task
    # Když zmáčknete Ctrl+C, server tuto smyčku bezpečně odstřelí
    if cleanup_task:
        cleanup_task.cancel()

@app.get("/", tags=["Health Check"])
async def root():
    return {"status": "online", "version": settings.VERSION}