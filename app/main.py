from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.api.api import api_router
from app.core.logging import setup_logging
from app.core.exceptions import AppBaseException, app_exception_handler

setup_logging()

app = FastAPI(
    title="Little Chemik API",
    description="Profesionální backend pro analýzu nukleových kyselin.",
    version="1.2.0"
)

# CORS nastavení
app.add_middleware(
    CORSMiddleware,
    allow_origins=["*"],
    allow_methods=["*"],
    allow_headers=["*"],
)

app.add_exception_handler(AppBaseException, app_exception_handler)

app.include_router(api_router)

@app.get("/", tags=["Health Check"])
async def root():
    return {"status": "online", "version": "1.2.0"}