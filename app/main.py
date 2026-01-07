# app/main.py
from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from app.api.api import api_router
from app.core.logging import setup_logging
from app.core.exceptions import AppBaseException, app_exception_handler
from app.core.config import settings

setup_logging()

app = FastAPI(
    title=settings.PROJECT_NAME,
    version=settings.VERSION
)

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
    return {"status": "online", "version": settings.VERSION}