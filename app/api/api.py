from fastapi import APIRouter
from app.api.v1.endpoints import molecules

api_router = APIRouter()
api_router.include_router(molecules.router, prefix="/molecules", tags=["molecules"])