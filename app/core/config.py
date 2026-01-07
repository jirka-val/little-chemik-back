from pydantic_settings import BaseSettings
from pathlib import Path
from typing import List


class Settings(BaseSettings):
    # Základní metadata aplikace
    PROJECT_NAME: str = "Little Chemik API"
    VERSION: str = "1.2.0"
    API_V1_STR: str = "/api/v1"

    # CORS
    BACKEND_CORS_ORIGINS: List[str] = ["*"]

    # Cesty k datům
    BASE_DIR: Path = Path(__file__).resolve().parent.parent.parent
    PDB_DATA_DIR: Path = BASE_DIR / "data" / "pdb_files"

    # Externí zdroje
    RCSB_PDB_URL: str = "https://files.rcsb.org/download"

    class Config:
        env_file = ".env"
        case_sensitive = True


settings = Settings()