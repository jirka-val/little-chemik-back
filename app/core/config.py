from pydantic import field_validator
from pydantic_settings import BaseSettings, SettingsConfigDict
from pathlib import Path
from typing import List


class Settings(BaseSettings):
    PROJECT_NAME: str = "Little Chemik API"
    VERSION: str = "1.2.0"
    API_V1_STR: str = "/api/v1"

    # Kdo smí volat API z prohlížeče z JINÉHO originu (CORS). V produkci
    # frontend i API běží pod stejnou adresou (nginx proxy /api), takže CORS
    # není potřeba vůbec a nezáleží na tom, na jaké doméně/IP app běží.
    # Lokální vývoj jde přes Vite proxy (vite.config.ts), taky same-origin;
    # localhost porty tu zůstávají pro případ, kdy frontend volá backend
    # napřímo (VITE_API_URL=http://localhost:8000). Další originy přes .env:
    # BACKEND_CORS_ORIGINS=https://a.cz,https://b.cz
    BACKEND_CORS_ORIGINS: List[str] = [
        "http://localhost:5173",   # Standardní Vite port
        "http://localhost:5174",   # Alternativní Vite port
    ]

    @field_validator("BACKEND_CORS_ORIGINS", mode="before")
    @classmethod
    def split_cors_origins(cls, v):
        if isinstance(v, str) and not v.startswith("["):
            return [origin.strip() for origin in v.split(",") if origin.strip()]
        return v

    BASE_DIR: Path = Path(__file__).resolve().parent.parent.parent
    PDB_DATA_DIR: Path = BASE_DIR / "data" / "pdb_files"

    # Pracovní složky uživatelských relací (maže je garbage collector).
    WORKSPACE_DIR: Path = BASE_DIR / "temp_workspaces"
    # Soubory silových polí rozbalené z katalogu (viz ForceFieldService).
    FF_CACHE_DIR: Path = BASE_DIR / "data" / "ff_cache"
    FF_FORGE_CACHE_DIR: Path = BASE_DIR / "data" / "ff_cache_forge"

    RCSB_PDB_URL: str = "https://files.rcsb.org/download"

    FORCE_FIELDS_CLASSIFICATION_FILE: Path = BASE_DIR / "data" / "force_fields.json"
    FF_CATALOG_SNAPSHOT_FILE: Path = BASE_DIR / "data" / "ff_catalog.json"
    # Jak často (v sekundách) se má katalog FF automaticky obnovovat z IDA na
    # pozadí (viz app/workspaces/tasks/ff_catalog_refresher.py). Výchozí 24h
    # odpovídá Pavlovu "aktualizace jednou denně v noci".
    FF_CATALOG_REFRESH_INTERVAL_SECONDS: int = 24 * 60 * 60
    # Sdílený token pro admin operace (přeřazování FF mezi tiery). Prázdné =
    # endpoint je zamčený úplně, dokud si ho nasazení nenastaví v .env -
    # bezpečnější výchozí stav než "otevřeno pro každého".
    ADMIN_TOKEN: str = ""

    # Hlášení chyb (viz app/services/incidents/). Vypnutí přes .env
    # INCIDENT_REPORTS_ENABLED=false vypne deník akcí, snímky PDB i příjem
    # reportů. Reporty se ukládají jen se souhlasem uživatele, mimo
    # temp_workspaces (ty maže garbage collector) - v Dockeru musí být
    # INCIDENTS_DIR připojený jako volume, jinak zmizí s každým nasazením.
    INCIDENT_REPORTS_ENABLED: bool = True
    INCIDENTS_DIR: Path = BASE_DIR / "data" / "incidents"
    # Strop pro všechny reporty dohromady; po jeho dosažení se nové reporty
    # odmítají (nic se nemaže automaticky).
    INCIDENTS_MAX_TOTAL_MB: int = 5000
    # Strop pro jeden report - snímky PDB nad něj se vynechají (od nejstarších).
    INCIDENT_MAX_MB: int = 300

    # extra="ignore": klíč v .env, který Settings nezná (např. FORGE_URL pro
    # scripts/incidents.py), nesmí shodit start backendu.
    model_config = SettingsConfigDict(env_file=".env", case_sensitive=True, extra="ignore")


settings = Settings()