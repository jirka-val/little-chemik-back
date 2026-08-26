import asyncio
import logging

from app.core.config import settings
from app.services.ff_catalog_service import catalog_service

logger = logging.getLogger("api")


async def refresh_ff_catalog_periodically():
    """
    Nekonečná smyčka, stejný tvar jako garbage_collector.cleanup_old_workspaces:
    jednou za FF_CATALOG_REFRESH_INTERVAL_SECONDS (výchozí 24h) stáhne čerstvý
    seznam FF z IDA a uloží ho jako lokální snapshot (FFCatalogService).
    Odpovídá Pavlovu "aktualizaci by dělal např. každý den v noci" - přesný
    čas dne se neřeší, důležitá je pravidelnost bez ručního zásahu.
    """
    logger.info("FF catalog background refresher has been started.")

    while True:
        try:
            snapshot = await asyncio.to_thread(catalog_service.refresh_catalog)
            logger.info(f"FF catalog auto-refresh finished: {len(snapshot['forcefields'])} force field(s).")
        except Exception as e:
            logger.error(f"FF catalog auto-refresh failed: {e}")

        await asyncio.sleep(settings.FF_CATALOG_REFRESH_INTERVAL_SECONDS)
