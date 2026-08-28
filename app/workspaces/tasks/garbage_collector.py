import os
import time
import asyncio
import logging
import shutil  # Nutné pro mazání celých složek
from app.workspaces.manager import WORKSPACE_DIR

logger = logging.getLogger("api")

# Jak dlouho může složka žít (2 hodiny)
MAX_AGE_SECONDS = 2 * 60 * 60

# Jak často se má proces spouštět (každou hodinu)
CLEANUP_INTERVAL = 3600


def _cleanup_old_workspaces_sync() -> int:
    """
    Synchronní tělo úklidu (listdir/getmtime/rmtree jsou blokující syscally).
    Voláno přes asyncio.to_thread, ať po dobu mazání velkých/starých
    workspace adresářů nezůstane stát event loop pro všechny ostatní
    současně obsluhované requesty.
    """
    now = time.time()
    deleted_count = 0

    if os.path.exists(WORKSPACE_DIR):
        # Procházíme vše v hlavní složce temp_workspaces
        for item in os.listdir(WORKSPACE_DIR):
            item_path = os.path.join(WORKSPACE_DIR, item)

            # Zjistíme čas poslední změny složky/souboru
            item_age = now - os.path.getmtime(item_path)

            if item_age > MAX_AGE_SECONDS:
                try:
                    if os.path.isdir(item_path):
                        # Smaže složku i se všemi soubory uvnitř (topologie, pdb...)
                        shutil.rmtree(item_path)
                    else:
                        # Pro jistotu, kdyby tam zůstal zapomenutý soubor mimo složku
                        os.remove(item_path)

                    deleted_count += 1
                except Exception as sub_e:
                    logger.error(f"Failed to delete {item_path}: {sub_e}")

    return deleted_count


async def cleanup_old_workspaces():
    """
    Nekonečná smyčka pro promazávání starých a opuštěných adresářů ve workspace.
    Teď už korektně maže celé složky (UUID složky).
    """
    logger.info("Garbage Collector has been successfully started in the background.")

    while True:
        try:
            deleted_count = await asyncio.to_thread(_cleanup_old_workspaces_sync)

            if deleted_count > 0:
                logger.info(f"Garbage Collector cleanup finished: Removed {deleted_count} old workspace directories.")

        except Exception as e:
            logger.error(f"Error during Garbage Collector run: {e}")

        # Čekání na další interval
        await asyncio.sleep(CLEANUP_INTERVAL)