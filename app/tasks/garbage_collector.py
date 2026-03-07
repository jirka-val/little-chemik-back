import os
import time
import asyncio
import logging
from app.workspaces.manager import WORKSPACE_DIR

logger = logging.getLogger("api")

# Nastavení času: Jak dlouho může soubor žít v sekundách (2 hodiny = 7200 sekund)
MAX_AGE_SECONDS = 2 * 60 * 60

# Nastavení frekvence úklidu: Jak často se má uklízečka probudit (1 hodina = 3600 sekund)
CLEANUP_INTERVAL = 3600


async def cleanup_old_workspaces():
    """
    Nekonečná smyčka, která běží na pozadí po celou dobu životnosti serveru.
    Promazává staré a opuštěné PDB soubory.
    """
    logger.info("Garbage Collector byl úspěšně spuštěn na pozadí.")

    while True:
        try:
            now = time.time()
            deleted_count = 0

            # Zkontrolujeme, zda složka vůbec existuje
            if os.path.exists(WORKSPACE_DIR):
                for filename in os.listdir(WORKSPACE_DIR):
                    file_path = os.path.join(WORKSPACE_DIR, filename)

                    # Ignorujeme složky, zajímají nás jen soubory
                    if os.path.isfile(file_path):
                        # Zjistíme, kdy byl soubor naposledy upraven (nebo vytvořen)
                        file_age = now - os.path.getmtime(file_path)

                        # Pokud je starší než náš limit, smažeme ho
                        if file_age > MAX_AGE_SECONDS:
                            os.remove(file_path)
                            deleted_count += 1

            if deleted_count > 0:
                logger.info(f"Garbage Collector úřadoval: Smazáno {deleted_count} starých workspace souborů.")

        except Exception as e:
            logger.error(f"Chyba při běhu Garbage Collectoru: {e}")

        # Uspíme proces, aby nežral výkon procesoru. Vzbudí se zase za hodinu.
        await asyncio.sleep(CLEANUP_INTERVAL)