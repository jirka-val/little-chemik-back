import logging
import sys
import threading
from collections import deque
from typing import Any, Dict, List


class InMemoryLogHandler(logging.Handler):
    """
    Drží posledních `capacity` log záznamů v paměti procesu, ať je frontend
    (Console panel v sidebaru) může pollovat přes /api/system/logs.

    Připojený POUZE na `console_logger` (viz níže), ne na root logger - ten
    zůstává jen ve stdoutu/terminálu. Bez tohohle oddělení by se do panelu
    (viditelného v prohlížeči) dostalo úplně všechno včetně workspace ID a
    dalších detailů z běžných per-request logů, což už jednou byl reálný
    problém (unik citlivých identifikátorů + zahlcení nedůležitými zprávami).
    Sem smí jen krátké, obecné stavové zprávy bez identifikátorů.

    Každý záznam dostane rostoucí `id`, takže frontend může pollovat
    přírůstkově (`since_id`) místo opakovaného posílání celého bufferu.
    """

    def __init__(self, capacity: int = 100):
        super().__init__()
        self._buffer: deque[Dict[str, Any]] = deque(maxlen=capacity)
        self._lock = threading.Lock()
        self._next_id = 1

    def emit(self, record: logging.LogRecord) -> None:
        try:
            message = record.getMessage()
        except Exception:
            message = str(record.msg)

        with self._lock:
            entry = {
                "id": self._next_id,
                "timestamp": record.created,
                "level": record.levelname,
                "message": message,
            }
            self._next_id += 1
            self._buffer.append(entry)

    def get_since(self, since_id: int = 0, limit: int = 200) -> List[Dict[str, Any]]:
        with self._lock:
            entries = [e for e in self._buffer if e["id"] > since_id]
        return entries[-limit:]


console_log_buffer = InMemoryLogHandler()

# Vyhrazený kanál pro krátké, obecné stavové zprávy určené pro Console panel
# na frontendu (builder začal/skončil/spadl) - NIKDY sem nedávat workspace_id,
# cesty k souborům ani jiné identifikátory. Běžné per-request logování
# (co dělá který endpoint, s jakými parametry) zůstává na modulových
# loggerech (`logging.getLogger(__name__)`) a jde jen do stdoutu/terminálu.
console_logger = logging.getLogger("app.console")
console_logger.addHandler(console_log_buffer)
console_logger.setLevel(logging.INFO)


def setup_logging():
    logging.basicConfig(
        level=logging.INFO,
        format="%(asctime)s [%(levelname)s] %(name)s: %(message)s",
        handlers=[logging.StreamHandler(sys.stdout)]
    )


logger = logging.getLogger(__name__)
