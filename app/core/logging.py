import logging
import sys
import threading
from collections import OrderedDict, deque
from contextvars import ContextVar
from typing import Any, Dict, List, Optional


class InMemoryLogHandler(logging.Handler):
    """
    Drží poslední log záznamy v paměti procesu, ať je frontend
    (Console panel v sidebaru) může pollovat přes /api/system/logs.

    Připojený POUZE na `console_logger` (viz níže), ne na root logger - ten
    zůstává jen ve stdoutu/terminálu. Bez tohohle oddělení by se do panelu
    (viditelného v prohlížeči) dostalo úplně všechno včetně workspace ID a
    dalších detailů z běžných per-request logů, což už jednou byl reálný
    problém (unik citlivých identifikátorů + zahlcení nedůležitými zprávami).
    Sem smí jen krátké, obecné stavové zprávy bez identifikátorů.

    Každý záznam dostane rostoucí `id`, takže frontend může pollovat
    přírůstkově (`since_id`) místo opakovaného posílání celého bufferu.

    Záznam si pamatuje workspace, v jehož requestu vznikl (console_workspace,
    nastavuje ConsoleWorkspaceMiddleware), aby jeden uživatel neviděl zprávy
    z přípravy jiného. Samotné ID se ven nikdy neposílá.
    """

    def __init__(self, general_capacity: int = 200, workspace_capacity: int = 1000, max_workspaces: int = 500):
        super().__init__()
        # Každý workspace má vlastní buffer ("vlastní terminál") - dřív byl
        # jeden společný na 100 záznamů a zprávy jednoho uživatele vytlačovaly
        # zprávy ostatních. Obecné zprávy (bez workspace) mají svůj.
        self._general: deque[Dict[str, Any]] = deque(maxlen=general_capacity)
        self._by_workspace: "OrderedDict[str, deque[Dict[str, Any]]]" = OrderedDict()
        self._workspace_capacity = workspace_capacity
        self._max_workspaces = max_workspaces
        self._lock = threading.Lock()
        self._next_id = 1

    def emit(self, record: logging.LogRecord) -> None:
        try:
            message = record.getMessage()
        except Exception:
            message = str(record.msg)

        workspace_id = console_workspace.get()
        with self._lock:
            entry = {
                "id": self._next_id,
                "timestamp": record.created,
                "level": record.levelname,
                "message": message,
            }
            self._next_id += 1
            if workspace_id is None:
                self._general.append(entry)
                return
            buffer = self._by_workspace.get(workspace_id)
            if buffer is None:
                buffer = self._by_workspace[workspace_id] = deque(maxlen=self._workspace_capacity)
                # Nejdéle nepoužité workspace zahodíme, ať paměť neroste.
                while len(self._by_workspace) > self._max_workspaces:
                    self._by_workspace.popitem(last=False)
            else:
                self._by_workspace.move_to_end(workspace_id)
            buffer.append(entry)

    def get_since(self, since_id: int = 0, limit: int = 1000, workspace_id: Optional[str] = None) -> List[Dict[str, Any]]:
        """Záznamy novější než since_id: obecné (bez workspace) + ty z daného workspace."""
        with self._lock:
            entries = [dict(e) for e in self._general if e["id"] > since_id]
            if workspace_id is not None:
                entries += [dict(e) for e in self._by_workspace.get(workspace_id, ()) if e["id"] > since_id]
        entries.sort(key=lambda e: e["id"])
        return entries[-limit:]


# Workspace aktuálního requestu (viz ConsoleWorkspaceMiddleware). Přenáší se
# i do run_in_threadpool, protože ten kopíruje kontext.
console_workspace: ContextVar[Optional[str]] = ContextVar("console_workspace", default=None)

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
