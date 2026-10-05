"""
Console panel ("System output"): každý workspace má vlastní buffer, takže
zprávy jednoho uživatele nevytlačí zprávy jiného (app/core/logging.py).
"""

import logging

import pytest

from app.core.logging import InMemoryLogHandler, console_workspace

pytestmark = pytest.mark.unit


def _emit(handler: InMemoryLogHandler, message: str, workspace_id=None) -> None:
    token = console_workspace.set(workspace_id)
    try:
        handler.emit(logging.LogRecord("app.console", logging.INFO, __file__, 0, message, None, None))
    finally:
        console_workspace.reset(token)


def test_workspaces_do_not_evict_each_other():
    handler = InMemoryLogHandler(general_capacity=5, workspace_capacity=3)
    _emit(handler, "mine", "ws-a")
    for i in range(10):
        _emit(handler, f"other {i}", "ws-b")
    assert [e["message"] for e in handler.get_since(workspace_id="ws-a")] == ["mine"]
    assert [e["message"] for e in handler.get_since(workspace_id="ws-b")] == ["other 7", "other 8", "other 9"]


def test_general_and_own_messages_in_order():
    handler = InMemoryLogHandler()
    _emit(handler, "general 1")
    _emit(handler, "own", "ws-a")
    _emit(handler, "foreign", "ws-b")
    _emit(handler, "general 2")
    entries = handler.get_since(workspace_id="ws-a")
    assert [e["message"] for e in entries] == ["general 1", "own", "general 2"]
    assert "_workspace" not in entries[0]
    since = handler.get_since(since_id=entries[1]["id"], workspace_id="ws-a")
    assert [e["message"] for e in since] == ["general 2"]
    assert [e["message"] for e in handler.get_since()] == ["general 1", "general 2"]


def test_oldest_workspace_dropped_beyond_limit():
    handler = InMemoryLogHandler(max_workspaces=2)
    _emit(handler, "a", "ws-a")
    _emit(handler, "b", "ws-b")
    _emit(handler, "a again", "ws-a")
    _emit(handler, "c", "ws-c")
    assert handler.get_since(workspace_id="ws-b") == []
    assert [e["message"] for e in handler.get_since(workspace_id="ws-a")] == ["a", "a again"]
