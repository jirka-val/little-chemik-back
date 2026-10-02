"""Ties Console panel messages to the workspace of the request that produced them."""

import re

from app.core.logging import console_workspace

_WORKSPACE_IN_PATH = re.compile(r"/([0-9a-f]{8}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{4}-[0-9a-f]{12})(?:/|$)")


class ConsoleWorkspaceMiddleware:
    """
    Pure ASGI middleware: sets `console_workspace` from a workspace UUID in
    the request path, so console_logger messages emitted while serving it
    (also in run_in_threadpool) are shown only to that workspace. Requests
    without one in the path leave the messages untagged (visible to everyone).
    """

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        match = _WORKSPACE_IN_PATH.search(scope.get("path", ""))
        token = console_workspace.set(match.group(1) if match else None)
        try:
            await self.app(scope, receive, send)
        finally:
            console_workspace.reset(token)
