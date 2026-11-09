"""
Kdo poslal právě obsluhovaný požadavek - pro limity na klienta (job_queue.py).

uvicorn běží s --proxy-headers, takže scope["client"] je skutečná adresa
klienta z X-Forwarded-For od nginx, ne adresa nginx kontejneru.
"""

from contextvars import ContextVar
from typing import Optional

client_address: ContextVar[Optional[str]] = ContextVar("client_address", default=None)


class ClientAddressMiddleware:
    """Pure ASGI middleware: nastaví client_address pro dobu požadavku."""

    def __init__(self, app):
        self.app = app

    async def __call__(self, scope, receive, send):
        if scope["type"] != "http":
            await self.app(scope, receive, send)
            return
        client = scope.get("client")
        token = client_address.set(client[0] if client else None)
        try:
            await self.app(scope, receive, send)
        finally:
            client_address.reset(token)
