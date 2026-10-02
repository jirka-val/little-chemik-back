"""Shared access checks for the HTTP layer."""

import secrets

from app.core.config import settings
from app.core.exceptions import ForbiddenError


def require_admin(token: str) -> None:
    """
    Raise ForbiddenError unless `token` equals Settings.ADMIN_TOKEN. An empty
    ADMIN_TOKEN keeps every admin endpoint locked. The comparison takes the
    same time whatever the token, so it cannot be guessed character by character.
    """
    if not settings.ADMIN_TOKEN or not secrets.compare_digest(token.encode(), settings.ADMIN_TOKEN.encode()):
        raise ForbiddenError()
