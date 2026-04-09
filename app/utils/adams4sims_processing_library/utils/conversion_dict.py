from __future__ import annotations

import json
from functools import lru_cache
from importlib.resources import files
from typing import Any

_RESOURCE_NAME = "converting_dictionary.json"


@lru_cache(maxsize=1)
def get_conversion_dict() -> dict[str, Any]:
    """
    Return the bundled conversion dictionary as a cached mapping.
    """
    resource = files(__package__).joinpath(_RESOURCE_NAME)
    with resource.open("r", encoding="utf-8") as f:
        try:
            return json.load(f)
        except json.JSONDecodeError as e:
            raise RuntimeError(
                "Bundled converting_dictionary.json is malformed. Please report this issue."
            ) from e


__all__ = ["get_conversion_dict"]
