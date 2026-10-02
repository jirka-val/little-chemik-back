"""
Structure preparation services built on the FORGE builder (app/builder).

The builder is vendored as a flat set of modules that import each other by
bare name (`from forge_molecule_parser import ...`), so its directory has to
be on sys.path. This package puts it there once; every module in it can then
import builder modules directly.
"""

import sys
from pathlib import Path

BUILDER_DIR = Path(__file__).resolve().parents[2] / "builder"
if str(BUILDER_DIR) not in sys.path:
    sys.path.insert(0, str(BUILDER_DIR))
