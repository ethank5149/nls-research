"""Centralized configuration for the figures package."""

import os
from pathlib import Path

# Resolve output directory relative to this project's figures/ directory,
# regardless of the current working directory. This makes saves work the
# same whether running locally on unraid or over SMB from another machine.
_PROJECT_ROOT = Path(__file__).resolve().parent  # figures/
OUTPUT_DIR = Path(os.environ.get("FIGURES_OUTPUT_DIR", _PROJECT_ROOT / "output"))
OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
