"""Compatibility entry point for the canonical portable bandwidth checker."""
from pathlib import Path
import runpy
import sys

script = Path(__file__).resolve().parents[1] / "portable_checks" / "rev" / "bandwidth_check.py"
sys.path.insert(0, str(script.parent))
runpy.run_path(str(script), run_name="__main__")
