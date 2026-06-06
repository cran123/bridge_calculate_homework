#!/usr/bin/env python3

from pathlib import Path
import runpy


if __name__ == "__main__":
    script = Path(__file__).resolve().parents[1] / "others" / "tools" / "abaqus_to_stappp.py"
    runpy.run_path(str(script), run_name="__main__")
