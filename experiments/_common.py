"""Shared helpers: provenance stamping and result writing for simulation experiments."""
from __future__ import annotations

import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "src"))


def stamp(config: dict) -> dict:
    return {
        "simulation_only": True,
        "warning": "Simulated data with known ground truth. Not biological evidence.",
        "utc_time": time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime()),
        "python": platform.python_version(),
        "numpy": np.__version__,
        "config": config,
    }


def write_results(name: str, df: pd.DataFrame, config: dict, out_dir: Path | None = None) -> Path:
    out = (out_dir or ROOT / "results") / name
    out.mkdir(parents=True, exist_ok=True)
    df.to_csv(out / "summary.csv", index=False)
    (out / "provenance.json").write_text(json.dumps(stamp(config), indent=2, default=str))
    (out / "summary.md").write_text(df.to_markdown(index=False, floatfmt=".3f")
                                    if _has_tabulate() else df.to_string(index=False))
    return out


def _has_tabulate() -> bool:
    try:
        import tabulate  # noqa: F401
        return True
    except ImportError:
        return False
