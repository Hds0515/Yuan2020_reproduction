"""Audit regenerated image-only curves against V2 CSVs after extraction."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

PROJECT_ROOT = Path(__file__).resolve().parents[1]


def compare(regenerated: Path, baseline: Path) -> dict[str, object]:
    new = pd.read_csv(regenerated)
    old = pd.read_csv(baseline)
    common = [c for c in new.columns if c != "time_s" and c in old.columns]
    metrics: dict[str, object] = {}
    for column in common:
        reference = np.interp(new["time_s"], old["time_s"], old[column])
        residual = new[column].to_numpy() - reference
        metrics[column] = {
            "rmse": float(np.sqrt(np.mean(residual**2))),
            "mae": float(np.mean(np.abs(residual))),
            "max_abs": float(np.max(np.abs(residual))),
        }
    return metrics


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, default=PROJECT_ROOT)
    args = parser.parse_args()
    result = {}
    for name in ("Fig14", "Fig15", "Fig16"):
        result[name] = compare(
            args.output_root / "digitization" / "regenerated" / f"{name}_digitized.csv",
            PROJECT_ROOT / "digitization" / f"{name}_cleaned_1s.csv",
        )
    path = args.output_root / "digitization" / "v2_comparison.json"
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
