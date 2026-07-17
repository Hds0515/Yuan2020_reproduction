"""Digitize PWM40/60/80 average-temperature paths from Guo Fig. 11."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs_v14" / "digitization"
EXPECTED_SIZE = (691, 442)
EXPECTED_SHA = "B0DF5C5A0E33D42DF6421B801A98234F88CA56A4594FAD845CA45525856C87D8"


def digest(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("figure11", type=Path)
    args = parser.parse_args()
    if digest(args.figure11) != EXPECTED_SHA:
        raise RuntimeError("unexpected publisher image hash")
    image = Image.open(args.figure11).convert("RGB")
    if image.size != EXPECTED_SIZE:
        raise RuntimeError(f"expected {EXPECTED_SIZE}, got {image.size}")
    pixels = np.asarray(image)
    masks = {
        40: lambda q: (q[:, 0] > 140) & (q[:, 0] > q[:, 1] * 1.25) & (q[:, 0] > q[:, 2] * 1.15),
        60: lambda q: (q[:, 2] > 100) & (q[:, 2] > q[:, 0] * 1.15) & (q[:, 2] > q[:, 1] * 1.05),
        80: lambda q: (q[:, 1] > 70) & (q[:, 1] > q[:, 0] * 1.10) & (q[:, 1] > q[:, 2] * 1.05),
    }
    expected = {
        40: lambda t: 34.0 + (t - 150.0) * (25.0 / 600.0),
        60: lambda t: 33.7 + (t - 150.0) * (24.6 / 800.0),
        80: lambda t: 32.5 + (t - 150.0) * (19.5 / 800.0),
    }
    purposes = {40: "calibration", 60: "validation", 80: "final_confirmation"}
    x_left, x_right = 48, 630
    y_top, y_bottom = 214, 383
    rows = []
    for pwm in (40, 60, 80):
        end = 750 if pwm == 40 else 950
        for time_s in range(150, end + 1, 100):
            x = round(x_left + time_s / 1100.0 * (x_right - x_left))
            candidate_y: list[int] = []
            for xx in range(x - 2, x + 3):
                segment = pixels[205:390, xx, :]
                candidate_y.extend((np.where(masks[pwm](segment))[0] + 205).tolist())
            temperatures = 20.0 + (y_bottom - np.asarray(candidate_y)) * 40.0 / (y_bottom - y_top)
            temperatures = temperatures[(temperatures >= 20.0) & (temperatures <= 65.0)]
            if len(temperatures) == 0:
                raise RuntimeError(f"no path pixels for PWM{pwm} at {time_s}s")
            value = float(temperatures[np.argmin(np.abs(temperatures - expected[pwm](time_s)))])
            y = int(round(y_bottom - (value - 20.0) / 40.0 * (y_bottom - y_top)))
            rows.append(
                {
                    "PWM_percent": pwm,
                    "purpose": purposes[pwm],
                    "time_s": time_s,
                    "average_temperature_C": value,
                    "pixel_x": x,
                    "pixel_y": y,
                    "uncertainty_C": 1.0,
                    "source": "Guo2024 Fig.11",
                }
            )
    OUT.mkdir(parents=True, exist_ok=True)
    csv_path = OUT / "guo2024_multi_pwm_average_temperature.csv"
    pd.DataFrame(rows).to_csv(csv_path, index=False)
    audit = {
        "source_sha256": digest(args.figure11),
        "source_size_pixels": list(EXPECTED_SIZE),
        "points": len(rows),
        "output_sha256": digest(csv_path),
        "excluded_PWM100_reason": "paper reports flooding; thermal-only scope",
    }
    (OUT / "digitization_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
