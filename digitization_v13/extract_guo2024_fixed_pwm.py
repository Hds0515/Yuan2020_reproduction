"""Digitize Guo et al. (2024) Fig. 5(e,f) from the publisher image."""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
OUT = ROOT / "outputs_v13" / "digitization"
EXPECTED_SIZE = (811, 984)
EXPECTED_SHA256 = "A50F12D8C0653DB83A0CAFC1D2F7961689CFB1BF0E65FA42831263712B4EEA48"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().upper()


def extract(path: Path) -> pd.DataFrame:
    if sha256(path) != EXPECTED_SHA256:
        raise RuntimeError("publisher image hash does not match the preregistration")
    image = Image.open(path).convert("RGB")
    if image.size != EXPECTED_SIZE:
        raise RuntimeError(f"expected {EXPECTED_SIZE}, got {image.size}")
    pixels = np.asarray(image)
    masks = {
        "outlet": lambda q: (q[:, 0] > 150) & (q[:, 0] > q[:, 1] * 1.25) & (q[:, 0] > q[:, 2] * 1.25),
        "intermediate": lambda q: (q[:, 2] > 100) & (q[:, 2] > q[:, 0] * 1.15) & (q[:, 2] > q[:, 1] * 1.05),
        "inlet": lambda q: (q[:, 1] > 70) & (q[:, 1] > q[:, 0] * 1.15) & (q[:, 1] > q[:, 2] * 1.05),
    }
    peaks = {"outlet": 64.0, "intermediate": 60.0, "inlet": 52.0}
    panels = {
        "ramp_2A": ("calibration", 69, 391),
        "ramp_4A": ("final_confirmation", 472, 795),
    }
    y_top, y_bottom = 679, 950
    rows: list[dict[str, object]] = []
    for scenario, (purpose, x_left, x_right) in panels.items():
        for time_s in range(50, 2000, 100):
            x_pixel = round(x_left + time_s / 2000.0 * (x_right - x_left))
            triangular_fraction = 1.0 - abs(time_s - 1000.0) / 1000.0
            for region, predicate in masks.items():
                candidate_y: list[int] = []
                for x in range(x_pixel - 2, x_pixel + 3):
                    segment = pixels[676:952, x, :]
                    candidate_y.extend((np.where(predicate(segment))[0] + 676).tolist())
                temperatures = 20.0 + (y_bottom - np.asarray(candidate_y)) * 48.0 / (y_bottom - y_top)
                expected = 23.0 + (peaks[region] - 23.0) * triangular_fraction
                close = temperatures[np.abs(temperatures - expected) < 5.0]
                if len(close) == 0:
                    raise RuntimeError(f"no colour-path match for {scenario} {time_s} {region}")
                temperature = float(np.median(close))
                y_pixel = int(round(y_bottom - (temperature - 20.0) / 48.0 * (y_bottom - y_top)))
                rows.append(
                    {
                        "scenario": scenario,
                        "purpose": purpose,
                        "time_s": time_s,
                        "flow_region": region,
                        "temperature_C": temperature,
                        "pixel_x": x_pixel,
                        "pixel_y": y_pixel,
                        "uncertainty_C": 0.75,
                        "source": "Guo2024 Fig.5(e)" if scenario == "ramp_2A" else "Guo2024 Fig.5(f)",
                    }
                )
    return pd.DataFrame(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("figure5", type=Path)
    args = parser.parse_args()
    OUT.mkdir(parents=True, exist_ok=True)
    frame = extract(args.figure5)
    csv_path = OUT / "guo2024_fixed_pwm_regional_temperatures.csv"
    frame.to_csv(csv_path, index=False)
    audit = {
        "source_url": "https://ars.els-cdn.com/content/image/1-s2.0-S0360319924000478-gr5.jpg",
        "source_sha256": sha256(args.figure5),
        "source_size_pixels": list(EXPECTED_SIZE),
        "points": len(frame),
        "calibration_points": int((frame.purpose == "calibration").sum()),
        "final_confirmation_points": int((frame.purpose == "final_confirmation").sum()),
        "output_sha256": sha256(csv_path),
    }
    (OUT / "digitization_audit.json").write_text(json.dumps(audit, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
