"""Extract the Gong et al. (2024) 75-point temperature field.

The script operates on page 11 rendered at 180 dpi with pdftoppm.  Fixed pixel
coordinates are part of the frozen audit.  Temperatures are obtained by nearest
colour matching against the figure's own 40--65 degC colour bar.  Dynamic
anchors are manual readings from Fig. 9 and carry a wider uncertainty field.
"""

from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image


ROOT = Path(__file__).resolve().parents[1]
OUTPUT = ROOT / "outputs_v11" / "digitization"
EXPECTED_PAGE_SIZE = (1489, 1985)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def extract_static(page_path: Path) -> pd.DataFrame:
    image = Image.open(page_path).convert("RGB")
    if image.size != EXPECTED_PAGE_SIZE:
        raise RuntimeError(
            f"Expected a 180-dpi page of size {EXPECTED_PAGE_SIZE}, got {image.size}"
        )
    pixels = np.asarray(image)
    colourbar_y = np.arange(180, 517)
    colourbar_rgb = pixels[colourbar_y, 260].astype(float)
    colourbar_temperature = 65.0 - (colourbar_y - 180.0) / (516.0 - 180.0) * 25.0
    rectangles = ((206, 353), (358, 506), (508, 659), (663, 812), (818, 968))
    rows: list[dict[str, object]] = []
    # Source rectangles run from hot outlet at the top to cool inlet at bottom.
    row_fractions = (5.0 / 6.0, 3.0 / 6.0, 1.0 / 6.0)
    column_fractions = (0.1, 0.3, 0.5, 0.7, 0.9)
    for group, (x0, x1) in enumerate(rectangles, start=1):
        for flow_row, row_fraction in enumerate(row_fractions, start=1):
            y_pixel = round(861 + row_fraction * (943 - 861))
            for column, column_fraction in enumerate(column_fractions, start=1):
                x_pixel = round(x0 + column_fraction * (x1 - x0))
                rgb = pixels[y_pixel, x_pixel].astype(float)
                distance = np.sum(np.square(colourbar_rgb - rgb), axis=1)
                match = int(np.argmin(distance))
                rows.append(
                    {
                        "cell_group": group,
                        "flow_row": flow_row,
                        "in_plane_column": column,
                        "temperature_C": float(colourbar_temperature[match]),
                        "x_pixel": x_pixel,
                        "y_pixel": y_pixel,
                        "rgb_r": int(rgb[0]),
                        "rgb_g": int(rgb[1]),
                        "rgb_b": int(rgb[2]),
                        "colour_distance": float(np.sqrt(distance[match])),
                        "digitization_uncertainty_C": 0.75,
                        "source": "Gong2024 Fig.8(b)",
                    }
                )
    return pd.DataFrame(rows)


def dynamic_anchors() -> pd.DataFrame:
    # Manual curve readings locked before model calibration. The larger stated
    # uncertainty covers line thickness, axes interpolation, and rasterisation.
    uav = (
        (0, 16.0, 2.0),
        (100, 16.0, 2.0),
        (200, 20.0, 2.0),
        (280, 23.0, 30.0),
        (300, 25.0, 26.0),
        (330, 35.0, 26.0),
        (360, 44.0, 26.0),
        (400, 50.0, 26.0),
        (450, 54.0, 26.0),
        (500, 55.5, 26.0),
        (580, 56.0, 26.0),
        (610, 53.0, 20.0),
        (640, 41.0, 20.0),
        (680, 32.0, 0.0),
    )
    random = (
        (0, 30.0, 0.0),
        (100, 27.5, 1.0),
        (200, 26.0, 1.0),
        (320, 25.5, 8.0),
        (360, 31.0, 10.0),
        (400, 39.0, 20.0),
        (450, 46.0, 15.0),
        (500, 47.0, 24.0),
        (550, 52.0, 30.0),
        (600, 57.0, 30.0),
        (650, 54.0, 16.0),
        (700, 52.0, 8.0),
        (750, 50.5, 10.0),
        (850, 51.0, 12.0),
        (900, 44.0, 1.0),
        (970, 39.5, 0.0),
    )
    records: list[dict[str, object]] = []
    for scenario, values, purpose in (
        ("uav_like", uav, "calibration"),
        ("random_load", random, "final_confirmation"),
    ):
        for time_s, temperature_C, current_A in values:
            records.append(
                {
                    "scenario": scenario,
                    "purpose": purpose,
                    "time_s": time_s,
                    "center_temperature_C": temperature_C,
                    "load_current_A": current_A,
                    "digitization_uncertainty_C": 1.25,
                    "source": f"Gong2024 Fig.9 {'(a)' if scenario == 'uav_like' else '(b)'}",
                }
            )
    return pd.DataFrame(records)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("page11", type=Path)
    args = parser.parse_args()
    OUTPUT.mkdir(parents=True, exist_ok=True)
    static = extract_static(args.page11)
    dynamic = dynamic_anchors()
    static_path = OUTPUT / "gong2024_fig8_75_temperatures.csv"
    dynamic_path = OUTPUT / "gong2024_fig9_dynamic_anchors.csv"
    static.to_csv(static_path, index=False)
    dynamic.to_csv(dynamic_path, index=False)
    audit = {
        "rendered_page_path": str(args.page11),
        "rendered_page_sha256": sha256(args.page11),
        "render_dpi": 180,
        "static_points": len(static),
        "dynamic_anchors": len(dynamic),
        "static_temperature_range_C": [
            float(static.temperature_C.min()),
            float(static.temperature_C.max()),
        ],
        "outputs": {
            static_path.name: sha256(static_path),
            dynamic_path.name: sha256(dynamic_path),
        },
    }
    (OUTPUT / "digitization_audit.json").write_text(
        json.dumps(audit, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(audit, indent=2))


if __name__ == "__main__":
    main()
