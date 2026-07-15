"""Digitize Fig. 2, 7, 9, 10 and 17 as independent external validation data.

The routine preserves raw matched pixels, data coordinates, overlays and an
explicit uncertainty budget.  Fig. 7 is decoded by matching plate pixels to the
printed color bar of each panel; its uncertainty is necessarily larger than a
line plot because of JPEG compression, antialiasing and black channel ribs.
"""

from __future__ import annotations

import argparse
from dataclasses import dataclass
import json
from pathlib import Path
import sys

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image
from scipy.ndimage import median_filter
from scipy.spatial import cKDTree

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from digitization.digitization_utils import CurveSpec, clean_trace, extract_curve, pixel_to_data, render_overlay, resolution_uncertainty  # noqa: E402


@dataclass(frozen=True)
class LineFigure:
    name: str
    rectangle: list[int]
    x_range: list[float]
    y_range: list[float]
    retained_y: list[int]
    curves: list[CurveSpec]
    exclusions: list[list[int]]
    x_column: str = "time_s"


LINE_FIGURES = {
    "Fig2": LineFigure(
        "Fig2", [139, 59, 1595, 1094], [0.0, 100.0], [0.0, 4.0e-4], [59, 1094],
        [CurveSpec("corotation_mass_flow_kg_s", (0, 114, 189), 62), CurveSpec("reversion_mass_flow_kg_s", (0, 180, 180), 62)],
        [[250, 80, 690, 270]], "duty_ratio_percent",
    ),
    "Fig9": LineFigure(
        "Fig9", [137, 23, 1608, 1080], [0.0, 60.0], [20.0, 45.0], [23, 1080],
        [CurveSpec("surface_1_air_C", (126, 47, 142), 58), CurveSpec("surface_2_air_C", (217, 83, 25), 58)],
        [[650, 30, 1550, 200]], "time_s",
    ),
    "Fig10": LineFigure(
        "Fig10", [132, 20, 1608, 1084], [0.0, 60.0], [45.0, 60.0], [20, 1084],
        [CurveSpec("surface_1_C", (0, 114, 189), 60), CurveSpec("surface_2_C", (217, 83, 25), 60)],
        [[900, 35, 1550, 200]], "time_s",
    ),
}


def masked_image(image: np.ndarray, exclusions: list[list[int]]) -> np.ndarray:
    result = image.copy()
    for x0, y0, x1, y1 in exclusions:
        result[y0:y1 + 1, x0:x1 + 1] = 255
    return result


def digitize_line_figure(specification: LineFigure, output: Path) -> dict[str, object]:
    image_path = PROJECT_ROOT / "source_images" / f"{specification.name}.jpg"
    image = np.asarray(Image.open(image_path).convert("RGB"))
    working = masked_image(image, specification.exclusions)
    raw_frames: list[pd.DataFrame] = []
    coordinate_frames: list[pd.DataFrame] = []
    cleaned: dict[str, pd.DataFrame] = {}
    metrics: dict[str, object] = {}
    x0, _, x1, _ = specification.rectangle
    sample_x = np.linspace(specification.x_range[0], specification.x_range[1], x1 - x0 + 1)
    final = pd.DataFrame({specification.x_column: sample_x})
    for curve in specification.curves:
        points, trace = extract_curve(working, curve, specification.rectangle, specification.retained_y)
        trace = clean_trace(trace, max_isolated_deviation_px=13.0)
        cleaned[curve.column] = trace
        data_x, data_y = pixel_to_data(
            trace.pixel_x.to_numpy(), trace.pixel_y_clean.to_numpy(), specification.rectangle,
            specification.x_range, specification.y_range,
        )
        final[curve.column] = np.interp(sample_x, data_x, data_y)
        px_x, px_y = pixel_to_data(
            points.pixel_x.to_numpy(), points.pixel_y.to_numpy(), specification.rectangle,
            specification.x_range, specification.y_range,
        )
        points[specification.x_column] = px_x
        points["data_value"] = px_y
        raw_frames.append(points)
        coordinate_frames.append(pd.DataFrame({"curve": curve.column, specification.x_column: data_x, "data_value": data_y}))
        metrics[curve.column] = {
            "matched_pixels": int(len(points)),
            "direct_valid_column_ratio": float(trace.observed_before_interpolation.mean()),
            "target_rgb": list(curve.rgb),
            "threshold": curve.threshold,
        }
    generated = output / "digitization" / "external"
    generated.mkdir(parents=True, exist_ok=True)
    final.to_csv(generated / f"{specification.name}_digitized.csv", index=False)
    pd.concat(raw_frames, ignore_index=True).to_csv(generated / f"{specification.name}_raw_pixels.csv", index=False)
    pd.concat(coordinate_frames, ignore_index=True).to_csv(generated / f"{specification.name}_coordinates.csv", index=False)
    render_overlay(image_path, generated / f"{specification.name}_overlay.png", cleaned, specification.curves, specification.rectangle)
    calibration = {
        "plot_rectangle_pixels": specification.rectangle,
        "x_range": specification.x_range,
        "y_range": specification.y_range,
    }
    return {
        "figure": specification.name,
        "axis_calibration": calibration,
        "resolution_uncertainty": resolution_uncertainty(calibration),
        "exclusion_rectangles": specification.exclusions,
        "curves": metrics,
    }


FIG7_PANELS = {
    # velocity: (plate crop, colour-bar crop, minimum K, maximum K)
    4: ([126, 2, 1400, 259], [31, 286, 1394, 311], 295.0, 338.0),
    5: ([127, 440, 1400, 696], [31, 727, 1394, 753], 295.0, 338.0),
    6: ([127, 881, 1400, 1138], [31, 1169, 1394, 1194], 295.0, 338.0),
    8: ([1705, 6, 2976, 262], [1608, 293, 2973, 318], 295.0, 335.0),
    10: ([1705, 447, 2976, 703], [1608, 734, 2973, 759], 295.0, 335.0),
    12: ([1705, 888, 2976, 1145], [1608, 1175, 2973, 1201], 295.0, 335.0),
}


def colorbar_palette(image: np.ndarray, crop: list[int], t_min: float, t_max: float) -> tuple[np.ndarray, np.ndarray]:
    x0, y0, x1, y1 = crop
    bar = image[y0:y1 + 1, x0:x1 + 1]
    # Median across the bar thickness suppresses border and tick artefacts.
    palette = np.median(bar, axis=0)
    temperature = np.linspace(t_min, t_max, palette.shape[0])
    valid = (palette.max(axis=1) - palette.min(axis=1) > 35) & (palette.mean(axis=1) < 245)
    return palette[valid], temperature[valid]


def decode_fig7_panel(image: np.ndarray, plate: list[int], bar: list[int], t_min: float, t_max: float, velocity: float) -> tuple[pd.DataFrame, dict[str, object], np.ndarray]:
    x0, y0, x1, y1 = plate
    region = image[y0:y1 + 1, x0:x1 + 1].astype(float)
    palette, palette_temperature = colorbar_palette(image, bar, t_min, t_max)
    # Downsample nearly identical adjacent colour-bar pixels and query with a KD-tree.
    select = np.linspace(0, len(palette) - 1, min(384, len(palette))).astype(int)
    palette = palette[select]
    palette_temperature = palette_temperature[select]
    flat = region.reshape(-1, 3)
    tree = cKDTree(palette)
    nearest_distance, index = tree.query(flat, k=1, workers=-1)
    decoded = palette_temperature[index]
    decoded = decoded.reshape(region.shape[:2])
    nearest_distance = nearest_distance.reshape(region.shape[:2])
    saturation = region.max(axis=2) - region.min(axis=2)
    # Exclude black ribs, white background, antialiased borders and poor palette matches.
    valid = (region.mean(axis=2) > 35) & (region.mean(axis=2) < 245) & (saturation > 25) & (nearest_distance < 55)
    decoded[~valid] = np.nan
    # Flow direction is bottom -> top in Fig. 7.  Use robust row medians across the central plate.
    row_temperature = np.nanmedian(decoded, axis=1)
    row_count = np.sum(np.isfinite(decoded), axis=1)
    row_temperature = pd.Series(row_temperature).interpolate(limit_direction="both").to_numpy()
    row_temperature = median_filter(row_temperature, size=5, mode="nearest")
    flow_coordinate = np.linspace(1.0, 0.0, len(row_temperature))  # raw top->bottom
    order = np.argsort(flow_coordinate)
    profile = pd.DataFrame({
        "velocity_m_s": velocity,
        "flow_coordinate_inlet_0_outlet_1": flow_coordinate[order],
        "temperature_K": row_temperature[order],
        "direct_pixel_count": row_count[order],
        "palette_uncertainty_K": np.full(len(row_temperature), 1.5),
    })
    # The first/last few rows contain sloped plate borders.  Exclude 3% at each end
    # from regional and extremum metrics while preserving raw decoded rows in the CSV.
    interior = profile[(profile.flow_coordinate_inlet_0_outlet_1 >= 0.03) & (profile.flow_coordinate_inlet_0_outlet_1 <= 0.97)].copy()
    third = np.array_split(interior.temperature_K.to_numpy(), 3)
    coordinate = interior.flow_coordinate_inlet_0_outlet_1.to_numpy()
    robust_temperature = interior.temperature_K.to_numpy()
    hotspot_index = int(np.argmax(robust_temperature))
    robust_max = float(np.quantile(robust_temperature, 0.98))
    robust_min = float(np.quantile(robust_temperature, 0.02))
    summary = {
        "velocity_m_s": velocity,
        "T_inlet_region_K": float(np.mean(third[0])),
        "T_middle_region_K": float(np.mean(third[1])),
        "T_outlet_region_K": float(np.mean(third[2])),
        "Tmax_K": robust_max,
        "Tmin_K": robust_min,
        "DeltaT_K": robust_max - robust_min,
        "hotspot_location_normalized": float(coordinate[hotspot_index]),
        "valid_pixel_fraction": float(np.mean(valid)),
        "reported_palette_uncertainty_K": 1.5,
    }
    return profile, summary, decoded


def digitize_fig7(output: Path) -> dict[str, object]:
    image_path = PROJECT_ROOT / "source_images" / "Fig7.jpg"
    image = np.asarray(Image.open(image_path).convert("RGB"))
    generated = output / "digitization" / "external"
    generated.mkdir(parents=True, exist_ok=True)
    profiles: list[pd.DataFrame] = []
    summaries: list[dict[str, object]] = []
    figure, axes = plt.subplots(2, 3, figsize=(13, 7), sharex=True, sharey=True)
    for axis, (velocity, (plate, bar, t_min, t_max)) in zip(axes.flat, FIG7_PANELS.items(), strict=True):
        profile, summary, _ = decode_fig7_panel(image, plate, bar, t_min, t_max, float(velocity))
        profiles.append(profile)
        summaries.append(summary)
        axis.plot(profile.flow_coordinate_inlet_0_outlet_1, profile.temperature_K)
        axis.fill_between(profile.flow_coordinate_inlet_0_outlet_1, profile.temperature_K - 1.5, profile.temperature_K + 1.5, alpha=0.2)
        axis.set_title(f"{velocity} m/s")
        axis.grid(alpha=0.25)
    for axis in axes[-1]:
        axis.set_xlabel("Normalized flow coordinate")
    for axis in axes[:, 0]:
        axis.set_ylabel("Temperature (K)")
    figure.suptitle("Fig. 7 palette-decoded flow-direction profiles")
    figure.tight_layout()
    figure.savefig(generated / "Fig7_profiles.png", dpi=200)
    plt.close(figure)
    pd.concat(profiles, ignore_index=True).to_csv(generated / "Fig7_spatial_profiles.csv", index=False)
    pd.DataFrame(summaries).to_csv(generated / "Fig7_region_summary.csv", index=False)
    return {
        "figure": "Fig7",
        "method": "nearest printed-colorbar palette decoding; robust row median; flow is bottom-to-top",
        "panels": summaries,
        "limitations": [
            "JPEG compression, antialiasing and black channel ribs introduce colour uncertainty.",
            "The temperature profile is a row-wise plate statistic, not original CFD node data.",
            "A conservative ±1.5 K palette uncertainty is attached to every profile point.",
        ],
    }


def digitize_fig17(output: Path) -> dict[str, object]:
    image_path = PROJECT_ROOT / "source_images" / "Fig17.jpg"
    image = np.asarray(Image.open(image_path).convert("RGB"))
    rectangle = [170, 17, 1723, 1188]
    x_range = [0.0, 1200.0]
    y_range = [295.0, 315.0]
    x0, y0, x1, y1 = rectangle
    crop = image[y0:y1 + 1, x0:x1 + 1].astype(float)
    targets = {"blue_signal": np.array((0, 114, 189), float), "orange_signal": np.array((217, 83, 25), float)}
    traces: dict[str, np.ndarray] = {}
    raw_frames: list[pd.DataFrame] = []
    pixel_x = np.arange(x0, x1 + 1)
    for name, target in targets.items():
        distance = np.linalg.norm(crop - target, axis=2)
        saturation = crop.max(axis=2) - crop.min(axis=2)
        mask = (distance < 65) & (saturation > 35)
        rows, columns = np.nonzero(mask)
        raw = pd.DataFrame({"signal": name, "pixel_x": columns + x0, "pixel_y": rows + y0, "color_distance": distance[rows, columns]})
        data_x, data_y = pixel_to_data(raw.pixel_x.to_numpy(), raw.pixel_y.to_numpy(), rectangle, x_range, y_range)
        raw["time_s"] = data_x
        raw["signal_value"] = data_y
        raw_frames.append(raw)
        trace = np.full(mask.shape[1], np.nan)
        for column in np.unique(columns):
            candidate = rows[columns == column]
            trace[column] = float(np.median(candidate))
        valid = np.isfinite(trace)
        trace = np.interp(np.arange(len(trace)), np.flatnonzero(valid), trace[valid])
        traces[name] = trace
    time_s, _ = pixel_to_data(pixel_x, np.full(len(pixel_x), y0), rectangle, x_range, y_range)
    _, blue_value = pixel_to_data(pixel_x, traces["blue_signal"] + y0, rectangle, x_range, y_range)
    _, orange_value = pixel_to_data(pixel_x, traces["orange_signal"] + y0, rectangle, x_range, y_range)
    difference = blue_value - orange_value
    sign = np.sign(median_filter(difference, size=5, mode="nearest"))
    changes = np.flatnonzero(sign[1:] * sign[:-1] < 0) + 1
    events: list[float] = []
    for index in changes:
        value = float(time_s[index])
        if value < 70.0:
            continue
        if not events or value - events[-1] > 5.0:
            events.append(value)
    generated = output / "digitization" / "external"
    generated.mkdir(parents=True, exist_ok=True)
    pd.concat(raw_frames, ignore_index=True).to_csv(generated / "Fig17_raw_pixels.csv", index=False)
    pd.DataFrame({"time_s": time_s, "blue_signal": blue_value, "orange_signal": orange_value}).to_csv(generated / "Fig17_digitized.csv", index=False)
    event_frame = pd.DataFrame({"event_index": np.arange(1, len(events) + 1), "switch_time_s": events})
    if len(events) > 1:
        event_frame["interval_from_previous_s"] = np.r_[np.nan, np.diff(events)]
    event_frame.to_csv(generated / "Fig17_switch_events.csv", index=False)
    fig, ax = plt.subplots(figsize=(10, 4))
    ax.plot(time_s, blue_value, label="blue signal")
    ax.plot(time_s, orange_value, label="orange signal")
    for event in events:
        ax.axvline(event, alpha=0.15)
    ax.set(xlabel="Time (s)", ylabel="State-machine plotted value", title="Fig. 17 digitized complementary switching signals")
    ax.legend(); ax.grid(alpha=0.25); fig.tight_layout(); fig.savefig(generated / "Fig17_overlay_reconstruction.png", dpi=180); plt.close(fig)
    return {
        "figure": "Fig17",
        "switch_event_count": len(events),
        "median_switch_interval_s": float(np.median(np.diff(events))) if len(events) > 1 else None,
        "uncertainty_s": 0.5 * (x_range[1] - x_range[0]) / (x1 - x0),
        "note": "Switches are sign changes of the two complementary colour traces after 70 s.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = (args.output_dir or root / "outputs_v5").resolve()
    audits: dict[str, object] = {}
    for specification in LINE_FIGURES.values():
        audits[specification.name] = digitize_line_figure(specification, output)
    audits["Fig7"] = digitize_fig7(output)
    audits["Fig17"] = digitize_fig17(output)
    target = output / "digitization" / "external_digitization_audit.json"
    target.parent.mkdir(parents=True, exist_ok=True)
    target.write_text(json.dumps(audits, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps({name: "completed" for name in audits}, indent=2))


if __name__ == "__main__":
    main()
