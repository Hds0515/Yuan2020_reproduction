"""Regenerate Fig. 14-16 CSV files directly from the source JPG pixels."""

from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from PIL import Image

try:
    from .digitization_utils import (
        CurveSpec,
        clean_trace,
        extract_curve,
        pixel_to_data,
        render_overlay,
        resolution_uncertainty,
    )
except ImportError:  # Direct script execution.
    from digitization_utils import (
        CurveSpec,
        clean_trace,
        extract_curve,
        pixel_to_data,
        render_overlay,
        resolution_uncertainty,
    )

PROJECT_ROOT = Path(__file__).resolve().parents[1]

FIGURE_CURVES = {
    "Fig14": [
        CurveSpec("experiment_T1_C", (120, 48, 136), 46),
        CurveSpec("experiment_T2_C", (216, 88, 32), 46),
        CurveSpec("experiment_T3_C", (232, 180, 40), 46),
        CurveSpec("paper_simulation_T1_C", (160, 24, 56), 42),
        CurveSpec("paper_simulation_T2_C", (120, 168, 48), 42),
        CurveSpec("paper_simulation_T3_C", (88, 184, 224), 46),
        CurveSpec("reference_C", (8, 248, 8), 52),
    ],
    "Fig15": [
        CurveSpec("experiment_current_A", (8, 112, 184), 52),
        CurveSpec("paper_simulation_current_A", (248, 8, 248), 48),
    ],
    "Fig16": [
        CurveSpec("experiment_middle_C", (8, 112, 184), 50),
        CurveSpec("experiment_T1_C", (216, 88, 32), 46),
        CurveSpec("experiment_T2_C", (232, 180, 40), 48),
        CurveSpec("paper_simulation_middle_C", (120, 48, 136), 46),
        CurveSpec("paper_simulation_T1_C", (120, 168, 48), 44),
        CurveSpec("paper_simulation_T2_C", (88, 184, 224), 48),
        CurveSpec("reference_C", (248, 8, 248), 48),
    ],
}

FIG16_DETAIL_CURVES = [
    CurveSpec("experiment_middle_C", (8, 112, 184), 50),
    CurveSpec("experiment_T1_C", (216, 88, 32), 46),
    CurveSpec("experiment_T2_C", (232, 180, 40), 48),
    CurveSpec("paper_simulation_middle_C", (120, 48, 136), 46),
    CurveSpec("paper_simulation_T1_C", (120, 168, 48), 44),
    CurveSpec("paper_simulation_T2_C", (88, 184, 224), 48),
]


def digitize_figure(name: str, output_root: Path) -> dict[str, object]:
    calibration_path = PROJECT_ROOT / "digitization" / "axis_calibration.json"
    calibration = json.loads(calibration_path.read_text(encoding="utf-8"))[name]
    image_path = PROJECT_ROOT / "source_images" / f"{name}.jpg"
    image = np.asarray(Image.open(image_path).convert("RGB"))
    specs = FIGURE_CURVES[name]

    regenerated = output_root / "digitization" / "regenerated"
    intermediate = regenerated / "intermediate"
    overlay_dir = output_root / "digitization"
    regenerated.mkdir(parents=True, exist_ok=True)
    intermediate.mkdir(parents=True, exist_ok=True)

    all_points: list[pd.DataFrame] = []
    all_preclean: list[pd.DataFrame] = []
    all_coordinates: list[pd.DataFrame] = []
    cleaned: dict[str, pd.DataFrame] = {}
    final = pd.DataFrame({"time_s": np.arange(0.0, 1201.0, 1.0)})
    curve_metrics: dict[str, object] = {}

    for spec in specs:
        points, trace = extract_curve(
            image,
            spec,
            calibration["plot_rectangle_pixels"],
            calibration["retained_y_pixels"],
        )
        trace = clean_trace(trace)
        cleaned[spec.column] = trace
        point_x, point_y = pixel_to_data(
            points["pixel_x"].to_numpy(),
            points["pixel_y"].to_numpy(),
            calibration["plot_rectangle_pixels"],
            calibration["x_range"],
            calibration["y_range"],
        )
        points["data_x"] = point_x
        points["data_y"] = point_y
        all_points.append(points)

        raw_x, raw_y = pixel_to_data(
            trace["pixel_x"].to_numpy(),
            trace["pixel_y_raw"].to_numpy(),
            calibration["plot_rectangle_pixels"],
            calibration["x_range"],
            calibration["y_range"],
        )
        preclean = trace.copy()
        preclean["data_x"] = raw_x
        preclean["data_y_raw"] = raw_y
        all_preclean.append(preclean)

        clean_x, clean_y = pixel_to_data(
            trace["pixel_x"].to_numpy(),
            trace["pixel_y_clean"].to_numpy(),
            calibration["plot_rectangle_pixels"],
            calibration["x_range"],
            calibration["y_range"],
        )
        coordinate = pd.DataFrame({"curve": spec.column, "data_x": clean_x, "data_y": clean_y})
        all_coordinates.append(coordinate)
        final[spec.column] = np.interp(final["time_s"], clean_x, clean_y)
        observed = trace["observed_before_interpolation"].to_numpy(bool)
        curve_metrics[spec.column] = {
            "valid_column_ratio": float(np.mean(observed)),
            "matched_pixel_count": int(len(points)),
            "rejected_column_count": int(trace["rejected_as_isolated"].sum()),
            "target_rgb": list(spec.rgb),
            "color_distance_threshold": spec.threshold,
        }

    pd.concat(all_points, ignore_index=True).to_csv(
        intermediate / f"{name}_pixel_points.csv", index=False
    )
    pd.concat(all_preclean, ignore_index=True).to_csv(
        intermediate / f"{name}_preclean.csv", index=False
    )
    pd.concat(all_coordinates, ignore_index=True).to_csv(
        intermediate / f"{name}_data_coordinates.csv", index=False
    )
    final.to_csv(regenerated / f"{name}_digitized.csv", index=False)
    render_overlay(
        image_path,
        overlay_dir / f"{name}_pixel_overlay.png",
        cleaned,
        specs,
        calibration["plot_rectangle_pixels"],
    )
    return {
        "figure": name,
        "source_image": str(image_path),
        "source_resolution_pixels": list(Image.open(image_path).size),
        "plot_rectangle_pixels": calibration["plot_rectangle_pixels"],
        "uncertainty": resolution_uncertainty(calibration),
        "curves": curve_metrics,
    }


def digitize_fig16_detail(panel: str, output_root: Path) -> dict[str, object]:
    """Digitize Fig. 16(b)/(c) directly instead of cropping the main-panel trace."""
    calibration = json.loads(
        (PROJECT_ROOT / "digitization" / "axis_calibration.json").read_text(encoding="utf-8")
    )[panel]
    image_path = PROJECT_ROOT / "source_images" / "Fig16.jpg"
    image = np.asarray(Image.open(image_path).convert("RGB"))
    regenerated = output_root / "digitization" / "regenerated"
    intermediate = regenerated / "intermediate"
    overlay_dir = output_root / "digitization"
    regenerated.mkdir(parents=True, exist_ok=True)
    intermediate.mkdir(parents=True, exist_ok=True)
    start, stop = map(int, calibration["x_range"])
    final = pd.DataFrame({"time_s": np.arange(start, stop + 1, 1.0)})
    all_points: list[pd.DataFrame] = []
    all_preclean: list[pd.DataFrame] = []
    all_coordinates: list[pd.DataFrame] = []
    cleaned: dict[str, pd.DataFrame] = {}
    curve_metrics: dict[str, object] = {}
    for spec in FIG16_DETAIL_CURVES:
        points, trace = extract_curve(
            image,
            spec,
            calibration["plot_rectangle_pixels"],
            calibration["retained_y_pixels"],
        )
        trace = clean_trace(trace)
        cleaned[spec.column] = trace
        point_x, point_y = pixel_to_data(
            points["pixel_x"].to_numpy(),
            points["pixel_y"].to_numpy(),
            calibration["plot_rectangle_pixels"],
            calibration["x_range"],
            calibration["y_range"],
        )
        points["data_x"] = point_x
        points["data_y"] = point_y
        all_points.append(points)
        raw_x, raw_y = pixel_to_data(
            trace["pixel_x"].to_numpy(),
            trace["pixel_y_raw"].to_numpy(),
            calibration["plot_rectangle_pixels"],
            calibration["x_range"],
            calibration["y_range"],
        )
        preclean = trace.copy()
        preclean["data_x"] = raw_x
        preclean["data_y_raw"] = raw_y
        all_preclean.append(preclean)
        clean_x, clean_y = pixel_to_data(
            trace["pixel_x"].to_numpy(),
            trace["pixel_y_clean"].to_numpy(),
            calibration["plot_rectangle_pixels"],
            calibration["x_range"],
            calibration["y_range"],
        )
        all_coordinates.append(
            pd.DataFrame({"curve": spec.column, "data_x": clean_x, "data_y": clean_y})
        )
        final[spec.column] = np.interp(final["time_s"], clean_x, clean_y)
        observed = trace["observed_before_interpolation"].to_numpy(bool)
        curve_metrics[spec.column] = {
            "direct_valid_column_ratio": float(np.mean(observed)),
            "matched_pixel_count": int(len(points)),
            "rejected_column_count": int(trace["rejected_as_isolated"].sum()),
        }
    pd.concat(all_points, ignore_index=True).to_csv(
        intermediate / f"{panel}_pixel_points.csv", index=False
    )
    pd.concat(all_preclean, ignore_index=True).to_csv(
        intermediate / f"{panel}_preclean.csv", index=False
    )
    pd.concat(all_coordinates, ignore_index=True).to_csv(
        intermediate / f"{panel}_data_coordinates.csv", index=False
    )
    final.to_csv(regenerated / f"{panel}_digitized.csv", index=False)
    render_overlay(
        image_path,
        overlay_dir / f"{panel}_pixel_overlay.png",
        cleaned,
        FIG16_DETAIL_CURVES,
        calibration["plot_rectangle_pixels"],
    )
    return {
        "figure": panel,
        "source_image": str(image_path),
        "plot_rectangle_pixels": calibration["plot_rectangle_pixels"],
        "x_range": calibration["x_range"],
        "y_range": calibration["y_range"],
        "uncertainty": resolution_uncertainty(calibration),
        "curves": curve_metrics,
    }


def compare_detail_to_main(output_root: Path, details: dict[str, object]) -> dict[str, object]:
    regenerated = output_root / "digitization" / "regenerated"
    main = pd.read_csv(regenerated / "Fig16_digitized.csv")
    result: dict[str, object] = {}
    for panel in ("Fig16b", "Fig16c"):
        local = pd.read_csv(regenerated / f"{panel}_digitized.csv")
        panel_metrics: dict[str, object] = {}
        detail_metrics = details[panel]["curves"]  # type: ignore[index]
        for column in [spec.column for spec in FIG16_DETAIL_CURVES]:
            main_value = np.interp(local["time_s"], main["time_s"], main[column])
            residual = local[column].to_numpy() - main_value
            panel_metrics[column] = {
                "local_vs_main_overlap_rmse_C": float(np.sqrt(np.mean(residual**2))),
                "direct_valid_column_ratio": detail_metrics[column]["direct_valid_column_ratio"],  # type: ignore[index]
            }
        result[panel] = panel_metrics
    path = output_root / "digitization" / "detail_panel_comparison.json"
    path.write_text(json.dumps(result, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    return result


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--output-root", type=Path, default=PROJECT_ROOT)
    args = parser.parse_args()
    output_root = args.output_root.resolve()
    metrics = {name: digitize_figure(name, output_root) for name in FIGURE_CURVES}
    details = {panel: digitize_fig16_detail(panel, output_root) for panel in ("Fig16b", "Fig16c")}
    comparison = compare_detail_to_main(output_root, details)
    metrics.update(details)
    metrics["detail_vs_main"] = comparison
    report_path = output_root / "digitization" / "digitization_metrics.json"
    report_path.parent.mkdir(parents=True, exist_ok=True)
    report_path.write_text(json.dumps(metrics, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(metrics, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
