"""Image-only curve extraction utilities for Yuan et al. (2020) figures."""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from PIL import Image


@dataclass(frozen=True)
class CurveSpec:
    column: str
    rgb: tuple[int, int, int]
    threshold: float = 48.0


def pixel_to_data(
    pixel_x: np.ndarray,
    pixel_y: np.ndarray,
    rectangle: list[int],
    x_range: list[float],
    y_range: list[float],
) -> tuple[np.ndarray, np.ndarray]:
    x0, y0, x1, y1 = rectangle
    x = x_range[0] + (pixel_x - x0) * (x_range[1] - x_range[0]) / (x1 - x0)
    y = y_range[1] - (pixel_y - y0) * (y_range[1] - y_range[0]) / (y1 - y0)
    return x, y


def data_to_pixel(
    x: np.ndarray,
    y: np.ndarray,
    rectangle: list[int],
    x_range: list[float],
    y_range: list[float],
) -> tuple[np.ndarray, np.ndarray]:
    x0, y0, x1, y1 = rectangle
    px = x0 + (x - x_range[0]) * (x1 - x0) / (x_range[1] - x_range[0])
    py = y0 + (y_range[1] - y) * (y1 - y0) / (y_range[1] - y_range[0])
    return px, py


def _weighted_row(rows: np.ndarray, distances: np.ndarray) -> float:
    """Return a robust color-core row, retaining line thickness as evidence."""
    if len(rows) == 1:
        return float(rows[0])
    weights = 1.0 / np.maximum(distances, 1.0)
    order = np.argsort(rows)
    ordered_rows = rows[order]
    cumulative = np.cumsum(weights[order])
    return float(ordered_rows[np.searchsorted(cumulative, cumulative[-1] / 2.0)])


def extract_curve(
    image_rgb: np.ndarray,
    spec: CurveSpec,
    rectangle: list[int],
    retained_y_pixels: list[int],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Extract all color-matched pixels and a robust per-column raw trace."""
    x0, _, x1, _ = rectangle
    keep_y0, keep_y1 = retained_y_pixels
    crop = image_rgb[keep_y0 : keep_y1 + 1, x0 : x1 + 1].astype(float)
    target = np.asarray(spec.rgb, dtype=float)
    distance = np.linalg.norm(crop - target, axis=2)
    saturation = crop.max(axis=2) - crop.min(axis=2)
    mask = (distance <= spec.threshold) & (saturation >= 35.0)

    point_rows, point_cols = np.nonzero(mask)
    points = pd.DataFrame(
        {
            "curve": spec.column,
            "pixel_x": point_cols + x0,
            "pixel_y": point_rows + keep_y0,
            "color_distance": distance[point_rows, point_cols],
        }
    )

    trace_rows: list[dict[str, float | int | str]] = []
    for local_x in range(mask.shape[1]):
        rows = np.flatnonzero(mask[:, local_x])
        if rows.size == 0:
            trace_rows.append(
                {
                    "curve": spec.column,
                    "pixel_x": local_x + x0,
                    "pixel_y_raw": np.nan,
                    "matched_pixel_count": 0,
                    "minimum_color_distance": np.nan,
                }
            )
            continue
        dist = distance[rows, local_x]
        trace_rows.append(
            {
                "curve": spec.column,
                "pixel_x": local_x + x0,
                "pixel_y_raw": _weighted_row(rows + keep_y0, dist),
                "matched_pixel_count": int(rows.size),
                "minimum_color_distance": float(np.min(dist)),
            }
        )
    return points, pd.DataFrame(trace_rows)


def clean_trace(trace: pd.DataFrame, max_isolated_deviation_px: float = 9.0) -> pd.DataFrame:
    """Reject isolated color artifacts, then interpolate only missing columns."""
    result = trace.copy()
    y = result["pixel_y_raw"].to_numpy(dtype=float)
    series = pd.Series(y)
    local = series.rolling(9, center=True, min_periods=3).median().to_numpy()
    finite = np.isfinite(y)
    isolated = finite & np.isfinite(local) & (np.abs(y - local) > max_isolated_deviation_px)
    # Preserve sustained ramps/steps: reject only when both immediate neighbors follow the local trace.
    prev_ok = np.r_[False, np.abs(y[:-1] - local[:-1]) <= max_isolated_deviation_px]
    next_ok = np.r_[np.abs(y[1:] - local[1:]) <= max_isolated_deviation_px, False]
    rejected = isolated & prev_ok & next_ok
    y[rejected] = np.nan
    valid = np.isfinite(y)
    if valid.sum() < 2:
        raise RuntimeError(f"Curve {trace['curve'].iloc[0]} has fewer than two valid columns")
    x = result["pixel_x"].to_numpy(dtype=float)
    result["pixel_y_clean"] = np.interp(x, x[valid], y[valid])
    result["rejected_as_isolated"] = rejected
    result["observed_before_interpolation"] = valid
    return result


def render_overlay(
    image_path: Path,
    output_path: Path,
    cleaned: dict[str, pd.DataFrame],
    specs: list[CurveSpec],
    rectangle: list[int],
) -> None:
    image = np.asarray(Image.open(image_path).convert("RGB"))
    fig, ax = plt.subplots(figsize=(12, 8), dpi=160)
    ax.imshow(image)
    for spec in specs:
        trace = cleaned[spec.column]
        ax.plot(
            trace["pixel_x"],
            trace["pixel_y_clean"],
            color=np.asarray(spec.rgb) / 255.0,
            linewidth=1.0,
            alpha=0.9,
            label=spec.column,
        )
    x0, y0, x1, y1 = rectangle
    ax.plot([x0, x1, x1, x0, x0], [y0, y0, y1, y1, y0], "k--", lw=0.8)
    ax.set_xlim(max(0, x0 - 20), min(image.shape[1], x1 + 20))
    ax.set_ylim(min(image.shape[0], y1 + 20), max(0, y0 - 20))
    ax.set_title(f"Image-only digitization overlay: {image_path.stem}")
    ax.legend(fontsize=7, ncol=2, loc="lower left")
    ax.set_axis_off()
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fig.tight_layout()
    fig.savefig(output_path, bbox_inches="tight")
    plt.close(fig)


def resolution_uncertainty(config: dict[str, object]) -> dict[str, float]:
    x0, y0, x1, y1 = config["plot_rectangle_pixels"]
    xmin, xmax = config["x_range"]
    ymin, ymax = config["y_range"]
    return {
        "x_units_per_pixel": (xmax - xmin) / (x1 - x0),
        "y_units_per_pixel": (ymax - ymin) / (y1 - y0),
        "x_half_pixel_uncertainty": 0.5 * (xmax - xmin) / (x1 - x0),
        "y_half_pixel_uncertainty": 0.5 * (ymax - ymin) / (y1 - y0),
    }
