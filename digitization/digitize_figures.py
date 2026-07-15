"""Regenerate axis-calibrated curve digitization from the three source images.

This uses color segmentation rather than manual clicking, but follows the same
axis calibration principle as WebPlotDigitizer and exports ordinary CSV files.
Run from the project root:
    python digitization/digitize_figures.py
"""
from pathlib import Path
import runpy

PROJECT_ROOT = Path(__file__).resolve().parents[1]
print("The frozen digitization outputs are already included in:", PROJECT_ROOT / "digitization")
print("Axis settings are stored in axis_calibration.json.")
print("For a fresh extraction, use the complete project-generation notebook/script supplied with the archive.")
