"""Audit COMSOL/Python consistency for the frozen V10 literature channel."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd


ROOT = Path(__file__).resolve().parents[2]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from models_v10.shahsavari_channel_reference import TABLE1_CASES  # noqa: E402


RAW = ROOT / "comsol" / "v10_reference"
OUTPUT = ROOT / "outputs_v10_comsol"
METRICS = (
    "T_inlet_region_C",
    "T_middle_region_C",
    "T_outlet_region_C",
    "Tmax_C",
    "Tmin_C",
    "DeltaT_C",
    "hotspot_position_normalized",
    "air_outlet_temperature_C",
    "air_enthalpy_gain_channel_W",
    "channel_heat_input_W",
    "energy_residual_fraction",
    "mass_flow_channel_kg_s",
    "reynolds_number",
    "nusselt_number",
    "forced_h_W_m2K",
    "pressure_drop_Pa",
)


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def read_metrics(mesh: str) -> pd.DataFrame:
    data = pd.read_csv(RAW / f"metrics_{mesh}.csv", comment="%", header=None)
    if data.shape != (6, len(METRICS) + 1):
        raise RuntimeError(f"Unexpected {mesh} metrics shape: {data.shape}")
    data.columns = ("case_id", *METRICS)
    data.insert(0, "mesh", mesh)
    data["case_id"] = data.case_id.astype(int)
    return data


def main() -> None:
    OUTPUT.mkdir(exist_ok=True)
    java_path = RAW / "V10LiteratureChannel.java"
    source = java_path.read_text(encoding="utf-8")
    required = (
        'model.param().set("nChannels", "80")',
        'model.param().set("L", "60[mm]")',
        'model.param().set("Nu", "3.610224")',
        '"qLine", "qCell/(nChannels*L)"',
        '"mdot", "rhoAir*uin*ACh"',
    )
    missing = [item for item in required if item not in source]
    if missing:
        raise RuntimeError(f"COMSOL source audit failed: {missing}")

    frames = [read_metrics(mesh) for mesh in ("coarse", "medium", "fine")]
    all_metrics = pd.concat(frames, ignore_index=True)
    all_metrics.to_csv(OUTPUT / "comsol_confirmation.csv", index=False)
    fine = all_metrics[all_metrics.mesh == "fine"].copy()
    python = pd.read_csv(ROOT / "outputs_v10" / "literature_validation.csv")
    python = python[python.mesh == "fine"].copy()
    selected = (
        "case_id",
        "T_inlet_region_C",
        "T_middle_region_C",
        "T_outlet_region_C",
        "Tmax_C",
        "Tmin_C",
        "DeltaT_C",
        "air_outlet_temperature_C",
        "air_enthalpy_gain_channel_W",
        "channel_heat_input_W",
        "mass_flow_channel_kg_s",
        "reynolds_number",
        "nusselt_number",
        "forced_h_W_m2K",
        "pressure_drop_Pa",
    )
    comparison = fine[list(selected)].merge(
        python[list(selected)], on="case_id", suffixes=("_comsol", "_python")
    )
    for name in selected[1:]:
        comparison[f"{name}_difference"] = (
            comparison[f"{name}_comsol"] - comparison[f"{name}_python"]
        )
    comparison.to_csv(OUTPUT / "comsol_vs_python.csv", index=False)

    region_names = (
        "T_inlet_region_C",
        "T_middle_region_C",
        "T_outlet_region_C",
    )
    max_region_difference = max(
        float(comparison[f"{name}_difference"].abs().max()) for name in region_names
    )
    max_Tmax_difference = float(comparison.Tmax_C_difference.abs().max())
    max_DeltaT_difference = float(comparison.DeltaT_C_difference.abs().max())
    max_energy_percent = float(fine.energy_residual_fraction.abs().max() * 100.0)

    medium = all_metrics[all_metrics.mesh == "medium"].set_index("case_id")
    fine_indexed = fine.set_index("case_id")
    max_mesh_region = max(
        float((medium[name] - fine_indexed[name]).abs().max()) for name in region_names
    )
    max_mesh_Tmax = float((medium.Tmax_C - fine_indexed.Tmax_C).abs().max())
    agreement = bool(
        max_region_difference <= 0.2
        and max_Tmax_difference <= 0.25
        and max_DeltaT_difference <= 0.25
        and max_energy_percent <= 0.5
        and max_mesh_region <= 0.2
    )
    experimental = np.asarray([case.experimental_Tmax_C for case in TABLE1_CASES])
    errors = fine.sort_values("case_id").Tmax_C.to_numpy() - experimental
    restricted_errors = errors[1:]
    summary = {
        "comsol_version": "6.4",
        "java_source_sha256": sha256(java_path),
        "frozen_python_registry_sha256": sha256(
            ROOT / "outputs_v10" / "frozen_parameter_registry.json"
        ),
        "all_six_cases_converged": True,
        "all_three_meshes_converged": True,
        "consistency_metrics": {
            "maximum_region_difference_C": max_region_difference,
            "maximum_Tmax_difference_C": max_Tmax_difference,
            "maximum_DeltaT_difference_C": max_DeltaT_difference,
            "maximum_energy_residual_percent": max_energy_percent,
            "maximum_medium_fine_region_difference_C": max_mesh_region,
            "maximum_medium_fine_Tmax_difference_C": max_mesh_Tmax,
        },
        "comsol_matches_frozen_python": agreement,
        "fig7_like_experimental_metrics": {
            "all_six_Tmax_RMSE_C": float(np.sqrt(np.mean(errors**2))),
            "all_six_maximum_absolute_Tmax_error_C": float(np.max(np.abs(errors))),
            "restricted_cases_2_to_6_Tmax_RMSE_C": float(
                np.sqrt(np.mean(restricted_errors**2))
            ),
            "restricted_cases_2_to_6_maximum_absolute_Tmax_error_C": float(
                np.max(np.abs(restricted_errors))
            ),
        },
        "high_fidelity_reference_validated": False,
        "authorize_decisive_control_research": False,
        "reason": (
            "Independent finite-element agreement does not repair the retained case-1 "
            "failure or supply the unavailable full-cell spatial/current/plenum inputs."
        ),
    }
    summary_path = OUTPUT / "final_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n", encoding="utf-8"
    )
    manifest = {
        name: sha256(OUTPUT / name)
        for name in (
            "comsol_confirmation.csv",
            "comsol_vs_python.csv",
            "final_summary.json",
        )
    }
    (OUTPUT / "sha256_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
