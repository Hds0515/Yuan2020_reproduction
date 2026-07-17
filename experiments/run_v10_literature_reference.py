"""Run the frozen V10 literature-constrained representative-channel audit."""

from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy.stats import spearmanr


ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from models_v10.shahsavari_channel_reference import (  # noqa: E402
    TABLE1_CASES,
    frozen_parameter_record,
    solve_channel,
)


OUTPUT = ROOT / "outputs_v10"
MESHES = {"coarse": 50, "medium": 100, "fine": 200}
SOURCE_URL = (
    "https://www.sfu.ca/~mbahrami/pdf/2012/"
    "S.%20Shahsavari%2C%20A.%20Desouza%2C%20M.%20Bahrami%2C%20E.%20Kjeang%20-%20"
    "Thermal%20analysis%20of%20air-cooled%20PEM%20fuel%20cells.pdf"
)
SOURCE_PDF_SHA256 = "3E5CC71356AEAC3F726DCCB212D74842B59E7EB0DCF2BD21E20012539C47E768"


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def rmse(values: np.ndarray) -> float:
    return float(np.sqrt(np.mean(np.square(values))))


def run() -> dict[str, object]:
    OUTPUT.mkdir(parents=True, exist_ok=True)
    parameter_registry = frozen_parameter_record()
    parameter_registry["primary_source"] = {
        "citation": (
            "Shahsavari S, Desouza A, Bahrami M, Kjeang E. Thermal analysis of "
            "air-cooled PEM fuel cells. Int J Hydrogen Energy. 2012;37:18261-18271."
        ),
        "doi": "10.1016/j.ijhydene.2012.09.075",
        "url": SOURCE_URL,
        "verified_pdf_sha256": SOURCE_PDF_SHA256,
        "geometry_source": "Table 2",
        "validation_source": "Table 1",
    }
    parameter_registry["validation_protocol"] = {
        "parameters_fit_to_table1": False,
        "all_cases_full_cell_gate": [1, 2, 3, 4, 5, 6],
        "restricted_representative_channel_check": [2, 3, 4, 5, 6],
        "reason_case1_not_excluded_from_full_gate": (
            "Case 1 is retained in the full gate. The restricted check is reported "
            "separately because the paper states that its cell-level model used local "
            "current-density data and inlet/outlet plenums that are absent here."
        ),
        "full_gate_thresholds": {
            "Tmax_RMSE_C": 1.5,
            "maximum_absolute_Tmax_error_C": 2.0,
            "maximum_energy_residual_percent": 0.5,
            "medium_fine_Tmax_relative_percent": 1.0,
            "medium_fine_region_difference_C": 0.2,
            "full_cell_spatial_map_required": True,
        },
    }
    registry_path = OUTPUT / "frozen_parameter_registry.json"
    registry_path.write_text(
        json.dumps(parameter_registry, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    rows: list[dict[str, object]] = []
    profile_rows: list[dict[str, object]] = []
    for mesh_name, nodes in MESHES.items():
        for case in TABLE1_CASES:
            solution = solve_channel(case, nodes=nodes)
            row = solution.summary()
            row["mesh"] = mesh_name
            row["Tmax_error_vs_experiment_C"] = (
                solution.Tmax_C - case.experimental_Tmax_C
            )
            row["Tmax_error_vs_paper_simulation_C"] = (
                solution.Tmax_C - case.paper_simulated_Tmax_C
            )
            rows.append(row)
            if mesh_name == "fine":
                for x_m, solid_C, air_C in zip(
                    solution.x_m,
                    solution.solid_temperature_C,
                    solution.air_temperature_C,
                    strict=True,
                ):
                    profile_rows.append(
                        {
                            "case_id": case.case_id,
                            "x_normalized": x_m / solution.x_m[-1],
                            "solid_temperature_C": solid_C,
                            "air_temperature_C": air_C,
                        }
                    )

    results = pd.DataFrame(rows)
    results.to_csv(OUTPUT / "literature_validation.csv", index=False)
    profiles = pd.DataFrame(profile_rows)
    profiles.to_csv(OUTPUT / "temperature_profiles.csv", index=False)

    mesh_rows: list[dict[str, object]] = []
    region_columns = [
        "T_inlet_region_C",
        "T_middle_region_C",
        "T_outlet_region_C",
    ]
    for case in TABLE1_CASES:
        medium = results[(results.case_id == case.case_id) & (results.mesh == "medium")].iloc[0]
        fine = results[(results.case_id == case.case_id) & (results.mesh == "fine")].iloc[0]
        mesh_rows.append(
            {
                "case_id": case.case_id,
                "medium_fine_Tmax_difference_C": abs(medium.Tmax_C - fine.Tmax_C),
                "medium_fine_Tmax_relative_percent": (
                    abs(medium.Tmax_C - fine.Tmax_C) / fine.Tmax_C * 100.0
                ),
                "maximum_medium_fine_region_difference_C": max(
                    abs(float(medium[column]) - float(fine[column]))
                    for column in region_columns
                ),
            }
        )
    mesh = pd.DataFrame(mesh_rows)
    mesh.to_csv(OUTPUT / "mesh_convergence.csv", index=False)

    fine = results[results.mesh == "fine"].sort_values("case_id").copy()
    errors_all = fine.Tmax_error_vs_experiment_C.to_numpy(float)
    restricted = fine[fine.case_id.isin([2, 3, 4, 5, 6])]
    errors_restricted = restricted.Tmax_error_vs_experiment_C.to_numpy(float)
    paper_errors = (
        fine.paper_simulated_Tmax_C.to_numpy(float)
        - fine.experimental_Tmax_C.to_numpy(float)
    )
    trend = spearmanr(fine.experimental_Tmax_C, fine.Tmax_C)
    maximum_energy_residual_percent = float(
        fine.energy_residual_fraction.abs().max() * 100.0
    )
    maximum_mesh_Tmax_percent = float(mesh.medium_fine_Tmax_relative_percent.max())
    maximum_mesh_region_C = float(mesh.maximum_medium_fine_region_difference_C.max())
    full_metrics = {
        "Tmax_RMSE_C": rmse(errors_all),
        "maximum_absolute_Tmax_error_C": float(np.max(np.abs(errors_all))),
        "mean_absolute_Tmax_error_C": float(np.mean(np.abs(errors_all))),
        "Tmax_trend_Spearman": float(trend.statistic),
        "Tmax_trend_p_value": float(trend.pvalue),
    }
    restricted_metrics = {
        "case_ids": [2, 3, 4, 5, 6],
        "Tmax_RMSE_C": rmse(errors_restricted),
        "maximum_absolute_Tmax_error_C": float(np.max(np.abs(errors_restricted))),
        "mean_absolute_Tmax_error_C": float(np.mean(np.abs(errors_restricted))),
    }
    numerical_gate = bool(
        maximum_energy_residual_percent <= 0.5
        and maximum_mesh_Tmax_percent <= 1.0
        and maximum_mesh_region_C <= 0.2
    )
    full_temperature_gate = bool(
        full_metrics["Tmax_RMSE_C"] <= 1.5
        and full_metrics["maximum_absolute_Tmax_error_C"] <= 2.0
    )
    restricted_temperature_gate = bool(
        restricted_metrics["Tmax_RMSE_C"] <= 1.5
        and restricted_metrics["maximum_absolute_Tmax_error_C"] <= 2.0
    )

    summary: dict[str, object] = {
        "stage": "V10 literature-constrained representative-channel reference",
        "model_parameters_fitted_to_validation_cases": False,
        "numerical_metrics": {
            "maximum_energy_residual_percent": maximum_energy_residual_percent,
            "maximum_medium_fine_Tmax_relative_percent": maximum_mesh_Tmax_percent,
            "maximum_medium_fine_region_difference_C": maximum_mesh_region_C,
        },
        "all_six_case_full_cell_metrics": full_metrics,
        "restricted_representative_channel_metrics": restricted_metrics,
        "paper_reported_simulation_Tmax_RMSE_C": rmse(paper_errors),
        "gates": {
            "numerical_conservation_and_mesh": numerical_gate,
            "all_six_case_Tmax": full_temperature_gate,
            "restricted_cases_2_to_6_Tmax": restricted_temperature_gate,
            "full_cell_spatial_map": False,
            "full_cell_energy_partition": False,
            "full_cell_plenum_and_local_current_inputs_available": False,
            "high_fidelity_reference_validated": False,
            "authorize_decisive_ROM_MMEKF_MPC": False,
        },
        "decision": (
            "The frozen model is a numerically credible representative-channel core, "
            "but it is not a validated full-cell/manifold reference. Case 1 and the "
            "missing 24-point spatial map prevent the high-fidelity gate from passing."
        ),
        "case1_is_retained_as_failure": True,
        "parameter_registry_sha256": sha256(registry_path),
    }
    summary_path = OUTPUT / "final_summary.json"
    summary_path.write_text(
        json.dumps(summary, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    figure, axes = plt.subplots(1, 2, figsize=(11.5, 4.4))
    axes[0].plot(
        fine.experimental_Tmax_C,
        fine.Tmax_C,
        "o",
        color="#1f77b4",
        label="Frozen channel model",
    )
    limits = [30.0, 100.0]
    axes[0].plot(limits, limits, "--", color="0.35", label="1:1")
    for record in fine.itertuples():
        axes[0].annotate(str(record.case_id), (record.experimental_Tmax_C, record.Tmax_C))
    axes[0].set(xlim=limits, ylim=limits, xlabel="Experimental Tmax (°C)", ylabel="Model Tmax (°C)")
    axes[0].legend(frameon=False)
    axes[0].grid(alpha=0.25)

    for case_id in fine.case_id:
        subset = profiles[profiles.case_id == case_id]
        axes[1].plot(
            subset.x_normalized,
            subset.solid_temperature_C,
            label=f"Case {case_id}",
        )
    axes[1].set(
        xlabel="Normalized channel coordinate",
        ylabel="Solid temperature (°C)",
        xlim=(0.0, 1.0),
    )
    axes[1].grid(alpha=0.25)
    axes[1].legend(frameon=False, ncol=2)
    figure.tight_layout()
    figure.savefig(OUTPUT / "validation_overview.png", dpi=180)
    plt.close(figure)

    manifest_files = [
        "frozen_parameter_registry.json",
        "literature_validation.csv",
        "temperature_profiles.csv",
        "mesh_convergence.csv",
        "final_summary.json",
        "validation_overview.png",
    ]
    manifest = {name: sha256(OUTPUT / name) for name in manifest_files}
    (OUTPUT / "sha256_manifest.json").write_text(
        json.dumps(manifest, indent=2) + "\n",
        encoding="utf-8",
    )
    return summary


if __name__ == "__main__":
    print(json.dumps(run(), indent=2, ensure_ascii=False))
