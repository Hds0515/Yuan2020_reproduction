"""Build the required machine-readable final status from actual artifacts."""

from __future__ import annotations

import csv
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent


def main() -> None:
    identification = json.loads((ROOT / "identification_v3" / "identification_metrics.json").read_text(encoding="utf-8"))
    five = json.loads((ROOT / "outputs" / "five_node" / "consistency_metrics.json").read_text(encoding="utf-8"))
    observer = json.loads((ROOT / "outputs" / "observer" / "observer_summary.json").read_text(encoding="utf-8"))
    control = json.loads((ROOT / "outputs" / "control" / "control_summary.json").read_text(encoding="utf-8"))
    comsol = json.loads((ROOT / "comsol" / "status.json").read_text(encoding="utf-8"))
    digitization_files = [ROOT / "digitization" / "regenerated" / f"Fig{number}_digitized.csv" for number in (14, 15, 16)]
    summary = {
        "curve_digitization_reproducible": all(path.exists() and path.stat().st_size > 0 for path in digitization_files),
        "three_node_reproduction_passed": bool(identification["fig14_calibration_rmse_C"] <= 0.35 and identification["fig16_independent_validation_rmse_C"] <= 0.40 and identification["maximum_relative_energy_residual"] < 1e-6),
        "independent_validation_rmse_C": float(identification["fig16_independent_validation_rmse_C"]),
        "comsol_compiled": bool(comsol["java_compile_success"]),
        "comsol_solved": bool(comsol["stationary_solution_converged"]),
        "mesh_convergence_passed": bool(comsol["mesh_independence_completed"]),
        "energy_balance_passed": bool(comsol["energy_balance_completed"]),
        "comsol_three_node_cross_validation_passed": bool(comsol["three_node_cross_validation_completed"]),
        "five_node_energy_consistency_passed": bool(five["maximum_absolute_energy_residual_W"] < 1e-6 and abs(five["area_preservation_ratio"] - 1.0) < 1e-12),
        "ekf_completed": (ROOT / "outputs" / "observer" / "sensor_layout_metrics.csv").exists(),
        "multiple_model_observer_completed": bool(observer.get("multiple_model_improves_mean") is not None),
        "mpc_completed": bool(control.get("scenario_count", 0) >= 5),
        "remaining_limitations": [
            "COMSOL Java compiles, but the stationary solve did not converge; no MPH result, mesh convergence, CFD energy balance, or cross-validation may be claimed.",
            "The paper does not publish full CAD/channel geometry; the COMSOL model is explicitly equivalent, not exact.",
            "The polarization curve is a documented proxy rather than the authors' full electrochemical model.",
            "K_cool, Kp, and Ki are strongly entangled without measured fan duty; product parameters are more defensible than individual values.",
            "The current bootstrap count is a reproducibility smoke run (8); use at least 200 resamples for publication-grade intervals.",
            "MPC conclusions are simulation-only and require experimental or converged CFD corroboration.",
        ],
    }
    output = ROOT / "outputs"; output.mkdir(exist_ok=True)
    (output / "final_summary.json").write_text(json.dumps(summary, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    with (output / "final_summary.csv").open("w", newline="", encoding="utf-8-sig") as handle:
        writer = csv.writer(handle); writer.writerow(["item", "value"])
        for key, value in summary.items():
            writer.writerow([key, json.dumps(value, ensure_ascii=False) if isinstance(value, (list, dict)) else value])
    print(json.dumps(summary, ensure_ascii=False, indent=2))


if __name__ == "__main__":
    main()
