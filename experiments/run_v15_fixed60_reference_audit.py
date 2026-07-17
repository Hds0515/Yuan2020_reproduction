"""Audit COMSOL/Python consistency and package the V15 scoped reference."""

from __future__ import annotations

import hashlib
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
if str(ROOT) not in sys.path:
    sys.path.insert(0, str(ROOT))

from experiments.run_v13_fixed_pwm_confirmation import load, metrics, prediction  # noqa: E402
from models_v15 import Fixed60HighResolutionReference  # noqa: E402


OUT = ROOT / "outputs_v15"


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest().lower()


def read_comsol() -> pd.DataFrame:
    names = [
        "time_s",
        "inlet_C",
        "intermediate_C",
        "outlet_C",
        "regional_max_C",
        "heat_input_W",
        "heat_rejection_W",
        "stored_heat_rate_W",
        "energy_residual_percent",
    ]
    return pd.read_csv(
        OUT / "comsol_fixed60_regional_timeseries.csv",
        comment="%",
        header=None,
        names=names,
    )


def main() -> None:
    comsol = read_comsol()
    times, target, python_regions = prediction(
        load(), "ramp_4A", flow_rows=3, columns=5
    )
    comsol_regions = (
        comsol.set_index("time_s")
        .loc[times, ["inlet_C", "intermediate_C", "outlet_C"]]
        .to_numpy(float)
    )
    difference = comsol_regions - python_regions

    rows: list[dict[str, object]] = []
    region_names = ("inlet", "intermediate", "outlet")
    for ti, time_s in enumerate(times):
        for ri, region in enumerate(region_names):
            rows.append(
                {
                    "time_s": time_s,
                    "region": region,
                    "COMSOL_C": comsol_regions[ti, ri],
                    "Python_C": python_regions[ti, ri],
                    "difference_C": difference[ti, ri],
                }
            )
    pd.DataFrame(rows).to_csv(OUT / "comsol_vs_python.csv", index=False)

    delivered = Fixed60HighResolutionReference()
    final_confirmation = prediction(load(), "ramp_4A", flow_rows=6, columns=10)
    experimental = metrics(*final_confirmation)
    voltage_40A = float(np.clip(33.5 - 0.34 * 40.0, 19.5, 33.5))
    surface, core = delivered.steady_state(
        40.0, voltage_40A, duty=0.60, ambient_C=23.0
    )
    python_energy = delivered.energy_audit(
        surface, core, 40.0, voltage_40A, 0.60, 23.0
    )

    # Exact event times have a discontinuous source and two valid one-sided
    # derivatives. Energy scoring therefore uses the mid-hold samples (50 mod
    # 100 s), which are also the digitized experimental comparison times.
    energy_scoring = comsol.loc[(comsol.time_s % 100.0) == 50.0]
    comsol_energy_max = float(energy_scoring.energy_residual_percent.abs().max())

    record = delivered.validation_record()
    record.update(
        {
            "created_utc": datetime.now(timezone.utc).isoformat(),
            "python_model_source_sha256": sha256(
                ROOT / "models_v15/fixed60_high_resolution_reference.py"
            ),
            "comsol_model_source_sha256": sha256(
                ROOT / "comsol/v15_reference/V15Fixed60NetworkConfirmation.java"
            ),
            "comsol_runtime_model_sha256": sha256(
                ROOT / "comsol/v15_reference/v15_fixed60_network_confirmation.mph"
            ),
            "comsol_version": "6.4.0.293",
            "parameter_changes_from_V13": 0,
        }
    )
    frozen = OUT / "frozen_high_resolution_reference.json"
    frozen.write_text(json.dumps(record, indent=2) + "\n", encoding="utf-8")
    frozen_sha = sha256(frozen)

    consistency_checks = {
        "regional_maximum_difference": float(np.max(np.abs(difference))) <= 0.05,
        "regional_RMSE": float(np.sqrt(np.mean(np.square(difference)))) <= 0.05,
        "COMSOL_energy_residual": comsol_energy_max <= 0.5,
    }
    result = {
        "frozen_reference_sha256": frozen_sha,
        "parent_V13_sha256": record["frozen_parent_sha256"],
        "COMSOL": {
            "version": "6.4.0.293",
            "state_count": 1200,
            "transient_solver": "BDF",
            "converged": True,
            "reported_core_solve_time_s": 7.09266290,
            "COMSOL_vs_Python_maximum_difference_C": float(np.max(np.abs(difference))),
            "COMSOL_vs_Python_RMSE_C": float(np.sqrt(np.mean(np.square(difference)))),
            "COMSOL_energy_residual_max_mid_hold_percent": comsol_energy_max,
            "checks": consistency_checks,
            "implementation_consistency_passed": bool(all(consistency_checks.values())),
        },
        "delivered_model": {
            "entry_point": "models_v15.Fixed60HighResolutionReference",
            "state_count": delivered.state_count,
            "fine_comparison_state_count": 19200,
            "scope": record["scope"],
        },
        "experimental_evidence": {
            "independent_spatial_RMSE_C": 1.1972718732972596,
            "fixed60_regional_dynamic": experimental,
            "medium_fine_max_sensor_difference_C": 0.12070005118104632,
            "steady_energy_audit": python_energy,
        },
        "decisions": {
            "scoped_high_resolution_reference_available": True,
            "COMSOL_implementation_confirmed": bool(all(consistency_checks.values())),
            "strict_experimental_high_fidelity_gate_passed": False,
            "variable_PWM_high_fidelity_plant_available": False,
            "ROM_MM_EKF_oracle_MPC_authorized": False,
            "reason": "V13 maximum-error and regional-order thresholds and V14 variable-PWM confirmation remain failed.",
        },
        "reference_temperature_C": 55.0,
        "reference_temperature_is_not_a_safety_limit": True,
    }
    (OUT / "validation_metrics.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    (OUT / "final_summary.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
