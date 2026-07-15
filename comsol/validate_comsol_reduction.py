"""Validate the frozen three-node reduction against COMSOL regional temperatures.

1. Run Yuan2020Equivalent3D.java in COMSOL 6.4.
2. Convert the COMSOL table to:
   uin_m_s,T_zone1_C,T_zone2_C,T_zone3_C,Tmax_C
3. Save it as comsol_zone_temperatures.csv in this folder.
4. Run:
   python validate_comsol_reduction.py
"""
from pathlib import Path
import json
import numpy as np
import pandas as pd

HERE = Path(__file__).resolve().parent
input_path = HERE / "comsol_zone_temperatures.csv"
if not input_path.exists():
    raise FileNotFoundError(
        "COMSOL result is not present. Copy the completed template to "
        "comsol_zone_temperatures.csv after running COMSOL."
    )

data = pd.read_csv(input_path)
required = ["uin_m_s", "T_zone1_C", "T_zone2_C", "T_zone3_C", "Tmax_C"]
missing = [column for column in required if column not in data.columns]
if missing:
    raise ValueError(f"Missing columns: {missing}")

# A first-order normalized mapping between CFD inlet speed and fan command.
# Replace this with the measured Fig. 2 fan map when raw fan data become available.
data["normalized_fan_command"] = np.clip(data["uin_m_s"] / 12.0, 0.0, 1.0)

zone = data[["T_zone1_C", "T_zone2_C", "T_zone3_C"]].to_numpy()
gradient = np.max(zone, axis=1) - np.min(zone, axis=1)
summary = {
    "number_of_operating_points": int(len(data)),
    "maximum_COMSOL_zone_gradient_C": float(np.max(gradient)),
    "mean_COMSOL_zone_gradient_C": float(np.mean(gradient)),
    "note": (
        "This script performs the COMSOL-side regional data audit. "
        "Dynamic three-node comparison should use the exported regional "
        "temperature transients if a time-dependent COMSOL study is added."
    ),
}
(HERE.parent / "outputs" / "comsol_reduction_validation.json").write_text(
    json.dumps(summary, indent=2), encoding="utf-8"
)
print(json.dumps(summary, indent=2))
