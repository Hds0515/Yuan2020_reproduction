"""Select and freeze the single primary heat model using complexity-aware criteria."""

from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
import sys

import numpy as np
import pandas as pd
import yaml

PROJECT_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PROJECT_ROOT))

from identification_v3.run_identification import load_cases, simulate_both  # noqa: E402
from models.three_node_model import ThermalParameters  # noqa: E402


def canonical_hash(data: dict[str, object]) -> str:
    payload = json.dumps(data, sort_keys=True, ensure_ascii=False, separators=(",", ":"))
    return hashlib.sha256(payload.encode("utf-8")).hexdigest()


def load_candidate(path: Path) -> dict[str, object]:
    return json.loads(path.read_text(encoding="utf-8"))


def vector(candidate: dict[str, object]) -> np.ndarray:
    values = candidate["values"]
    names = ["C_th_J_K", "K_cool_W_K", "K_node_W_K", "Kp_1_K", "Ki_1_Ks"]
    if candidate["heat_model"] == "B_quadratic_current":
        names += ["heat_a0_W", "heat_a1_W_A", "heat_a2_W_A2"]
    return np.asarray([float(values[name]) for name in names])


def metrics(candidate: dict[str, object], cases: dict[str, object], config: dict[str, object]) -> dict[str, float | int | bool]:
    heat_model = str(candidate["heat_model"])
    result14, result16 = simulate_both(vector(candidate), heat_model, cases, config)
    mask14 = np.asarray(cases["fig14_time"]) <= 1200.0
    mask16 = np.asarray(cases["fig16_time"]) <= 700.0
    validation = np.asarray(cases["fig16_time"]) >= 700.0
    residual = np.concatenate(
        [
            (result14.temperature_C[mask14] - np.asarray(cases["fig14_target"])[mask14]).ravel(),
            (result16.temperature_C[mask16] - np.asarray(cases["fig16_target"])[mask16]).ravel(),
        ]
    )
    validation_residual = (
        result16.temperature_C[validation] - np.asarray(cases["fig16_target"])[validation]
    ).ravel()
    n = int(residual.size)
    k = 5 if heat_model == "A_thermoneutral_proxy" else 8
    rss = float(np.sum(residual**2))
    mse = max(rss / n, np.finfo(float).tiny)
    aic = n * math.log(mse) + 2 * k
    aicc = aic + (2 * k * (k + 1)) / max(n - k - 1, 1)
    bic = n * math.log(mse) + k * math.log(n)
    return {
        "physically_admissible": bool(np.all(vector(candidate) > 0)),
        "parameter_count": k,
        "calibration_observation_count": n,
        "calibration_rss_C2": rss,
        "calibration_rmse_C": float(np.sqrt(mse)),
        "AIC": float(aic),
        "AICc": float(aicc),
        "BIC": float(bic),
        "independent_validation_rmse_C": float(np.sqrt(np.mean(validation_residual**2))),
        "weighted_jacobian_condition_number": float(candidate.get("metrics", {}).get("weighted_jacobian_condition_number", np.nan)),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--root", type=Path, default=PROJECT_ROOT)
    parser.add_argument("--output-dir", type=Path)
    args = parser.parse_args()
    root = args.root.resolve()
    output = (args.output_dir or root / "identification_v4_primary").resolve()
    output.mkdir(parents=True, exist_ok=True)
    config = yaml.safe_load((root / "configs" / "reproduction_config.yaml").read_text(encoding="utf-8"))
    digitized = root / "outputs_v4" / "digitization" / "regenerated"
    if not digitized.exists():
        digitized = root / "digitization" / "regenerated"
    cases = load_cases(root, digitized)
    candidates = {
        "A_thermoneutral_proxy": load_candidate(root / "identification_v3" / "frozen_parameters_v3.json"),
        "B_quadratic_current": load_candidate(root / "outputs_v4" / "identification" / "identified_parameters_v4.json"),
    }
    rows = []
    for name, candidate in candidates.items():
        row = {"heat_model": name, **metrics(candidate, cases, config)}
        rows.append(row)
    frame = pd.DataFrame(rows).sort_values("AICc")
    best_aicc = float(frame.iloc[0].AICc)
    frame["delta_AICc"] = frame.AICc - best_aicc
    frame["delta_BIC_from_best"] = frame.BIC - float(frame.BIC.min())

    admissible = frame[frame.physically_admissible]
    if admissible.empty:
        raise RuntimeError("No physically admissible candidate model")
    top = admissible.sort_values("AICc").iloc[0]
    near = admissible[admissible.AICc <= float(top.AICc) + 2.0]
    if len(near) > 1:
        selected_row = near.sort_values(["parameter_count", "BIC"]).iloc[0]
        rule = "delta_AICc<2: selected the lower-parameter candidate"
    else:
        selected_row = top
        rule = "selected minimum AICc among physically admissible candidates"
    # Validation is not used to fit; it is an external sanity check.  If the more complex
    # candidate is worse on validation and only marginally better in calibration, keep A.
    model_a = frame[frame.heat_model == "A_thermoneutral_proxy"].iloc[0]
    model_b = frame[frame.heat_model == "B_quadratic_current"].iloc[0]
    if (
        float(model_b.independent_validation_rmse_C) > float(model_a.independent_validation_rmse_C)
        and float(model_b.calibration_rmse_C) >= float(model_a.calibration_rmse_C) - 0.005
    ):
        selected_row = model_a
        rule += "; external validation and parsimony override selected model A"

    selected_name = str(selected_row.heat_model)
    selected = json.loads(json.dumps(candidates[selected_name]))
    selected["model"] = "Yuan2020_three_node_primary_v5"
    selected["selection"] = {
        "selected_heat_model": selected_name,
        "selection_rule": rule,
        "validation_not_used_for_parameter_fitting": True,
        "candidate_metrics": frame.to_dict(orient="records"),
        "sensitivity_model": "B_quadratic_current" if selected_name == "A_thermoneutral_proxy" else "A_thermoneutral_proxy",
    }
    selected.pop("sha256_without_hash_field", None)
    selected["canonical_content_sha256"] = canonical_hash(selected)
    target = output / "frozen_primary_model.json"
    target.write_text(json.dumps(selected, ensure_ascii=False, indent=2) + "\n", encoding="utf-8")
    file_hash = hashlib.sha256(target.read_bytes()).hexdigest()
    frame["selected_primary"] = frame.heat_model == selected_name
    frame.to_csv(output / "model_selection_metrics.csv", index=False)
    report = f"""# Primary model selection\n\n- Selected model: **{selected_name}**\n- Rule: {rule}\n- File SHA-256: `{file_hash}`\n- Validation data were not used for fitting.\n\n{frame.to_markdown(index=False)}\n\nThe quadratic-current model remains a heat-source sensitivity analysis, not the downstream primary model.\n"""
    (output / "model_selection_report.md").write_text(report, encoding="utf-8")
    (output / "primary_model_file_sha256.txt").write_text(file_hash + "\n", encoding="utf-8")
    print(json.dumps({"selected": selected_name, "sha256": file_hash, "rule": rule}, indent=2))


if __name__ == "__main__":
    main()
