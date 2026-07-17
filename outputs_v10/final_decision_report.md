# V10 machine-linked decision report

- Numerical conservation and mesh gate: **PASS**.
- Frozen Python / independent COMSOL consistency: **PASS**.
- Restricted representative-channel cases 2–6: **PASS** (`Tmax RMSE = 0.940°C` in COMSOL).
- All-six-case full-cell temperature gate: **FAIL** (`Tmax RMSE = 7.968°C`; retained case 1 error `+19.403°C`).
- Full-cell spatial/hotspot gate: **NOT VALIDATED**.
- High-fidelity reference gate: **FAIL**.
- Decisive ROM, MM-EKF, and Oracle MPC authorization: **NO**.

See `docs_v10/final_decision_report.md` for the scientific interpretation and `outputs_v10_comsol/final_summary.json` for independent finite-element metrics.
