# V15 decision

COMSOL 6.4 independently reproduced the frozen fixed-60% thermal equations.
The 1200-state BDF solve converged, COMSOL–Python maximum regional difference
was 0.00739 °C, and the mid-hold energy residual was 0.00996%.

The delivered Python entry point is
`models_v15.Fixed60HighResolutionReference` with 4800 states.  Its usable scope
is 23 °C ambient, 2–40 A, and fixed 60% PWM with measured current and voltage
inputs.  Independent spatial RMSE is 1.197 °C and fixed-duty regional dynamic
RMSE is 1.485 °C.

This does not overturn the negative frozen results: V13 missed the maximum
error and ordering gates, and V14 failed variable-PWM confirmation.  Therefore
the artifact is retained as an experimentally constrained high-resolution
reference.  It does not authorize ROM, MM-EKF, or Oracle MPC conclusions and
must not be described as a validated variable-duty high-fidelity plant.
