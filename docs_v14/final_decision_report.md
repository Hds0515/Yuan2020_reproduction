# V14 multi-PWM model decision

## Frozen test

V14 estimated only `fan_flow_exponent_relative_to_60pct` from the PWM40 trace.
The resulting frozen JSON has SHA-256
`7cd86ccb0e23b0b12f409c0a9bc9c910fc04cd878f115b1a5cf2a1ebd1ec447d`.
PWM60 was an independent validation trajectory and PWM80 was opened once as the
final confirmation trajectory.  No parameter was changed afterward.

## Outcome

The model reproduced the monotonic temperature trend, mesh convergence, and
the numerical energy ledger, but it missed the independent PWM60 and PWM80
temperature trajectories.  The final-confirmation RMSE was 3.055 °C and its
maximum absolute error was 6.111 °C.  Therefore the model is not a validated
variable-PWM high-fidelity plant.

The discrepancy is consistent with omitted duty-dependent fan behavior,
pressure-flow interaction, and water-management effects.  Those mechanisms
cannot be identified uniquely from average-temperature curves alone.  They are
recorded as a model boundary instead of being fitted post hoc.

## Scope decision

The fixed-60% PWM experimental evidence and the independent spatial field are
retained for a narrower confirmation study.  Until that independent
implementation check is complete, the model must not be used to authorize ROM,
MM-EKF, or decisive MPC experiments.  The 55 °C quantity is a reference, not a
safety limit.
