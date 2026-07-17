# V12 gate report

## Outcome

V12 resolves the steady spatial field but fails the dynamic portability gate.
Its frozen parameters must not be retuned.

The delivered finite-volume discretization has 4800 thermal states and was
checked against a 19200-state mesh.  The held-out central spatial group passed
all registered temperature gates:

- spatial RMSE: 1.197 °C;
- Tmax absolute error: 1.336 °C;
- DeltaT absolute error: 1.983 °C;
- hotspot-row accuracy: 100%;
- medium-fine maximum sensor difference: 0.068 °C;
- steady energy residual: below 1e-12%.

The independent Fig. 9(b) random-load trace failed with RMSE 8.265 °C and
maximum error 14.257 °C.  The model used a fan-duty history reconstructed from
current because neither Gong et al. (2024) article publishes PWM telemetry for
the two Fig. 9 experiments.  The reconstruction calibrated on Fig. 9(a) did
not transfer to Fig. 9(b).  Early in Fig. 9(b), the measured stack cools from
30.5 °C to about 25.5 °C while the inferred minimum fan command leaves the
model near 31 °C.  This is direct evidence that current alone does not identify
the cooling input.

## Scientific decision

The spatial/steady finite-volume structure is retained as a validated
component.  V12 as a complete transient reference model is rejected.  The
successor must validate dynamics against an experiment with explicit fan PWM
or mass flow.  Guo et al. (2024), DOI 10.1016/j.ijhydene.2024.01.045, uses the
same 75-thermocouple stack platform and reports ramp experiments at fixed 60%
PWM; those experiments remove the V12 input ambiguity.

ROM, observer, and MPC conclusions remain unauthorized at this checkpoint.

