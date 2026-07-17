# V13 final gate report

V13 is retained as a near-pass negative checkpoint.  It corrected the energy
ledger to use 1.48 V/cell thermoneutral voltage and independently reproduced
the fixed-60%-PWM 4 A ramp with overall regional RMSE 1.485 °C.  Each regional
RMSE was below 1.63 °C, peak-time error was 0 s, the 4800-to-19200-state
difference was 0.121 °C, and the energy residual was below 1e-12%.

Two preregistered gates narrowly failed:

- maximum confirmation error was 4.271 °C (limit 4.0 °C);
- regional ordering accuracy was 94.74% (limit 95%).

The total heat capacity approached its upper bound and surface-capacity
fraction approached its lower bound.  This indicates that the two-capacitance
structure is compensating for unresolved distributed stack/fixture thermal
mass.  The model is not retuned after the confirmation result and is not yet
the final high-fidelity deliverable.

The next version must use the independent multi-PWM trajectories in Guo et al.
Fig. 11 to test whether fan-command changes transfer.  A result that passes
only fixed PWM is not sufficient for variable-duty control research.

