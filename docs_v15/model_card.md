# V15 fixed-60% high-resolution reference model

## Delivered entry point

```python
from models_v15 import Fixed60HighResolutionReference

plant = Fixed60HighResolutionReference()
plant.assert_in_scope(ambient_C=23.0, duty=0.60, current_A=30.0)
```

The delivered mesh contains 40 stack positions, 6 cathode-flow subdivisions,
10 in-plane subdivisions, and two thermal masses at every finite volume.  It
therefore has 4800 temperature states.  A 19200-state refinement was used for
the mesh comparison.

## Validity domain

- ambient temperature: 23 °C;
- current: 2–40 A;
- fan command: fixed at 60% PWM;
- measured current and stack voltage are exogenous inputs;
- forward cathode flow only;
- thermal states only—membrane water, liquid flooding, voltage, degradation,
  fan pressure-flow dynamics, and reversal are not predicted.

Calling `assert_in_scope` with a different fan duty or ambient temperature
raises an error.  This is intentional: V14 showed that a single PWM-to-flow
law does not transfer across PWM40, PWM60, and PWM80.

## Evidence

| Check | Result |
|---|---:|
| Independent spatial RMSE | 1.197 °C |
| Fixed-60% regional dynamic RMSE | 1.485 °C |
| Worst regional RMSE | 1.621 °C |
| Dynamic maximum absolute error | 4.271 °C |
| Regional-order accuracy | 94.737% |
| Medium–fine sensor difference | 0.121 °C |
| COMSOL–Python maximum difference | 0.00739 °C |
| COMSOL energy residual, mid-hold maximum | 0.00996% |

The spatial, aggregate dynamic, mesh, energy, and implementation-consistency
checks pass.  The preregistered V13 maximum-error limit (4.0 °C) and ordering
limit (95%) were narrowly missed.  Consequently this artifact is an
**experimentally constrained high-resolution reference**, not a fully accepted
variable-duty high-fidelity control plant.

## Provenance

The 40-cell geometry and 75-thermocouple layout come from Gong et al.,
*International Journal of Hydrogen Energy* 60 (2024), 1134–1146,
DOI `10.1016/j.ijhydene.2024.02.270`.  The fixed-60% ramp evidence comes from
Guo et al., *International Journal of Hydrogen Energy* 57 (2024), 601–612,
DOI `10.1016/j.ijhydene.2024.01.045`.  The thermoneutral heat ledger uses
1.48 V per cell.  All fitted values are inherited unchanged from the frozen
V13 record.

The frozen V15 reference JSON SHA-256 is
`05566ac5903af63fdda10c1dd941a3bfbd6b7f607409007d341e320debf81ae7`.
The 55 °C quantity remains a reference temperature, not a safety limit.
