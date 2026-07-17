# V11 40-cell reference-model gate report

## Decision

V11 is **not accepted as a high-fidelity reference model**.  The first locked
confirmation run was retained as a negative result and its parameters must not
be retuned.

The model is a conservative 600-state thermal network (40 cells x 3 cathode
flow regions x 5 in-plane regions) based on the geometry and 75-thermocouple
experiment reported by Gong et al. (2024).  Fig. 8 groups 2/4 and Fig. 9(a)
were used for calibration.  Fig. 8 groups 1/5, Fig. 8 group 3, and Fig. 9(b)
were excluded from calibration.

## Locked results

| Gate | Threshold | Result | Pass |
|---|---:|---:|:---:|
| Independent spatial RMSE | <=2.0 °C | 3.370 °C | No |
| Independent Tmax absolute error | <=2.5 °C | 4.034 °C | No |
| Independent DeltaT absolute error | <=2.5 °C | 3.461 °C | No |
| Final spatial RMSE | <=2.0 °C | 2.652 °C | No |
| Final dynamic RMSE | <=2.5 °C | 4.570 °C | No |
| Final dynamic maximum error | <=4.0 °C | 12.034 °C | No |
| Hotspot-row accuracy | >=90% | 100% | Yes |
| Medium-fine sensor difference | <=0.2 °C | 0.800 °C | No |
| Energy residual | <=0.5% | <1e-12% | Yes |

## Structural diagnosis

The locked solution drives the stack-edge heat-transfer coefficient close to
its upper bound, intercell conductance close to its lower bound, and the fixed
source nonuniformity close to its upper bound.  It overstates the centre-edge
temperature difference: the two outer validation groups are underpredicted by
about 2.5--3.5 °C while the central held-out group is overpredicted by 2.3 °C
on average.

The 0.800 °C medium-fine discrepancy additionally shows that conductances were
parameterized per region instead of exclusively by face area and distance.
That formulation cannot be promoted to a mesh-independent reference model.

Dynamic errors reveal a second structural limitation: total cathode flow was
assumed linear in PWM duty although the experiment supplies no such mapping,
and stack voltage was represented by a smooth closure rather than the measured
time trace.  These effects cannot be repaired by changing the frozen V11
parameters.

## Required next model structure

The successor must use face-area-consistent finite-volume conductances,
explicit fan/plenum mixing, measured voltage for heat generation, a bounded
nonlinear PWM-to-flow relation, and separately auditable solid/contact thermal
mass.  The V11 confirmation data may be used only as a documented structural
failure signal; V11 itself remains frozen.

No ROM, MM-EKF, or MPC result is authorized from V11.

