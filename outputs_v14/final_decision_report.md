# V14 final decision

The preregistered multi-PWM high-fidelity gate **failed**.  The one fitted
fan-flow exponent was frozen after PWM40 calibration and was not adjusted after
opening PWM60 or PWM80.

| Evidence block | Result |
|---|---:|
| PWM40 calibration RMSE | 1.562 °C |
| PWM60 independent validation RMSE | 2.607 °C |
| PWM80 final-confirmation RMSE | 3.055 °C |
| PWM80 maximum absolute error | 6.111 °C |
| Independent spatial RMSE retained from V12 | 1.197 °C |
| Fixed-60% regional dynamic RMSE retained from V13 | 1.485 °C |
| Medium–fine maximum sensor difference | 0.099 °C |
| Energy residual | 7.63e-14% |

The PWM60 and PWM80 RMSE thresholds, and the PWM80 maximum-error threshold,
failed.  This is evidence that a single static PWM-to-flow power law does not
represent the experimental fan/stack interaction across duty ratios.  Adding
parameters after observing the confirmation case would invalidate the frozen
test, so no further fan-curve fitting is performed in V14.

The result does not erase the useful spatial and fixed-60% PWM evidence.  The
next admissible step is an independent COMSOL implementation check of that
narrower thermal-plant scope.  Variable-duty ROM, observer, and MPC conclusions
remain unauthorized.  The 55 °C value remains a reference temperature only.
