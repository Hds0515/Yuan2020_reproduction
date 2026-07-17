# V14 multi-PWM preregistration

V14 does not retune the V12 spatial parameters or the V13 fixed-60% thermal
parameters.  It identifies one fan-flow exponent from the independent PWM40
trajectory in Guo et al. Fig. 11, checks cross-experiment PWM60, and opens the
PWM80 trace once as final confirmation.

PWM100 is excluded before fitting because Guo et al. explicitly identify
flooding in that experiment.  V14 contains no membrane-water or flooding state,
so treating PWM100 as a thermal-model validation case would violate model
scope.  PWM40 data after the paper's declared overheating termination are also
excluded.

Passing V14 combines three independent layers of evidence: V12 spatial field,
V13 regional dynamics, and V14 variable-fan average temperature.  It does not
make 55 °C a safety limit.
