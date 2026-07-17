# V12 model-structure preregistration

V12 targets a high-resolution **thermal digital twin**, not recovery of an
unpublished CAD model and not an electrochemical multiphysics model.  The
reference experiment is the 40-cell, 45-channel open-cathode stack of Gong et
al. (2024), which reports the active area, channel dimensions, fan model, 75
thermocouple locations, a steady 30 A temperature field, and two dynamic load
programs.

The V11 negative result is retained under tag
`v11_high_fidelity_gate_failed`.  It showed that a conductance-per-region
network is mesh dependent and that linear PWM-to-flow mapping cannot reproduce
the dynamic experiment.  V12 therefore changes structure before any new
optimization.

## Physical structure

Each of the 40 cells is divided into three cathode-flow regions and five
in-plane regions.  Every region has a surface/plate temperature and an internal
assembly temperature, giving 1200 states.  Surface nodes exchange heat through
area-and-distance finite-volume conductances.  Internal nodes exchange heat
only with their colocated surface node, which provides a physically
interpretable fast/slow thermal response without adding an unmeasured heat
source.

Cathode air is marched from inlet to outlet with an exact heat-exchanger
effectiveness relation, so its enthalpy gain equals surface heat loss.  The fan
footprint includes an axial hub deficit followed by plenum mixing.  Total flow
above the reported 12% PWM start point follows a bounded power law.  PWM fan
speed curves are model-specific, so the exponent is calibrated and explicitly
bounded rather than asserted as a manufacturer value for the unpublished fan
controller.

Heat generation is

`I * (40 * 1.253 - V_stack)`.

For dynamic experiments, `I` and `V_stack` are digitized inputs.  The model
does not predict voltage.  Natural convection and radiation are external loss
paths; graphite emissivity remains fixed at 0.82.

## Calibration isolation

Static groups 1, 2, 4, and 5 plus Fig. 9(a) are calibration data.  Static group
3 and Fig. 9(b) are final confirmation data and cannot enter the objective.
The predecessor's failure means this is not a fully blinded replication; that
limitation will remain explicit in the final report.  Nevertheless, V12
parameters are frozen before its confirmation calculation and are not adjusted
afterward.

COMSOL 6.4 is an independent implementation gate, not a second fitting
environment.  Its geometry, coefficients, inputs, and frozen parameters must
match the Python finite-volume model within 0.20 °C at the aggregated sensors.

