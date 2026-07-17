# Validation result and research scope

## What is now usable

The repository now contains a reproducible 4800-state stack thermal reference
for fixed-60% PWM experiments.  It is suitable for:

- reproducing the published spatial temperature field at the frozen operating
  point;
- studying heat-load transients within the fixed-duty evidence domain;
- testing aggregation and numerical integration methods without changing the
  plant parameters;
- developing software interfaces that will later accept a broader plant.

## What remains prohibited

It is not suitable for decisive variable-duty ROM, MM-EKF, or Oracle MPC
claims.  The V14 confirmation errors increased from 1.562 °C at PWM40 to
2.607 °C at PWM60 and 3.055 °C at PWM80.  This establishes a missing physical
input model rather than a numerical implementation problem.  Likely omitted
states include fan pressure-flow interaction, duty-dependent plenum flow,
membrane water, and flooding.  Average-temperature curves cannot uniquely
identify those mechanisms.

The V13 fixed-duty evidence is strong but not a clean strict-gate pass: its
maximum error was 4.271 °C against a 4.0 °C threshold, and regional ordering
was 94.737% against 95%.  Those thresholds are retained unchanged.

## Decision

The delivered artifact may be called an experimentally constrained,
high-resolution fixed-duty reference model.  It must not be called a validated
variable-PWM high-fidelity control plant.  A future variable-duty plant needs
new measured fan speed/flow or pressure telemetry and an independent
confirmation experiment reserved before fitting.
