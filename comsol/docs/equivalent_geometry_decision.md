# Equivalent geometry decision

## Selected option: B — representative parallel channel

The active COMSOL model solves one declared representative channel and scales every reported extensive flow quantity by

`areaScale = Ain / (W * Hair)`.

The declared values are `Ain = 3962 mm²`, `W = 82.333333 mm`, and `Hair = 1 mm`, giving `areaScale ≈ 48.12`. The Stage 1 report evaluates both the scaled surface integral and the independent target

`m_dot,total = rho * u * Ain`.

Their agreement is an explicit acceptance condition; `areaScale` is therefore active rather than a dead parameter. The same multiplier must be applied to heat input, air enthalpy rise, and any other extensive balance in Stage 3.

This option preserves the representative channel's declared hydraulic diameter instead of replacing it with the much larger hydraulic diameter implied by a single full-area duct. The channel count and exact CAD remain unreported by Yuan et al.; consequently this is an equivalent parallel-channel model, not a reconstruction of the unpublished manifold or channel network.

`h_ext` and `Tamb` were removed from Stage 1 because heat transfer is disabled. Stage 3 must declare any external-convection assumption explicitly before those terms are introduced.
