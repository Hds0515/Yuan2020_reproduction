"""Audited entry point for the V15 fixed-60%-PWM reference plant.

The numerical equations and frozen parameters are inherited from V13.  V15
adds no fit parameters; it packages the 4800-state mesh with an explicit scope
guard and links it to the independent COMSOL implementation check.
"""

from __future__ import annotations

from dataclasses import asdict

from models_v13.fixed_pwm_plant import DynamicParameters, FixedPWM60Plant


FROZEN_PARENT_SHA256 = "223c4a77df0c44dfa401e014348db87d0eeb74a358da8f4b516b2835a5dffb1d"
FROZEN_DYNAMIC = DynamicParameters(
    effective_total_heat_capacity_J_K=3999.999989994387,
    surface_heat_capacity_fraction=0.010010000000056184,
    surface_core_exchange_W_K_total=8.065428145396455,
    fixed60_airflow_multiplier_relative_to_v12=2.072349298596786,
    current_dependent_axial_heat_skew=0.4128642328339601,
)


class Fixed60HighResolutionReference(FixedPWM60Plant):
    """4800-state thermal reference with a deliberately narrow validity domain."""

    ambient_C = 23.0
    fixed_duty = 0.60
    minimum_current_A = 2.0
    maximum_current_A = 40.0

    def __init__(self) -> None:
        super().__init__(FROZEN_DYNAMIC, flow_rows=6, columns=10)

    def assert_in_scope(self, *, ambient_C: float, duty: float, current_A: float) -> None:
        if abs(float(ambient_C) - self.ambient_C) > 1e-9:
            raise ValueError("V15 is validated only at 23 degC ambient")
        if abs(float(duty) - self.fixed_duty) > 1e-9:
            raise ValueError("V15 is not a variable-PWM plant; duty must remain 0.60")
        if not self.minimum_current_A <= float(current_A) <= self.maximum_current_A:
            raise ValueError("current is outside the 2-40 A evidence range")

    def validation_record(self) -> dict[str, object]:
        return {
            "model_class": type(self).__name__,
            "frozen_parent_sha256": FROZEN_PARENT_SHA256,
            "dynamic_parameters": asdict(FROZEN_DYNAMIC),
            "state_count": self.state_count,
            "scope": {
                "ambient_C": self.ambient_C,
                "duty": self.fixed_duty,
                "current_A": [self.minimum_current_A, self.maximum_current_A],
                "measured_current_and_voltage_are_exogenous": True,
            },
            "limitations": [
                "not validated for variable PWM",
                "does not resolve membrane water or flooding",
                "does not reconstruct the original fan plenum or manifold CAD",
                "V13 strict maximum-error and ordering thresholds were narrowly missed",
            ],
        }
