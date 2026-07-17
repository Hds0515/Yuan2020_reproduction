"""Variable-PWM extension of the frozen V13 thermal plant."""

from __future__ import annotations

import numpy as np

from models_v13.fixed_pwm_plant import DynamicParameters, FixedPWM60Plant


V13_DYNAMIC = DynamicParameters(
    effective_total_heat_capacity_J_K=3999.999989994387,
    surface_heat_capacity_fraction=0.010010000000056184,
    surface_core_exchange_W_K_total=8.065428145396455,
    fixed60_airflow_multiplier_relative_to_v12=2.072349298596786,
    current_dependent_axial_heat_skew=0.4128642328339601,
)
V13_FROZEN_SHA256 = "223c4a77df0c44dfa401e014348db87d0eeb74a358da8f4b516b2835a5dffb1d"


class MultiPWMPlant(FixedPWM60Plant):
    def __init__(
        self,
        fan_flow_exponent_relative_to_60pct: float,
        *,
        flow_rows: int = 3,
        columns: int = 5,
    ) -> None:
        self.fan_flow_exponent_relative_to_60pct = float(
            fan_flow_exponent_relative_to_60pct
        )
        super().__init__(V13_DYNAMIC, flow_rows=flow_rows, columns=columns)

    def _v13_flow_at_60pct(self) -> float:
        fixed = self.physics
        normalized = (
            (0.60 - fixed.minimum_commanded_duty)
            / (fixed.reference_duty - fixed.minimum_commanded_duty)
        )
        factor = fixed.fan_start_flow_fraction + (
            1.0 - fixed.fan_start_flow_fraction
        ) * normalized**self.parameters.pwm_to_flow_exponent
        return self.parameters.effective_total_airflow_at_57pct_kg_s * factor

    def total_airflow_kg_s(self, duty: float) -> float:
        fixed = self.physics
        ratio = np.clip(
            (float(duty) - fixed.minimum_commanded_duty)
            / (0.60 - fixed.minimum_commanded_duty),
            0.0,
            None,
        )
        factor = fixed.fan_start_flow_fraction + (
            1.0 - fixed.fan_start_flow_fraction
        ) * ratio**self.fan_flow_exponent_relative_to_60pct
        return float(self._v13_flow_at_60pct() * factor)

    def record(self) -> dict[str, object]:
        return {
            **self.v13_record(),
            "v13_frozen_sha256": V13_FROZEN_SHA256,
            "fan_flow_exponent_relative_to_60pct": self.fan_flow_exponent_relative_to_60pct,
            "validated_PWM_range_percent": [40, 80],
        }
