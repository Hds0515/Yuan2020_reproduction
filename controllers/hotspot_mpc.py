"""Small deterministic nonlinear MPC for hotspot and gradient regulation."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from models.five_node_model import FiveNodeParameters, step_five_node


@dataclass
class HotspotMPC:
    parameters: FiveNodeParameters
    heat_model: str
    reference_C: float = 55.0
    horizon_steps: int = 12
    prediction_dt_s: float = 3.0
    hotspot_weight: float = 8.0
    gradient_weight: float = 2.0
    fan_weight: float = 0.22
    duty_move_weight: float = 0.8
    direction_switch_weight: float = 2.0
    minimum_direction_dwell_steps: int = 10
    maximum_duty_step: float = 0.1
    maximum_duty: float = 1.0
    previous_duty: float = 0.0
    previous_direction: int = -1
    steps_since_direction_change: int = 10

    def command(self, state_estimate_C: np.ndarray, current_preview_A: np.ndarray) -> tuple[float, int, float]:
        preview = np.asarray(current_preview_A, dtype=float)
        if preview.size == 0:
            raise ValueError("current_preview_A cannot be empty")
        all_candidates = np.unique(np.r_[np.linspace(0.0, self.maximum_duty, 11), self.previous_duty])
        candidates = all_candidates[np.abs(all_candidates - self.previous_duty) <= self.maximum_duty_step + 1e-12]
        best: tuple[float, float, int] | None = None
        allowed_directions = (-1, 1) if self.steps_since_direction_change >= self.minimum_direction_dwell_steps else (self.previous_direction,)
        for direction in allowed_directions:
            for duty in candidates:
                state = np.asarray(state_estimate_C, dtype=float).copy()
                cost = self.duty_move_weight * (duty - self.previous_duty) ** 2
                cost += self.direction_switch_weight * float(direction != self.previous_direction)
                for j in range(self.horizon_steps):
                    current = float(preview[min(j, len(preview) - 1)])
                    state = step_five_node(
                        state, current, float(duty), direction, self.prediction_dt_s,
                        self.parameters, self.heat_model,
                    ).next_temperature_C
                    hotspot_violation = max(float(np.max(state) - self.reference_C), 0.0)
                    gradient = float(np.max(state) - np.min(state))
                    cost += self.hotspot_weight * hotspot_violation**2
                    cost += self.gradient_weight * gradient**2
                    cost += self.fan_weight * float(duty) ** 2
                proposal = (float(cost), float(duty), int(direction))
                if best is None or proposal < best:
                    best = proposal
        assert best is not None
        cost, duty, direction = best
        if direction != self.previous_direction:
            self.steps_since_direction_change = 0
        else:
            self.steps_since_direction_change += 1
        self.previous_duty = duty
        self.previous_direction = direction
        return duty, direction, cost
