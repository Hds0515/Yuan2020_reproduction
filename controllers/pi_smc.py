"""Paper-inspired PI fan command with discrete airflow-direction switching."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np


@dataclass
class PISMCController:
    kp_1_K: float
    ki_1_Ks: float
    reference_C: float = 55.0
    deadband_C: float = 0.5
    update_period_s: float = 10.0
    direction: int = -1
    feedback_mode: str = "middle"
    maximum_duty_step: float = 0.1
    maximum_duty: float = 1.0
    previous_duty: float = 0.0
    integral_error_Ks: float = 0.0
    next_direction_update_s: float = 0.0

    def command(self, state_estimate_C: np.ndarray, time_s: float, dt_s: float) -> tuple[float, int]:
        state = np.asarray(state_estimate_C, dtype=float)
        if self.feedback_mode == "middle":
            feedback = float(state[len(state) // 2])
        elif self.feedback_mode == "mean":
            feedback = float(np.mean(state))
        elif self.feedback_mode == "hotspot":
            feedback = float(np.max(state))
        else:
            raise ValueError(f"Unknown feedback_mode: {self.feedback_mode}")
        error = feedback - self.reference_C
        unsaturated = self.kp_1_K * error + self.ki_1_Ks * self.integral_error_Ks
        duty = float(np.clip(unsaturated, 0.0, self.maximum_duty))
        if 0.0 < duty < self.maximum_duty or (duty >= self.maximum_duty and error < 0.0) or (duty <= 0.0 and error > 0.0):
            self.integral_error_Ks += error * dt_s
        duty = float(np.clip(self.kp_1_K * error + self.ki_1_Ks * self.integral_error_Ks, 0.0, self.maximum_duty))
        duty = float(np.clip(duty, self.previous_duty - self.maximum_duty_step,
                             self.previous_duty + self.maximum_duty_step))
        if time_s + 1e-12 >= self.next_direction_update_s:
            end_difference = float(state[0] - state[-1])
            if end_difference >= self.deadband_C:
                self.direction = 1
            elif end_difference <= -self.deadband_C:
                self.direction = -1
            self.next_direction_update_s += self.update_period_s
        self.previous_duty = duty
        return duty, self.direction
