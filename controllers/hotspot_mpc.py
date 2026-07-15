"""Deterministic control-blocking MPC for hotspot and gradient regulation."""

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
    mean_temperature_weight: float = 0.0
    gradient_weight: float = 2.0
    fan_weight: float = 0.22
    duty_move_weight: float = 0.8
    direction_switch_weight: float = 2.0
    minimum_direction_dwell_s: float = 10.0
    reversal_dead_time_s: float = 2.0
    control_dt_s: float = 1.0
    maximum_duty_step: float = 0.1
    maximum_duty: float = 1.0
    previous_duty: float = 0.0
    previous_direction: int = -1
    time_since_direction_change_s: float = 10.0
    reversal_dead_time_remaining_s: float = 0.0
    last_planned_duty: tuple[float, ...] = ()
    last_planned_direction: tuple[int, ...] = ()

    @property
    def steps_since_direction_change(self) -> int:
        """Backward-compatible view used by the original tests and callers."""
        return int(round(self.time_since_direction_change_s / self.control_dt_s))

    @steps_since_direction_change.setter
    def steps_since_direction_change(self, value: int) -> None:
        self.time_since_direction_change_s = float(value) * self.control_dt_s

    def _block_lengths(self) -> tuple[int, int, int]:
        if self.horizon_steps < 3:
            raise ValueError("control-blocking MPC requires a horizon of at least three steps")
        quotient, remainder = divmod(self.horizon_steps, 3)
        return tuple(quotient + (1 if i < remainder else 0) for i in range(3))  # type: ignore[return-value]

    def _duty_sequences(self) -> list[tuple[float, float, float]]:
        step = self.maximum_duty_step
        clip = lambda value: float(np.clip(value, 0.0, self.maximum_duty))
        previous = self.previous_duty
        candidates = {(previous, previous, previous)}
        for sign in (-1.0, 1.0):
            first = clip(previous + sign * step)
            second = clip(first + sign * step)
            third = clip(second + sign * step)
            candidates.add((first, first, first))
            candidates.add((first, second, third))
            candidates.add((previous, first, second))
            candidates.add((first, previous, previous))
        return sorted(candidates)

    def _direction_sequences(self, blocks: tuple[int, int, int]) -> list[tuple[int, int, int]]:
        valid: list[tuple[int, int, int]] = []
        other = -self.previous_direction
        candidates = (
            (self.previous_direction, self.previous_direction, self.previous_direction),
            (other, other, other),
            (self.previous_direction, other, other),
            (self.previous_direction, self.previous_direction, other),
        )
        for sequence in candidates:
            previous = self.previous_direction
            elapsed = self.time_since_direction_change_s
            feasible = True
            for direction, length in zip(sequence, blocks, strict=True):
                if direction != previous:
                    if elapsed + 1e-12 < self.minimum_direction_dwell_s:
                        feasible = False
                        break
                    elapsed = 0.0
                    previous = direction
                elapsed += length * self.prediction_dt_s
            if feasible:
                valid.append(sequence)
        return valid

    def command(self, state_estimate_C: np.ndarray, current_preview_A: np.ndarray) -> tuple[float, int, float]:
        preview = np.asarray(current_preview_A, dtype=float)
        if preview.size == 0:
            raise ValueError("current_preview_A cannot be empty")
        blocks = self._block_lengths()
        duty_sequences = self._duty_sequences()
        direction_sequences = self._direction_sequences(blocks)
        best: tuple[float, tuple[float, float, float], tuple[int, int, int]] | None = None
        for directions in direction_sequences:
            for duties in duty_sequences:
                state = np.asarray(state_estimate_C, dtype=float).copy()
                cost = 0.0
                previous_duty = self.previous_duty
                previous_direction = self.previous_direction
                preview_index = 0
                for duty, direction, length in zip(duties, directions, blocks, strict=True):
                    cost += self.duty_move_weight * (duty - previous_duty) ** 2
                    changed = direction != previous_direction
                    cost += self.direction_switch_weight * float(changed)
                    dead_steps = int(np.ceil(self.reversal_dead_time_s / self.prediction_dt_s)) if changed else 0
                    for local_step in range(length):
                        current = float(preview[min(preview_index, len(preview) - 1)])
                        effective_duty = 0.0 if local_step < dead_steps else duty
                        state = step_five_node(
                            state, current, effective_duty, direction, self.prediction_dt_s,
                            self.parameters, self.heat_model,
                        ).next_temperature_C
                        hotspot_violation = max(float(np.max(state) - self.reference_C), 0.0)
                        mean_error = float(np.mean(state) - self.reference_C)
                        gradient = float(np.max(state) - np.min(state))
                        cost += self.hotspot_weight * hotspot_violation**2
                        cost += self.mean_temperature_weight * mean_error**2
                        cost += self.gradient_weight * gradient**2
                        cost += self.fan_weight * effective_duty**2
                        preview_index += 1
                    previous_duty = duty
                    previous_direction = direction
                proposal = (float(cost), duties, directions)
                if best is None or proposal[0] < best[0]:
                    best = proposal
        assert best is not None
        cost, duties, directions = best
        duty, direction = duties[0], directions[0]
        self.last_planned_duty = duties
        self.last_planned_direction = directions
        if direction != self.previous_direction:
            self.time_since_direction_change_s = 0.0
            self.reversal_dead_time_remaining_s = self.reversal_dead_time_s
        else:
            self.time_since_direction_change_s += self.control_dt_s
        if self.reversal_dead_time_remaining_s > 0.0:
            duty = 0.0
            self.reversal_dead_time_remaining_s = max(
                0.0, self.reversal_dead_time_remaining_s - self.control_dt_s
            )
        self.previous_duty = float(duty)
        self.previous_direction = direction
        return float(duty), int(direction), float(cost)
