"""Frozen Yuan-2020 three-node thermal reproduction."""

from __future__ import annotations
import json
import math
from pathlib import Path
import numpy as np


def load_frozen_parameters(path: str | Path) -> np.ndarray:
    data = json.loads(Path(path).read_text(encoding="utf-8"))
    values = data["values"]
    return np.asarray([
        values["C_th_J_K"],
        values["K_cool_W_K"],
        values["K_node_W_K"],
        values["Kp_1_K"],
        values["Ki_1_Ks"],
    ], dtype=float)


def stack_voltage(current_A: float, temperature_C: float) -> float:
    current = max(float(current_A), 0.0)
    limiting_current = 45.0
    effective_current = min(current, limiting_current - 1e-8)
    voltage = (
        38.4
        - 0.8 * math.log1p(effective_current)
        - 0.055 * effective_current
        - 0.0005 * effective_current**2
        - math.log(1.0 / (1.0 - effective_current / limiting_current))
        - 0.002 * (temperature_C - 55.0) ** 2
    )
    return max(voltage, 1.0)


def simulate(
    parameters: np.ndarray,
    current_A: np.ndarray,
    initial_temperature_C: np.ndarray,
    *,
    bidirectional: bool,
    fan_initially_enabled: bool,
    initial_direction: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    C_th, K_cool, K_node, Kp, Ki = parameters
    count = len(current_A)
    temperature = np.empty((count, 3))
    temperature[0] = initial_temperature_C
    duty = np.zeros(count)
    direction = np.ones(count, dtype=int)

    integral_error = 0.0
    enabled = fan_initially_enabled
    flow_direction = initial_direction
    coolant_forward = np.asarray([24.8, 30.05, 35.3])

    for index in range(count - 1):
        state = temperature[index]
        current = current_A[index]

        if not enabled and state[1] >= 53.0:
            enabled = True

        error = state[1] - 55.0
        if enabled:
            command = np.clip(Kp * error + Ki * integral_error, 0.0, 1.0)
            can_integrate = (
                0.0 < command < 1.0
                or (command <= 0.0 and error > 0.0)
                or (command >= 1.0 and error < 0.0)
            )
            if can_integrate:
                integral_error += error
            command = np.clip(Kp * error + Ki * integral_error, 0.0, 1.0)
        else:
            command = 0.0

        if bidirectional and index % 10 == 0:
            difference = state[0] - state[2]
            if difference >= 0.5:
                flow_direction = 1
            elif difference <= -0.5:
                flow_direction = -1
        elif not bidirectional:
            flow_direction = 1

        coolant = coolant_forward if flow_direction == 1 else coolant_forward[::-1]
        voltage = stack_voltage(current, float(np.mean(state)))
        heat = max(40.0 * current * 1.253 - voltage * current, 0.0)
        source = heat / 3.0
        cooling = K_cool * command * np.maximum(state - coolant, 0.0)

        q12 = K_node * (state[0] - state[1])
        q23 = K_node * (state[1] - state[2])
        conduction = np.asarray([-q12, q12 - q23, q23])

        temperature[index + 1] = state + (
            source - cooling + conduction
        ) / (C_th / 3.0)
        duty[index] = command
        direction[index] = flow_direction

    duty[-1] = duty[-2]
    direction[-1] = direction[-2]
    return temperature, duty, direction
