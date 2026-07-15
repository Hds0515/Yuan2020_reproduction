"""Auditable three-node thermal model for the Yuan-2020 reproduction.

The electrochemical voltage is a documented proxy because the paper does not
publish every parameter needed to reconstruct its full electrochemical model.
"""

from __future__ import annotations

from dataclasses import dataclass
import json
import math
from pathlib import Path

import numpy as np


@dataclass(frozen=True)
class ThermalParameters:
    C_th_J_K: float
    K_cool_W_K: float
    K_node_W_K: float
    Kp_1_K: float
    Ki_1_Ks: float
    heat_a0_W: float = 0.0
    heat_a1_W_A: float = 0.0
    heat_a2_W_A2: float = 0.0


@dataclass(frozen=True)
class SimulationResult:
    temperature_C: np.ndarray
    fan_duty: np.ndarray
    airflow_direction: np.ndarray
    node_heat_generation_W: np.ndarray
    node_cooling_heat_W: np.ndarray
    edge_conduction_W: np.ndarray
    node_conduction_net_W: np.ndarray
    node_energy_residual_W: np.ndarray
    total_energy_residual_W: np.ndarray

    @property
    def T1_C(self) -> np.ndarray:
        return self.temperature_C[:, 0]

    @property
    def T2_C(self) -> np.ndarray:
        return self.temperature_C[:, 1]

    @property
    def T3_C(self) -> np.ndarray:
        return self.temperature_C[:, 2]


def load_frozen_parameters(path: str | Path) -> ThermalParameters:
    values = json.loads(Path(path).read_text(encoding="utf-8"))["values"]
    return ThermalParameters(
        C_th_J_K=float(values["C_th_J_K"]),
        K_cool_W_K=float(values["K_cool_W_K"]),
        K_node_W_K=float(values["K_node_W_K"]),
        Kp_1_K=float(values["Kp_1_K"]),
        Ki_1_Ks=float(values["Ki_1_Ks"]),
        heat_a0_W=float(values.get("heat_a0_W", 0.0)),
        heat_a1_W_A=float(values.get("heat_a1_W_A", 0.0)),
        heat_a2_W_A2=float(values.get("heat_a2_W_A2", 0.0)),
    )


def stack_voltage_proxy_V(current_A: float, temperature_C: float) -> float:
    """Reproducible polarization proxy, not the authors' full model."""
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


def heat_generation_W(
    current_A: float,
    mean_temperature_C: float,
    parameters: ThermalParameters,
    heat_model: str,
) -> float:
    current = max(float(current_A), 0.0)
    if heat_model == "A_thermoneutral_proxy":
        # 40 cells * I * thermoneutral voltage - electrical stack power.
        voltage = stack_voltage_proxy_V(current, mean_temperature_C)
        return max(40.0 * current * 1.253 - voltage * current, 0.0)
    if heat_model == "B_quadratic_current":
        heat = parameters.heat_a0_W + parameters.heat_a1_W_A * current + parameters.heat_a2_W_A2 * current**2
        return max(float(heat), 0.0)
    raise ValueError(f"Unknown heat model: {heat_model}")


def simulate(
    parameters: ThermalParameters,
    current_A: np.ndarray,
    initial_temperature_C: np.ndarray,
    *,
    dt_s: float,
    bidirectional: bool,
    initial_airflow_direction: int,
    fan_initially_enabled: bool,
    fan_enable_temperature_C: float,
    reference_temperature_C: float,
    smc_deadband_C: float,
    smc_update_period_s: float,
    heat_model: str,
    coolant_node_temperatures_forward_C: np.ndarray,
) -> SimulationResult:
    """Simulate with explicit Euler integration and conditional anti-windup."""
    if dt_s <= 0:
        raise ValueError("dt_s must be positive")
    if initial_airflow_direction not in (-1, 1):
        raise ValueError("initial_airflow_direction must be -1 or 1")
    current = np.asarray(current_A, dtype=float)
    initial = np.asarray(initial_temperature_C, dtype=float)
    coolant_forward = np.asarray(coolant_node_temperatures_forward_C, dtype=float)
    if initial.shape != (3,) or coolant_forward.shape != (3,):
        raise ValueError("initial temperature and coolant arrays must have three values")
    if current.ndim != 1 or len(current) < 2:
        raise ValueError("current_A must be a one-dimensional array with at least two samples")

    n = len(current)
    temperature = np.empty((n, 3), dtype=float)
    temperature[0] = initial
    duty = np.zeros(n, dtype=float)
    direction = np.full(n, initial_airflow_direction, dtype=int)
    heat_nodes = np.zeros((n, 3), dtype=float)
    cooling_nodes = np.zeros((n, 3), dtype=float)
    edge_conduction = np.zeros((n, 2), dtype=float)
    conduction_net = np.zeros((n, 3), dtype=float)
    node_residual = np.zeros((n, 3), dtype=float)
    total_residual = np.zeros(n, dtype=float)

    integral_error_Ks = 0.0
    enabled = bool(fan_initially_enabled)
    flow_direction = int(initial_airflow_direction)
    next_smc_update_s = 0.0
    node_capacity_J_K = parameters.C_th_J_K / 3.0

    for k in range(n - 1):
        state = temperature[k]
        time_s = k * dt_s
        if not enabled and state[1] >= fan_enable_temperature_C:
            enabled = True

        error = state[1] - reference_temperature_C
        if enabled:
            unsaturated = parameters.Kp_1_K * error + parameters.Ki_1_Ks * integral_error_Ks
            command = float(np.clip(unsaturated, 0.0, 1.0))
            drives_out_of_upper = command >= 1.0 and error < 0.0
            drives_out_of_lower = command <= 0.0 and error > 0.0
            if 0.0 < command < 1.0 or drives_out_of_upper or drives_out_of_lower:
                integral_error_Ks += error * dt_s
            command = float(
                np.clip(
                    parameters.Kp_1_K * error + parameters.Ki_1_Ks * integral_error_Ks,
                    0.0,
                    1.0,
                )
            )
        else:
            command = 0.0

        if bidirectional and time_s + 1e-12 >= next_smc_update_s:
            difference = state[0] - state[2]
            if difference >= smc_deadband_C:
                flow_direction = 1
            elif difference <= -smc_deadband_C:
                flow_direction = -1
            next_smc_update_s += smc_update_period_s
        elif not bidirectional:
            flow_direction = initial_airflow_direction

        coolant = coolant_forward if flow_direction == 1 else coolant_forward[::-1]
        heat_total = heat_generation_W(current[k], float(np.mean(state)), parameters, heat_model)
        generation = np.full(3, heat_total / 3.0)
        cooling = parameters.K_cool_W_K * command * np.maximum(state - coolant, 0.0)
        q12 = parameters.K_node_W_K * (state[0] - state[1])
        q23 = parameters.K_node_W_K * (state[1] - state[2])
        conduction = np.asarray([-q12, q12 - q23, q23])
        rhs_W = generation - cooling + conduction
        delta = dt_s * rhs_W / node_capacity_J_K
        temperature[k + 1] = state + delta

        duty[k] = command
        direction[k] = flow_direction
        heat_nodes[k] = generation
        cooling_nodes[k] = cooling
        edge_conduction[k] = [q12, q23]
        conduction_net[k] = conduction
        storage_W = node_capacity_J_K * (temperature[k + 1] - state) / dt_s
        node_residual[k] = storage_W - rhs_W
        total_residual[k] = float(np.sum(node_residual[k]))

    duty[-1] = duty[-2]
    direction[-1] = direction[-2]
    heat_nodes[-1] = heat_nodes[-2]
    cooling_nodes[-1] = cooling_nodes[-2]
    edge_conduction[-1] = edge_conduction[-2]
    conduction_net[-1] = conduction_net[-2]
    node_residual[-1] = node_residual[-2]
    total_residual[-1] = total_residual[-2]
    return SimulationResult(
        temperature_C=temperature,
        fan_duty=duty,
        airflow_direction=direction,
        node_heat_generation_W=heat_nodes,
        node_cooling_heat_W=cooling_nodes,
        edge_conduction_W=edge_conduction,
        node_conduction_net_W=conduction_net,
        node_energy_residual_W=node_residual,
        total_energy_residual_W=total_residual,
    )
