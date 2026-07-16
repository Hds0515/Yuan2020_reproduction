"""Energy-consistent five-node extension of the frozen three-node model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from models.three_node_model import ThermalParameters, heat_generation_W


@dataclass(frozen=True)
class FiveNodeParameters:
    total_heat_capacity_J_K: float
    cooling_coefficient_per_node_W_K: float
    conduction_edge_W_K: float
    controller_parameters: ThermalParameters
    coolant_forward_C: np.ndarray

    @classmethod
    def from_three_node(
        cls,
        parameters: ThermalParameters,
        coolant_three_C: np.ndarray,
        node_count: int = 5,
    ) -> "FiveNodeParameters":
        if node_count != 5:
            raise ValueError("This research extension is defined for five nodes")
        coolant = np.interp(np.linspace(0, 1, node_count), np.linspace(0, 1, 3), coolant_three_C)
        # Preserve total cooled area: 3*Kcool in the original becomes 5*Kcool5.
        cooling = parameters.K_cool_W_K * 3.0 / node_count
        # Conductance kA/dx scales inversely with grid spacing.
        conduction = parameters.K_node_W_K * (node_count - 1) / (3 - 1)
        return cls(parameters.C_th_J_K, cooling, conduction, parameters, coolant)


@dataclass(frozen=True)
class FiveNodeStep:
    next_temperature_C: np.ndarray
    generation_W: np.ndarray
    cooling_W: np.ndarray
    edge_conduction_W: np.ndarray
    energy_residual_W: float


def interpolate_three_to_five(temperature_three_C: np.ndarray) -> np.ndarray:
    return np.interp(np.linspace(0, 1, 5), np.linspace(0, 1, 3), np.asarray(temperature_three_C))


def step_five_node(
    state_C: np.ndarray,
    current_A: float,
    duty: float,
    direction: int,
    dt_s: float,
    parameters: FiveNodeParameters,
    heat_model: str,
) -> FiveNodeStep:
    state = np.asarray(state_C, dtype=float)
    if state.shape != (5,):
        raise ValueError("state_C must contain five node temperatures")
    coolant = parameters.coolant_forward_C if direction == 1 else parameters.coolant_forward_C[::-1]
    heat = heat_generation_W(current_A, float(np.mean(state)), parameters.controller_parameters, heat_model)
    generation = np.full(5, heat / 5.0)
    cooling = parameters.cooling_coefficient_per_node_W_K * float(np.clip(duty, 0, 1)) * np.maximum(state - coolant, 0.0)
    # Positive edge heat is defined from node i+1 toward node i.  This sign
    # convention makes conduction dissipative: a hotter neighbour heats the
    # cooler node and loses the same amount of energy.
    edge = parameters.conduction_edge_W_K * np.diff(state)
    conduction = np.zeros(5)
    conduction[:-1] += edge
    conduction[1:] -= edge
    rhs = generation - cooling + conduction
    capacity = parameters.total_heat_capacity_J_K / 5.0
    next_state = state + dt_s * rhs / capacity
    storage = capacity * np.sum(next_state - state) / dt_s
    residual = float(storage - np.sum(generation - cooling))
    return FiveNodeStep(next_state, generation, cooling, edge, residual)


def simulate_open_loop(
    parameters: FiveNodeParameters,
    current_A: np.ndarray,
    duty: np.ndarray,
    direction: np.ndarray,
    initial_temperature_C: np.ndarray,
    dt_s: float,
    heat_model: str,
) -> tuple[np.ndarray, np.ndarray]:
    count = len(current_A)
    temperature = np.empty((count, 5))
    temperature[0] = initial_temperature_C
    residual = np.zeros(count)
    for k in range(count - 1):
        result = step_five_node(temperature[k], current_A[k], duty[k], int(direction[k]), dt_s, parameters, heat_model)
        temperature[k + 1] = result.next_temperature_C
        residual[k] = result.energy_residual_W
    residual[-1] = residual[-2]
    return temperature, residual
