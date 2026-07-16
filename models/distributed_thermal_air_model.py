"""High-resolution one-dimensional coupled solid-air thermal model.

This model is independent of the five-node controller discretization.  It uses
paper material/heat-flux information and Fig. 2's digitized fan-flow curve, with
four explicitly identified equivalent parameters: total solid-air UA, external
loss conductance, axial solid conductance and effective heat-flux fraction.
"""

from __future__ import annotations

from dataclasses import dataclass
import numpy as np


@dataclass(frozen=True)
class DistributedParameters:
    node_count: int
    total_UA_W_K: float
    total_external_loss_W_K: float
    axial_edge_conductance_W_K: float
    effective_heat_fraction: float
    total_heat_capacity_J_K: float = 10.5 * 460.0
    heat_flux_W_m2: float = 2425.5
    active_area_m2: float = 0.0247
    air_cp_J_kgK: float = 1007.0
    inlet_air_temperature_C: float = 25.0
    ambient_temperature_C: float = 25.0

    @property
    def total_heat_W(self) -> float:
        return self.heat_flux_W_m2 * self.active_area_m2 * self.effective_heat_fraction


@dataclass(frozen=True)
class DistributedSteadyResult:
    solid_temperature_C: np.ndarray
    air_node_inlet_C: np.ndarray
    air_node_outlet_C: np.ndarray
    air_heat_W: np.ndarray
    external_loss_W: np.ndarray
    energy_residual_W: float


def _air_effectiveness(parameters: DistributedParameters, mass_flow_kg_s: float) -> float:
    if mass_flow_kg_s <= 0:
        raise ValueError("mass_flow_kg_s must be positive")
    ua_node = parameters.total_UA_W_K / parameters.node_count
    return float(1.0 - np.exp(-ua_node / (mass_flow_kg_s * parameters.air_cp_J_kgK)))


def steady_state(
    parameters: DistributedParameters,
    mass_flow_kg_s: float,
    direction: int = 1,
) -> DistributedSteadyResult:
    """Solve the linear steady coupled problem exactly.

    Direction +1 means air traverses node 0 -> N-1; -1 reverses the traversal.
    """
    if direction not in (-1, 1):
        raise ValueError("direction must be -1 or 1")
    n = parameters.node_count
    q_node = parameters.total_heat_W / n
    h_loss_node = parameters.total_external_loss_W_K / n
    k_edge = parameters.axial_edge_conductance_W_K
    effectiveness = _air_effectiveness(parameters, mass_flow_kg_s)
    traversal = np.arange(n) if direction == 1 else np.arange(n - 1, -1, -1)

    def residual(temperature_C: np.ndarray) -> np.ndarray:
        air = parameters.inlet_air_temperature_C
        r = np.zeros(n)
        for node in traversal:
            q_air = mass_flow_kg_s * parameters.air_cp_J_kgK * effectiveness * (
                temperature_C[node] - air
            )
            conduction = 0.0
            if node > 0:
                conduction += k_edge * (temperature_C[node - 1] - temperature_C[node])
            if node < n - 1:
                conduction += k_edge * (temperature_C[node + 1] - temperature_C[node])
            r[node] = (
                q_node
                - q_air
                - h_loss_node * (temperature_C[node] - parameters.ambient_temperature_C)
                + conduction
            )
            air += q_air / (mass_flow_kg_s * parameters.air_cp_J_kgK)
        return r

    zero = np.zeros(n)
    offset = residual(zero)
    matrix = np.column_stack([residual(np.eye(n)[i]) - offset for i in range(n)])
    solid = np.linalg.solve(matrix, -offset)
    air_in = np.empty(n)
    air_out = np.empty(n)
    q_air_vector = np.zeros(n)
    air = parameters.inlet_air_temperature_C
    for node in traversal:
        air_in[node] = air
        q_air = mass_flow_kg_s * parameters.air_cp_J_kgK * effectiveness * (solid[node] - air)
        air += q_air / (mass_flow_kg_s * parameters.air_cp_J_kgK)
        air_out[node] = air
        q_air_vector[node] = q_air
    external = h_loss_node * (solid - parameters.ambient_temperature_C)
    residual_W = float(parameters.total_heat_W - np.sum(q_air_vector) - np.sum(external))
    return DistributedSteadyResult(solid, air_in, air_out, q_air_vector, external, residual_W)


def dynamic_reversal(
    parameters: DistributedParameters,
    mass_flow_kg_s: float,
    initial_temperature_C: np.ndarray,
    switch_time_s: float,
    final_time_s: float,
    dt_s: float,
    initial_direction: int = 1,
) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
    """Transient solid dynamics with quasi-steady sequential air heating."""
    n = parameters.node_count
    state = np.asarray(initial_temperature_C, dtype=float).copy()
    if state.shape != (n,):
        raise ValueError("initial_temperature_C has the wrong shape")
    times = np.arange(0.0, final_time_s + 0.5 * dt_s, dt_s)
    history = np.empty((len(times), n))
    outlet = np.empty(len(times))
    residual = np.empty(len(times))
    capacity = parameters.total_heat_capacity_J_K / n
    q_node = parameters.total_heat_W / n
    h_loss_node = parameters.total_external_loss_W_K / n
    k_edge = parameters.axial_edge_conductance_W_K
    effectiveness = _air_effectiveness(parameters, mass_flow_kg_s)
    for index, time_s in enumerate(times):
        history[index] = state
        direction = initial_direction if time_s < switch_time_s else -initial_direction
        traversal = np.arange(n) if direction == 1 else np.arange(n - 1, -1, -1)
        air = parameters.inlet_air_temperature_C
        q_air_vector = np.zeros(n)
        for node in traversal:
            q_air_vector[node] = mass_flow_kg_s * parameters.air_cp_J_kgK * effectiveness * (state[node] - air)
            air += q_air_vector[node] / (mass_flow_kg_s * parameters.air_cp_J_kgK)
        outlet[index] = air
        conduction = np.zeros(n)
        edge = k_edge * np.diff(state)
        conduction[:-1] += edge
        conduction[1:] -= edge
        external = h_loss_node * (state - parameters.ambient_temperature_C)
        rhs = q_node - q_air_vector - external + conduction
        residual[index] = float(parameters.total_heat_W - np.sum(q_air_vector) - np.sum(external) - capacity * np.sum(rhs / capacity))
        if index < len(times) - 1:
            state = state + dt_s * rhs / capacity
    return times, history, outlet, residual
