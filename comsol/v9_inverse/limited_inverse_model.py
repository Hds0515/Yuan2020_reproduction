"""Limited equivalent inverse model for the Yuan et al. Fig. 7 audit.

This module deliberately identifies only four *equivalent* parameters.  It is
not a CAD reconstruction and must not be described as the authors' original
geometry.  The extensive quantities (total inlet area, heat input, material
properties and external heat loss law) remain frozen from V8.

Identified parameters
---------------------
h0_W_m2K_at_8m_s
    Effective internal heat-transfer coefficient at 8 m/s.  It lumps the
    unresolved hydraulic diameter and effective flow-solid contact area.
velocity_exponent_m
    Velocity sensitivity in h ~ u**m.
downstream_cooling_bias_beta
    Spatial cooling bias.  Positive beta weakens cooling toward the outlet.
heat_source_skew_eta
    Linear heat-source skew.  Positive eta shifts heat generation toward the
    outlet; the total heat input remains exactly conserved.
"""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np

from comsol.v8_runtime.cooling_structure_models import (
    PHYSICS,
    SIGMA,
    churchill_chu_h_W_m2K,
)


PARAMETER_NAMES = (
    "h0_W_m2K_at_8m_s",
    "velocity_exponent_m",
    "downstream_cooling_bias_beta",
    "heat_source_skew_eta",
)

LOWER_BOUNDS = np.array([1.0, 0.10, -8.0, -0.95], dtype=float)
UPPER_BOUNDS = np.array([300.0, 1.50, 8.0, 0.95], dtype=float)
DEFAULT_INITIAL = np.array([90.0, 0.65, 1.0, -0.25], dtype=float)


@dataclass(frozen=True)
class InverseModelResult:
    velocity_m_s: float
    flow_coordinate: np.ndarray
    surface_C: np.ndarray
    air_C: np.ndarray
    T_inlet_region_C: float
    T_middle_region_C: float
    T_outlet_region_C: float
    Tmax_C: float
    Tmin_C: float
    DeltaT_C: float
    hotspot_location_normalized: float
    air_outlet_temperature_C: float
    air_enthalpy_gain_W: float
    natural_convection_loss_W: float
    radiation_loss_W: float
    total_heat_input_W: float
    energy_residual_relative: float
    runtime_s: float


def _air_heat_matrix(effectiveness: np.ndarray, mcp_W_K: float) -> tuple[np.ndarray, np.ndarray]:
    """Return q_air = C @ T_surface + d for sequential plug-flow heating."""
    count = len(effectiveness)
    matrix = np.zeros((count, count), dtype=float)
    offset = np.zeros(count, dtype=float)
    air_constant = PHYSICS.ambient_C
    air_coefficients = np.zeros(count, dtype=float)
    for index, epsilon in enumerate(effectiveness):
        matrix[index, index] += mcp_W_K * epsilon
        matrix[index, :] -= mcp_W_K * epsilon * air_coefficients
        offset[index] = -mcp_W_K * epsilon * air_constant
        air_constant = (1.0 - epsilon) * air_constant
        air_coefficients = (1.0 - epsilon) * air_coefficients
        air_coefficients[index] += epsilon
    return matrix, offset


def _conduction_matrix(nodes: int, conductance_W_K: float) -> np.ndarray:
    matrix = np.zeros((nodes, nodes), dtype=float)
    for index in range(nodes):
        if index > 0:
            matrix[index, index - 1] += conductance_W_K
            matrix[index, index] -= conductance_W_K
        if index < nodes - 1:
            matrix[index, index + 1] += conductance_W_K
            matrix[index, index] -= conductance_W_K
    return matrix


def solve_inverse_model(
    velocity_m_s: float,
    parameters: np.ndarray,
    *,
    nodes: int = 120,
    max_iterations: int = 100,
    tolerance_C: float = 1.0e-9,
) -> InverseModelResult:
    """Solve the conservative steady 1-D equivalent conjugate model."""
    start = perf_counter()
    parameters = np.asarray(parameters, dtype=float)
    if parameters.shape != (4,):
        raise ValueError("parameters must contain exactly four values")
    if np.any(parameters < LOWER_BOUNDS) or np.any(parameters > UPPER_BOUNDS):
        raise ValueError("parameters are outside the frozen bounds")
    if velocity_m_s <= 0.0:
        raise ValueError("velocity_m_s must be positive")
    if nodes < 12:
        raise ValueError("nodes must be at least 12")

    h0, exponent_m, cooling_bias, heat_skew = map(float, parameters)
    coordinate = (np.arange(nodes, dtype=float) + 0.5) / nodes
    dx = PHYSICS.length_m / nodes

    mdot = PHYSICS.air_rho_kg_m3 * velocity_m_s * PHYSICS.stack_inlet_area_m2
    mcp = mdot * PHYSICS.air_cp_J_kgK

    h_base = h0 * (velocity_m_s / 8.0) ** exponent_m
    # The exponential field is not normalised: beta represents a true
    # downstream cooling imbalance, not a redistribution with fixed mean h.
    h_internal = h_base * np.exp(-cooling_bias * (coordinate - 0.5))
    ua_segment = h_internal * PHYSICS.width_m * dx
    effectiveness = 1.0 - np.exp(-ua_segment / mcp)
    q_air_matrix, q_air_offset = _air_heat_matrix(effectiveness, mcp)

    heat_shape = 1.0 + heat_skew * (2.0 * coordinate - 1.0)
    if np.any(heat_shape <= 0.0):
        raise ValueError("heat-source skew produced a non-positive local heat source")
    q_heat = PHYSICS.heat_input_cell_W * heat_shape / heat_shape.sum()

    axial_conductance = (
        PHYSICS.graphite_k_W_mK * PHYSICS.solid_cross_section_m2 / dx
    )
    conduction_matrix = _conduction_matrix(nodes, axial_conductance)
    external_area_segment = PHYSICS.external_area_cell_m2 / nodes

    surface_C = np.linspace(45.0, 60.0, nodes)
    ambient_C = PHYSICS.ambient_C
    for _ in range(max_iterations):
        h_natural = churchill_chu_h_W_m2K(surface_C)
        surface_K = surface_C + 273.15
        ambient_K = ambient_C + 273.15
        h_radiation = PHYSICS.emissivity * SIGMA * (
            surface_K**3
            + surface_K**2 * ambient_K
            + surface_K * ambient_K**2
            + ambient_K**3
        )
        external_conductance = (h_natural + h_radiation) * external_area_segment
        system_matrix = (
            conduction_matrix
            - q_air_matrix
            - np.diag(external_conductance)
        )
        right_hand_side = (
            q_air_offset
            - q_heat
            - external_conductance * ambient_C
        )
        updated = np.linalg.solve(system_matrix, right_hand_side)
        if np.max(np.abs(updated - surface_C)) < tolerance_C:
            surface_C = updated
            break
        surface_C = 0.5 * surface_C + 0.5 * updated
    else:
        raise RuntimeError("inverse model nonlinear heat-loss iteration did not converge")

    air_C = np.empty(nodes + 1, dtype=float)
    air_C[0] = ambient_C
    q_air = np.empty(nodes, dtype=float)
    for index, epsilon in enumerate(effectiveness):
        air_C[index + 1] = air_C[index] + epsilon * (
            surface_C[index] - air_C[index]
        )
        q_air[index] = mcp * (air_C[index + 1] - air_C[index])

    h_natural = churchill_chu_h_W_m2K(surface_C)
    surface_K = surface_C + 273.15
    ambient_K = ambient_C + 273.15
    q_natural = h_natural * external_area_segment * (surface_C - ambient_C)
    q_radiation = (
        PHYSICS.emissivity
        * SIGMA
        * external_area_segment
        * (surface_K**4 - ambient_K**4)
    )
    q_air_total = float(q_air.sum())
    q_natural_total = float(q_natural.sum())
    q_radiation_total = float(q_radiation.sum())
    energy_residual = abs(
        PHYSICS.heat_input_cell_W
        - q_air_total
        - q_natural_total
        - q_radiation_total
    ) / PHYSICS.heat_input_cell_W

    thirds = np.array_split(surface_C, 3)
    hotspot_index = int(np.argmax(surface_C))
    return InverseModelResult(
        velocity_m_s=float(velocity_m_s),
        flow_coordinate=coordinate,
        surface_C=surface_C,
        air_C=air_C,
        T_inlet_region_C=float(np.mean(thirds[0])),
        T_middle_region_C=float(np.mean(thirds[1])),
        T_outlet_region_C=float(np.mean(thirds[2])),
        Tmax_C=float(np.max(surface_C)),
        Tmin_C=float(np.min(surface_C)),
        DeltaT_C=float(np.ptp(surface_C)),
        hotspot_location_normalized=float(coordinate[hotspot_index]),
        air_outlet_temperature_C=float(air_C[-1]),
        air_enthalpy_gain_W=q_air_total,
        natural_convection_loss_W=q_natural_total,
        radiation_loss_W=q_radiation_total,
        total_heat_input_W=float(PHYSICS.heat_input_cell_W),
        energy_residual_relative=float(energy_residual),
        runtime_s=float(perf_counter() - start),
    )


def result_row(result: InverseModelResult) -> dict[str, float]:
    """Convert a result into a compact machine-readable row."""
    return {
        "velocity_m_s": result.velocity_m_s,
        "T_inlet_region_C": result.T_inlet_region_C,
        "T_middle_region_C": result.T_middle_region_C,
        "T_outlet_region_C": result.T_outlet_region_C,
        "Tmax_C": result.Tmax_C,
        "Tmin_C": result.Tmin_C,
        "DeltaT_C": result.DeltaT_C,
        "hotspot_location_normalized": result.hotspot_location_normalized,
        "air_outlet_temperature_C": result.air_outlet_temperature_C,
        "air_enthalpy_gain_W": result.air_enthalpy_gain_W,
        "natural_convection_loss_W": result.natural_convection_loss_W,
        "radiation_loss_W": result.radiation_loss_W,
        "total_heat_input_W": result.total_heat_input_W,
        "energy_residual_relative": result.energy_residual_relative,
        "runtime_s": result.runtime_s,
    }
