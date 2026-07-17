"""High-resolution thermal reference for the Gong et al. (2024) stack.

The model resolves 40 cells, three cathode-flow segments and five in-plane
regions per cell (600 thermal states).  Air enthalpy is marched conservatively
through each cathode path.  Heat conduction is resolved in all three structured
directions, and the fan footprint, stack-edge losses and serpentine heat-source
pattern are explicit.

It is not an electrochemical model.  Stack voltage is an exogenous/empirical
input used only to close the measured electrical-to-thermal energy balance.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from time import perf_counter

import numpy as np
from scipy.sparse import bmat, csr_matrix, diags, lil_matrix
from scipy.sparse.linalg import expm_multiply, spsolve


@dataclass(frozen=True)
class GongStackGeometry:
    cell_count: int = 40
    active_area_m2: float = 0.0051
    cathode_channel_count: int = 45
    cathode_channel_depth_m: float = 0.0014
    cathode_channel_width_m: float = 0.0010
    cathode_channel_length_m: float = 0.0608
    cathode_ridge_width_m: float = 0.0014
    bipolar_plate_thickness_m: float = 0.0025
    mea_thickness_m: float = 0.00062
    stack_length_m: float = 0.13642
    total_stack_mass_kg: float = 1.871
    fan_diameter_m: float = 0.092
    maximum_fan_flow_CFM: float = 173.46
    maximum_fan_static_pressure_mmAq: float = 65.3

    @property
    def active_width_m(self) -> float:
        return self.active_area_m2 / self.cathode_channel_length_m

    @property
    def repeat_thickness_m(self) -> float:
        return self.stack_length_m / self.cell_count


@dataclass(frozen=True)
class FixedPhysics:
    graphite_in_plane_k_W_mK: float = 60.0
    effective_through_stack_k_W_mK: float = 1.0
    graphite_heat_capacity_J_kgK: float = 710.0
    air_heat_capacity_J_kgK: float = 1006.0
    emissivity: float = 0.82
    stefan_boltzmann_W_m2K4: float = 5.670374419e-8
    reversible_voltage_V_per_cell: float = 1.253
    reference_duty: float = 0.57
    reference_current_A: float = 30.0
    ambient_static_C: float = 25.0


@dataclass(frozen=True)
class CalibratedParameters:
    effective_total_thermal_mass_fraction: float = 0.72
    forced_convection_multiplier: float = 1.0
    effective_total_airflow_at_reference_duty_kg_s: float = 0.050
    stack_edge_heat_transfer_coefficient_W_m2K: float = 8.0
    intercell_contact_conductance_W_K_per_region: float = 0.18
    serpentine_heat_source_nonuniformity: float = 0.18


PARAMETER_NAMES = tuple(CalibratedParameters.__dataclass_fields__)
LOWER_BOUNDS = np.asarray((0.40, 0.45, 0.020, 3.0, 0.04, 0.0), dtype=float)
UPPER_BOUNDS = np.asarray((0.95, 2.20, 0.090, 18.0, 0.80, 0.35), dtype=float)
PRIOR = np.asarray((0.72, 1.0, 0.050, 8.0, 0.18, 0.18), dtype=float)
PRIOR_SCALE = np.asarray((0.18, 0.45, 0.020, 4.0, 0.20, 0.12), dtype=float)


def vector_to_parameters(values: np.ndarray) -> CalibratedParameters:
    return CalibratedParameters(**dict(zip(PARAMETER_NAMES, map(float, values), strict=True)))


def measured_stack_voltage_V(current_A: float) -> float:
    """Smooth voltage closure consistent with Gong Fig. 7/9 operating range."""

    current = max(float(current_A), 0.0)
    return max(20.0, 41.0 - 0.50 * current - 0.002 * current**2)


def fan_duty_from_current(current_A: float) -> float:
    """Published experimental rule: +3 percentage points for each +2 A."""

    current = max(float(current_A), 0.0)
    if current <= 2.0:
        return 0.15
    return float(np.clip(0.15 + 0.015 * (current - 2.0), 0.15, 0.80))


class Gong40CellReference:
    """Structured 3-D stack thermal model with conservative air paths."""

    def __init__(
        self,
        parameters: CalibratedParameters | None = None,
        *,
        cells: int = 40,
        flow_rows: int = 3,
        in_plane_columns: int = 5,
    ) -> None:
        if cells % 40 != 0 or flow_rows % 3 != 0 or in_plane_columns % 5 != 0:
            raise ValueError("mesh dimensions must be integer refinements of 40x3x5")
        self.geometry = GongStackGeometry()
        self.physics = FixedPhysics()
        self.parameters = parameters or CalibratedParameters()
        self.cells = cells
        self.flow_rows = flow_rows
        self.in_plane_columns = in_plane_columns
        self.state_count = cells * flow_rows * in_plane_columns
        self._matrix_cache: dict[tuple[float, float], csr_matrix] = {}

    def index(self, cell: int, row: int, column: int) -> int:
        return (cell * self.flow_rows + row) * self.in_plane_columns + column

    def _fan_weights(self) -> np.ndarray:
        """Return the fixed axial-fan footprint over cells and channel columns.

        The fan hub suppresses flow near the stack centre while the outer ring
        receives more flow.  A smooth U-shaped profile is used because Gong et
        al. report the fan envelope but do not publish a measured face-velocity
        map.  Its mean is one, so it cannot alter the calibrated total airflow.
        """

        z = (np.arange(self.cells) + 0.5) / self.cells
        x = (np.arange(self.in_plane_columns) + 0.5) / self.in_plane_columns
        rz = np.abs(z - 0.5) / 0.5
        rx = np.abs(x - 0.5) / 0.5
        axial_ring = 0.75 + 0.50 * rz**1.8
        lateral_ring = 0.95 + 0.10 * rx**1.5
        footprint = axial_ring[:, None] * lateral_ring[None, :]
        return footprint / np.mean(footprint)

    def _flow_distribution(self, duty: float) -> np.ndarray:
        reference_flow = self.parameters.effective_total_airflow_at_reference_duty_kg_s
        total = reference_flow * max(duty, 0.02) / self.physics.reference_duty
        weights = self._fan_weights()
        return total * weights / np.sum(weights)

    def _conductance_matrix(self, duty: float, ambient_C: float) -> csr_matrix:
        key = (round(float(duty), 8), round(float(ambient_C), 6))
        cached = self._matrix_cache.get(key)
        if cached is not None:
            return cached
        g = self.geometry
        p = self.parameters
        rows = self.flow_rows
        columns = self.in_plane_columns
        dx = g.active_width_m / columns
        dy = g.cathode_channel_length_m / rows
        dz = g.stack_length_m / self.cells
        matrix = lil_matrix((self.state_count, self.state_count), dtype=float)

        # In-plane graphite conduction.
        gx = self.physics.graphite_in_plane_k_W_mK * g.bipolar_plate_thickness_m * dy / dx
        gy = self.physics.graphite_in_plane_k_W_mK * g.bipolar_plate_thickness_m * dx / dy
        gz = p.intercell_contact_conductance_W_K_per_region * 40.0 / self.cells

        def connect(first: int, second: int, conductance: float) -> None:
            matrix[first, first] -= conductance
            matrix[second, second] -= conductance
            matrix[first, second] += conductance
            matrix[second, first] += conductance

        for cell in range(self.cells):
            for row in range(rows):
                for column in range(columns):
                    node = self.index(cell, row, column)
                    if column + 1 < columns:
                        connect(node, self.index(cell, row, column + 1), gx)
                    if row + 1 < rows:
                        connect(node, self.index(cell, row + 1, column), gy)
                    if cell + 1 < self.cells:
                        connect(node, self.index(cell + 1, row, column), gz)

        # Conservative cathode air enthalpy paths.
        flow = self._flow_distribution(duty)
        channels_per_column = g.cathode_channel_count / columns
        wall_area = (
            channels_per_column
            * (g.cathode_channel_width_m + 2.0 * g.cathode_channel_depth_m)
            * dy
        )
        reference_h = 85.0 * p.forced_convection_multiplier
        mean_path_flow = np.mean(flow)
        for cell in range(self.cells):
            for column in range(columns):
                mdot = max(flow[cell, column], 1e-9)
                h_local = reference_h * (mdot / mean_path_flow) ** 0.33
                ua = h_local * wall_area
                mcp = mdot * self.physics.air_heat_capacity_J_kgK
                effectiveness = 1.0 - np.exp(-ua / mcp)
                air_weights: dict[int, float] = {}
                for row in range(rows):
                    node = self.index(cell, row, column)
                    removal = mcp * effectiveness
                    matrix[node, node] -= removal
                    for upstream, weight in air_weights.items():
                        matrix[node, upstream] += removal * weight
                    air_weights = {
                        upstream: (1.0 - effectiveness) * weight
                        for upstream, weight in air_weights.items()
                    }
                    air_weights[node] = air_weights.get(node, 0.0) + effectiveness

        # Stack outer edges: natural convection plus fixed-emissivity radiation.
        ambient_K = ambient_C + 273.15
        reference_surface_K = 55.0 + 273.15
        h_radiation = (
            self.physics.emissivity
            * self.physics.stefan_boltzmann_W_m2K4
            * (reference_surface_K**2 + ambient_K**2)
            * (reference_surface_K + ambient_K)
        )
        h_edge = p.stack_edge_heat_transfer_coefficient_W_m2K + h_radiation
        for cell in range(self.cells):
            for row in range(rows):
                for column in range(columns):
                    area = 0.0
                    if column in (0, columns - 1):
                        area += dz * dy
                    if row in (0, rows - 1):
                        area += dz * dx
                    if cell in (0, self.cells - 1):
                        area += dx * dy
                    matrix[self.index(cell, row, column), self.index(cell, row, column)] -= h_edge * area

        result = matrix.tocsr()
        self._matrix_cache[key] = result
        return result

    def heat_generation_W(self, current_A: float, voltage_V: float | None = None) -> np.ndarray:
        current = max(float(current_A), 0.0)
        voltage = measured_stack_voltage_V(current) if voltage_V is None else float(voltage_V)
        total = max(current * (self.physics.reversible_voltage_V_per_cell * 40.0 - voltage), 0.0)
        column_base = np.asarray((-0.55, 0.30, 0.12, 0.48, -0.35), dtype=float)
        fine_x = (np.arange(self.in_plane_columns) + 0.5) / self.in_plane_columns
        base_x = (np.arange(5) + 0.5) / 5.0
        column_pattern = np.interp(fine_x, base_x, column_base)
        shape = np.ones((self.cells, self.flow_rows, self.in_plane_columns), dtype=float)
        shape *= 1.0 + self.parameters.serpentine_heat_source_nonuniformity * column_pattern[None, None, :]
        shape /= np.sum(shape)
        return total * shape.ravel()

    def capacity_J_K(self) -> np.ndarray:
        total = (
            self.geometry.total_stack_mass_kg
            * self.parameters.effective_total_thermal_mass_fraction
            * self.physics.graphite_heat_capacity_J_kgK
        )
        return np.full(self.state_count, total / self.state_count, dtype=float)

    def steady_state(
        self,
        current_A: float,
        *,
        duty: float | None = None,
        ambient_C: float | None = None,
        voltage_V: float | None = None,
    ) -> np.ndarray:
        duty = fan_duty_from_current(current_A) if duty is None else float(duty)
        ambient = self.physics.ambient_static_C if ambient_C is None else float(ambient_C)
        matrix = self._conductance_matrix(duty, ambient)
        theta = spsolve(-matrix, self.heat_generation_W(current_A, voltage_V))
        return np.asarray(theta + ambient).reshape(
            self.cells, self.flow_rows, self.in_plane_columns
        )

    def simulate_piecewise(
        self,
        times_s: np.ndarray,
        currents_A: np.ndarray,
        *,
        ambient_C: float,
        initial_temperature_C: float,
        duties: np.ndarray | None = None,
        sample_step_s: float = 1.0,
    ) -> tuple[np.ndarray, np.ndarray]:
        times = np.asarray(times_s, dtype=float)
        currents = np.asarray(currents_A, dtype=float)
        if len(times) != len(currents) or np.any(np.diff(times) <= 0):
            raise ValueError("times and currents must have equal length and increasing times")
        if duties is None:
            duties = np.asarray([fan_duty_from_current(value) for value in currents])
        else:
            duties = np.asarray(duties, dtype=float)
        capacity = self.capacity_J_K()
        inverse_capacity = diags(1.0 / capacity, format="csr")
        theta = np.full(self.state_count, initial_temperature_C - ambient_C, dtype=float)
        sample_times = np.arange(times[0], times[-1] + sample_step_s * 0.5, sample_step_s)
        output = np.empty((len(sample_times), self.state_count), dtype=float)
        output[0] = theta + ambient_C
        sample_index = 1
        start = perf_counter()
        for interval in range(len(times) - 1):
            t0, t1 = times[interval], times[interval + 1]
            matrix = self._conductance_matrix(float(duties[interval]), ambient_C)
            dynamic = inverse_capacity @ matrix
            source = self.heat_generation_W(float(currents[interval])) / capacity
            augmented = bmat(
                [
                    [dynamic, csr_matrix(source.reshape(-1, 1))],
                    [csr_matrix((1, self.state_count)), csr_matrix((1, 1))],
                ],
                format="csr",
            )
            augmented_state = np.concatenate((theta, (1.0,)))
            current_time = t0
            while sample_index < len(sample_times) and sample_times[sample_index] <= t1 + 1e-9:
                delta = sample_times[sample_index] - current_time
                if delta > 0.0:
                    augmented_state = expm_multiply(augmented * delta, augmented_state)
                output[sample_index] = augmented_state[:-1] + ambient_C
                current_time = sample_times[sample_index]
                sample_index += 1
            remaining = t1 - current_time
            if remaining > 1e-12:
                augmented_state = expm_multiply(augmented * remaining, augmented_state)
            theta = augmented_state[:-1]
        self.last_runtime_s = perf_counter() - start
        return sample_times, output.reshape(
            len(sample_times), self.cells, self.flow_rows, self.in_plane_columns
        )

    def sensor_field(self, temperature: np.ndarray) -> np.ndarray:
        """Return the five experimental group representatives as 5x3x5."""

        if temperature.shape != (self.cells, self.flow_rows, self.in_plane_columns):
            raise ValueError("unexpected temperature shape")
        scale = self.cells // 40
        base_indices = np.asarray((3, 11, 19, 27, 35))
        cell_indices = np.clip((base_indices + 0.5) * scale - 0.5, 0, self.cells - 1).astype(int)
        selected = temperature[cell_indices]
        if self.flow_rows == 3 and self.in_plane_columns == 5:
            return selected
        # Aggregate refinements to the physical 3x5 thermocouple regions.
        return selected.reshape(5, 3, self.flow_rows // 3, 5, self.in_plane_columns // 5).mean(axis=(2, 4))

    def energy_audit(
        self,
        temperature_C: np.ndarray,
        current_A: float,
        duty: float,
        ambient_C: float,
    ) -> dict[str, float]:
        theta = temperature_C.ravel() - ambient_C
        matrix = self._conductance_matrix(duty, ambient_C)
        q = self.heat_generation_W(current_A)
        residual_vector = q + matrix @ theta
        q_input = float(np.sum(q))
        return {
            "total_heat_input_W": q_input,
            "total_heat_rejection_W": float(-np.sum(matrix @ theta)),
            "maximum_nodal_residual_W": float(np.max(np.abs(residual_vector))),
            "energy_residual_percent": float(abs(np.sum(residual_vector)) / q_input * 100.0),
        }

    def frozen_record(self) -> dict[str, object]:
        return {
            "geometry": asdict(self.geometry),
            "fixed_physics": asdict(self.physics),
            "calibrated_parameters": asdict(self.parameters),
            "state_layout": {
                "cells": self.cells,
                "flow_rows": self.flow_rows,
                "in_plane_columns": self.in_plane_columns,
                "state_count": self.state_count,
            },
            "scope": "40-cell high-resolution thermal reference; electrochemistry reduced to measured voltage heat balance",
        }
