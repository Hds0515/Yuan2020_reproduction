"""Mesh-consistent thermal digital twin of the Gong et al. 40-cell stack.

The plant resolves 40 stack positions, cathode-flow and in-plane directions,
and two thermal masses per finite volume.  It is deliberately a thermal model:
measured current and stack voltage close the heat balance, while voltage and
water state are not predicted.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from time import perf_counter

import numpy as np
from scipy.sparse import bmat, csr_matrix, diags, lil_matrix
from scipy.sparse.linalg import expm_multiply, spsolve


@dataclass(frozen=True)
class Geometry:
    cell_count: int = 40
    active_area_m2: float = 0.0051
    cathode_channel_count: int = 45
    channel_depth_m: float = 0.0014
    channel_width_m: float = 0.0010
    channel_length_m: float = 0.0608
    ridge_width_m: float = 0.0014
    bipolar_plate_thickness_m: float = 0.0025
    mea_thickness_m: float = 0.00062
    stack_length_m: float = 0.13642

    @property
    def active_width_m(self) -> float:
        return self.active_area_m2 / self.channel_length_m


@dataclass(frozen=True)
class FixedPhysics:
    graphite_in_plane_k_W_mK: float = 60.0
    air_cp_J_kgK: float = 1006.0
    emissivity: float = 0.82
    stefan_boltzmann_W_m2K4: float = 5.670374419e-8
    reversible_voltage_V_per_cell: float = 1.253
    reference_duty: float = 0.57
    fan_start_duty: float = 0.12
    minimum_commanded_duty: float = 0.15
    fan_start_flow_fraction: float = 0.08
    reference_internal_h_W_m2K: float = 85.0
    plenum_mixing_fraction: float = 0.55


@dataclass(frozen=True)
class Parameters:
    effective_total_heat_capacity_J_K: float = 1350.0
    surface_heat_capacity_fraction: float = 0.32
    surface_core_exchange_W_K_total: float = 65.0
    effective_total_airflow_at_57pct_kg_s: float = 0.030
    pwm_to_flow_exponent: float = 2.0
    channel_heat_transfer_multiplier: float = 1.2
    stack_edge_heat_transfer_coefficient_W_m2K: float = 7.0
    effective_through_stack_conductivity_W_mK: float = 1.2
    heat_source_nonuniformity: float = 0.35
    fan_hub_deficit: float = 0.22
    fan_release_time_constant_s: float = 35.0


PARAMETER_NAMES = tuple(Parameters.__dataclass_fields__)
LOWER = np.asarray((700.0, 0.10, 12.0, 0.012, 0.8, 0.55, 3.0, 0.15, 0.0, 0.0, 5.0))
UPPER = np.asarray((2600.0, 0.65, 180.0, 0.070, 3.5, 2.4, 14.0, 4.0, 0.75, 0.45, 120.0))
PRIOR = np.asarray((1350.0, 0.32, 65.0, 0.030, 2.0, 1.2, 7.0, 1.2, 0.35, 0.22, 35.0))
PRIOR_SCALE = np.asarray((500.0, 0.16, 45.0, 0.018, 0.8, 0.5, 3.5, 1.0, 0.24, 0.15, 30.0))


def vector_to_parameters(values: np.ndarray) -> Parameters:
    return Parameters(**dict(zip(PARAMETER_NAMES, map(float, values), strict=True)))


def fan_duty_from_current(current_A: float) -> float:
    current = max(float(current_A), 0.0)
    return float(np.clip(0.15 + 0.015 * max(current - 2.0, 0.0), 0.15, 0.80))


def fan_duty_with_release_lag(
    times_s: np.ndarray,
    currents_A: np.ndarray,
    release_time_constant_s: float,
) -> np.ndarray:
    """Reconstruct duty with immediate rise and delayed release.

    Gong et al. identify fan/current timing as the cause of the observed
    dynamic undershoot and recommend lagging fan adjustment behind current.
    PWM telemetry is not published, so this bounded first-order release is an
    explicitly identified actuator-input model rather than a measurement.
    """

    times = np.asarray(times_s, dtype=float)
    targets = np.asarray([fan_duty_from_current(value) for value in currents_A])
    duty = np.empty_like(targets)
    duty[0] = targets[0]
    tau = max(float(release_time_constant_s), 1e-6)
    for index in range(1, len(times)):
        if targets[index] >= duty[index - 1]:
            duty[index] = targets[index]
        else:
            decay = np.exp(-(times[index] - times[index - 1]) / tau)
            duty[index] = targets[index] + (duty[index - 1] - targets[index]) * decay
    return duty


class Gong40CellDigitalTwin:
    """Finite-volume thermal plant with exact air-path energy removal."""

    def __init__(
        self,
        parameters: Parameters | None = None,
        *,
        cells: int = 40,
        flow_rows: int = 3,
        columns: int = 5,
    ) -> None:
        if cells != 40 or flow_rows % 3 or columns % 5:
            raise ValueError("cells must be 40; in-plane mesh must refine 3x5")
        self.geometry = Geometry()
        self.physics = FixedPhysics()
        self.parameters = parameters or Parameters()
        self.cells = cells
        self.flow_rows = flow_rows
        self.columns = columns
        self.surface_count = cells * flow_rows * columns
        self.state_count = 2 * self.surface_count
        self._surface_cache: dict[tuple[float, float], csr_matrix] = {}
        self._system_cache: dict[tuple[float, float], csr_matrix] = {}

    def index(self, cell: int, row: int, column: int) -> int:
        return (cell * self.flow_rows + row) * self.columns + column

    @property
    def cell_sizes(self) -> tuple[float, float, float]:
        return (
            self.geometry.active_width_m / self.columns,
            self.geometry.channel_length_m / self.flow_rows,
            self.geometry.stack_length_m / self.cells,
        )

    def _fan_weights(self) -> np.ndarray:
        p = self.parameters
        z = (np.arange(self.cells) + 0.5) / self.cells
        x = (np.arange(self.columns) + 0.5) / self.columns
        radial_z = np.abs(z - 0.5) / 0.5
        radial_x = np.abs(x - 0.5) / 0.5
        radius = np.sqrt(np.square(radial_z[:, None]) + 0.35 * np.square(radial_x[None, :]))
        raw = 1.0 - p.fan_hub_deficit * np.exp(-np.square(radius / 0.38))
        mixed = 1.0 + (raw / np.mean(raw) - 1.0) * (1.0 - self.physics.plenum_mixing_fraction)
        return mixed / np.mean(mixed)

    def total_airflow_kg_s(self, duty: float) -> float:
        p = self.parameters
        fixed = self.physics
        # The paper distinguishes the 12% electrical start threshold from the
        # 15% initial commanded operating point.  The latter is the fan-curve
        # knee; below it the speed/flow is not reported as stable.
        normalized = np.clip(
            (float(duty) - fixed.minimum_commanded_duty)
            / (fixed.reference_duty - fixed.minimum_commanded_duty),
            0.0,
            None,
        )
        factor = fixed.fan_start_flow_fraction + (1.0 - fixed.fan_start_flow_fraction) * normalized**p.pwm_to_flow_exponent
        return float(p.effective_total_airflow_at_57pct_kg_s * factor)

    def _flow_distribution(self, duty: float) -> np.ndarray:
        weights = self._fan_weights()
        return self.total_airflow_kg_s(duty) * weights / np.sum(weights)

    def _surface_matrix(self, duty: float, ambient_C: float) -> csr_matrix:
        key = (round(float(duty), 7), round(float(ambient_C), 5))
        if key in self._surface_cache:
            return self._surface_cache[key]
        g, f, p = self.geometry, self.physics, self.parameters
        dx, dy, dz = self.cell_sizes
        matrix = lil_matrix((self.surface_count, self.surface_count), dtype=float)

        # Mesh-invariant face conductances.
        gx = f.graphite_in_plane_k_W_mK * g.bipolar_plate_thickness_m * dy / dx
        gy = f.graphite_in_plane_k_W_mK * g.bipolar_plate_thickness_m * dx / dy
        gz = p.effective_through_stack_conductivity_W_mK * dx * dy / dz

        def connect(a: int, b: int, conductance: float) -> None:
            matrix[a, a] -= conductance
            matrix[b, b] -= conductance
            matrix[a, b] += conductance
            matrix[b, a] += conductance

        for cell in range(self.cells):
            for row in range(self.flow_rows):
                for column in range(self.columns):
                    node = self.index(cell, row, column)
                    if column + 1 < self.columns:
                        connect(node, self.index(cell, row, column + 1), gx)
                    if row + 1 < self.flow_rows:
                        connect(node, self.index(cell, row + 1, column), gy)
                    if cell + 1 < self.cells:
                        connect(node, self.index(cell + 1, row, column), gz)

        # Air enthalpy is marched conservatively through each cell/column path.
        path_flow = self._flow_distribution(duty)
        mean_flow = float(np.mean(path_flow))
        channels_per_path = g.cathode_channel_count / self.columns
        wetted_area = channels_per_path * (g.channel_width_m + 2.0 * g.channel_depth_m) * dy
        for cell in range(self.cells):
            for column in range(self.columns):
                mdot = max(float(path_flow[cell, column]), 1e-10)
                h = f.reference_internal_h_W_m2K * p.channel_heat_transfer_multiplier * (mdot / mean_flow) ** 0.33
                effectiveness = 1.0 - np.exp(-h * wetted_area / (mdot * f.air_cp_J_kgK))
                mcp = mdot * f.air_cp_J_kgK
                upstream_weights: dict[int, float] = {}
                for row in range(self.flow_rows):
                    node = self.index(cell, row, column)
                    removal = mcp * effectiveness
                    matrix[node, node] -= removal
                    for upstream, weight in upstream_weights.items():
                        matrix[node, upstream] += removal * weight
                    upstream_weights = {
                        upstream: (1.0 - effectiveness) * weight
                        for upstream, weight in upstream_weights.items()
                    }
                    upstream_weights[node] = upstream_weights.get(node, 0.0) + effectiveness

        ambient_K = ambient_C + 273.15
        reference_K = 55.0 + 273.15
        h_rad = f.emissivity * f.stefan_boltzmann_W_m2K4 * (
            reference_K**2 + ambient_K**2
        ) * (reference_K + ambient_K)
        h_external = p.stack_edge_heat_transfer_coefficient_W_m2K + h_rad
        for cell in range(self.cells):
            for row in range(self.flow_rows):
                for column in range(self.columns):
                    area = 0.0
                    if column in (0, self.columns - 1):
                        area += dz * dy
                    if row in (0, self.flow_rows - 1):
                        area += dz * dx
                    if cell in (0, self.cells - 1):
                        area += dx * dy
                    node = self.index(cell, row, column)
                    matrix[node, node] -= h_external * area

        result = matrix.tocsr()
        self._surface_cache[key] = result
        return result

    def heat_generation_W(self, current_A: float, voltage_V: float) -> np.ndarray:
        current = max(float(current_A), 0.0)
        total = max(current * (self.physics.reversible_voltage_V_per_cell * 40.0 - float(voltage_V)), 0.0)
        base_pattern = np.asarray((-0.65, 0.10, 0.28, 0.52, -0.25), dtype=float)
        base_x = (np.arange(5) + 0.5) / 5.0
        fine_x = (np.arange(self.columns) + 0.5) / self.columns
        pattern = np.interp(fine_x, base_x, base_pattern)
        shape = np.ones((self.cells, self.flow_rows, self.columns), dtype=float)
        shape *= 1.0 + self.parameters.heat_source_nonuniformity * pattern[None, None, :]
        shape /= np.sum(shape)
        return total * shape.ravel()

    def capacities_J_K(self) -> tuple[np.ndarray, np.ndarray]:
        total = self.parameters.effective_total_heat_capacity_J_K
        fraction = self.parameters.surface_heat_capacity_fraction
        return (
            np.full(self.surface_count, total * fraction / self.surface_count),
            np.full(self.surface_count, total * (1.0 - fraction) / self.surface_count),
        )

    def _system_matrix(self, duty: float, ambient_C: float) -> csr_matrix:
        key = (round(float(duty), 7), round(float(ambient_C), 5))
        if key in self._system_cache:
            return self._system_cache[key]
        surface = self._surface_matrix(duty, ambient_C)
        per_node = self.parameters.surface_core_exchange_W_K_total / self.surface_count
        exchange = diags(np.full(self.surface_count, per_node), format="csr")
        system = bmat(
            [[surface - exchange, exchange], [exchange, -exchange]],
            format="csr",
        )
        self._system_cache[key] = system
        return system

    def steady_state(
        self,
        current_A: float,
        voltage_V: float,
        *,
        duty: float,
        ambient_C: float,
    ) -> tuple[np.ndarray, np.ndarray]:
        system = self._system_matrix(duty, ambient_C)
        source = np.concatenate((self.heat_generation_W(current_A, voltage_V), np.zeros(self.surface_count)))
        theta = spsolve(-system, source)
        surface = theta[: self.surface_count].reshape(self.cells, self.flow_rows, self.columns) + ambient_C
        core = theta[self.surface_count :].reshape(self.cells, self.flow_rows, self.columns) + ambient_C
        return surface, core

    def simulate_anchors(
        self,
        times_s: np.ndarray,
        currents_A: np.ndarray,
        voltages_V: np.ndarray,
        *,
        ambient_C: float,
        initial_temperature_C: float,
        duties: np.ndarray | None = None,
    ) -> np.ndarray:
        times = np.asarray(times_s, dtype=float)
        currents = np.asarray(currents_A, dtype=float)
        voltages = np.asarray(voltages_V, dtype=float)
        if not (len(times) == len(currents) == len(voltages)) or np.any(np.diff(times) <= 0):
            raise ValueError("times, currents and voltages must align on an increasing grid")
        if duties is None:
            duties = np.asarray([fan_duty_from_current(value) for value in currents])
        else:
            duties = np.asarray(duties, dtype=float)
        cs, cc = self.capacities_J_K()
        inverse_capacity = diags(1.0 / np.concatenate((cs, cc)), format="csr")
        theta = np.full(self.state_count, initial_temperature_C - ambient_C)
        output = np.empty((len(times), self.state_count), dtype=float)
        output[0] = theta + ambient_C
        start = perf_counter()
        for interval in range(len(times) - 1):
            system = self._system_matrix(float(duties[interval]), ambient_C)
            dynamic = inverse_capacity @ system
            q = self.heat_generation_W(currents[interval], voltages[interval])
            source = np.concatenate((q / cs, np.zeros(self.surface_count)))
            augmented = bmat(
                [
                    [dynamic, csr_matrix(source.reshape(-1, 1))],
                    [csr_matrix((1, self.state_count)), csr_matrix((1, 1))],
                ],
                format="csr",
            )
            augmented_state = expm_multiply(
                augmented * (times[interval + 1] - times[interval]),
                np.concatenate((theta, (1.0,))),
            )
            theta = augmented_state[:-1]
            output[interval + 1] = theta + ambient_C
        self.last_runtime_s = perf_counter() - start
        return output

    def surface_from_state(self, state: np.ndarray) -> np.ndarray:
        return state[..., : self.surface_count].reshape(
            *state.shape[:-1], self.cells, self.flow_rows, self.columns
        )

    def sensor_field(self, surface_temperature: np.ndarray) -> np.ndarray:
        expected = (self.cells, self.flow_rows, self.columns)
        if surface_temperature.shape[-3:] != expected:
            raise ValueError(f"expected trailing dimensions {expected}")
        base_indices = np.asarray((3, 11, 19, 27, 35))
        selected = surface_temperature[..., base_indices, :, :]
        if self.flow_rows == 3 and self.columns == 5:
            return selected
        leading = selected.shape[:-3]
        reshaped = selected.reshape(
            *leading,
            5,
            3,
            self.flow_rows // 3,
            5,
            self.columns // 5,
        )
        return reshaped.mean(axis=(-3, -1))

    def energy_audit(
        self,
        surface_C: np.ndarray,
        core_C: np.ndarray,
        current_A: float,
        voltage_V: float,
        duty: float,
        ambient_C: float,
    ) -> dict[str, float]:
        theta = np.concatenate((surface_C.ravel() - ambient_C, core_C.ravel() - ambient_C))
        system = self._system_matrix(duty, ambient_C)
        q = self.heat_generation_W(current_A, voltage_V)
        residual = np.concatenate((q, np.zeros(self.surface_count))) + system @ theta
        q_total = float(np.sum(q))
        return {
            "total_heat_input_W": q_total,
            "total_heat_rejection_W": float(-np.sum(system @ theta)),
            "maximum_nodal_residual_W": float(np.max(np.abs(residual))),
            "energy_residual_percent": float(abs(np.sum(residual)) / max(q_total, 1e-12) * 100.0),
        }

    def record(self) -> dict[str, object]:
        return {
            "geometry": asdict(self.geometry),
            "fixed_physics": asdict(self.physics),
            "parameters": asdict(self.parameters),
            "mesh": {
                "cells": self.cells,
                "flow_rows": self.flow_rows,
                "columns": self.columns,
                "surface_states": self.surface_count,
                "total_states": self.state_count,
            },
            "scope": "thermal digital twin; measured current and voltage are exogenous inputs",
        }
