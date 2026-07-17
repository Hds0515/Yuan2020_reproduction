"""Galerkin-projected five-zone, two-mass ROM for the frozen V15 plant."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm
from scipy.sparse import csr_matrix

from models_v15 import Fixed60HighResolutionReference
from models_v16.physics_constrained_five_zone_rom import (
    AXIAL_HEAT_SKEW_AT_40A,
    CELL_COUNT,
    THERMONEUTRAL_VOLTAGE_V_PER_CELL,
    TOTAL_HEAT_CAPACITY_J_K,
    ZONE_COUNT,
)


@dataclass(frozen=True)
class ProjectedParameters:
    capacity_J_K: np.ndarray
    thermal_operator_W_K: np.ndarray
    tmax_offset_C: np.ndarray
    tmin_offset_C: np.ndarray

    def __post_init__(self) -> None:
        capacity = np.asarray(self.capacity_J_K, dtype=float)
        operator = np.asarray(self.thermal_operator_W_K, dtype=float)
        if capacity.shape != (10,) or operator.shape != (10, 10):
            raise ValueError("projected model requires 10 capacities and a 10x10 operator")
        if np.any(capacity <= 0.0) or not np.isclose(
            capacity.sum(), TOTAL_HEAT_CAPACITY_J_K, atol=1e-5
        ):
            raise ValueError("projected capacities must preserve the V15 total")
        off_diagonal = operator.copy()
        off_diagonal[np.diag_indices(10)] = 0.0
        if np.any(off_diagonal < -1e-10):
            raise ValueError("projected thermal operator must be Metzler")
        if np.asarray(self.tmax_offset_C).shape != (2,) or np.asarray(
            self.tmin_offset_C
        ).shape != (2,):
            raise ValueError("extreme reconstruction needs two coefficients")

    def to_dict(self) -> dict[str, object]:
        return {
            "capacity_J_K": self.capacity_J_K.tolist(),
            "thermal_operator_W_K": self.thermal_operator_W_K.tolist(),
            "tmax_offset_C_intercept_and_current40": self.tmax_offset_C.tolist(),
            "tmin_offset_C_intercept_and_current40": self.tmin_offset_C.tolist(),
        }


class ProjectedTwoMassFiveZoneROM:
    """Ten-state positive ROM: surface and core state for each of five zones."""

    ambient_C = 23.0
    fixed_duty = 0.60
    state_count = 10

    def __init__(self, parameters: ProjectedParameters) -> None:
        self.parameters = parameters
        self._continuous = parameters.thermal_operator_W_K / parameters.capacity_J_K[:, None]

    @classmethod
    def project_from_v15(
        cls,
        plant: Fixed60HighResolutionReference,
        *,
        tmax_offset_C: np.ndarray | None = None,
        tmin_offset_C: np.ndarray | None = None,
    ) -> "ProjectedTwoMassFiveZoneROM":
        if plant.flow_rows % ZONE_COUNT:
            raise ValueError("projection mesh flow rows must be divisible by five")
        surface_count = plant.surface_count
        full_count = plant.state_count
        rows_per_zone = plant.flow_rows // ZONE_COUNT
        row_indices: list[int] = []
        column_indices: list[int] = []
        for cell in range(plant.cells):
            for row in range(plant.flow_rows):
                zone = row // rows_per_zone
                for column in range(plant.columns):
                    node = plant.index(cell, row, column)
                    row_indices.extend((node, surface_count + node))
                    column_indices.extend((zone, ZONE_COUNT + zone))
        prolongation = csr_matrix(
            (
                np.ones(len(row_indices), dtype=float),
                (np.asarray(row_indices), np.asarray(column_indices)),
            ),
            shape=(full_count, 10),
        )
        cs, cc = plant.capacities_J_K()
        full_capacity = np.concatenate((cs, cc))
        capacity = np.asarray(prolongation.T @ full_capacity).ravel()
        full_operator = plant._system_matrix(0.60, 23.0)
        operator = np.asarray((prolongation.T @ full_operator @ prolongation).todense())
        parameters = ProjectedParameters(
            capacity,
            operator,
            np.zeros(2) if tmax_offset_C is None else np.asarray(tmax_offset_C),
            np.zeros(2) if tmin_offset_C is None else np.asarray(tmin_offset_C),
        )
        model = cls(parameters)
        model.assert_physical()
        return model

    @staticmethod
    def heat_input_W(current_A: float, stack_voltage_V: float) -> float:
        return max(
            float(current_A)
            * (THERMONEUTRAL_VOLTAGE_V_PER_CELL * CELL_COUNT - float(stack_voltage_V)),
            0.0,
        )

    @staticmethod
    def heat_distribution(current_A: float) -> np.ndarray:
        centers = (np.arange(ZONE_COUNT, dtype=float) + 0.5) / ZONE_COUNT
        strength = AXIAL_HEAT_SKEW_AT_40A * max(float(current_A), 0.0) / 40.0
        raw = 1.0 + strength * (2.0 * centers - 1.0)
        return raw / raw.sum()

    @property
    def continuous_matrix(self) -> np.ndarray:
        return self._continuous.copy()

    def assert_physical(self) -> None:
        off_diagonal = self._continuous.copy()
        off_diagonal[np.diag_indices(10)] = 0.0
        if np.any(off_diagonal < -1e-12):
            raise RuntimeError("projected dynamics lost positivity")
        if np.max(np.real(np.linalg.eigvals(self._continuous))) >= 0.0:
            raise RuntimeError("projected dynamics are not asymptotically stable")

    def exact_discrete_transition(
        self, current_A: float, stack_voltage_V: float, dt_s: float
    ) -> tuple[np.ndarray, np.ndarray]:
        heat = self.heat_input_W(current_A, stack_voltage_V)
        source_W = np.r_[heat * self.heat_distribution(current_A), np.zeros(5)]
        source_rate = source_W / self.parameters.capacity_J_K
        augmented = np.zeros((11, 11), dtype=float)
        augmented[:10, :10] = self._continuous
        augmented[:10, 10] = source_rate
        discrete = expm(augmented * float(dt_s))
        return discrete[:10, :10], discrete[:10, 10]

    def simulate(
        self,
        times_s: np.ndarray,
        currents_A: np.ndarray,
        voltages_V: np.ndarray,
        initial_state_C: np.ndarray,
    ) -> np.ndarray:
        times = np.asarray(times_s, dtype=float)
        current = np.asarray(currents_A, dtype=float)
        voltage = np.asarray(voltages_V, dtype=float)
        initial = np.asarray(initial_state_C, dtype=float)
        if initial.shape != (10,) or not (
            len(times) == len(current) == len(voltage)
        ) or np.any(np.diff(times) <= 0):
            raise ValueError("invalid projected-ROM simulation inputs")
        output = np.empty((len(times), 10), dtype=float)
        output[0] = initial
        for index, dt_s in enumerate(np.diff(times)):
            transition, forcing = self.exact_discrete_transition(
                current[index], voltage[index], float(dt_s)
            )
            output[index + 1] = self.ambient_C + transition @ (
                output[index] - self.ambient_C
            ) + forcing
        return output

    def surface_zones(self, state_C: np.ndarray) -> np.ndarray:
        state = np.asarray(state_C, dtype=float)
        return state[..., :5]

    def extreme_temperatures(
        self, state_C: np.ndarray, current_A: float
    ) -> tuple[float, float, float]:
        surface = self.surface_zones(state_C)
        basis = np.asarray((1.0, np.clip(float(current_A) / 40.0, 0.0, 1.0)))
        maximum = float(np.max(surface) + self.parameters.tmax_offset_C @ basis)
        minimum = float(np.min(surface) - self.parameters.tmin_offset_C @ basis)
        return maximum, minimum, maximum - minimum

    def instantaneous_energy_audit(
        self, state_C: np.ndarray, current_A: float, stack_voltage_V: float
    ) -> dict[str, float]:
        state = np.asarray(state_C, dtype=float)
        theta = state - self.ambient_C
        heat = self.heat_input_W(current_A, stack_voltage_V)
        source = np.r_[heat * self.heat_distribution(current_A), np.zeros(5)]
        rate = self._continuous @ theta + source / self.parameters.capacity_J_K
        storage = float(self.parameters.capacity_J_K @ rate)
        rejection = float(-np.sum(self.parameters.thermal_operator_W_K @ theta))
        residual = storage - (heat - rejection)
        return {
            "heat_input_W": heat,
            "heat_rejection_W": rejection,
            "stored_heat_rate_W": storage,
            "energy_residual_W": residual,
            "energy_residual_percent": abs(residual) / max(heat, 1e-12) * 100.0,
        }


def project_full_state_to_ten_zones(
    plant: Fixed60HighResolutionReference, state_C: np.ndarray
) -> np.ndarray:
    """Average a divisible flow mesh into five surface and five core zones."""

    if plant.flow_rows % ZONE_COUNT:
        raise ValueError("projection mesh flow rows must be divisible by five")
    state = np.asarray(state_C, dtype=float)
    surface = state[..., : plant.surface_count].reshape(
        *state.shape[:-1], plant.cells, plant.flow_rows, plant.columns
    )
    core = state[..., plant.surface_count :].reshape(
        *state.shape[:-1], plant.cells, plant.flow_rows, plant.columns
    )
    rows_per_zone = plant.flow_rows // ZONE_COUNT

    def zones(field: np.ndarray) -> np.ndarray:
        leading = field.shape[:-3]
        reshaped = field.reshape(
            *leading,
            plant.cells,
            ZONE_COUNT,
            rows_per_zone,
            plant.columns,
        )
        return reshaped.mean(axis=(-4, -2, -1))

    return np.concatenate((zones(surface), zones(core)), axis=-1)
