"""Physics-constrained five-zone ROM for the V15 fixed-60%-PWM reference.

The model is deliberately narrow: it represents five equal-volume zones along
the cathode-flow coordinate at 23 degC ambient and fixed 60 percent PWM.  It is
an energy-conserving control-oriented approximation of the V15 thermal
reference, not a replacement for the experimental high-fidelity gate.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.linalg import expm


ZONE_COUNT = 5
TOTAL_HEAT_CAPACITY_J_K = 3999.999989994387
THERMONEUTRAL_VOLTAGE_V_PER_CELL = 1.48
CELL_COUNT = 40
AXIAL_HEAT_SKEW_AT_40A = 0.4128642328339601


def _as_vector(values: np.ndarray, size: int, name: str) -> np.ndarray:
    array = np.asarray(values, dtype=float)
    if array.shape != (size,):
        raise ValueError(f"{name} must have shape ({size},)")
    return array


def conservative_five_zone_average(surface_temperature_C: np.ndarray) -> np.ndarray:
    """Conservatively remap an equal-width flow mesh to five zone averages.

    The trailing dimensions must be ``(cell, flow_row, in_plane_column)``.
    Overlap integration, rather than point interpolation, preserves the global
    surface-temperature integral for any number of flow rows.
    """

    surface = np.asarray(surface_temperature_C, dtype=float)
    if surface.ndim < 3:
        raise ValueError("surface temperature needs cell, row and column dimensions")
    rows = surface.shape[-2]
    row_mean = surface.mean(axis=(-3, -1))
    overlap = np.zeros((ZONE_COUNT, rows), dtype=float)
    for zone in range(ZONE_COUNT):
        zone_left, zone_right = zone / ZONE_COUNT, (zone + 1) / ZONE_COUNT
        for row in range(rows):
            row_left, row_right = row / rows, (row + 1) / rows
            overlap[zone, row] = max(
                0.0, min(zone_right, row_right) - max(zone_left, row_left)
            )
    weights = overlap * ZONE_COUNT
    result = np.einsum("...r,zr->...z", row_mean, weights)
    global_from_zones = result.mean(axis=-1)
    global_from_mesh = row_mean.mean(axis=-1)
    if not np.allclose(global_from_zones, global_from_mesh, atol=1e-10, rtol=1e-10):
        raise RuntimeError("five-zone remap did not preserve the global mean")
    return result


@dataclass(frozen=True)
class FiveZoneParameters:
    capacity_J_K: np.ndarray
    adjacent_conductance_W_K: np.ndarray
    ambient_loss_W_K: np.ndarray
    tmax_offset_C: np.ndarray
    tmin_offset_C: np.ndarray

    def __post_init__(self) -> None:
        capacity = _as_vector(self.capacity_J_K, 5, "capacity_J_K")
        conductance = _as_vector(
            self.adjacent_conductance_W_K, 4, "adjacent_conductance_W_K"
        )
        loss = _as_vector(self.ambient_loss_W_K, 5, "ambient_loss_W_K")
        maximum = _as_vector(self.tmax_offset_C, 2, "tmax_offset_C")
        minimum = _as_vector(self.tmin_offset_C, 2, "tmin_offset_C")
        if np.any(capacity <= 0) or not np.isclose(
            capacity.sum(), TOTAL_HEAT_CAPACITY_J_K, atol=1e-5
        ):
            raise ValueError("capacities must be positive and preserve V15 total heat capacity")
        if np.any(conductance < 0) or np.any(loss <= 0):
            raise ValueError("conductance must be nonnegative and ambient loss positive")
        if np.any(maximum < 0) or np.any(minimum < 0):
            raise ValueError("extreme-temperature reconstruction offsets must be nonnegative")

    def to_dict(self) -> dict[str, list[float]]:
        return {
            "capacity_J_K": self.capacity_J_K.tolist(),
            "adjacent_conductance_W_K": self.adjacent_conductance_W_K.tolist(),
            "ambient_loss_W_K": self.ambient_loss_W_K.tolist(),
            "tmax_offset_C_intercept_and_current40": self.tmax_offset_C.tolist(),
            "tmin_offset_C_intercept_and_current40": self.tmin_offset_C.tolist(),
        }


class PhysicsConstrainedFiveZoneROM:
    """Stable RC-chain ROM with conservative heat generation."""

    ambient_C = 23.0
    fixed_duty = 0.60

    def __init__(self, parameters: FiveZoneParameters) -> None:
        self.parameters = parameters
        self._continuous_matrix = self._build_continuous_matrix()

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
        if np.any(raw <= 0):
            raise ValueError("heat distribution became nonphysical")
        return raw / raw.sum()

    def _build_continuous_matrix(self) -> np.ndarray:
        p = self.parameters
        thermal = np.zeros((ZONE_COUNT, ZONE_COUNT), dtype=float)
        for edge, conductance in enumerate(p.adjacent_conductance_W_K):
            thermal[edge, edge] -= conductance
            thermal[edge + 1, edge + 1] -= conductance
            thermal[edge, edge + 1] += conductance
            thermal[edge + 1, edge] += conductance
        thermal[np.diag_indices(ZONE_COUNT)] -= p.ambient_loss_W_K
        return thermal / p.capacity_J_K[:, None]

    @property
    def continuous_matrix(self) -> np.ndarray:
        return self._continuous_matrix.copy()

    def stability_eigenvalues(self) -> np.ndarray:
        return np.linalg.eigvals(self._continuous_matrix)

    def assert_physical(self) -> None:
        matrix = self._continuous_matrix
        off_diagonal = matrix.copy()
        off_diagonal[np.diag_indices(ZONE_COUNT)] = 0.0
        if np.any(off_diagonal < -1e-12):
            raise RuntimeError("ROM is not a positive/Metzler thermal system")
        if np.max(np.real(self.stability_eigenvalues())) >= 0.0:
            raise RuntimeError("ROM is not asymptotically stable")

    def exact_discrete_transition(
        self, current_A: float, stack_voltage_V: float, dt_s: float
    ) -> tuple[np.ndarray, np.ndarray]:
        source_W = self.heat_input_W(current_A, stack_voltage_V) * self.heat_distribution(
            current_A
        )
        source_rate = source_W / self.parameters.capacity_J_K
        augmented = np.zeros((ZONE_COUNT + 1, ZONE_COUNT + 1), dtype=float)
        augmented[:ZONE_COUNT, :ZONE_COUNT] = self._continuous_matrix
        augmented[:ZONE_COUNT, ZONE_COUNT] = source_rate
        discrete = expm(augmented * float(dt_s))
        return discrete[:ZONE_COUNT, :ZONE_COUNT], discrete[:ZONE_COUNT, ZONE_COUNT]

    def step(
        self,
        state_C: np.ndarray,
        current_A: float,
        stack_voltage_V: float,
        dt_s: float,
    ) -> np.ndarray:
        state = _as_vector(state_C, 5, "state_C")
        transition, forcing = self.exact_discrete_transition(
            current_A, stack_voltage_V, dt_s
        )
        theta = state - self.ambient_C
        return self.ambient_C + transition @ theta + forcing

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
        if not (len(times) == len(current) == len(voltage)) or np.any(np.diff(times) <= 0):
            raise ValueError("times/current/voltage must align on a strictly increasing grid")
        output = np.empty((len(times), ZONE_COUNT), dtype=float)
        output[0] = _as_vector(initial_state_C, 5, "initial_state_C")
        for index, dt_s in enumerate(np.diff(times)):
            output[index + 1] = self.step(
                output[index], current[index], voltage[index], float(dt_s)
            )
        return output

    def extreme_temperatures(
        self, state_C: np.ndarray, current_A: float
    ) -> tuple[float, float, float]:
        state = _as_vector(state_C, 5, "state_C")
        normalized_current = np.clip(float(current_A) / 40.0, 0.0, 1.0)
        basis = np.asarray((1.0, normalized_current))
        maximum = float(np.max(state) + self.parameters.tmax_offset_C @ basis)
        minimum = float(np.min(state) - self.parameters.tmin_offset_C @ basis)
        return maximum, minimum, maximum - minimum

    def hotspot_zone(self, state_C: np.ndarray) -> int:
        return int(np.argmax(_as_vector(state_C, 5, "state_C")))

    def instantaneous_energy_audit(
        self, state_C: np.ndarray, current_A: float, stack_voltage_V: float
    ) -> dict[str, float]:
        state = _as_vector(state_C, 5, "state_C")
        theta = state - self.ambient_C
        heat = self.heat_input_W(current_A, stack_voltage_V)
        source = heat * self.heat_distribution(current_A)
        rate = self._continuous_matrix @ theta + source / self.parameters.capacity_J_K
        storage = float(self.parameters.capacity_J_K @ rate)
        rejection = float(self.parameters.ambient_loss_W_K @ theta)
        residual = storage - (heat - rejection)
        return {
            "heat_input_W": heat,
            "ambient_rejection_W": rejection,
            "stored_heat_rate_W": storage,
            "energy_residual_W": residual,
            "energy_residual_percent": abs(residual) / max(heat, 1e-12) * 100.0,
        }
