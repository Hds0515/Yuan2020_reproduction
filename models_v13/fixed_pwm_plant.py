"""V13 fixed-60%-PWM thermal plant built on frozen V12 spatial physics."""

from __future__ import annotations

from dataclasses import asdict, dataclass

import numpy as np
from scipy.sparse import bmat, csr_matrix, diags
from scipy.sparse.linalg import expm_multiply

from models_v12.gong40cell_digital_twin import (
    PARAMETER_NAMES as V12_PARAMETER_NAMES,
    Gong40CellDigitalTwin,
    Parameters as V12Parameters,
)


V12_FROZEN_SHA256 = "35cc2b8edce569ee18e01322baf88f2b21d74a42044ff40150fb4cf2c5c06f95"
V12_SPATIAL = {
    "effective_total_airflow_at_57pct_kg_s": 0.021934345212215003,
    "pwm_to_flow_exponent": 2.5917841202975773,
    "channel_heat_transfer_multiplier": 1.5434528081728003,
    "stack_edge_heat_transfer_coefficient_W_m2K": 13.9999972210056,
    "effective_through_stack_conductivity_W_mK": 2.8745487024696486,
    "heat_source_nonuniformity": 0.3721976849771562,
    "fan_hub_deficit": 0.449999,
}


@dataclass(frozen=True)
class DynamicParameters:
    effective_total_heat_capacity_J_K: float = 1800.0
    surface_heat_capacity_fraction: float = 0.12
    surface_core_exchange_W_K_total: float = 40.0
    fixed60_airflow_multiplier_relative_to_v12: float = 1.5
    current_dependent_axial_heat_skew: float = 0.25


NAMES = tuple(DynamicParameters.__dataclass_fields__)
LOWER = np.asarray((700.0, 0.01, 5.0, 0.5, 0.0))
UPPER = np.asarray((4000.0, 0.65, 180.0, 3.0, 0.60))
PRIOR = np.asarray((1800.0, 0.12, 40.0, 1.5, 0.25))
SCALE = np.asarray((900.0, 0.16, 45.0, 0.6, 0.20))


def vector_to_dynamic(values: np.ndarray) -> DynamicParameters:
    return DynamicParameters(**dict(zip(NAMES, map(float, values), strict=True)))


class FixedPWM60Plant(Gong40CellDigitalTwin):
    thermoneutral_voltage_V_per_cell = 1.48

    def __init__(
        self,
        dynamic: DynamicParameters | None = None,
        *,
        flow_rows: int = 3,
        columns: int = 5,
    ) -> None:
        self.dynamic_parameters = dynamic or DynamicParameters()
        values = {name: getattr(V12Parameters(), name) for name in V12_PARAMETER_NAMES}
        values.update(V12_SPATIAL)
        values.update(
            {
                "effective_total_heat_capacity_J_K": self.dynamic_parameters.effective_total_heat_capacity_J_K,
                "surface_heat_capacity_fraction": self.dynamic_parameters.surface_heat_capacity_fraction,
                "surface_core_exchange_W_K_total": self.dynamic_parameters.surface_core_exchange_W_K_total,
                "effective_total_airflow_at_57pct_kg_s": V12_SPATIAL["effective_total_airflow_at_57pct_kg_s"]
                * self.dynamic_parameters.fixed60_airflow_multiplier_relative_to_v12,
            }
        )
        super().__init__(V12Parameters(**values), flow_rows=flow_rows, columns=columns)

    def heat_generation_W(self, current_A: float, voltage_V: float) -> np.ndarray:
        current = max(float(current_A), 0.0)
        total = max(
            current * (self.thermoneutral_voltage_V_per_cell * 40.0 - float(voltage_V)),
            0.0,
        )
        base_column = np.asarray((-0.65, 0.10, 0.28, 0.52, -0.25))
        column_pattern = np.interp(
            (np.arange(self.columns) + 0.5) / self.columns,
            (np.arange(5) + 0.5) / 5.0,
            base_column,
        )
        row_coordinate = 2.0 * (np.arange(self.flow_rows) + 0.5) / self.flow_rows - 1.0
        axial_strength = self.dynamic_parameters.current_dependent_axial_heat_skew * current / 40.0
        shape = np.ones((self.cells, self.flow_rows, self.columns))
        shape *= 1.0 + self.parameters.heat_source_nonuniformity * column_pattern[None, None, :]
        shape *= 1.0 + axial_strength * row_coordinate[None, :, None]
        shape /= np.sum(shape)
        return total * shape.ravel()

    def simulate_from_regional_initial_state(
        self,
        times_s: np.ndarray,
        currents_A: np.ndarray,
        voltages_V: np.ndarray,
        initial_regions_C: np.ndarray,
        *,
        ambient_C: float = 23.0,
        duty: float = 0.60,
    ) -> np.ndarray:
        times = np.asarray(times_s, dtype=float)
        current = np.asarray(currents_A, dtype=float)
        voltage = np.asarray(voltages_V, dtype=float)
        regions = np.asarray(initial_regions_C, dtype=float)
        if regions.shape != (3,):
            raise ValueError("initial regions must be [inlet, intermediate, outlet]")
        surface = np.empty((self.cells, self.flow_rows, self.columns))
        refinement = self.flow_rows // 3
        for physical_row in range(3):
            surface[:, physical_row * refinement : (physical_row + 1) * refinement, :] = regions[physical_row]
        theta = np.concatenate((surface.ravel(), surface.ravel())) - ambient_C
        cs, cc = self.capacities_J_K()
        inverse_capacity = diags(1.0 / np.concatenate((cs, cc)), format="csr")
        output = np.empty((len(times), self.state_count))
        output[0] = theta + ambient_C
        dynamic_matrix = inverse_capacity @ self._system_matrix(duty, ambient_C)
        for interval in range(len(times) - 1):
            q = self.heat_generation_W(current[interval], voltage[interval])
            source = np.concatenate((q / cs, np.zeros(self.surface_count)))
            augmented = bmat(
                [
                    [dynamic_matrix, csr_matrix(source.reshape(-1, 1))],
                    [csr_matrix((1, self.state_count)), csr_matrix((1, 1))],
                ],
                format="csr",
            )
            theta = expm_multiply(
                augmented * (times[interval + 1] - times[interval]),
                np.concatenate((theta, (1.0,))),
            )[:-1]
            output[interval + 1] = theta + ambient_C
        return output

    def v13_record(self) -> dict[str, object]:
        return {
            "v12_frozen_sha256": V12_FROZEN_SHA256,
            "v12_spatial_parameters": V12_SPATIAL,
            "dynamic_parameters": asdict(self.dynamic_parameters),
            "thermoneutral_voltage_V_per_cell": self.thermoneutral_voltage_V_per_cell,
            "state_count": self.state_count,
            "scope": "23 C ambient, 2-40 A, fixed 60 percent PWM thermal plant",
        }
