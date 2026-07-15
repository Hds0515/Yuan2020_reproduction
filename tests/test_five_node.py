from __future__ import annotations

import numpy as np

from models.five_node_model import FiveNodeParameters, step_five_node
from models.three_node_model import ThermalParameters


def _parameters() -> FiveNodeParameters:
    three = ThermalParameters(5000.0, 30.0, 7.0, 0.07, 0.001)
    return FiveNodeParameters.from_three_node(three, np.array([25.0, 30.0, 35.0]))


def test_conduction_is_dissipative_and_conservative() -> None:
    initial = np.array([50.0, 50.0, 60.0, 50.0, 50.0])
    result = step_five_node(initial, 0.0, 0.0, 1, 1.0, _parameters(), "A_thermoneutral_proxy")
    assert result.next_temperature_C[2] < initial[2]
    assert result.next_temperature_C[1] > initial[1]
    assert result.next_temperature_C[3] > initial[3]
    assert abs(result.energy_residual_W) < 1e-9


def test_five_node_scaling_preserves_capacity_and_cooled_area() -> None:
    parameters = _parameters()
    assert parameters.total_heat_capacity_J_K == 5000.0
    assert np.isclose(5 * parameters.cooling_coefficient_per_node_W_K, 3 * 30.0)
