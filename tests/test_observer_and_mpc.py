from __future__ import annotations

import numpy as np

from controllers.hotspot_mpc import HotspotMPC
from models.five_node_model import FiveNodeParameters
from models.three_node_model import ThermalParameters
from observers.ekf import ExtendedKalmanFilter


def _parameters() -> FiveNodeParameters:
    three = ThermalParameters(5000.0, 30.0, 7.0, 0.07, 0.001)
    return FiveNodeParameters.from_three_node(three, np.array([25.0, 30.0, 35.0]))


def test_ekf_measurement_reduces_measured_node_error() -> None:
    observer = ExtendedKalmanFilter(
        state_C=np.full(5, 50.0), covariance=np.eye(5), process_covariance=np.eye(5) * 0.001,
        measurement_variance_C2=0.01, sensor_indices=(1, 4), parameters=_parameters(),
    )
    before = np.linalg.norm(observer.state_C[[1, 4]] - np.array([52.0, 48.0]))
    observer.update(np.array([52.0, 48.0]))
    after = np.linalg.norm(observer.state_C[[1, 4]] - np.array([52.0, 48.0]))
    assert after < before
    assert np.all(np.linalg.eigvalsh(observer.covariance) >= -1e-12)


def test_mpc_returns_feasible_command_and_obeys_dwell() -> None:
    controller = HotspotMPC(_parameters(), "A_thermoneutral_proxy", horizon_steps=3)
    duty, direction, cost = controller.command(np.array([56.0, 55.5, 55.0, 54.5, 54.0]), np.full(3, 20.0))
    assert 0.0 <= duty <= 1.0
    assert direction in (-1, 1)
    assert np.isfinite(cost)
    controller.steps_since_direction_change = 0
    controller.previous_direction = direction
    _, next_direction, _ = controller.command(np.array([54.0, 54.5, 55.0, 55.5, 56.0]), np.full(3, 20.0))
    assert next_direction == direction
