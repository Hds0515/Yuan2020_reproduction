"""Extended Kalman filter for the nonlinear five-node thermal model."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from models.five_node_model import FiveNodeParameters, step_five_node


@dataclass
class ExtendedKalmanFilter:
    state_C: np.ndarray
    covariance: np.ndarray
    process_covariance: np.ndarray
    measurement_variance_C2: float
    sensor_indices: tuple[int, ...]
    parameters: FiveNodeParameters
    heat_model: str = "A_thermoneutral_proxy"

    def _transition(self, state: np.ndarray, current_A: float, duty: float, direction: int, dt_s: float) -> np.ndarray:
        return step_five_node(state, current_A, duty, direction, dt_s, self.parameters, self.heat_model).next_temperature_C

    def predict(self, current_A: float, duty: float, direction: int, dt_s: float) -> None:
        prior = self.state_C.copy()
        prediction = self._transition(prior, current_A, duty, direction, dt_s)
        jacobian = np.empty((5, 5))
        epsilon = 1e-4
        for column in range(5):
            perturbed = prior.copy(); perturbed[column] += epsilon
            jacobian[:, column] = (self._transition(perturbed, current_A, duty, direction, dt_s) - prediction) / epsilon
        self.state_C = prediction
        self.covariance = jacobian @ self.covariance @ jacobian.T + self.process_covariance

    def update(self, measurement_C: np.ndarray) -> float:
        indices = np.asarray(self.sensor_indices, dtype=int)
        measurement = np.asarray(measurement_C, dtype=float)
        innovation = measurement - self.state_C[indices]
        h = np.zeros((len(indices), 5)); h[np.arange(len(indices)), indices] = 1.0
        r = np.eye(len(indices)) * self.measurement_variance_C2
        innovation_covariance = h @ self.covariance @ h.T + r
        gain = self.covariance @ h.T @ np.linalg.inv(innovation_covariance)
        self.state_C = self.state_C + gain @ innovation
        identity = np.eye(5)
        self.covariance = (identity - gain @ h) @ self.covariance @ (identity - gain @ h).T + gain @ r @ gain.T
        sign, logdet = np.linalg.slogdet(innovation_covariance)
        return float(-0.5 * (innovation @ np.linalg.solve(innovation_covariance, innovation) + logdet + len(indices) * np.log(2 * np.pi)))
