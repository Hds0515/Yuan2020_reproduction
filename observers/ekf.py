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

@dataclass
class FastExtendedKalmanFilter(ExtendedKalmanFilter):
    """Five-node EKF with an analytic state-transition Jacobian.

    It is mathematically equivalent to the finite-difference implementation for
    the current five-node model but is substantially faster for Monte Carlo
    validation.
    """

    def _analytic_jacobian(
        self, state: np.ndarray, current_A: float, duty: float, direction: int, dt_s: float
    ) -> np.ndarray:
        p = self.parameters
        node_capacity = p.total_heat_capacity_J_K / 5.0
        jac_rhs = np.zeros((5, 5), dtype=float)
        # Heat-generation derivative for the thermoneutral/polarization proxy.
        if self.heat_model == "A_thermoneutral_proxy":
            mean_temperature = float(np.mean(state))
            dq_d_each_state = 0.004 * max(float(current_A), 0.0) * (mean_temperature - 55.0) / 5.0
            # Each node receives one fifth of total generation.
            jac_rhs += dq_d_each_state / 5.0
        coolant = p.coolant_forward_C if direction == 1 else p.coolant_forward_C[::-1]
        active = state > coolant
        jac_rhs[np.arange(5), np.arange(5)] -= (
            p.cooling_coefficient_per_node_W_K * float(np.clip(duty, 0.0, 1.0)) * active
        )
        k = p.conduction_edge_W_K
        for edge in range(4):
            jac_rhs[edge, edge] -= k
            jac_rhs[edge, edge + 1] += k
            jac_rhs[edge + 1, edge] += k
            jac_rhs[edge + 1, edge + 1] -= k
        return np.eye(5) + dt_s * jac_rhs / node_capacity

    def predict(self, current_A: float, duty: float, direction: int, dt_s: float) -> None:
        prior = self.state_C.copy()
        prediction = self._transition(prior, current_A, duty, direction, dt_s)
        jacobian = self._analytic_jacobian(prior, current_A, duty, direction, dt_s)
        self.state_C = prediction
        self.covariance = jacobian @ self.covariance @ jacobian.T + self.process_covariance
