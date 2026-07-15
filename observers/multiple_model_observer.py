"""Likelihood-weighted bank of five-node EKFs."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from observers.ekf import ExtendedKalmanFilter


@dataclass
class MultipleModelObserver:
    filters: list[ExtendedKalmanFilter]
    probabilities: np.ndarray
    probability_floor: float = 1e-4

    def predict(self, current_A: float, duty: float, direction: int, dt_s: float) -> None:
        for observer in self.filters:
            observer.predict(current_A, duty, direction, dt_s)

    def update(self, measurement_C: np.ndarray) -> None:
        log_likelihood = np.asarray([observer.update(measurement_C) for observer in self.filters])
        log_weight = np.log(np.maximum(self.probabilities, self.probability_floor)) + log_likelihood
        log_weight -= np.max(log_weight)
        weight = np.exp(log_weight)
        self.probabilities = np.maximum(weight / np.sum(weight), self.probability_floor)
        self.probabilities /= np.sum(self.probabilities)

    @property
    def state_C(self) -> np.ndarray:
        return np.sum([weight * observer.state_C for weight, observer in zip(self.probabilities, self.filters, strict=True)], axis=0)
