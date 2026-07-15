"""Likelihood-weighted bank of five-node EKFs."""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

from observers.ekf import ExtendedKalmanFilter


@dataclass
class MultipleModelObserver:
    filters: list[ExtendedKalmanFilter]
    probabilities: np.ndarray
    model_names: tuple[str, ...] | None = None
    input_scales: np.ndarray | None = None
    probability_floor: float = 1e-4

    def __post_init__(self) -> None:
        count = len(self.filters)
        if self.probabilities.shape != (count,):
            raise ValueError("probabilities must have one entry per filter")
        if self.model_names is None:
            self.model_names = tuple(f"model_{index}" for index in range(count))
        if len(self.model_names) != count:
            raise ValueError("model_names must have one entry per filter")
        if self.input_scales is None:
            self.input_scales = np.ones(count)
        self.input_scales = np.asarray(self.input_scales, dtype=float)
        if self.input_scales.shape != (count,):
            raise ValueError("input_scales must have one entry per filter")

    def predict(self, current_A: float, duty: float, direction: int, dt_s: float) -> None:
        assert self.input_scales is not None
        for observer, scale in zip(self.filters, self.input_scales, strict=True):
            observer.predict(current_A * float(scale), duty, direction, dt_s)

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
