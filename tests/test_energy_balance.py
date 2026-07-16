import numpy as np

from test_three_node import run


def test_discrete_energy_balance() -> None:
    result = run(0.5)
    scale = np.max(np.sum(result.node_heat_generation_W, axis=1))
    relative = np.max(np.abs(result.total_energy_residual_W)) / scale
    assert relative < 1e-12
