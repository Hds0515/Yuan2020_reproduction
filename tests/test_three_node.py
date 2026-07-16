import numpy as np

from models.three_node_model import ThermalParameters, simulate


PARAMETERS = ThermalParameters(4381.600240346621, 31.722471972119184, 7.1833293414246695, 0.0727170185330114, 0.0014992061134433243)


def run(dt: float, duration: float = 100.0):
    count = int(duration / dt) + 1
    return simulate(
        PARAMETERS,
        np.full(count, 10.0),
        np.array([53.0, 53.2, 53.8]),
        dt_s=dt,
        bidirectional=False,
        initial_airflow_direction=1,
        fan_initially_enabled=True,
        fan_enable_temperature_C=53.0,
        reference_temperature_C=55.0,
        smc_deadband_C=0.5,
        smc_update_period_s=10.0,
        heat_model="A_thermoneutral_proxy",
        coolant_node_temperatures_forward_C=np.array([24.8, 30.05, 35.3]),
    )


def test_dt_is_explicit_and_convergent() -> None:
    coarse = run(1.0)
    fine = run(0.5)
    assert np.max(np.abs(coarse.temperature_C[-1] - fine.temperature_C[-1])) < 0.08


def test_outputs_have_requested_shapes() -> None:
    result = run(1.0)
    assert result.temperature_C.shape == (101, 3)
    assert result.node_heat_generation_W.shape == (101, 3)
    assert result.node_cooling_heat_W.shape == (101, 3)
    assert result.edge_conduction_W.shape == (101, 2)
    assert set(np.unique(result.airflow_direction)).issubset({-1, 1})
