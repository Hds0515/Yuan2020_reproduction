from __future__ import annotations

import numpy as np

from models_v10.shahsavari_channel_reference import (
    TABLE1_CASES,
    LiteratureGeometry,
    shah_london_constant_heat_flux_nusselt,
    solve_channel,
)


def test_reported_trapezoid_geometry_is_positive_and_auditable() -> None:
    geometry = LiteratureGeometry()
    assert geometry.channel_count == 80
    assert geometry.channel_short_width_m > 0.0
    assert geometry.channel_area_m2 > 0.0
    assert np.isclose(
        geometry.hydraulic_diameter_m,
        4.0 * geometry.channel_area_m2 / geometry.channel_wetted_perimeter_m,
    )


def test_square_duct_nusselt_is_standard_value() -> None:
    assert np.isclose(shah_london_constant_heat_flux_nusselt(1.0), 3.610224, atol=1e-6)


def test_channel_energy_is_conserved() -> None:
    result = solve_channel(TABLE1_CASES[2], nodes=80)
    assert abs(result.energy_residual_fraction) < 1e-10
    assert np.isclose(result.air_enthalpy_gain_channel_W, result.channel_heat_input_W)


def test_temperature_and_hotspot_are_physically_ordered() -> None:
    result = solve_channel(TABLE1_CASES[3], nodes=80)
    assert result.Tmax_C > result.Tmin_C > TABLE1_CASES[3].inlet_temperature_C
    assert result.air_outlet_temperature_C > TABLE1_CASES[3].inlet_temperature_C
    assert result.hotspot_position_normalized > 0.9
    assert result.reynolds_number < 800.0
