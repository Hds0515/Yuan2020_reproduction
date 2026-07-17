"""Conservative representative-channel model based on Shahsavari et al. (2012).

This is a literature-constrained *reference channel*, not a reconstruction of
Yuan et al.'s undisclosed cooling CAD.  No parameter in this module is fitted
to the six validation temperatures.  The supplied cell heat loads and inlet
velocities are the experimental inputs reported in Shahsavari et al., Table 1.

The solid is represented by an axially conducting laminate.  Cathode air gains
exactly the heat removed by a fully developed laminar forced-convection term.
The trapezoidal channel is reduced to its reported area and wetted perimeter;
the square-duct Shah--London Nusselt number is an explicit closure assumption.
"""

from __future__ import annotations

from dataclasses import asdict, dataclass
from time import perf_counter

import numpy as np


@dataclass(frozen=True)
class LiteratureGeometry:
    """Dimensions and solid properties from Shahsavari et al., Table 2."""

    flow_length_m: float = 0.060
    cell_width_m: float = 0.280
    channel_count: int = 80
    cathode_channel_height_m: float = 0.0025
    cathode_channel_long_width_m: float = 0.0025
    cathode_wall_angle_deg: float = 80.0
    bipolar_plate_thickness_m: float = 0.0055
    gdl_anode_thickness_m: float = 0.0002
    gdl_cathode_thickness_m: float = 0.0002
    ccm_thickness_m: float = 0.00005
    bipolar_plate_in_plane_k_W_mK: float = 60.0
    gdl_in_plane_k_W_mK: float = 10.0
    ccm_k_W_mK: float = 1.5

    @property
    def channel_pitch_m(self) -> float:
        return self.cell_width_m / self.channel_count

    @property
    def channel_short_width_m(self) -> float:
        side_reduction = 2.0 * self.cathode_channel_height_m / np.tan(
            np.deg2rad(self.cathode_wall_angle_deg)
        )
        return self.cathode_channel_long_width_m - side_reduction

    @property
    def channel_area_m2(self) -> float:
        return (
            0.5
            * (self.cathode_channel_long_width_m + self.channel_short_width_m)
            * self.cathode_channel_height_m
        )

    @property
    def channel_wetted_perimeter_m(self) -> float:
        side_length = self.cathode_channel_height_m / np.sin(
            np.deg2rad(self.cathode_wall_angle_deg)
        )
        return (
            self.cathode_channel_long_width_m
            + self.channel_short_width_m
            + 2.0 * side_length
        )

    @property
    def hydraulic_diameter_m(self) -> float:
        return 4.0 * self.channel_area_m2 / self.channel_wetted_perimeter_m

    @property
    def axial_kA_W_m_K(self) -> float:
        pitch = self.channel_pitch_m
        return pitch * (
            self.bipolar_plate_in_plane_k_W_mK * self.bipolar_plate_thickness_m
            + self.gdl_in_plane_k_W_mK
            * (self.gdl_anode_thickness_m + self.gdl_cathode_thickness_m)
            + self.ccm_k_W_mK * self.ccm_thickness_m
        )


@dataclass(frozen=True)
class AirProperties:
    """Dry-air constants near the 21 degC experiment inlet condition."""

    density_kg_m3: float = 1.204
    heat_capacity_J_kgK: float = 1006.0
    thermal_conductivity_W_mK: float = 0.02514
    dynamic_viscosity_Pa_s: float = 1.825e-5


@dataclass(frozen=True)
class ExperimentalCase:
    case_id: int
    inlet_velocity_m_s: float
    total_cell_heat_W: float
    experimental_Tmax_C: float
    paper_simulated_Tmax_C: float
    inlet_temperature_C: float = 21.0


TABLE1_CASES: tuple[ExperimentalCase, ...] = (
    ExperimentalCase(1, 1.76, 55.7, 76.0, 77.0),
    ExperimentalCase(2, 2.19, 37.7, 65.0, 68.0),
    ExperimentalCase(3, 1.63, 27.6, 61.0, 64.0),
    ExperimentalCase(4, 1.45, 29.1, 67.0, 69.0),
    ExperimentalCase(5, 0.95, 6.3, 35.0, 36.0),
    ExperimentalCase(6, 0.94, 13.7, 50.0, 54.0),
)


@dataclass
class ChannelSolution:
    case_id: int
    nodes: int
    x_m: np.ndarray
    solid_temperature_C: np.ndarray
    air_temperature_C: np.ndarray
    inlet_velocity_m_s: float
    total_cell_heat_W: float
    experimental_Tmax_C: float
    paper_simulated_Tmax_C: float
    reynolds_number: float
    nusselt_number: float
    forced_h_W_m2K: float
    mass_flow_channel_kg_s: float
    air_enthalpy_gain_channel_W: float
    channel_heat_input_W: float
    energy_residual_fraction: float
    pressure_drop_Pa: float
    solve_time_s: float

    @property
    def Tmax_C(self) -> float:
        return float(np.max(self.solid_temperature_C))

    @property
    def Tmin_C(self) -> float:
        return float(np.min(self.solid_temperature_C))

    @property
    def DeltaT_C(self) -> float:
        return self.Tmax_C - self.Tmin_C

    @property
    def hotspot_position_normalized(self) -> float:
        index = int(np.argmax(self.solid_temperature_C))
        return float(self.x_m[index] / self.x_m[-1])

    @property
    def air_outlet_temperature_C(self) -> float:
        return float(self.air_temperature_C[-1])

    def region_mean_C(self, start: float, stop: float) -> float:
        xi = self.x_m / self.x_m[-1]
        mask = (xi >= start) & (xi <= stop)
        return float(np.mean(self.solid_temperature_C[mask]))

    def summary(self) -> dict[str, float | int]:
        return {
            "case_id": self.case_id,
            "nodes": self.nodes,
            "inlet_velocity_m_s": self.inlet_velocity_m_s,
            "total_cell_heat_W": self.total_cell_heat_W,
            "experimental_Tmax_C": self.experimental_Tmax_C,
            "paper_simulated_Tmax_C": self.paper_simulated_Tmax_C,
            "T_inlet_region_C": self.region_mean_C(0.0, 1.0 / 3.0),
            "T_middle_region_C": self.region_mean_C(1.0 / 3.0, 2.0 / 3.0),
            "T_outlet_region_C": self.region_mean_C(2.0 / 3.0, 1.0),
            "Tmax_C": self.Tmax_C,
            "Tmin_C": self.Tmin_C,
            "DeltaT_C": self.DeltaT_C,
            "hotspot_position_normalized": self.hotspot_position_normalized,
            "air_outlet_temperature_C": self.air_outlet_temperature_C,
            "air_enthalpy_gain_channel_W": self.air_enthalpy_gain_channel_W,
            "channel_heat_input_W": self.channel_heat_input_W,
            "energy_residual_fraction": self.energy_residual_fraction,
            "reynolds_number": self.reynolds_number,
            "nusselt_number": self.nusselt_number,
            "forced_h_W_m2K": self.forced_h_W_m2K,
            "mass_flow_channel_kg_s": self.mass_flow_channel_kg_s,
            "pressure_drop_Pa": self.pressure_drop_Pa,
            "solve_time_s": self.solve_time_s,
        }


def shah_london_constant_heat_flux_nusselt(aspect_ratio: float = 1.0) -> float:
    """Fully developed rectangular-duct Nu for uniform heat flux.

    The polynomial expects 0 < aspect_ratio <= 1.  A value of one is retained
    as the declared equivalent-rectangle approximation because the reported
    cathode long width and height are both 2.5 mm.
    """

    alpha = float(aspect_ratio)
    if not 0.0 < alpha <= 1.0:
        raise ValueError("aspect_ratio must be in (0, 1]")
    return 8.235 * (
        1.0
        - 2.0421 * alpha
        + 3.0853 * alpha**2
        - 2.4765 * alpha**3
        + 1.0578 * alpha**4
        - 0.1861 * alpha**5
    )


def solve_channel(
    case: ExperimentalCase,
    *,
    nodes: int = 200,
    geometry: LiteratureGeometry | None = None,
    air: AirProperties | None = None,
) -> ChannelSolution:
    """Solve the conservative 1-D solid/air channel equations."""

    if nodes < 10:
        raise ValueError("nodes must be at least 10")
    geometry = geometry or LiteratureGeometry()
    air = air or AirProperties()
    start = perf_counter()

    length = geometry.flow_length_m
    dx = length / nodes
    x_m = (np.arange(nodes, dtype=float) + 0.5) * dx
    nusselt = shah_london_constant_heat_flux_nusselt(1.0)
    forced_h = nusselt * air.thermal_conductivity_W_mK / geometry.hydraulic_diameter_m
    mdot = air.density_kg_m3 * case.inlet_velocity_m_s * geometry.channel_area_m2
    mcp = mdot * air.heat_capacity_J_kgK
    heat_channel = case.total_cell_heat_W / geometry.channel_count
    q_node = heat_channel / nodes
    convectance_node = forced_h * geometry.channel_wetted_perimeter_m * dx
    axial_conductance = geometry.axial_kA_W_m_K / dx

    # Ta[i] = c[i] + sum_j D[i,j] * Ts[j].  The upwind energy equation is
    # exact for the discrete control volumes, so convective heat removed from
    # the solid is identical to air enthalpy gain to round-off.
    transfer = convectance_node / mcp
    c = np.empty(nodes, dtype=float)
    D = np.zeros((nodes, nodes), dtype=float)
    inlet = case.inlet_temperature_C
    c_prev = inlet
    d_prev = np.zeros(nodes, dtype=float)
    for index in range(nodes):
        c[index] = c_prev
        D[index, :] = d_prev
        c_prev = c_prev * (1.0 - transfer)
        d_prev = d_prev * (1.0 - transfer)
        d_prev[index] += transfer

    matrix = convectance_node * (np.eye(nodes) - D)
    rhs = np.full(nodes, q_node, dtype=float) + convectance_node * c
    for index in range(nodes - 1):
        matrix[index, index] += axial_conductance
        matrix[index, index + 1] -= axial_conductance
        matrix[index + 1, index + 1] += axial_conductance
        matrix[index + 1, index] -= axial_conductance

    solid_C = np.linalg.solve(matrix, rhs)
    air_inlet_to_cells = c + D @ solid_C
    air_outlet = c_prev + d_prev @ solid_C
    air_C = np.concatenate([air_inlet_to_cells[1:], [air_outlet]])
    air_gain = mcp * (air_outlet - inlet)
    residual = (heat_channel - air_gain) / heat_channel

    reynolds = (
        air.density_kg_m3
        * case.inlet_velocity_m_s
        * geometry.hydraulic_diameter_m
        / air.dynamic_viscosity_Pa_s
    )
    # Square-duct approximation: Darcy f*Re = 56.91.
    darcy_friction = 56.91 / reynolds
    pressure_drop = (
        darcy_friction
        * length
        / geometry.hydraulic_diameter_m
        * 0.5
        * air.density_kg_m3
        * case.inlet_velocity_m_s**2
    )

    return ChannelSolution(
        case_id=case.case_id,
        nodes=nodes,
        x_m=x_m,
        solid_temperature_C=solid_C,
        air_temperature_C=air_C,
        inlet_velocity_m_s=case.inlet_velocity_m_s,
        total_cell_heat_W=case.total_cell_heat_W,
        experimental_Tmax_C=case.experimental_Tmax_C,
        paper_simulated_Tmax_C=case.paper_simulated_Tmax_C,
        reynolds_number=float(reynolds),
        nusselt_number=float(nusselt),
        forced_h_W_m2K=float(forced_h),
        mass_flow_channel_kg_s=float(mdot),
        air_enthalpy_gain_channel_W=float(air_gain),
        channel_heat_input_W=float(heat_channel),
        energy_residual_fraction=float(residual),
        pressure_drop_Pa=float(pressure_drop),
        solve_time_s=float(perf_counter() - start),
    )


def frozen_parameter_record() -> dict[str, object]:
    """Return a serializable record of every model constant and assumption."""

    geometry = LiteratureGeometry()
    air = AirProperties()
    return {
        "geometry": asdict(geometry),
        "derived_geometry": {
            "channel_pitch_m": geometry.channel_pitch_m,
            "channel_short_width_m": geometry.channel_short_width_m,
            "channel_area_m2": geometry.channel_area_m2,
            "channel_wetted_perimeter_m": geometry.channel_wetted_perimeter_m,
            "hydraulic_diameter_m": geometry.hydraulic_diameter_m,
            "axial_kA_W_m_K": geometry.axial_kA_W_m_K,
        },
        "air_properties": asdict(air),
        "heat_transfer_closure": {
            "correlation": "Shah-London fully developed rectangular duct, uniform heat flux",
            "declared_equivalent_aspect_ratio": 1.0,
            "nusselt_number": shah_london_constant_heat_flux_nusselt(1.0),
            "fitted_to_validation_data": False,
        },
        "scope": "representative cathode channel; not a recovered full-cell manifold/CAD",
    }
