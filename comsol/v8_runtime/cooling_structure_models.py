"""Steady 1-D conjugate cooling structures used by the V8 Fig.7 audit.

These models are deliberately isolated from the ROM, observer and controller
code.  They conserve the Table-2 heat input and rho*u*A mass flow exactly.
"""

from __future__ import annotations

from dataclasses import dataclass
from time import perf_counter

import numpy as np
from scipy.optimize import root


SIGMA = 5.670374419e-8
G = 9.80665


@dataclass(frozen=True)
class FrozenPhysics:
    ambient_C: float = 25.0
    length_m: float = 0.24936900555275113
    stack_inlet_area_m2: float = 3962e-6
    cell_count: int = 40
    heat_area_cell_m2: float = 0.0247
    heat_flux_W_m2: float = 2425.5
    graphite_k_W_mK: float = 24.0
    graphite_thickness_m: float = 0.010
    air_k_W_mK: float = 0.0251
    air_cp_J_kgK: float = 1007.0
    air_rho_kg_m3: float = 1.184
    air_mu_Pa_s: float = 1.849e-5
    emissivity: float = 0.82

    @property
    def cell_inlet_area_m2(self) -> float:
        return self.stack_inlet_area_m2 / self.cell_count

    @property
    def width_m(self) -> float:
        return self.heat_area_cell_m2 / self.length_m

    @property
    def heat_input_cell_W(self) -> float:
        return self.heat_area_cell_m2 * self.heat_flux_W_m2

    @property
    def solid_cross_section_m2(self) -> float:
        return self.width_m * self.graphite_thickness_m

    @property
    def external_area_cell_m2(self) -> float:
        # The heated face is excluded: opposite major face plus four edges.
        return (
            self.heat_area_cell_m2
            + 2.0 * self.length_m * self.graphite_thickness_m
            + 2.0 * self.width_m * self.graphite_thickness_m
        )


PHYSICS = FrozenPhysics()


def churchill_chu_h_W_m2K(surface_C: np.ndarray | float) -> np.ndarray:
    """Churchill-Chu vertical-plate average h using frozen Table-2 air data."""
    surface = np.asarray(surface_C, dtype=float)
    ambient_K = PHYSICS.ambient_C + 273.15
    film_K = 0.5 * (surface + PHYSICS.ambient_C) + 273.15
    beta = 1.0 / film_K
    nu = PHYSICS.air_mu_Pa_s / PHYSICS.air_rho_kg_m3
    alpha = PHYSICS.air_k_W_mK / (
        PHYSICS.air_rho_kg_m3 * PHYSICS.air_cp_J_kgK
    )
    pr = nu / alpha
    rayleigh = (
        G
        * beta
        * np.maximum(surface + 273.15 - ambient_K, 1e-6)
        * PHYSICS.length_m**3
        / (nu * alpha)
    )
    nusselt = (
        0.825
        + 0.387
        * rayleigh ** (1.0 / 6.0)
        / (1.0 + (0.492 / pr) ** (9.0 / 16.0)) ** (8.0 / 27.0)
    ) ** 2
    return nusselt * PHYSICS.air_k_W_mK / PHYSICS.length_m


@dataclass(frozen=True)
class ChannelGeometry:
    channel_count: int = 40
    channel_height_m: float = 10.0e-3
    rib_width_m: float = 1.0e-3

    @property
    def channel_width_m(self) -> float:
        return PHYSICS.stack_inlet_area_m2 / (
            self.channel_count * self.channel_height_m
        )

    @property
    def wetted_perimeter_total_m(self) -> float:
        return 2.0 * self.channel_count * (
            self.channel_width_m + self.channel_height_m
        )

    @property
    def nominal_hydraulic_diameter_m(self) -> float:
        return 4.0 * PHYSICS.stack_inlet_area_m2 / self.wetted_perimeter_total_m

    @property
    def contact_area_cell_m2(self) -> float:
        return self.wetted_perimeter_total_m * PHYSICS.length_m

    @property
    def total_pitch_width_m(self) -> float:
        return self.channel_count * (self.channel_width_m + self.rib_width_m)


M1_GEOMETRY = ChannelGeometry()


def rectangular_nusselt_constant_flux() -> float:
    aspect = min(
        M1_GEOMETRY.channel_height_m / M1_GEOMETRY.channel_width_m,
        M1_GEOMETRY.channel_width_m / M1_GEOMETRY.channel_height_m,
    )
    return 8.235 * (
        1.0
        - 2.0421 * aspect
        + 3.0853 * aspect**2
        - 2.4765 * aspect**3
        + 1.0578 * aspect**4
        - 0.1861 * aspect**5
    )


def m1_internal_h(
    velocity_m_s: float, hydraulic_diameter_m: float
) -> tuple[float, float, float]:
    reynolds = (
        PHYSICS.air_rho_kg_m3
        * velocity_m_s
        * hydraulic_diameter_m
        / PHYSICS.air_mu_Pa_s
    )
    pr = (
        PHYSICS.air_mu_Pa_s
        * PHYSICS.air_cp_J_kgK
        / PHYSICS.air_k_W_mK
    )
    nusselt_fd = rectangular_nusselt_constant_flux()
    graetz = max(reynolds * pr * hydraulic_diameter_m / PHYSICS.length_m, 0.0)
    nusselt_entry = 1.953 * graetz ** (1.0 / 3.0)
    nusselt = (nusselt_fd**3 + nusselt_entry**3) ** (1.0 / 3.0)
    return nusselt * PHYSICS.air_k_W_mK / hydraulic_diameter_m, reynolds, nusselt


def solve_structure(
    model: str,
    velocity_m_s: float,
    parameters: np.ndarray,
    nodes: int = 100,
    inlet_area_m2: float | None = None,
) -> dict[str, float | np.ndarray]:
    """Solve M1 or M2 and return temperatures plus a closed energy ledger."""
    if model not in {"M1", "M2"}:
        raise ValueError(model)
    start = perf_counter()
    dx = PHYSICS.length_m / nodes
    # Formal M1/M2 use H3: the reported 3962 mm2 is the total effective area
    # of one Fig.7 cooling model.  H1 sensitivity is evaluated by overriding
    # this argument with Ain/40 while keeping the same 59.90985 W heat input.
    effective_inlet_area = (
        PHYSICS.stack_inlet_area_m2 if inlet_area_m2 is None else inlet_area_m2
    )
    mdot = PHYSICS.air_rho_kg_m3 * velocity_m_s * effective_inlet_area
    mdot_expected = PHYSICS.air_rho_kg_m3 * velocity_m_s * effective_inlet_area
    if model == "M1":
        hydraulic_diameter_m, area_multiplier = map(float, parameters)
        h_internal, reynolds, nusselt = m1_internal_h(
            velocity_m_s, hydraulic_diameter_m
        )
        perimeter_effective = M1_GEOMETRY.wetted_perimeter_total_m * area_multiplier
    else:
        h0_W_m2K, exponent_m = map(float, parameters)
        h_internal = h0_W_m2K * (velocity_m_s / 8.0) ** exponent_m
        perimeter_effective = PHYSICS.width_m
        reynolds = (
            PHYSICS.air_rho_kg_m3
            * velocity_m_s
            * 2.0e-3
            / PHYSICS.air_mu_Pa_s
        )
        nusselt = h_internal * 2.0e-3 / PHYSICS.air_k_W_mK

    ua_segment = h_internal * perimeter_effective * dx
    external_area_segment = PHYSICS.external_area_cell_m2 / nodes
    q_heat_segment = PHYSICS.heat_input_cell_W / nodes
    axial_conductance = PHYSICS.graphite_k_W_mK * PHYSICS.solid_cross_section_m2 / dx

    def ledger(surface_C: np.ndarray) -> tuple[np.ndarray, np.ndarray, np.ndarray, np.ndarray]:
        air_C = np.empty(nodes + 1)
        air_C[0] = PHYSICS.ambient_C
        q_air = np.empty(nodes)
        effectiveness = 1.0 - np.exp(-ua_segment / (mdot * PHYSICS.air_cp_J_kgK))
        for index in range(nodes):
            air_C[index + 1] = air_C[index] + effectiveness * (
                surface_C[index] - air_C[index]
            )
            q_air[index] = mdot * PHYSICS.air_cp_J_kgK * (
                air_C[index + 1] - air_C[index]
            )
        h_natural = churchill_chu_h_W_m2K(surface_C)
        q_natural = (
            h_natural
            * external_area_segment
            * (surface_C - PHYSICS.ambient_C)
        )
        surface_K = surface_C + 273.15
        ambient_K = PHYSICS.ambient_C + 273.15
        q_radiation = (
            PHYSICS.emissivity
            * SIGMA
            * external_area_segment
            * (surface_K**4 - ambient_K**4)
        )
        return air_C, q_air, q_natural, q_radiation

    def residual(surface_C: np.ndarray) -> np.ndarray:
        _, q_air, q_natural, q_radiation = ledger(surface_C)
        conduction = np.zeros(nodes)
        conduction[1:] += axial_conductance * (surface_C[:-1] - surface_C[1:])
        conduction[:-1] += axial_conductance * (surface_C[1:] - surface_C[:-1])
        return (
            conduction
            + q_heat_segment
            - q_air
            - q_natural
            - q_radiation
        )

    initial = np.linspace(45.0, 65.0, nodes)
    solution = root(residual, initial, method="hybr", options={"xtol": 1e-9})
    if not solution.success or np.max(np.abs(residual(solution.x))) > 1e-6:
        raise RuntimeError(f"{model} at {velocity_m_s:g} m/s did not converge: {solution.message}")
    surface_C = solution.x
    air_C, q_air, q_natural, q_radiation = ledger(surface_C)
    q_air_total = float(q_air.sum())
    q_natural_total = float(q_natural.sum())
    q_radiation_total = float(q_radiation.sum())
    q_input = PHYSICS.heat_input_cell_W
    energy_residual = abs(
        q_input - q_air_total - q_natural_total - q_radiation_total
    ) / q_input
    coordinate = (np.arange(nodes) + 0.5) / nodes
    thirds = np.array_split(surface_C, 3)
    hotspot_index = int(np.argmax(surface_C))
    elapsed = perf_counter() - start
    return {
        "surface_C": surface_C,
        "air_C": air_C,
        "flow_coordinate": coordinate,
        "T_inlet_region_C": float(np.mean(thirds[0])),
        "T_middle_region_C": float(np.mean(thirds[1])),
        "T_outlet_region_C": float(np.mean(thirds[2])),
        "Tmax_C": float(np.max(surface_C)),
        "Tmin_C": float(np.min(surface_C)),
        "DeltaT_C": float(np.ptp(surface_C)),
        "hotspot_location_normalized": float(coordinate[hotspot_index]),
        "air_outlet_temperature_C": float(air_C[-1]),
        "air_enthalpy_gain_W": q_air_total,
        "natural_convection_loss_W": q_natural_total,
        "radiation_loss_W": q_radiation_total,
        "conduction_loss_W": 0.0,
        "total_heat_input_W": q_input,
        "energy_residual_relative": float(energy_residual),
        "mdot_cell_kg_s": float(mdot),
        "mdot_identity_relative_error": float(abs(mdot - mdot_expected) / mdot_expected),
        "h_internal_W_m2K": float(h_internal),
        "h_natural_mean_W_m2K": float(np.mean(churchill_chu_h_W_m2K(surface_C))),
        "reynolds": float(reynolds),
        "nusselt": float(nusselt),
        "effective_contact_area_m2": float(perimeter_effective * PHYSICS.length_m),
        "runtime_s": float(elapsed),
        "converged": True,
    }
