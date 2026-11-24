"""One simulator for every energy route in docs/energy.md.

Every route turns a few parameters into the same energy account:

    electricity --(driver efficiency)--> driver energy (laser light, accelerator beam, pulsed current)
                --(coupling)--> energy in the fuel --(the route's physics model)--> fusion energy
                --(blanket)--> heat --(thermal efficiency)--> electricity out

"electric gain" = electricity out / electricity in; net energy means > 1.

Every design also lists what it demands of technology (laser energy and
power, magnetic field, density, current, pinch stability, muon supply ...),
each compared with the best achieved so far. Its "stretch" is the largest of
those ratios: 1 means today's best would do, 10 means some demand is ten times
beyond it. Unproven physics (spin-polarised fuel surviving in a plasma) counts
as infinite stretch. The search finds the best electric gain of each route
within a stretch budget.

    python -m nanofusion sim list
    python -m nanofusion sim run tamped_ball radius_mm=0.3 t_i_kev=20 b_tesla=3000
    python -m nanofusion sim search --samples 300

The physics is the project's own models (energy.py, confinement.py,
alternatives.py, degenerate.py, walls.py); this module only wires them into one
account. Sections refer to docs/energy.md.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from dataclasses import dataclass, field, replace
from functools import lru_cache
from typing import Any, Callable

import numpy as np

from . import alternatives as ALT
from . import confinement as C
from . import degenerate as D
from . import energy as EN
from . import physics as P
from . import walls as W

KEV = P.J_PER_KEV
MEV = 1e3 * KEV
Q_DT_MEV = P.Q_VALUE_MEV["DT"]
MU0 = 4e-7 * math.pi
INF = math.inf


def _num(x: float) -> str:
    """Two significant digits without an exponent below 10^5: 0.42, 31, 310."""
    return f"{float(f'{x:.2g}'):g}" if abs(x) < 1e5 else f"{x:.1e}"


def deplete(raw: float, most: float) -> float:
    """Fusion energy once the fuel burns down: a fuel burning at a fixed rate would release `raw`, but at most
    `most` is there to burn. For two equal reactant densities n' = -k n^2, the burnt share is x / (1 + x)."""
    if most <= 0:
        return 0.0
    x = raw / most
    return most * x / (1 + x)


@dataclass
class Settings:
    """Efficiencies shared by every route."""

    laser_efficiency: float = 0.10  # electricity -> light, fusion-laser designs (NIF itself: ~0.5%)
    pulsed_power_efficiency: float = 0.15  # electricity -> pinch, about what Z delivers to its load
    accelerator_efficiency: float = 0.40  # electricity -> beam, accelerator-driven reactor designs (PSI: 18%)
    thermal_efficiency: float = 0.40  # heat -> electricity
    laser_to_ions: float = 0.15  # light -> crystal explosion energy (or fuel heat), best case (sections 7, 12)
    implosion_efficiency: float = 0.01  # light -> energy of squeezed fuel, NIF-class implosions
    muon_cost_gev: float = 4.7  # beam energy per muon, optimised accelerator design (section 9)
    blanket: str = "lithium"  # "none", "lithium" (breeds tritium: +4.78 MeV per D-T fusion) or "uranium"
    k_eff: float = 0.0  # uranium blanket multiplication (subcritical), used when blanket == "uranium"

    def efficiency(self, driver: str) -> float:
        return {"laser": self.laser_efficiency, "pulsed power": self.pulsed_power_efficiency,
                "accelerator": self.accelerator_efficiency}[driver]


# Best achieved so far, for each demand (value, what achieved it).
TODAY = {
    "laser energy per shot (J)": (2.08e6, "NIF"),
    "picosecond laser energy (J)": (4.4e3, "NIF ARC, four 30 ps beamlets"),
    "laser peak power (W)": (1.02e16, "ELI-NP, 10 PW"),
    "laser intensity (W/cm²)": (1.1e23, "2021 record"),
    "magnetic field on the target (T)": (800.0, "laser-driven coil, claimed"),
    "cold squeezed density (g/cm³)": (1000.0, "NIF design density"),
    "radial squeeze of the fuel": (35.0, "NIF-class implosions, 30–40"),
    "hot-spot areal density (g/cm²)": (0.44, "NIF's best shots"),
    "hot-spot pressure (Gbar)": (600.0, "NIF's best shots"),
    "pinch current (A)": (26e6, "Sandia's Z"),
    "energy into the pinch (J)": (3e6, "Z (2.7 MJ of X-rays from ~20 MJ stored)"),
    "stable Alfvén times of a dense hot pinch": (1.0, "fibre pinches broke up at peak current"),
    "muons per second": (1e8, "PSI, J-PARC"),
    "liquid or solid D-T density (× liquid hydrogen)": (1.45, "solid D-T"),
    "projectile speed (m/s)": (8e4, "dust accelerators"),
    "crystals' distance from the core (m)": (100e-6, "closer than ~0.1 mm is swept into the implosion"),
    "blanket k_eff": (0.98, "accelerator-driven reactor designs"),
}
LOWER_IS_HARDER = {"crystals' distance from the core (m)"}
UNPROVEN = INF


@dataclass
class Requirement:
    name: str
    needed: float

    @property
    def best(self) -> float:
        return TODAY[self.name][0] if self.name in TODAY else math.nan

    @property
    def stretch(self) -> float:
        if self.name not in TODAY:
            return UNPROVEN
        if self.name == "blanket k_eff":
            return max(1.0, (1 - self.best) / max(1 - self.needed, 1e-12))
        if self.name in LOWER_IS_HARDER:
            return max(1.0, self.best / max(self.needed, 1e-300))
        return max(1.0, self.needed / self.best)


@dataclass
class Outcome:
    route: str
    params: dict
    fusion_J: float
    fuel_J: float
    drivers: dict  # driver kind -> energy it must deliver (J)
    dt_fusions: float = 0.0  # D-T fusions, whose neutrons can feed a blanket
    requirements: list = field(default_factory=list)
    notes: str = ""
    settings: Settings = field(default_factory=Settings)

    @property
    def driver_J(self) -> float:
        return sum(self.drivers.values())

    @property
    def electricity_in_J(self) -> float:
        return sum(e / self.settings.efficiency(k) for k, e in self.drivers.items())

    @property
    def blanket_extra_mev(self) -> float:
        s = self.settings
        if s.blanket == "none":
            return 0.0
        return ALT.blanket_energy_per_fusion_mev(s.k_eff if s.blanket == "uranium" else 0.0,
                                                  uranium=s.blanket == "uranium") - Q_DT_MEV

    @property
    def heat_J(self) -> float:
        return self.fusion_J + self.dt_fusions * self.blanket_extra_mev * MEV

    @property
    def target_gain(self) -> float:
        return self.fusion_J / self.driver_J if self.driver_J > 0 else 0.0

    @property
    def electric_gain(self) -> float:
        e_in = self.electricity_in_J
        return self.heat_J * self.settings.thermal_efficiency / e_in if e_in > 0 else 0.0

    @property
    def all_requirements(self) -> list:
        reqs = list(self.requirements)
        if self.settings.blanket == "uranium" and self.dt_fusions > 0:
            reqs.append(Requirement("blanket k_eff", self.settings.k_eff))
        if self.params.get("spin"):
            reqs.append(Requirement("spin-polarised D-T that survives in the fuel", 1.0))
        return reqs

    @property
    def stretch(self) -> float:
        return max([r.stretch for r in self.all_requirements], default=1.0)

    @property
    def limiting(self) -> str:
        reqs = self.all_requirements
        if not reqs:
            return ""
        worst = max(reqs, key=lambda r: r.stretch)
        return worst.name if worst.stretch > 1 else ""

    def summary(self) -> dict:
        return {"route": self.route, "params": self.params, "fusion_J": self.fusion_J, "fuel_J": self.fuel_J,
                "drivers_J": self.drivers, "electricity_in_J": self.electricity_in_J, "heat_J": self.heat_J,
                "target_gain": self.target_gain, "electric_gain": self.electric_gain, "stretch": self.stretch,
                "limiting": self.limiting, "notes": self.notes,
                "requirements": [{"name": r.name, "needed": r.needed, "best": r.best, "stretch": r.stretch}
                                 for r in self.all_requirements]}


@dataclass
class Param:
    name: str
    low: float | None = None
    high: float | None = None
    choices: tuple | None = None
    log: bool = True
    default: Any = None
    help: str = ""

    def sample(self, rng: np.random.Generator) -> Any:
        if self.choices is not None:
            return self.choices[int(rng.integers(len(self.choices)))]
        edge = rng.random()
        if edge < 0.1:  # the ends of a range are often where today's technology sits
            return float(self.low)
        if edge < 0.2:
            return float(self.high)
        if self.log:
            return float(math.exp(rng.uniform(math.log(self.low), math.log(self.high))))
        return float(rng.uniform(self.low, self.high))

    def perturb(self, value: Any, rng: np.random.Generator, scale: float) -> Any:
        if self.choices is not None:
            return self.sample(rng) if rng.random() < 0.2 else value
        if self.log:
            out = value * math.exp(rng.normal(0.0, scale * math.log(self.high / self.low)))
        else:
            out = value + rng.normal(0.0, scale * (self.high - self.low))
        return float(min(max(out, self.low), self.high))

    def parse(self, text: str) -> Any:
        if self.choices is not None:
            for c in self.choices:
                if str(c).lower() == text.lower():
                    return c
            if all(isinstance(c, bool) for c in self.choices):
                return text.lower() in ("1", "true", "yes", "on")
            raise ValueError(f"{self.name} must be one of {self.choices}")
        return float(text)


@dataclass
class Route:
    name: str
    section: str
    title: str
    params: list
    run: Callable
    dt: bool = True  # D-T fuel: spin polarisation and a blanket can apply

    def defaults(self) -> dict:
        return {p.name: p.default for p in self.params}


SPIN = Param("spin", choices=(False, True), default=False, help="spin-polarised D-T (cross section x1.5, untested)")


# --- the routes -------------------------------------------------------------------------------

def _cold_catcher(p: dict, s: Settings) -> Outcome:
    n = P.lid_molecule_density()
    targets = {"T": n} if p["catcher"] == "tritide" else {"D": n, "Li6": n}
    gain = float(EN.thick_target(p["ion_kev"], "D", "LiD", targets)["gain"][0])
    fusion = gain  # per joule of ion energy
    dt = fusion / (Q_DT_MEV * MEV) if p["catcher"] == "tritide" else 0.0
    return Outcome("cold_catcher", p, fusion, 1.0, {"laser": 1.0 / s.laser_to_ions}, dt,
                   notes="per joule of crystal ions; section 1")


def _crystal_for_ion_energy(mean_ion_kev: float, fuel: str = "DT") -> tuple[float, float]:
    """Crystal radius (m) whose Coulomb explosion gives that mean ion energy, and its stripping intensity (W/cm^2)."""
    rho = EN.charge_density(fuel)
    per_r2 = 0.9 * 0.6 * P.E_CHARGE * rho / (3 * P.EPS0) / KEV
    radius = math.sqrt(mean_ion_kev / per_r2)
    return radius, P.intensity_w_cm2(rho * radius / (3 * P.EPS0))


def _crystal_ball(p: dict, s: Settings) -> Outcome:
    r = EN.array_burn(p["fuel"], p["crystal_nm"], p["fill"], p["ball_mm"] * 1e-3,
                      fast_in_final_plasma=p["bound"] == "generous")
    parts, _ = EN.FUELS[p["fuel"]]
    mean_mass = sum(f * EN.SPECIES[sp][0] for sp, f in parts.items()) / sum(parts.values()) * P.AMU
    q = Q_DT_MEV if p["fuel"] == "DT" else P.Q_VALUE_MEV["Li6"]
    most = r["mass_g"] * 1e-3 / mean_mass / 2 * q * MEV  # every pair of nuclei fused once
    fusion = deplete(max(r["gain_after_radiation"], 0.0) * r["energy_in_J"], most)
    laser = r["energy_in_J"] / s.laser_to_ions
    reqs = [Requirement("laser energy per shot (J)", laser), Requirement("picosecond laser energy (J)", laser),
            Requirement("laser peak power (W)", laser / max(r["burn_time_s"], 1e-15)),
            Requirement("laser intensity (W/cm²)", r["stripping_intensity_W_cm2"])]
    dt = fusion / (Q_DT_MEV * MEV) if p["fuel"] == "DT" else 0.0
    return Outcome("crystal_ball", p, fusion, r["energy_in_J"], {"laser": laser}, dt, reqs,
                   notes="untamped ball of exploding crystals; section 3")


def _tamped_ball(p: dict, s: Settings) -> Outcome:
    ball = C.Ball(radius_m=p["radius_mm"] * 1e-3, fuel=p["fuel"], t_i_kev=p["t_i_kev"],
                  t_e_kev=1.0 if p["ions_first"] else None, fill=p["fill"], b_tesla=p["b_tesla"], tamper="gold",
                  hole_fraction=0.01)
    r = C.burn(ball)
    laser = r["energy_in_J"] / s.laser_to_ions
    mean_ion = 1.5 * p["t_i_kev"] if p["ions_first"] else 3 * p["t_i_kev"]
    _, intensity = _crystal_for_ion_energy(mean_ion)
    power = intensity * math.pi * (p["radius_mm"] * 0.1) ** 2
    reqs = [Requirement("laser energy per shot (J)", laser), Requirement("picosecond laser energy (J)", laser),
            Requirement("laser peak power (W)", power), Requirement("laser intensity (W/cm²)", intensity),
            Requirement("magnetic field on the target (T)", p["b_tesla"])]
    dt = r["fusion_J"] / (Q_DT_MEV * MEV) if p["fuel"] == "DT" else 0.0
    return Outcome("tamped_ball", p, r["fusion_J"], r["energy_in_J"], {"laser": laser}, dt, reqs,
                   notes=f"gold tamper, 1% entrance hole; {r['stopped']}; section 12")


@lru_cache(maxsize=1)
def _fast_ignition_source() -> dict:
    """Section 7: deuterons that stop in 1.2 g/cm2 of 400 g/cm3 D-T, from Coulomb-exploding 6LiD crystals."""
    rho = 400.0
    n = rho * 1e3 / (2.5 * P.AMU)
    plasma = EN.Plasma({"D": n / 2, "T": n / 2}, t_e=5.0, t_i=5.0)
    energies = (2e3, 5e3, 1e4, 2e4, 4e4)
    ranges = [EN.beam_in_plasma(e, "D", plasma)["path_m"] * rho * 100 for e in energies]
    e_needed = float(np.interp(1.2, ranges, energies))
    rho_q = EN.charge_density("6LiD")
    radius = math.sqrt(e_needed / (0.9 * P.E_CHARGE * rho_q / (3 * P.EPS0) / KEV))
    e = np.linspace(e_needed / 400, e_needed, 400)
    w = np.sqrt(e)
    window = float(np.sum((w * e)[e >= 0.75 * e_needed]) / (np.sum(w * e) + 3 * np.sum(w * e)))
    e_strip = rho_q * radius / (3 * P.EPS0)
    pulse = P.Pulse(e0=e_strip, wavelength=0.8e-6, fwhm=6e-15)
    q = rho_q * 4 / 3 * math.pi * radius**3
    explosion = 0.6 * P.K_COULOMB * q * q / radius
    efficiency = min(1.0, explosion / (pulse.fluence() * math.pi * radius**2))
    return {"deuteron_keV": e_needed, "window": window, "efficiency": efficiency,
            "intensity_W_cm2": P.intensity_w_cm2(e_strip)}


FAST_IGNITION_YIELD_J = 30e6  # design yield of ion fast ignition (section 7)
FAST_IGNITION_COMPRESSION_J = 250e3  # compression laser of HiPER-class designs
FAST_IGNITION_SPOT_J = 10e3  # ions needed in the 20 um spot
SHAPED_CRYSTAL_BOOST = 6.0  # double-layer crystals beam 5-10 times better (section 7)


def _fast_ignition(p: dict, s: Settings) -> Outcome:
    src = _fast_ignition_source()
    d = p["distance_um"] * 1e-6
    theta = math.atan(10e-6 / d)
    geometric = (1 - math.cos(theta)) / 2 * (SHAPED_CRYSTAL_BOOST if p["shaped"] else 1.0)
    ignitor = FAST_IGNITION_SPOT_J / (src["efficiency"] * src["window"] * min(geometric, 1.0))
    laser = FAST_IGNITION_COMPRESSION_J + ignitor
    reqs = [Requirement("laser energy per shot (J)", laser), Requirement("picosecond laser energy (J)", ignitor),
            Requirement("laser peak power (W)", ignitor / 6e-15),
            Requirement("laser intensity (W/cm²)", src["intensity_W_cm2"]),
            Requirement("crystals' distance from the core (m)", d)]
    fusion = FAST_IGNITION_YIELD_J
    return Outcome("fast_ignition", p, fusion, FAST_IGNITION_SPOT_J,
                   {"laser": laser}, fusion / (Q_DT_MEV * MEV), reqs,
                   notes="crystals as the ion spark for 400 g/cm3 D-T; yield fixed at the 30 MJ design; section 7")


STICKING = {"PSI 1987 (0.45%)": 0.0045, "PSI direct (0.56%)": 0.0056, "theory (0.65%)": 0.0065}
CYCLE_RATE_PER_PHI = 2.1e8 / 1.3  # 1/s per unit of liquid-hydrogen density, from the Los Alamos record


def _muon(p: dict, s: Settings) -> Outcome:
    y = ALT.fusions_per_muon(CYCLE_RATE_PER_PHI * p["phi"], STICKING[p["sticking"]])
    fusion = y * Q_DT_MEV * MEV
    beam = s.muon_cost_gev * 1e3 * MEV
    reqs = [Requirement("liquid or solid D-T density (× liquid hydrogen)", p["phi"])]
    return Outcome("muon_catalysis", p, fusion, beam, {"accelerator": beam}, y, reqs,
                   notes=f"{y:.0f} fusions per muon; section 9")


def _muon_squeezed(p: dict, s: Settings) -> Outcome:
    rho, core = p["rho_g_cc"], p["core_um"] * 1e-6
    hold = D.squeeze_lifetime(rho, core)
    sticking = D.muon_reactivation(rho)["sticking"]
    y = min(1 / sticking if sticking > 0 else INF, D.DTMU_FUSION_RATE * hold)
    cc = ALT.cold_compression(rho)
    electrons = 4 / 3 * math.pi * core**3 * rho * 1e3 / (0.5 * (P.MASS_U["D"] + P.MASS_U["T"]) * P.AMU)
    squeeze = electrons * cc["electron_energy_keV"] * KEV
    muons = p["muons"]
    fusion = muons * y * Q_DT_MEV * MEV
    beam = muons * s.muon_cost_gev * 1e3 * MEV
    laser = squeeze / s.implosion_efficiency
    reqs = [Requirement("cold squeezed density (g/cm³)", rho),
            Requirement("radial squeeze of the fuel", (rho / 0.25) ** (1 / 3)),
            Requirement("laser energy per shot (J)", laser),
            Requirement("muons per second", muons * ALT.MUON_DECAY_RATE)]  # all made within one muon lifetime
    return Outcome("muon_squeezed", p, fusion, squeeze, {"laser": laser, "accelerator": beam}, muons * y, reqs,
                   notes=f"{y:.0f} fusions per muon in a {hold * 1e12:.0f} ps squeeze; section 17")


def _pycnonuclear(p: dict, s: Settings) -> Outcome:
    cc = ALT.cold_compression(p["rho_g_cc"])
    laser = cc["compression_J"] / s.implosion_efficiency
    reqs = [Requirement("cold squeezed density (g/cm³)", p["rho_g_cc"]),
            Requirement("radial squeeze of the fuel", cc["convergence_from_solid"]),
            Requirement("laser energy per shot (J)", laser)]
    fusion = cc["yield_if_burnt_J"]
    return Outcome("pycnonuclear", p, fusion, cc["compression_J"], {"laser": laser}, fusion / (Q_DT_MEV * MEV), reqs,
                   notes="cold fuel squeezed until it fuses, burnt completely (generous); section 10")


def _projectiles(p: dict, s: Settings) -> Outcome:
    shot = ALT.crystal_projectile(p["radius_nm"], p["kev_per_deuteron"])
    gain = float(EN.thick_target(p["kev_per_deuteron"], "D", "TiD2", {"T": 1e29})["gain"][0])
    fusion = gain * shot["kinetic_J"]
    reqs = [Requirement("projectile speed (m/s)", shot["speed_m_s"])]
    return Outcome("projectiles", p, fusion, shot["kinetic_J"], {"accelerator": shot["kinetic_J"]},
                   fusion / (Q_DT_MEV * MEV), reqs, notes="whole crystals fired into a tritide target; section 16")


def _fuel_volume(radius_um: float) -> float:
    return 4 / 3 * math.pi * (radius_um * 1e-6) ** 3


def _pauli_beam(p: dict, s: Settings) -> Outcome:
    old = D.ELECTRON_DRAG_FACTOR
    D.ELECTRON_DRAG_FACTOR = p["drag"]
    try:
        best = D.best_in_flight("D", "DT", p["rho_g_cc"])
    finally:
        D.ELECTRON_DRAG_FACTOR = old
    _, n_e = D.fuel_ions("DT", p["rho_g_cc"])
    volume = _fuel_volume(p["fuel_um"])
    squeeze = n_e * D.squeeze_energy_kev(n_e) * KEV * volume
    beam = squeeze  # the trick lasts while the fuel stays cold (section 17)
    fusion = best["gain"] * beam
    laser = squeeze / s.implosion_efficiency + beam / s.laser_to_ions
    reqs = [Requirement("cold squeezed density (g/cm³)", p["rho_g_cc"]),
            Requirement("radial squeeze of the fuel", (p["rho_g_cc"] / 0.25) ** (1 / 3)),
            Requirement("laser energy per shot (J)", laser)]
    return Outcome("pauli_beam", p, fusion, squeeze + beam, {"laser": laser}, fusion / (Q_DT_MEV * MEV), reqs,
                   notes=f"deuteron returns {_num(best['gain'])} of its energy in flight; section 17")


def _son_fisch(p: dict, s: Settings) -> Outcome:
    h = D.hot_ion_balance("pB", p["rho_g_cc"], p["t_i_kev"])
    ions, n_e = D.fuel_ions(D.pb_fuel(0.15), p["rho_g_cc"])
    pressure = sum(ions.values()) * p["t_i_kev"] * KEV + 2 / 3 * h["squeeze_J_m3"]
    hold = p["fuel_um"] * 1e-6 / (4 * math.sqrt(5 / 3 * pressure / (p["rho_g_cc"] * 1e3)))
    volume = _fuel_volume(p["fuel_um"])
    boron = ions["B11"] * volume
    fusion = deplete(h["fusion_W_m3"] * hold * min(1.0, h["self_heating"]) * volume,
                     boron * P.Q_VALUE_MEV["pB"] * MEV)
    squeeze = h["squeeze_J_m3"] * volume
    heat = h["ion_heat_J_m3"] * volume
    laser = squeeze / s.implosion_efficiency + heat / s.laser_to_ions
    reqs = [Requirement("cold squeezed density (g/cm³)", p["rho_g_cc"]),
            Requirement("radial squeeze of the fuel", (p["rho_g_cc"] / 1.0) ** (1 / 3)),
            Requirement("laser energy per shot (J)", laser)]
    return Outcome("son_fisch_pB", p, fusion, squeeze + heat, {"laser": laser}, 0.0, reqs,
                   notes=f"ions self-heating ratio {_num(h['self_heating'])}; section 17")


def _opaque_pb(p: dict, s: Settings) -> Outcome:
    ball = W.opaque_ball(p["rho_g_cc"], p["t_i_kev"], p["radius_um"] * 1e-6)
    _, n_e = D.fuel_ions(D.pb_fuel(0.15), p["rho_g_cc"])
    squeeze = n_e * D.squeeze_energy_kev(n_e) * KEV * _fuel_volume(p["radius_um"])
    heat = max(ball["energy_J"] - squeeze, 0.0)
    ions, _ = D.fuel_ions(D.pb_fuel(0.15), p["rho_g_cc"])
    fusion = deplete(ball["gain"] * ball["energy_J"] * min(1.0, ball["burning"]),
                     ions["B11"] * _fuel_volume(p["radius_um"]) * P.Q_VALUE_MEV["pB"] * MEV)
    laser = squeeze / s.implosion_efficiency + heat / s.laser_to_ions
    reqs = [Requirement("cold squeezed density (g/cm³)", p["rho_g_cc"]),
            Requirement("hot-spot areal density (g/cm²)", ball["rhoR_g_cm2"]),
            Requirement("hot-spot pressure (Gbar)", ball["pressure_Pa"] / 1e14),
            Requirement("laser energy per shot (J)", laser)]
    return Outcome("opaque_pB", p, fusion, ball["energy_J"], {"laser": laser}, 0.0, reqs,
                   notes=f"{ball['trapped']:.0%} of the X-rays trapped; section 18")


PINCH_LENGTH_M = 0.01
PINCH_WALL_M = 0.01


def _pinch(p: dict, s: Settings) -> Outcome:
    current = p["current_ma"] * 1e6
    one = W.fibre_pinch(current, t_kev=p["t_kev"], alfven_times=p["alfven_times"])
    radius = one["fibre_radius_m"] / p["squeeze"]
    ratio = one["radius_m"] / radius  # fusion / heat grows as 1 / radius, before the fuel burns down
    heat = one["heat_J_m"] * PINCH_LENGTH_M
    raw = ratio * heat
    fusion = deplete(raw, one["ions_per_m"] / 2 * Q_DT_MEV * MEV * PINCH_LENGTH_M)
    # X-rays: bremsstrahlung over the same n^2 x lifetime as the fusion. Whatever the alphas don't make up (they
    # escape below ~1.35 MA) the generator must resupply, or the pinch cools.
    t_kev = p["t_kev"]
    brem_over_fusion = 5.35e-37 * math.sqrt(t_kev) * 4 / (EN.reactivity("D", "T", t_kev) * Q_DT_MEV * MEV)
    alphas = D.ALPHA_KEV["DT"] / (Q_DT_MEV * 1e3) if one["alphas_held"] else 0.0
    resupply = max(0.0, brem_over_fusion * raw - alphas * fusion)
    field_j = MU0 * current**2 / (4 * math.pi) * (math.log(max(PINCH_WALL_M / radius, 1.0)) + 0.25) * PINCH_LENGTH_M
    pressure_gbar = MU0 * current**2 / (8 * math.pi**2 * radius**2) / 1e14  # Bennett: mean plasma pressure
    drive = field_j + heat + resupply
    reqs = [Requirement("pinch current (A)", current), Requirement("energy into the pinch (J)", drive),
            Requirement("stable Alfvén times of a dense hot pinch", p["alfven_times"]),
            Requirement("radial squeeze of the fuel", p["squeeze"]),
            Requirement("hot-spot pressure (Gbar)", pressure_gbar)]
    return Outcome("fibre_pinch", p, fusion, heat + resupply, {"pulsed power": drive}, fusion / (Q_DT_MEV * MEV),
                   reqs, notes=f"1 cm of D-T fibre; fusion / plasma heat {_num(fusion / heat)}; X-rays to resupply "
                   f"{_num(resupply / heat)} x the heat; section 19")


CHANNEL_VALUES = {"measured": ((2.0, 0.9), (1.4, 0.95)), "generous": ((3.0, 0.5), (2.0, 0.9))}


def _channelling(p: dict, s: Settings) -> Outcome:
    e = p["ion_kev"]
    random = float(EN.thick_target(e, "D", "TiD2", {"T": 1e29})["gain"][0])
    axial, planar = CHANNEL_VALUES[p["values"]]
    if p["source"] == "aligned beam":
        ret = W.channelled_return(random, [(1.0, *axial, 1.0)])
        drivers = {"accelerator": 1.0}
    else:
        psi = W.critical_angle(e, 1, 22, 4.45e-10 / math.sqrt(2))
        ret = W.channelled_return(random, [(W.axial_acceptance(psi, 26), *axial, 0.5),
                                           (W.planar_acceptance(psi / 2.5, 13), *planar, 0.5)])
        drivers = {"laser": 1.0 / s.laser_to_ions}
    return Outcome("channelling", p, ret, 1.0, drivers, ret / (Q_DT_MEV * MEV),
                   notes="per joule of deuterons into a titanium tritide crystal; section 20")


ROUTES = {r.name: r for r in [
    Route("cold_catcher", "1", "Crystal ions stopping in a cold catcher", [
        Param("ion_kev", 10.0, 5000.0, default=300.0, help="deuteron energy (keV)"),
        Param("catcher", choices=("tritide", "LiD"), default="tritide")], _cold_catcher),
    Route("crystal_ball", "3", "Ball of exploding crystals, untamped", [
        Param("fuel", choices=("DT", "6LiDT", "6LiD"), default="DT"),
        Param("crystal_nm", 1.5, 50.0, default=15.0, help="crystal radius (nm)"),
        Param("fill", 0.1, 1.0, default=1.0, help="share of the ball filled with crystals"),
        Param("ball_mm", 0.01, 100.0, default=3.0, help="ball radius (mm)"),
        Param("bound", choices=("generous", "conservative"), default="generous")], _crystal_ball),
    Route("tamped_ball", "12", "Ball held by a gold tamper and a magnetic field", [
        Param("radius_mm", 0.01, 10.0, default=0.3, help="fuel radius (mm)"),
        Param("t_i_kev", 5.0, 80.0, default=20.0, help="starting ion temperature (keV)"),
        Param("b_tesla", 1.0, 3e4, default=3000.0, help="field along the wall (T)"),
        Param("fill", 0.01, 1.0, default=1.0, help="fuel density / solid"),
        Param("ions_first", choices=(False, True), default=False, help="crystals heat the ions first"),
        Param("fuel", choices=("DT", "DD"), default="DT"), SPIN], _tamped_ball),
    Route("fast_ignition", "7", "Exploding crystals as the spark for compressed D-T", [
        Param("distance_um", 30.0, 3000.0, default=100.0, help="crystals' distance from the spark (um)"),
        Param("shaped", choices=(False, True), default=True, help="double-layer crystals that beam")],
        _fast_ignition),
    Route("muon_catalysis", "9", "Muon-catalysed fusion in liquid or solid D-T", [
        Param("phi", 0.1, 10.0, default=1.2, help="density (x liquid hydrogen)"),
        Param("sticking", choices=tuple(STICKING), default="PSI 1987 (0.45%)")], _muon),
    Route("muon_squeezed", "17", "Muon catalysis in squeezed fuel (stuck muons freed)", [
        Param("rho_g_cc", 10.0, 2000.0, default=100.0, help="squeezed D-T density (g/cm3)"),
        Param("core_um", 10.0, 3000.0, default=100.0, help="squeezed core radius (um)"),
        Param("muons", 1e8, 1e18, default=1e12, help="muons in the core during the squeeze")], _muon_squeezed),
    Route("pycnonuclear", "10", "Cold fuel squeezed until it fuses", [
        Param("rho_g_cc", 1e3, 1e8, default=1e7, help="cold D-T density (g/cm3)")], _pycnonuclear),
    Route("projectiles", "16", "Whole crystals fired by an accelerator", [
        Param("radius_nm", 2.0, 50.0, default=15.0, help="crystal radius (nm)"),
        Param("kev_per_deuteron", 10.0, 3000.0, default=100.0)], _projectiles),
    Route("pauli_beam", "17", "Crystal ions into cold fuel squeezed until it is degenerate", [
        Param("rho_g_cc", 25.0, 2.5e5, default=2500.0, help="squeezed D-T density (g/cm3)"),
        Param("fuel_um", 10.0, 1000.0, default=100.0, help="squeezed fuel radius (um)"),
        Param("drag", 1.0, 1.6, log=False, default=1.0, help="electron drag x linear response"), SPIN],
        _pauli_beam),
    Route("son_fisch_pB", "17", "Hot proton-boron ions in degenerate electrons", [
        Param("rho_g_cc", 1e4, 1e7, default=3.5e5, help="fuel density (g/cm3)"),
        Param("t_i_kev", 50.0, 300.0, default=100.0, help="ion temperature (keV)"),
        Param("fuel_um", 0.1, 100.0, default=5.0, help="fuel radius (um)")], _son_fisch, dt=False),
    Route("opaque_pB", "18", "Proton-boron ball that traps its own X-rays", [
        Param("rho_g_cc", 100.0, 1e5, default=1e4, help="fuel density (g/cm3)"),
        Param("t_i_kev", 150.0, 500.0, default=200.0, help="ion temperature (keV)"),
        Param("radius_um", 0.5, 1e4, default=20.0, help="ball radius (um)")], _opaque_pb, dt=False),
    Route("fibre_pinch", "19", "D-T fibre pinched by its own current", [
        Param("current_ma", 0.3, 100.0, default=26.0, help="current (MA)"),
        Param("t_kev", 3.0, 30.0, default=10.0, help="Bennett temperature (keV)"),
        Param("alfven_times", 1.0, 1e4, default=1.0, help="how long the pinch stays stable"),
        Param("squeeze", 1.0, 300.0, default=6.0, help="radial squeeze of the fibre"), SPIN], _pinch),
    Route("channelling", "20", "Ions channelled through a tritide crystal", [
        Param("ion_kev", 50.0, 3000.0, default=300.0),
        Param("source", choices=("aligned beam", "exploding crystal"), default="aligned beam"),
        Param("values", choices=("measured", "generous"), default="measured")], _channelling),
]}


def evaluate(route: str, params: dict | None = None, settings: Settings | None = None) -> Outcome:
    """Run one route with the given parameters (missing ones take their defaults)."""
    r = ROUTES[route]
    p = r.defaults()
    p.update(params or {})
    s = settings or Settings()
    old_spin = EN.DT_SPIN_FACTOR
    EN.DT_SPIN_FACTOR = 1.5 if p.get("spin") else 1.0
    try:
        out = r.run(p, s)
    finally:
        EN.DT_SPIN_FACTOR = old_spin
    out.params = dict(p)
    out.settings = s
    return out


# --- the search ---------------------------------------------------------------------------------

BUDGETS = (1.0, 10.0, 100.0, INF)


def _safe(route: str, params: dict, settings: Settings) -> Outcome | None:
    try:
        out = evaluate(route, params, settings)
    except (ValueError, ZeroDivisionError, OverflowError, FloatingPointError):
        return None
    if not all(math.isfinite(x) for x in (out.fusion_J, out.electricity_in_J)) or out.electricity_in_J <= 0:
        return None
    return out


def search(routes=None, samples: int = 200, refine: int = 80, seed: int = 1, settings: Settings | None = None,
           budgets=BUDGETS) -> dict:
    """Best electric gain of each route within each stretch budget.

    Random sampling of every parameter (log-uniform for ranges), then a
    hill-climb from the best design within each budget. Returns
    {route: {budget: Outcome or None}}.
    """
    settings = settings or Settings()
    rng = np.random.default_rng(seed)
    results = {}
    with np.errstate(all="ignore"):
        for name in routes or ROUTES:
            route = ROUTES[name]
            pool = [o for o in (_safe(name, {q.name: q.sample(rng) for q in route.params}, settings)
                                for _ in range(samples)) if o is not None]
            default = _safe(name, {}, settings)
            if default is not None:
                pool.append(default)
            best_by_budget = {}
            for budget in budgets:
                inside = [o for o in pool if o.stretch <= budget]
                if not inside:
                    best_by_budget[budget] = None
                    continue
                best = max(inside, key=lambda o: o.electric_gain)
                for step in range(refine):
                    scale = 0.15 * (1 - step / refine) + 0.01
                    trial = {q.name: q.perturb(best.params[q.name], rng, scale) for q in route.params}
                    out = _safe(name, trial, settings)
                    if out is not None and out.stretch <= budget and out.electric_gain > best.electric_gain:
                        best = out
                        pool.append(out)
                best_by_budget[budget] = best
            results[name] = best_by_budget
    return results


def budget_label(budget: float) -> str:
    return "unlimited" if math.isinf(budget) else ("today's technology" if budget == 1 else f"≤ {budget:g}× today")


# --- command line ---------------------------------------------------------------------------------

def _fmt(x: float) -> str:
    if x == 0 or not math.isfinite(x):
        return str(x)
    return f"{x:.3g}"


def _print_outcome(o: Outcome) -> None:
    print(f"{o.route}: {ROUTES[o.route].title} (docs/energy.md section {ROUTES[o.route].section})")
    print("  parameters: " + ", ".join(f"{k}={_fmt(v) if isinstance(v, float) else v}" for k, v in o.params.items()))
    print(f"  fusion {_fmt(o.fusion_J)} J, energy in the fuel {_fmt(o.fuel_J)} J, "
          + ", ".join(f"{k} {_fmt(v)} J" for k, v in o.drivers.items()))
    print(f"  fusion / driver energy {_fmt(o.target_gain)}; heat with {o.settings.blanket} blanket {_fmt(o.heat_J)} J")
    print(f"  electricity in {_fmt(o.electricity_in_J)} J -> electric gain {_fmt(o.electric_gain)}"
          + ("  (NET ENERGY)" if o.electric_gain > 1 else ""))
    for r in o.all_requirements:
        best = f"{_fmt(r.best)} ({TODAY[r.name][1]})" if r.name in TODAY else "never demonstrated"
        print(f"  needs {r.name}: {_fmt(r.needed)}  vs best {best}  -> stretch {_fmt(r.stretch)}")
    print(f"  stretch {_fmt(o.stretch)}" + (f", limited by {o.limiting}" if o.limiting else "") + f"; {o.notes}")


def _settings_from(args: argparse.Namespace) -> Settings:
    return replace(Settings(), laser_efficiency=args.laser_efficiency, thermal_efficiency=args.thermal_efficiency,
                   blanket=args.blanket, k_eff=args.k_eff)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(prog="python -m nanofusion sim", description=__doc__.split("\n\n")[0])
    sub = parser.add_subparsers(dest="command", required=True)
    sub.add_parser("list", help="list the routes and their parameters")
    run = sub.add_parser("run", help="run one route: ROUTE name=value ...")
    run.add_argument("route", choices=sorted(ROUTES))
    run.add_argument("values", nargs="*", help="parameter=value")
    find = sub.add_parser("search", help="search every route for the best electric gain")
    find.add_argument("--routes", default="", help="comma-separated route names (default: all)")
    find.add_argument("--samples", type=int, default=200)
    find.add_argument("--refine", type=int, default=80)
    find.add_argument("--seed", type=int, default=1)
    find.add_argument("-o", "--output", default="", help="write the results as JSON")
    for p in (run, find):
        p.add_argument("--blanket", choices=("none", "lithium", "uranium"), default="lithium")
        p.add_argument("--k-eff", type=float, default=0.0, help="uranium blanket k_eff")
        p.add_argument("--laser-efficiency", type=float, default=0.10)
        p.add_argument("--thermal-efficiency", type=float, default=0.40)
    args = parser.parse_args(argv)
    if args.command == "list":
        for r in ROUTES.values():
            print(f"{r.name} (section {r.section}): {r.title}")
            for q in r.params:
                span = "/".join(str(c) for c in q.choices) if q.choices else f"{_fmt(q.low)}–{_fmt(q.high)}"
                print(f"    {q.name} = {q.default}  [{span}]  {q.help}")
        return 0
    settings = _settings_from(args)
    if args.command == "run":
        route = ROUTES[args.route]
        params = {}
        known = {q.name: q for q in route.params}
        for item in args.values:
            key, _, value = item.partition("=")
            if key not in known or not value:
                print(f"unknown parameter {item!r}; parameters: {', '.join(known)}", file=sys.stderr)
                return 2
            params[key] = known[key].parse(value)
        _print_outcome(evaluate(args.route, params, settings))
        return 0
    routes = [r for r in args.routes.split(",") if r] or None
    if routes and any(r not in ROUTES for r in routes):
        print(f"unknown route; routes: {', '.join(ROUTES)}", file=sys.stderr)
        return 2
    results = search(routes, args.samples, args.refine, args.seed, settings)
    for budget in BUDGETS:
        print(f"\n== Best electric gain, {budget_label(budget)} ==")
        ranked = sorted(((name, by[budget]) for name, by in results.items() if by[budget] is not None),
                        key=lambda t: -t[1].electric_gain)
        for name, o in ranked:
            print(f"  {name:16s} electric gain {_fmt(o.electric_gain):>9s}  stretch {_fmt(o.stretch):>8s}"
                  + (f"  limited by {o.limiting}" if o.limiting else ""))
    if args.output:
        with open(args.output, "w") as f:
            json.dump({name: {budget_label(b): (o.summary() if o else None) for b, o in by.items()}
                       for name, by in results.items()}, f, indent=1, default=float)
    return 0


if __name__ == "__main__":
    sys.exit(main())
