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
from collections.abc import Callable
from dataclasses import dataclass, field, replace
from functools import lru_cache
from typing import Any

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
    # Pulsed power delivered to a pinch -> fuel of an X-ray-driven capsule: ~2/3 becomes X-rays, a capsule absorbs
    # ~6-11% of the hohlraum's energy (Z dynamic hohlraum; Bailey et al. 2004, Vesey 2007), and ~6% of that ends up
    # in the fuel, as at NIF.
    xray_drive_efficiency: float = 0.0035
    muon_cost_gev: float = 4.7  # beam energy per muon, optimised accelerator design (section 9)
    blanket: str = "lithium"  # "none", "lithium" (breeds tritium: +4.78 MeV per D-T fusion) or "uranium"
    k_eff: float = 0.0  # uranium blanket multiplication (subcritical), used when blanket == "uranium"

    def efficiency(self, driver: str) -> float:
        return {"laser": self.laser_efficiency, "pulsed power": self.pulsed_power_efficiency,
                "accelerator": self.accelerator_efficiency, "electricity": 1.0}[driver]


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
    "fuel areal density (g/cm²)": (0.6, "NIF's best shots"),
    "fuel adiabat (lower is harder)": (1.5, "NIF's 2012 low-foot designs"),
    "compressed magnetic field (T)": (7.1e3, "NIF magnetised implosion (Sio et al. 2023)"),
    "electron screening potential (eV)": (1500.0, "6Li(d,a) in palladium, the largest measured"),
    "antiprotons per shot": (6e12, "the most ever stored, Fermilab's Recycler"),
}
LOWER_IS_HARDER = {"crystals' distance from the core (m)", "fuel adiabat (lower is harder)"}
# Untested theory: a demand not in TODAY is an idea never shown to work. It adds no technology stretch, but the
# design is flagged, and the search leaves it out unless asked for theory (search(..., theory=True)).
# Route parameters that switch on untested theory, with their "off" values:
THEORY_OFF = {"spin": False, "alpha_channelling": 0.0, "direct_conversion": False, "rf_stripping": "none",
              "stability": "breaks"}
DIRECT_EFFICIENCY = 0.6  # charged fusion products turned straight into electricity (design claims)


@dataclass
class Requirement:
    name: str
    needed: float

    @property
    def best(self) -> float:
        return TODAY[self.name][0] if self.name in TODAY else math.nan

    @property
    def untested(self) -> bool:
        return self.name not in TODAY

    @property
    def stretch(self) -> float:
        if self.untested:
            return 1.0
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
    direct_share: float = 0.0  # share of the fusion energy turned straight into electricity (direct conversion)

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
    def electricity_out_J(self) -> float:
        direct = self.fusion_J * self.direct_share
        return (self.heat_J - direct) * self.settings.thermal_efficiency + direct * DIRECT_EFFICIENCY

    @property
    def electric_gain(self) -> float:
        e_in = self.electricity_in_J
        return self.electricity_out_J / e_in if e_in > 0 else 0.0

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
        """How far beyond today's best technology this design goes (untested theory not counted)."""
        return max([r.stretch for r in self.all_requirements], default=1.0)

    @property
    def untested(self) -> list:
        """The untested theoretical ideas this design relies on."""
        return [r.name for r in self.all_requirements if r.untested]

    @property
    def limiting(self) -> str:
        reqs = [r for r in self.all_requirements if not r.untested]
        if not reqs:
            return ""
        worst = max(reqs, key=lambda r: r.stretch)
        return worst.name if worst.stretch > 1 else ""

    def summary(self) -> dict:
        return {"route": self.route, "params": self.params, "fusion_J": self.fusion_J, "fuel_J": self.fuel_J,
                "drivers_J": self.drivers, "electricity_in_J": self.electricity_in_J, "heat_J": self.heat_J,
                "target_gain": self.target_gain, "electric_gain": self.electric_gain, "stretch": self.stretch,
                "limiting": self.limiting, "untested": self.untested, "notes": self.notes,
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
    effort: float = 1.0  # share of the search's samples spent on this route (slow routes get fewer)

    def defaults(self) -> dict:
        return {p.name: p.default for p in self.params}


SPIN = Param("spin", choices=(False, True), default=False, help="spin-polarised D-T (cross section x1.5, untested)")
RF_STRIPPING = Param("rf_stripping", choices=("none", "Kulsrud/TEK"), default="none",
                     help="free stuck muons by field re-acceleration (untested; costs cryogenic cooling)")
PB_THEORY = [Param("alpha_channelling", 0.0, 0.9, log=False, default=0.0,
                   help="share of the alphas' electron heating handed to ions by waves (untested)"),
             Param("direct_conversion", choices=(False, True), default=False,
                   help="charged products turned straight into electricity (untested at scale)")]


# --- the routes -------------------------------------------------------------------------------

def _ion_source(p: dict, s: Settings) -> dict:
    """Driver energy for one joule of fast ions: exploding crystals (laser) or an accelerator."""
    return {"laser": 1.0 / s.laser_to_ions} if p["source"] == "crystals" else {"accelerator": 1.0}


def _screening(p: dict) -> list:
    return [Requirement("electron screening potential (eV)", p["screening_ev"])] if p["screening_ev"] > 0 else []


def _cold_catcher(p: dict, s: Settings) -> Outcome:
    n = P.lid_molecule_density()
    targets = {"T": n} if p["catcher"] == "tritide" else {"D": n, "Li6": n}
    gain = float(EN.thick_target(p["ion_kev"], "D", "LiD", targets, screening_ev=p["screening_ev"])["gain"][0])
    fusion = gain  # per joule of ion energy
    dt = fusion / (Q_DT_MEV * MEV) if p["catcher"] == "tritide" else 0.0
    return Outcome("cold_catcher", p, fusion, 1.0, _ion_source(p, s), dt, _screening(p),
                   notes="per joule of fast deuterons; sections 1 and 5")


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


# Untested theory: re-accelerate stuck (mu-alpha)+ ions with fields so they are stripped (Kulsrud, PPPL-2681; the
# TEK concept, Tajima, Eliezer & Kulsrud 1988-90: effective sticking below 1e-3). Both papers find it needs the fuel
# split into thin rods kept below ~25 K, and removing each fusion's 3.5 MeV of alpha heat at that temperature costs
# ~30-35 MeV of electricity, more than the 17.6 MeV the fusion releases.
TEK_STICKING = 1e-3
TEK_COOLING_MEV = 30.0


def _muon_cycle(p: dict) -> tuple[float, float, list]:
    """Fusions per muon, cooling electricity per muon (J), and untested demands, with or without RF stripping."""
    sticking = STICKING[p["sticking"]]
    reqs, cooling = [], 0.0
    if p.get("rf_stripping", "none") != "none":
        sticking = TEK_STICKING
        reqs.append(Requirement("stuck muons re-accelerated and freed by fields (Kulsrud, TEK)", TEK_STICKING))
    y = ALT.fusions_per_muon(CYCLE_RATE_PER_PHI * p["phi"], sticking)
    if reqs:
        cooling = y * TEK_COOLING_MEV * MEV
    return y, cooling, reqs


def _muon(p: dict, s: Settings) -> Outcome:
    y, cooling, theory = _muon_cycle(p)
    fusion = y * Q_DT_MEV * MEV
    beam = s.muon_cost_gev * 1e3 * MEV
    reqs = [Requirement("liquid or solid D-T density (× liquid hydrogen)", p["phi"])] + theory
    drivers = {"accelerator": beam, "electricity": cooling} if cooling else {"accelerator": beam}
    return Outcome("muon_catalysis", p, fusion, beam, drivers, y, reqs,
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


def _theory_pb(p: dict) -> tuple[list, float]:
    """Untested theory on the proton-boron routes: alpha channelling and direct conversion."""
    reqs = []
    if p.get("alpha_channelling", 0.0) > 0:
        reqs.append(Requirement("alpha channelling (waves hand alpha energy to the ions)", p["alpha_channelling"]))
    if p.get("direct_conversion"):
        reqs.append(Requirement(f"direct conversion of fusion products at {DIRECT_EFFICIENCY:.0%}", DIRECT_EFFICIENCY))
    return reqs, (1.0 if p.get("direct_conversion") else 0.0)


def _son_fisch(p: dict, s: Settings) -> Outcome:
    h = D.hot_ion_balance("pB", p["rho_g_cc"], p["t_i_kev"], alpha_channelling=p["alpha_channelling"])
    ions, _ = D.fuel_ions(D.pb_fuel(0.15), p["rho_g_cc"])
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
    theory, direct = _theory_pb(p)
    return Outcome("son_fisch_pB", p, fusion, squeeze + heat, {"laser": laser}, 0.0, reqs + theory,
                   notes=f"ions self-heating ratio {_num(h['self_heating'])}; section 17", direct_share=direct)


def _opaque_pb(p: dict, s: Settings) -> Outcome:
    ball = W.opaque_ball(p["rho_g_cc"], p["t_i_kev"], p["radius_um"] * 1e-6, alpha_channelling=p["alpha_channelling"])
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
    theory, direct = _theory_pb(p)
    return Outcome("opaque_pB", p, fusion, ball["energy_J"], {"laser": laser}, 0.0, reqs + theory,
                   notes=f"{ball['trapped']:.0%} of the X-rays trapped; section 18", direct_share=direct)


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



# --- the dense pinch in time (zpinch.py) ---------------------------------------------------------

PINCH_STABILITY = {  # untested ideas that would keep a dense pinch stable longer
    "breaks": None,
    "rising current": "a dense pinch stays stable while its current rises (NRL 1987; disputed by 2-D MHD)",
    "sheared flow": "sheared-flow stabilisation of a pinch at fibre density (shown only below 0.6 MA, 5e24 m⁻³)",
}
PREHEAT_COUPLING = 0.5  # share of the preheat laser's energy that heats the fibre (MagLIF-like)


def _dynamic_pinch(p: dict, s: Settings) -> Outcome:
    """A D-T fibre on a pulsed-power generator, followed in time (zpinch.run).

    The generator is a capacitor bank behind 15 nH, sized to drive current_ma into the fibre at its starting
    radius in 100 ns. A laser preheats the fibre so it can carry the current. The pinch heats ohmically and by
    compression, radiates, collapses above the Pease-Braginskii current, and burns, until the m = 0 mode has
    grown for alfven_times e-folds (or never, with an untested stability idea). Each shot pays for the whole bank.
    """
    from . import zpinch as Z  # local: zpinch imports burn1d, which is only needed here

    pinch = Z.Pinch(fuel="DT", radius_um=p["radius_um"], current_ma=p["current_ma"], start_kev=p["preheat_kev"],
                    growth_limit=p["alfven_times"], rise_stabilised=p["stability"] == "rising current",
                    shear=0.1 if p["stability"] == "sheared flow" else 0.0)
    r = Z.run(pinch)
    preheat = 3 * r["line_density"] * p["preheat_kev"] * KEV * pinch.length_mm * 1e-3 / PREHEAT_COUPLING
    reqs = [Requirement("pinch current (A)", p["current_ma"] * 1e6),
            Requirement("energy into the pinch (J)", r["delivered_J"]),
            Requirement("stable Alfvén times of a dense hot pinch", p["alfven_times"]),
            Requirement("radial squeeze of the fuel", p["radius_um"] / max(r["min_radius_um"], 1e-3)),
            Requirement("laser energy per shot (J)", preheat)]
    if PINCH_STABILITY[p["stability"]]:
        reqs.append(Requirement(PINCH_STABILITY[p["stability"]], 1.0))
    return Outcome("dynamic_pinch", p, r["yield_J"], r["delivered_J"], {"pulsed power": r["bank_J"], "laser": preheat},
                   r["yield_J"] / (Q_DT_MEV * MEV), reqs,
                   notes=f"{r['stopped']} at {r['t_stop_ns']:.0f} ns; {r['peak_current_MA']:.1f} MA reached; "
                   f"{r['burn_up']:.1%} burnt; hottest {r['peak_T_keV']:.2g} keV; narrowest {r['min_radius_um']:.2g} µm; "
                   "zpinch.py")


# --- routes added later: the reference route, and combinations of ideas ---------------------------

# Burnt share of compressed D-T = rhoR / (rhoR + H_B): H_B = 6.3 g/cm2 at 20-70 keV (Fraley et al. 1974), 7.3 at its
# minimum near 40 keV (Atzeni & Meyer-ter-Vehn 2004); larger at lower temperatures.
BURN_PARAMETER = 7.0


HOT_SPOT_RHOR = 0.3  # g/cm2: an isobaric D-T hot spot ignites at rhoR ~ 0.3 g/cm2 and ~10 keV (Lindl's criterion)
HOT_SPOT_KEV = 10.0


def _hot_spot(pressure_pa: float) -> dict:
    """The igniting hot spot at a given pressure: 10 keV, rhoR = 0.3 g/cm2, density from P = 2 n kT."""
    mean_mass = 0.5 * (P.MASS_U["D"] + P.MASS_U["T"]) * P.AMU
    rho = pressure_pa * mean_mass / (2 * HOT_SPOT_KEV * KEV)  # kg/m^3
    radius = HOT_SPOT_RHOR * 10 / rho  # m (rhoR in kg/m^2 = 10 x g/cm^2)
    volume = 4 / 3 * math.pi * radius**3
    return {"radius_m": radius, "mass_kg": rho * volume, "energy_J": 1.5 * pressure_pa * volume,
            "rhoR_g_cm2": HOT_SPOT_RHOR, "T_keV": HOT_SPOT_KEV}


def _icf_target(p: dict) -> dict:
    """The compressed D-T target shared by the inertial-fusion routes: an isobaric hot spot in a dense shell."""
    rho, rhor = p["rho_g_cc"], p["rho_r"]
    _, n_e = D.fuel_ions("DT", rho)
    pressure = p["adiabat"] * ALT._pressure(n_e) * KEV  # Pa
    hot = _hot_spot(pressure)
    r_hot_cm = hot["radius_m"] * 100
    shell_cm = rhor / rho
    shell_kg = 4 / 3 * math.pi * rho * ((r_hot_cm + shell_cm) ** 3 - r_hot_cm**3) * 1e-3
    mass_kg = shell_kg + hot["mass_kg"]
    pairs = mass_kg / ((P.MASS_U["D"] + P.MASS_U["T"]) * P.AMU)
    total_rhor = rhor + hot["rhoR_g_cm2"]
    burnt = total_rhor / (total_rhor + BURN_PARAMETER / EN.DT_SPIN_FACTOR)
    cold = p["adiabat"] * D.squeeze_energy_kev(n_e) * KEV * shell_kg / (rho * 1e3) * n_e
    reqs = [Requirement("fuel areal density (g/cm²)", rhor), Requirement("cold squeezed density (g/cm³)", rho),
            Requirement("fuel adiabat (lower is harder)", p["adiabat"]),
            Requirement("hot-spot pressure (Gbar)", pressure / 1e14)]
    return {"hot": hot, "cold_J": cold, "fusion_J": burnt * pairs * Q_DT_MEV * MEV, "pressure_Pa": pressure,
            "mass_kg": mass_kg, "burnt": burnt, "requirements": reqs}


def _compressed_icf(p: dict, s: Settings) -> Outcome:
    """Ordinary inertial fusion: a hot spot inside a dense shell of D-T at the same pressure (isobaric).

    The shell's density and adiabat set the pressure (adiabat x the degenerate pressure); at that pressure the
    hot spot is at 10 keV with rhoR = 0.3 g/cm2, the usual ignition condition.
    The shell (areal density rho_r) surrounds it. Fusion: the burnt share rhoR/(rhoR + H_B) of all the fuel
    (H_B scaled by the spin factor). Energy in the fuel: the hot spot's heat plus the shell's squeeze. Driver: an
    implosion that puts 1% of its light into the fuel (laser), or 0.35% of the pulsed power delivered (X-rays
    from a pinch: see Settings); or
    (ignitor = "crystal spark") the implosion squeezes the shell and exploding shaped crystals 0.1 mm away light
    the hot spot (section 7). NIF-like inputs (950 g/cm3, adiabat 3, rho_r 0.6: ~600 Gbar, ~0.15 mg) give a
    target gain of ~2; NIF's record is 4.1.
    """
    t = _icf_target(p)
    hot, cold, fusion, pressure = t["hot"], t["cold_J"], t["fusion_J"], t["pressure_Pa"]
    mass_kg, burnt, reqs = t["mass_kg"], t["burnt"], t["requirements"]
    drive_efficiency = s.implosion_efficiency if p["driver"] == "laser" else s.xray_drive_efficiency
    if p["ignitor"] == "crystal spark":
        src = _fast_ignition_source()
        theta = math.atan(10e-6 / 100e-6)
        share = src["efficiency"] * src["window"] * (1 - math.cos(theta)) / 2 * SHAPED_CRYSTAL_BOOST
        spark = hot["energy_J"] / share
        implosion = cold / drive_efficiency
        reqs += [Requirement("picosecond laser energy (J)", spark), Requirement("laser peak power (W)", spark / 6e-15),
                 Requirement("laser intensity (W/cm²)", src["intensity_W_cm2"])]
    else:
        spark = 0.0
        implosion = (hot["energy_J"] + cold) / drive_efficiency
    if p["driver"] == "laser":
        drivers = {"laser": implosion + spark}
        reqs.append(Requirement("laser energy per shot (J)", implosion + spark))
    else:
        drivers = {"pulsed power": implosion, "laser": spark} if spark else {"pulsed power": implosion}
        reqs.append(Requirement("energy into the pinch (J)", implosion))
    return Outcome("compressed_icf", p, fusion, hot["energy_J"] + cold, drivers, fusion / (Q_DT_MEV * MEV), reqs,
                   notes=f"{mass_kg * 1e6:.2g} mg of D-T, {burnt:.0%} burnt, hot spot {hot['radius_m'] * 1e6:.0f} µm at "
                   f"{pressure / 1e14:.0f} Gbar; sections 7 and 11")


def icf_target(rho_g_cc: float, rho_r: float, adiabat: float) -> dict:
    """The compressed D-T target of compressed_icf: hot spot, cold shell, fusion if it burns (see _icf_target)."""
    return _icf_target({"rho_g_cc": rho_g_cc, "rho_r": rho_r, "adiabat": adiabat})


def icf_shell_for_mass(pressure_gbar: float, adiabat: float, mass_kg: float) -> float:
    """Cold-shell areal density (g/cm2) that makes a target at this pressure and adiabat hold mass_kg of D-T."""
    rho = icf_shell_density(pressure_gbar, adiabat)
    lo, hi = 1e-4, 20.0
    for _ in range(80):
        mid = 0.5 * (lo + hi)
        lo, hi = (mid, hi) if icf_target(rho, mid, adiabat)["mass_kg"] < mass_kg else (lo, mid)
    return 0.5 * (lo + hi)


def icf_shell_density(pressure_gbar: float, adiabat: float) -> float:
    """Cold-shell density (g/cm3) at which adiabat x the degenerate pressure equals pressure_gbar."""
    lo, hi = 1.0, 1e6
    for _ in range(100):
        mid = math.sqrt(lo * hi)
        _, n_e = D.fuel_ions("DT", mid)
        if adiabat * ALT._pressure(n_e) * KEV / 1e14 < pressure_gbar:
            lo = mid
        else:
            hi = mid
    return math.sqrt(lo * hi)


def icf_design(fuel_budget_J: float, max_pressure_gbar: float, min_adiabat: float, max_rho_r: float = 10.0,
               adiabats=None, pressures=None) -> dict | None:
    """The compressed D-T target (as compressed_icf) that releases the most fusion with at most fuel_budget_J in
    its fuel, a hot-spot pressure at most max_pressure_gbar and an adiabat at least min_adiabat.

    At a given pressure and adiabat the hot spot is fixed, and more cold-shell areal density adds both fuel energy
    and fusion; so for each (pressure, adiabat) the shell is made as thick as the budget allows.
    """
    adiabats = adiabats or [min_adiabat * f for f in (1.0, 1.25, 1.5, 2.0)]
    pressures = pressures or [max_pressure_gbar * f for f in (0.5, 0.65, 0.8, 0.9, 1.0)]
    best = None
    for adiabat in adiabats:
        for pressure in pressures:
            rho = icf_shell_density(pressure, adiabat)

            def target(rho_r: float, rho: float = rho, adiabat: float = adiabat) -> dict:
                return _icf_target({"rho_g_cc": rho, "rho_r": rho_r, "adiabat": adiabat})

            if target(1e-4)["hot"]["energy_J"] > fuel_budget_J:
                continue
            lo, hi = 1e-4, max_rho_r
            if target(hi)["hot"]["energy_J"] + target(hi)["cold_J"] <= fuel_budget_J:
                lo = hi
            else:
                for _ in range(60):
                    mid = 0.5 * (lo + hi)
                    t = target(mid)
                    lo, hi = (mid, hi) if t["hot"]["energy_J"] + t["cold_J"] <= fuel_budget_J else (lo, mid)
            t = target(lo)
            design = {"rho_g_cc": rho, "rho_r": lo, "adiabat": adiabat, "pressure_gbar": pressure,
                      "fusion_J": t["fusion_J"], "fuel_J": t["hot"]["energy_J"] + t["cold_J"],
                      "hot_spot_J": t["hot"]["energy_J"], "mass_mg": t["mass_kg"] * 1e6, "burnt": t["burnt"],
                      "hot_spot_um": t["hot"]["radius_m"] * 1e6}
            if best is None or design["fusion_J"] > best["fusion_J"]:
                best = design
    return best


# MagLIF on Z (2020 best shot, z3289): 13.5 kJ of hot-spot heat from ~20 MJ stored, 0.07%. With 15% of the stored
# energy counted as delivered (Settings.pulsed_power_efficiency), that is ~0.5% of the delivered energy.
LINER_COUPLING = 0.005


def _magnetised_liner(p: dict, s: Settings) -> Outcome:
    """Section 12's magnetised, tamped fuel, but as a column squeezed by a pulsed-power liner (as in MagLIF)."""
    ball = C.Ball(radius_m=p["radius_mm"] * 1e-3, fuel="DT", t_i_kev=p["t_i_kev"], fill=p["fill"],
                  b_tesla=p["b_tesla"], tamper="gold", field="axial", aspect=p["aspect"])
    r = C.burn(ball)
    drive = r["energy_in_J"] / LINER_COUPLING
    reqs = [Requirement("energy into the pinch (J)", drive), Requirement("compressed magnetic field (T)", p["b_tesla"])]
    return Outcome("magnetised_liner", p, r["fusion_J"], r["energy_in_J"], {"pulsed power": drive},
                   r["fusion_J"] / (Q_DT_MEV * MEV), reqs,
                   notes=f"{LINER_COUPLING:.1%} of the pulsed power reaches the fuel; {r['stopped']}; sections 12, 19")


def _muon_factory(p: dict, s: Settings) -> Outcome:
    """Exploding micrometre crystals make the pions that make muons (section 10), which then catalyse (section 9)."""
    f = ALT.crystal_muon_factory(p["e_max_gev_u"], s.muon_cost_gev, s.laser_to_ions)
    y, cooling, theory = _muon_cycle(p)
    laser = float(f["laser_GeV_per_muon"]) * 1e3 * MEV
    reqs = [Requirement("laser intensity (W/cm²)", float(f["strip_intensity_W_cm2"])),
            Requirement("liquid or solid D-T density (× liquid hydrogen)", p["phi"])] + theory
    drivers = {"laser": laser, "electricity": cooling} if cooling else {"laser": laser}
    return Outcome("muon_factory", p, y * Q_DT_MEV * MEV, laser * s.laser_to_ions, drivers, y, reqs,
                   notes=f"{float(f['laser_GeV_per_muon']):.2g} GeV of light per muon, {y:.0f} fusions each; "
                   "sections 9-10")


def _hollow_shell(p: dict, s: Settings) -> Outcome:
    """Section 8: ions of a hollow D-T shell imploding onto a focus of radius r_c. The focus returns at most
    0.5 pm / r_c of the ions' energy; ions springing back into the cold wall add the cold-catcher return."""
    focus = 0.5e-12 / (p["focus_nm"] * 1e-9)
    wall = float(EN.thick_target(p["ion_kev"], "D", "DT", {"T": 2.25e28})["gain"][0])
    fusion = focus + wall
    reqs = [Requirement("a microbubble implosion (simulated only)", 1.0)]
    return Outcome("hollow_shell", p, fusion, 1.0, {"laser": 1.0 / s.laser_to_ions}, fusion / (Q_DT_MEV * MEV),
                   reqs, notes="per joule of imploding ions; section 8")


LATTICE_RETURN = {"measured (NASA)": 4e-14, "perfect photon source": 3e-7}


def _lattice(p: dict, s: Settings) -> Outcome:
    """Section 9: NASA's lattice confinement fusion, per joule of electron beam."""
    fusion = LATTICE_RETURN[p["case"]]
    return Outcome("lattice_confinement", p, fusion, 1.0, {"accelerator": 1.0}, 0.0,
                   notes="gamma rays split deuterons in erbium deuteride; D-D fusion; section 9")


PB_AVALANCHE_K = 0.02  # secondary fusions per fusion, at most (section 4)


def _pb_catcher(p: dict, s: Settings) -> Outcome:
    """Section 4: fast protons stopping in cold boron or decaborane, with screening (section 5) and the avalanche."""
    material, n_b = ("boron", 1.3e29) if p["target"] == "boron" else ("decaborane", 4.56e28)
    gain = float(EN.thick_target(p["proton_kev"], "p", material, {"B11": n_b}, screening_ev=p["screening_ev"])
                 ["gain"][0])
    if p["avalanche"]:
        gain /= 1 - PB_AVALANCHE_K
    return Outcome("pB_catcher", p, gain, 1.0, _ion_source(p, s), 0.0, _screening(p),
                   notes="per joule of fast protons; sections 4-5")



# --- untested theoretical ideas ------------------------------------------------------------------

# Antiproton-initiated fusion (ICAN-II, AIMStar: Lewis, Gaidos, Smith et al., Penn State, ~1997-2001): antiprotons
# annihilate on uranium seeded into compressed D-T and the fission fragments heat a hot spot. Each annihilation
# fissions a nucleus whose fragments leave ~0.17 GeV locally; the annihilation pions mostly escape.
ANTIPROTON_DEPOSIT_J = 0.17e3 * MEV
# Electricity per antiproton: CERN's AD spends ~2 mJ of 26 GeV proton beam per delivered antiproton (5e5 protons
# each); at 40% beam efficiency that is ~5 mJ. Forward's engineered factory: 2 m c^2 / 1e-4 = ~3 uJ.
ANTIPROTON_ELECTRICITY_J = {"CERN today": 2.1e-3 / 0.4, "optimised factory (Forward)": 1.88e3 * MEV / 1e-4}


def _antimatter(p: dict, s: Settings) -> Outcome:
    """Compressed D-T (as compressed_icf) whose hot spot is heated by antiproton-induced fission."""
    t = _icf_target(p)
    antiprotons = t["hot"]["energy_J"] / ANTIPROTON_DEPOSIT_J
    beam = antiprotons * ANTIPROTON_ELECTRICITY_J[p["antiprotons"]] * s.accelerator_efficiency
    laser = t["cold_J"] / s.implosion_efficiency
    reqs = t["requirements"] + [Requirement("laser energy per shot (J)", laser),
                                Requirement("antiprotons per shot", antiprotons),
                                Requirement("antiprotons light a hot spot by fission (ICAN, AIMStar)", 1.0)]
    if p["antiprotons"] != "CERN today":
        reqs.append(Requirement("an antiproton factory 10,000 times more efficient than CERN's", 1e-4))
    return Outcome("antimatter_spark", p, t["fusion_J"], t["hot"]["energy_J"] + t["cold_J"],
                   {"laser": laser, "accelerator": beam}, t["fusion_J"] / (Q_DT_MEV * MEV), reqs,
                   notes=f"{antiprotons:.2g} antiprotons per shot; section 22")


# Dynamically assisted tunnelling (Queisser & Schuetzhold, PRC 100, 041601 (2019); Wang, PRC 102, 011601 (2020)):
# an intense light field raises the D-T tunnelling probability. Wang's numbers for 800 nm light at 5e21 W/cm2: over
# 1e9 at 1 keV, over 100 at 5 keV, over 10 at 10 keV (taken as equal here, generous); a 2024 calculation (Thomson et
# al., PRC 110, 034614) finds only +7-70% even at 1e27-1e29 W/cm2.
ASSISTED_INTENSITY_W_CM2 = 5e21
ASSISTED_BOOST = {"Wang 2020 (optimistic)": ((1.0, 5.0, 10.0, 30.0), (9.0, 2.0, 1.0, 0.0)),
                  "Thomson et al. 2024 (pessimistic)": ((1.0, 30.0), (math.log10(1.7), 0.0))}


@lru_cache(maxsize=256)
def _assisted_reactivity(t_kev: float, model: str) -> float:
    """Maxwellian D-T <sigma v> (m^3/s) with the cross section multiplied by the light field's boost."""
    energies, logs = ASSISTED_BOOST[model]
    mu = P.reduced_mass_u(P.MASS_U["D"], P.MASS_U["T"]) * P.AMU
    e = np.geomspace(1e-2, max(60.0 * t_kev, 50.0), 4000)
    boost = 10 ** np.interp(np.log10(e), np.log10(energies), logs, left=logs[0], right=0.0)
    integrand = EN.sigma_m2("D", "T", e) * boost * e * np.exp(-e / t_kev)
    return float(math.sqrt(8 / (math.pi * mu)) * (t_kev * KEV) ** -1.5 * np.trapezoid(integrand, e) * KEV**2)


def _laser_assisted(p: dict, s: Settings) -> Outcome:
    """A cool D-T ball (0.5-5 keV: the point is fusion without full heating) bathed in intense light for as long
    as it holds together, its fusion boosted by assisted tunnelling. The light passes through (none of it is
    counted as absorbed), so the laser must supply intensity x cross-section x time; heating the ball to T costs
    its heat / laser_to_ions, and squeezing it to rho costs its electrons' degenerate energy / implosion_efficiency."""
    rho, t_kev, radius = p["rho_g_cc"], p["t_kev"], p["radius_mm"] * 1e-3
    n = rho * 1e3 / (0.5 * (P.MASS_U["D"] + P.MASS_U["T"]) * P.AMU)
    volume = 4 / 3 * math.pi * radius**3
    sound = math.sqrt(5 / 3 * 2 * n * t_kev * KEV / (rho * 1e3))
    hold = radius / (4 * sound)
    raw = (n / 2) ** 2 * _assisted_reactivity(t_kev, p["model"]) * Q_DT_MEV * MEV * volume * hold
    fusion = deplete(raw, n / 2 * volume * Q_DT_MEV * MEV)
    heat = 3 * n * t_kev * KEV * volume
    squeeze = n * D.squeeze_energy_kev(n) * KEV * volume  # cold squeeze to rho, by implosion (1% of the light)
    light = ASSISTED_INTENSITY_W_CM2 * 1e4 * math.pi * radius**2 * hold
    laser = light + heat / s.laser_to_ions + squeeze / s.implosion_efficiency
    reqs = [Requirement("cold squeezed density (g/cm³)", rho), Requirement("radial squeeze of the fuel",
                                                                          max(1.0, (rho / 0.25) ** (1 / 3))),
            Requirement("laser energy per shot (J)", laser), Requirement("laser peak power (W)", light / hold),
            Requirement("laser intensity (W/cm²)", ASSISTED_INTENSITY_W_CM2),
            Requirement("light-assisted tunnelling (Queisser & Schützhold)", 1.0)]
    return Outcome("laser_assisted", p, fusion, heat, {"laser": laser}, fusion / (Q_DT_MEV * MEV), reqs,
                   notes=f"light held {hold * 1e12:.2g} ps; section 22")


ROUTES = {r.name: r for r in [
    Route("cold_catcher", "1, 5", "Fast deuterons stopping in a cold catcher", [
        Param("ion_kev", 10.0, 5000.0, default=300.0, help="deuteron energy (keV)"),
        Param("catcher", choices=("tritide", "LiD"), default="tritide"),
        Param("source", choices=("crystals", "accelerator"), default="crystals"),
        Param("screening_ev", 0.0, 3000.0, log=False, default=0.0, help="electron screening potential (eV)")],
        _cold_catcher),
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
        Param("sticking", choices=tuple(STICKING), default="PSI 1987 (0.45%)"), RF_STRIPPING], _muon),
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
        Param("fuel_um", 0.1, 100.0, default=5.0, help="fuel radius (um)"), *PB_THEORY], _son_fisch, dt=False),
    Route("opaque_pB", "18", "Proton-boron ball that traps its own X-rays", [
        Param("rho_g_cc", 100.0, 1e5, default=1e4, help="fuel density (g/cm3)"),
        Param("t_i_kev", 150.0, 500.0, default=200.0, help="ion temperature (keV)"),
        Param("radius_um", 0.5, 1e4, default=20.0, help="ball radius (um)"), *PB_THEORY], _opaque_pb, dt=False),
    Route("fibre_pinch", "19", "D-T fibre pinched by its own current", [
        Param("current_ma", 0.3, 100.0, default=26.0, help="current (MA)"),
        Param("t_kev", 3.0, 30.0, default=10.0, help="Bennett temperature (keV)"),
        Param("alfven_times", 1.0, 1e4, default=1.0, help="how long the pinch stays stable"),
        Param("squeeze", 1.0, 300.0, default=6.0, help="radial squeeze of the fibre"), SPIN], _pinch),
    Route("dynamic_pinch", "24", "D-T fibre pinch followed in time, driven by a generator circuit", [
        Param("current_ma", 0.5, 100.0, default=26.0, help="peak current into a fixed load (MA)"),
        Param("radius_um", 20.0, 2000.0, default=240.0, help="fibre radius (um)"),
        Param("preheat_kev", 0.01, 0.5, default=0.1, help="laser preheat of the fibre (keV)"),
        Param("alfven_times", 1.0, 1e4, default=1.0, help="e-folds of the m = 0 mode before the pinch breaks"),
        Param("stability", choices=tuple(PINCH_STABILITY), default="breaks",
              help="untested ways to stay stable longer")], _dynamic_pinch, effort=0.25),
    Route("channelling", "20", "Ions channelled through a tritide crystal", [
        Param("ion_kev", 50.0, 3000.0, default=300.0),
        Param("source", choices=("aligned beam", "exploding crystal"), default="aligned beam"),
        Param("values", choices=("measured", "generous"), default="measured")], _channelling),
    Route("compressed_icf", "7, 11", "Ordinary inertial fusion: compressed D-T lit by a hot spot (the reference)", [
        Param("rho_g_cc", 50.0, 5000.0, default=950.0, help="density of the cold shell (g/cm3)"),
        Param("rho_r", 0.2, 10.0, default=0.6, help="areal density of the cold shell (g/cm2)"),
        Param("adiabat", 1.0, 5.0, log=False, default=3.0, help="cold fuel pressure / Fermi pressure"),
        Param("ignitor", choices=("hot spot", "crystal spark"), default="hot spot"),
        Param("driver", choices=("laser", "pulsed-power X-rays"), default="laser"), SPIN], _compressed_icf),
    Route("magnetised_liner", "12, 19", "Magnetised fuel column squeezed by a pulsed-power liner (MagLIF-like)", [
        Param("radius_mm", 0.02, 3.0, default=0.1, help="fuel radius at stagnation (mm)"),
        Param("aspect", 2.0, 100.0, default=30.0, help="column length / radius"),
        Param("t_i_kev", 3.0, 40.0, default=10.0, help="starting temperature (keV)"),
        Param("fill", 0.01, 1.0, default=1.0, help="fuel density / solid"),
        Param("b_tesla", 100.0, 1e5, default=1e4, help="compressed axial field (T)"), SPIN], _magnetised_liner),
    Route("muon_factory", "9, 10", "Exploding crystals make the muons for muon catalysis", [
        Param("e_max_gev_u", 0.5, 5.0, default=2.0, help="top ion energy of the explosion (GeV per nucleon)"),
        Param("phi", 0.1, 10.0, default=1.2, help="fuel density (x liquid hydrogen)"),
        Param("sticking", choices=tuple(STICKING), default="PSI 1987 (0.45%)"), RF_STRIPPING], _muon_factory),
    Route("hollow_shell", "8", "Hollow shell imploding onto its centre (microbubble implosion)", [
        Param("ion_kev", 10.0, 1000.0, default=100.0),
        Param("focus_nm", 0.01, 10.0, default=0.5, help="radius of the focus (nm)")], _hollow_shell),
    Route("lattice_confinement", "9", "NASA's lattice confinement fusion", [
        Param("case", choices=tuple(LATTICE_RETURN), default="measured (NASA)")], _lattice),
    Route("pB_catcher", "4, 5", "Fast protons stopping in cold boron fuel", [
        Param("proton_kev", 100.0, 10000.0, default=700.0),
        Param("target", choices=("boron", "decaborane"), default="boron"),
        Param("avalanche", choices=(False, True), default=True),
        Param("source", choices=("crystals", "accelerator"), default="crystals"),
        Param("screening_ev", 0.0, 3000.0, log=False, default=0.0, help="electron screening potential (eV)")],
        _pb_catcher),
    Route("antimatter_spark", "22", "Compressed D-T lit by antiproton-induced fission (ICAN, AIMStar)", [
        Param("rho_g_cc", 50.0, 5000.0, default=950.0, help="density of the cold shell (g/cm3)"),
        Param("rho_r", 0.2, 10.0, default=0.6, help="areal density of the cold shell (g/cm2)"),
        Param("adiabat", 1.0, 5.0, log=False, default=3.0, help="cold fuel pressure / Fermi pressure"),
        Param("antiprotons", choices=tuple(ANTIPROTON_ELECTRICITY_J), default="CERN today")], _antimatter),
    Route("laser_assisted", "22", "D-T fusion boosted by light-assisted tunnelling", [
        Param("rho_g_cc", 0.25, 1000.0, default=0.25, help="fuel density (g/cm3)"),
        Param("t_kev", 0.5, 5.0, default=2.0, help="fuel temperature (keV)"),
        Param("radius_mm", 0.01, 10.0, default=0.5, help="ball radius (mm)"),
        Param("model", choices=tuple(ASSISTED_BOOST), default="Wang 2020 (optimistic)")], _laser_assisted),
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
           budgets=BUDGETS, fixed: dict | None = None, start=None, theory: bool = False) -> dict:
    """Best electric gain of each route within each stretch budget.

    Random sampling of every parameter (log-uniform for ranges), then a
    hill-climb from the best design within each budget. `fixed` holds
    parameters at given values (for routes that have them); `start` adds known designs (parameter dicts) to
    the sample. With theory=False (the default) untested theoretical ideas are held off and designs that rely
    on one are left out; theory=True lets them in (their stretch counts only the technology they need).
    Returns {route: {budget: Outcome or None}}.
    """
    settings = settings or Settings()
    rng = np.random.default_rng(seed)
    results = {}
    with np.errstate(all="ignore"):
        for name in routes or ROUTES:
            route = ROUTES[name]
            names = {q.name for q in route.params}
            held = {k: v for k, v in (fixed or {}).items() if k in names}
            if not theory:
                held.update({k: v for k, v in THEORY_OFF.items() if k in names and k not in held})
            free = [q for q in route.params if q.name not in held]
            n_samples, n_refine = max(1, int(samples * route.effort)), max(1, int(refine * route.effort))
            pool = [o for o in (_safe(name, {**{q.name: q.sample(rng) for q in free}, **held}, settings)
                                for _ in range(n_samples)) if o is not None]
            for seed_design in [dict(held)] + [dict(d) for d in (start or [])]:
                seed_design.update(held)
                known = {q.name for q in route.params}
                out = _safe(name, {k: v for k, v in seed_design.items() if k in known}, settings)
                if out is not None:
                    pool.append(out)
            best_by_budget = {}
            for budget in budgets:
                inside = [o for o in pool if o.stretch <= budget and (theory or not o.untested)]
                if not inside:
                    best_by_budget[budget] = None
                    continue
                best = max(inside, key=lambda o: o.electric_gain)
                for step in range(n_refine):
                    scale = 0.15 * (1 - step / n_refine) + 0.01
                    trial = {**{q.name: q.perturb(best.params[q.name], rng, scale) for q in free}, **held}
                    out = _safe(name, trial, settings)
                    if (out is not None and out.stretch <= budget and (theory or not out.untested)
                            and out.electric_gain > best.electric_gain):
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
        if r.untested:
            print(f"  relies on untested theory: {r.name}")
            continue
        best = f"{_fmt(r.best)} ({TODAY[r.name][1]})"
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
    find.add_argument("--theory", action="store_true", help="let in untested theoretical ideas")
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
    results = search(routes, args.samples, args.refine, args.seed, settings, theory=args.theory)
    for budget in BUDGETS:
        print(f"\n== Best electric gain, {budget_label(budget)} ==")
        ranked = sorted(((name, by[budget]) for name, by in results.items() if by[budget] is not None),
                        key=lambda t: -t[1].electric_gain)
        for name, o in ranked:
            print(f"  {name:16s} electric gain {_fmt(o.electric_gain):>9s}  stretch {_fmt(o.stretch):>8s}"
                  + (f"  limited by {o.limiting}" if o.limiting else "")
                  + (f"  [untested: {'; '.join(o.untested)}]" if o.untested else ""))
    if args.output:
        with open(args.output, "w") as f:
            json.dump({name: {budget_label(b): (o.summary() if o else None) for b, o in by.items()}
                       for name, by in results.items()}, f, indent=1, default=float)
    return 0


if __name__ == "__main__":
    sys.exit(main())
