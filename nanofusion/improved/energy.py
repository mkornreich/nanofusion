"""Energy balance: can laser-exploded nanocrystals give more fusion energy than goes in?

Three levels of optimism:

1. Cold catcher, perfect laser: every joule becomes ion kinetic energy, and the
   ions stop in cold material. Gain = Q x reactions per ion / ion energy.
2. Hot catcher: the ions slow down in a plasma that is already hot, where
   electrons drag far less (Spitzer/Stix slowing). This is the "beam into hot
   plasma" energy multiplication of neutral-beam-driven fusion.
3. Crystal-array burn: a ball of crystals is exploded at once. The ion energy
   plus the electron heating form a plasma that burns until it flies apart
   (0-D model with bremsstrahlung). Gain = fusion energy / energy put in.

Formulas
--------
* Maxwellian reactivity: <sigma v> = sqrt(8/(pi mu)) (kT)^-3/2 Int sigma(E) E exp(-E/kT) dE.
* Slowing of a fast ion (Stix 1972; Heidbrink, Nucl. Fusion review):
  dE/dt = -(2E/tau_se) [1 + (E_c/E)^(3/2)],
  tau_se = 6.27e8 A T_e[eV]^(3/2) / (Z^2 n_e[cm^-3] lnL) s,
  E_c = 14.8 A T_e (sum_j n_j Z_j^2 / (n_e A_j))^(2/3).
* Bremsstrahlung: P = 5.35e-37 n_e sum_j n_j Z_j^2 sqrt(T_e[keV]) W/m^3 (NRL formulary).
* Disassembly (burn) time of a sphere of radius L: tau = L / (4 c_s), the
  usual inertial-confinement estimate, with c_s = sqrt((n_e T_e + sum n_j T_i) / rho).
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np

from . import physics as P

KEV = P.J_PER_KEV
# Spin polarisation of D-T fuel: parallel spins raise the D-T cross section by 1.5 (theory, e.g. Kulsrud et al.,
# PRL 49, 1248 (1982)). Experiments only; leave at 1.0 for unpolarised fuel.
DT_SPIN_FACTOR = 1.0
# Above 550 keV (centre of mass) the D-T cross section is held constant by default, which overstates it.
# True switches to the measured fall-off (Bosch & Hale's 550-4700 keV fit). It matters only for MeV ions.
DT_HIGH_ENERGY = False
SPECIES = {  # mass (u), charge
    "p": (P.MASS_U["p"], 1),
    "He4": (P.MASS_U["He4"], 2),
    "B11": (P.MASS_U["B11"], 5),
    "D": (P.MASS_U["D"], 1),
    "T": (P.MASS_U["T"], 1),
    "Li6": (P.MASS_U["Li6"], 3),
    "Li7": (P.MASS_U["Li7"], 3),
}


def sigma_m2(a: str, b: str, e_cm_kev):
    """Fusion cross section (m^2) for ion species a+b at CM energy (keV). 0 if not modelled."""
    pair = {a, b}
    e = np.asarray(e_cm_kev, dtype=float)
    if pair == {"p", "B11"}:
        return P.sigma_pb11_mb(e) * P.MB_TO_M2
    if pair == {"D"}:
        return P.sigma_dd_mb(e) * P.MB_TO_M2
    if pair == {"D", "T"}:
        return P.sigma_dt_mb(e, high_energy=DT_HIGH_ENERGY) * P.MB_TO_M2 * DT_SPIN_FACTOR
    for li in ("Li6", "Li7"):
        if pair == {"D", li}:
            mu = P.reduced_mass_u(P.MASS_U["D"], P.MASS_U[li])
            return P.sigma_dli_mb(e, P.S_DLI_KEV_MB[li], P.gamow_constant(1, 3, mu)) * P.MB_TO_M2
    return np.zeros_like(e)


def q_value(a: str, b: str) -> float:
    pair = {a, b}
    if pair == {"p", "B11"}:
        return P.Q_VALUE_MEV["pB"]
    if pair == {"D"}:
        return P.Q_VALUE_MEV["DD"]
    if pair == {"D", "T"}:
        return P.Q_VALUE_MEV["DT"]
    for li in ("Li6", "Li7"):
        if pair == {"D", li}:
            return P.Q_VALUE_MEV[li]
    return 0.0


def reactivity(a: str, b: str, t_kev: float) -> float:
    """Maxwellian <sigma v> (m^3/s) at temperature t_kev."""
    mu = P.reduced_mass_u(SPECIES[a][0], SPECIES[b][0]) * P.AMU
    e = np.geomspace(1e-3, max(60.0 * t_kev, 50.0), 4000)  # keV
    integrand = sigma_m2(a, b, e) * e * np.exp(-e / t_kev)
    integral = np.trapezoid(integrand, e) * KEV**2  # m^2 J^2
    return float(math.sqrt(8 / (math.pi * mu)) * (t_kev * KEV) ** -1.5 * integral)


@dataclass
class Plasma:
    """Fully ionised plasma: ion densities (m^-3) by species, electron/ion temperatures (keV)."""

    ions: dict = field(default_factory=dict)
    t_e: float = 10.0
    t_i: float = 10.0

    @property
    def n_e(self) -> float:
        return sum(n * SPECIES[s][1] for s, n in self.ions.items())

    @property
    def mass_density(self) -> float:
        return sum(n * SPECIES[s][0] * P.AMU for s, n in self.ions.items())

    def coulomb_log(self) -> float:
        n_cm3 = self.n_e * 1e-6
        t_ev = self.t_e * 1e3
        return max(2.0, 24.0 - math.log(math.sqrt(n_cm3) / t_ev))

    def slowing_time(self, projectile: str) -> float:
        a, z = SPECIES[projectile]
        return 6.27e8 * a * (self.t_e * 1e3) ** 1.5 / (z**2 * self.n_e * 1e-6 * self.coulomb_log())

    def critical_energy(self, projectile: str) -> float:
        a = SPECIES[projectile][0]
        s = sum(n * SPECIES[sp][1] ** 2 / SPECIES[sp][0] for sp, n in self.ions.items()) / self.n_e
        return 14.8 * a * self.t_e * s ** (2 / 3)

    def sound_speed(self) -> float:
        pressure = (self.n_e * self.t_e + sum(self.ions.values()) * self.t_i) * KEV
        return math.sqrt(pressure / self.mass_density)

    def bremsstrahlung(self) -> float:
        """W/m^3."""
        zsum = sum(n * SPECIES[s][1] ** 2 for s, n in self.ions.items())
        return 5.35e-37 * self.n_e * zsum * math.sqrt(self.t_e)

    def fusion_power(self) -> tuple[float, dict]:
        """Thermonuclear fusion power density (W/m^3) and its channels."""
        names = list(self.ions)
        channels = {}
        for i, a in enumerate(names):
            for b in names[i:]:
                q = q_value(a, b)
                if q == 0:
                    continue
                rate = self.ions[a] * self.ions[b] * reactivity(a, b, self.t_i)
                if a == b:
                    rate /= 2
                channels[f"{a}+{b}"] = rate * q * 1.602176634e-13
        return sum(channels.values()), channels


def beam_in_plasma(e0_kev: float, projectile: str, plasma: Plasma, max_path_m: float = math.inf,
                   points: int = 3000) -> dict:
    """Fusion by one fast ion slowing down in a hot plasma.

    Returns reactions by channel, fusion energy (MeV), the ion's energy
    multiplication (fusion energy / its starting energy), and the path length
    it needs to thermalise. ``max_path_m`` stops the ion when it leaves a
    plasma of that size. Target ions are treated as at rest (beam >> thermal).
    """
    m_p = SPECIES[projectile][0]
    tau = plasma.slowing_time(projectile)
    e_c = plasma.critical_energy(projectile)
    e_end = max(1.5 * plasma.t_i, 1e-3)
    if e0_kev <= e_end:
        return {"channels": {}, "fusion_MeV": 0.0, "gain": 0.0, "path_m": 0.0}
    e = np.geomspace(e0_kev, e_end, points)
    v = np.sqrt(2 * e * KEV / (m_p * P.AMU))
    dedt = (2 * e / tau) * (1 + (e_c / e) ** 1.5)  # keV/s, magnitude
    dt = np.concatenate([[0.0], np.cumsum(0.5 * (1 / dedt[1:] + 1 / dedt[:-1]) * -np.diff(e))])
    path = np.concatenate([[0.0], np.cumsum(0.5 * (v[1:] + v[:-1]) * np.diff(dt))])
    inside = path <= max_path_m
    channels = {}
    for target, n_t in plasma.ions.items():
        q = q_value(projectile, target)
        if q == 0:
            continue
        m_t = SPECIES[target][0]
        e_cm = e * m_t / (m_p + m_t)
        rate = n_t * sigma_m2(projectile, target, e_cm) * v  # 1/s
        p = np.concatenate([[0.0], np.cumsum(0.5 * (rate[1:] + rate[:-1]) * np.diff(dt))])
        channels[f"{projectile}+{target}"] = float(p[inside][-1])
    fusion = sum(channels[k] * q_value(*k.split("+")) for k in channels)
    return {"channels": channels, "fusion_MeV": fusion, "gain": fusion / (e0_kev / 1e3),
            "path_m": float(path[-1]), "slowing_time_s": float(dt[-1])}


# --- fuels and crystal arrays -------------------------------------------------

FUELS = {
    # name: (ion fractions per formula unit, formula units per m^3 in the solid)
    "6LiD": ({"Li6": 1, "D": 1}, P.lid_molecule_density()),
    "7LiD": ({"Li7": 1, "D": 1}, P.lid_molecule_density()),
    "6LiDT": ({"Li6": 1, "D": 0.5, "T": 0.5}, P.lid_molecule_density()),  # lithium deuteride-tritide
    "DT": ({"D": 1, "T": 1}, 2.25e28),  # solid D-T (about 0.2 g/cm^3): molecules per m^3
    "B10H14": ({"B11": 10, "p": 14}, 4.56e27),  # decaborane-11 at 0.94 g/cm^3
    "D2": ({"D": 2}, 3.0e28),  # solid deuterium (~0.2 g/cm^3) as in laser-cluster experiments
}


def fuel_plasma(fuel: str, fill: float, t_e: float, t_i: float) -> Plasma:
    parts, n_units = FUELS[fuel]
    return Plasma({s: f * n_units * fill for s, f in parts.items()}, t_e=t_e, t_i=t_i)


def charge_density(fuel: str) -> float:
    """Charge density (C/m^3) of the fully ionised solid fuel."""
    parts, n_units = FUELS[fuel]
    return n_units * sum(f * SPECIES[s][1] for s, f in parts.items()) * P.E_CHARGE


def explosion_energy_per_ion_kev(fuel: str, radius_m: float) -> dict:
    """Mean kinetic energy per ion species after a uniform-sphere Coulomb explosion.

    E(r) = z e rho r^2 / (3 eps0), averaged over a uniform sphere: 3/5 of the
    surface value. The MD runs reach about 90% of this; we use 0.9.
    """
    rho = charge_density(fuel)
    out = {}
    for s in FUELS[fuel][0]:
        z = SPECIES[s][1]
        out[s] = 0.9 * 0.6 * z * P.E_CHARGE * rho * radius_m**2 / (3 * P.EPS0) / KEV
    return out


def array_burn(fuel: str, crystal_radius_nm: float, fill: float, ball_radius_m: float,
               t_e_heat_kev: float = 1.0, fast_in_final_plasma: bool = True) -> dict:
    """A ball of crystals (filling a fraction ``fill`` of its volume) exploded at once.

    Energy in: the Coulomb explosion energy of every crystal (which the laser
    must at least supply) plus heating the electrons to t_e_heat_kev.

    Phase 1, fast ions: every ion slows down in the plasma, fusing along the
    way, and escapes if its path is longer than the ball. With
    fast_in_final_plasma=True the plasma is taken to be at its final
    temperature at once (generous). With False, the ions slow down on the
    initial t_e_heat_kev electrons, where they lose energy much faster
    (conservative).
    Phase 2, thermal burn: the deposited energy is shared by all particles
    (temperature T_eq). The ball burns for tau = L/(4 c_s), with bremsstrahlung
    losses. No alpha self-heating, which slightly understates a burning case.
    """
    parts, n_units = FUELS[fuel]
    n_u = n_units * fill
    volume = 4 / 3 * math.pi * ball_radius_m**3
    e_ion = explosion_energy_per_ion_kev(fuel, crystal_radius_nm * 1e-9)
    n_ions = {s: f * n_u for s, f in parts.items()}
    n_e = sum(n * SPECIES[s][1] for s, n in n_ions.items())
    energy_in = (sum(n_ions[s] * e_ion[s] for s in n_ions) + 1.5 * n_e * t_e_heat_kev) * KEV * volume  # J
    n_particles = n_e + sum(n_ions.values())
    t_eq = energy_in / volume / (1.5 * n_particles * KEV)
    plasma = Plasma(n_ions, t_e=t_eq, t_i=t_eq)
    # Phase 1: fast ions
    fast = 0.0
    slowing_plasma = plasma if fast_in_final_plasma else Plasma(n_ions, t_e=t_e_heat_kev, t_i=t_e_heat_kev)
    for s, n in n_ions.items():
        res = beam_in_plasma(e_ion[s], s, slowing_plasma, max_path_m=ball_radius_m)
        fast += n * volume * res["fusion_MeV"] * 1.602176634e-13
    # Phase 2: thermal burn
    p_fus, channels = plasma.fusion_power()
    p_brem = plasma.bremsstrahlung()
    tau = ball_radius_m / (4 * plasma.sound_speed())
    thermal = p_fus * volume * tau
    brem = p_brem * volume * tau
    fusion = fast + thermal
    return {
        "fuel": fuel, "crystal_radius_nm": crystal_radius_nm, "fill": fill, "ball_radius_m": ball_radius_m,
        "ion_energy_keV": e_ion, "T_eq_keV": t_eq, "burn_time_s": tau,
        "energy_in_J": energy_in, "fusion_fast_J": fast, "fusion_thermal_J": thermal,
        "bremsstrahlung_J": brem, "fusion_to_brem": p_fus / p_brem if p_brem else math.inf,
        "gain": fusion / energy_in, "gain_after_radiation": (fusion - brem) / energy_in,
        "areal_density_g_cm2": plasma.mass_density * ball_radius_m * 0.1,
        "mass_g": plasma.mass_density * volume * 1e3,
        "stripping_intensity_W_cm2": P.intensity_w_cm2(charge_density(fuel) * crystal_radius_nm * 1e-9
                                                       / (3 * P.EPS0)),
    }


# --- generic cold-target yields and the proton-boron "avalanche" ------------------

def _channel_sigma(projectile: str, target: str):
    return lambda e_cm: sigma_m2(projectile, target, e_cm)


def thick_target(e0_kev, projectile: str, material: str, targets: dict, points: int = 2000,
                 screening_ev: float = 0.0) -> dict:
    """Reactions per ion of lab energy e0 stopping in a cold material.

    targets: {target species: number density (m^-3)} to react with.
    screening_ev: electron-screening potential U_e. Bound electrons lower the
    Coulomb barrier by about U_e, so the cross section is taken at E + U_e
    (Assenbaum et al.). That equals the usual exp(pi eta U_e / E) factor when
    U_e << E, and stays finite at low E.
    Returns {"reactions": array, "fusion_MeV": array, "gain": array}.
    """
    e0 = np.atleast_1d(np.asarray(e0_kev, dtype=float))
    m_p, z_p = SPECIES[projectile]
    grid = np.concatenate([[0.0], np.geomspace(0.05, max(float(e0.max()), 1.0), points)])
    stop = P.stopping_kev_per_m(np.maximum(grid, 1e-3), z_p, m_p, material)
    total = np.zeros_like(grid)
    energy = np.zeros_like(grid)
    for target, n_t in targets.items():
        q = q_value(projectile, target)
        if q == 0:
            continue
        m_t = SPECIES[target][0]
        e_cm = grid * m_t / (m_p + m_t)
        sigma = sigma_m2(projectile, target, e_cm + screening_ev * 1e-3)
        dydE = n_t * sigma / stop
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (dydE[1:] + dydE[:-1]) * np.diff(grid))])
        total += cum
        energy += cum * q
    y = np.interp(e0, grid, total)
    mev = np.interp(e0, grid, energy)
    return {"reactions": y, "fusion_MeV": mev, "gain": mev / (e0 / 1e3)}


def pb11_avalanche(material: str = "decaborane", alpha_kev: float = 2900.0, nuclear_boost: float = 1.0,
                   hot_plasma: "Plasma | None" = None, points: int = 600) -> dict:
    """Secondary p+11B fusions triggered per fusion, via alpha -> proton knock-on.

    Each fusion gives 3 alphas (about 2.9 MeV each). While an alpha slows down
    it hits protons elastically (Rutherford: dsigma/dT = 2 pi (z_a z_p k e^2)^2 /
    (m_p v_a^2 T^2), T <= 0.64 E_a). Each knocked-on proton then fuses with
    probability Y_p(T) before stopping. ``nuclear_boost`` multiplies the elastic
    cross section to allow for nuclear (non-Rutherford) scattering, which is
    stronger at MeV energies. An avalanche needs k >= 1.
    """
    atoms, _ = P._material(material) if hot_plasma is None else ({}, 0.0)
    if hot_plasma is None:
        n_p = atoms.get(1, 0.0)
        n_b = atoms.get(5, 0.0)
    else:
        n_p = hot_plasma.ions.get("p", 0.0)
        n_b = hot_plasma.ions.get("B11", 0.0)
    m_a, m_pr = SPECIES["He4"][0] * P.AMU, SPECIES["p"][0] * P.AMU
    e_a = np.geomspace(alpha_kev, 20.0, points)  # alpha energy along its path, keV
    if hot_plasma is None:
        stop_a = P.stopping_kev_per_m(e_a, 2, SPECIES["He4"][0], material)  # keV/m
    else:
        tau = hot_plasma.slowing_time("He4")
        e_c = hot_plasma.critical_energy("He4")
        v_tmp = np.sqrt(2 * e_a * KEV / m_a)
        stop_a = (2 * e_a / tau) * (1 + (e_c / e_a) ** 1.5) / v_tmp
    # proton yields as a function of knock-on energy
    t_grid = np.geomspace(1.0, 0.64 * alpha_kev, 400)
    if hot_plasma is None:
        y_p = thick_target(t_grid, "p", material, {"B11": n_b})["reactions"]
    else:
        y_p = np.array([beam_in_plasma(t, "p", hot_plasma)["channels"].get("p+B11", 0.0) for t in t_grid])
    k_e2 = P.K_COULOMB * P.E_CHARGE**2  # J m
    secondary = 0.0
    for i in range(len(e_a) - 1):
        e_mid = 0.5 * (e_a[i] + e_a[i + 1])
        ds = (e_a[i] - e_a[i + 1]) / np.interp(e_mid, e_a[::-1], stop_a[::-1])  # m
        v2 = 2 * e_mid * KEV / m_a
        k_r = 2 * math.pi * (2 * 1 * k_e2) ** 2 / (m_pr * v2)  # m^2 J
        t_max = 0.64 * e_mid
        t = t_grid[t_grid <= t_max]
        if len(t) < 2:
            continue
        integrand = k_r / (t * KEV) ** 2 * np.interp(t, t_grid, y_p) * KEV  # per keV of T
        secondary += n_p * ds * nuclear_boost * np.trapezoid(integrand, t)
    k = 3 * secondary
    return {"k": k, "per_alpha": secondary, "material": material, "nuclear_boost": nuclear_boost}


# --- chain reaction in D-T carried by the fusion products ----------------------------

N_ELASTIC_14MEV_B = {"D": 0.639, "T": 0.933}  # neutron elastic cross sections at 14.1 MeV (ENDF/B-VIII.0), barn
# Elastic scattering at 14 MeV is peaked forwards: mean cosine of the scattering angle in the centre-of-mass
# frame (ENDF/B-VIII.0). The struck nucleus then gets less energy on average: 4.2 MeV for D, not 6.3.
N_MEAN_COSINE_14MEV = {"D": 0.330, "T": 0.464}
E_NEUTRON_KEV = 14070.0
E_ALPHA_KEV = 3520.0
NUCLEAR_KNOCK_ON_KEV = 500.0


def dt_chain(hot_plasma: "Plasma | None" = None, material: str = "DT", nuclear_boost: float = 1.0,
             first_collision_probability: float | None = None, neutrons: int = 20000, seed: int = 3) -> dict:
    """Secondary D-T fusions per D-T fusion, carried by its 14.1 MeV neutron and its 3.5 MeV alpha.

    Both products knock deuterons and tritons on. A knock-on of energy T then
    fuses with probability Y(T) before it stops: in cold D-T (thick-target
    yield) or in a hot plasma (``beam_in_plasma``). A chain reaction needs
    k >= 1.

    * Neutron: each elastic collision gives the struck nucleus up to
      T_max = 8/9 of the neutron's energy for D and 3/4 for T. Scattering is
      peaked forwards (distribution exp(b cos theta), with the measured mean
      cosine at 14 MeV, falling in proportion to the neutron's energy), so
      the mean share is less than half of T_max. With
      first_collision_probability=None the neutron is followed down to
      20 keV in an infinite medium (the best case). A number counts only the
      first collision, weighted by that probability, for a body thinner than
      the 14 MeV mean free path (about 4.7 g/cm^2).
    * Alpha: Rutherford knock-ons on D and T while it slows down, as in
      ``pb11_avalanche``. ``nuclear_boost`` multiplies the cross section for
      knock-ons above 0.5 MeV, where nuclear scattering outweighs Coulomb
      scattering; softer collisions are the alpha's ordinary heating of ions.
    """
    if hot_plasma is None:
        n_h = P._material(material)[0][1]
        n_d = n_t = 0.5 * n_h
    else:
        n_d, n_t = hot_plasma.ions["D"], hot_plasma.ions["T"]
    t_grid = np.geomspace(1.0, 0.9 * E_NEUTRON_KEV, 260)
    if hot_plasma is None:
        y_d = thick_target(t_grid, "D", material, {"T": n_t})["reactions"]
        y_t = thick_target(t_grid, "T", material, {"D": n_d})["reactions"]
    else:
        y_d = np.array([beam_in_plasma(t, "D", hot_plasma)["channels"].get("D+T", 0.0) for t in t_grid])
        y_t = np.array([beam_in_plasma(t, "T", hot_plasma)["channels"].get("T+D", 0.0) for t in t_grid])
    yields = {"D": y_d, "T": y_t}
    density = {"D": n_d, "T": n_t}

    # alpha knock-ons
    m_a = SPECIES["He4"][0] * P.AMU
    e_a = np.geomspace(E_ALPHA_KEV, 20.0, 500)
    v2 = 2 * e_a * KEV / m_a
    if hot_plasma is None:
        stop_a = P.stopping_kev_per_m(e_a, 2, SPECIES["He4"][0], material)
    else:
        tau = hot_plasma.slowing_time("He4")
        e_c = hot_plasma.critical_energy("He4")
        stop_a = (2 * e_a / tau) * (1 + (e_c / e_a) ** 1.5) / np.sqrt(v2)
    k_e2 = P.K_COULOMB * P.E_CHARGE**2
    k_alpha = 0.0
    for i in range(len(e_a) - 1):
        e_mid = 0.5 * (e_a[i] + e_a[i + 1])
        ds = (e_a[i] - e_a[i + 1]) / (0.5 * (stop_a[i] + stop_a[i + 1]))
        for s in ("D", "T"):
            m_s = SPECIES[s][0] * P.AMU
            t_max = 4 * m_a * m_s / (m_a + m_s) ** 2 * e_mid
            t = t_grid[t_grid <= t_max]
            if len(t) < 2:
                continue
            k_r = 2 * math.pi * (2 * k_e2) ** 2 / (m_s * 2 * e_mid * KEV / m_a)  # m^2 J
            boost = np.where(t >= NUCLEAR_KNOCK_ON_KEV, nuclear_boost, 1.0)
            integrand = boost * k_r / (t * KEV) ** 2 * np.interp(t, t_grid, yields[s]) * KEV
            k_alpha += density[s] * ds * float(np.trapezoid(integrand, t))

    # neutron knock-ons
    gen = np.random.default_rng(seed)
    w_d = density["D"] * N_ELASTIC_14MEV_B["D"]
    p_d = w_d / (w_d + density["T"] * N_ELASTIC_14MEV_B["T"])
    t_max_frac = {"D": 4 * 2.014 / (1 + 2.014) ** 2, "T": 4 * 3.016 / (1 + 3.016) ** 2}
    e = np.full(neutrons, E_NEUTRON_KEV)
    total = np.zeros(neutrons)
    hits = np.zeros(neutrons)
    first = None
    alive = e > 20.0
    while alive.any():
        n_alive = int(alive.sum())
        on_d = gen.uniform(size=n_alive) < p_d
        frac = np.where(on_d, t_max_frac["D"], t_max_frac["T"])
        mean_cos = np.where(on_d, N_MEAN_COSINE_14MEV["D"], N_MEAN_COSINE_14MEV["T"]) * e[alive] / E_NEUTRON_KEV
        b = np.maximum(mean_cos * (3 - mean_cos**2) / (1 - mean_cos**2), 1e-6)  # inverse Langevin function
        cos = np.log(np.exp(-b) + gen.uniform(size=n_alive) * (np.exp(b) - np.exp(-b))) / b
        t = 0.5 * (1 - cos) * frac * e[alive]
        y = np.where(on_d, np.interp(t, t_grid, y_d, left=0.0), np.interp(t, t_grid, y_t, left=0.0))
        if first is None:
            first = float(y.mean())
        total[alive] += y
        hits[alive] += 1
        e[alive] -= t
        alive = e > 20.0
    k_neutron = float(total.mean()) if first_collision_probability is None else first_collision_probability * first
    return {"k": k_alpha + k_neutron, "k_alpha": k_alpha, "k_neutron": k_neutron,
            "first_knock_on_yield": first, "collisions_per_neutron": float(hits.mean()),
            "yield_at_5MeV": {s: float(np.interp(5000.0, t_grid, yields[s])) for s in ("D", "T")}}
