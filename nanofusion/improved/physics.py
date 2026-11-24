"""Physics inputs for the improved model: constants, fusion cross sections,
stopping power, and the laser pulse.

Sources
-------
* D+D cross sections: Bosch & Hale, Nucl. Fusion 32, 611 (1992), with S-factor
  coefficients for D(d,p)T and D(d,n)3He.
* D+6Li -> 2 alpha: bare S(0) = 21.7 MeV b (Grineviciute et al., Phys. Rev. C
  91, 014601 (2015), arXiv:1406.5735). This is the (d,alpha) channel only.
* D+T: Bosch & Hale (1992) T(d,n)4He fit (valid 0.5-550 keV; held constant above).
* p+11B -> 3 alpha: piecewise S-factor of Nevins & Swain (2000) as refitted in
  arXiv:2601.00241 (Table 1), E_G = 22.589 MeV, valid to 10 MeV.
* D+7Li -> n + 8Be (n0+n1): S(0) = 5.4 MeV b (Sabourov et al., Phys. Rev. C
  73, 015801 (2006)). This is a lower bound for all 7Li+d channels; a 2023
  laser-plasma measurement gave 24 +- 8 MeV b near 173 keV.
* Stopping of D in LiD: Lindhard-Scharff electronic stopping at low energy,
  Bethe at high energy, combined as 1/S = 1/S_low + 1/S_high. Mean
  excitation energy I = 36.5 eV (LiH, PDG). Expect factor-of-2 accuracy.
* LiD crystal: rock salt, Li-D distance 2.04 A (lattice constant about 4.07 A).
"""

from __future__ import annotations

import math
from dataclasses import dataclass

import numpy as np

# --- constants (SI) ----------------------------------------------------------
E_CHARGE = 1.602176634e-19
K_COULOMB = 8.9875517873681764e9
EPS0 = 8.8541878128e-12
C_LIGHT = 2.99792458e8
AMU = 1.66053906660e-27
AMU_KEV = 931494.10242  # u c^2 in keV
M_E = 9.1093837015e-31
M_E_KEV = 510.99895
J_PER_KEV = 1.602176634e-16
MB_TO_M2 = 1e-31
ALPHA_FS = 7.2973525693e-3
BOHR_RADIUS_NM = 0.0529177
BOHR_VELOCITY = 2.18769126e6
E2_EV_NM = 1.439964547  # e^2 / (4 pi eps0) in eV nm

MASS_U = {"p": 1.00727646688, "D": 2.01410177812, "T": 3.01604928, "He4": 4.00260325,
          "Li6": 6.0151228874, "Li7": 7.016003434, "B11": 11.0093054}
CHARGE = {"D": 1, "Li6": 3, "Li7": 3}
SPACING = 2.04e-10  # Li-D nearest-neighbour distance (m)

S_DLI_KEV_MB = {"Li6": 21.7e3 * 1e3, "Li7": 5.4e3 * 1e3}  # MeV b -> keV mb
E_CAP_KEV = 4900.0  # Bosch-Hale validity limit; cross sections held constant above it
Q_VALUE_MEV = {"DD": 3.65, "Li6": 22.37, "Li7": 15.03, "DT": 17.59, "pB": 8.68}  # D+D: average of its two branches
E_CAP_DT_KEV = 550.0  # upper limit of the Bosch-Hale D-T fit used here
E_CAP_DT_HIGH_KEV = 4700.0  # upper limit of Bosch & Hale's high-energy D-T fit (optional, see sigma_dt_mb)
I_LID_EV = 36.5


def reduced_mass_u(a: float, b: float) -> float:
    return a * b / (a + b)


def gamow_constant(z1: int, z2: int, mu_u: float) -> float:
    """B_G in keV^0.5 for sigma = S/E exp(-B_G/sqrt(E)), E = centre-of-mass energy in keV."""
    return math.pi * ALPHA_FS * z1 * z2 * math.sqrt(2.0 * mu_u * AMU_KEV)


def sigma_dd_mb(e_kev):
    """Total D+D cross section (both branches), millibarn; E = CM energy in keV."""
    e = np.clip(np.asarray(e_kev, dtype=float), 1e-3, E_CAP_KEV)
    with np.errstate(over="ignore"):
        return _sigma_dd(e, e_kev)


def _sigma_dd(e, e_kev):
    s_p = 5.5576e4 + e * (2.1054e2 + e * (-3.2638e-2 + e * (1.4987e-6 + e * 1.8181e-10)))
    s_n = 5.3701e4 + e * (3.3027e2 + e * (-1.2706e-1 + e * (2.9327e-5 + e * -2.5151e-9)))
    out = (s_p + s_n) / (e * np.exp(31.3970 / np.sqrt(e)))
    return np.where(np.asarray(e_kev) < 1e-3, 0.0, out)


def sigma_dt_mb(e_kev, high_energy: bool = False):
    """D+T -> n + alpha (Bosch-Hale), millibarn; E = CM energy in keV.

    By default the 0.5-550 keV fit is used and held constant above 550 keV,
    which overstates the cross section there. ``high_energy=True`` switches
    to Bosch & Hale's second fit from 550 to 4700 keV (it joins the first
    within 1%), where the cross section falls from 0.26 b to below 0.1 b, and
    holds it constant above 4700 keV.
    """
    raw = np.asarray(e_kev, dtype=float)
    e = np.clip(raw, 1e-3, E_CAP_DT_KEV)
    num = 6.927e4 + e * (7.454e8 + e * (2.050e6 + e * (5.2000e4 + e * 0.0)))
    den = 1.0 + e * (6.380e1 + e * (-9.950e-1 + e * (6.981e-5 + e * 1.728e-4)))
    with np.errstate(over="ignore"):
        out = (num / den) / (e * np.exp(34.3827 / np.sqrt(e)))
    if high_energy:
        eh = np.clip(raw, E_CAP_DT_KEV, E_CAP_DT_HIGH_KEV)
        s_high = -1.4714e6 / (1.0 + eh * (-8.4127e-3 + eh * (4.7983e-6 + eh * (-1.0748e-9 + eh * 8.5184e-14))))
        out = np.where(raw > E_CAP_DT_KEV, s_high / (eh * np.exp(34.3827 / np.sqrt(eh))), out)
    return np.where(raw < 1e-3, 0.0, out)


_PB_A = (2.0235e6, 4.0102e6, 1.3220e6, 4.9451e6, 4.3430e5)
_PB_E = (622.2, 1388.4, 2492.4, 3528.6, 4703.6)  # keV
_PB_W = (99.6, 449.9, 238.6, 398.5, 152.5)  # keV


def s_factor_pb11_mev_b(e_kev):
    """p-11B S-factor (MeV b), CM energy in keV (arXiv:2601.00241, Table 1)."""
    e = np.asarray(e_kev, dtype=float)
    s1 = 197 + 0.240 * e + 2.31e-4 * e**2 + 1.82e4 / ((e - 148.0) ** 2 + 2.35**2)
    x = (e - 400.0) / 100.0
    s2 = 330.2 + 102.436 * x - 58.481 * x**2 + 0.0933 * x**5
    s3 = 0.209689 + sum(a / ((e - c) ** 2 + w**2) for a, c, w in zip(_PB_A, _PB_E, _PB_W))
    return np.where(e <= 400, s1, np.where(e <= 700, s2, s3))


def sigma_pb11_mb(e_kev):
    """p + 11B -> 3 alpha cross section, millibarn; held constant above 10 MeV."""
    e = np.clip(np.asarray(e_kev, dtype=float), 1e-3, 1e4)
    with np.errstate(over="ignore"):
        out = s_factor_pb11_mev_b(e) * 1e6 / e * np.exp(-np.sqrt(22589.0 / e))  # MeV b -> keV mb
    return np.where(np.asarray(e_kev) < 1e-3, 0.0, out)


def sigma_dli_mb(e_kev, s_kev_mb: float, bg: float):
    e = np.clip(np.asarray(e_kev, dtype=float), 1e-3, E_CAP_KEV)
    with np.errstate(over="ignore"):
        out = s_kev_mb / (e * np.exp(bg / np.sqrt(e)))
    return np.where(np.asarray(e_kev) < 1e-3, 0.0, out)


# --- stopping power of deuterons in LiD ----------------------------------------------

def lid_molecule_density() -> float:
    """LiD molecules per m^3 (4 per cubic cell of side 2 x spacing)."""
    return 4.0 / (2 * SPACING) ** 3


def stopping_in_lid_kev_per_m(e_lab_kev, z1: int = 1, mass_u: float = MASS_U["D"]):
    """Electronic stopping of an ion (charge z1, mass in u) in cold LiD, keV per metre (+-~50%).

    Lindhard-Scharff at low speed, Bethe at high speed, combined as
    1/S = 1/S_low + 1/S_high. Heavier ions (z1 > 1) carry a speed-dependent
    effective charge z_eff = z1 (1 - exp(-0.92 v / (v0 z1^(2/3)))).
    """
    e = np.maximum(np.asarray(e_lab_kev, dtype=float), 1e-6)
    n_mol_nm3 = lid_molecule_density() * 1e-27
    v = np.sqrt(2 * e * J_PER_KEV / (mass_u * AMU))

    def lss(z2):  # Lindhard-Scharff electronic stopping cross section per target atom, eV nm^2 per (v/v0)
        z = (z1 ** (2 / 3) + z2 ** (2 / 3)) ** 1.5
        return 8 * math.pi * E2_EV_NM * BOHR_RADIUS_NM * z1 ** (7 / 6) * z2 / z

    s_low = (lss(3) + lss(1)) * (v / BOHR_VELOCITY) * n_mol_nm3  # eV/nm
    z_eff = z1 if z1 == 1 else z1 * (1 - np.exp(-0.92 * v / (BOHR_VELOCITY * z1 ** (2 / 3))))
    me_v2_ev = M_E_KEV * 1e3 * (v / C_LIGHT) ** 2
    n_e = 4 * n_mol_nm3
    # ln(1 + x) instead of ln(x): equal at high energy, and stays positive at low
    # energy, where the Lindhard term then dominates.
    log = np.log1p(2 * me_v2_ev / I_LID_EV)
    s_high = 4 * math.pi * E2_EV_NM**2 * z_eff**2 * n_e / me_v2_ev * log  # eV/nm
    s = 1.0 / (1.0 / s_low + 1.0 / np.maximum(s_high, 1e-300))
    return s * 1e-3 / 1e-9  # eV/nm -> keV/m


MATERIALS = {
    # name: ({Z: atoms per m^3}, mean excitation energy I in eV)
    "LiD": ({3: None, 1: None}, I_LID_EV),  # filled in lazily from the lattice (one Li and one D per molecule)
    "decaborane": ({5: 4.56e28, 1: 6.38e28}, 0.0),  # B10H14 at 0.94 g/cm^3; I from Bragg's rule
    "boron": ({5: 1.30e29}, 76.0),  # 2.34 g/cm^3
    "D2": ({1: 6.0e28}, 0.0),  # solid deuterium, 0.2 g/cm^3
    "DT": ({1: 4.5e28}, 0.0),  # solid deuterium-tritium, 0.22 g/cm^3
    "ErD3": ({68: 2.67e28, 1: 8.0e28}, 0.0),  # erbium deuteride, 8e22 D/cm^3 (NASA lattice-confinement tests)
    "TiD2": ({22: 5.0e28, 1: 1.0e29}, 0.0),  # titanium deuteride, 1e23 D/cm^3
    "CD2": ({6: 4.0e28, 1: 8.0e28}, 0.0),  # deuterated polyethylene, about 1.06 g/cm^3
}
_I_ATOM = {1: 19.2, 3: 40.0, 5: 76.0, 6: 78.0, 22: 233.0, 68: 720.0}  # eV (ICRU/PDG elemental values, approximate)


def _material(name: str) -> tuple[dict, float]:
    atoms, i_ev = MATERIALS[name]
    if name == "LiD":
        n = lid_molecule_density()
        atoms = {3: n, 1: n}
    if not i_ev:
        n_e = sum(z * n for z, n in atoms.items())
        i_ev = math.exp(sum(z * n * math.log(_I_ATOM[z]) for z, n in atoms.items()) / n_e)
    return atoms, i_ev


def stopping_kev_per_m(e_lab_kev, z1: int, mass_u: float, material: str):
    """Electronic stopping of an ion (charge z1, mass in u) in a cold material, keV/m (+-~50%).

    Same model as stopping_in_lid_kev_per_m, for any material in MATERIALS.
    """
    atoms, i_ev = _material(material)
    e = np.maximum(np.asarray(e_lab_kev, dtype=float), 1e-6)
    v = np.sqrt(2 * e * J_PER_KEV / (mass_u * AMU))
    s_low = 0.0
    for z2, n in atoms.items():
        z = (z1 ** (2 / 3) + z2 ** (2 / 3)) ** 1.5
        s_low = s_low + 8 * math.pi * E2_EV_NM * BOHR_RADIUS_NM * z1 ** (7 / 6) * z2 / z * (v / BOHR_VELOCITY) * n * 1e-27
    z_eff = z1 if z1 == 1 else z1 * (1 - np.exp(-0.92 * v / (BOHR_VELOCITY * z1 ** (2 / 3))))
    me_v2_ev = M_E_KEV * 1e3 * (v / C_LIGHT) ** 2
    n_e = sum(z * n for z, n in atoms.items()) * 1e-27
    s_high = 4 * math.pi * E2_EV_NM**2 * z_eff**2 * n_e / me_v2_ev * np.log1p(2 * me_v2_ev / i_ev)
    s = 1.0 / (1.0 / s_low + 1.0 / np.maximum(s_high, 1e-300))
    return s * 1e6  # eV/nm -> keV/m


def stopping_d_in_lid_kev_per_m(e_lab_kev):
    """Electronic stopping of deuterons in cold LiD, keV per metre (+-~50%)."""
    return stopping_in_lid_kev_per_m(e_lab_kev, 1, MASS_U["D"])


def beam_target(e_lab_kev, projectile: str, catcher: str, lithium: str, points: int = 2000) -> dict:
    """Reactions per ion stopping in a cold catcher, for each initial lab energy.

    projectile: "D" or "Li" (the crystal's lithium isotope). catcher: "LiD"
    (lithium deuteride with the same isotope) or "LiT" (lithium tritide).
    Returns {channel: array}. Each ion reacts at most once.
    """
    e0 = np.atleast_1d(np.asarray(e_lab_kev, dtype=float))
    n = lid_molecule_density()
    m_d, m_t, m_li = MASS_U["D"], MASS_U["T"], MASS_U[lithium]
    bg_dli = gamow_constant(1, 3, reduced_mass_u(m_d, m_li))
    s_dli = S_DLI_KEV_MB[lithium]
    top = max(float(e0.max()) if e0.size else 1.0, 1.0)
    grid = np.concatenate([[0.0], np.geomspace(0.05, top, points)])
    if projectile == "D":
        stop = stopping_in_lid_kev_per_m(np.maximum(grid, 1e-3), 1, m_d)
        channels = {"D+Li": sigma_dli_mb(grid * m_li / (m_d + m_li), s_dli, bg_dli)}
        if catcher == "LiD":
            channels["D+D"] = sigma_dd_mb(grid * 0.5)
        else:
            channels["D+T"] = sigma_dt_mb(grid * m_t / (m_d + m_t))
    else:  # lithium projectile: only the deuterium in an LiD catcher (no Li+T data used)
        stop = stopping_in_lid_kev_per_m(np.maximum(grid, 1e-3), 3, m_li)
        channels = {"Li+D": sigma_dli_mb(grid * m_d / (m_d + m_li), s_dli, bg_dli)} if catcher == "LiD" else {}
    out = {}
    de = np.diff(grid)
    for name, sigma in channels.items():
        dydE = n * sigma * MB_TO_M2 / stop
        cum = np.concatenate([[0.0], np.cumsum(0.5 * (dydE[1:] + dydE[:-1]) * de)])
        out[name] = np.interp(e0, grid, cum)
    total = sum(out.values()) if out else np.zeros_like(e0)
    p_any = -np.expm1(-total)
    with np.errstate(invalid="ignore", divide="ignore"):
        share = np.where(total > 0, p_any / np.where(total > 0, total, 1), 1.0)
    return {k: v * share for k, v in out.items()}


def thick_target_yield(e_lab_kev, lithium: str, points: int = 2000) -> tuple[np.ndarray, np.ndarray]:
    """Reactions per deuteron stopping in cold LiD: (D+Li, D+D), for each initial lab energy."""
    y = beam_target(e_lab_kev, "D", "LiD", lithium, points)
    return y["D+Li"], y["D+D"]


Q_BY_CHANNEL = {"D+D": Q_VALUE_MEV["DD"], "D+T": Q_VALUE_MEV["DT"]}


def channel_q(channel: str, lithium: str) -> float:
    return Q_BY_CHANNEL.get(channel, Q_VALUE_MEV[lithium])


# --- laser pulse -------------------------------------------------------------------

@dataclass(frozen=True)
class Pulse:
    """Linearly polarised (x) Gaussian pulse travelling along +z, from a vector potential.

    E(t) = -dA/dt with A(t) = -(E0/omega) f(t) sin(omega (t - t_c) + phase),
    f = exp(-(t - t_c)^2 / tau_g^2). The field has zero net area, as a real
    propagating pulse must. ``phase`` = 0 puts a field crest at the pulse
    centre; 90 degrees a zero crossing.
    """

    e0: float  # peak field, V/m
    wavelength: float  # m
    fwhm: float  # intensity FWHM, s
    centre: float = 0.0  # s
    phase: float = 0.0  # rad

    @property
    def omega(self) -> float:
        return 2 * math.pi * C_LIGHT / self.wavelength

    @property
    def tau_g(self) -> float:
        return self.fwhm / math.sqrt(2 * math.log(2))

    def field(self, t: float) -> float:
        s = t - self.centre
        f = math.exp(-(s / self.tau_g) ** 2)
        theta = self.omega * s + self.phase
        return self.e0 * f * (math.cos(theta) - 2 * s / (self.omega * self.tau_g**2) * math.sin(theta))

    def fluence(self) -> float:
        """Energy per area, J/m^2 (cycle-averaged intensity integrated over time)."""
        return 0.5 * C_LIGHT * EPS0 * self.e0**2 * self.tau_g * math.sqrt(math.pi / 2)


def intensity_w_cm2(e0: float) -> float:
    return 0.5 * C_LIGHT * EPS0 * e0**2 * 1e-4


def field_from_intensity(i_w_cm2: float) -> float:
    return math.sqrt(2 * i_w_cm2 * 1e4 / (C_LIGHT * EPS0))


def electron_a0(e0: float, wavelength: float) -> float:
    """Normalised vector potential; above ~1 electrons are relativistic."""
    omega = 2 * math.pi * C_LIGHT / wavelength
    return E_CHARGE * e0 / (M_E * omega * C_LIGHT)


def critical_density(wavelength: float) -> float:
    omega = 2 * math.pi * C_LIGHT / wavelength
    return EPS0 * M_E * omega**2 / E_CHARGE**2


def stripped_charge_density() -> float:
    """Charge density of fully ionised LiD, C/m^3 (4 Li3+ and 4 D+ per cubic cell)."""
    return 16 * E_CHARGE / (2 * SPACING) ** 3


def coulomb_explosion_max_kev(radius: float, z: int) -> float:
    """Final energy of an ion of charge z at the surface of a uniformly charged sphere."""
    return z * E_CHARGE * stripped_charge_density() * radius**2 / (3 * EPS0) / J_PER_KEV


def stripping_field(radius: float) -> float:
    """Field at the surface of a fully ionised LiD sphere; the laser must exceed it to pull electrons out."""
    return stripped_charge_density() * radius / (3 * EPS0)
