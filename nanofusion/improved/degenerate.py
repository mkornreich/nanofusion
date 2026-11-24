"""Pauli blocking: fusion in fuel squeezed until its electrons are degenerate.

Every route in docs/energy.md fails at the same wall: electrons drag fast ions
to a stop long before they fuse. In fuel squeezed until its electrons are
degenerate, an ion slower than the electrons' Fermi velocity can only kick the
few electrons near the Fermi surface, because every lower state is full. The
drag per metre then stops growing with density (Fermi & Teller, Phys. Rev. 72,
399 (1947)), while fusion per metre keeps growing. Three tests:

1. ``in_flight``: a fast ion (from exploding crystals) slowing down in cold,
   squeezed fuel. In cold matter at normal density it returns at most ~0.4%
   of its energy (section 1 of docs/energy.md).
2. ``hot_ion_balance``: hot ions in cold degenerate electrons, proposed for
   proton-boron by Son & Fisch (Phys. Lett. A 329, 76 (2004)).
3. ``muon_reactivation``: a muon stuck to a fusion alpha is knocked free by
   collisions with fuel nuclei. Less electron drag means more knocks before
   the stuck ion slows down.

Physics
-------
* Electron stopping at zero temperature: Lindhard's random-phase dielectric
  function, integrated numerically over the particle-hole continuum and the
  plasmon,
      dE/dx = 4 pi Z^2 e^4 n_e L / (m_e v^2),
      L = 6/(pi chi^2) Int_0^{v/v_F} u du Int_0^inf z dz Im(-1/eps(z, u)),
  with z = k/(2 k_F), u = omega/(k v_F), chi^2 = 1/(pi k_F a_0) (Lindhard &
  Winther, Mat. Fys. Medd. Dan. Vid. Selsk. 34, no. 4 (1964), eqs. 2 and 8).
  L tends to Bethe's ln(2 m v^2 / hbar w_p) for v >> v_F, and to the
  Lindhard-Winther friction L = C1 (v/v_F)^3 for v << v_F (their eqs. 12-13;
  ~[ln(1 + 1/chi^2) - 1/(1 + chi^2)] / 2 with Thomas-Fermi screening). The
  electrons are free (pressure-ionised fuel). Bound electrons, which matter at
  a few times solid density and below, are not modelled. Linear response
  underestimates the drag where the electrons are strongly coupled: for slow
  protons in aluminium (r_s = 2.07) it gives ~60% of the measured friction.
  Squeezed fuel is far more weakly coupled (r_s ~ 0.3 at 1,000 times solid
  density), but ELECTRON_DRAG_FACTOR lets the drag be raised to test this.
* Slow ions in relativistic electrons (Fermi energy >~ 50 keV): the friction
  scales with the square of the electrons' mass at the Fermi surface, so it
  is multiplied by 1 + (p_F / m_e c)^2. This is applied to the hot-ion drain,
  the only test that reaches such densities.
* Ion stopping: binary Coulomb collisions with the cold fuel ions,
      dE/dx = sum_j 4 pi Z^2 Z_j^2 e^4 n_j L_j / (m_j v^2),
  L_j = ln(1 + (lambda / b_min)^2) / 2, with lambda the Thomas-Fermi screening
  length and the ion-sphere radius added in quadrature, and b_min the
  classical and quantum closest approaches added in quadrature.
* Fusion: the Bosch-Hale cross sections (D-T with its measured fall-off above
  550 keV) and the p-11B S-factor used in energy.py.
"""

from __future__ import annotations

import math
from functools import lru_cache

import numpy as np

from . import alternatives as ALT
from . import energy as EN
from . import physics as P

KEV = P.J_PER_KEV
HBAR = 1.054571817e-34
A0 = 5.29177210903e-11  # Bohr radius, m
E2 = P.E_CHARGE**2 / (4 * math.pi * P.EPS0)  # e^2 / (4 pi eps0), J m

# Multiplies the electrons' drag on fast ions (in_flight and the stuck muon). 1.0 = linear response (Lindhard);
# ~1.6 would match the measured friction on slow protons in aluminium. Experiments only.
ELECTRON_DRAG_FACTOR = 1.0
MU_MASS_U = 0.1134289259  # negative muon
SPECIES = dict(EN.SPECIES)
SPECIES["muHe"] = (P.MASS_U["He4"] + MU_MASS_U, 1)  # the stuck (mu-alpha)+ ion: charge +1 seen from outside

FUELS = {  # nuclei per formula unit
    "DT": {"D": 0.5, "T": 0.5},
    "boron": {"B11": 1.0},
}
ALPHA_KEV = {"DT": 3520.0, "pB": P.Q_VALUE_MEV["pB"] * 1e3 / 3}  # one alpha (D-T), or the mean of three (p-B)
CHARGED_KEV = {"DT": 3520.0, "pB": P.Q_VALUE_MEV["pB"] * 1e3}  # charged-particle energy per reaction


def fuel_ions(fuel, rho_g_cc: float) -> tuple[dict, float]:
    """Ion densities (m^-3) by species and the electron density of fully ionised fuel at rho (g/cm^3).

    ``fuel`` is a key of FUELS or a dict {species: share of nuclei}.
    """
    parts = FUELS[fuel] if isinstance(fuel, str) else fuel
    total = sum(parts.values())
    mean_mass = sum(f * SPECIES[s][0] for s, f in parts.items()) / total
    n_nuclei = rho_g_cc * 1e3 / (mean_mass * P.AMU)
    ions = {s: n_nuclei * f / total for s, f in parts.items()}
    return ions, sum(n * SPECIES[s][1] for s, n in ions.items())


def fermi_gas(n_e: float) -> dict:
    """Zero-temperature electron gas: k_F (1/m), v_F (m/s, non-relativistic), Fermi kinetic energy (keV,
    relativistic), chi^2 = 1/(pi k_F a0), x_F = p_F/(m_e c), and the plasma frequency (1/s)."""
    k_f = (3 * math.pi**2 * n_e) ** (1 / 3)
    p_f = HBAR * k_f
    x_f = p_f / (P.M_E * P.C_LIGHT)
    return {"k_F": k_f, "v_F": p_f / P.M_E, "E_F_keV": P.M_E_KEV * (math.sqrt(1 + x_f**2) - 1),
            "chi2": 1 / (math.pi * k_f * A0), "x_F": x_f,
            "omega_p": math.sqrt(n_e * P.E_CHARGE**2 / (P.EPS0 * P.M_E))}


def lindhard_winther(chi2: float) -> float:
    """Thomas-Fermi approximation to the Lindhard-Winther low-speed coefficient C1 (f1(z, 0) set to 1;
    Ferrell & Ritchie 1977). low_speed_coefficient gives the exact RPA value."""
    return 0.5 * (math.log1p(1 / chi2) - 1 / (1 + chi2))


# --- Lindhard (RPA) stopping number of a zero-temperature electron gas ------------------------

V_RATIO_MIN = 0.02  # below this, L = J (v/v_F)^3 exactly enough
V_RATIO_MAX = 300.0  # above this, L grows as Bethe's 2 ln v


def _f1(z, u):
    with np.errstate(divide="ignore", invalid="ignore"):
        a = (1 - (z - u) ** 2) * np.log(np.abs((z - u + 1) / (z - u - 1)))
        b = (1 - (z + u) ** 2) * np.log(np.abs((z + u + 1) / (z + u - 1)))
    a = np.nan_to_num(a, nan=0.0, posinf=0.0, neginf=0.0)  # 0 x log(0) at the edges of the continuum
    b = np.nan_to_num(b, nan=0.0, posinf=0.0, neginf=0.0)
    return 0.5 + (a + b) / (8 * z)


def _f2(z, u):
    inner = z + u < 1
    edge = (np.abs(z - u) < 1) & ~inner
    return np.where(inner, 0.5 * math.pi * u, np.where(edge, math.pi / (8 * z) * (1 - (z - u) ** 2), 0.0))


def low_speed_coefficient(chi2: float) -> float:
    """C1 with L = C1 (v/v_F)^3 for v -> 0: Int_0^1 z^3 dz / (z^2 + chi^2 f1(z, 0))^2 (Lindhard & Winther
    eq. 13, the exact RPA limit; their Table 1 gives 0.787 at chi^2 = 0.1)."""
    return _low_speed(float(f"{chi2:.6g}"))


@lru_cache(maxsize=512)
def _low_speed(chi2: float) -> float:
    z = np.geomspace(1e-7, 1.0, 20001)[:-1]
    f = _f1(z, np.zeros_like(z))
    return float(np.trapezoid(z**3 / (z**2 + chi2 * f) ** 2, z))


@lru_cache(maxsize=128)
def _stopping_table(chi2: float) -> tuple[np.ndarray, np.ndarray]:
    """(v/v_F grid, stopping number L) for one electron density, by numerical integration."""
    # particle-hole continuum: g(u) = Int z Im(-1/eps) dz over |z - u| < 1
    u = np.concatenate([np.linspace(0.0, 2.0, 1001)[1:], np.geomspace(2.0, V_RATIO_MAX, 801)[1:]])
    s = np.geomspace(1e-8, 1.0, 900)
    lo = np.maximum(u - 1.0, 0.0)[:, None]
    z = lo + ((u + 1.0)[:, None] - lo) * s[None, :]
    uu = u[:, None]
    e1 = 1 + chi2 / z**2 * _f1(z, uu)
    e2 = chi2 / z**2 * _f2(z, uu)
    g = np.trapezoid(z * e2 / (e1**2 + e2**2), z, axis=1)
    u0 = np.concatenate([[0.0], u])
    ug = np.concatenate([[0.0], u * g])
    continuum = np.concatenate([[0.0], np.cumsum(0.5 * (ug[1:] + ug[:-1]) * np.diff(u0))])
    # plasmon: Im(-1/eps) = pi delta(eps1) where eps1(z, u) = 0 above the continuum (u > z + 1)
    # Long-wavelength plasmons (u_p ~ chi / (sqrt(3) z)) are faster than any ion in the table; starting
    # above them also avoids f1's rounding errors at tiny z.
    zp = np.geomspace(math.sqrt(chi2) / (2 * V_RATIO_MAX), 3.0, 4000)
    lo_u = zp + 1 + 1e-9
    has = (1 + chi2 / zp**2 * _f1(zp, lo_u)) < 0
    zp, lo_u = zp[has], lo_u[has]
    hi_u = lo_u + 2 * math.sqrt(chi2) / zp + 1.0
    while True:
        neg = (1 + chi2 / zp**2 * _f1(zp, hi_u)) < 0
        if not neg.any():
            break
        hi_u = np.where(neg, 2 * hi_u, hi_u)
    for _ in range(70):
        mid = 0.5 * (lo_u + hi_u)
        neg = (1 + chi2 / zp**2 * _f1(zp, mid)) < 0
        lo_u, hi_u = np.where(neg, mid, lo_u), np.where(neg, hi_u, mid)
    u_p = 0.5 * (lo_u + hi_u)
    keep = u_p <= V_RATIO_MAX  # plasmons faster than any ion in the table never get excited
    zp, u_p = zp[keep], u_p[keep]
    h = 1e-6 * u_p
    slope = np.abs(chi2 / zp**2 * (_f1(zp, u_p + h) - _f1(zp, u_p - h)) / (2 * h))
    with np.errstate(divide="ignore", invalid="ignore"):
        weight = np.where(slope > 0, zp * u_p * math.pi / slope, 0.0)
    vr = np.geomspace(V_RATIO_MIN, V_RATIO_MAX, 600)
    plasmon = np.zeros_like(vr)
    if zp.size > 1:
        dz = np.gradient(zp)
        plasmon = ((u_p[None, :] <= vr[:, None]) * (weight * dz)[None, :]).sum(axis=1)
    total = 6 / (math.pi * chi2) * (np.interp(vr, u0, continuum) + plasmon)
    return vr, total


def stopping_number(v_ratio, chi2: float) -> np.ndarray:
    """Lindhard stopping number L at v/v_F for a gas with chi^2 = 1/(pi k_F a0) (zero temperature)."""
    key = float(f"{chi2:.4g}")
    vr, table = _stopping_table(key)
    x = np.asarray(v_ratio, dtype=float)
    mid = np.interp(np.log(np.clip(x, vr[0], vr[-1])), np.log(vr), table)
    low = table[0] * (x / vr[0]) ** 3
    high = table[-1] + 2 * np.log(np.maximum(x, vr[-1]) / vr[-1])
    return np.where(x < vr[0], low, np.where(x > vr[-1], high, mid))


def electron_stopping(e_kev, mass_u: float, z: float, n_e: float) -> np.ndarray:
    """Stopping of an ion by a zero-temperature free electron gas (Lindhard RPA), keV per metre."""
    gas = fermi_gas(n_e)
    v = np.sqrt(2 * np.asarray(e_kev, dtype=float) * KEV / (mass_u * P.AMU))
    return 4 * math.pi * z**2 * E2**2 * n_e * stopping_number(v / gas["v_F"], gas["chi2"]) / (P.M_E * v**2) / KEV


def screening_length(ions: dict, n_e: float) -> float:
    """Thomas-Fermi screening length and ion-sphere radius, added in quadrature (m)."""
    e_f = HBAR**2 * (3 * math.pi**2 * n_e) ** (2 / 3) / (2 * P.M_E)
    lam_tf = math.sqrt(2 * P.EPS0 * e_f / (3 * n_e * P.E_CHARGE**2))
    a_i = (3 / (4 * math.pi * sum(ions.values()))) ** (1 / 3)
    return math.hypot(lam_tf, a_i)


def ion_stopping(e_kev, mass_u: float, z: float, ions: dict, n_e: float) -> np.ndarray:
    """Stopping of an ion by Coulomb collisions with cold fuel ions, keV per metre."""
    e = np.asarray(e_kev, dtype=float)
    v = np.sqrt(2 * e * KEV / (mass_u * P.AMU))
    lam = screening_length(ions, n_e)
    total = np.zeros_like(e)
    for species, n in ions.items():
        m_t, z_t = SPECIES[species]
        mu = mass_u * m_t / (mass_u + m_t) * P.AMU
        b_min = np.hypot(z * z_t * E2 / (mu * v**2), HBAR / (2 * mu * v))
        log = 0.5 * np.log1p((lam / b_min) ** 2)
        total = total + 4 * math.pi * z**2 * z_t**2 * E2**2 * n * log / (m_t * P.AMU * v**2)
    return total / KEV


def _sigma_m2(projectile: str, target: str, e_cm_kev):
    if {projectile, target} == {"D", "T"}:  # with the measured fall-off above 550 keV
        return P.sigma_dt_mb(e_cm_kev, high_energy=True) * P.MB_TO_M2 * EN.DT_SPIN_FACTOR
    return EN.sigma_m2(projectile, target, e_cm_kev)


def _cumulative(y: np.ndarray, x: np.ndarray) -> np.ndarray:
    return np.concatenate([[0.0], np.cumsum(0.5 * (y[1:] + y[:-1]) * np.diff(x))])


# --- 1. a fast ion slowing down in cold squeezed fuel -----------------------------------------

def in_flight(e0_kev, projectile: str, fuel, rho_g_cc: float, points: int = 500) -> dict:
    """Fusion by one fast ion slowing down in cold fuel squeezed to rho (g/cm^3).

    Returns, for each starting energy: reactions, fusion energy / ion energy
    ("gain", full Q including neutrons, as in sections 1-2), the share of the
    ion's energy given to electrons, and its path length.
    """
    ions, n_e = fuel_ions(fuel, rho_g_cc)
    e0 = np.atleast_1d(np.asarray(e0_kev, dtype=float))
    m_p, z_p = SPECIES[projectile]
    grid = np.concatenate([[0.0], np.geomspace(0.1, max(float(e0.max()), 1.0), points)])
    e = np.maximum(grid, 0.05)
    s_e = electron_stopping(e, m_p, z_p, n_e) * ELECTRON_DRAG_FACTOR
    s_i = ion_stopping(e, m_p, z_p, ions, n_e)
    stop = s_e + s_i
    reactions = np.zeros_like(grid)
    energy = np.zeros_like(grid)
    for target, n_t in ions.items():
        q = EN.q_value(projectile, target)
        if q == 0:
            continue
        m_t = SPECIES[target][0]
        cum = _cumulative(n_t * _sigma_m2(projectile, target, grid * m_t / (m_p + m_t)) / stop, grid)
        reactions += cum
        energy += cum * q
    path = _cumulative(1 / stop, grid)
    to_electrons = _cumulative(s_e / stop, grid)
    mev = np.interp(e0, grid, energy)
    return {"rho_g_cc": rho_g_cc, "e0_keV": e0, "reactions": np.interp(e0, grid, reactions),
            "gain": mev / (e0 / 1e3), "electron_share": np.interp(e0, grid, to_electrons) / e0,
            "range_m": np.interp(e0, grid, path), "fermi": fermi_gas(n_e), "n_e": n_e}


def best_in_flight(projectile: str, fuel, rho_g_cc: float, e_lo: float = 30.0, e_hi: float = 5000.0) -> dict:
    """The starting energy that gives the largest return, and that return."""
    e0 = np.geomspace(e_lo, e_hi, 120)
    res = in_flight(e0, projectile, fuel, rho_g_cc)
    k = int(np.argmax(res["gain"]))
    return {"e0_keV": float(e0[k]), "gain": float(res["gain"][k]), "electron_share": float(res["electron_share"][k]),
            "range_m": float(res["range_m"][k]), "E_F_keV": res["fermi"]["E_F_keV"]}


def squeeze_energy_kev(n_e: float) -> float:
    """Mean kinetic energy per electron of the cold squeezed fuel (keV): the least work to squeeze it."""
    return float(ALT.degenerate_electron_energy_kev(n_e))


# --- 2. hot ions in cold degenerate electrons (Son & Fisch) ---------------------------------

def friction_coefficient(n_e: float, z: float = 1.0, t_e_kev: float = 0.0) -> float:
    """Drag per unit speed (kg/s) on an ion much slower than the electrons.

    Zero temperature: Lindhard-Winther, F = (4 Z^2 e^4 m_e^2 J / (3 pi hbar^3)) v, times
    1 + (p_F/m_e c)^2 for relativistic electrons. If t_e_kev > 0 it is combined
    harmonically with Spitzer's classical friction at that temperature, which
    takes over once the electrons are hotter than their Fermi energy. The
    harmonic sum is generous (it underestimates the drag near T_F, by up to ~2).
    """
    gas = fermi_gas(n_e)
    eta = 4 * z**2 * E2**2 * P.M_E**2 * low_speed_coefficient(gas["chi2"]) / (3 * math.pi * HBAR**3)
    eta *= 1 + gas["x_F"] ** 2
    if t_e_kev > 0:
        t = t_e_kev * KEV
        log = max(2.0, 24.0 - math.log(math.sqrt(n_e * 1e-6) / (t_e_kev * 1e3)))
        classical = 4 / 3 * math.sqrt(2 * math.pi) * z**2 * E2**2 * n_e * math.sqrt(P.M_E) * log / t**1.5
        eta = 1 / (1 / eta + 1 / classical)
    return eta


def bremsstrahlung_w_m3(n_e: float, ions: dict, t_e_kev: float, fermi_kev: float = 0.0) -> float:
    """Bremsstrahlung (W/m^3): NRL's classical formula with the relativistic and electron-electron
    corrections used in confinement.py, times the Pauli suppression for degenerate electrons.

    Suppression: only electrons within ~T of the Fermi surface can emit (the
    final state must be empty). With Kramers' spectrum this gives
    S = (pi^2.5 / 8) (T/T_F)^1.5 ~ 2.2 (T/T_F)^1.5 for T << T_F (derived here,
    not from a reference), capped at 1. fermi_kev = 0 means no suppression.
    """
    z2 = sum(n * SPECIES[s][1] ** 2 for s, n in ions.items())
    theta = t_e_kev / P.M_E_KEV
    relativistic = 1 + 0.7936 * theta + 1.874 * theta**2 + 3 / math.sqrt(2) * theta * n_e / z2
    classical = 5.35e-37 * n_e * z2 * math.sqrt(t_e_kev) * relativistic
    if fermi_kev <= 0:
        return classical
    return classical * min(1.0, math.pi**2.5 / 8 * (t_e_kev / fermi_kev) ** 1.5)


def pb_fuel(boron_fraction: float) -> dict:
    return {"p": 1 - boron_fraction, "B11": boron_fraction}


def alpha_share_to_ions(fuel: str, rho_g_cc: float, boron_fraction: float = 0.15) -> float:
    """Share of a fusion alpha's energy given to ions while it slows down (zero-temperature electrons,
    target ions at rest; generous, since hot target ions take less)."""
    parts = FUELS["DT"] if fuel == "DT" else pb_fuel(boron_fraction)
    res = in_flight(ALPHA_KEV[fuel], "He4", parts, rho_g_cc, points=300)
    return 1.0 - float(res["electron_share"][0])


def _drain(ions: dict, n_e: float, t_i_kev: float, t_e_kev: float) -> float:
    """Power (W/m^3) flowing from ions at t_i to electrons at t_e: each ion loses (2 eta / M)(3/2)(T_i - T_e)."""
    return sum(n * 2 * friction_coefficient(n_e, SPECIES[s][1], t_e_kev) / (SPECIES[s][0] * P.AMU)
               * 1.5 * (t_i_kev - t_e_kev) * KEV for s, n in ions.items())


def hot_ion_balance(fuel: str, rho_g_cc: float, t_i_kev: float, boron_fraction: float = 0.15,
                    alpha_channelling: float = 0.0) -> dict:
    """Steady power balance of hot ions in degenerate electrons (W/m^3).

    * fusion: Maxwellian reactivity at t_i_kev. Its charged products give a
      share f_i of their energy to the ions, the rest to the electrons;
    * drain: every ion is slower than the electrons, so it loses energy at
      dE/dt = -(2 eta / M)(E - 3/2 T_e). While the electrons are degenerate
      eta does not depend on density;
    * the electrons must radiate what they receive (drain + (1 - f_i) of the
      charged power) as bremsstrahlung, with the Pauli suppression. That fixes
      their temperature T_e.
    "self_heating" = f_i P_charged / drain at that T_e (equivalently, charged
    fusion power / radiation). Above 1 the ions keep themselves hot.
    alpha_channelling (untested theory, Fisch & Rax 1992): waves pass that share of what the alphas would give
    the electrons to the ions instead.
    "self_heating_cold" = f_i P_charged / drain with the electrons held at zero
    temperature instead: they drain the ions at the full degenerate rate, and
    whatever they receive is taken to vanish. "time_to_heat_s": how long
    zero-temperature electrons would take to hold their Fermi energy as heat if
    they could not radiate at all.
    """
    parts = FUELS["DT"] if fuel == "DT" else pb_fuel(boron_fraction)
    ions, n_e = fuel_ions(parts, rho_g_cc)
    gas = fermi_gas(n_e)
    if fuel == "DT":
        rate = ions["D"] * ions["T"] * EN.reactivity("D", "T", t_i_kev)
        per_fuel_ion = rate / ions["D"]
    else:
        rate = ions["p"] * ions["B11"] * EN.reactivity("p", "B11", t_i_kev)
        per_fuel_ion = rate / ions["B11"]
    p_fusion = rate * EN.q_value(*(("D", "T") if fuel == "DT" else ("p", "B11"))) * 1e3 * KEV
    p_charged = rate * CHARGED_KEV[fuel] * KEV
    f_i = alpha_share_to_ions(fuel, rho_g_cc, boron_fraction)
    f_i += alpha_channelling * (1 - f_i)  # untested: waves hand this share of the electrons' part to the ions
    drain_cold = _drain(ions, n_e, t_i_kev, 0.0)

    def excess(t_e: float) -> float:  # heat into the electrons minus what they radiate; falls as t_e rises
        return (_drain(ions, n_e, t_i_kev, t_e) + (1 - f_i) * p_charged
                - bremsstrahlung_w_m3(n_e, ions, t_e, gas["E_F_keV"]))

    lo, hi = 1e-3, 10 * t_i_kev
    if excess(hi) > 0:
        lo = hi
    for _ in range(60):
        mid = math.sqrt(lo * hi)
        if hi / lo < 1.0005:
            break
        lo, hi = (mid, hi) if excess(mid) > 0 else (lo, mid)
    t_e = math.sqrt(lo * hi)
    radiated = bremsstrahlung_w_m3(n_e, ions, t_e, gas["E_F_keV"])
    eps = squeeze_energy_kev(n_e)
    heat_to_fermi = n_e * (math.hypot(eps, 1.5 * gas["E_F_keV"]) - eps) * KEV  # J/m^3, as in alternatives.ignition_ball
    return {"rho_g_cc": rho_g_cc, "T_i_keV": t_i_kev, "E_F_keV": gas["E_F_keV"], "fusion_W_m3": p_fusion,
            "charged_W_m3": p_charged, "alpha_to_ions": f_i, "drain_cold_W_m3": drain_cold,
            "self_heating_cold": f_i * p_charged / drain_cold, "T_e_keV": t_e, "T_e_over_T_F": t_e / gas["E_F_keV"],
            "radiated_W_m3": radiated, "self_heating": p_charged / radiated,
            "time_to_heat_s": heat_to_fermi / (drain_cold + (1 - f_i) * p_charged),
            "burn_rate_per_s": per_fuel_ion, "squeeze_J_m3": n_e * eps * KEV,
            "ion_heat_J_m3": 1.5 * sum(ions.values()) * t_i_kev * KEV}


def self_heating_density(fuel: str, t_i_kev: float, boron_fraction: float = 0.15, cold: bool = False,
                         rho_lo: float = 1.0, rho_hi: float = 1e9) -> float:
    """Density (g/cm^3) above which hot ions keep themselves hot (inf if none up to rho_hi).

    cold=True holds the electrons at zero temperature (see hot_ion_balance)."""
    key = "self_heating_cold" if cold else "self_heating"

    def ok(rho: float) -> bool:
        return hot_ion_balance(fuel, rho, t_i_kev, boron_fraction)[key] >= 1

    if not ok(rho_hi):
        return math.inf
    if ok(rho_lo):
        return rho_lo
    lo, hi = rho_lo, rho_hi
    while hi / lo > 1.01:
        mid = math.sqrt(lo * hi)
        lo, hi = (lo, mid) if ok(mid) else (mid, hi)
    return math.sqrt(lo * hi)


# --- 3. stripping a stuck muon (muon-catalysed fusion) -------------------------------------

MU_HE_REDUCED_MASS = (MU_MASS_U * P.MASS_U["He4"] / (MU_MASS_U + P.MASS_U["He4"])) / (P.M_E / P.AMU)  # in m_e
MU_HE_BINDING_KEV = 4 * 13.605693e-3 * MU_HE_REDUCED_MASS  # 1s binding, about 10.9 keV
MU_HE_KEV = 3460.0  # starting kinetic energy of the stuck (mu-alpha)+ ion (17.59 MeV shared with the neutron)
# Rudd et al., Rev. Mod. Phys. 57, 965 (1985): ionisation of atomic hydrogen by protons,
# sigma = 1/(1/s_low + 1/s_high), s_low = 4 pi a0^2 C x^D, s_high = 4 pi a0^2 (A ln(1+x) + B)/x,
# x = (m_e/m_p) E_p / I, I = 13.6 eV.
RUDD_H = (0.28, 1.15, 0.44, 0.907)
# Stripping fraction at liquid density (phi ~ 1.2-1.5): theory, e.g. Stodden et al., PRA 41, 1281 (1990), as used
# by Kamimura, Kino & Yamashita, PRC 107, 034607 (2023). It rises only weakly with density in that range.
REACTIVATION_LIQUID = 0.35
LIQUID_NUCLEI_M3 = 4.25e28  # liquid hydrogen density
STICKING_INITIAL = 0.009  # initial sticking, before stripping (theory 0.86-0.93%)
DTMU_FUSION_RATE = 1.0e12  # 1/s: fusion inside the d-t-mu molecule (1.0-1.15e12), the fastest any cycle can go


def hydrogen_ionisation_m2(e_p_kev):
    """Ionisation of H(1s) by protons of energy e_p_kev (Rudd's fit), m^2."""
    a, b, c, d = RUDD_H
    x = np.asarray(e_p_kev, dtype=float) * 1e3 * (P.M_E / (P.MASS_U["p"] * P.AMU)) / 13.6057
    unit = 4 * math.pi * A0**2
    low = unit * c * x**d
    high = unit * (a * np.log1p(x) + b) / x
    return 1 / (1 / low + 1 / high)


def stripping_m2(e_kev):
    """Ionisation of the (mu-alpha)+ 1s state by a singly charged nucleus, m^2 (lab energy of the ion, keV).

    Born scaling of hydrogen: a one-electron (here one-muon) ion of nuclear
    charge Z and reduced mass m_r has sigma(v) = sigma_H(v / Z) (m_e/m_r)^2 / Z^4,
    with the same projectile.
    """
    m = SPECIES["muHe"][0]
    v = np.sqrt(2 * np.asarray(e_kev, dtype=float) * KEV / (m * P.AMU))
    e_h = 0.5 * P.MASS_U["p"] * P.AMU * (v / 2) ** 2 / KEV
    return hydrogen_ionisation_m2(e_h) / (MU_HE_REDUCED_MASS**2 * 16)


def orbit_matching_kev() -> float:
    """Energy (keV) at which the stuck ion moves as fast as the muon in its orbit (2 alpha c), ~410 keV.
    Below it, stripping falls off steeply and the Born-scaled cross section is least reliable."""
    return 0.5 * SPECIES["muHe"][0] * P.AMU * (2 * P.ALPHA_FS * P.C_LIGHT) ** 2 / KEV


def _strip_integral(rho_g_cc: float, cold: bool, min_kev: float = 1.0) -> float:
    """Int n sigma_strip dE / S over the stuck ion's slowing down (expected number of strippings, uncalibrated)."""
    m = SPECIES["muHe"][0]
    e = np.geomspace(min_kev, MU_HE_KEV, 600)
    if cold:  # liquid or solid D-T: bound electrons (the project's cold-matter stopping), scaled to density
        n = LIQUID_NUCLEI_M3 * rho_g_cc / liquid_density_g_cc()
        s_e = P.stopping_kev_per_m(e, 1, m, "DT") * n / P.MATERIALS["DT"][0][1]
        ions = {"D": n / 2, "T": n / 2}
        n_e = n
    else:
        ions, n_e = fuel_ions("DT", rho_g_cc)
        n = sum(ions.values())
        s_e = electron_stopping(e, m, 1, n_e) * ELECTRON_DRAG_FACTOR
    s_i = ion_stopping(e, m, 1, ions, n_e)
    y = n * stripping_m2(e) / (s_e + s_i)
    return float(np.trapezoid(y, e))


def liquid_density_g_cc() -> float:
    return LIQUID_NUCLEI_M3 * 0.5 * (P.MASS_U["D"] + P.MASS_U["T"]) * P.AMU / 1e3


@lru_cache(maxsize=1)
def stripping_calibration() -> float:
    """Factor on the stripping integral that reproduces the reactivation measured at liquid density.

    The single-state model leaves out excitation and de-excitation of the
    stuck muon, so it is tuned to the measured value once and then used at
    every density.
    """
    return -math.log(1 - REACTIVATION_LIQUID) / _strip_integral(liquid_density_g_cc(), cold=True)


def muon_reactivation(rho_g_cc: float, fast_only: bool = False) -> dict:
    """Stripping fraction R, effective sticking, and the fusions per muon sticking alone would allow,
    in squeezed (pressure-ionised) D-T. fast_only counts only strippings above orbit_matching_kev()."""
    k = stripping_calibration()
    lowest = orbit_matching_kev() if fast_only else 1.0
    survive = math.exp(-k * _strip_integral(rho_g_cc, cold=False, min_kev=lowest))  # 1 - R, without rounding
    sticking = STICKING_INITIAL * survive
    return {"rho_g_cc": rho_g_cc, "reactivation": 1 - survive, "sticking": sticking,
            "fusions_limit": 1 / sticking if sticking > 0 else math.inf}


def squeeze_lifetime(rho_g_cc: float, radius_m: float) -> float:
    """How long cold fuel squeezed to rho stays squeezed: R / (4 c_s), with the degenerate sound speed."""
    return radius_m / (4 * ALT.cold_compression(rho_g_cc)["sound_speed_m_s"])
