"""Three more walls, three new attacks.

docs/energy.md keeps meeting the same few walls. Each test here attacks one
that no earlier section touched:

1. X-rays. Proton-boron fails because its bremsstrahlung carries off more
   than fusion makes. ``opaque_ball``: make the fuel thick enough to trap its
   own X-rays.
2. Holding the fuel. Section 12 used a tamper and a field. ``fibre_pinch``:
   run megaamperes through a fibre of solid fuel, so that the current's own
   field squeezes and holds it (a dense Z-pinch), with no laser at all.
3. Electron drag. ``channelled_return``: steer the ions along the open
   channels of a tritide crystal, where there are fewer electrons and the
   ions' flux peaks between the rows of atoms.

Physics
-------
* X-ray trapping, frequency by frequency: free-free absorption
  alpha_nu = 3.7e8 T^-1/2 n_e sum(n_i Z^2) nu^-3 (1 - exp(-h nu / kT)) cm^-1 (cgs,
  Gaunt factor 1; its Rosseland mean is Kramers' 2.06e-25 T^-3.5 n_e sum(n_i Z^2)),
  with Thomson scattering folded in as an effective optical depth: the
  lengthened path tau_a (1 + tau_s) while that is small, the thermalisation
  depth sqrt(tau_a (tau_a + tau_s)) beyond. Bremsstrahlung is emitted with a flat
  spectrum ~exp(-h nu / kT) and escapes a uniform sphere with Osterbrock's
  probability P(tau) = 3/(4 tau) [1 - 1/(2 tau^2) + (1/tau + 1/(2 tau^2)) exp(-2 tau)],
  which tends to blackbody emission from the surface when tau >> 1. The
  trapped frequencies fill with Planck radiation, which the ball must pay
  for, and so do thermal electron-positron pairs (relativistic equilibrium
  density) once photons near m_e c^2 are trapped.
* Hot ions and electrons: the ion-to-electron drain and the Pauli-suppressed
  bremsstrahlung of degenerate.py (the drain merges into Spitzer's for hot
  electrons); the alpha share to ions from Stix's critical energy.
* Z-pinch: Bennett's equilibrium mu0 I^2 / 8 pi = (N_i + N_e) k T; the pinch
  lives k Alfven times, a / v_A, with v_A = I sqrt(mu0 / (4 pi M)) for a mass M
  per metre (independent of the radius); fusion alphas are held by the
  pinch's own field once a B(a) > 0.27 T m, i.e. above ~1.35 MA.
* Channelling: Lindhard's critical angle psi_1 = sqrt(2 Z1 Z2 e^2 / (E d)).
"""

from __future__ import annotations

import math

import numpy as np

from . import degenerate as D
from . import energy as EN
from . import physics as P

KEV = P.J_PER_KEV
K_B = 1.380649e-23
A_RAD = 7.565723e-16  # radiation constant, J m^-3 K^-4
SIGMA_T = 6.6524587e-29  # Thomson cross section, m^2
MU0 = 4e-7 * math.pi
ALPHA_LARMOR_TM = 0.27  # 3.5 MeV alpha: Larmor radius x field, T m (as in confinement.py)


# --- 1. trapping the X-rays: an opaque proton-boron ball -------------------------------------

def _bessel_k2(z: float) -> float:
    """Modified Bessel function K_2(z) = Int_0^inf exp(-z cosh t) cosh(2t) dt, for z >~ 0.1."""
    t = np.linspace(0.0, math.acosh(1 + 60.0 / z), 4000)
    return float(np.trapezoid(np.exp(-z * (np.cosh(t) - 1)) * np.cosh(2 * t), t)) * math.exp(-z)


def pair_density(n_e0: float, t_kev: float, photon_share: float = 1.0) -> float:
    """Thermal positrons per m^3 in equilibrium with radiation at t_kev, for n_e0 electrons from the fuel.

    photon_share: the photons near m_e c^2 as a share of their full Planck number. Pairs are made and
    destroyed in pairs of photons, so n+ n- scales as its square.

    n+ n- = n_eq^2, with the relativistic (Maxwell-Juttner) density of one species at zero chemical potential,
    n_eq = (m c / hbar)^3 theta K_2(1/theta) / pi^2, theta = kT / m c^2. For kT << m c^2 this is the familiar
    2 (m k T / 2 pi hbar^2)^1.5 exp(-m c^2 / kT).
    """
    theta = t_kev / P.M_E_KEV
    n_eq = photon_share * (P.M_E * P.C_LIGHT / D.HBAR) ** 3 * theta * _bessel_k2(1 / theta) / math.pi**2
    return 0.5 * (math.sqrt(n_e0**2 + 4 * n_eq**2) - n_e0)


def free_free_opacity(ions: dict, n_e: float, t_kev: float) -> float:
    """Rosseland-mean free-free absorption coefficient (1/m), Gaunt factor 1 (Kramers' law)."""
    t = t_kev * KEV / K_B
    z2 = sum(n * D.SPECIES[s][1] ** 2 for s, n in ions.items())
    return 2.06e-25 * t**-3.5 * (z2 * 1e-6) * (n_e * 1e-6) * 100.0


X_GRID = np.geomspace(1e-5, 40.0, 800)  # photon energy / kT
H_PLANCK = 6.62607015e-34


def escape_probability(tau) -> np.ndarray:
    """Escape probability of photons made uniformly in a sphere of radial optical depth tau (Osterbrock)."""
    tau = np.asarray(tau, dtype=float)
    t = np.maximum(tau, 1e-3)
    p = 3 / (4 * t) * (1 - 1 / (2 * t**2) + (1 / t + 1 / (2 * t**2)) * np.exp(-2 * t))
    return np.where(tau < 1e-3, 1 - 0.75 * tau, p)


def trapping(ions: dict, n_e: float, t_kev: float, radius_m: float) -> dict:
    """How much of a ball's bremsstrahlung escapes, and how much Planck radiation it traps (J/m^3)."""
    t_k = t_kev * KEV / K_B
    z2 = sum(n * D.SPECIES[s][1] ** 2 for s, n in ions.items())
    nu = X_GRID * t_kev * KEV / H_PLANCK
    alpha = 3.7e8 * t_k**-0.5 * (z2 * 1e-6) * (n_e * 1e-6) * nu**-3 * -np.expm1(-X_GRID) * 100.0  # 1/m
    tau_a = alpha * radius_m
    tau_s = n_e * SIGMA_T * radius_m
    # Scattering lengthens each photon's path: absorption depth ~ tau_a (1 + tau_s) while that is small, and
    # the thermalisation depth sqrt(tau_a (tau_a + tau_s)) once photons random-walk; the two meet at tau_a tau_s ~ 1.
    escape = escape_probability(np.minimum(tau_a * (1 + tau_s), np.sqrt(tau_a * (tau_a + tau_s))))
    emitted = np.exp(-X_GRID)
    escaped = float(np.trapezoid(emitted * escape, X_GRID) / np.trapezoid(emitted, X_GRID))
    planck = 15 / math.pi**4 * X_GRID**3 / np.expm1(X_GRID)
    trapped_planck = float(np.trapezoid(planck * (1 - escape), X_GRID))
    pair_photons = 1 - float(np.interp(P.M_E_KEV / t_kev, X_GRID, escape))
    return {"escaped": escaped, "radiation_J_m3": A_RAD * t_k**4 * trapped_planck,
            "pair_photons_trapped": pair_photons, "thomson_depth": tau_s}


def alpha_share_classical(ions: dict, n_e: float, t_e_kev: float, e0_kev: float = D.ALPHA_KEV["pB"]) -> float:
    """Share of an alpha's energy given to ions in a classical plasma (Stix): (1/y) Int_0^y du / (1 + u^1.5),
    y = E0 / E_c, E_c = 14.8 A_alpha T_e (sum n_j Z_j^2 / (n_e A_j))^(2/3)."""
    s = sum(n * D.SPECIES[sp][1] ** 2 / D.SPECIES[sp][0] for sp, n in ions.items()) / n_e
    y = e0_kev / (14.8 * D.SPECIES["He4"][0] * t_e_kev * s ** (2 / 3))
    u = np.linspace(0.0, y, 400)
    return float(np.trapezoid(1 / (1 + u**1.5), u) / y)


def opaque_ball(rho_g_cc: float, t_i_kev: float, radius_m: float, boron_fraction: float = 0.15,
                alpha_channelling: float = 0.0) -> dict:
    """Steady burn of a free (untamped) proton-boron ball, with its X-rays trapped if it is thick enough.

    The ions sit at t_i_kev. The electrons settle where what they receive
    (the ions' drain plus the alphas' share) equals what they lose: thin-plasma
    bremsstrahlung, less the part the ball reabsorbs (see trapping). Returns:
    * "trapped": the share of the bremsstrahlung that is reabsorbed ("opaque" when over half);
    * "burning": charged fusion power / radiation lost (>= 1: the ions keep themselves hot);
    * "gain": fusion energy released while the ball holds together (R / 4 c_s) / the energy it holds
      (ions, electrons, and for an opaque ball its radiation and thermal pairs);
    * the electron temperature, optical depth, energy, mass and pressure.
    alpha_channelling (untested theory, Fisch & Rax 1992): waves pass that share of what the alphas would give
    the electrons to the ions instead.
    """
    ions, n_e0 = D.fuel_ions(D.pb_fuel(boron_fraction), rho_g_cc)
    fermi = D.fermi_gas(n_e0)["E_F_keV"]
    rate = ions["p"] * ions["B11"] * EN.reactivity("p", "B11", t_i_kev)
    charged = rate * D.CHARGED_KEV["pB"] * KEV

    def state(t_e: float) -> dict:
        trap = trapping(ions, n_e0, t_e, radius_m)
        positrons = pair_density(n_e0, t_e, trap["pair_photons_trapped"])
        loss = D.bremsstrahlung_w_m3(n_e0, ions, t_e, fermi) * trap["escaped"]
        return {"trapped": 1 - trap["escaped"], "loss": loss, "n_e": n_e0 + 2 * positrons,
                "radiation": trap["radiation_J_m3"], "pairs": positrons}

    def excess(t_e: float) -> float:
        s = state(t_e)
        f_i = alpha_share_classical(ions, s["n_e"], t_e)
        f_i += alpha_channelling * (1 - f_i)  # untested: waves hand this share of the electrons' part to the ions
        return D._drain(ions, s["n_e"], t_i_kev, t_e) + (1 - f_i) * charged - s["loss"]

    lo, hi = 1e-2, 3 * t_i_kev
    if excess(hi) > 0:
        lo = hi
    while hi / lo > 1.001:
        mid = math.sqrt(lo * hi)
        lo, hi = (mid, hi) if excess(mid) > 0 else (lo, mid)
    t_e = math.sqrt(lo * hi)
    s = state(t_e)
    n_ions = sum(ions.values())
    eps = D.squeeze_energy_kev(n_e0)
    heat = 1.5 * n_ions * t_i_kev * KEV + n_e0 * math.hypot(eps, 1.5 * t_e) * KEV
    pairs_j = s["pairs"] * (2 * P.M_E_KEV + 3 * t_e) * KEV
    energy_density = heat + s["radiation"] + pairs_j
    pressure = n_ions * t_i_kev * KEV + s["n_e"] * t_e * KEV + s["radiation"] / 3
    sound = math.sqrt(5 / 3 * pressure / (rho_g_cc * 1e3))
    hold = radius_m / (4 * sound)
    volume = 4 / 3 * math.pi * radius_m**3
    fusion_density = rate * P.Q_VALUE_MEV["pB"] * 1e3 * KEV * hold
    return {"rho_g_cc": rho_g_cc, "T_i_keV": t_i_kev, "radius_m": radius_m, "T_e_keV": t_e,
            "trapped": s["trapped"], "opaque": s["trapped"] >= 0.5, "burning": charged / s["loss"],
            "gain": fusion_density / energy_density, "energy_J": energy_density * volume,
            "mass_kg": rho_g_cc * 1e3 * volume, "rhoR_g_cm2": rho_g_cc * radius_m * 100,
            "pressure_Pa": pressure, "radiation_share": (s["radiation"] + pairs_j) / energy_density,
            "positrons_per_electron": s["pairs"] / n_e0, "hold_s": hold}


def smallest_burning_ball(rho_g_cc: float, temperatures=(150.0, 200.0, 300.0, 500.0), boron_fraction: float = 0.15,
                          r_lo: float = 1e-7, r_hi: float = 30.0, points: int = 161) -> dict | None:
    """The proton-boron ball holding the least energy that keeps its ions hot and releases its own energy
    before flying apart (burning >= 1, gain >= 1). None if no ball up to r_hi does."""
    best = None
    for t_i in temperatures:
        for radius in np.geomspace(r_lo, r_hi, points):
            ball = opaque_ball(rho_g_cc, t_i, float(radius), boron_fraction)
            if ball["burning"] >= 1 and ball["gain"] >= 1:
                if best is None or ball["energy_J"] < best["energy_J"]:
                    best = ball
                break
    return best


def smallest_opaque_ball(rho_g_cc: float, t_i_kev: float = 300.0, boron_fraction: float = 0.15,
                         r_lo: float = 1e-7, r_hi: float = 1e3, points: int = 241) -> dict | None:
    """The smallest proton-boron ball at rho and t_i that traps at least half of its X-rays."""
    for radius in np.geomspace(r_lo, r_hi, points):
        ball = opaque_ball(rho_g_cc, t_i_kev, float(radius), boron_fraction)
        if ball["opaque"]:
            return ball
    return None


# --- 2. holding the fuel with a current: a fibre Z-pinch -------------------------------------

def alpha_confining_current() -> float:
    """Current (A) above which a pinch's own field holds D-T alphas within its radius: B(a) a = mu0 I / 2 pi."""
    return 2 * math.pi * ALPHA_LARMOR_TM / MU0


def bennett_temperature(current_a: float, fuel: str, radius_m: float) -> float:
    """Temperature (keV) at which a fibre of solid `fuel` (a key of energy.FUELS) of the given radius is in
    Bennett equilibrium with current_a (T_e = T_i)."""
    parts, units = EN.FUELS[fuel]
    ions_per_m = units * sum(parts.values()) * math.pi * radius_m**2
    charge = sum(f * EN.SPECIES[s][1] for s, f in parts.items()) / sum(parts.values())
    return MU0 * current_a**2 / (8 * math.pi * (1 + charge) * ions_per_m) / KEV


def fibre_pinch(current_a: float, t_kev: float = 10.0, alfven_times: float = 1.0, fuel: str = "DT",
                fusion_over_heat: float = 1.0, wall_m: float = 1e-2) -> dict:
    """A Bennett pinch of solid fuel at t_kev carrying current_a.

    The fibre's line density is chosen so that Bennett's relation gives
    t_kev. The pinch lives `alfven_times` x a / v_A. Fusion per metre in that
    time grows as 1/a (denser fuel burns faster), so there is a radius at which
    fusion equals `fusion_over_heat` x the plasma's heat. Returns that radius,
    the squeeze from the solid fibre, the density, the lifetime, the heat and
    the field energy per metre (out to a return conductor at wall_m).
    fusion_over_heat = 5 means the alphas alone resupply the heat (ignition-like).
    """
    parts, units = EN.FUELS[fuel]
    n_solid = units * sum(parts.values())
    charge = sum(f * EN.SPECIES[s][1] for s, f in parts.items()) / sum(parts.values())
    mass = sum(f * EN.SPECIES[s][0] for s, f in parts.items()) / sum(parts.values()) * P.AMU
    ions_per_m = MU0 * current_a**2 / (8 * math.pi * (1 + charge) * t_kev * KEV)
    a_solid = math.sqrt(ions_per_m / (math.pi * n_solid))
    v_alfven = current_a * math.sqrt(MU0 / (4 * math.pi * ions_per_m * mass))
    shares = {s: f / sum(parts.values()) for s, f in parts.items()}
    # fusion energy per metre per second at radius a = (sum over pairs of N_a N_b <sv> Q) / (pi a^2)
    rate_q = 0.0
    names = list(shares)
    for i, a in enumerate(names):
        for b in names[i:]:
            q = EN.q_value(a, b)
            if q == 0:
                continue
            pair = shares[a] * shares[b] * ions_per_m**2 * EN.reactivity(a, b, t_kev) * (0.5 if a == b else 1.0)
            rate_q += pair * q * 1e3 * KEV
    heat = 1.5 * (1 + charge) * ions_per_m * t_kev * KEV
    radius = rate_q * alfven_times / (math.pi * v_alfven * fusion_over_heat * heat)
    field = MU0 * current_a**2 / (4 * math.pi) * (math.log(max(wall_m / radius, 1.0)) + 0.25)
    return {"current_A": current_a, "T_keV": t_kev, "ions_per_m": ions_per_m, "fibre_radius_m": a_solid,
            "radius_m": radius, "squeeze": a_solid / radius,
            "rho_g_cc": ions_per_m * mass / (math.pi * radius**2) / 1e3, "lifetime_s": alfven_times * radius / v_alfven,
            "alfven_speed": v_alfven, "heat_J_m": heat, "field_J_m": field, "fusion_J_m": fusion_over_heat * heat,
            "alphas_held": current_a >= alpha_confining_current()}


# --- 3. steering the ions through a crystal: channelling ----------------------------------

def critical_angle(e_kev: float, z1: int, z2: int, spacing_m: float) -> float:
    """Lindhard's axial critical angle (rad): psi_1 = sqrt(2 Z1 Z2 e^2 / (E d))."""
    return math.sqrt(2 * z1 * z2 * D.E2 / (e_kev * KEV * spacing_m))


def axial_acceptance(psi: float, directions: int) -> float:
    """Share of an isotropic source within psi of one of `directions` axes (cones of half-angle psi)."""
    return min(1.0, directions * (1 - math.cos(psi)) / 2)


def planar_acceptance(psi: float, planes: int) -> float:
    """Share of an isotropic source within psi of one of `planes` crystal planes (bands of half-width psi
    around great circles; overlaps ignored, so generous)."""
    return min(1.0, planes * math.sin(psi))


def channelled_return(random_return: float, channels) -> float:
    """Fusion return of a crystal catcher when some ions channel.

    channels: iterable of (accepted, flux_peaking, stopping_ratio, kept). A
    share `accepted` of the ions channels; such an ion meets interstitial
    targets `flux_peaking` times as often and loses energy at `stopping_ratio`
    of the random rate, for the share `kept` of its path before it
    dechannels. The rest behave as in a random (amorphous) catcher.
    """
    total = 1.0
    for accepted, flux_peaking, stopping_ratio, kept in channels:
        total += accepted * kept * (flux_peaking / stopping_ratio - 1)
    return random_return * total


def magnetised_electron_field(n_e: float, t_e_kev: float) -> float:
    """Field (T) at which electrons gyrate inside a Debye length, the onset of magnetised (reduced) drag:
    m v_e / (e B) = lambda_D."""
    debye = math.sqrt(P.EPS0 * t_e_kev * KEV / (n_e * P.E_CHARGE**2))
    v_e = math.sqrt(t_e_kev * KEV / P.M_E)
    return P.M_E * v_e / (P.E_CHARGE * debye)
