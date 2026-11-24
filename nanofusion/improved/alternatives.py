"""Two "cold" fusion routes in crystals that need no hot plasma.

1. Muon-catalysed fusion in solid (or compressed) deuterium-tritium.
   A negative muon replaces the electron of a hydrogen atom. Its orbit is 207
   times smaller, so a d-t-muon molecule brings the nuclei close enough to fuse
   within a picosecond. The muon is then usually freed to catalyse again,
   until it decays (2.2 us) or sticks to the alpha particle.

   Fusions per muon (the standard cycle formula, e.g. Jones, Nature 321, 127
   (1986); arXiv:2605.26432):
       Y = 1 / (lambda_0 / lambda_c + omega_s)
   lambda_0 = 4.55e5 /s is the muon decay rate, lambda_c the cycle rate (about
   proportional to density), and omega_s the effective sticking per fusion.

2. Lattice confinement fusion (NASA Glenn; Pines et al. and Steinetz et al.,
   Phys. Rev. C 101, 044609 and 044610 (2020)).
   Gamma rays split deuterons in erbium or titanium deuteride. The
   photoneutrons knock other deuterons to tens of keV, which then fuse with
   cold, electron-screened deuterons in the lattice.
"""

from __future__ import annotations

import numpy as np

from . import energy as EN

# --- muon-catalysed fusion -----------------------------------------------------

MUON_DECAY_RATE = 1 / 2.197e-6  # 1/s
DT_FUSION_MEV = 17.59


def fusions_per_muon(cycle_rate: float, sticking: float) -> float:
    """Mean D-T fusions catalysed by one muon."""
    return 1.0 / (MUON_DECAY_RATE / cycle_rate + sticking)


def muon_gain(fusions: float, muon_cost_gev: float, energy_per_fusion_mev: float = DT_FUSION_MEV) -> float:
    """Fusion energy per muon divided by the energy it took to make the muon."""
    return fusions * energy_per_fusion_mev / (muon_cost_gev * 1e3)


def fusions_needed(muon_cost_gev: float, energy_per_fusion_mev: float = DT_FUSION_MEV) -> float:
    return muon_cost_gev * 1e3 / energy_per_fusion_mev


# --- lattice confinement fusion -----------------------------------------------------

DEUTERON_BINDING_MEV = 2.2246
DD_FUSION_MEV = 3.65  # mean of the two branches


def photoneutron_energies(endpoint_mev: float = 2.9, n: int = 20000, seed: int = 1) -> np.ndarray:
    """Photoneutron energies (keV) from a bremsstrahlung beam just above the D(gamma,n)p threshold.

    The neutron gets about half of (E_gamma - 2.2246 MeV). Taking the photon
    spectrum to fall linearly to its endpoint and the cross section to rise as
    the square root of the excess energy gives p(E_n) ~ E_n^0.5 (E_max - E_n).
    For a 2.9 MeV endpoint the mean is 145 keV, which is the value NASA
    calculated for its tests.
    """
    e_max = (endpoint_mev - DEUTERON_BINDING_MEV) / 2 * 1e3
    gen = np.random.default_rng(seed)
    out = np.empty(0)
    while len(out) < n:  # rejection sampling
        x = gen.uniform(0, e_max, 4 * n)
        keep = gen.uniform(0, 1, 4 * n) * (e_max**1.5 * 0.385) < x**0.5 * (e_max - x)
        out = np.concatenate([out, x[keep]])
    return out[:n]


def lattice_fusions_per_neutron(neutrons_kev: np.ndarray, material: str, d_density: float, screening_ev: float,
                                collisions: float | None = None, seed: int = 2) -> dict:
    """D-D fusions per photoneutron in a metal deuteride.

    Each collision with a deuteron gives it an energy T uniform in [0, 8/9 E_n]
    (isotropic elastic scattering in the centre-of-mass frame). The knocked-on
    deuteron then fuses with cold, screened lattice deuterons as it stops
    (thick-target yield with the cross section at E + U_e).

    collisions=None follows each neutron down to 1 keV in an infinite sample
    (every possible collision with a deuteron, the best case). A number gives
    the mean count of deuteron collisions per neutron in a finite sample, and
    only first collisions are counted.
    """
    gen = np.random.default_rng(seed)
    grid = np.geomspace(0.1, max(float(neutrons_kev.max()), 1.0), 400)
    y_grid = EN.thick_target(grid, "D", material, {"D": d_density}, screening_ev=screening_ev)["reactions"]

    def yield_of(t):
        return np.interp(t, grid, y_grid, left=0.0)

    if collisions is not None:
        t = gen.uniform(0, 8 / 9, len(neutrons_kev)) * neutrons_kev
        fusions = collisions * float(np.mean(yield_of(t)))
        return {"fusions_per_neutron": fusions, "deuteron_hits_per_neutron": collisions,
                "mean_deuteron_keV": float(np.mean(t))}
    e = neutrons_kev.copy()
    total = np.zeros_like(e)
    hits = np.zeros_like(e)
    alive = e > 1.0
    while alive.any():
        t = gen.uniform(0, 8 / 9, alive.sum()) * e[alive]
        total[alive] += yield_of(t)
        hits[alive] += 1
        e[alive] -= t
        alive = e > 1.0
    return {"fusions_per_neutron": float(total.mean()), "deuteron_hits_per_neutron": float(hits.mean()),
            "mean_deuteron_keV": None}


# --- squeezing a cold crystal until it fuses (pycnonuclear fusion) --------------------

HBAR_C_MEV_FM = 197.3269804
M_E_MEV = 0.51099895
S_ZERO_MEV_B = {"DD": 0.11, "DT": 11.5}  # astrophysical S-factors near zero energy (Bosch & Hale fits)


def pycnonuclear_rate(n_cm3: float, z: int, a: float, s_mev_b: float) -> float:
    """Zero-temperature reactions per cm^3 per s in a bcc lattice of identical nuclei.

    Salpeter & Van Horn, ApJ 155, 183 (1969), static lattice:
        R = 3.90e46 Z^4 A rho S lambda^(7/4) exp(-2.638 / sqrt(lambda)),
        lambda = (rho / 1.3574e11 g/cm^3)^(1/3) / (Z^2 A^(4/3)),
    with rho = n A m_u in g/cm^3 and S in MeV b. Model uncertainty is
    orders of magnitude for light nuclei, whose zero-point motion is large.
    """
    rho = n_cm3 * a * 1.66053906660e-24
    lam = (rho / 1.3574e11) ** (1 / 3) / (z**2 * a ** (4 / 3))
    return 3.90e46 * z**4 * a * rho * s_mev_b * lam**1.75 * np.exp(-2.638 / np.sqrt(lam))


def pycnonuclear_rate_per_nucleus(rho_g_cc: float, fuel: str) -> float:
    """Reactions per second per nucleus in cold D (fuel 'DD') or equimolar D-T (fuel 'DT').

    D-T is treated as one component with the D-T reduced mass (A = 2 mu = 2.4),
    with half the pairs being D-T.
    """
    mean_a = 2.014 if fuel == "DD" else 2.515
    n = rho_g_cc / (mean_a * 1.66053906660e-24)
    if fuel == "DD":
        return pycnonuclear_rate(n, 1, 2.014, S_ZERO_MEV_B["DD"]) / n
    a_eff = 2 * 2.014 * 3.016 / (2.014 + 3.016)
    return 0.5 * pycnonuclear_rate(n, 1, a_eff, S_ZERO_MEV_B["DT"]) / n


def degenerate_electron_energy_kev(n_e_m3: float) -> float:
    """Mean kinetic energy per electron (keV) of a zero-temperature Fermi gas, relativistic."""
    p_f = HBAR_C_MEV_FM * 1e-15 * (3 * np.pi**2 * n_e_m3) ** (1 / 3)  # MeV (as pc)
    x = p_f / M_E_MEV
    if x < 1e-3:
        return 0.6 * M_E_MEV * x**2 / 2 * 1e3
    total = 3 * M_E_MEV / (8 * x**3) * (x * (2 * x**2 + 1) * np.sqrt(1 + x**2) - np.arcsinh(x))
    return (total - M_E_MEV) * 1e3


def cold_compression(rho_g_cc: float, fuel: str = "DT", ignition_burn: float = 0.01) -> dict:
    """What it takes to fuse a cold crystal of fuel by compression alone.

    The fuel is compressed without heating, so its pressure is that of the
    degenerate electrons. It flies apart in tau = R / (4 c_s). Ignition: the
    pycnonuclear burn must consume `ignition_burn` of the fuel within tau; the
    alpha particles then heat it to thermonuclear temperatures (0.4% would
    heat the ions to 5 keV). The fuel radius R that achieves this sets the
    mass and the (ideal, shock-free) compression energy.
    """
    mean_a = 2.014 if fuel == "DD" else 2.515
    n = rho_g_cc * 1e3 / (mean_a * 1.66053906660e-27)  # m^-3, one electron per nucleus
    eps = degenerate_electron_energy_kev(n)
    h = 1e-3
    pressure = n**2 * (degenerate_electron_energy_kev(n * (1 + h)) - degenerate_electron_energy_kev(n * (1 - h))) \
        * 1.602176634e-16 / (2 * h * n)
    gamma = np.log(_pressure(n * (1 + h)) / _pressure(n * (1 - h))) / np.log((1 + h) / (1 - h))
    sound = np.sqrt(gamma * pressure / (rho_g_cc * 1e3))
    rate = pycnonuclear_rate_per_nucleus(rho_g_cc, fuel)
    radius = ignition_burn * 4 * sound / rate if rate > 0 else np.inf
    radius = min(radius, 1e30)  # beyond any meaning; keeps the arithmetic finite
    mass = 4 / 3 * np.pi * radius**3 * rho_g_cc * 1e3
    n_nuclei = mass / (mean_a * 1.66053906660e-27)
    energy = n_nuclei * eps * 1.602176634e-16
    q = DT_FUSION_MEV if fuel == "DT" else DD_FUSION_MEV
    return {"rho_g_cc": rho_g_cc, "fuel": fuel, "rate_per_nucleus_s": rate, "electron_energy_keV": eps,
            "pressure_Gbar": pressure / 1e14, "sound_speed_m_s": sound,
            "convergence_from_solid": (rho_g_cc / (0.21 if fuel == "DT" else 0.2)) ** (1 / 3),
            "ignition_radius_m": radius, "mass_kg": mass, "compression_J": energy,
            "yield_if_burnt_J": n_nuclei / 2 * q * 1.602176634e-13,
            "gain_ideal": q * 1e3 / (2 * eps) if eps > 0 else np.inf}


def _pressure(n: float) -> float:
    h = 1e-4
    return n**2 * (degenerate_electron_energy_kev(n * (1 + h)) - degenerate_electron_energy_kev(n * (1 - h))) / (2 * h * n)


# --- exploding crystals as a muon factory -------------------------------------------

PION_THRESHOLD_GEV_U = 0.3  # little pion production below about 0.3 GeV per nucleon
PION_SATURATION_GEV_U = 1.0  # muon production per GeV of beam is best from about 1 GeV per nucleon


def crystal_muon_factory(e_max_gev_u: float, best_cost_gev: float, laser_to_ions: float) -> dict:
    """Laser energy per muon if exploding 6LiD crystals drive the pion target.

    Every ion of a Coulomb-exploded 6LiD crystal has the same energy per
    nucleon at the same radius (D and 6Li3+ share one charge-to-mass ratio),
    rising as r^2 to e_max_gev_u, so nucleon energies are spread as
    dN/dE ~ sqrt(E). Muons per GeV of beam rise linearly from the pion
    threshold to their best value (1 / best_cost_gev) at 1 GeV per nucleon.
    The source is isotropic, which is fine here: the pion target can surround it.
    """
    e = np.linspace(e_max_gev_u / 2000, e_max_gev_u, 2000)
    weight = np.sqrt(e)
    efficiency = np.clip((e - PION_THRESHOLD_GEV_U) / (PION_SATURATION_GEV_U - PION_THRESHOLD_GEV_U), 0, 1) / best_cost_gev
    beam = float(np.sum(weight * e))
    muons = float(np.sum(weight * e * efficiency))
    beam_per_muon = beam / muons if muons > 0 else np.inf
    rho_q = EN.charge_density("6LiD")
    per_nucleon_j = 0.9 * 1.602176634e-19 * rho_q / (3 * 8.8541878128e-12) / 2.014  # J per m^2 of R^2, per u
    radius = np.sqrt(e_max_gev_u * 1.602176634e-10 / per_nucleon_j)
    field = rho_q * radius / (3 * 8.8541878128e-12)
    intensity = 0.5 * 2.99792458e8 * 8.8541878128e-12 * field**2 / 1e4
    return {"e_max_GeV_u": e_max_gev_u, "crystal_radius_m": radius, "strip_intensity_W_cm2": intensity,
            "beam_GeV_per_muon": beam_per_muon, "laser_GeV_per_muon": beam_per_muon / laser_to_ions}


# --- fusion in dense matter, all regimes (Gasques et al., Phys. Rev. C 72, 025806 (2005)) -------------

E2_EV_CM = 1.439964547e-7
HBAR_C_EV_CM = 1.973269804e-5
HBAR_ERG_S = 1.054571817e-27
E2_ERG_CM = 2.307077552e-19
AMU_G = 1.66053906660e-24
AMU_EV = 931.49410242e6
C_SC = 1.0754
DENSE_MODELS = {
    # name: (C_exp, C_pyc, C_pl, C_T, Lambda)  -- Table I of Gasques et al.
    "optimal (static bcc lattice)": (2.638, 3.90, 1.25, 0.724, 0.5),
    "maximum rate": (2.450, 50.0, 1.25, 0.904, 0.35),
    "minimum rate": (2.650, 0.5, 1.25, 0.711, 0.65),
    "Ogata-Iyetomi-Ichimaru (Monte Carlo)": (2.460, 0.00181, 1.809, 0.893, 0.5),
}


def _pair(fuel: str) -> tuple[float, float, float]:
    """(barrier mass A, mean ion mass, fraction of pairs that react) for 'DD' or equimolar 'DT'.

    For D-T the barrier mass is twice the reduced mass, 2.41 u, which is also
    the effective mass of the mixture's ion plasma frequency.
    """
    if fuel == "DD":
        return 2.0141, 2.0141, 1.0
    return 2 * 2.0141 * 3.0160 / (2.0141 + 3.0160), 2.515, 0.5


def dense_fusion_rate(rho_g_cc: float, t_ev: float, fuel: str = "DT",
                      model: str = "optimal (static bcc lattice)") -> float:
    """Reactions per nucleus per second in dense fuel, in all five burning regimes.

    Gasques et al. 2005, Eq. (28): R = R_pyc(T=0) + dR(T), interpolating
    between classical thermonuclear burning, strong plasma screening (the
    Monte Carlo h0 ~ 1.0754 Gamma), thermally enhanced and zero-temperature
    pycnonuclear burning. The fuel is treated as one component with the pair's
    barrier mass. The constants were tuned for carbon, so for hydrogen isotopes
    take the spread between models as a lower bound on the uncertainty.
    """
    c_exp, c_pyc, c_pl, c_t, lam_c = DENSE_MODELS[model]
    a_mass, mean_mass, pair_fraction = _pair(fuel)
    s_key = "DD" if fuel == "DD" else "DT"
    n = rho_g_cc / (mean_mass * AMU_G)  # cm^-3
    a = (3 / (4 * np.pi * n)) ** (1 / 3)
    mc2 = a_mass * AMU_EV
    t_p = np.sqrt(4 * np.pi * n * E2_EV_CM * HBAR_C_EV_CM**2 / mc2)  # eV
    e_a = mc2 * (E2_EV_CM / HBAR_C_EV_CM) ** 2  # eV
    # zero-temperature pycnonuclear part, Eq. (25), with rho for the barrier mass
    rho_b = n * a_mass * AMU_G
    lam = (rho_b / 1.3574e11) ** (1 / 3) / a_mass ** (4 / 3)
    r_pyc = rho_b * a_mass * _s_factor(s_key, t_p) * c_pyc * 1e46 * lam ** (3 - c_pl) * np.exp(-c_exp / np.sqrt(lam))
    # temperature-dependent part
    if t_ev > 0:
        t_tilde = np.sqrt(t_ev**2 + c_t**2 * t_p**2)
        gamma_t = E2_EV_CM / (a * t_tilde)
        tau_t = 3 * (np.pi / 2) ** (2 / 3) * (e_a / t_tilde) ** (1 / 3)
        phi = np.sqrt(gamma_t) / ((C_SC**4 / 9) + gamma_t**2) ** 0.25
        boltz = np.exp(-lam_c * t_p / t_ev)
        exponent = -tau_t + C_SC * gamma_t * phi * boltz - lam_c * t_p / t_ev
        g = (t_ev**2 * (2 / 3) + t_p**2 * (2 / 3) * (c_pl + 0.5)) / (t_ev**2 + t_p**2)
        pre = 8 * np.pi ** (1 / 3) / (np.sqrt(3) * 2 ** (1 / 3)) * (e_a / t_tilde) ** g
        tau = 3 * (np.pi / 2) ** (2 / 3) * (e_a / t_ev) ** (1 / 3)
        e_pk = t_p + (E2_EV_CM / a + t_ev * tau / 3) * boltz
        s_erg_cm2 = _s_factor(s_key, e_pk) * 1.602176634e-6 * 1e-24
        d_r = n**2 / 2 * s_erg_cm2 * HBAR_ERG_S / (a_mass * AMU_G * E2_ERG_CM) * pre * np.exp(exponent)
    else:
        d_r = 0.0
    return pair_fraction * (r_pyc + d_r) / n


def _s_factor(key: str, e_ev: float) -> float:
    """Astrophysical S-factor (MeV b) at a centre-of-mass energy in eV, from the Bosch-Hale cross sections."""
    from . import physics as P

    e_kev = float(np.clip(e_ev / 1e3, 0.5, 500.0))
    if key == "DD":
        sigma_b = float(P.sigma_dd_mb(e_kev)) * 1e-3
        b_g = P.gamow_constant(1, 1, P.reduced_mass_u(P.MASS_U["D"], P.MASS_U["D"]))
    else:
        sigma_b = float(P.sigma_dt_mb(e_kev)) * 1e-3
        b_g = P.gamow_constant(1, 1, P.reduced_mass_u(P.MASS_U["D"], P.MASS_U["T"]))
    return sigma_b * e_kev * np.exp(b_g / np.sqrt(e_kev)) / 1e3


def ignition_ball(rho_g_cc: float, t_ev: float, fuel: str = "DT", model: str = "optimal (static bcc lattice)",
                  burn: float = 0.01, alpha_rho_r: float = 0.3) -> dict:
    """Smallest ball of fuel at (rho, T) that ignites, and the energy it holds.

    Two conditions: it burns `burn` of itself before flying apart,
    rate x R / (4 c_s) >= burn (as in cold_compression), and it is thick
    enough to stop its own alpha particles, rho R >= alpha_rho_r g/cm^2 (the
    usual ~0.3 g/cm^2 of inertial fusion). Energy = electron energy
    (degenerate plus thermal, joined smoothly) + ion thermal energy, per
    nucleus, times the number of nuclei. It excludes the driver's losses.
    """
    _, mean_mass, _ = _pair(fuel)
    n_m3 = rho_g_cc * 1e3 / (mean_mass * 1.66053906660e-27)
    eps_deg = degenerate_electron_energy_kev(n_m3) * 1e3  # eV
    eps_e = np.hypot(eps_deg, 1.5 * t_ev)
    p_deg = _pressure(n_m3) * 1e3  # eV per m^3
    pressure = (np.hypot(p_deg, n_m3 * t_ev) + n_m3 * t_ev) * 1.602176634e-19  # Pa
    sound = np.sqrt(5 / 3 * pressure / (rho_g_cc * 1e3))
    rate = dense_fusion_rate(rho_g_cc, t_ev, fuel, model)
    radius = min(burn * 4 * sound / rate, 1e30) if rate > 0 else 1e30
    radius = max(radius, alpha_rho_r / (rho_g_cc * 100))
    nuclei = 4 / 3 * np.pi * radius**3 * n_m3
    energy = nuclei * (eps_e + 1.5 * t_ev) * 1.602176634e-19
    gamma = E2_EV_CM / ((3 / (4 * np.pi * n_m3 * 1e-6)) ** (1 / 3) * t_ev) if t_ev > 0 else np.inf
    return {"rho_g_cc": rho_g_cc, "T_eV": t_ev, "rate_per_nucleus_s": rate, "radius_m": radius,
            "mass_kg": nuclei * mean_mass * 1.66053906660e-27, "energy_J": energy,
            "rhoR_g_cm2": rho_g_cc * radius * 100, "coupling_Gamma": gamma}


# --- a fission blanket around a fusion neutron source (fusion-fission hybrid) ------------------

FISSION_MEV = 200.0  # energy per fission
NEUTRONS_PER_FISSION = 2.5
BREEDING_CAPTURE_MEV = 4.78  # 6Li(n,t)4He, the capture that replaces the burnt triton
# A 14 MeV neutron in thick uranium splits uranium-238 directly and multiplies through (n,2n) and (n,3n).
# Measured in a natural-uranium pile (McEnhill & Weale 1962): 1.18 fissions of U-238 and 0.28 of U-235 per
# neutron, about 290 MeV. Calculated for pure U-238 (J. D. Lee): 233 MeV and 4.4 captures.
FAST_FISSIONS = 1.2  # direct fissions per 14 MeV neutron
NEUTRON_WORTH = 3.0  # chain-starting neutrons one 14 MeV neutron turns into (Bethe 1979: "two to four")


def blanket_energy_per_fusion_mev(k_eff: float, uranium: bool = True, nu: float = NEUTRONS_PER_FISSION) -> float:
    """Heat released per D-T fusion when its neutron enters a blanket.

    uranium=False: lithium only, the fusion energy plus one tritium-breeding
    capture (about 22 MeV). uranium=True: the neutron also causes
    FAST_FISSIONS fissions directly and leaves NEUTRON_WORTH neutrons, each
    starting chains that produce k/(1-k) further neutrons in a blanket of
    multiplication k_eff, so k / (nu (1 - k)) fissions. k_eff = 0 is
    depleted uranium with no chain: about 260 MeV per fusion, close to the
    measured 290 MeV of a natural-uranium pile.
    """
    if not 0 <= k_eff < 1:
        raise ValueError("k_eff must be in [0, 1)")
    heat = DT_FUSION_MEV + BREEDING_CAPTURE_MEV
    if uranium:
        heat += FISSION_MEV * (FAST_FISSIONS + NEUTRON_WORTH * k_eff / (nu * (1 - k_eff)))
    return heat


def hybrid_plant_gain(target_gain: float, driver_efficiency: float, k_eff: float, uranium: bool = True,
                      thermal_efficiency: float = 0.4) -> float:
    """Electricity out / electricity in for a fusion source inside a blanket.

    target_gain is fusion energy per unit of driver energy (laser light or
    beam) and driver_efficiency is driver energy per unit of electricity.
    """
    heat = blanket_energy_per_fusion_mev(k_eff, uranium)
    return thermal_efficiency * driver_efficiency * target_gain * heat / DT_FUSION_MEV


def hybrid_k_needed(target_gain: float, driver_efficiency: float, plant_gain: float = 1.0,
                    thermal_efficiency: float = 0.4) -> float | None:
    """Uranium-blanket k_eff at which the plant's electricity out / in reaches `plant_gain`.

    0.0 means uranium with no chain reaction is enough; None means a lithium
    blanket with no uranium at all is enough.
    """
    per_fusion = plant_gain * DT_FUSION_MEV / (thermal_efficiency * driver_efficiency * target_gain)
    if per_fusion <= blanket_energy_per_fusion_mev(0.0, uranium=False):
        return None
    x = (per_fusion - blanket_energy_per_fusion_mev(0.0)) * NEUTRONS_PER_FISSION / (FISSION_MEV * NEUTRON_WORTH)
    return x / (1 + x) if x > 0 else 0.0  # x = k / (1 - k)


def fission_share(k_eff: float) -> float:
    """Share of a uranium blanket's heat that comes from fission."""
    total = blanket_energy_per_fusion_mev(k_eff)
    return (total - DT_FUSION_MEV - BREEDING_CAPTURE_MEV) / total


# --- bare ions in a magnetic trap (no electrons to drag on them or to radiate) ------------------

def brillouin_density(b_tesla: float, mass_u: float = 2.515) -> float:
    """Highest density (m^-3) of a single-sign plasma that a field can hold: n = eps0 B^2 / (2 m)."""
    return 8.8541878128e-12 * b_tesla**2 / (2 * mass_u * 1.66053906660e-27)


def ion_trap(b_tesla: float, radius_m: float, length_m: float, t_kev: float = 60.0) -> dict:
    """D-T fusion from a column of bare ions held at the Brillouin density limit.

    With no electrons there is no electron drag and no bremsstrahlung, and
    like-particle collisions do not carry ions across the field, so on paper
    the energy balance is unlimited. The price is density: the ions' own
    charge caps it at the Brillouin limit, and the voltage across the column
    grows as n a^2.
    """
    n = brillouin_density(b_tesla)
    volume = np.pi * radius_m**2 * length_m
    power_density = 0.25 * n * n * EN.reactivity("D", "T", t_kev) * DT_FUSION_MEV * 1.602176634e-13
    thermal = 1.5 * n * t_kev * 1.602176634e-16
    field_energy = b_tesla**2 / (2 * 4e-7 * np.pi)
    return {"b_tesla": b_tesla, "density_m3": n, "power_density_W_m3": power_density,
            "power_W": power_density * volume,
            "voltage_V": 1.602176634e-19 * n * radius_m**2 / (4 * 8.8541878128e-12),
            "thermal_J": thermal * volume, "field_J": field_energy * volume,
            "confinement_for_breakeven_s": thermal / power_density,
            "field_payback_s": field_energy / power_density}


# --- throwing whole crystals at a target (cluster-impact fusion) ------------------------------

# A positively charged grain starts to shed ions above a surface field of about 10 V/nm (the field ion emission
# limit used for charged dust, Dzhanoev et al., A&A 591, A147 (2016)); metals hold 19-57 V/nm, hydrogen ice far less.
FIELD_EVAPORATION_V_M = 1.0e10


def crystal_projectile(radius_nm: float, kev_per_deuteron: float, fuel: str = "D2",
                       surface_field: float = FIELD_EVAPORATION_V_M) -> dict:
    """A fuel crystal charged to the field-evaporation limit and accelerated as one projectile.

    Returns its charge, the voltage needed for the given energy per deuteron
    (every nucleus moves at the same speed, so heavier ones carry
    proportionally more energy), and the speed.
    """
    r = radius_nm * 1e-9
    parts, n_units = EN.FUELS[fuel]
    units = n_units * 4 / 3 * np.pi * r**3
    mass_u = sum(f * EN.SPECIES[s][0] for s, f in parts.items())
    charge_e = 4 * np.pi * 8.8541878128e-12 * r * r * surface_field / 1.602176634e-19
    kinetic_kev = units * mass_u / EN.SPECIES["D"][0] * kev_per_deuteron
    speed = np.sqrt(2 * kev_per_deuteron * 1.602176634e-16 / (EN.SPECIES["D"][0] * 1.66053906660e-27))
    return {"nuclei": units * sum(parts.values()), "charge_e": charge_e,
            "charge_per_mass_e_u": charge_e / (units * mass_u),
            "voltage_V": kinetic_kev * 1e3 / charge_e, "speed_m_s": float(speed),
            "kinetic_J": kinetic_kev * 1.602176634e-16}


def cluster_impact_yield(molecules: int, kev: float, molecule_mass_u: float = 20.03, deuterons_per_molecule: int = 2,
                         target: str = "TiD2", d_density: float = 1.0e29) -> dict:
    """D-D fusions per cluster for a (D2O)_n cluster of total energy `kev` stopping in a deuterated target.

    Each deuteron carries its mass share of the energy and is treated as an
    independent beam ion (thick-target yield), which is what the molecular
    binding energy of a few eV against keV of motion implies.
    """
    e_d = kev * EN.SPECIES["D"][0] / (molecules * molecule_mass_u)
    y = float(EN.thick_target([e_d], "D", target, {"D": d_density})["reactions"][0])
    return {"keV_per_deuteron": e_d, "fusions_per_deuteron": y, "fusions_per_cluster": y * molecules * deuterons_per_molecule}


# --- ions with a spread of energies stopping in a cold target (nanowire-array experiments) -------------

DEUTERON_DENSITY = {"CD2": 8.0e28}  # deuterons per m^3


def spectrum_yield(temperature_kev: float, e_max_kev: float, e_min_kev: float = 0.0, projectile: str = "D",
                   material: str = "CD2", targets: dict | None = None) -> dict:
    """Reactions per ion for an exponential spectrum, dN/dE ~ exp(-E / T) from e_min to e_max, stopping in a cold solid.

    Laser-accelerated ions usually have such a spectrum. Default target:
    deuterated polyethylene, as in the nanowire-array experiments.
    """
    targets = {"D": DEUTERON_DENSITY[material]} if targets is None else targets
    e = np.linspace(max(e_min_kev, e_max_kev / 4000), e_max_kev, 4000)
    w = np.exp(-(e - e[0]) / temperature_kev)
    w /= w.sum()
    res = EN.thick_target(e, projectile, material, targets)
    return {"mean_keV": float(np.sum(w * e)), "reactions_per_ion": float(np.sum(w * res["reactions"])),
            "fusion_MeV_per_ion": float(np.sum(w * res["fusion_MeV"]))}
