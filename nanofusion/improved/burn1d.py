"""One-dimensional implosion and burn of a D-T capsule: what the static model leaves out.

The simulator's compressed_icf route (simulator.py) is static. It takes an
isobaric hot spot at 10 keV inside a cold shell and burns the share
rhoR / (rhoR + 7 g/cm2). This module follows the capsule in time instead:

* The D-T shell, with the ablator still left on it, coasts inward at the
  implosion velocity.
* It stagnates on its vapour, which it heats by compression. Heat conducted
  out of the vapour ablates the shell's inner surface into the hot spot.
* The hot spot ignites (or doesn't) by alpha heating. The burn spreads into the
  dense shell, and the capsule blows itself apart.

What it adds to the static model:
* alpha heating and burn propagation, including alphas that escape into the
  shell;
* conduction losses and bremsstrahlung;
* disassembly: how long the fuel stays together;
* the cost of pressure: the hot spot's pressure comes from the shell's kinetic
  energy, which the laser must supply;
* mixing: carbon from the ablator mixed into the fuel that becomes the hot spot
  (it radiates and dilutes);
* asymmetry: a share of the kinetic energy that never reaches the hot spot
  (residual kinetic energy), modelled in 1-D as a slower shell.

Physics
-------
* Lagrangian hydrodynamics in spherical geometry, with von Neumann-Richtmyer
  artificial viscosity.
* One temperature for ions and electrons.
* Equation of state: ideal ions plus electrons whose pressure goes from ideal
  to degenerate, P = n_i T + n_e sqrt(T^2 + (2/5 E_F)^2), and internal energy
  3/2 P V.
* D-T reactivity of Bosch & Hale (1992). The fuel burns up into helium ash.
* Alphas: deposited along straight radial paths, in units of their range in
  each zone. The range uses Spitzer slowing (the same as confinement.py), with
  the electrons' mean energy floored by degeneracy. Inward-going alphas cross
  the centre and continue outward.
* Electron heat conduction: Spitzer-Härm, flux-limited (f = 0.1), implicit.
* Bremsstrahlung escapes from zones above ~1 keV. Cooler, dense zones are
  taken to be opaque, so their X-rays are re-absorbed where they are made.

Not included: two temperatures, radiation transport, the ablator's partial
ionisation and tungsten dopant, the drive phase (the run starts at peak
velocity), and real 2-D/3-D asymmetry.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field

import numpy as np
from scipy.linalg import solve_banded

from . import physics as P
from .confinement import ALPHA_V0, E_ALPHA_KEV

KEV = P.J_PER_KEV
AMU = P.AMU
M_E = P.M_E
HBAR = 1.054571817e-34
A_DT, A_HE, A_C = 2.5, 4.0026, 12.011
E_FUSION_J = 17.59e3 * KEV
E_ALPHA_J = E_ALPHA_KEV * KEV
BREMS = 5.35e-37  # W m^3 keV^-1/2, per n_e * sum(n Z^2)


@dataclass
class Capsule:
    """A capsule at peak implosion velocity, coasting inward."""

    fuel_mg: float = 0.22  # D-T ice
    velocity_kms: float = 394.0  # implosion velocity (N221204, LLNL post-shot simulation)
    adiabat: float = 2.45  # fuel pressure / degenerate pressure
    ablator_mg: float = 0.24  # ablator left on the fuel (carbon): 5.7% of N221204's 4.25 mg
    # where the run starts, and the fuel density then: not published; chosen so that the clean yields match
    # LLNL's 1-D simulations (N210808: 4.2-6.2 MJ, here 5.1; N221204-class: ~15 MJ, here 20)
    inner_radius_um: float = 300.0
    flight_density: float = 4.5  # g/cm3
    vapour_ug: float = 1.8  # D-T gas inside: 0.45 mg/cm3 in the 0.985 mm hollow of a NIF capsule
    mix_ng: float = 0.0  # carbon mixed into the vapour and the inner fuel
    mix_depth: float = 0.1  # share of the fuel mass (from the inside) the mix is spread through
    rke: float = 0.0  # share of the kinetic energy that never reaches the hot spot
    outer_pressure: float = 1.0  # ablation pressure left on the outside, in units of the shell's starting pressure
    alpha_heating: bool = True
    zones: tuple = (20, 110, 16)  # vapour, fuel, ablator

    @property
    def kinetic_J(self) -> float:
        return 0.5 * (self.fuel_mg + self.ablator_mg) * 1e-6 * (self.velocity_kms * 1e3) ** 2


@dataclass
class State:
    r: np.ndarray  # node radii (m)
    u: np.ndarray  # node velocities (m/s)
    m: np.ndarray  # zone masses (kg)
    e: np.ndarray  # zone specific internal energy (J/kg)
    fuel: np.ndarray  # D + T ions per zone
    ash: np.ndarray  # helium ions per zone
    carbon: np.ndarray  # carbon ions per zone
    history: dict = field(default_factory=dict)


# --- material properties ---------------------------------------------------------------------------------------------

def bosch_hale_dt(t_kev: np.ndarray) -> np.ndarray:
    """Maxwellian D-T <sigma v> (m^3/s), Bosch & Hale, Nucl. Fusion 32, 611 (1992); 0.2-100 keV."""
    t = np.clip(t_kev, 0.2, 100.0)
    c1, c2, c3, c4, c5, c6, c7 = 1.17302e-9, 1.51361e-2, 7.51886e-2, 4.60643e-3, 1.35e-2, -1.0675e-4, 1.366e-5
    theta = t / (1 - t * (c2 + t * (c4 + t * c6)) / (1 + t * (c3 + t * (c5 + t * c7))))
    xi = (34.3827**2 / (4 * theta)) ** (1 / 3)
    sv = c1 * theta * np.sqrt(xi / (1124656.0 * t**3)) * np.exp(-3 * xi) * 1e-6
    return np.where(t_kev < 0.2, 0.0, sv)


def fermi_energy(n_e: np.ndarray) -> np.ndarray:
    return HBAR**2 * (3 * math.pi**2 * n_e) ** (2 / 3) / (2 * M_E)


def temperature(pressure: np.ndarray, n_i: np.ndarray, n_e: np.ndarray) -> np.ndarray:
    """Invert P = n_i T + n_e sqrt(T^2 + c^2) for T (J); c = 2/5 E_F."""
    c = 0.4 * fermi_energy(n_e)
    a = n_i**2 - n_e**2
    b = -2 * pressure * n_i
    k = pressure**2 - (n_e * c) ** 2
    lin = k / np.maximum(-b, 1e-300)
    disc = np.sqrt(np.maximum(b * b - 4 * a * k, 0.0))
    quad = np.where(np.abs(a) > 1e-9 * n_e**2, 2 * k / np.maximum(-b + disc, 1e-300), lin)
    return np.maximum(np.where(np.abs(a) > 1e-9 * n_e**2, quad, lin), 0.0)


def pressure_of(t: np.ndarray, n_i: np.ndarray, n_e: np.ndarray) -> np.ndarray:
    return n_i * t + n_e * np.sqrt(t * t + (0.4 * fermi_energy(n_e)) ** 2)


def alpha_range_factor(u_c: np.ndarray) -> np.ndarray:
    """Vectorised confinement.alpha_range_factor."""
    a = np.maximum(u_c, 1e-9)
    log = np.log((1 + a) ** 2 / (1 - a + a * a)) / 6
    atan = (np.arctan((2 - a) / (a * math.sqrt(3))) + math.pi / 6) / math.sqrt(3)
    return 1.0 - a * (log + atan)


def coulomb_log(n_e: np.ndarray, t_kev: np.ndarray) -> np.ndarray:
    return np.maximum(2.0, 24.0 - np.log(np.sqrt(n_e * 1e-6) / (t_kev * 1e3)))


# Bremsstrahlung in four photon groups, x = h nu / T: share of the emission, and the group's mean x.
GROUPS = ((1 - math.exp(-1), 0.418), (math.exp(-1) - math.exp(-2), 1.418), (math.exp(-2) - math.exp(-4), 2.69),
          (math.exp(-4), 5.0))
FREE_FREE = 7.66e-58  # m^-1: absorption = this x n_e sum(n Z^2) / (sqrt(T) E^3) x (1 - exp(-E/T)); T, E in keV
EMITTING_KEV = 0.2  # cooler zones emit too little to matter, and absorb what they emit


def radiation(t_kev: np.ndarray, n_e: np.ndarray, z2: np.ndarray, dr: np.ndarray, vol: np.ndarray):
    """Bremsstrahlung emitted and absorbed by each zone (W), and the power that escapes.

    Photons from each hot zone travel along straight radial paths, half outward and half inward through the
    centre, and are absorbed by free-free absorption (Kramers, Gaunt factor 1) in the zones they cross.
    """
    nz = len(t_kev)
    t = np.maximum(t_kev, 0.01)
    emitted = BREMS * n_e * z2 * np.sqrt(t) * vol * (t_kev > EMITTING_KEV)
    absorbed = np.zeros(nz)
    idx = np.nonzero(emitted > 0)[0]
    if len(idx) == 0:
        return emitted, absorbed, 0.0
    strength = FREE_FREE * n_e * z2 / np.sqrt(t) * dr  # optical depth x E^3 / (1 - exp(-E/T))
    rows = np.arange(len(idx))
    later = np.arange(nz)[None, :] > idx[:, None]
    before = np.arange(nz)[None, :] < idx[:, None]
    same = ~later & ~before
    escaped = 0.0
    for share, x in GROUPS:
        energy = x * t[idx]  # keV
        k = strength[None, :] / energy[:, None] ** 3 * -np.expm1(-energy[:, None] / t[None, :])
        k = np.minimum(k, 1e3)
        c = np.cumsum(k, axis=1)
        cprev = c - k
        kii, cii, cprev_ii = k[rows, idx], c[rows, idx], cprev[rows, idx]
        # outward half
        t_in = np.where(later, np.maximum(cprev - cii[:, None] + 0.5 * kii[:, None], 0.0), 0.0)
        t_out = np.maximum(c - cii[:, None] + 0.5 * kii[:, None], 0.0)  # (only used where j >= i)
        out = np.where(later | same, np.exp(-t_in) - np.exp(-t_out), 0.0)
        esc_out = np.exp(-(c[:, -1] - cii + 0.5 * kii))
        # inward half: in through the zones inside, then out through all of them
        t_enter = np.maximum(0.5 * kii[:, None] + cprev_ii[:, None] - c, 0.0)  # (only used where j < i)
        inward = np.where(before, np.exp(-t_enter) - np.exp(-(t_enter + k)), 0.0)
        inward += np.where(same, -np.expm1(-0.5 * kii)[:, None], 0.0)
        t_c = (0.5 * kii + cprev_ii)[:, None]
        inward += np.exp(-(t_c + cprev)) - np.exp(-(t_c + c))
        esc_in = np.exp(-(t_c[:, 0] + c[:, -1]))
        power = emitted[idx] * share
        absorbed += power @ (0.5 * (out + inward))
        escaped += float(power @ (0.5 * (esc_out + esc_in)))
    return emitted, absorbed, escaped


# --- set-up ----------------------------------------------------------------------------------------------------------

def _graded(total: float, n: int, ratio: float) -> np.ndarray:
    """n masses summing to total, each `ratio` times the one before."""
    w = ratio ** np.arange(n)
    return total * w / w.sum()


def build(c: Capsule) -> State:
    nv, nf, na = c.zones
    rho_f = c.flight_density * 1e3
    m_fuel, m_abl = c.fuel_mg * 1e-6, c.ablator_mg * 1e-6
    r_in = c.inner_radius_um * 1e-6
    m_vap = c.vapour_ug * 1e-9
    masses = np.concatenate([_graded(m_vap, nv, 0.93)[::-1] if nv else [], _graded(m_fuel, nf, 1.035),
                             _graded(m_abl, na, 1.1) if na else []])
    kind = np.array(["vapour"] * nv + ["fuel"] * nf + ["ablator"] * na)
    rho_vap = m_vap / (4 / 3 * math.pi * r_in**3)
    dens = np.where(kind == "vapour", rho_vap, np.where(kind == "fuel", rho_f, 0.8 * rho_f))
    vol = masses / dens
    r = np.concatenate([[0.0], np.cbrt(np.cumsum(vol) * 3 / (4 * math.pi))])
    r *= r_in / r[nv] if nv else 1.0  # exact inner radius
    # carbon mixed into the vapour and the inner mix_depth of the fuel, as a uniform atom share
    fuel_ions = np.where(kind == "ablator", 0.0, masses / (A_DT * AMU))
    carbon = np.where(kind == "ablator", masses / (A_C * AMU), 0.0)
    if c.mix_ng > 0:
        inner = np.cumsum(np.where(kind == "fuel", masses, 0.0)) <= c.mix_depth * m_fuel + 1e-30
        region = (kind == "vapour") | ((kind == "fuel") & inner)
        share = masses * region / (masses * region).sum()
        n_c = c.mix_ng * 1e-12 * share / (A_C * AMU)  # replaces fuel of the same mass
        carbon = carbon + n_c
        fuel_ions = fuel_ions - n_c * A_C / A_DT
    vol = 4 / 3 * math.pi * (r[1:] ** 3 - r[:-1] ** 3)
    n_i = (fuel_ions + carbon) / vol
    n_e = (fuel_ions + 6 * carbon) / vol
    p_cold = c.adiabat * 0.4 * n_e * fermi_energy(n_e)
    t = np.where(kind == "vapour", 0.03 * KEV, temperature(p_cold, n_i, n_e))
    e = 1.5 * pressure_of(t, n_i, n_e) * vol / masses
    v = c.velocity_kms * 1e3 * math.sqrt(1 - c.rke)
    u = np.where(np.arange(len(r)) > nv, -v, 0.0)
    if nv:
        u[nv] = -v
    return State(r=r, u=u, m=masses, e=e, fuel=fuel_ions, ash=np.zeros_like(masses), carbon=carbon)


# --- one run ---------------------------------------------------------------------------------------------------------

def run(c: Capsule, t_max: float = 3e-9, cfl: float = 0.25, snapshot_every: int = 0) -> dict:
    """Implode, burn and disassemble the capsule; yield and burn-averaged diagnostics."""
    s = build(c)
    r, u, m, e = s.r, s.u, s.m, s.e
    nz = len(m)
    snapshots = []
    lost = np.zeros(nz)
    t = 0.0
    dt = 1e-13
    yield_j = 0.0
    kinetic0 = float(0.5 * np.sum(0.5 * (m[:-1] + m[1:]) * u[1:-1] ** 2) + 0.25 * m[-1] * u[-1] ** 2)
    # the ablation pressure still on the outside: the outer zone's starting pressure, held until stagnation
    p_out = c.outer_pressure * 2 / 3 * m[-1] / (4 / 3 * math.pi * (r[-1] ** 3 - r[-2] ** 3)) * e[-1]
    drive_work = 0.0
    energy0 = kinetic0 + float(np.sum(m * e))
    alpha_total = removed = 0.0
    times, powers, rho_rs, p_centre = [], [], [], []
    w_sum = w_t = w_p = w_rhor = 0.0
    peak_rhor = 0.0
    radiated = 0.0
    steps = 0
    while t < t_max:
        steps += 1
        vol = 4 / 3 * math.pi * (r[1:] ** 3 - r[:-1] ** 3)
        rho = m / vol
        ions = s.fuel + s.ash + s.carbon
        electrons = s.fuel + 2 * s.ash + 6 * s.carbon
        n_i, n_e = ions / vol, electrons / vol
        p = 2 / 3 * rho * e
        temp = temperature(p, n_i, n_e)
        t_kev = temp / KEV
        cs = np.sqrt(5 / 3 * p / rho)
        du = u[1:] - u[:-1]
        q = np.where(du < 0, rho * (2.0 * du * du + 0.1 * cs * np.abs(du)), 0.0)
        dr = r[1:] - r[:-1]

        # fusion, alphas, bremsstrahlung (rates)
        sv = bosch_hale_dt(t_kev)
        rate = (s.fuel / 2) ** 2 * sv / vol  # reactions per second per zone
        power = rate * E_FUSION_J
        z2 = (s.fuel + 4 * s.ash + 36 * s.carbon) / vol
        emitted, absorbed, escaping = radiation(t_kev, n_e, z2, dr, vol)
        net_rad = emitted - absorbed  # W lost by each zone

        # time step: sound and flow, and no zone's energy changing by more than 10% from sources
        dt_hydro = cfl * float(np.min(dr / (cs + np.abs(du) + 1e-3)))
        heat_rate = (rate * E_ALPHA_J + np.abs(net_rad)) / (m * e)
        dt_src = 0.1 / max(float(heat_rate.max()), 1e-30)
        dt = min(dt_hydro, dt_src, 1.3 * dt, t_max - t + 1e-18)

        # hydrodynamics: symplectic Euler on the staggered grid
        ptot = p + q
        force = np.zeros_like(r)
        area = 4 * math.pi * r * r
        force[1:-1] = -area[1:-1] * (ptot[1:] - ptot[:-1])
        pushing = p_out if u[-1] < 0 else 0.0  # released once the outside stops moving in
        force[-1] = area[-1] * (ptot[-1] - pushing)
        m_node = np.zeros_like(r)
        m_node[1:-1] = 0.5 * (m[:-1] + m[1:])
        m_node[-1] = 0.5 * m[-1]
        u[1:] += dt * force[1:] / m_node[1:]
        u[0] = 0.0
        r_new = r + dt * u
        r_new[0] = 0.0
        drive_work += pushing * 4 / 3 * math.pi * (r[-1] ** 3 - r_new[-1] ** 3)
        if np.any(np.diff(r_new) <= 0):
            raise FloatingPointError(f"tangled mesh at t = {t:.3e} s")
        vol_new = 4 / 3 * math.pi * (r_new[1:] ** 3 - r_new[:-1] ** 3)
        e = e - ptot * (vol_new - vol) / m
        r = r_new
        vol = vol_new
        rho = m / vol

        # burn
        burnt = np.minimum(rate * dt, 0.25 * s.fuel)
        s.fuel -= 2 * burnt
        s.ash += burnt
        fusion_now = float(burnt.sum()) * E_FUSION_J
        yield_j += fusion_now

        # alphas: deposit along radial paths, measured in ranges
        if c.alpha_heating and fusion_now > 0:
            e_alpha = burnt * E_ALPHA_J
            electrons = s.fuel + 2 * s.ash + 6 * s.carbon
            n_e = electrons / vol
            t_eff = np.sqrt(t_kev**2 + (0.6 * fermi_energy(n_e) / KEV) ** 2)
            tau_se = 6.27e8 * A_HE * (t_eff * 1e3) ** 1.5 / (4 * n_e * 1e-6 * coulomb_log(n_e, t_eff))
            s_za = (s.fuel / A_DT + s.ash * 4 / A_HE + s.carbon * 36 / A_C) / electrons
            e_c = 14.8 * A_HE * t_eff * s_za ** (2 / 3)
            rng = ALPHA_V0 * tau_se * alpha_range_factor(np.sqrt(e_c / E_ALPHA_KEV))  # m
            depth = np.concatenate([[0.0], np.cumsum((r[1:] - r[:-1]) / rng)])  # in ranges
            mid = 0.5 * (depth[:-1] + depth[1:])
            lo, hi = depth[:-1][None, :], depth[1:][None, :]

            def overlap(a, b, lo=lo, hi=hi):
                return np.clip(np.minimum(b[:, None], hi) - np.maximum(a[:, None], lo), 0.0, None)

            out = overlap(mid, mid + 1.0)
            inward = overlap(np.maximum(mid - 1.0, 0.0), mid) + overlap(np.zeros_like(mid), np.maximum(1.0 - mid, 0.0))
            share = 0.5 * (out + inward)
            dep = e_alpha @ share
            alpha_total += float(dep.sum())
            e = e + dep / m
        else:
            dep = np.zeros_like(m)

        # bremsstrahlung: emitted, partly re-absorbed in the shell, the rest escapes
        loss = np.clip(net_rad * dt, -10 * e * m, 0.5 * e * m)
        e = e - loss / m
        radiated += escaping * dt
        lost += loss
        removed += float(loss.sum())

        # electron heat conduction, implicit in T
        ions = s.fuel + s.ash + s.carbon
        electrons = s.fuel + 2 * s.ash + 6 * s.carbon
        n_i, n_e = ions / vol, electrons / vol
        p = 2 / 3 * rho * e
        temp = temperature(p, n_i, n_e)
        t_kev = np.maximum(temp / KEV, 1e-4)
        z_eff = (s.fuel + 4 * s.ash + 36 * s.carbon) / electrons
        tau_e = 3.44e5 * (t_kev * 1e3) ** 1.5 / (n_e * 1e-6 * coulomb_log(n_e, t_kev) * z_eff)
        kappa = 3.16 * n_e * temp * tau_e / M_E  # W / (m J)
        centres = 0.5 * (r[1:] + r[:-1])
        gap = centres[1:] - centres[:-1]
        k_face = 0.5 * (kappa[1:] + kappa[:-1])
        grad = np.abs(temp[1:] - temp[:-1]) / gap
        n_face = 0.5 * (n_e[1:] + n_e[:-1])
        t_face = 0.5 * (temp[1:] + temp[:-1])
        q_max = 0.1 * n_face * t_face * np.sqrt(t_face / M_E)
        k_face = k_face / (1 + k_face * grad / np.maximum(q_max, 1e-300))
        g = k_face * 4 * math.pi * r[1:-1] ** 2 / gap * dt  # J per J of temperature difference
        dp_dt = n_i + n_e * temp / np.sqrt(temp**2 + (0.4 * fermi_energy(n_e)) ** 2)
        cap = 1.5 * dp_dt * vol  # J per J of temperature
        ab = np.zeros((3, nz))
        ab[1] = cap
        ab[1, :-1] += g
        ab[1, 1:] += g
        ab[0, 1:] = -g
        ab[2, :-1] = -g
        t_new = np.maximum(solve_banded((1, 1), ab, cap * temp), 0.0)
        e = e + cap * (t_new - temp) / m
        e = np.maximum(e, 1e-6)

        t += dt
        # diagnostics
        rhor = float(np.sum(rho * (r[1:] - r[:-1]))) / 10  # g/cm2
        fuel_power = fusion_now / dt
        peak_rhor = max(peak_rhor, rhor)
        if fusion_now > 0:
            w_sum += fusion_now
            w_t += fusion_now * float(np.sum(power * t_kev) / max(power.sum(), 1e-300))
            w_p += fusion_now * float(np.sum(power * (2 / 3 * rho * e)) / max(power.sum(), 1e-300))
            w_rhor += fusion_now * rhor
        if snapshot_every and steps % snapshot_every == 0:
            snapshots.append({"t": t, "r": r.copy(), "u": u.copy(), "rho": rho.copy(), "T_keV": t_new / KEV,
                              "p": 2 / 3 * rho * e, "dep": dep.copy()})
        if steps % 5 == 0:
            times.append(t)
            powers.append(fuel_power)
            rho_rs.append(rhor)
            p_centre.append(float(2 / 3 * rho[0] * e[0]))
        if len(powers) > 50 and max(powers) > 0 and powers[-1] < 1e-3 * max(powers) and u[-1] > 0 and u[1] > 0:
            break
    times, powers = np.array(times), np.array(powers)
    kinetic_end = float(0.5 * np.sum(0.5 * (m[:-1] + m[1:]) * u[1:-1] ** 2) + 0.25 * m[-1] * u[-1] ** 2)
    supplied = energy0 + drive_work + alpha_total
    energy_error = (kinetic_end + float(np.sum(m * e)) + removed - supplied) / supplied
    k = int(np.argmax(powers)) if len(powers) else 0
    half = powers >= 0.5 * powers[k] if len(powers) else np.array([])
    fwhm = float(times[half][-1] - times[half][0]) if half.any() else 0.0
    fuel_mass = c.fuel_mg * 1e-6
    pairs = fuel_mass / (A_DT * AMU) / 2
    return {"yield_J": yield_j, "burn_up": yield_j / (pairs * E_FUSION_J),
            "bang_time_ps": float(times[k]) * 1e12 if len(times) else 0.0, "burn_width_ps": fwhm * 1e12,
            "burn_T_keV": w_t / w_sum if w_sum else 0.0, "burn_pressure_gbar": w_p / w_sum / 1e14 if w_sum else 0.0,
            "burn_rhoR": w_rhor / w_sum if w_sum else 0.0, "peak_rhoR": peak_rhor,
            "peak_centre_pressure_gbar": max(p_centre) / 1e14 if p_centre else 0.0,
            "kinetic_J": kinetic0, "drive_work_J": float(drive_work), "energy_error": energy_error, "radiated_J": radiated, "steps": steps,
            "radiated_by_zone": lost, "snapshots": snapshots}
