"""A dense Z-pinch in time: a current-driven column of fuel, from cold to burning or broken.

The simulator's fibre_pinch route (walls.fibre_pinch) is static. It assumes a
Bennett pinch already at its temperature, held for a chosen number of Alfvén
times. This model follows the pinch from the start:
* a generator (a capacitor bank behind the machine's inductance, sized to give
  the stated peak current into a fixed load) drives the current through the
  column's own changing inductance and resistance, so a load whose impedance
  climbs draws less;
* its field squeezes the column, which heats by compression and by its own
  resistance (ohmic heating);
* it loses heat by X-rays (with the column's own opacity) and by conduction
  along the field lines to the electrodes;
* it gains heat from the alphas its self-field holds in.

Where the X-rays beat ohmic heating, above the Pease-Braginskii current, the
column collapses radiatively until it heats by compression, its X-rays are
trapped, or degeneracy stops it.

Instabilities: the sausage (m = 0) mode grows at about v_A / a. The pinch
breaks once the growth adds up to `growth_limit` e-folds (about that many
Alfvén times), counted from when it is hot (> 50 eV). Two ways to slow it:
* `rise_stabilised`: no growth while the current is still rising. NRL's
  fibres (Sethian et al. 1987) looked stable during the rise and broke once it
  stopped; Lindemuth's 2-D simulations put that down to a cold core under an
  unstable corona. Untested at high current.
* `shear`: a sheared axial flow, as in the sheared-flow-stabilised Z-pinch
  (Shumlak et al.), slows the growth and stops it at 0.1 k v_A. Shown for
  ~10^3 Alfvén times at <= 0.6 MA and <= 5e24 per m3, in deuterium; untested at
  fibre-pinch density.

Single zone (uniform column), one temperature, ideal gas plus degenerate
electrons. Not included: the cold fibre's ionisation and the low-density
corona that carried the current in 1980s experiments, beam-target neutrons
from instabilities, and impurity radiation from the electrodes.
"""

from __future__ import annotations

import itertools
import math
from dataclasses import dataclass
from types import SimpleNamespace

import numpy as np
from scipy.integrate import solve_ivp

from . import physics as P
from .burn1d import BREMS, FREE_FREE, bosch_hale_dt, fermi_energy
from .confinement import (
    ALPHA_V0,
    E_ALPHA_KEV,
    alpha_range_factor,
    kappa_perp_electron,
    magnetic_thickness,
    sphere_deposition,
)

KEV = P.J_PER_KEV
MU0 = 4e-7 * math.pi
M_E = P.M_E
A_HE = 4.0026
E_DT_J = 17.59e3 * KEV
E_DD_J = 0.5 * (3.27e3 + 4.03e3) * KEV  # the two D-D branches, equally likely
SOLID_DT_IONS = 0.25e3 / (2.5 * P.AMU)  # per m3 (solid D-T, 0.25 g/cm3)
SOLID_D2_IONS = 0.2e3 / (2.014 * P.AMU)  # solid deuterium, 0.2 g/cm3
SIGMA_SB = 5.670374419e-8
K_B = 1.380649e-23


def dd_reactivity(t_kev: float) -> float:
    """Maxwellian D-D <sigma v> (m^3/s), both branches, Bosch & Hale (1992)."""
    t = min(max(t_kev, 0.2), 100.0)
    total = 0.0
    for bg, mrc2, c in ((31.3970, 937814.0, (5.65718e-12, 3.41267e-3, 1.99167e-3, 0.0, 1.05060e-5, 0.0, 0.0)),
                        (31.3970, 937814.0, (5.43360e-12, 5.85778e-3, 7.68222e-3, 0.0, -2.96400e-6, 0.0, 0.0))):
        c1, c2, c3, c4, c5, c6, c7 = c
        theta = t / (1 - t * (c2 + t * (c4 + t * c6)) / (1 + t * (c3 + t * (c5 + t * c7))))
        xi = (bg**2 / (4 * theta)) ** (1 / 3)
        total += c1 * theta * math.sqrt(xi / (mrc2 * t**3)) * math.exp(-3 * xi) * 1e-6
    return total if t_kev >= 0.2 else 0.0


@dataclass
class Pinch:
    """A Z-pinch load and its driver."""

    fuel: str = "DT"  # "DT" (50/50) or "DD"; `tritium` overrides the mix
    tritium: float | None = None  # share of the hydrogen that is tritium (0.5 for "DT", 0 for "DD")
    radius_um: float = 240.0  # starting radius (a solid fibre, or the gas column)
    ions_per_m: float = 0.0  # line density; 0 = a solid fibre of radius_um
    length_mm: float = 10.0
    current_ma: float = 26.0  # peak current
    rise_ns: float = 100.0  # 0 to peak into a fixed load
    hold_ns: float = 50.0  # run this long after the nominal peak
    driver_nh: float = 15.0  # the machine's own inductance, seen from the load
    start_kev: float = 0.01
    growth_limit: float = 1.0  # e-folds of the m = 0 mode before the pinch breaks (~1 Alfvén time)
    rise_stabilised: bool = False  # no growth while the current rises
    shear: float = 0.0  # axial-flow shear in units of k v_A; 0.1 stabilises fully
    wall_mm: float = 10.0  # return conductor
    alpha_heating: bool = True
    # a liner (MagLIF-like): a metal tube carrying the current, imploding onto preheated, magnetised fuel gas
    liner_um: float = 0.0  # outer radius; 0 = no liner
    liner_aspect: float = 6.0  # outer radius / wall thickness
    liner_density: float = 1850.0  # kg/m3 (beryllium)
    fill_mg_cc: float = 1.0  # fuel gas density inside the liner
    bz_tesla: float = 0.0  # axial field frozen into the fuel at the start
    ice_um: float = 0.0  # frozen fuel on the liner's inner wall (high-gain MagLIF); the gas fills the rest
    ice_adiabat: float = 2.0  # the ice's pressure / its degenerate pressure once compressed
    end_losses: bool = True  # heat conducted along the field out of the ends (1-D design codes leave this out)
    # The axial field leaks out of hot fuel (Nernst advection, diffusion into the cold liner): it grows as
    # (convergence)^this, not ^2. 1.57 matches MagLIF's measured B R of 0.25-0.34 T m at convergence ~30-45 from
    # 10-16 T (Lewis et al. 2023).
    flux_exponent: float = 1.57
    # The liner breaks up (magneto-Rayleigh-Taylor) once the fuel has converged this many times; 0 = never. MagLIF
    # has reached ~30-40; the model alone would go much further.
    max_convergence: float = 0.0

    @property
    def lined(self) -> bool:
        return self.liner_um > 0

    @property
    def liner_inner_um(self) -> float:
        return self.liner_um * (1 - 1 / self.liner_aspect)

    @property
    def fuel_radius_um(self) -> float:
        """The hot fuel's starting radius: inside the ice, inside the liner, or the column itself."""
        return self.liner_inner_um - self.ice_um if self.lined else self.radius_um

    @property
    def ice_ions(self) -> float:
        """Ions per metre in the ice layer."""
        if not (self.lined and self.ice_um > 0):
            return 0.0
        r_in, r_out = self.fuel_radius_um * 1e-6, self.liner_inner_um * 1e-6
        return SOLID_DT_IONS * math.pi * (r_out**2 - r_in**2)

    @property
    def outer_radius_um(self) -> float:
        return self.liner_um if self.lined else self.radius_um

    @property
    def tritium_share(self) -> float:
        return self.tritium if self.tritium is not None else (0.5 if self.fuel == "DT" else 0.0)

    @property
    def line_density(self) -> float:
        if self.lined:
            return self.fill_mg_cc * math.pi * (self.fuel_radius_um * 1e-6) ** 2 / self.ion_mass
        if self.ions_per_m > 0:
            return self.ions_per_m
        solid = SOLID_DT_IONS if self.tritium_share > 0.25 else SOLID_D2_IONS
        return solid * math.pi * (self.radius_um * 1e-6) ** 2

    @property
    def ion_mass(self) -> float:
        return (2.014 + 1.002 * self.tritium_share) * P.AMU


def circuit(p: Pinch) -> dict:
    """The generator as a capacitor bank behind the machine's inductance (driver_nh).

    Sized so that into the pinch at its starting radius (treated as fixed) the current would reach current_ma at
    rise_ns, as a quarter sine. A load whose impedance grows (a collapsing, resistive or swelling pinch) draws less.
    """
    length = p.length_mm * 1e-3
    l_load = MU0 * length / (2 * math.pi) * (math.log(p.wall_mm / (p.outer_radius_um * 1e-3)) + 0.25)
    l_sc = p.driver_nh * 1e-9 + l_load
    cap = (2 * p.rise_ns * 1e-9 / math.pi) ** 2 / l_sc
    volts = p.current_ma * 1e6 * math.sqrt(l_sc / cap)
    return {"C": cap, "V0": volts, "L0": p.driver_nh * 1e-9, "stored_J": 0.5 * cap * volts**2}


def bennett_kev(current_a: float, ions_per_m: float) -> float:
    """Bennett temperature (keV) of a hydrogen-isotope pinch: N (1 + Z) T = mu0 I^2 / 8 pi."""
    return MU0 * current_a**2 / (16 * math.pi * ions_per_m) / KEV


def pease_braginskii_ma(coulomb_log: float = 10.0) -> float:
    """Current at which bremsstrahlung equals ohmic heating in a Bennett pinch (hydrogen, optically thin).

    Ohmic I^2 eta / (pi a^2) against brems C n^2 sqrt(T) pi a^2, with Bennett's N T = mu0 I^2 / 16 pi:
    the ratio depends on I alone.
    """
    eta_coeff = 1.03e-4 * coulomb_log  # Spitzer (perpendicular): eta = this / T_eV^1.5 ohm m
    # I^2 = 256 pi^2 eta T_J^2 / (C mu0^2 sqrt(T_keV)), T-independent: evaluate at 1 keV
    t_ev, t_j = 1e3, KEV
    eta = eta_coeff / t_ev**1.5
    return math.sqrt(256 * math.pi**2 * eta * t_j**2 / (BREMS * MU0**2 * 1.0)) / 1e6


def _rates(p: Pinch, t: float, y, drive: dict) -> dict:
    a, adot, heat, deut, ash, i_now, volts, trit = y[0], y[1], y[2], y[3], y[4], y[10], y[11], y[12]
    deut, trit = max(deut, 0.0), max(trit, 0.0)
    fuel = deut + trit
    a = max(a, 1e-9)
    area = math.pi * a * a
    ions = fuel + ash
    electrons = fuel + 2 * ash
    n_e = electrons / area
    particles = ions + electrons
    t_j = max(2 / 3 * heat / particles, 1e-3 * KEV)
    t_kev = t_j / KEV
    b = MU0 * i_now / (2 * math.pi * a)
    # degenerate electrons resist compression beyond solid density (below it, the solid's cohesion balances them)
    n_solid = SOLID_DT_IONS
    p_deg = max(0.0, 0.4 * n_e * fermi_energy(n_e) - 0.4 * n_solid * fermi_energy(n_solid))
    mass = (2.014 * deut + 3.016 * trit + 3.5 * ash) * P.AMU
    iced = p.lined and p.ice_um > 0
    accel_liner = 0.0
    if p.lined:
        # An incompressible thick liner, outer radius r_out, inner radius r_in (the fuel's edge, or the ice's).
        # Lagrangian of the shell: T = pi rho r_in^2 rdot^2 ln(r_out / r_in); forces: inside pressure, field outside.
        a_fuel0 = p.fuel_radius_um * 1e-6
        bz = p.bz_tesla * (a_fuel0 / a) ** p.flux_exponent  # the axial field, partly leaked out
        p_fuel = (particles * t_j + p_deg * area) / area + bz * bz / (2 * MU0)
        r_in, rin_dot = (y[15], y[16]) if iced else (a, adot)
        shell = (p.liner_um * 1e-6) ** 2 - (p.liner_inner_um * 1e-6) ** 2
        r_out = math.sqrt(r_in * r_in + shell)
        p_field = MU0 * i_now**2 / (8 * math.pi**2 * r_out**2)
        rho_l = p.liner_density
        m_eff = 2 * math.pi * rho_l * r_in * r_in * math.log(r_out / r_in)
        inertia = math.pi * rho_l * rin_dot**2 * (2 * r_in * math.log(r_out / r_in) + r_in**3 / r_out**2 - r_in)
        if iced:
            # the ice between the hot fuel (radius a) and the liner (r_in): cold, compressed on its adiabat; its mass
            # is shared between the two surfaces it moves with
            ice = max(y[17], 1e-6 * p.ice_ions)
            n_ice = ice / max(math.pi * (r_in * r_in - a * a), 1e-30)
            p_ice = p.ice_adiabat * max(0.0, 0.4 * n_ice * fermi_energy(n_ice) - 0.4 * n_solid * fermi_energy(n_solid))
            m_half = 0.5 * ice * 2.5 * P.AMU
            accel = 2 * math.pi * a * (p_fuel - p_ice) / m_half
            accel_liner = (2 * math.pi * r_in * (p_ice - p_field) - inertia) / (m_eff + m_half)
        else:
            accel = (2 * math.pi * a * (p_fuel - p_field) - inertia) / m_eff
        b = bz  # the axial field is what holds the alphas and the heat in
    else:
        # radial motion of the column (mass per metre ~ all its ions)
        r_out = a
        accel = (2 * (particles * t_j + p_deg * area) / a - MU0 * i_now**2 / (4 * math.pi * a)) / mass
    # heating and losses, per metre
    lnl = max(2.0, 24.0 - math.log(math.sqrt(n_e * 1e-6) / (t_kev * 1e3)))
    eta = 1.03e-4 * lnl / (t_kev * 1e3) ** 1.5 * (1 + ash / max(electrons, 1.0))
    ohmic = 0.0 if p.lined else i_now**2 * eta / area  # inside a liner the fuel carries no current
    z2 = (fuel + 4 * ash) / area
    brems_thin = BREMS * n_e * z2 * math.sqrt(t_kev) * area
    tau = FREE_FREE * n_e * z2 / math.sqrt(t_kev) / t_kev**3 * a  # at h nu = T
    escape = (1 - math.exp(-tau)) / tau if tau > 1e-6 else 1.0
    surface = 2 * math.pi * a * SIGMA_SB * (t_j / K_B) ** 4  # black body
    xrays = min(brems_thin * escape, surface)
    # D-T, and D-D in its two equally likely branches: T + p (which breeds tritium in place) and 3He + n
    r_dt = deut * trit * bosch_hale_dt(t_kev) / area
    r_dd = 0.5 * deut**2 * dd_reactivity(t_kev) / area
    reactions = r_dt + r_dd
    fusion = r_dt * E_DT_J + r_dd * E_DD_J
    alpha_j = (r_dt * E_ALPHA_KEV + r_dd * 0.5 * (4.03e3 + 0.82e3)) * KEV  # the charged products
    deposited = 0.0
    escaping = 0.0
    if p.alpha_heating and alpha_j > 0:
        t_e_eff = math.sqrt(t_kev**2 + (0.6 * fermi_energy(n_e) / KEV) ** 2)
        tau_se = 6.27e8 * A_HE * (t_e_eff * 1e3) ** 1.5 / (4 * n_e * 1e-6 * lnl)
        e_c = 14.8 * A_HE * t_e_eff
        rng = ALPHA_V0 * tau_se * alpha_range_factor(math.sqrt(e_c / E_ALPHA_KEV))
        deposited = alpha_j * sphere_deposition(a / rng + magnetic_thickness(a, b))
    escaping = alpha_j - deposited
    # heat conducted along the field to the electrodes (flux-limited Spitzer)
    tau_e = 3.44e5 * (t_kev * 1e3) ** 1.5 / (n_e * 1e-6 * lnl)
    kappa = 3.16 * n_e * t_j * tau_e / M_E
    half = 0.5 * p.length_mm * 1e-3
    flux = min(kappa * t_j / half, 0.1 * n_e * t_j * math.sqrt(t_j / M_E))
    ends = 2 * flux * area / (2 * half) if p.end_losses else 0.0
    wall_loss = 0.0
    if p.lined:
        # heat conducted across the axial field into the cold liner (Braginskii, magnetised), flux-limited
        x = P.E_CHARGE * b / M_E * tau_e
        kappa_perp = n_e * t_j * tau_e / M_E * kappa_perp_electron(x)
        wall_loss = 2 * math.pi * a * min(kappa_perp * t_j / (0.5 * a), 0.1 * n_e * t_j * math.sqrt(t_j / M_E))
    compression = -(2 * particles * t_j / a) * adot  # -p dA/dt
    dheat = ohmic + compression + deposited - xrays - ends - wall_loss
    ablation = 0.0
    if iced and y[17] > 0.05 * p.ice_ions:
        # heat conducted into the ice and the alphas that escape into it ablate ice into the hot fuel: the burn
        # wave. That energy comes back with the ablated fuel, which joins the hot fuel at its temperature.
        into_ice = wall_loss + escaping
        ablation = into_ice / (3 * t_j)
        dheat += into_ice
    # the circuit: bank voltage across the machine inductance and the pinch's inductance and resistance
    length = p.length_mm * 1e-3
    wall = p.wall_mm * 1e-3
    inside = r_out < wall
    if iced:
        rdot = y[15] * y[16] / r_out
    else:
        rdot = a * adot / r_out  # the current flows at the outer surface (the liner's, or the column's)
    l_p = MU0 * length / (2 * math.pi) * ((math.log(wall / r_out) if inside else 0.0) + 0.25)
    dl_p = -MU0 * length / (2 * math.pi) * rdot / r_out if inside else 0.0
    r_p = 0.0 if p.lined else eta * length / area
    di = (volts - i_now * r_p - i_now * dl_p) / (drive["L0"] + l_p)
    dv = -i_now / drive["C"]
    v_alfven = abs(i_now) * math.sqrt(MU0 / (4 * math.pi * mass))
    suppress = max(0.0, 1 - p.shear / 0.1)
    rising = p.rise_stabilised and di > 0
    growth = suppress * v_alfven / a if t_kev > 0.05 and not rising and not p.lined else 0.0
    return {"r_dt": r_dt, "r_dd": r_dd, "di": di, "dv": dv, "l_p": l_p, "r_p": r_p, "dl_p": dl_p, "a": a, "accel": accel, "dheat": dheat, "reactions": reactions, "fusion": fusion, "ohmic": ohmic,
            "xrays": xrays, "deposited": deposited, "ends": ends, "growth": growth, "T_keV": t_kev, "n_e": n_e,
            "current": i_now, "b": b, "v_alfven": v_alfven, "compression": compression, "volts": volts,
            "wall_loss": wall_loss, "r_out": r_out, "accel_liner": accel_liner, "ablation": ablation}


# Memory and time guards: the run is integrated in CHUNKS, keeping at most KEEP_PER_CHUNK solver steps of each
# (evenly spaced in step number, so busy stretches stay well sampled); a run that needs more than MAX_CALLS
# evaluations stops early as "solver gave up" (a stiff, pathological design) instead of growing without bound.
CHUNKS = 40
KEEP_PER_CHUNK = 1000
MAX_CALLS = 300_000


class _SolverBudget(Exception):
    pass


def _integrate(rhs, t_end: float, y0: list, atol: list, max_step: float, events: tuple) -> SimpleNamespace:
    calls = [0]

    def counted(t, y):
        calls[0] += 1
        if calls[0] > MAX_CALLS:
            raise _SolverBudget
        return rhs(t, y)

    ts, ys = [0.0], [np.asarray(y0, dtype=float)]
    status, hit, gave_up = 0, [[] for _ in events], False
    y_now = np.asarray(y0, dtype=float)
    edges = np.linspace(0.0, t_end, CHUNKS + 1)
    for t0, t1 in itertools.pairwise(edges):
        try:
            seg = solve_ivp(counted, (t0, t1), y_now, method="LSODA", max_step=max_step, rtol=1e-6, atol=atol,
                            events=events)
        except _SolverBudget:
            gave_up = True
            break
        n = len(seg.t)
        if n > 1:
            keep = np.unique(np.linspace(1, n - 1, min(n - 1, KEEP_PER_CHUNK)).astype(int))
            ts.extend(seg.t[keep])
            ys.extend(seg.y[:, keep].T)
        y_now = seg.y[:, -1]
        if seg.status == 1:
            status, hit = 1, [list(e) for e in seg.t_events]
            break
        if seg.status < 0:
            gave_up = True
            break
    return SimpleNamespace(t=np.array(ts), y=np.array(ys).T, status=status, t_events=hit, gave_up=gave_up)


def run(p: Pinch, t_end_ns: float | None = None, max_step_ns: float = 0.05) -> dict:
    """Drive the pinch; returns its yield per shot and the history."""
    n0 = p.line_density
    a0 = p.fuel_radius_um * 1e-6
    heat0 = 1.5 * 2 * n0 * p.start_kev * KEV
    wall = p.wall_mm * 1e-3
    length = p.length_mm * 1e-3
    t_end = (t_end_ns or p.rise_ns + p.hold_ns) * 1e-9
    drive = circuit(p)

    def rhs(t, y):
        r = _rates(p, t, y, drive)
        # energy into the load: ohmic heating plus the work of its inductance changing (per metre)
        work = r["ohmic"] + 0.5 * y[10] ** 2 * r["dl_p"] / length
        r_dt, r_dd = r["r_dt"], r["r_dd"]
        abl = r["ablation"]  # ice ablated into the hot fuel (ions per metre per second), with the ice's mix
        return [y[1], r["accel"], r["dheat"], -r_dt - 2 * r_dd + abl * (1 - f_t), r_dt + r_dd, r["growth"],
                r["fusion"], work, r["xrays"], r["deposited"], r["di"], r["dv"], -r_dt + 0.5 * r_dd + abl * f_t,
                r_dt, 0.5 * r_dd, y[16], r["accel_liner"], -abl]

    def broken(t, y):
        return y[5] - p.growth_limit
    broken.terminal = True

    def collapsed(t, y):
        return y[0] - 1e-7
    collapsed.terminal = True

    def blown_apart(t, y):
        return y[0] - wall
    blown_apart.terminal = True

    def liner_broken(t, y):
        return y[0] - a0 / p.max_convergence if p.max_convergence > 0 else 1.0
    liner_broken.terminal = True

    f_t = p.tritium_share
    ice0 = p.ice_ions
    y0 = [a0, 0.0, heat0, float(n0) * (1 - f_t), 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, drive["V0"], float(n0) * f_t,
          0.0, 0.0, p.liner_inner_um * 1e-6 if p.lined else 0.0, 0.0, ice0]
    atol = [1e-10, 1e-2, heat0 * 1e-6, n0 * 1e-12, n0 * 1e-12, 1e-6, 1e-6, 1e-3, 1e-3, 1e-6, 1.0, 1.0, n0 * 1e-14,
            1.0, 1.0, 1e-10, 1e-2, max(ice0, n0) * 1e-12]
    sol = _integrate(rhs, t_end, y0, atol, max_step_ns * 1e-9, (broken, collapsed, blown_apart, liner_broken))
    y = sol.y[:, -1]
    t_stop = float(sol.t[-1])
    end = _rates(p, t_stop, y, drive)
    field = 0.5 * end["l_p"] * y[10] ** 2 / length
    history = [_rates(p, t, sol.y[:, k], drive) | {"t_ns": t * 1e9, "work": sol.y[7, k], "fuel": sol.y[3, k],
                                                    "gas": sol.y[3, k] + sol.y[12, k], "r_in": sol.y[15, k],
                                                    "ice": sol.y[17, k]}
               for k, t in enumerate(sol.t)]
    # An upper bound for the ice: once the hot fuel ignites, the compressed ice burns like an ICF fuel shell, the
    # share rhoR / (rhoR + 7 g/cm2) of the gas and ice together (the 0-D ablation above is the lower bound).
    ice_bound = 0.0
    if ice0 > 0:
        best_share, ice_then = 0.0, 0.0
        for h in history:
            # the hot fuel has ignited when its alphas out-heat its losses
            ignited = h["T_keV"] > 5 and h["deposited"] > h["xrays"] + h["ends"] + h["wall_loss"]
            if not ignited or h["r_in"] <= h["a"]:
                continue
            ion = 2.5 * P.AMU
            rho_gas = h["gas"] * ion / (math.pi * h["a"] ** 2)
            rho_ice = h["ice"] * ion / (math.pi * (h["r_in"] ** 2 - h["a"] ** 2))
            rho_r = (rho_gas * h["a"] + rho_ice * (h["r_in"] - h["a"])) / 10  # g/cm2
            share = rho_r / (rho_r + 7.0)
            if share > best_share:
                best_share, ice_then = share, h["ice"]
        ice_bound = best_share * ice_then / 2 * E_DT_J * length
    # what the generator puts into the load: the field's energy out to the return conductor plus the work done on
    # the plasma, at its largest (energy an expanding, burnt column pushes back is not counted as recovered)
    delivered = max(0.5 * h["l_p"] * h["current"] ** 2 / length + h["work"] for h in history)
    drawn = max(0.5 * drive["C"] * (drive["V0"] ** 2 - h["volts"] ** 2) for h in history)
    # a trace of tritium (1%) in a deuterium pinch: D-T neutrons from the same history, an ion-temperature gauge
    trace = 0.0
    if p.tritium_share == 0:
        for h0, h1 in itertools.pairwise(history):
            dt = (h1["t_ns"] - h0["t_ns"]) * 1e-9
            rate = sum(0.99 * 0.01 * h["fuel"] ** 2 * bosch_hale_dt(h["T_keV"]) / (math.pi * h["a"] ** 2)
                       for h in (h0, h1)) / 2
            trace += rate * dt * length
    status = "solver gave up" if sol.gave_up else "end of pulse"
    if sol.status == 1:
        status = ("broke up" if len(sol.t_events[0]) else "collapsed" if len(sol.t_events[1]) else
                  "liner broke up" if len(sol.t_events[3]) else "blown apart")
    hot = [h for h in history if h["T_keV"] > 0.05]
    peak = max(history, key=lambda h: h["fusion"]) if history else end
    burning = [h["t_ns"] for h in history if h["fusion"] >= 0.5 * peak["fusion"]] if peak["fusion"] > 0 else []
    burn_width = (burning[-1] - burning[0]) if burning else 0.0
    # the Alfvén times the pinch had held when 95% of its yield was out: what it needed, not what it was allowed
    burnt_by = next((k for k in range(len(sol.t)) if sol.y[6, k] >= 0.95 * y[6]), len(sol.t) - 1)
    alfven_to_burn = float(sol.y[5, burnt_by]) if y[6] > 0 else float(y[5])
    # what is left in the circuit at the end (bank and inductances): energy a recovering driver could take back
    circuit_end = 0.5 * drive["C"] * y[11] ** 2 + 0.5 * (drive["L0"] + end["l_p"]) * y[10] ** 2
    return {"yield_J": y[6] * length, "fusion_per_m_J": y[6], "neutrons": (y[13] + y[14]) * length,
            "neutrons_dt": y[13] * length, "neutrons_dd": y[14] * length, "circuit_end_J": circuit_end,
            "stopped": status, "trace_dt_neutrons": trace,
            "t_stop_ns": t_stop * 1e9, "min_radius_um": min(h["a"] for h in history) * 1e6,
            "peak_T_keV": max(h["T_keV"] for h in history), "T_at_peak_keV": peak["T_keV"],
            "radius_at_peak_um": peak["a"] * 1e6, "density_at_peak_g_cc": peak["n_e"] * p.ion_mass / 1e3,
            "burn_up": 1 - (y[3] + y[12] + y[17]) / (n0 + ice0), "delivered_J": delivered * length,
            "ice_ions": ice0, "ice_left": y[17] / ice0 if ice0 else 0.0, "ice_burn_bound_J": ice_bound,
            "burn_width_ns": burn_width, "bang_time_ns": peak.get("t_ns", 0.0),
            "field_J": field * length, "ohmic_and_work_J": y[7] * length, "xrays_J": y[8] * length,
            "alpha_heat_J": y[9] * length, "hot_alfven_times": y[5], "alfven_times_to_burn": alfven_to_burn,
            "line_density": n0,
            "bennett_keV_at_peak_current": bennett_kev(p.current_ma * 1e6, n0), "hot_ns": (hot[-1]["t_ns"] - hot[0][
                "t_ns"]) if hot else 0.0, "peak_current_MA": max(abs(h["current"]) for h in history) / 1e6,
            "bank_J": drive["stored_J"], "drawn_J": drawn, "history": history}
