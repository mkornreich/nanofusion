"""Holding a hot ball of fuel together: a tamper, a magnetic field and self-heating.

``energy.array_burn`` lets a ball of exploded crystals fly apart at its own
sound speed and burn at a fixed temperature. That is why it needs megajoules
(docs/energy.md, section 3). Magnetised target fusion and its relatives
(MagLIF, magnetised inertial fusion) use three known ways to do better
without compressing the fuel any further:

* A heavy tamper. Hot fuel inside a thick gold shell can only expand as fast
  as it pushes the gold, which is about a tenth of its own sound speed.
* A magnetic field. Fuel touching a cold wall loses its heat to it within
  picoseconds. A field parallel to the wall slows heat conduction across it
  by orders of magnitude, and it traps the 3.5 MeV alpha particles once the
  ball is a few of their gyro-radii across (B R of order 1 T m).
* Alpha self-heating. Trapped alphas heat the fuel further, so a ball that
  starts near 10 keV can run away to higher temperatures ("ignition").

Model (0-D, two temperatures)
-----------------------------
State: thermal energy of ions and of electrons, energy still carried by fast
alphas, particle numbers, and the ball radius. Terms:

* Fusion: n_D n_T <sigma v>(T_i), with burn-up. ``fuel="DD"`` is deuterium
  alone in its most generous form: every triton and helium-3 it makes is
  taken to burn at once (fully catalysed D-D, 21.6 MeV per D-D reaction, of
  which 13.4 MeV in charged particles that are deposited like alphas).
* Alphas: a fraction f of each alpha's 3.52 MeV is deposited in the fuel.
  Unmagnetised sphere: Krokhin & Rozanov's f(R / l_alpha). With a field the
  fuel is taken to be thicker by b^2 / sqrt(9 b^2 + 1000), b = R / r_Larmor,
  the effective thickness of Basko, Kemp & Meyer-ter-Vehn (Nucl. Fusion 40,
  59 (2000)), as used in MagLIF modelling. The range l_alpha comes from the
  same Stix slowing-down formula as ``energy.beam_in_plasma``. At solid
  density an alpha takes about a nanosecond to slow down, so the deposited
  energy passes through a reservoir that empties on the slowing-down time.
  A fraction 1 / (1 + 32 keV / T_e) goes to the ions (Fraley et al. 1974).
* Electron-ion energy exchange (Braginskii).
* Bremsstrahlung, with relativistic corrections in the form used by Rider
  (Phys. Plasmas 4, 1039 (1997)). All of it is lost: the fuel is transparent
  to its own X-rays. Cyclotron emission is ignored because it lies far below
  the plasma frequency at these densities and cannot leave.
* Heat conduction to the wall (tamped balls only): Braginskii's
  perpendicular conductivities for electrons and ions, a temperature fall
  to the wall over R / 4, and free-streaming limits of 0.1 n T v_th for
  electrons and 0.5 n T v_th for ions.
  ``bohm=True`` replaces the electron conductivity by Bohm diffusion,
  T / (16 e B), when that is larger: the usual pessimistic bound.
* Expansion.
  - ``tamper="none"``: a rarefaction eats into the ball at the sound speed
    and the untouched core keeps burning. At fixed temperature that gives
    the burn time R / (4 c_s) of ``energy.array_burn``.
  - ``tamper="gold"``: the wall moves at the piston speed u of the shock
    the fuel pressure p drives into cold gold, p = rho0 U_s u, using gold's
    measured shock fit and its higher compression at gigabar pressures. That
    holds while the tamper is thicker than the distance the shock travels
    (reported as ``tamper_thickness_needed_m``). For a shell of finite
    thickness (``tamper_mm``) the wall also moves at least as fast as the
    pressure has accelerated the whole shell, so a thin shell holds less.
  The fuel does p dV work on the wall and cools adiabatically.
* The field is frozen into the fuel (resistive decay takes microseconds at
  these temperatures): B ~ R^-2. Its energy counts as energy put in.
* Shape and field direction. By default the fuel is a sphere and the field
  lies along its wall everywhere (``field="ring"``). That is the ideal case:
  it needs a ring-shaped field around an axis, driven by megaamperes through
  the fuel. ``aspect`` makes the fuel a column instead (length / diameter),
  and ``field="axial"`` points the field along it, as a coil would. Heat and
  alphas then also escape along the field to the two end walls, at the
  unmagnetised rate, as in MagLIF.

Left out, all of which would make things worse: wall material mixing into
the fuel (``mix`` gives a rough bremsstrahlung-only estimate), the Nernst
effect carrying field towards the wall, instabilities, how the energy gets
in, and the entrance hole it needs (``hole_fraction`` gives a rough
estimate of the outflow).
"""

from __future__ import annotations

import math
from dataclasses import dataclass
from typing import NamedTuple

import numpy as np

from . import energy as EN
from . import physics as P

KEV = P.J_PER_KEV
MU0 = 4e-7 * math.pi
E_ALPHA_KEV = 3520.0
Q_DT_KEV = P.Q_VALUE_MEV["DT"] * 1e3
A_D, A_T, A_HE = P.MASS_U["D"], P.MASS_U["T"], P.MASS_U["He4"]
A_FUEL = 0.5 * (A_D + A_T)
ALPHA_V0 = math.sqrt(2 * E_ALPHA_KEV * KEV / (A_HE * P.AMU))  # m/s
ALPHA_LARMOR_TM = A_HE * P.AMU * ALPHA_V0 / (2 * P.E_CHARGE)  # gyro-radius times field: 0.27 T m
SOLID_DT_IONS = 2 * EN.FUELS["DT"][1]  # nuclei per m^3 in solid D-T, as elsewhere in this project
GOLD = (19240.0, 3056.0, 1.572)  # density (kg/m^3) and shock fit U_s = c0 + s u_p (LASL Shock Hugoniot Data, 1980)
HIGH_PRESSURE_COMPRESSION = 4.5  # shock compression of gold at gigabar pressures
IONISING_PRESSURE = 2.0e13  # Pa (200 Mbar): where the compression is halfway from the measured fit to that value
GOLD_CHARGE_SQUARED = 60.0**2  # mean ion charge squared of gold mixed into keV fuel (roughly nickel-like)
GRADIENT_FACTOR = 4.0  # wall heat flux = GRADIENT_FACTOR x conductivity x T / R
GROWTH_LIMIT = 5.0  # the burn is followed until the ball's radius has grown this many times
FLUX_LIMIT = 0.1  # electrons: heat flux at most 0.1 n T v_th, the usual limiter of laser-plasma codes
FLUX_LIMIT_ION = 0.5  # ions run almost freely onto a wall they can reach (one-sided flux is 0.8 n T v_th)


class Fuel(NamedTuple):
    ions_per_m3: float  # nuclei per m^3 in the solid
    mass_u: float  # mean mass of a fuel ion
    z2_over_a: float  # mean Z^2 / A of a fuel ion
    pair: tuple  # the reacting species
    pair_factor: float  # reactions per m^3 per s = pair_factor n^2 <sigma v>
    q_kev: float  # energy per reaction
    charged_kev: float  # of which carried by charged products
    burnt: float  # fuel ions used up per reaction
    ash_ions: float  # ions left behind per reaction
    ash_z2: float  # sum of their Z^2
    ash_mass_u: float  # their total mass
    ash_z2_over_a: float  # sum of their Z^2 / A


_A_P = P.MASS_U["p"]
FUELS = {
    "DT": Fuel(SOLID_DT_IONS, A_FUEL, 0.5 / A_D + 0.5 / A_T, ("D", "T"), 0.25, Q_DT_KEV, E_ALPHA_KEV,
               2.0, 1.0, 4.0, A_HE, 4 / A_HE),
    # Fully catalysed D-D: D(d,p)T and D(d,n)3He in equal shares, then T(d,n)4He and 3He(d,p)4He at once.
    # Per D-D reaction: (4.03 + 3.27 + 17.59 + 18.35) / 2 MeV, three deuterons used, one alpha and one proton left.
    "DD": Fuel(2 * EN.FUELS["D2"][1], A_D, 1 / A_D, ("D", "D"), 0.5, 21620.0, 13360.0,
               3.0, 2.0, 5.0, A_HE + _A_P, 4 / A_HE + 1 / _A_P),
}


@dataclass
class Ball:
    """A ball of fuel at the start of the burn."""

    radius_m: float
    fuel: str = "DT"  # "DT", or "DD" for deuterium alone (fully catalysed, a generous bound)
    t_i_kev: float = 10.0
    t_e_kev: float | None = None  # None: same as the ions
    fill: float = 1.0  # density as a fraction of the solid's
    b_tesla: float = 0.0
    tamper: str = "none"  # "none" or "gold"
    tamper_mm: float | None = None  # thickness of the gold shell; None: thicker than the shock travels
    alpha_heating: bool = True
    bohm: bool = False
    hole_fraction: float = 0.0  # share of the wall left open (laser entrance)
    mix: float = 0.0  # gold atoms mixed in per fuel ion
    wall_speed_factor: float = 1.0  # multiplies the wall speed, to test how much the answer depends on it
    aspect: float | None = None  # None: a sphere; a number: a column, length / diameter
    field: str = "ring"  # "ring": along the wall everywhere; "axial": along a column, meeting its end walls

    def __post_init__(self) -> None:
        if self.fuel not in FUELS:
            raise ValueError(f"fuel must be one of {', '.join(FUELS)}")
        if self.tamper not in ("none", "gold"):
            raise ValueError("tamper must be 'none' or 'gold'")
        if self.tamper_mm is not None and self.tamper_mm <= 0:
            raise ValueError("tamper_mm must be positive")
        if self.field not in ("ring", "axial"):
            raise ValueError("field must be 'ring' or 'axial'")
        if self.aspect is not None and (self.aspect <= 0 or self.tamper == "none"):
            raise ValueError("a column needs aspect > 0 and a tamper")
        if self.field == "axial" and self.aspect is None:
            raise ValueError("an axial field needs a column (set aspect)")
        if self.t_e_kev is None:
            self.t_e_kev = self.t_i_kev

    @property
    def volume_factor(self) -> float:
        """Volume / radius^3: a sphere, or a column of length 2 x aspect x radius."""
        return 4 / 3 * math.pi if self.aspect is None else 2 * self.aspect * math.pi

    @property
    def ions(self) -> float:
        return FUELS[self.fuel].ions_per_m3 * self.fill * self.volume_factor * self.radius_m**3

    @property
    def field_energy_J(self) -> float:
        return self.b_tesla**2 / (2 * MU0) * self.volume_factor * self.radius_m**3

    @property
    def energy_in_J(self) -> float:
        """Thermal energy of ions and electrons plus the magnetic field energy inside the fuel."""
        return 1.5 * self.ions * (self.t_i_kev + self.t_e_kev) * KEV + self.field_energy_J


# --- pieces of physics ---------------------------------------------------------------

_SV_CACHE: dict = {}


def dt_reactivity(t_kev: float, pair: tuple = ("D", "T")) -> float:
    """Maxwellian reactivity (m^3/s) of D-T (or another pair), interpolated from energy.reactivity."""
    key = (pair, EN.DT_SPIN_FACTOR, EN.DT_HIGH_ENERGY)
    if key not in _SV_CACHE:
        t = np.geomspace(0.1, 3000.0, 160)
        sv = np.array([EN.reactivity(*pair, float(x)) for x in t])
        _SV_CACHE[key] = (np.log(t), np.log(np.maximum(sv, 1e-300)))
    log_t, log_sv = _SV_CACHE[key]
    if t_kev < 0.1:
        return 0.0
    return math.exp(float(np.interp(math.log(min(t_kev, 3000.0)), log_t, log_sv)))


def coulomb_log_ei(n_e: float, t_e_kev: float) -> float:
    return max(2.0, 24.0 - math.log(math.sqrt(n_e * 1e-6) / (t_e_kev * 1e3)))


def coulomb_log_ii(n_i: float, t_i_kev: float) -> float:
    return max(2.0, 23.0 - math.log(math.sqrt(2 * n_i * 1e-6) / (t_i_kev * 1e3) ** 1.5))


def tau_e(n_e: float, t_e_kev: float) -> float:
    """Braginskii electron collision time (s) for singly charged ions (NRL Plasma Formulary)."""
    return 3.44e5 * (t_e_kev * 1e3) ** 1.5 / (n_e * 1e-6 * coulomb_log_ei(n_e, t_e_kev))


def tau_i(n_i: float, t_i_kev: float, mass_u: float = A_FUEL) -> float:
    """Braginskii ion collision time (s) for singly charged ions."""
    return 2.09e7 * (t_i_kev * 1e3) ** 1.5 * math.sqrt(mass_u) / (n_i * 1e-6 * coulomb_log_ii(n_i, t_i_kev))


def kappa_perp_electron(x: float) -> float:
    """Braginskii electron heat conductivity across a field, in units of n T tau_e / m_e; x = omega_ce tau_e."""
    return (4.664 * x * x + 11.92) / (x**4 + 14.79 * x * x + 3.7703)


def kappa_perp_ion(x: float) -> float:
    """Braginskii ion heat conductivity across a field, in units of n T tau_i / m_i; x = omega_ci tau_i."""
    return (2.0 * x * x + 2.645) / (x**4 + 2.70 * x * x + 0.677)


def alpha_range_factor(u_c: float) -> float:
    """Path length of a slowing alpha in units of v0 tau_se: integral of u^3 / (u^3 + u_c^3) from 0 to 1.

    u_c = sqrt(E_c / E_0): the speed, in units of the starting speed, below
    which ions rather than electrons do most of the slowing.
    """
    if u_c < 1e-9:
        return 1.0
    a = u_c
    log = math.log((1 + a) ** 2 / (1 - a + a * a)) / 6
    atan = (math.atan((2 - a) / (a * math.sqrt(3))) + math.pi / 6) / math.sqrt(3)
    return 1.0 - a * (log + atan)


def sphere_deposition(thickness: float) -> float:
    """Fraction of alpha energy left in a uniform sphere of radius `thickness` alpha ranges (Krokhin & Rozanov)."""
    x = thickness
    if x <= 0:
        return 0.0
    if x < 0.5:
        return 1.5 * x - 0.8 * x * x
    return 1.0 - 1 / (4 * x) + 1 / (160 * x**3)


def magnetic_thickness(radius_m: float, b_tesla: float) -> float:
    """Extra effective thickness (in alpha ranges) from the field: b^2 / sqrt(9 b^2 + 1000), b = R / r_Larmor."""
    b = radius_m * abs(b_tesla) / ALPHA_LARMOR_TM
    return b * b / math.sqrt(9 * b * b + 1000)


def shock_compression(pressure: float, tamper: tuple = GOLD) -> float:
    """Density ratio across the shock that `pressure` (Pa) drives into the cold tamper.

    Up to about 10 Mbar this is the measured fit U_s = c0 + s u_p. That fit can
    never compress gold more than s / (s - 1) = 2.75-fold, but above about
    100 Mbar the shocked gold ionises and compresses further, to 4-5 fold at
    gigabar pressures. The compression is taken to rise smoothly from the fit
    to 4.5 around 200 Mbar, which should be good to about 10% in wall speed.
    """
    rho0, c0, s = tamper
    u = (-c0 + math.sqrt(c0 * c0 + 4 * s * pressure / rho0)) / (2 * s)
    fit = (c0 + s * u) / (c0 + (s - 1) * u)
    return fit + (HIGH_PRESSURE_COMPRESSION - s / (s - 1)) * pressure / (pressure + IONISING_PRESSURE)


def piston_speed(pressure: float, tamper: tuple = GOLD) -> float:
    """Speed (m/s) of the wall of a thick tamper pushed by `pressure` (Pa): p = rho0 U_s u, U_s = u eta / (eta - 1)."""
    eta = shock_compression(pressure, tamper)
    return math.sqrt(pressure * (1 - 1 / eta) / tamper[0])


# --- the burn ------------------------------------------------------------------------

# state vector
_EI, _EE, _EA, _NDT, _NA, _NE, _R, _VR, _XI, _FUS, _BREM, _COND, _SHOCK, _WORK, _LEFT, _ALPHA = range(16)


def _derivatives(y: np.ndarray, ball: Ball, r0: float, tamper_mass: float) -> tuple[np.ndarray, dict]:
    fuel = FUELS[ball.fuel]
    e_i, e_e, e_a = max(y[_EI], 1e-300), max(y[_EE], 1e-300), max(y[_EA], 0.0)
    n_dt, n_a, n_e = max(y[_NDT], 0.0), max(y[_NA], 0.0), max(y[_NE], 1e-300)  # fuel ions, reactions' ash, electrons
    tamped = ball.tamper != "none"
    r = y[_R] if tamped else r0 * max(1.0 - y[_XI], 1e-6)
    vol = ball.volume_factor * r**3
    n_ions = n_dt + fuel.ash_ions * n_a
    t_i = max(e_i / (1.5 * n_ions * KEV), 1e-3)
    t_e = max(e_e / (1.5 * n_e * KEV), 1e-3)
    ne_d = n_e / vol
    ni_d = n_ions / vol
    b = ball.b_tesla * (r0 / r) ** 2 if tamped else 0.0
    # sum over ions of n Z^2 / A per electron (for slowing and energy exchange), and of n Z^2 (bremsstrahlung)
    s_za = (fuel.z2_over_a * n_dt + fuel.ash_z2_over_a * n_a) / n_e
    z2 = (n_dt + fuel.ash_z2 * n_a + ball.mix * n_dt * GOLD_CHARGE_SQUARED) / vol

    rate = fuel.pair_factor * n_dt**2 / vol * dt_reactivity(t_i, fuel.pair)  # reactions per second
    t_el = tau_e(ne_d, t_e)

    # alpha particles
    source = deposit = f_dep = 0.0
    if ball.alpha_heating:
        tau_se = 6.27e8 * A_HE * (t_e * 1e3) ** 1.5 / (4 * ne_d * 1e-6 * coulomb_log_ei(ne_d, t_e))
        e_c = 14.8 * A_HE * t_e * s_za ** (2 / 3)
        alpha_range = ALPHA_V0 * tau_se * alpha_range_factor(math.sqrt(e_c / E_ALPHA_KEV))
        half_length = r if ball.aspect is None else ball.aspect * r
        if ball.field == "axial":  # held in sideways by the field, but free to leave along it through the ends
            thickness = min(r / alpha_range + magnetic_thickness(r, b), half_length / alpha_range)
        else:
            short = min(r, half_length)
            thickness = short / alpha_range + magnetic_thickness(short, b)
        f_dep = sphere_deposition(thickness)
        tau_energy = tau_se / (2 * (1 + (2 * e_c / E_ALPHA_KEV) ** 1.5))
        source = rate * fuel.charged_kev * KEV * f_dep
        deposit = e_a / (max(f_dep, 1e-3) * tau_energy)
    to_ions = 1 / (1 + 32.0 / t_e)

    exchange = 3 * P.M_E / P.AMU * s_za * n_e * (t_e - t_i) * KEV / t_el  # W, to the ions
    theta = t_e / P.M_E_KEV
    relativistic = 1 + 0.7936 * theta + 1.874 * theta**2 + 3 / math.sqrt(2) * theta * ne_d / z2
    brem = 5.35e-37 * ne_d * z2 * math.sqrt(t_e) * relativistic * vol

    p_i = ni_d * t_i * KEV
    p_e = ne_d * t_e * KEV
    p_a = 2 / 3 * e_a / vol
    p_th = p_i + p_e + p_a
    rho = (n_dt * fuel.mass_u + n_a * fuel.ash_mass_u) * P.AMU / vol

    d = np.zeros_like(y)
    cond_e = cond_i = 0.0
    loss_rate = 0.0  # share of the fuel leaving per second
    enthalpy = 1.0
    dvol = 0.0
    if tamped:
        # wall areas: all of a sphere, or the side and the two ends of a column
        side = 4 * math.pi * r * r if ball.aspect is None else 4 * ball.aspect * math.pi * r * r
        ends = 0.0 if ball.aspect is None else 2 * math.pi * r * r
        area = side + ends
        m_i = fuel.mass_u * P.AMU
        t_io = tau_i(ni_d, t_i, fuel.mass_u)
        x_e = P.E_CHARGE * b / P.M_E * t_el
        x_i = P.E_CHARGE * b / m_i * t_io
        limit_e = FLUX_LIMIT * ne_d * t_e * KEV * math.sqrt(t_e * KEV / P.M_E)
        limit_i = FLUX_LIMIT_ION * ni_d * t_i * KEV * math.sqrt(t_i * KEV / m_i)

        def fluxes(distance: float, magnetised: bool) -> tuple[float, float]:
            """Electron and ion heat flux (W/m^2) to a wall `distance` away, across the field or along it."""
            q_e = GRADIENT_FACTOR * ne_d * t_e * KEV * t_el / P.M_E * t_e * KEV / distance
            q_e *= kappa_perp_electron(x_e if magnetised else 0.0)
            if ball.bohm and b > 0 and magnetised:
                q_e = max(q_e, GRADIENT_FACTOR * ne_d * t_e * KEV / (16 * P.E_CHARGE * b) * t_e * KEV / distance)
            q_i = GRADIENT_FACTOR * ni_d * t_i * KEV * t_io / m_i * t_i * KEV / distance
            q_i *= kappa_perp_ion(x_i if magnetised else 0.0)
            return min(q_e, limit_e), min(q_i, limit_i)

        q_e, q_i = fluxes(r, True)
        cond_e, cond_i = side * q_e, side * q_i
        if ends:
            q_e, q_i = fluxes(ball.aspect * r, ball.field != "axial")
            cond_e, cond_i = cond_e + ends * q_e, cond_i + ends * q_i
        p_mag = b * b / (6 * MU0)  # a tangled field pushes like a gas with gamma = 4/3
        eta = shock_compression(p_th + p_mag)
        speed = ball.wall_speed_factor * piston_speed(p_th + p_mag)
        d[_SHOCK] = speed * eta / (eta - 1)  # shock speed in the gold: how thick the tamper must be
        if tamper_mass > 0:  # a finite shell is also accelerated as a whole
            d[_VR] = area * (p_th + p_mag) / tamper_mass
            speed = max(speed, y[_VR])
        dvol = area * speed
        d[_R] = dvol * r / (3 * vol)  # the shape is kept: every wall recedes at `speed`
        if ball.hole_fraction > 0:
            loss_rate = ball.hole_fraction * 0.5625 * math.sqrt(5 / 3 * p_th / rho) * area / vol
            enthalpy = 5 / 3
    else:
        d[_XI] = math.sqrt((p_i + p_e) / rho) / r0
        loss_rate = 3 * d[_XI] / max(1.0 - y[_XI], 1e-6)

    out = enthalpy * loss_rate
    gone = rate * 1.5 * t_i * KEV * (fuel.burnt - fuel.ash_ions)  # heat of the ions that react, net of their ash
    d[_EI] = to_ions * deposit + exchange - cond_i - p_i * dvol - out * e_i - gone
    d[_EE] = (1 - to_ions) * deposit - exchange - brem - cond_e - p_e * dvol - out * e_e
    d[_EA] = source - deposit - p_a * dvol - out * e_a
    d[_NDT] = -fuel.burnt * rate - loss_rate * n_dt
    d[_NA] = rate - loss_rate * n_a
    d[_NE] = -loss_rate * n_e
    d[_FUS] = rate * fuel.q_kev * KEV
    d[_BREM] = brem
    d[_COND] = cond_e + cond_i
    d[_WORK] = p_th * dvol  # work done on the wall
    d[_LEFT] = out * (e_i + e_e + e_a) + gone  # heat carried off by fuel that leaves or reacts
    d[_ALPHA] = source  # alpha energy kept in the fuel
    info = {"t_i": t_i, "t_e": t_e, "radius": r, "fusion_W": d[_FUS], "f_dep": f_dep, "pressure": p_th,
            "field": b, "sound": math.sqrt((p_i + p_e) / rho)}
    return d, info


# Dormand-Prince 5(4) coefficients
_A = ((), (1 / 5,), (3 / 40, 9 / 40), (44 / 45, -56 / 15, 32 / 9),
      (19372 / 6561, -25360 / 2187, 64448 / 6561, -212 / 729),
      (9017 / 3168, -355 / 33, 46732 / 5247, 49 / 176, -5103 / 18656),
      (35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84))
_B5 = np.array([35 / 384, 0.0, 500 / 1113, 125 / 192, -2187 / 6784, 11 / 84, 0.0])
_B4 = np.array([5179 / 57600, 0.0, 7571 / 16695, 393 / 640, -92097 / 339200, 187 / 2100, 1 / 40])


def burn(ball: Ball, rtol: float = 1e-5, history: bool = False) -> dict:
    """Follow the ball until it stops burning. Returns the fusion energy, the gain and where the energy went.

    A tamped ball is followed until its radius has grown GROWTH_LIMIT-fold (its density is then under 1% of the
    start) or its fusion power has become negligible, whichever comes first; "stopped" says which. The wall model
    is not trusted beyond that growth. A ball that has ignited is still burning there, so its return is
    understated (by 20-50% in the cases tried).
    """
    r0 = ball.radius_m
    n0 = ball.ions
    y = np.zeros(16)
    y[_EI] = 1.5 * n0 * ball.t_i_kev * KEV
    y[_EE] = 1.5 * n0 * ball.t_e_kev * KEV
    y[_NDT] = n0
    y[_NE] = n0
    y[_R] = r0
    tamper_mass = 0.0
    if ball.tamper_mm is not None:
        shell = ball.tamper_mm * 1e-3
        if ball.aspect is None:
            outer = 4 / 3 * math.pi * (r0 + shell) ** 3
        else:
            outer = math.pi * (r0 + shell) ** 2 * 2 * (ball.aspect * r0 + shell)
        tamper_mass = GOLD[0] * (outer - ball.volume_factor * r0**3)
    thermal0 = y[_EI] + y[_EE]
    scale = np.array([thermal0] * 3 + [n0] * 3 + [r0, 1e6, 1.0] + [thermal0] * 3 + [r0] + [thermal0] * 3)
    d, info = _derivatives(y, ball, r0, tamper_mass)
    crossing = r0 / info["sound"]
    t, h = 0.0, 1e-5 * crossing
    stopped = "step limit"
    peak_ti, peak_te, t_double = info["t_i"], info["t_e"], None
    track = [(0.0, info["t_i"], info["t_e"], r0, 0.0)] if history else None
    steps = 0
    while steps < 200000:
        steps += 1
        k = [d]
        for row in _A[1:]:
            k.append(_derivatives(y + h * sum(c * ki for c, ki in zip(row, k)), ball, r0, tamper_mass)[0])
        k = np.array(k)
        y5 = y + h * (_B5 @ k)
        error = float(np.max(np.abs(h * ((_B5 - _B4) @ k)) / (1e-9 * scale + rtol * np.abs(y5))))
        if error <= 1.0 or h < 1e-22:
            if y5[_R] > GROWTH_LIMIT * r0:  # stop exactly where the ball has grown that far, not a step beyond
                part = (GROWTH_LIMIT * r0 - y[_R]) / (y5[_R] - y[_R])
                y5 = y + part * (y5 - y)
                h *= part
            t += h
            y = y5
            d, info = _derivatives(y, ball, r0, tamper_mass)
            peak_ti, peak_te = max(peak_ti, info["t_i"]), max(peak_te, info["t_e"])
            if t_double is None and info["radius"] >= 2 ** (1 / 3) * r0:
                t_double = t
            if history:
                track.append((t, info["t_i"], info["t_e"], info["radius"], float(y[_FUS])))
            if ball.tamper == "none" and (y[_XI] > 0.995 or y[_NDT] < 1e-4 * n0):
                stopped = "flew apart"
            elif y[_R] >= GROWTH_LIMIT * r0 * (1 - 1e-12):
                stopped = "grew to the limit"
            elif y[_NDT] < 1e-4 * n0:
                stopped = "fuel burnt"
            elif info["t_i"] < 0.2 or (t > 3 * crossing
                                       and info["fusion_W"] * t < 1e-5 * max(y[_FUS], 1e-6 * thermal0)):
                stopped = "burn died out"
            if stopped != "step limit":
                break
        h *= min(5.0, max(0.2, 0.9 * error ** -0.2)) if error > 0 else 5.0
    e_in = ball.energy_in_J
    fuel = FUELS[ball.fuel]
    burnt = fuel.burnt * y[_FUS] / (fuel.q_kev * KEV) / n0
    out = {"radius_m": r0, "fuel": ball.fuel, "t_i_keV": ball.t_i_kev, "t_e_keV": ball.t_e_kev, "fill": ball.fill,
           "b_tesla": ball.b_tesla, "tamper": ball.tamper, "energy_in_J": e_in, "field_energy_J": ball.field_energy_J,
           "fusion_J": float(y[_FUS]), "gain": float(y[_FUS]) / e_in, "burn_fraction": float(burnt),
           "bremsstrahlung_J": float(y[_BREM]), "conduction_J": float(y[_COND]), "work_on_wall_J": float(y[_WORK]),
           "carried_off_J": float(y[_LEFT]), "alpha_heating_J": float(y[_ALPHA]),
           "heat_left_J": float(y[_EI] + y[_EE] + y[_EA]), "heat_at_start_J": thermal0,
           "peak_t_i_keV": peak_ti, "peak_t_e_keV": peak_te, "ignited": peak_ti > 1.5 * max(ball.t_i_kev, ball.t_e_kev),
           "burn_time_s": t, "time_to_double_volume_s": t_double, "tamper_thickness_needed_m": float(y[_SHOCK]),
           "mass_g": n0 * fuel.mass_u * P.AMU * 1e3,
           "areal_density_g_cm2": ball.fill * fuel.ions_per_m3 * fuel.mass_u * P.AMU * r0 * 0.1,
           "b_r_tesla_m": ball.b_tesla * r0, "steps": steps, "stopped": stopped}
    if history:
        out["history"] = track
    return out


def smallest_ball(target_gain: float, temperatures=(5.0, 7.0, 10.0, 14.0, 20.0, 28.0, 40.0),
                  r_lo: float = 2e-6, r_hi: float = 0.2, rel: float = 0.03, **settings) -> dict | None:
    """The cheapest ball (least energy put in) that reaches `target_gain`, over starting temperature and radius.

    For each starting temperature the smallest radius that reaches the gain is
    found by bisection (gain rises with size). Returns the burn() result of the
    cheapest one, or None if no ball up to r_hi gets there.
    """
    best = None
    for t0 in temperatures:
        make = lambda r, t0=t0: burn(Ball(radius_m=r, t_i_kev=t0, **settings))  # noqa: E731
        hi = make(r_hi)
        if hi["gain"] < target_gain:
            continue
        lo_r, hi_r = r_lo, r_hi
        if make(lo_r)["gain"] >= target_gain:
            hi_r, hi = lo_r, make(lo_r)
        while hi_r / lo_r > 1 + rel:
            mid = math.sqrt(lo_r * hi_r)
            res = make(mid)
            if res["gain"] >= target_gain:
                hi_r, hi = mid, res
            else:
                lo_r = mid
        if best is None or hi["energy_in_J"] < best["energy_in_J"]:
            best = hi
    return best


def main(argv: list[str] | None = None) -> int:
    import argparse
    import json
    import sys

    p = argparse.ArgumentParser(prog="python -m nanofusion ball",
                                description="Burn of a hot fuel ball held by a gold tamper and a magnetic field.")
    p.add_argument("--radius", type=float, default=0.28, help="ball radius, mm")
    p.add_argument("--temperature", type=float, default=20.0, help="starting ion temperature, keV")
    p.add_argument("--electron-temperature", type=float, help="starting electron temperature, keV (default: same)")
    p.add_argument("--fuel", choices=list(FUELS), default="DT")
    p.add_argument("--fill", type=float, default=1.0, help="density as a fraction of the solid's")
    p.add_argument("--field", type=float, default=3000.0, help="magnetic field, tesla")
    p.add_argument("--tamper", choices=("none", "gold"), default="gold")
    p.add_argument("--tamper-mm", type=float, help="gold shell thickness, mm (default: thick)")
    p.add_argument("--aspect", type=float, help="make the fuel a column: length / diameter (default: a sphere)")
    p.add_argument("--axial-field", action="store_true",
                   help="field along the column, so heat also escapes to its ends (default: along the wall everywhere)")
    p.add_argument("--hole", type=float, default=0.0, help="share of the wall left open")
    p.add_argument("--mix", type=float, default=0.0, help="gold atoms mixed in per fuel ion")
    p.add_argument("--bohm", action="store_true", help="Bohm heat loss across the field")
    p.add_argument("--no-alpha-heating", action="store_true")
    p.add_argument("-o", "--output", help="write the result as JSON")
    args = p.parse_args(argv)
    if args.radius <= 0 or args.temperature <= 0 or not 0 < args.fill <= 1 or not 0 <= args.hole < 1:
        print("radius and temperature must be positive, fill in (0, 1], hole in [0, 1).", file=sys.stderr)
        return 2
    try:
        ball = Ball(radius_m=args.radius * 1e-3, t_i_kev=args.temperature, t_e_kev=args.electron_temperature,
                    fuel=args.fuel, fill=args.fill, b_tesla=args.field, tamper=args.tamper, tamper_mm=args.tamper_mm,
                    alpha_heating=not args.no_alpha_heating, bohm=args.bohm, hole_fraction=args.hole, mix=args.mix,
                    aspect=args.aspect, field="axial" if args.axial_field else "ring")
    except ValueError as exc:
        print(exc, file=sys.stderr)
        return 2
    res = burn(ball)
    if args.output:
        with open(args.output, "w") as fh:
            json.dump(res, fh, indent=1, default=float)
    e_in = res["energy_in_J"]
    print(f"{res['mass_g'] * 1e6:.3g} ug of {args.fuel}, energy in {e_in:.3g} J (field {res['field_energy_J']:.3g} J), "
          f"B R = {res['b_r_tesla_m']:.3g} T m")
    print(f"fusion {res['fusion_J']:.3g} J, gain {res['gain']:.3g}, burnt {res['burn_fraction']:.2%}, "
          f"peak ion temperature {res['peak_t_i_keV']:.3g} keV{' (ignited)' if res['ignited'] else ''}; "
          f"stopped after {res['burn_time_s'] * 1e9:.3g} ns: {res['stopped']}")
    heat = res["heat_at_start_J"] + res["alpha_heating_J"]
    print(f"of its heat (start + alphas, {heat:.3g} J): {res['conduction_J'] / heat:.0%} conducted to the wall, "
          f"{res['bremsstrahlung_J'] / heat:.0%} radiated, {res['work_on_wall_J'] / heat:.0%} work on the wall")
    if res["time_to_double_volume_s"]:
        print(f"volume doubles after {res['time_to_double_volume_s'] * 1e9:.3g} ns; "
              f"tamper must be at least {res['tamper_thickness_needed_m'] * 1e3:.3g} mm thick")
    return 0
