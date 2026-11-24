"""Improved physics model of laser-driven LiD nanocrystal fusion.

What changes compared with the original (Java) model:

* Every ion feels every other ion's Coulomb field. The original only used
  neighbours within 2 lattice steps, which removes the Coulomb explosion of the
  stripped crystal.
* The laser is a real pulse: a Gaussian envelope from a vector potential
  (zero net field area), so there is no unphysical DC push. The magnetic
  force q v x B uses B = E/c, which the original's "Lorentz" term omitted.
* Velocity-Verlet integration in double precision with an adaptive time step,
  instead of 6-digit-rounded explicit Euler.
* Fusion probability comes from measured cross sections, which include
  quantum tunnelling, instead of "pairs closer than 2 fm". It is counted
  inside the crystal during the run, and afterwards for ions that fly out and
  stop in surrounding cold LiD (beam-target).
* Checks whether the laser can really strip all electrons, plus an energy balance.

Kept from the proposal: all electrons are removed at t = 0 and the nuclei
start at rest on LiD lattice sites. Electrons are not simulated.
"""

from __future__ import annotations

import math
import time
from dataclasses import asdict, dataclass, field
from typing import Callable, Optional

import numpy as np

from ..model import Frame, SimulationResult, SimulatorProperties
from . import physics as P
from .native import HIST_EDGES_KEV, NBINS, coulomb_step, native_available


class Cancelled(Exception):
    pass


@dataclass
class ImprovedProperties:
    name: str = "Improved model"
    shape: str = "sphere"  # "sphere", "cube" or "disk"
    size: float = 2.0  # sphere and disk: radius in nm; cube: ions per edge
    lithium: str = "Li6"  # "Li6" (enriched, best case) or "Li7" (as in the original)
    intensity_w_cm2: float = 1e18
    wavelength_um: float = 3.6
    fwhm_fs: float = 6.0
    phase_deg: float = 0.0  # carrier phase at the pulse centre (0 = field crest)
    centre_fs: float = 0.0  # pulse centre relative to the stripping time t = 0
    duration_fs: float = 40.0
    frames: int = 200  # frames kept for the viewer
    softening_nm: float = 0.001
    max_step_fs: float = 0.02
    target_radius_nm: float = 15.0  # crystal size to extrapolate yields to (30 nm crystal)
    profile: str = "uniform"  # "uniform" crystal, or "expanded": pre-expanded by a prepulse
    core_density: float = 1.0  # expanded profile: central density / solid LiD density
    expansion: float = 3.0  # expanded profile: how far the edge is stretched (edge radius / original)
    geometry: str = "single"  # "single" crystal, or "pair": two crystals side by side along x
    gap_nm: float = 1.0  # pair: surface-to-surface gap
    thickness_nm: float = 1.0  # disk: thickness along z (the disk's axis)
    dot_nm: float = 0.0  # disk: if > 0, a lithium-only disk with a deuterium layer of this radius on its +z face

    def validate(self) -> None:
        problems = []
        if self.shape not in ("sphere", "cube", "disk"):
            problems.append("Shape must be 'sphere', 'cube' or 'disk'.")
        if self.lithium not in ("Li6", "Li7"):
            problems.append("Lithium isotope must be Li6 or Li7.")
        if self.shape in ("sphere", "disk") and not 0.2 <= self.size <= 20:
            problems.append("Sphere or disk radius must be between 0.2 and 20 nm.")
        if self.shape == "disk" and not 0.1 <= self.thickness_nm <= 40:
            problems.append("Disk thickness must be between 0.1 and 40 nm.")
        if not 0 <= self.dot_nm <= self.size:
            problems.append("Deuterium dot radius must be between 0 and the disk radius.")
        if self.shape == "cube" and not (2 <= self.size <= 60 and float(self.size).is_integer()):
            problems.append("Cube edge must be a whole number of ions between 2 and 60.")
        for label, value, lo, hi in (
            ("Intensity (W/cm²)", self.intensity_w_cm2, 0, 1e26),
            ("Wavelength (µm)", self.wavelength_um, 0.1, 100),
            ("Pulse FWHM (fs)", self.fwhm_fs, 0.5, 1000),
            ("Duration (fs)", self.duration_fs, 0.1, 5000),
            ("Target radius (nm)", self.target_radius_nm, 0.1, 1e4),
        ):
            if not lo <= value <= hi:
                problems.append(f"{label} must be between {lo:g} and {hi:g}.")
        if not 2 <= self.frames <= 2000:
            problems.append("Frames must be between 2 and 2000.")
        if self.profile not in ("uniform", "expanded"):
            problems.append("Profile must be 'uniform' or 'expanded'.")
        if not 1.0 <= self.expansion <= 20.0:
            problems.append("Expansion must be between 1 and 20.")
        if not 0.01 <= self.core_density <= 1.0:
            problems.append("Core density must be between 0.01 and 1 (fraction of solid LiD).")
        if self.geometry not in ("single", "pair"):
            problems.append("Geometry must be 'single' or 'pair'.")
        if not 0 <= self.gap_nm <= 1000:
            problems.append("Gap must be between 0 and 1000 nm.")
        if problems:
            raise ValueError("\n".join(problems))


def build_crystal(props: ImprovedProperties) -> tuple[np.ndarray, np.ndarray]:
    """Rock-salt LiD: positions (m, centred on the origin) and is_lithium flags."""
    a = P.SPACING
    if props.shape == "cube":
        n = int(props.size)
        idx = np.indices((n, n, n)).reshape(3, -1).T
    elif props.shape == "disk":
        r_sites = props.size * 1e-9 / a
        h_sites = props.thickness_nm * 1e-9 / (2 * a)
        m = int(math.ceil(r_sites)) + 1
        mz = int(math.ceil(h_sites)) + 1
        rng, rz = np.arange(-m, m + 1), np.arange(-mz, mz + 1)
        idx = np.stack(np.meshgrid(rng, rng, rz, indexing="ij"), axis=-1).reshape(-1, 3)
        keep = (idx[:, 0].astype(float) ** 2 + idx[:, 1].astype(float) ** 2 <= r_sites**2) & \
            (np.abs(idx[:, 2]) <= h_sites + 1e-9)
        idx = idx[keep]
        if props.dot_nm > 0:
            # Double-layer target: lithium-only disk, plus one deuterium layer on its +z face.
            idx = idx[(idx.sum(axis=1) % 2) == 1]
            top = int(idx[:, 2].max()) + 1
            d_sites = props.dot_nm * 1e-9 / a
            md = int(math.ceil(d_sites)) + 1
            rd = np.arange(-md, md + 1)
            dot = np.stack(np.meshgrid(rd, rd, [top], indexing="ij"), axis=-1).reshape(-1, 3)
            dot = dot[dot[:, 0].astype(float) ** 2 + dot[:, 1].astype(float) ** 2 <= d_sites**2]
            pos = np.concatenate([idx, dot]).astype(float) * a
            is_li = np.concatenate([np.ones(len(idx), bool), np.zeros(len(dot), bool)])
            order = np.lexsort((pos[:, 0], pos[:, 1], pos[:, 2]))
            return pos[order], is_li[order]
    else:
        r_sites = props.size * 1e-9 / a
        m = int(math.ceil(r_sites)) + 1
        rng = np.arange(-m, m + 1)
        idx = np.stack(np.meshgrid(rng, rng, rng, indexing="ij"), axis=-1).reshape(-1, 3)
        idx = idx[np.sum(idx.astype(float) ** 2, axis=1) <= r_sites**2]
    is_li = (idx.sum(axis=1) % 2) == 1
    pos = idx.astype(float) * a
    pos -= pos.mean(axis=0)
    if props.profile == "expanded":
        pos = _expand(pos, props.core_density, props.expansion)
    if props.geometry == "pair":
        half = np.max(np.linalg.norm(pos, axis=1)) + props.gap_nm * 1e-9 / 2 + a / 2
        pos = np.concatenate([pos - [half, 0, 0], pos + [half, 0, 0]])
        is_li = np.concatenate([is_li, is_li])
    order = np.lexsort((pos[:, 0], pos[:, 1], pos[:, 2]))
    return pos[order], is_li[order]


def _expand(pos: np.ndarray, core_density: float, expansion: float) -> np.ndarray:
    """Stretch the crystal radially so its density falls off towards the edge.

    Models a crystal that a weak prepulse has let expand before the main
    pulse. Each ion moves out along its own radius:
        r' = core^(-1/3) * r * (1 + (E - 1) (r/R)^2),
    which never brings ions closer together (both r'/r and dr'/dr are >= 1).
    The centre keeps ``core_density`` x solid density; at E = 3 the edge is 63
    times thinner than the centre.

    Kaplan, Dubetsky & Shkolnikov (PRL 91, 143401) showed that profiles like
    this make a Coulomb explosion form shock shells: inner ions overtake and
    collide with outer ones wherever the local density is below a third of the
    mean density inside that radius.
    """
    r = np.linalg.norm(pos, axis=1)
    r_max = r.max()
    if r_max == 0:
        return pos
    factor = core_density ** (-1 / 3) * (1 + (expansion - 1) * (r / r_max) ** 2)
    return pos * factor[:, None]


@dataclass
class ImprovedReport:
    inputs: dict
    particles: int
    deuterons: int
    radius_nm: float
    native_kernel: bool
    steps: int
    wall_seconds: float
    simulated_fs: float
    laser: dict = field(default_factory=dict)
    validity: dict = field(default_factory=dict)
    dynamics: dict = field(default_factory=dict)
    fusion: dict = field(default_factory=dict)
    energy: dict = field(default_factory=dict)
    scaled: dict = field(default_factory=dict)
    spectra: dict = field(default_factory=dict)
    size_table: list = field(default_factory=list)
    pair_exposure: dict = field(default_factory=dict)


def run_improved(
    props: ImprovedProperties,
    progress: Optional[Callable[[int, int], None]] = None,
    cancelled: Optional[Callable[[], bool]] = None,
) -> tuple[SimulationResult, ImprovedReport]:
    props.validate()
    started = time.monotonic()
    pos, is_li = build_crystal(props)
    kind = is_li.astype(np.int32)
    li_mass = P.MASS_U[props.lithium]
    mass = np.where(is_li, li_mass, P.MASS_U["D"]) * P.AMU
    q = np.where(is_li, 3.0, 1.0) * P.E_CHARGE
    q_over_m = q / mass
    vel = np.zeros_like(pos)

    e0 = P.field_from_intensity(props.intensity_w_cm2)
    pulse = P.Pulse(e0=e0, wavelength=props.wavelength_um * 1e-6, fwhm=props.fwhm_fs * 1e-15,
                    centre=props.centre_fs * 1e-15, phase=math.radians(props.phase_deg))
    bg_dli = P.gamow_constant(1, 3, P.reduced_mass_u(P.MASS_U["D"], li_mass))
    s_dli = P.S_DLI_KEV_MB[props.lithium]
    soft2 = (props.softening_nm * 1e-9) ** 2
    rc = 2 * P.SPACING  # reaction volume radius
    rc2, inv_vc = rc * rc, 1.0 / (4 / 3 * math.pi * rc**3)

    exposure_dd = np.zeros(NBINS)
    exposure_dli = np.zeros(NBINS)
    last_hist = [None, None]

    def accel(x, v, t):
        f, stats, h_dd, h_dli = coulomb_step(x, v, q, mass, kind, soft2, rc2, inv_vc, s_dli, bg_dli)
        last_hist[0], last_hist[1] = h_dd, h_dli
        ex = pulse.field(t)
        a = f / mass[:, None]
        # Laser: E along x, B = E/c along y (propagation +z): F = q (E + v x B)
        b = ex / P.C_LIGHT
        a[:, 0] += q_over_m * (ex - v[:, 2] * b)
        a[:, 2] += q_over_m * (v[:, 0] * b)
        return a, stats

    t_end = props.duration_fs * 1e-15
    period = props.wavelength_um * 1e-6 / P.C_LIGHT
    dt_max = min(props.max_step_fs * 1e-15, period / 200)
    dt_min = dt_max * 1e-4
    t = 0.0
    a, stats = accel(pos, vel, t)
    potential0 = stats[5]
    yield_dli = yield_dd = 0.0
    min_dli = stats[2]
    min_dli_t = 0.0
    max_e_dli = stats[3]
    frame_times = [t_end * k / (props.frames - 1) for k in range(props.frames)]
    next_index = 0
    frames: list[Frame] = []
    steps = 0
    work_by_laser = 0.0

    def snapshot():
        frames.append(Frame(time=_decimal_time(t), execution_time_ms=int((time.monotonic() - started) * 1000),
                            positions=pos.copy()))

    while True:
        if next_index < len(frame_times) and t >= frame_times[next_index] * (1 - 1e-9):
            snapshot()
            next_index += 1
        if t >= t_end * (1 - 1e-12):
            if next_index < len(frame_times):  # final state, if rounding skipped it
                snapshot()
            break
        next_frame = frame_times[min(next_index, len(frame_times) - 1)]
        if cancelled is not None and cancelled():
            raise Cancelled()
        amax = float(np.sqrt(np.max(np.einsum("ij,ij->i", a, a))))
        # Longer steps once the pulse is over (the remaining motion is smooth expansion).
        dt = dt_max if abs(t - pulse.centre) < 4 * pulse.tau_g + 0.5 * period else 10 * dt_max
        if stats[4] < 1e299:
            dt = min(dt, 0.02 * stats[4])  # resolve approaching neighbours
        if amax > 0:
            dt = min(dt, math.sqrt(0.002 * min(stats[6], P.SPACING) / amax))  # smooth force changes
        dt = max(dt_min, min(dt, next_frame - t, t_end - t))
        v_half = vel + 0.5 * dt * a
        pos += dt * v_half
        t += dt
        work_by_laser += float(np.sum(q * pulse.field(t - 0.5 * dt) * v_half[:, 0])) * dt
        a, stats = accel(pos, v_half, t)
        vel = v_half + 0.5 * dt * a
        yield_dli += stats[0] * dt
        yield_dd += stats[1] * dt
        exposure_dd += last_hist[0] * dt
        exposure_dli += last_hist[1] * dt
        if stats[2] < min_dli:
            min_dli, min_dli_t = stats[2], t
        max_e_dli = max(max_e_dli, stats[3])
        steps += 1
        if progress is not None:
            progress(int(1000 * t / t_end), 1000)

    report = _report(props, pos, vel, is_li, mass, q, e0, pulse, steps, started, t,
                     yield_dli, yield_dd, min_dli, min_dli_t, max_e_dli, potential0, stats[5], work_by_laser)
    report.pair_exposure = {
        "edges_keV": HIST_EDGES_KEV.tolist(),
        "D_D": exposure_dd.tolist(),
        "D_Li": exposure_dli.tolist(),
        "note": "sum over close pairs of |v_rel| dt / V_c (1/m^2), binned by relative kinetic energy; "
                "times a cross section (m^2) gives expected reactions",
    }
    report.final_state = {"velocities": vel, "is_lithium": is_li, "mass": mass}  # for tools; not serialised
    result = _as_result(props, frames, is_li, report)
    return result, report


def _decimal_time(t: float):
    from decimal import Decimal

    return Decimal(f"{t:.4e}")


def _report(props, pos, vel, is_li, mass, q, e0, pulse, steps, started, t, yield_dli, yield_dd,
            min_dli, min_dli_t, max_e_dli, pot0, pot_end, work_by_laser) -> ImprovedReport:
    ke = 0.5 * mass * np.einsum("ij,ij->i", vel, vel) / P.J_PER_KEV  # keV
    d = ~is_li
    n = len(q)
    crystals = 2 if props.geometry == "pair" else 1
    radius = _equivalent_radius(n // crystals)  # of one crystal
    wavelength = props.wavelength_um * 1e-6
    e_strip = P.stripping_field(radius)
    density_e = 16 / (2 * P.SPACING) ** 3
    n_crit = P.critical_density(wavelength)
    skin = P.C_LIGHT / math.sqrt(density_e * P.E_CHARGE**2 / (P.EPS0 * P.M_E))

    # Beam-target: every ion that leaves stops in a cold catcher (LiD or LiT).
    y_bt_li, y_bt_dd = P.thick_target_yield(ke[d], props.lithium)
    q_li, q_dd = P.Q_VALUE_MEV[props.lithium], P.Q_VALUE_MEV["DD"]
    catchers = {c: _catcher_yields(ke[d], ke[~d], c, props.lithium) for c in ("LiD", "LiT")}

    # Extrapolate to the target crystal size: Coulomb-explosion velocities scale
    # as R at fixed density (energies as R^2), ion numbers as R^3.
    scale_r = props.target_radius_nm * 1e-9 / radius
    ke_all_scaled = _scaled_energies(vel, mass, d, scale_r)
    ke_scaled = ke_all_scaled[d]
    count_scale = scale_r**3
    n_target_d = int(round(d.sum() * count_scale / crystals))

    fluence = pulse.fluence()
    area = math.pi * radius**2
    area_t = math.pi * (props.target_radius_nm * 1e-9) ** 2
    mev = 1.602176634e-13

    hist_edges = np.linspace(0, max(float(ke.max()), 1e-6), 41)
    report = ImprovedReport(
        inputs=asdict(props),
        particles=n,
        deuterons=int(d.sum()),
        radius_nm=radius * 1e9,
        native_kernel=native_available(),
        steps=steps,
        wall_seconds=time.monotonic() - started,
        simulated_fs=t * 1e15,
    )
    report.laser = {
        "peak_field_V_per_m": e0,
        "intensity_W_per_cm2": props.intensity_w_cm2,
        "period_fs": wavelength / P.C_LIGHT * 1e15,
        "cycles_in_fwhm": props.fwhm_fs / (wavelength / P.C_LIGHT * 1e15),
        "electron_a0": P.electron_a0(e0, wavelength),
        "fluence_J_per_cm2": fluence * 1e-4,
    }
    report.validity = {
        "field_needed_to_strip_V_per_m": e_strip,
        "field_needed_to_strip_target_V_per_m": P.stripping_field(props.target_radius_nm * 1e-9),
        "full_stripping_plausible": e0 >= e_strip,
        "full_stripping_plausible_at_target_size": e0 >= P.stripping_field(props.target_radius_nm * 1e-9),
        "intensity_needed_to_strip_target_W_per_cm2": P.intensity_w_cm2(
            P.stripping_field(props.target_radius_nm * 1e-9)),
        "electron_density_over_critical": density_e / n_crit,
        "skin_depth_nm": skin * 1e9,
    }
    report.dynamics = {
        "max_D_keV": float(ke[d].max()),
        "mean_D_keV": float(ke[d].mean()),
        "max_Li_keV": float(ke[~d].max()) if (~d).any() else 0.0,
        "mean_Li_keV": float(ke[~d].mean()) if (~d).any() else 0.0,
        "coulomb_explosion_theory_max_D_keV": P.coulomb_explosion_max_kev(radius, 1),
        "coulomb_explosion_theory_max_Li_keV": P.coulomb_explosion_max_kev(radius, 3),
        "closest_D_Li_approach_A": min_dli * 1e10,
        "closest_D_Li_approach_time_fs": min_dli_t * 1e15,
        "max_D_Li_pair_energy_keV": max_e_dli,
        "D_drift_keV": float(0.5 * mass[d][0] * np.sum(vel[d].mean(axis=0) ** 2) / P.J_PER_KEV),
        "Li_drift_keV": float(0.5 * mass[~d][0] * np.sum(vel[~d].mean(axis=0) ** 2) / P.J_PER_KEV)
        if (~d).any() else 0.0,
        "energy_check_ion_KE_minus_potential_drop_keV": float(ke.sum() - (pot0 - pot_end) / P.J_PER_KEV),
        "laser_work_on_ions_keV": work_by_laser / P.J_PER_KEV,
    }
    lid = catchers["LiD"]
    report.fusion = {
        "in_crystal_D_Li": yield_dli,
        "in_crystal_D_D": yield_dd,
        "beam_target_D_Li": float(y_bt_li.sum()),
        "beam_target_D_D": float(y_bt_dd.sum()),
        "beam_target_Li_D": lid["channels"].get("Li+D", 0.0),
        "catchers": catchers,
        "total": yield_dli + yield_dd + lid["total"],
        "fusion_energy_MeV": yield_dli * q_li + yield_dd * q_dd + lid["energy_MeV"],
    }
    report.energy = {
        "laser_energy_through_crystal_MeV": fluence * area / mev,
        "gain": _ratio(report.fusion["fusion_energy_MeV"], fluence * area / mev),
    }
    scaled_lid = _catcher_yields(ke_scaled, ke_all_scaled[~d], "LiD", props.lithium, count_scale / crystals)
    report.scaled = {
        "target_radius_nm": props.target_radius_nm,
        "deuterons": n_target_d,
        "max_D_keV": float(ke_scaled.max()),
        "beam_target_D_Li": scaled_lid["channels"].get("D+Li", 0.0),
        "beam_target_D_D": scaled_lid["channels"].get("D+D", 0.0),
        "beam_target_Li_D": scaled_lid["channels"].get("Li+D", 0.0),
        "total": scaled_lid["total"],
        "fusion_energy_MeV": scaled_lid["energy_MeV"],
        "laser_energy_through_crystal_MeV": fluence * area_t / mev,
    }
    report.scaled["gain"] = _ratio(report.scaled["fusion_energy_MeV"], report.scaled["laser_energy_through_crystal_MeV"])
    report.size_table = [
        _size_row(r_nm, vel, mass, d, radius, crystals, props, pulse)
        for r_nm in sorted({5.0, 10.0, 15.0, 25.0, 50.0, 100.0, props.target_radius_nm})
    ]
    report.spectra = {
        "edges_keV": hist_edges.tolist(),
        "D": np.histogram(ke[d], hist_edges)[0].tolist(),
        "Li": np.histogram(ke[~d], hist_edges)[0].tolist(),
    }
    return report


def _scaled_energies(vel, mass, d, scale_r) -> np.ndarray:
    """Kinetic energies (keV) after scaling the crystal radius by ``scale_r``.

    Each species' mean velocity is the laser drift (the explosion is
    symmetric). Only the explosion part (velocity minus drift) grows with
    crystal size, linearly in R at fixed density, so its energy grows as R^2.
    """
    v_scaled = vel.copy()
    for species in (d, ~d):
        if species.any():
            drift = vel[species].mean(axis=0)
            v_scaled[species] = drift + (vel[species] - drift) * scale_r
    return 0.5 * mass * np.einsum("ij,ij->i", v_scaled, v_scaled) / P.J_PER_KEV


def _catcher_yields(ke_d, ke_li, catcher: str, lithium: str, weight: float = 1.0) -> dict:
    """Beam-target reactions when every D and Li ion stops in the catcher."""
    channels: dict[str, float] = {}
    for proj, ke in (("D", ke_d), ("Li", ke_li)):
        if len(ke) == 0:
            continue
        for name, y in P.beam_target(ke, proj, catcher, lithium).items():
            channels[name] = channels.get(name, 0.0) + float(y.sum()) * weight
    energy = sum(v * P.channel_q(k, lithium) for k, v in channels.items())
    return {"channels": channels, "total": sum(channels.values()), "energy_MeV": energy}


def _size_row(r_nm, vel, mass, d, radius, crystals, props, pulse) -> dict:
    """Yields and gain extrapolated to one crystal of radius r_nm, for both catchers.

    The gain uses the laser energy through the crystal at the run's intensity
    and, separately, at the lowest intensity that can strip a crystal of that size.
    """
    scale = r_nm * 1e-9 / radius
    ke = _scaled_energies(vel, mass, d, scale)
    weight = scale**3 / crystals  # ions in one target crystal per simulated ion
    mev = 1.602176634e-13
    area = math.pi * (r_nm * 1e-9) ** 2
    e_strip = P.stripping_field(r_nm * 1e-9)
    p_strip = P.Pulse(e0=e_strip, wavelength=pulse.wavelength, fwhm=pulse.fwhm)
    laser_run = pulse.fluence() * area / mev
    laser_strip = p_strip.fluence() * area / mev
    row = {
        "radius_nm": r_nm,
        "deuterons": int(round(d.sum() * weight)),
        "max_D_keV": float(ke[d].max()),
        "mean_D_keV": float(ke[d].mean()),
        "intensity_to_strip_W_per_cm2": P.intensity_w_cm2(e_strip),
        "stripped_at_run_intensity": pulse.e0 >= e_strip,
        "laser_MeV_at_run_intensity": laser_run,
        "laser_MeV_at_strip_threshold": laser_strip,
    }
    for catcher in ("LiD", "LiT"):
        y = _catcher_yields(ke[d], ke[~d], catcher, props.lithium, weight)
        row[catcher] = {**y, "gain_at_run_intensity": _ratio(y["energy_MeV"], laser_run),
                        "gain_at_strip_threshold": _ratio(y["energy_MeV"], laser_strip)}
    return row


def _ratio(a: float, b: float) -> Optional[float]:
    return a / b if b > 0 else None


def _equivalent_radius(n_ions: int) -> float:
    """Radius of a sphere holding n ions at LiD density."""
    volume = n_ions * P.SPACING**3
    return (3 * volume / (4 * math.pi)) ** (1 / 3)


def _as_result(props, frames, is_li, report) -> SimulationResult:
    n = len(is_li)
    base = SimulatorProperties(name=props.name, lattice_x_size=n, lattice_y_size=1, lattice_z_size=1,
                               iterations=len(frames) - 1)
    # Layer index for the 2D slice view: lattice planes along z.
    z0 = np.array([f for f in frames[0].positions[:, 2]])
    layer = np.round((z0 - z0.min()) / P.SPACING).astype(int)
    index = np.stack([np.zeros(n, int), np.zeros(n, int), layer], axis=1)
    return SimulationResult(
        properties=base,
        types=["Li" if li else "D" for li in is_li],
        lattice_index=index,
        frames=frames,
        fusion_events=None,
        saved=True,
    )
