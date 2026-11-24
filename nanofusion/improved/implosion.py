"""Spherical Coulomb implosion of a hollow fuel shell.

The original proposal hoped the laser would steer ions into head-on
collisions. The strongest form of that idea is spherical focusing: put fuel
ions on the inside wall of a spherical cavity and pull them all towards its
centre at once, so that they meet head-on in a tiny focus. Two published
schemes work this way:

* Microbubble implosion (Murakami et al., Sci. Rep. 8, 7537 (2018)). An
  intense laser fills a sub-micron bubble in solid hydrogen with hot
  electrons, and their charge pulls the bubble-wall ions inwards.
* The periodically oscillating plasma sphere (Barnes & Nebel, Phys. Plasmas
  5, 2498 (1998)). Ions oscillate in the harmonic well of a uniform electron
  cloud and all reach the centre together, every half period. It is a
  refinement of the inertial-electrostatic "fusor".

Model
-----
* The hot electrons are a fixed, uniform negative charge -Q_e filling the
  cavity (radius R). Inside it the field is linear, so an ion that starts at
  rest reaches the centre after pi / (2 omega), with
  omega^2 = z e k Q_e / (m R^3), whatever its starting radius: a perfect lens.
* The fully stripped ions sit in a layer on the cavity wall. Their total
  charge is Q_e / kappa, where kappa is the electron-to-ion charge ratio.
  Large kappa means few ions per electron.
* Every ion feels every other ion, using the same C kernel as the crystal
  model. A run therefore shows what limits the focus: the ions' own repulsion
  (space charge), the graininess of the layer, and their initial random motion.
* Fusion in flight uses the same binary-collision estimator as the crystal
  model. Its reaction volume is smaller here, because the focus can be smaller
  than an atom.

Known limits: the electrons do not move. Real hot electrons would partly follow
the ions into the focus, and would leak away and cool within about a picosecond.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from collections.abc import Callable
from dataclasses import asdict, dataclass

import numpy as np

from . import energy as EN
from . import physics as P
from .native import HIST_EDGES_KEV, NBINS, coulomb_step, native_available

FUELS = {
    # name: (species on even lattice sites, species on odd sites, nuclei per m^3 in the solid, wall material)
    "D2": ("D", "D", 6.0e28, "D2"),
    "DT": ("D", "T", 4.5e28, "DT"),
    "6LiD": ("D", "Li6", 2 * P.lid_molecule_density(), "LiD"),
}


class Cancelled(Exception):
    pass


@dataclass
class ImplosionProperties:
    cavity_nm: float = 6.0  # cavity (bubble) radius
    hot_density: float = 1e28  # hot-electron density inside the cavity, m^-3
    electron_ratio: float = 10.0  # kappa: hot-electron charge / imploding ion charge
    fuel: str = "D2"
    temperature_ev: float = 0.0  # initial random ion motion (thermal or zero-point), eV
    periods: float = 0.5  # simulated time in oscillation periods of the slowest ion (0.5 = one pass)
    softening_nm: float = 0.0005
    reaction_radius_nm: float = 0.05
    seed: int = 1

    def validate(self) -> None:
        problems = []
        for label, value, lo, hi in (
            ("Cavity radius (nm)", self.cavity_nm, 0.5, 40.0),
            ("Hot-electron density (m^-3)", self.hot_density, 1e25, 1e30),
            ("Electron-to-ion charge ratio", self.electron_ratio, 1.0, 1e5),
            ("Temperature (eV)", self.temperature_ev, 0.0, 1e4),
            ("Periods", self.periods, 0.05, 10.0),
            ("Softening (nm)", self.softening_nm, 1e-6, 0.1),
            ("Reaction radius (nm)", self.reaction_radius_nm, 0.005, 1.0),
        ):
            if not lo <= value <= hi:
                problems.append(f"{label} must be between {lo:g} and {hi:g}.")
        if self.fuel not in FUELS:
            problems.append(f"Fuel must be one of {', '.join(FUELS)}.")
        if problems:
            raise ValueError("\n".join(problems))


def build_shell(props: ImplosionProperties) -> dict:
    """Ions in a layer on the cavity wall carrying charge Q_e / kappa.

    Lattice sites (simple cubic at the solid's density) are taken from the
    innermost layers just outside radius R, thickening the layer until it
    holds enough charge, and a random subset (fixed seed) carries exactly the
    needed charge. A partial layer stands for a thinner charged sheath.
    """
    sp_even, sp_odd, n_nuc, _ = FUELS[props.fuel]
    a = n_nuc ** (-1 / 3)
    radius = props.cavity_nm * 1e-9
    electrons = props.hot_density * 4 / 3 * math.pi * radius**3
    ion_charge = electrons / props.electron_ratio
    z_even, z_odd = EN.SPECIES[sp_even][1], EN.SPECIES[sp_odd][1]
    if ion_charge < 2 * max(z_even, z_odd):
        raise ValueError(f"Only {ion_charge:.1f} ion charges: make the cavity larger or the charge ratio smaller.")
    thickness = a
    while True:
        m = int(math.ceil((radius + thickness) / a)) + 1
        rng = np.arange(-m, m + 1)
        idx = np.stack(np.meshgrid(rng, rng, rng, indexing="ij"), axis=-1).reshape(-1, 3)
        r = np.linalg.norm(idx * a, axis=1)
        sel = (r >= radius) & (r < radius + thickness)
        odd = (idx[sel].sum(axis=1) % 2) == 1
        charges = np.where(odd, z_odd, z_even)
        if charges.sum() >= ion_charge:
            break
        thickness += a
    gen = np.random.default_rng(props.seed)
    order = gen.permutation(int(sel.sum()))
    take = int(np.searchsorted(np.cumsum(charges[order]), ion_charge)) + 1
    keep = order[:take]
    pos = (idx[sel] * a).astype(float)[keep]
    species = np.where(odd[keep], sp_odd, sp_even)
    return {"pos": pos, "species": species, "electrons": electrons, "layer_thickness": thickness,
            "layer_fill": take / int(sel.sum()), "lattice_spacing": a}


def _bin_sigma(a: str, b: str) -> np.ndarray:
    """Cross section (m^2) averaged over each pair-exposure bin (uniform in log E)."""
    lo, hi = np.log10(HIST_EDGES_KEV[:-1]), np.log10(HIST_EDGES_KEV[1:])
    pts = lo[:, None] + (hi - lo)[:, None] * (np.arange(16) + 0.5)[None, :] / 16
    return EN.sigma_m2(a, b, 10.0**pts).mean(axis=1)


def run_implosion(props: ImplosionProperties,
                  progress: Callable[[int, int], None] | None = None,
                  cancelled: Callable[[], bool] | None = None) -> dict:
    props.validate()
    started = time.monotonic()
    shell = build_shell(props)
    pos = shell["pos"]
    species = shell["species"]
    n = len(pos)
    z = np.array([EN.SPECIES[s][1] for s in species], dtype=float)
    mass = np.array([P.MASS_U[s] for s in species]) * P.AMU
    q = z * P.E_CHARGE
    kind = (species != "D").astype(np.int32)  # 0 = deuterium, 1 = the other species
    other = FUELS[props.fuel][1]
    radius = props.cavity_nm * 1e-9
    electrons = shell["electrons"]
    k_qe = P.K_COULOMB * electrons * P.E_CHARGE  # V m
    ion_charge = float(z.sum())

    if other.startswith("Li"):
        s_dli = P.S_DLI_KEV_MB[other]
        bg_dli = P.gamow_constant(1, 3, P.reduced_mass_u(P.MASS_U["D"], P.MASS_U[other]))
    else:
        s_dli, bg_dli = 0.0, 1.0  # D+T is taken from the exposure histogram instead
    soft2 = (props.softening_nm * 1e-9) ** 2
    rc = props.reaction_radius_nm * 1e-9
    rc2, inv_vc = rc * rc, 1.0 / (4 / 3 * math.pi * rc**3)

    names = sorted(set(species.tolist()), key=lambda s: s != "D")
    omega = {s: math.sqrt(EN.SPECIES[s][1] * P.E_CHARGE * k_qe / (P.MASS_U[s] * P.AMU * radius**3)) for s in names}
    t_end = props.periods * 2 * math.pi / min(omega.values())
    dt_ref = 2 * math.pi / max(omega.values()) / 2000
    dt_min = dt_ref * 1e-6

    def potential_ext(r):
        return np.where(r < radius, -k_qe * (3 * radius**2 - r**2) / (2 * radius**3), -k_qe / np.maximum(r, 1e-30))

    def accel(x, v):
        f, stats, h_dd, h_x = coulomb_step(x, v, q, mass, kind, soft2, rc2, inv_vc, s_dli, bg_dli)
        r = np.sqrt(np.einsum("ij,ij->i", x, x))
        inv3 = np.where(r < radius, 1.0 / radius**3, 1.0 / np.maximum(r, 1e-30) ** 3)
        a = (f - (q * k_qe * inv3)[:, None] * x) / mass[:, None]
        return a, stats, h_dd, h_x, r

    vel = np.zeros_like(pos)
    if props.temperature_ev > 0:
        gen = np.random.default_rng(props.seed + 1)
        vel = gen.normal(size=pos.shape) * np.sqrt(props.temperature_ev * P.E_CHARGE / mass)[:, None]

    sig_dd = _bin_sigma("D", "D")
    sig_x = _bin_sigma("D", other) if other != "D" else np.zeros(NBINS)
    q_dd, q_x = EN.q_value("D", "D"), EN.q_value("D", other) if other != "D" else 0.0

    t = 0.0
    a, stats, h_dd, h_x, r = accel(pos, vel)
    ke_i = 0.5 * mass * np.einsum("ij,ij->i", vel, vel)
    e_total0 = float(ke_i.sum() + stats[5] + np.sum(q * potential_ext(r)))
    x_dd = np.zeros(NBINS)
    x_x = np.zeros(NBINS)
    y_dd_exact = y_dli_exact = 0.0
    max_ke, max_ke_t = float(ke_i.sum()), 0.0
    went_in = np.zeros(n, dtype=bool)
    wall_ke = np.full(n, np.nan)
    groups = {s: species == s for s in names}
    rec: dict[str, list] = {"t": [], "yield": [], "ke": []}
    for s in names:
        rec[f"median_{s}"] = []
        rec[f"rms_{s}"] = []
        rec[f"aim_{s}"] = []
    steps = 0

    def record():
        rec["t"].append(t)
        rec["yield"].append(float(x_dd @ sig_dd + x_x @ sig_x))
        rec["ke"].append(float(ke_i.sum()))
        speed = np.sqrt(np.einsum("ij,ij->i", vel, vel))
        moving = speed > 0
        aim = np.zeros(n)
        if moving.any():
            aim[moving] = np.linalg.norm(np.cross(pos[moving], vel[moving]), axis=1) / speed[moving]
        for s, g in groups.items():
            rec[f"median_{s}"].append(float(np.median(r[g])))
            rec[f"rms_{s}"].append(float(np.sqrt(np.mean(r[g] ** 2))))
            rec[f"aim_{s}"].append(float(np.median(aim[g])))

    record()
    while t < t_end * (1 - 1e-12):
        if cancelled is not None and cancelled():
            raise Cancelled()
        amax = float(np.sqrt(np.max(np.einsum("ij,ij->i", a, a))))
        dt = dt_ref
        if stats[4] < 1e299:
            dt = min(dt, 0.02 * stats[4])
        if amax > 0:
            dt = min(dt, math.sqrt(0.002 * min(stats[6], shell["lattice_spacing"]) / amax))
        dt = max(dt_min, min(dt, t_end - t))
        v_half = vel + 0.5 * dt * a
        pos = pos + dt * v_half
        t += dt
        a, stats, h_dd, h_x, r = accel(pos, v_half)
        vel = v_half + 0.5 * dt * a
        x_dd += h_dd * dt
        x_x += h_x * dt
        y_dd_exact += stats[1] * dt
        y_dli_exact += stats[0] * dt
        ke_i = 0.5 * mass * np.einsum("ij,ij->i", vel, vel)
        if ke_i.sum() > max_ke:
            max_ke, max_ke_t = float(ke_i.sum()), t
        went_in |= r < 0.9 * radius
        out = went_in & np.isnan(wall_ke) & (r >= radius) & (np.einsum("ij,ij->i", pos, vel) > 0)
        wall_ke[out] = ke_i[out]
        record()
        steps += 1
        if progress is not None:
            progress(int(1000 * t / t_end), 1000)

    e_total1 = float(ke_i.sum() + stats[5] + np.sum(q * potential_ext(r)))
    rec_np = {k: np.array(v) for k, v in rec.items()}
    kev = P.J_PER_KEV
    mev_j = 1.602176634e-13

    # Focus passes: the k-th focus of species s is due at (2k + 1) pi / (2 omega_s).
    foci = []
    for s in names:
        period = 2 * math.pi / omega[s]
        k = 0
        while (2 * k + 1) * period / 4 < t_end:
            t_due = (2 * k + 1) * period / 4
            win = (rec_np["t"] >= t_due - period / 4) & (rec_np["t"] <= min(t_due + period / 4, t_end))
            if win.sum() >= 3:
                med = rec_np[f"median_{s}"][win]
                i = int(np.argmin(med))
                y = rec_np["yield"][win]
                foci.append({
                    "species": s, "pass": k + 1, "due_fs": t_due * 1e15,
                    "time_fs": float(rec_np["t"][win][i] * 1e15),
                    "median_radius_nm": float(med[i] * 1e9),
                    "rms_radius_nm": float(rec_np[f"rms_{s}"][win].min() * 1e9),
                    "median_radius_over_R": float(med[i] / radius),
                    "reactions_in_window": float(y[-1] - y[0]),
                })
            k += 1

    # Impact-parameter spread on the way in, when the ions are halfway to the centre.
    aim = {}
    for s in names:
        half = np.nonzero(rec_np[f"median_{s}"] <= 0.5 * radius)[0]
        aim[s] = float(rec_np[f"aim_{s}"][half[0]] * 1e9) if len(half) else None

    channels = {"D+D": float(x_dd @ sig_dd)}
    if other != "D":
        channels[f"D+{other}"] = float(x_x @ sig_x)
    fusion_mev = channels["D+D"] * q_dd + (channels.get(f"D+{other}", 0.0) * q_x)
    n_d = int((species == "D").sum())
    column = {"D sees D": 2 * float(x_dd.sum()) / max(n_d, 1)}
    if other != "D":
        column[f"D sees {other}"] = float(x_x.sum()) / max(n_d, 1)

    # Ions that come back out to the wall with energy left hit cold fuel there (beam-target).
    wall_edges = HIST_EDGES_KEV
    wall_spectrum = {s: np.histogram(wall_ke[(species == s) & ~np.isnan(wall_ke)] / P.J_PER_KEV, wall_edges)[0].tolist()
                     for s in names}
    wall = FUELS[props.fuel][3]
    _, _, n_nuc, _ = FUELS[props.fuel]
    targets = {"D": n_nuc} if other == "D" else {"D": n_nuc / 2, other: n_nuc / 2}
    wall_reactions = 0.0
    wall_mev = 0.0
    for s in names:
        e = wall_ke[(species == s) & ~np.isnan(wall_ke)] / kev
        if len(e):
            res = EN.thick_target(e, s, wall, {k: v for k, v in targets.items() if EN.q_value(s, k) > 0})
            wall_reactions += float(res["reactions"].sum())
            wall_mev += float(res["fusion_MeV"].sum())

    field_min_j = P.K_COULOMB * (electrons * P.E_CHARGE) ** 2 / (10 * radius)
    total_mev = fusion_mev + wall_mev
    kappa = electrons / ion_charge
    z_mean = ion_charge / n
    bounce = (math.sqrt(1 + 4 / kappa) - 1) / 2  # thin uniform shell, no graininess
    every = max(1, len(rec_np["t"]) // 400)
    report = {
        "inputs": asdict(props),
        "ions": n,
        "species": {s: int(g.sum()) for s, g in groups.items()},
        "cavity_nm": props.cavity_nm,
        "hot_electrons": electrons,
        "electron_ratio": kappa,
        "layer_thickness_nm": shell["layer_thickness"] * 1e9,
        "layer_fill": shell["layer_fill"],
        "native_kernel": native_available(),
        "steps": steps,
        "wall_seconds": time.monotonic() - started,
        "simulated_fs": t * 1e15,
        "theory": {
            "focus_time_fs": {s: math.pi / (2 * omega[s]) * 1e15 for s in names},
            "period_fs": {s: 2 * math.pi / omega[s] * 1e15 for s in names},
            "centre_energy_keV": {s: EN.SPECIES[s][1] * P.E_CHARGE * k_qe / (2 * radius) / kev for s in names},
            "thin_shell_bounce_radius_over_R": bounce,
            "note": "centre_energy: kinetic energy at the centre for an ion starting at rest on the wall, "
                    "ignoring the other ions; bounce radius: where a smooth thin shell of the same charge "
                    "would stop (space-charge limit), mean charge " + f"{z_mean:.2f}",
        },
        "focus": foci,
        "dynamics": {
            "max_total_ion_KE_MeV": max_ke / mev_j,
            "max_total_ion_KE_time_fs": max_ke_t * 1e15,
            "mean_ion_KE_at_max_keV": max_ke / n / kev,
            "median_impact_parameter_halfway_nm": aim,
            "energy_drift_relative": (e_total1 - e_total0) / max(abs(e_total0), 1e-300),
        },
        "fusion": {
            "in_flight": channels,
            "in_flight_exact_D_D": y_dd_exact,
            "in_flight_exact_D_Li": y_dli_exact if other.startswith("Li") else None,
            "column_per_ion_m2": column,
            "wall_beam_target_reactions": wall_reactions,
            "wall_beam_target_MeV": wall_mev,
            "wall_spectrum": {"edges_keV": wall_edges.tolist(), **wall_spectrum,
                              "note": "kinetic energy of ions when they come back out through the cavity wall"},
            "fusion_MeV": total_mev,
        },
        "energy": {
            "ion_energy_MeV": max_ke / mev_j,
            "field_energy_min_MeV": field_min_j / mev_j,
            "gain_vs_ion_energy": total_mev / (max_ke / mev_j) if max_ke > 0 else None,
            "gain_vs_field_energy": total_mev / (field_min_j / mev_j),
            "note": "field_energy_min: electrostatic energy of the hot-electron cloud with its compensating "
                    "positive charge right at the wall, k Q_e^2 / (10 R); the laser must supply at least this",
        },
        "pair_exposure": {
            "edges_keV": HIST_EDGES_KEV.tolist(),
            "D_D": x_dd.tolist(),
            "D_X": x_x.tolist(),
            "X": other,
        },
        "series": {k: (v[::every] * (1e15 if k == "t" else 1e9 if k.startswith(("median", "rms", "aim"))
                                     else 1.0 / mev_j if k == "ke" else 1.0)).tolist()
                   for k, v in rec_np.items()},
        "series_units": {"t": "fs", "median_*/rms_*/aim_*": "nm", "ke": "MeV", "yield": "reactions (cumulative)"},
    }
    return report


def main(argv: list[str] | None = None) -> int:
    d = ImplosionProperties()
    p = argparse.ArgumentParser(prog="python -m nanofusion implosion",
                                description="Coulomb implosion of a hollow fuel shell onto its centre.")
    p.add_argument("--cavity", type=float, default=d.cavity_nm, help="cavity radius, nm")
    p.add_argument("--hot-density", type=float, default=d.hot_density, help="hot-electron density, m^-3")
    p.add_argument("--ratio", type=float, default=d.electron_ratio, help="hot-electron charge / ion charge")
    p.add_argument("--fuel", choices=list(FUELS), default=d.fuel)
    p.add_argument("--temperature", type=float, default=d.temperature_ev, help="initial ion temperature, eV")
    p.add_argument("--periods", type=float, default=d.periods, help="oscillation periods to simulate")
    p.add_argument("--reaction-radius", type=float, default=d.reaction_radius_nm, help="nm")
    p.add_argument("--seed", type=int, default=d.seed)
    p.add_argument("-o", "--output", help="write the report as JSON")
    args = p.parse_args(argv)
    props = ImplosionProperties(cavity_nm=args.cavity, hot_density=args.hot_density, electron_ratio=args.ratio,
                                fuel=args.fuel, temperature_ev=args.temperature, periods=args.periods,
                                reaction_radius_nm=args.reaction_radius, seed=args.seed)
    try:
        rep = run_implosion(props)
    except ValueError as exc:
        print(exc, file=sys.stderr)
        return 2
    if args.output:
        with open(args.output, "w") as fh:
            json.dump(rep, fh, indent=1)
    print(f"{rep['ions']} ions {rep['species']}, kappa {rep['electron_ratio']:.3g}, "
          f"{rep['steps']} steps, {rep['wall_seconds']:.1f} s, energy drift {rep['dynamics']['energy_drift_relative']:.1e}")
    th = rep["theory"]
    print(f"centre energy (no space charge) {th['centre_energy_keV']}, bounce radius/R {th['thin_shell_bounce_radius_over_R']:.3g}")
    for f in rep["focus"]:
        print(f"  focus {f['species']} #{f['pass']}: t={f['time_fs']:.1f} fs (due {f['due_fs']:.1f}), "
              f"median r/R={f['median_radius_over_R']:.3g}, reactions {f['reactions_in_window']:.3g}")
    fu, en = rep["fusion"], rep["energy"]
    print(f"in flight {fu['in_flight']} (exact D+D {fu['in_flight_exact_D_D']:.3g}); wall {fu['wall_beam_target_reactions']:.3g}")
    print(f"column per ion {fu['column_per_ion_m2']}")
    print(f"ion energy {en['ion_energy_MeV']:.3g} MeV, gain vs ions {en['gain_vs_ion_energy']:.3g}, "
          f"vs field {en['gain_vs_field_energy']:.3g}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
