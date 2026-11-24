"""Hollow-shell Coulomb implosions: can spherical focusing make fusion pay?

Runs nanofusion.improved.implosion over:
  A. the electron-to-ion charge ratio kappa (how hard the ions' own repulsion stops the focus);
  B. cavity size at fixed kappa, to test the scaling rule used to extrapolate;
  C. initial ion temperature;
  D. fuel (D2, DT, 6LiD);
  E. several oscillation passes;
  F. the reaction-volume radius (a numerical check);
then writes docs/results/implosion/summary.md.

    python3 tools/implosion_experiments.py [--analyse-only] [--force]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import energy as EN
from nanofusion.improved import physics as P
from nanofusion.improved.implosion import (
    FUELS,
    ImplosionProperties,
    run_implosion,
)

OUT = Path("docs/results/implosion")
HOT = 1e28  # hot-electron density, m^-3


def bounce(kappa: float) -> float:
    return (math.sqrt(1 + 4 / kappa) - 1) / 2


def rc_for(cavity: float, kappa: float) -> float:
    """Reaction-volume radius (nm): a quarter of the expected focus radius, 0.005-0.05 nm."""
    return min(0.05, max(0.005, 0.25 * bounce(kappa) * cavity))


def runs() -> dict[str, ImplosionProperties]:
    r: dict[str, ImplosionProperties] = {}
    for k in (2, 5, 10, 30, 100, 300, 1000):
        r[f"kappa{k}-R8-D2"] = ImplosionProperties(cavity_nm=8, electron_ratio=k, reaction_radius_nm=rc_for(8, k))
    for size in (4, 6, 9, 13):
        r[f"kappa10-R{size}-D2"] = ImplosionProperties(cavity_nm=size, electron_ratio=10,
                                                       reaction_radius_nm=rc_for(size, 10))
    for size in (6, 9, 13, 18):
        r[f"kappa100-R{size}-D2"] = ImplosionProperties(cavity_nm=size, electron_ratio=100,
                                                        reaction_radius_nm=rc_for(size, 100))
    for t in (0.01, 0.1, 1.0, 10.0):
        r[f"temp{t:g}-kappa100-R8-D2"] = ImplosionProperties(cavity_nm=8, electron_ratio=100, temperature_ev=t,
                                                             reaction_radius_nm=rc_for(8, 100))
    for fuel in ("DT", "6LiD"):
        r[f"kappa30-R8-{fuel}"] = ImplosionProperties(cavity_nm=8, electron_ratio=30, fuel=fuel,
                                                      reaction_radius_nm=rc_for(8, 30))
    r["passes-kappa100-R8-D2"] = ImplosionProperties(cavity_nm=8, electron_ratio=100, periods=2.0,
                                                     reaction_radius_nm=rc_for(8, 100))
    for rc in (0.01, 0.05):
        r[f"rc{rc:g}-kappa100-R8-D2"] = ImplosionProperties(cavity_nm=8, electron_ratio=100, reaction_radius_nm=rc)
    return r


def run_all(force: bool) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, props in runs().items():
        path = OUT / f"{name}.json"
        if path.exists() and not force:
            continue
        props.hot_density = HOT
        t0 = time.monotonic()
        rep = run_implosion(props)
        path.write_text(json.dumps(rep, indent=1))
        print(f"{name}: {rep['ions']} ions, {rep['steps']} steps, {time.monotonic() - t0:.0f} s", flush=True)


# --- analysis -------------------------------------------------------------------

def load(name: str) -> dict | None:
    path = OUT / f"{name}.json"
    return json.loads(path.read_text()) if path.exists() else None


def bin_sigma(a: str, b: str, edges: np.ndarray) -> np.ndarray:
    lo, hi = np.log10(edges[:-1]), np.log10(edges[1:])
    pts = lo[:, None] + (hi - lo)[:, None] * (np.arange(16) + 0.5)[None, :] / 16
    return EN.sigma_m2(a, b, 10.0**pts).mean(axis=1)


def scaled(rep: dict, s: float) -> dict:
    """The same implosion with every length times s at the same hot-electron density and kappa.

    Ion numbers grow as s^3, speeds as s (energies s^2), times are unchanged, so
    each pair-exposure bin moves to s^2 times its energy and its weight grows as
    s^4 (s^3 more ions, each crossing s times more matter). Field and ion
    energies grow as s^5. The wall ions' energies scale the same way.
    """
    pe = rep["pair_exposure"]
    edges = np.array(pe["edges_keV"]) * s**2
    other = pe["X"]
    x_dd = np.array(pe["D_D"]) * s**4
    x_x = np.array(pe["D_X"]) * s**4
    channels = {"D+D": float(x_dd @ bin_sigma("D", "D", edges))}
    fusion = channels["D+D"] * EN.q_value("D", "D")
    if other != "D":
        channels[f"D+{other}"] = float(x_x @ bin_sigma("D", other, edges))
        fusion += channels[f"D+{other}"] * EN.q_value("D", other)
    fuel = rep["inputs"]["fuel"]
    _, _, n_nuc, wall = FUELS[fuel]
    targets = {"D": n_nuc} if other == "D" else {"D": n_nuc / 2, other: n_nuc / 2}
    ws = rep["fusion"]["wall_spectrum"]
    wedges = np.array(ws["edges_keV"]) * s**2
    centres = np.sqrt(wedges[1:] * wedges[:-1])
    wall_mev = 0.0
    for sp, counts in ws.items():
        if sp in ("edges_keV", "note"):
            continue
        counts = np.array(counts, dtype=float) * s**3
        if counts.sum() == 0:
            continue
        res = EN.thick_target(centres, sp, wall, {k: v for k, v in targets.items() if EN.q_value(sp, k) > 0})
        wall_mev += float(counts @ res["fusion_MeV"])
    ion_mev = rep["energy"]["ion_energy_MeV"] * s**5
    field_mev = rep["energy"]["field_energy_min_MeV"] * s**5
    total = fusion + wall_mev
    return {"R_nm": rep["cavity_nm"] * s, "ions": rep["ions"] * s**3,
            "centre_keV": {k: v * s**2 for k, v in rep["theory"]["centre_energy_keV"].items()},
            "mean_keV": rep["dynamics"]["mean_ion_KE_at_max_keV"] * s**2,
            "channels": channels, "in_flight_MeV": fusion, "wall_MeV": wall_mev,
            "gain_ion": total / ion_mev, "gain_field": total / field_mev, "field_J": field_mev * 1.602176634e-13}


def ideal_dt(rep_d2: dict) -> dict:
    """A D-T shell whose D and T somehow reach the centre together (best case).

    Takes the pure-deuterium run's dynamics and makes half the ions tritium:
    D-T pairs are half of all pairs, D-D pairs a quarter, with the same
    relative energies (head-on D-T at equal ion energies has about the same
    centre-of-mass energy as head-on D-D).
    """
    rep = json.loads(json.dumps(rep_d2))
    x = np.array(rep_d2["pair_exposure"]["D_D"])
    rep["pair_exposure"]["D_D"] = (0.25 * x).tolist()
    rep["pair_exposure"]["D_X"] = (0.5 * x).tolist()
    rep["pair_exposure"]["X"] = "T"
    rep["inputs"]["fuel"] = "DT"
    return rep


def gain_bounds() -> list[str]:
    """Space-charge bound on the gain of a non-neutral focus, from the cross sections alone."""
    e = np.geomspace(1, 4900, 2000)
    k = P.EPS0 * 1.602176634e-13 / P.E_CHARGE**2  # m per (m^2 MeV)
    s_dt = float(np.max(EN.sigma_m2("D", "T", e)))
    s_dd = float(np.max(EN.sigma_m2("D", "D", e)))
    l_dt = k * s_dt * EN.q_value("D", "T")
    l_dd = 2 * k * s_dd * EN.q_value("D", "D")
    return ["## Space-charge bound", "",
            "Only so many ions fit in a focus of radius r_c before their own repulsion stops the rest: "
            "N e²/(4πε₀ r_c) ≤ E. So the matter an ion crosses is at most 4ε₀E/(e² r_c), and the fusion "
            "energy per unit of ion energy is at most",
            "",
            f"* D₂: 2ε₀σQ/(e² r_c) = {l_dd * 1e15:.3g} fm / r_c (σ at its maximum, {s_dd / 1e-28:.2g} b);",
            f"* D–T (50/50, arriving together): ε₀σQ/(e² r_c) = {l_dt * 1e12:.3g} pm / r_c "
            f"(σ = {s_dt / 1e-28:.2g} b at 64 keV).",
            "",
            "Breaking even on the ions' own energy needs a focus smaller than these lengths; an atom is about "
            "50,000 fm (50 pm) across.", ""]


def max_wall_kev(rep: dict) -> float:
    """Highest energy (keV) of any ion coming back out through the wall."""
    ws = rep["fusion"]["wall_spectrum"]
    edges = np.array(ws["edges_keV"])
    top = 0.0
    for sp, counts in ws.items():
        if sp in ("edges_keV", "note"):
            continue
        nz = np.nonzero(np.array(counts))[0]
        if len(nz):
            top = max(top, float(edges[nz.max() + 1]))
    return max(top, max(rep["theory"]["centre_energy_keV"].values()))


def g(x: float) -> str:
    return f"{x:.2g}" if x is not None else "–"


def focus_of(rep: dict, species: str = "D", k: int = 1) -> dict:
    for f in rep["focus"]:
        if f["species"] == species and f["pass"] == k:
            return f
    return {}


def column_formula(rep: dict) -> float:
    """Space-charge bound on the matter an ion crosses: 4 eps0 E / (e^2 r_c)."""
    e = rep["dynamics"]["mean_ion_KE_at_max_keV"] * P.J_PER_KEV
    rc = focus_of(rep)["median_radius_nm"] * 1e-9
    return 4 * P.EPS0 * e / (P.E_CHARGE**2 * rc)


def analyse() -> str:
    out = ["# Hollow-shell Coulomb implosions", "",
           f"Hot-electron density {HOT:.0e} m⁻³ throughout. Produced by `tools/implosion_experiments.py`.", ""]

    # A. kappa scan
    out += ["## A. Electron-to-ion charge ratio κ (D₂, cavity radius 8 nm)", "",
            "| κ | ions | energy at centre, no space charge (keV) | mean ion energy reached (keV) "
            "| focus radius / R (run) | smooth-shell bounce / R | column per ion (m⁻²) | column / space-charge bound "
            "| D+D in flight | D+D when the ions spring back into the wall | gain vs ion energy | gain vs field energy |",
            "|---|---|---|---|---|---|---|---|---|---|---|---|"]
    for k in (2, 5, 10, 30, 100, 300, 1000):
        rep = load(f"kappa{k}-R8-D2")
        if rep is None:
            continue
        f = focus_of(rep)
        col = rep["fusion"]["column_per_ion_m2"]["D sees D"]
        if col == 0:  # no pair ever came within the reaction radius: too few ions to sample the focus
            out.append(f"| {k} | {rep['ions']} | {rep['theory']['centre_energy_keV']['D']:.3g} "
                       f"| {rep['dynamics']['mean_ion_KE_at_max_keV']:.3g} "
                       f"| {focus_of(rep).get('median_radius_over_R', float('nan')):.3g} "
                       f"| {rep['theory']['thin_shell_bounce_radius_over_R']:.3g} | too few ions to sample | | | | | |")
            continue
        out.append(f"| {k} | {rep['ions']} | {rep['theory']['centre_energy_keV']['D']:.3g} "
                   f"| {rep['dynamics']['mean_ion_KE_at_max_keV']:.3g} | {f.get('median_radius_over_R', float('nan')):.3g} "
                   f"| {rep['theory']['thin_shell_bounce_radius_over_R']:.3g} | {col:.2g} | {col / column_formula(rep):.2f} "
                   f"| {rep['fusion']['in_flight_exact_D_D']:.2g} | {rep['fusion']['wall_beam_target_reactions']:.2g} "
                   f"| {g(rep['energy']['gain_vs_ion_energy'])} | {g(rep['energy']['gain_vs_field_energy'])} |")
    out.append("")

    # B. size scaling
    for kappa, sizes in ((10, (4, 6, 9, 13)), (100, (6, 9, 13, 18))):
        reps = {s: load(f"kappa{kappa}-R{s}-D2") for s in sizes}
        reps = {s: r for s, r in reps.items() if r is not None}
        if len(reps) < 2:
            continue
        out += [f"## B. Size scaling, κ = {kappa} (D₂)", "",
                "| R (nm) | ions | mean ion energy (keV) | focus radius / R | D+D in flight (run) "
                "| D+D predicted from the next smaller run |", "|---|---|---|---|---|---|"]
        rs = sorted(reps)
        for i, s in enumerate(rs):
            rep = reps[s]
            pred = "–"
            if i > 0:
                prev = reps[rs[i - 1]]
                pred = g(scaled(prev, s / rs[i - 1])["channels"]["D+D"])
            out.append(f"| {s} | {rep['ions']} | {rep['dynamics']['mean_ion_KE_at_max_keV']:.3g} "
                       f"| {focus_of(rep).get('median_radius_over_R', float('nan')):.3g} "
                       f"| {g(scaled(rep, 1.0)['channels']['D+D'])} | {pred} |")
        out.append("")
        base = reps[rs[-1]]
        out += [f"Extrapolated from the R = {rs[-1]} nm run (same κ and hot-electron density):", "",
                "| R (nm) | ions | energy at centre (keV) | D+D in flight | fusion energy from ions springing back into "
                "the wall | gain vs ion energy | gain vs field energy | field energy (J) |", "|---|---|---|---|---|---|---|---|"]
        for target in (20, 30, 50, 75, 100, 150, 200, 300):
            sc = scaled(base, target / rs[-1])
            share = sc["wall_MeV"] / (sc["wall_MeV"] + sc["in_flight_MeV"])
            out.append(f"| {target} | {sc['ions']:.2g} | {sc['centre_keV']['D']:.3g} | {g(sc['channels']['D+D'])} "
                       f"| {share:.0%} | {g(sc['gain_ion'])} | {g(sc['gain_field'])} | {sc['field_J']:.2g} |")
        out.append("")

    # C. temperature
    out += ["## C. Initial ion temperature (D₂, R = 8 nm, κ = 100)", "",
            "| temperature (eV) | focus radius / R | D+D in flight | gain vs ion energy |", "|---|---|---|---|"]
    for t, name in ((0.0, "kappa100-R8-D2"), (0.01, "temp0.01-kappa100-R8-D2"), (0.1, "temp0.1-kappa100-R8-D2"),
                    (1.0, "temp1-kappa100-R8-D2"), (10.0, "temp10-kappa100-R8-D2")):
        rep = load(name)
        if rep is None:
            continue
        out.append(f"| {t:g} | {focus_of(rep).get('median_radius_over_R', float('nan')):.3g} "
                   f"| {rep['fusion']['in_flight_exact_D_D']:.2g} | {g(rep['energy']['gain_vs_ion_energy'])} |")
    out.append("")

    # D. fuels
    out += ["## D. Fuel (R = 8 nm, κ = 30)", "",
            "| fuel | ions | D focus time (fs) | partner focus time (fs) | column of D seen by D (m⁻²) "
            "| column of partner seen by D (m⁻²) | reactions in flight | gain vs ion energy |",
            "|---|---|---|---|---|---|---|---|"]
    for fuel in ("D2", "DT", "6LiD"):
        rep = load("kappa30-R8-D2" if fuel == "D2" else f"kappa30-R8-{fuel}")
        if rep is None:
            continue
        other = rep["pair_exposure"]["X"]
        col = rep["fusion"]["column_per_ion_m2"]
        fd, fo = focus_of(rep, "D"), focus_of(rep, other)
        out.append(f"| {fuel} | {rep['ions']} | {fd.get('time_fs', float('nan')):.1f} "
                   f"| {fo.get('time_fs', float('nan')):.1f} | {col['D sees D']:.2g} "
                   f"| {g(col.get(f'D sees {other}'))} | {sum(rep['fusion']['in_flight'].values()):.2g} "
                   f"| {g(rep['energy']['gain_vs_ion_energy'])} |")
    out.append("")
    out += ["Fuels extrapolated from their R = 8 nm, κ = 30 runs:", "",
            "| fuel | R (nm) | energy at centre, D (keV) | reactions in flight | share from ions springing back "
            "into the wall | gain vs ion energy | gain vs field energy |",
            "|---|---|---|---|---|---|---|"]
    best = {}
    d2 = load("kappa30-R8-D2")
    for fuel in ("D2", "DT", "6LiD", "DT, arriving together"):
        if fuel == "DT, arriving together":
            rep = ideal_dt(d2) if d2 is not None else None
        else:
            rep = load("kappa30-R8-D2" if fuel == "D2" else f"kappa30-R8-{fuel}")
        if rep is None:
            continue
        for target in (20, 30, 50, 100, 200):
            sc = scaled(rep, target / 8)
            share = sc["wall_MeV"] / (sc["wall_MeV"] + sc["in_flight_MeV"])
            out.append(f"| {fuel} | {target} | {sc['centre_keV']['D']:.3g} | {g(sum(sc['channels'].values()))} "
                       f"| {share:.0%} | {g(sc['gain_ion'])} | {g(sc['gain_field'])} |")
        # Search sizes only while every ion stays within the measured cross-section range (4.9 MeV).
        r_max = 8 * math.sqrt(4900 / max_wall_kev(rep))
        grid = np.geomspace(8, r_max, 200)
        gains = [(scaled(rep, t / 8), t) for t in grid]
        sc, r_best = max(gains, key=lambda x: x[0]["gain_field"])
        best[fuel] = (sc, r_best)
    out.append("")
    if best:
        out.append("Best gain against the field energy, for cavities small enough that no ion exceeds 4.9 MeV "
                   "(the end of the measured cross sections):")
        out.append("")
        out.append("| fuel | best R (nm) | gain vs field energy | gain vs ion energy | of which from the focus "
                   "| share from ions springing back into the wall | fusion reactions per implosion |")
        out.append("|---|---|---|---|---|---|---|")
        for f, (sc, r) in best.items():
            share = sc["wall_MeV"] / (sc["wall_MeV"] + sc["in_flight_MeV"])
            out.append(f"| {f} | {r:.0f} | {sc['gain_field']:.2g} | {sc['gain_ion']:.2g} "
                       f"| {sc['gain_ion'] * (1 - share):.2g} | {share:.0%} "
                       f"| {g(sum(sc['channels'].values()))} in flight |")
        out.append("")

    out += gain_bounds()

    # E. passes
    rep = load("passes-kappa100-R8-D2")
    if rep is not None:
        out += ["## E. Repeated passes (D₂, R = 8 nm, κ = 100, two oscillation periods)", "",
                "| pass | time (fs) | focus radius / R | D+D in this pass |", "|---|---|---|---|"]
        for f in rep["focus"]:
            out.append(f"| {f['pass']} | {f['time_fs']:.1f} | {f['median_radius_over_R']:.3g} | {f['reactions_in_window']:.2g} |")
        out.append("")

    # F. reaction radius
    out += ["## F. Numerical check: reaction-volume radius (D₂, R = 8 nm, κ = 100)", "",
            "| radius (nm) | D+D in flight (exact estimator) | D+D from the exposure histogram |", "|---|---|---|"]
    for rc, name in ((0.01, "rc0.01-kappa100-R8-D2"), (rc_for(8, 100), "kappa100-R8-D2"), (0.05, "rc0.05-kappa100-R8-D2")):
        rep = load(name)
        if rep is None:
            continue
        out.append(f"| {rc:.3g} | {rep['fusion']['in_flight_exact_D_D']:.2g} | {rep['fusion']['in_flight']['D+D']:.2g} |")
    out.append("")
    return "\n".join(out) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--analyse-only", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if not args.analyse_only:
        run_all(args.force)
    text = analyse()
    (OUT / "summary.md").write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
