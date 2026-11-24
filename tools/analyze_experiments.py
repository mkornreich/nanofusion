"""Analyse tools/experiments.py output: does fusion *inside* the crystal grow with size?

For overtaking ions (shock shells, or deuterons overrunning lithium) the
relative speed of a collision grows in proportion to the crystal radius R,
because all explosion velocities do. Collisions between ordinary neighbours do
not grow. So the high-energy tail of the pair-exposure histogram should shift
as E -> s^2 E when the crystal grows by s. Its weight grows as s^4: s^3 more
ions, times s in relative speed.

This script:
1. measures how the hard tail moves between sizes (the energy above which the
   top 1% of exposure lies, against R);
2. checks the scaling rule by predicting the 3 nm in-crystal yield from the
   2 nm run;
3. if that works, extrapolates in-crystal yields to large crystals and puts
   them next to the beam-target yields.

    python3 tools/analyze_experiments.py [--dir docs/results/experiments]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import physics as P  # noqa: E402

SIZES = (1.0, 1.5, 2.0, 3.0)
TARGETS_NM = (5.0, 15.0, 30.0, 50.0, 100.0)


def sigma_m2(channel: str, lithium: str, e_kev: np.ndarray) -> np.ndarray:
    if channel == "D_D":
        return P.sigma_dd_mb(e_kev) * P.MB_TO_M2
    mu = P.reduced_mass_u(P.MASS_U["D"], P.MASS_U[lithium])
    return P.sigma_dli_mb(e_kev, P.S_DLI_KEV_MB[lithium], P.gamow_constant(1, 3, mu)) * P.MB_TO_M2


def tail_energy(edges: np.ndarray, h: np.ndarray, top: float = 0.01) -> float:
    """Energy above which the top fraction of (relative-speed-weighted) exposure lies."""
    c = np.cumsum(h[::-1])[::-1]
    if c[0] <= 0:
        return 0.0
    idx = np.nonzero(c >= top * c[0])[0][-1]
    return float(edges[idx])


def tail_yield(rep: dict, channel: str, lithium: str, scale: float, tail_from_kev: float) -> float:
    """In-crystal yield of a crystal `scale` times larger, scaling only bins above tail_from_kev."""
    edges = np.array(rep["pair_exposure"]["edges_keV"])
    centres = np.sqrt(edges[1:] * edges[:-1])
    h = np.array(rep["pair_exposure"][channel])
    crystals = 2 if rep["inputs"].get("geometry") == "pair" else 1
    tail = centres >= tail_from_kev
    scaled = scale**4 * np.sum(h[tail] * sigma_m2(channel, lithium, centres[tail] * scale**2)) / crystals
    rest = np.sum(h[~tail] * sigma_m2(channel, lithium, centres[~tail])) * scale**3 / crystals
    return float(scaled + rest)


def analyse_pairs(root: Path) -> list[str]:
    """Colliding crystals: two crystals of radius R, gap 2.5 R, Li-6, 1e20 W/cm^2."""
    files = {1.0: "pair-R1-gap2.5", 1.5: "pair-R1.5-gap3.75", 2.0: "pair-gap5"}
    reps = {r: json.loads((root / f"{n}.json").read_text()) for r, n in files.items() if (root / f"{n}.json").exists()}
    if len(reps) < 3:
        return []
    lithium = "Li6"
    edges = np.array(reps[2.0]["pair_exposure"]["edges_keV"])
    out = ["## Colliding crystals (two Li-6 crystals, gap = 2.5 R)", "",
           "| R (nm) | ions | hardest D–Li pair (keV) | top-1% D–D exposure above (keV) "
           "| top-1% D–Li exposure above (keV) | collision D+D | collision D+Li |", "|---|---|---|---|---|---|---|"]
    for r, rep in sorted(reps.items()):
        pe = rep["pair_exposure"]
        out.append(f"| {r:g} | {rep['particles']} | {rep['dynamics']['max_D_Li_pair_energy_keV']:.3g} "
                   f"| {tail_energy(edges, np.array(pe['D_D'])):.3g} | {tail_energy(edges, np.array(pe['D_Li'])):.3g} "
                   f"| {rep['fusion']['in_crystal_D_D']:.2g} | {rep['fusion']['in_crystal_D_Li']:.2g} |")
    alphas = {}
    for ch in ("D_D", "D_Li"):
        e1 = tail_energy(edges, np.array(reps[1.0]["pair_exposure"][ch]))
        e2 = tail_energy(edges, np.array(reps[2.0]["pair_exposure"][ch]))
        alphas[ch] = math.log(e2 / e1) / math.log(2.0) if e1 > 0 and e2 > 0 else float("nan")
    out += ["", f"Tail energy grows as R^{alphas['D_D']:.2f} (D–D) and R^{alphas['D_Li']:.2f} (D–Li) from 1 to 2 nm "
                f"(head-on collisions of the two explosions should give R^2).", ""]
    checks = []
    for ch in ("D_D", "D_Li"):
        thermal = tail_energy(edges, np.array(reps[1.0]["pair_exposure"][ch]), top=0.5)
        pred = tail_yield(reps[1.5], ch, lithium, 2.0 / 1.5, 10 * thermal) * 2  # tail_yield divides by 2 crystals
        sim = tail_yield(reps[2.0], ch, lithium, 1.0, 10 * thermal) * 2
        checks.append(f"{ch.replace('_', '+')}: predicted {pred:.2g}, simulated {sim:.2g}")
    out += ["Check, predicting the 2 nm pair from the 1.5 nm pair: " + "; ".join(checks) + ".", "",
            "| crystal radius (nm), gap 2.5 R | collision D+D per pair (extrapolated) | collision D+Li per pair |",
            "|---|---|---|"]
    radius_sim = reps[2.0]["radius_nm"]
    for target in TARGETS_NM:
        vals = []
        for ch in ("D_D", "D_Li"):
            thermal = tail_energy(edges, np.array(reps[1.0]["pair_exposure"][ch]), top=0.5)
            vals.append(tail_yield(reps[2.0], ch, lithium, target / radius_sim, 10 * thermal) * 2)
        out.append(f"| {target:g} | {vals[0]:.2g} | {vals[1]:.2g} |")
    out.append("")
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--dir", default="docs/results/experiments")
    args = ap.parse_args()
    root = Path(args.dir)
    out = ["# Fusion inside the crystal: size scaling", ""]
    for lithium in ("Li6", "Li7"):
        for profile in ("uniform", "expanded3"):
            reps = {}
            for r in SIZES:
                f = root / f"shock-{lithium}-{profile}-R{r:g}.json"
                if f.exists():
                    reps[r] = json.loads(f.read_text())
            if len(reps) < 2:
                continue
            out.append(f"## {lithium}, {profile}")
            out.append("")
            out.append("| R (nm) | ions | hardest D–Li pair (keV) | top-1% D–D exposure above (keV) "
                       "| top-1% D–Li exposure above (keV) | in-crystal D+D | in-crystal D+Li |")
            out.append("|---|---|---|---|---|---|---|")
            edges = None
            for r, rep in sorted(reps.items()):
                pe = rep["pair_exposure"]
                edges = np.array(pe["edges_keV"])
                out.append(f"| {r:g} | {rep['particles']} | {rep['dynamics']['max_D_Li_pair_energy_keV']:.3g} "
                           f"| {tail_energy(edges, np.array(pe['D_D'])):.3g} "
                           f"| {tail_energy(edges, np.array(pe['D_Li'])):.3g} "
                           f"| {rep['fusion']['in_crystal_D_D']:.2g} | {rep['fusion']['in_crystal_D_Li']:.2g} |")
            # Scaling exponent of the hard tail with R (between the smallest and largest size)
            rs = sorted(reps)
            lo, hi = reps[rs[0]], reps[rs[-1]]
            alphas = {}
            for ch in ("D_D", "D_Li"):
                e_lo = tail_energy(edges, np.array(lo["pair_exposure"][ch]))
                e_hi = tail_energy(edges, np.array(hi["pair_exposure"][ch]))
                alphas[ch] = math.log(e_hi / e_lo) / math.log(rs[-1] / rs[0]) if e_lo > 0 and e_hi > 0 else float("nan")
            out.append("")
            out.append(f"Tail energy grows as R^{alphas['D_D']:.2f} (D–D) and R^{alphas['D_Li']:.2f} (D–Li) "
                       f"from {rs[0]:g} to {rs[-1]:g} nm. Collisions that scale with the explosion give R^2; "
                       f"ordinary neighbours give R^0.")
            # Prediction check: 2 nm -> 3 nm, scaling the tail above 10x the thermal scale
            if 2.0 in reps and 3.0 in reps:
                base, big = reps[2.0], reps[3.0]
                rows = []
                passed = True
                for ch in ("D_D", "D_Li"):
                    thermal = tail_energy(edges, np.array(reps[rs[0]]["pair_exposure"][ch]), top=0.5)
                    pred = tail_yield(base, ch, lithium, 1.5, 10 * thermal)
                    actual_hist = tail_yield(big, ch, lithium, 1.0, 10 * thermal)
                    rows.append(f"{ch.replace('_', '+')}: predicted {pred:.2g}, simulated {actual_hist:.2g}")
                    if actual_hist > 0 and not 0.1 <= pred / actual_hist <= 10:
                        passed = False
                    if not alphas[ch] >= 1.0:
                        passed = False
                out.append("")
                out.append("Check, predicting the 3 nm crystal from the 2 nm one with the R² rule: "
                           + "; ".join(rows) + ".")
                out.append("")
                radius_sim = big["radius_nm"]
                if passed:
                    out.append("| target radius (nm) | in-crystal D+D (extrapolated) | in-crystal D+Li (extrapolated) |")
                    out.append("|---|---|---|")
                    for target in TARGETS_NM:
                        s = target / radius_sim
                        vals = []
                        for ch in ("D_D", "D_Li"):
                            thermal = tail_energy(edges, np.array(reps[rs[0]]["pair_exposure"][ch]), top=0.5)
                            vals.append(tail_yield(big, ch, lithium, s, 10 * thermal))
                        out.append(f"| {target:g} | {vals[0]:.2g} | {vals[1]:.2g} |")
                else:
                    out.append("The R² rule fails here: the hardest collisions do not grow with crystal size. "
                               "The in-crystal yield therefore grows only with the number of ions (R³):")
                    out.append("")
                    out.append("| target radius (nm) | in-crystal D+D (× R³) | in-crystal D+Li (× R³) |")
                    out.append("|---|---|---|")
                    for target in TARGETS_NM:
                        k = (target / radius_sim) ** 3
                        out.append(f"| {target:g} | {big['fusion']['in_crystal_D_D'] * k:.2g} "
                                   f"| {big['fusion']['in_crystal_D_Li'] * k:.2g} |")
            out.append("")
    out += analyse_pairs(root)
    text = "\n".join(out) + "\n"
    (root / "size_scaling.md").write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
