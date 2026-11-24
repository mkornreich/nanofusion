"""A practical dense Z-pinch experiment, designed with the time-dependent pinch model.

nanofusion/improved/zpinch.py follows a fibre (or gas column) driven by a
generator circuit: ohmic and compressional heating, X-rays with opacity,
radiative collapse above the Pease-Braginskii current, alpha heating, and
instability growth in Alfvén times. This script
* checks the model against known numbers (Pease-Braginskii current, Bennett
  temperature of the FuZE pinch, the 1980s fibre experiments);
* shows the three barriers the static fibre_pinch route left out;
* designs a deuterium-fibre experiment on Z, with a 1% tritium trace, whose
  outcomes separate the stability scenarios;
* and gives the D-T outlook on Z and on a next-generation 60 MA machine.

Writes docs/results/zpinch_experiment/summary.md and results.json.

    python3 tools/zpinch_experiment.py
"""

from __future__ import annotations

import json
import math
import multiprocessing as mp
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import simulator as S
from nanofusion.improved import zpinch as Z

OUT = Path("docs/results/zpinch_experiment")
SCENARIOS = {"breaks in ~1 Alfvén time (1980s fibres)": {},
             "stable while the current rises (NRL 1987)": {"rise_stabilised": True},
             "sheared flow, fully stable (untested at this density)": {"shear": 0.1}}
KEYS = ("yield_J", "neutrons", "trace_dt_neutrons", "stopped", "t_stop_ns", "min_radius_um", "peak_T_keV",
        "burn_up", "delivered_J", "bank_J", "drawn_J", "xrays_J", "peak_current_MA", "density_at_peak_g_cc",
        "hot_alfven_times", "line_density")
THERMAL, PULSED = 0.4, 0.15
BLANKET = S.ALT.blanket_energy_per_fusion_mev(0.0, uranium=False) / S.Q_DT_MEV


def sim(kw: dict) -> dict:
    warnings.simplefilter("ignore")
    r = Z.run(Z.Pinch(**kw))
    return {k: r[k] for k in KEYS} | {"inputs": kw}


def run_all(jobs: dict) -> dict:
    with ProcessPoolExecutor(8, mp_context=mp.get_context("fork")) as ex:
        return dict(zip(jobs, ex.map(sim, jobs.values())))


def compute() -> dict:
    jobs = {}
    # the 1980s-90s fibres (deuterium, cold start): NRL Poseidon, LANL HDZP-II, Imperial College MAGPIE
    fibres = {"NRL Poseidon (1987)": {"radius_um": 50, "current_ma": 0.64, "rise_ns": 130},
              "LANL HDZP-II": {"radius_um": 10, "current_ma": 0.7, "rise_ns": 100},
              "MAGPIE": {"radius_um": 50, "current_ma": 1.4, "rise_ns": 240}}
    for name, kw in fibres.items():
        for sc, extra in SCENARIOS.items():
            jobs[f"{name} | {sc}"] = dict(kw, fuel="DD", **extra)
    # the barriers on Z: D-T fibre, 240 um, 26 MA
    z = {"fuel": "DT", "radius_um": 240, "current_ma": 26}
    jobs["Z cold fibre | stable while the current rises"] = dict(z, rise_stabilised=True)
    for sc, extra in SCENARIOS.items():
        jobs[f"Z D-T preheated | {sc}"] = dict(z, start_kev=0.1, **extra)
    # the experiment: deuterium fibres on Z, preheated, with a 1% tritium trace (computed from the same history)
    for radius in (60, 120, 240):
        for sc, extra in SCENARIOS.items():
            jobs[f"exp D {radius} um | {sc}"] = dict(fuel="DD", radius_um=radius, current_ma=26, start_kev=0.1, **extra)
    jobs["exp D 120 um cold | stable while the current rises (NRL 1987)"] = {"fuel": "DD", "radius_um": 120,
                                                                                  "current_ma": 26, "rise_stabilised": True}
    # outlook: D-T on Z and on a 60 MA machine, sheared flow
    for cur in (26, 60):
        for radius in (120, 240, 480, 960, 1920):
            jobs[f"D-T {cur} MA {radius} um | sheared"] = {"fuel": "DT", "radius_um": radius, "current_ma": cur,
                                                               "start_kev": 0.1, "shear": 0.1}
    runs = run_all(jobs)
    # the static route, for comparison
    static = S.evaluate("fibre_pinch", {"current_ma": 26, "t_kev": 10, "alfven_times": 1, "squeeze": 6})
    return {"runs": runs, "static": {"fusion_J": static.fusion_J, "electric_gain": static.electric_gain,
                                     "drivers": static.drivers},
            "pease_braginskii_MA": {str(k): Z.pease_braginskii_ma(k) for k in (7, 10)}}


def fmt(x: float) -> str:
    return f"{x:.2g}" if x >= 1 else "< 1"


def energy(x: float) -> str:
    if x >= 1e6:
        return f"{x / 1e6:.3g} MJ"
    if x >= 1e3:
        return f"{x / 1e3:.3g} kJ"
    return f"{x:.2g} J"


def write(d: dict) -> list[str]:
    runs = d["runs"]
    pb = d["pease_braginskii_MA"]
    fuze_n = 1e23 * math.pi * 3e-3**2
    fuze_t = Z.bennett_kev(0.2e6, fuze_n)
    lines = ["# A practical dense Z-pinch experiment", "",
             ("Produced by `tools/zpinch_experiment.py` with the time-dependent pinch model "
              "`nanofusion/improved/zpinch.py`. The proposal is [docs/zpinch_experiment.md](../../zpinch_experiment.md)."),
             "", "## Checks", "",
             (f"* **Pease-Braginskii current** (X-rays beat ohmic heating): {pb['7']:.2f}–{pb['10']:.2f} MA for a "
             "Coulomb logarithm of 7–10, against Haines' 0.433 √lnΛ = 1.15–1.37 MA."),
             (f"* **Bennett temperature of FuZE (2018):** 200 kA, ~1e23 per m³, radius ~3 mm give {fuze_t:.1f} keV; "
             "measured 1–2 keV."),
             ("* **The 1980s–90s deuterium fibres**, cold start. Thermonuclear D-D neutrons from the model, against "
             "the measured neutrons, which were mostly beam-target from instabilities:"), "",
             ("| Fibre experiment | Breaks in ~1 Alfvén time | Stable while the current rises | Sheared flow | "
             "Measured (mostly beam-target) |"), "|---|---|---|---|---|"]
    measured = {"NRL Poseidon (1987)": "8e9 (1e12 claimed)", "LANL HDZP-II": "5e9", "MAGPIE": "1e8–5e9"}
    for name, seen in measured.items():
        cells = [fmt(runs[f"{name} | {sc}"]["neutrons"]) for sc in SCENARIOS]
        lines.append(f"| {name} | {' | '.join(cells)} | {seen} |")
    lines += ["",
              ("The model's thermonuclear yields are far below the measured ones in every scenario, which agrees "
              "with the finding that those neutrons came from beams made by the instabilities."), "",
              "## Three barriers the static model left out", "",
              ("A 240 µm D-T fibre, 1 cm long, on Z (26 MA into a fixed load in 100 ns; 15 nH of machine "
              "inductance):"), "",
              "| Case | Current reached | Narrowest | Hottest | Ended | Fusion | X-rays | From the bank |",
              "|---|---|---|---|---|---|---|---|"]
    rows = [("Cold fibre, stable while the current rises", "Z cold fibre | stable while the current rises")]
    rows += [(f"Preheated to 100 eV, {sc}", f"Z D-T preheated | {sc}") for sc in SCENARIOS]
    for label, key in rows:
        r = runs[key]
        lines.append(f"| {label} | {r['peak_current_MA']:.1f} MA | {r['min_radius_um']:.3g} µm | "
                     f"{r['peak_T_keV']:.2f} keV | {r['stopped']} at {r['t_stop_ns']:.0f} ns | {energy(r['yield_J'])} | "
                     f"{energy(r['xrays_J'])} | {energy(r['drawn_J'])} of {energy(r['bank_J'])} |")
    cold = runs["Z cold fibre | stable while the current rises"]
    sheared = runs[f"Z D-T preheated | {list(SCENARIOS)[2]}"]
    st = d["static"]
    lines += ["",
              (f"1. **A cold fibre can't take the current.** At ~10 eV its resistance is far above the generator's "
              f"impedance, so the current stalls at {cold['peak_current_MA']:.1f} MA. Real fibres grew a hot, "
              "thin corona that carried the current instead, and broke up. Preheating the fibre (100 eV costs ~5 kJ "
              "per cm, within Z-Beamlet's reach, as in MagLIF) makes the dense column conduct."),
              (f"2. **The collapse caps its own current.** Above the Pease-Braginskii current the column radiates "
              f"faster than it heats and collapses. The collapse's back-voltage holds Z to "
              f"{sheared['peak_current_MA']:.1f} MA, and most of the energy leaves as X-rays."),
              (f"3. **Stability decides everything, and even full stability isn't enough on Z.** Fully "
              f"flow-stabilised, the column makes {energy(sheared['yield_J'])} (burn-up {sheared['burn_up']:.1%}, "
              f"{sheared['peak_T_keV']:.1f} keV at most): no ignition."),
              (f"* The static route (26 MA, 10 keV, one Alfvén time, squeeze 6) assumed a Bennett pinch already hot: "
              f"{energy(st['fusion_J'])} of fusion for an electric gain of {st['electric_gain']:.2g}. The dynamic "
              "model shows the current never gets there."), "",
              "## The experiment: preheated deuterium fibres on Z, with a 1% tritium trace", "",
              ("Each row is one fibre radius; each cell is what the model predicts in one stability scenario: "
              "thermonuclear D-D neutrons / D-T neutrons from the tritium trace / narrowest radius."), "",
              f"| Fibre radius | {' | '.join(SCENARIOS)} |", "|---|---|---|---|"]
    for radius in (60, 120, 240):
        cells = []
        for sc in SCENARIOS:
            r = runs[f"exp D {radius} um | {sc}"]
            cells.append(f"{fmt(r['neutrons'])} / {fmt(r['trace_dt_neutrons'])} / {r['min_radius_um']:.2g} µm")
        lines.append(f"| {radius} µm | {' | '.join(cells)} |")
    c120 = runs["exp D 120 um cold | stable while the current rises (NRL 1987)"]
    lines += ["",
              (f"* Without preheat (120 µm, stable during the rise): {fmt(c120['neutrons'])} D-D neutrons, at "
              f"{c120['peak_current_MA']:.1f} MA."),
              ("* The scenarios differ by orders of magnitude in thermonuclear neutrons, and the D-T/D-D ratio "
              "measures the ion temperature. Beam-target neutrons, the 1980s background, are anisotropic and "
              "spectrally shifted; thermonuclear ones are not."), "",
              "## Outlook: D-T, fully flow-stabilised", "",
              ("| Machine | Fibre radius | Current reached | Fusion (1 cm) | Burn-up | Hottest | Bank | "
              "Electricity out ÷ in |"), "|---|---|---|---|---|---|---|---|"]
    best, best_at, z_best = 0.0, "", 0.0
    for cur in (26, 60):
        for radius in (120, 240, 480, 960, 1920):
            r = runs[f"D-T {cur} MA {radius} um | sheared"]
            gain = r["yield_J"] * BLANKET * THERMAL / (r["bank_J"] / PULSED)
            if gain > best:
                best, best_at = gain, f"{radius} µm at {cur} MA"
            if cur == 26:
                z_best = max(z_best, gain)
            lines.append(f"| {cur} MA | {radius} µm | {r['peak_current_MA']:.1f} MA | {energy(r['yield_J'])} | "
                         f"{r['burn_up']:.1%} | {r['peak_T_keV']:.1f} keV | {energy(r['bank_J'])} | {gain:.2g} |")
    lines += ["",
              ("* Each shot pays for the whole bank: once the column has burnt and blown apart, the energy left in it "
              "is wasted. With a lithium blanket, 40% of heat to electricity and 15% of electricity reaching the load."),
              f"* **On Z, even fully flow-stabilised, the best is {z_best:.2g} of breakeven.**",
              (f"* **On a 60 MA machine, wide fibres ignite** and burn up to a third of their fuel. The best, "
              f"{best_at}, returns {best:.2g}× the electricity: net energy, but only if the pinch stays stable for "
              "its whole burn, which no dense pinch has done."), ""]
    return lines


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / "results.json"
    if "--reuse" in sys.argv and cache.exists():
        d = json.loads(cache.read_text())
    else:
        d = compute()
        cache.write_text(json.dumps(d, indent=1, default=float))
    lines = write(d)
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
