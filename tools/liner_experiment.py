"""An experiment for the simulator's best liner design: a thin beryllium liner lined with D-T ice, at 34 MA.

The simulator's search (tools/simulator_search.py) found that a time-dependent pinch with a liner
(nanofusion/improved/zpinch.py) leads every route within ten times today's technology:
* a 4.5 mm beryllium liner, its wall 1/12 of its radius;
* 583 µm of D-T ice on its inner wall, 0.3 mg/cm3 of D-T gas, a 4.3 T field, almost no preheat;
* 34 MA.
This script asks what an experiment must show before anyone believes it:
* how the design depends on convergence (the liner breaking up early), current, wall, ice, field and gas;
* what a deuterium version fielded on Z now would show for each outcome;
* the D-T step, and the energy.

Writes docs/results/liner_experiment/summary.md and results.json (~1 minute on 8 cores).

    python3 tools/liner_experiment.py
"""

from __future__ import annotations

import json
import multiprocessing as mp
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import simulator as S
from nanofusion.improved import zpinch as Z

OUT = Path("docs/results/liner_experiment")
DESIGN = {"fuel": "DT", "current_ma": 34.06, "liner_um": 4495.0, "liner_aspect": 12.0, "fill_mg_cc": 0.3,
          "bz_tesla": 4.283, "ice_um": 583.3, "start_kev": 0.017, "hold_ns": 150.0}
KEYS = ("yield_J", "neutrons", "neutrons_dd", "neutrons_dt", "stopped", "t_stop_ns", "min_radius_um", "peak_T_keV",
        "T_at_peak_keV", "burn_up", "bank_J", "burn_width_ns", "bang_time_ns", "ice_left", "peak_current_MA",
        "line_density")
BLANKET = S.ALT.blanket_energy_per_fusion_mev(0.0, uranium=False) / S.Q_DT_MEV


def sim(kw: dict) -> dict:
    warnings.simplefilter("ignore")
    p = Z.Pinch(**kw)
    r = Z.run(p)
    out = {k: r[k] for k in KEYS}
    out["convergence"] = p.fuel_radius_um / r["min_radius_um"]
    out["inputs"] = kw
    return out


def gain(r: dict, pulsed: float = 0.15) -> float:
    heat = r["yield_J"] + r["neutrons_dt"] * (BLANKET - 1) * S.Q_DT_MEV * S.MEV
    preheat = 3 * r["line_density"] * r["inputs"]["start_kev"] * 1.602176634e-16 * 0.01 / 0.5
    return heat * 0.4 / (r["bank_J"] / pulsed + preheat / 0.1)


def compute() -> dict:
    jobs = {"design": dict(DESIGN)}
    for mc in (30, 40, 60, 100, 130, 150, 170):
        jobs[f"convergence|{mc}"] = dict(DESIGN, max_convergence=mc)
    for cur in (20, 26, 30, 34.06, 40, 50):
        jobs[f"current|{cur}"] = dict(DESIGN, current_ma=cur)
    for ar in (6, 9, 12):
        jobs[f"aspect|{ar}"] = dict(DESIGN, liner_aspect=ar)
    for ice in (0, 200, 400, 583.3):
        jobs[f"ice|{ice}"] = dict(DESIGN, ice_um=ice)
    for bz in (0.5, 4.283, 10, 20, 30):
        jobs[f"field|{bz}"] = dict(DESIGN, bz_tesla=bz)
    for fill in (0.1, 0.3, 1.0, 3.0):
        jobs[f"gas|{fill}"] = dict(DESIGN, fill_mg_cc=fill)
    # Stage 1: the same target in deuterium on Z, at the current Z can drive into it
    for cur in (20, 26):
        for ar in (6, 9, 12):
            for mc in (0, 40, 100):
                jobs[f"z|{cur}|{ar}|{mc}"] = dict(DESIGN, fuel="DD", current_ma=cur, liner_aspect=ar,
                                                  max_convergence=mc)
    names = list(jobs)
    with ProcessPoolExecutor(8, mp_context=mp.get_context("fork")) as ex:
        res = dict(zip(names, ex.map(sim, [jobs[n] for n in names])))
    return res


def fmt(x: float) -> str:
    if x == 0:
        return "0"
    if 1000 <= x < 1e6:
        return f"{float(f'{x:.2g}'):,.0f}"
    if x >= 1e6 or x < 0.01:
        return f"{x:.1e}"
    return f"{x:.3g}"


def energy(x: float) -> str:
    for unit, scale in (("GJ", 1e9), ("MJ", 1e6), ("kJ", 1e3)):
        if x >= scale:
            return f"{x / scale:.3g} {unit}"
    return f"{x:.2g} J"


def write(d: dict) -> list[str]:
    des = d["design"]
    p = Z.Pinch(**DESIGN)
    lines = ["# An experiment for the thin, ice-lined liner", "",
             ("Produced by `tools/liner_experiment.py` with `nanofusion/improved/zpinch.py`. The proposal is "
              "[docs/liner_experiment.md](../../liner_experiment.md)."), "",
             "## The design, as the simulator found it", "",
             "| | |", "|---|---|",
             (f"| Liner | beryllium, outer radius {p.liner_um / 1e3:.2f} mm, wall {p.liner_um / p.liner_aspect:.0f} µm "
             f"(1/{p.liner_aspect:.0f} of the radius), 1 cm long |"),
             f"| Fuel | {p.ice_um:.0f} µm of 50/50 D-T ice on the inner wall; {p.fill_mg_cc} mg/cm³ of D-T gas inside |",
             f"| Field, preheat | {p.bz_tesla:.1f} T axial; {p.start_kev * 1e3:.0f} eV (almost none) |",
             f"| Driver | {p.current_ma:.0f} MA in 100 ns into a fixed load; bank {energy(des['bank_J'])} |",
             (f"| Fusion | {energy(des['yield_J'])} ({des['burn_up']:.0%} of the fuel burnt), hottest "
             f"{des['peak_T_keV']:.0f} keV |"),
             (f"| Convergence | {des['convergence']:.0f}×: the fuel from {p.fuel_radius_um / 1e3:.2f} mm to "
             f"{des['min_radius_um']:.0f} µm |"),
             f"| Electricity out ÷ in | {fmt(gain(des))} (15% of electricity to the load, 40% heat to electricity) |",
             "", "## Its weak point: convergence", "",
             ("If the liner breaks up (magneto-Rayleigh-Taylor) before the fuel converges this far, compression "
              "stops there:"), "",
             "| Liner breaks up at | Fusion | Hottest | Electricity out ÷ in |", "|---|---|---|---|"]
    for mc in (30, 40, 60, 100, 130, 150, 170):
        r = d[f"convergence|{mc}"]
        lines.append(f"| {mc}× | {energy(r['yield_J'])} | {r['peak_T_keV']:.1f} keV | {fmt(gain(r))} |")
    lines.append(f"| never ({des['convergence']:.0f}×) | {energy(des['yield_J'])} | {des['peak_T_keV']:.0f} keV | "
                 f"{fmt(gain(des))} |")
    lines += ["",
              ("* **All or nothing.** The fuel only gets hot past ~100× and ignites near the full ~170×. At "
               "MagLIF's demonstrated 30–40×, the design makes tens of joules."), "",
              "## Other sensitivities, one at a time", "",
              "| Change | Fusion | Convergence | Electricity out ÷ in |", "|---|---|---|---|"]
    rows = [("current", "MA", (20, 26, 30, 34.06, 40, 50)), ("aspect", "wall ratio", (6, 9, 12)),
            ("ice", "µm of ice", (0, 200, 400, 583.3)), ("field", "T", (0.5, 4.283, 10, 20, 30)),
            ("gas", "mg/cm³ of gas", (0.1, 0.3, 1.0, 3.0))]
    for key, unit, values in rows:
        for v in values:
            r = d[f"{key}|{v}"]
            mark = " (design)" if r["inputs"] == DESIGN else ""
            lines.append(f"| {v:g} {unit}{mark} | {energy(r['yield_J'])} | {r['convergence']:.0f}× | "
                         f"{fmt(gain(r))} |")
    drop = {name: des["yield_J"] / max(d[key]["yield_J"], 1e-30)
            for name, key in (("30 MA", "current|30"), ("a wall ratio of 9", "aspect|9"),
                              ("1 mg/cm³ of gas", "gas|1.0"), ("20 T", "field|20"))}
    lines += ["",
              ("* **A knife-edge.** Small changes kill it: " + ", ".join(f"{k} cuts the fusion {fmt(v)}-fold"
                                                                      for k, v in drop.items())
               + ". An optimiser's narrow window like this is more often a model artefact than robust physics."),
              ("* **Below ~34 MA the fuel over-converges cold** (hundreds of times) without igniting: the model "
               "has nothing to stop it, a real liner would break up first."), "",
              "## Stage 1 on Z now: the same target in deuterium", "",
              ("Deuterium ice and gas, the same liner, at the current Z can drive. D-D fuel doesn't ignite, so this "
               "tests the implosion, not the burn. Each cell: D-D neutrons / stagnation radius / convergence."), "",
              "| Current, wall ratio | Liner intact | Breaks up at 40× | Breaks up at 100× |", "|---|---|---|---|"]
    for cur in (20, 26):
        for ar in (6, 9, 12):
            cells = []
            for mc in (0, 40, 100):
                r = d[f"z|{cur}|{ar}|{mc}"]
                cells.append(f"{fmt(r['neutrons_dd'])} / {r['min_radius_um']:.0f} µm / {r['convergence']:.0f}×")
            lines.append(f"| {cur} MA, 1/{ar} | {' | '.join(cells)} |")
    intact = [d[f"z|{c}|{a}|0"] for c in (20, 26) for a in (6, 9, 12)]
    broke = [d[f"z|{c}|{a}|40"] for c in (20, 26) for a in (6, 9, 12)]
    z12 = d["z|26|12|0"], d["z|26|12|40"]
    lines += ["",
              (f"* **The stagnation radius is the measurement.** With the liner intact the model squeezes the fuel to "
               f"{min(r['min_radius_um'] for r in intact):.0f}–{max(r['min_radius_um'] for r in intact):.0f} µm "
               f"({min(r['convergence'] for r in intact):.0f}–{max(r['convergence'] for r in intact):.0f}×, more "
               "than any implosion has reached: an upper bound). Breaking up at MagLIF's ~40×, it stops at "
               f"{min(r['min_radius_um'] for r in broke):.0f}–{max(r['min_radius_um'] for r in broke):.0f} µm. "
               "X-ray imaging of the stagnated column resolves the difference."),
              (f"* **The neutrons confirm it.** At 26 MA with the design's thin wall: {fmt(z12[0]['neutrons_dd'])} "
               f"D-D neutrons if it holds, {fmt(z12[1]['neutrons_dd'])} if it breaks at 40×. Secondary D-T neutrons "
               "(from D-D tritons) give the fuel's areal density."), ""]
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
