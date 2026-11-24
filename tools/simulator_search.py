"""Search every energy route in one simulator for the best electricity out / electricity in.

All the experiments of docs/energy.md share one energy account in
nanofusion/improved/simulator.py. This script samples each route's parameters
(and the boosters: spin-polarised fuel, blankets), refines the best designs,
and ranks them within budgets of how far beyond today's best technology they
may go.

Writes docs/results/simulator/summary.md and results.json.

    python3 tools/simulator_search.py
"""

from __future__ import annotations

import json
import math
import sys
from dataclasses import replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import simulator as S  # noqa: E402

OUT = Path("docs/results/simulator")
SAMPLES = 300
REFINE = 100
SEED = 1
URANIUM_K = 0.95  # a subcritical uranium blanket at MYRRHA's design k_eff (section 13)


def fmt(x: float) -> str:
    if x == 0:
        return "0"
    if not math.isfinite(x):
        return "∞"
    if 0.01 <= abs(x) < 10000:
        return f"{float(f'{x:.2g}'):,g}"
    mantissa, exponent = f"{x:.1e}".split("e")
    mantissa = mantissa.rstrip("0").rstrip(".")
    sup = str(int(exponent)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"10{sup}" if mantissa == "1" else f"{mantissa}×10{sup}"


def params_text(o: S.Outcome) -> str:
    parts = []
    for q in S.ROUTES[o.route].params:
        v = o.params[q.name]
        if isinstance(v, bool):
            if v:
                parts.append(q.name.replace("_", " "))
        elif isinstance(v, float):
            parts.append(f"{q.name} {fmt(v)}")
        else:
            parts.append(f"{q.name} {v}")
    return ", ".join(parts)


def leaderboard(results: dict) -> list[str]:
    head = "| Route (section) | " + " | ".join(S.budget_label(b) for b in S.BUDGETS) + " |"
    lines = [head, "|---|" + "---|" * len(S.BUDGETS)]
    order = sorted(results, key=lambda n: -max((o.electric_gain for o in results[n].values() if o), default=0))
    for name in order:
        cells = []
        for b in S.BUDGETS:
            o = results[name][b]
            cells.append("–" if o is None else (f"**{fmt(o.electric_gain)}**" if o.electric_gain >= 1
                                                else fmt(o.electric_gain)))
        lines.append(f"| {name} ({S.ROUTES[name].section}) | " + " | ".join(cells) + " |")
    return lines


def design(o: S.Outcome) -> list[str]:
    lines = [f"**{o.route}** ({S.ROUTES[o.route].title}): {params_text(o)}.", "",
             f"* Fusion {fmt(o.fusion_J)} J from " + ", ".join(f"{fmt(v)} J of {k}" for k, v in o.drivers.items())
             + f"; electricity out / in **{fmt(o.electric_gain)}**. {o.notes[:1].upper() + o.notes[1:]}."]
    for r in o.all_requirements:
        best = f"{fmt(r.best)} ({S.TODAY[r.name][1]})" if r.name in S.TODAY else "never demonstrated"
        mark = " ← the hardest" if r.name == o.limiting else ""
        lines.append(f"* {r.name}: {fmt(r.needed)} against {best}, ×{fmt(r.stretch)}{mark}")
    return lines + [""]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    fusion_only = S.search(samples=SAMPLES, refine=REFINE, seed=SEED)
    uranium = S.search(samples=SAMPLES // 2, refine=REFINE // 2, seed=SEED,
                       settings=replace(S.Settings(), blanket="uranium", k_eff=URANIUM_K), budgets=(1.0, 10.0))
    s = S.Settings()
    lines = ["# Every route in one simulator: the search for net energy", "",
             "Produced by `tools/simulator_search.py`. Code: `nanofusion/improved/simulator.py` "
             "(`python -m nanofusion sim list | run | search`).", "",
             "## The account", "",
             "Each route of docs/energy.md is wrapped in the same energy account:", "",
             "electricity → driver (laser light, accelerator beam or pulsed current) → energy in the fuel → "
             "fusion → heat (with a lithium blanket that breeds tritium: +4.78 MeV per D-T fusion) → electricity.",
             "",
             f"* Driver efficiencies: laser {s.laser_efficiency:.0%} (fusion-laser designs; NIF's own is ~0.5%), "
             f"pulsed power {s.pulsed_power_efficiency:.0%}, accelerators {s.accelerator_efficiency:.0%} (a design "
             f"figure). Light reaching crystal ions or fuel heat: {s.laser_to_ions:.0%}; light reaching squeezed "
             f"fuel: {s.implosion_efficiency:.0%}. Heat to electricity: {s.thermal_efficiency:.0%}.",
             "* **Electric gain** = electricity out / electricity in. Net energy means more than 1.",
             "* **Stretch**: each design lists what it demands (laser energy and power, field, density, "
             "pressure, current, how long a pinch stays stable, muon supply, ...) against the best ever achieved. "
             "The largest ratio is its stretch. Spin-polarised fuel, never shown to survive in a plasma, counts as "
             "infinite stretch.",
             f"* The search samples each route {SAMPLES} times (log-uniform over each range, including the ends, "
             f"with and without the boosters), then hill-climbs {REFINE} steps from the best design within each "
             "budget.", "",
             "## Best electric gain of each route, by how far beyond today's technology it may go", ""]
    lines += leaderboard(fusion_only)
    lines += ["", "Bold: more than 1, net electricity on paper.", "", "## The winning designs", ""]
    winners = {}
    for b in S.BUDGETS:
        ranked = sorted((o for by in fusion_only.values() if (o := by[b]) is not None),
                        key=lambda o: -o.electric_gain)
        winners[b] = ranked[0]
        lines += [f"### {S.budget_label(b).capitalize()}", ""] + design(ranked[0])
    pinch = fusion_only["fibre_pinch"]
    best10 = pinch[10.0]
    best100 = pinch[100.0]
    today = winners[1.0]
    lines += ["## What the search found", "",
              f"* **With today's technology nothing comes close.** The best is {today.route.replace('_', ' ')} at "
              f"{fmt(today.electric_gain)}: a muon returns about a tenth of the electricity it cost.",
              f"* **One route crosses 1 on paper: the dense D-T Z-pinch** (section 19). Within ten times today's "
              f"best on every demand at once it reaches {fmt(best10.electric_gain)}. Allowed a hundred times, it "
              f"reaches {fmt(best100.electric_gain)}, limited by {best100.limiting}. That needs the pinch held "
              f"stable about {fmt(best100.params['alfven_times'])} Alfvén times at {fmt(best100.params['current_ma'])} "
              "MA. No dense, hot pinch has stayed stable past its peak current; Zap Energy's sheared-flow pinches "
              "last 1,000–10,000 growth times, but at a millionth of solid density. The 1980 Los Alamos dense "
              "Z-pinch reactor study also found net power on paper, assuming stability.",
              "* **The crystal routes never get near.** Exploding crystals stay below a percent of breakeven "
              "within 100 times today's technology, held back by picosecond laser energy. The tamped, magnetised "
              "ball (section 12), the best of them, needs a picosecond laser far larger than any.",
              "* **Unlimited technology** lets several routes pass 1 on paper: pycnonuclear burning at "
              "white-dwarf density, a magnetised ball with spin-polarised fuel, muons in squeezed fuel. Each needs "
              "something many orders of magnitude beyond anything built.",
              "* **Corrections the search forced.** Its first runs found a pinch \"burning\" more than all its fuel, "
              "a pinch at 130,000 Gbar, and one at 3.6 keV whose X-rays would have cost 17 times its drive. Fuel "
              "burn-up, a pressure demand and X-ray resupply are now in the model; those designs are gone.", "",
              f"## With a uranium blanket (k_eff = {URANIUM_K})", "",
              "A subcritical uranium blanket turns each D-T neutron into fission heat: about "
              f"{fmt(S.ALT.blanket_energy_per_fusion_mev(URANIUM_K) / S.Q_DT_MEV)} times the fusion energy. That "
              "is a fission reactor driven by a neutron source (section 13).", "",
              "| Route | today's technology | ≤ 10× today |", "|---|---|---|"]
    for name in sorted(uranium, key=lambda n: -max((o.electric_gain for o in uranium[n].values() if o), default=0)):
        by = uranium[name]
        cells = ["–" if by[b] is None else (f"**{fmt(by[b].electric_gain)}**" if by[b].electric_gain >= 1
                                             else fmt(by[b].electric_gain)) for b in (1.0, 10.0)]
        lines.append(f"| {name} | " + " | ".join(cells) + " |")
    lines += ["", "* **Fission makes several routes \"net energy\"**, muon catalysis with today's accelerators "
              "among them, but over 99% of the heat is fission. That is the conclusion of section 13 again.", ""]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")
    data = {"fusion_only": {n: {S.budget_label(b): (o.summary() if o else None) for b, o in by.items()}
                            for n, by in fusion_only.items()},
            "uranium": {n: {S.budget_label(b): (o.summary() if o else None) for b, o in by.items()}
                        for n, by in uranium.items()},
            "settings": {"samples": SAMPLES, "refine": REFINE, "seed": SEED, "uranium_k_eff": URANIUM_K}}
    (OUT / "results.json").write_text(json.dumps(data, indent=1, default=float))
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
