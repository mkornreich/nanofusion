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

from nanofusion.improved import simulator as S

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


# Ideas that can be switched on or off in a route; the search re-optimises everything else around each choice.
IDEAS = [
    ("compressed_icf", "ignitor", ("hot spot", "crystal spark")),
    ("compressed_icf", "driver", ("laser", "pulsed-power X-rays")),
    ("compressed_icf", "spin", (False, True)),
    ("tamped_ball", "ions_first", (False, True)),
    ("tamped_ball", "fuel", ("DT", "DD")),
    ("tamped_ball", "spin", (False, True)),
    ("fibre_pinch", "spin", (False, True)),
    ("magnetised_liner", "spin", (False, True)),
    ("fast_ignition", "shaped", (False, True)),
    ("cold_catcher", "source", ("crystals", "accelerator")),
    ("cold_catcher", "screening_ev", (0.0, 1500.0)),
    ("pB_catcher", "avalanche", (False, True)),
    ("pB_catcher", "screening_ev", (0.0, 1500.0)),
    ("muon_catalysis", "sticking", ("PSI 1987 (0.45%)", "theory (0.65%)")),
    # untested theory
    ("muon_catalysis", "rf_stripping", ("none", "Kulsrud/TEK")),
    ("opaque_pB", "alpha_channelling", (0.0, 0.5, 0.9)),
    ("opaque_pB", "direct_conversion", (False, True)),
    ("son_fisch_pB", "alpha_channelling", (0.0, 0.5, 0.9)),
    ("son_fisch_pB", "direct_conversion", (False, True)),
    ("dynamic_pinch", "load", ("fibre", "gas column", "liner")),
    ("dynamic_pinch", "stability", ("breaks", "rising current", "sheared flow")),
    ("dynamic_pinch", "ice_burn", ("ablation", "as an ICF shell")),
]


def gain(o) -> float:
    return o.electric_gain if o is not None else 0.0


def pinch_design(o) -> str:
    """One line on a dynamic_pinch design: its load, current and the convergence it needs."""
    p = o.params
    squeeze = next(r.needed for r in o.all_requirements if r.name == "radial squeeze of the fuel")
    if p["load"] == "liner":
        what = (f"a {fmt(p['liner_mm'])} mm liner with a wall 1/{fmt(p['liner_aspect'])} of its radius, "
                f"{fmt(p['ice_um'])} µm of D-T ice and {fmt(p['bz_tesla'])} T")
    elif p["load"] == "gas column":
        what = f"a {fmt(p['column_mm'])} mm gas column of {fmt(p['line_mg_cm'])} mg/cm"
    else:
        what = f"a {fmt(p['radius_um'])} µm fibre"
    return f"{what} at {fmt(p['current_ma'])} MA, converging {fmt(squeeze)}× (NIF's best: 35)"


def lead_findings(icf: dict, pinch: dict, dyn: dict, dyn_theory: dict, lead10) -> list[str]:
    """Bullets on who leads beyond today's technology: inertial fusion, and the pinch, static and in time."""
    leads_pinch = lead10 is not None and lead10.route == "dynamic_pinch"
    icf_line = (f"Within ten times today's best on every demand it reaches {fmt(icf[10.0].electric_gain)} (target gain "
                f"{fmt(icf[10.0].target_gain)}); within a hundred, {fmt(icf[100.0].electric_gain)}. The ten-times design "
                f"needs a drive of {fmt(icf[10.0].driver_J / 1e6)} MJ, a cold shell of {fmt(icf[10.0].params['rho_r'])} "
                f"g/cm² ({fmt(icf[10.0].params['rho_r'] / S.TODAY['fuel areal density (g/cm²)'][0])} times NIF's) at an "
                f"adiabat of {fmt(icf[10.0].params['adiabat'])}. That matches the standard view that inertial fusion "
                "power needs target gains of tens to a hundred and an efficient driver.")
    out = [("* **Ordinary inertial fusion is the most proven route beyond that.** " if leads_pinch else
            "* **Ordinary inertial fusion leads everywhere beyond that.** ") + icf_line,
           (f"* **The dense Z-pinch, static**: {fmt(pinch[10.0].electric_gain)} within ten times today, "
            f"{fmt(pinch[100.0].electric_gain)} within a hundred, if it stays stable for many Alfvén times.")]
    if dyn[10.0] is not None:
        out.append(
            f"* **Followed in time, the pinch {'leads within ten times today' if leads_pinch else 'does differently'}"
             f"** (dynamic_pinch: {fmt(gain(dyn[10.0]))}, against the static {fmt(pinch[10.0].electric_gain)}; "
             f"{fmt(gain(dyn[100.0]))} within a hundred). Its best design is {pinch_design(dyn[10.0])}. That is in "
             "effect cylindrical inertial fusion, in a 0-D model with no liner instability, so its convergence is "
             "the number to doubt. A cold fibre can't take the current; a gas column or a liner can. With untested "
             f"ideas it reaches {fmt(gain(dyn_theory[10.0]))} within ten times today.")
    return out


def idea_effects(known: dict, found_designs: dict | None = None) -> tuple[list[str], list]:
    """Best electric gain of a route with one idea held at each setting, everything else re-optimised.
    `known` (the main search's results) seeds each search with that route's best designs, and is updated with
    any better design these searches find."""
    lines = ["| Route | Idea | Setting | " + " | ".join(S.budget_label(b) for b in S.BUDGETS) + " |",
             "|---|---|---|" + "---|" * len(S.BUDGETS)]
    rows = []
    for route, name, values in IDEAS:
        names = {q.name for q in S.ROUTES[route].params}
        for value in values:
            # other untested ideas stay off; this one is let in when it is switched on
            fixed = {k: v for k, v in S.THEORY_OFF.items() if k in names and k != name}
            fixed[name] = value
            theory = name in S.THEORY_OFF and value != S.THEORY_OFF[name]
            start = [o.params for o in known[route].values() if o is not None]
            found = S.search([route], samples=SAMPLES, refine=REFINE, seed=SEED, fixed=fixed, start=start,
                             theory=theory)[route]
            if found_designs is not None:
                found_designs.setdefault(route, []).extend(o.params for o in found.values() if o is not None)
            cells = ["–" if found[b] is None else fmt(found[b].electric_gain) for b in S.BUDGETS]
            label = f"{name} (untested)" if name in S.THEORY_OFF else name
            lines.append(f"| {route} | {label} | {value} | " + " | ".join(cells) + " |")
            rows.append({"route": route, "idea": name, "value": value,
                         **{S.budget_label(b): (found[b].electric_gain if found[b] else None) for b in S.BUDGETS}})
            for b in S.BUDGETS:  # a better design of the route (without untested theory) found here: keep it
                if found[b] is not None and not found[b].untested and (known[route][b] is None
                                             or found[b].electric_gain > known[route][b].electric_gain):
                    known[route][b] = found[b]
    return lines, rows


def theory_findings(plain: dict, theory: dict, idea_rows: list) -> list[str]:
    """Bullets on what the untested theoretical ideas change."""
    def best(results, budget, routes=None):
        o = best_of(results, budget, routes)
        return (o.electric_gain if o else 0.0), o

    out = []
    g10, o10 = best(theory, 10.0)
    g100, o100 = best(theory, 100.0)
    p10, l10 = best(plain, 10.0)
    if l10 is not None and o10.route != l10.route:
        opener = (f"* **Untested theory changes the leader within ten times today.** With it, the best is "
                  f"{o10.route.replace('_', ' ')} at {fmt(g10)}"
                  + (f" (with {'; '.join(o10.untested)})" if o10.untested else "")
                  + f", against {l10.route.replace('_', ' ')} at {fmt(p10)} without")
    else:
        opener = (f"* **Untested theory moves the numbers, not the verdict.** With it, the best within ten times "
                  f"today is {o10.route.replace('_', ' ')} at {fmt(g10)}"
                  + (f" (with {'; '.join(o10.untested)})" if o10.untested else ""))
    out.append(opener + f"; within a hundred, {o100.route.replace('_', ' ')} at {fmt(g100)}"
               + (f" (with {'; '.join(o100.untested)})" if o100.untested else "") + ".")
    pb = {"son_fisch_pB", "opaque_pB", "pB_catcher"}
    gp, _ = best(theory, S.INF, pb)
    gp0, _ = best(plain, S.INF, pb)
    out.append(f"* **Proton–boron with alpha channelling and direct conversion** reaches {fmt(gp)} with no limits on "
               f"technology, against {fmt(gp0)} without them. Channelling feeds the ions, but the X-rays still win.")
    ga, oa = best(theory, S.INF, {"antimatter_spark"})
    if oa:
        n = next(r.needed for r in oa.all_requirements if r.name == "antiprotons per shot")
        factory = any("factory" in u for u in oa.untested)
        out.append(f"* **Antimatter sparks** reach at most {fmt(ga)}, with {fmt(n)} antiprotons per shot, "
                   f"{fmt(n / S.TODAY['antiprotons per shot'][0])} times the most ever stored"
                   + (", from a factory 10,000 times more efficient than CERN's" if factory else "")
                   + ". Making antiprotons costs more than the fusion they spark.")
    gl, _ = best(theory, S.INF, {"laser_assisted"})
    out.append(f"* **Light-assisted tunnelling** reaches at most {fmt(gl)}, with no limits on technology: holding "
               "5×10²¹ W/cm² over the fuel costs more light than the extra fusion returns, even with the most "
               "optimistic tunnelling boost.")
    rf = {r["value"]: r[S.budget_label(1.0)] for r in idea_rows if r["idea"] == "rf_stripping"}
    if rf:
        out.append(f"* **Freeing stuck muons with fields** takes muon catalysis from {fmt(rf.get('none') or 0)} to "
                   f"{fmt(rf.get('Kulsrud/TEK') or 0)} with today's technology: twice the fusions per muon, but "
                   "keeping the fuel below 25 K costs ~30 MeV of electricity per 17.6 MeV fusion.")
    return out


def best_of(results: dict, budget: float, routes=None) -> S.Outcome | None:
    found = [o for name, by in results.items() if (routes is None or name in routes) and (o := by[budget]) is not None]
    return max(found, key=lambda o: o.electric_gain, default=None)


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    fusion_only = S.search(samples=SAMPLES, refine=REFINE, seed=SEED)
    uranium = S.search(samples=SAMPLES // 2, refine=REFINE // 2, seed=SEED,
                       settings=replace(S.Settings(), blanket="uranium", k_eff=URANIUM_K), budgets=(1.0, 10.0))
    found_designs: dict = {}
    idea_lines, idea_rows = idea_effects(fusion_only, found_designs)  # also improves fusion_only where it can
    # The same search with untested theory let in, starting from the best designs found so far.
    theory = {name: S.search([name], samples=SAMPLES, refine=REFINE, seed=SEED, theory=True,
                             start=[o.params for o in fusion_only[name].values() if o]
                             + found_designs.get(name, []))[name] for name in S.ROUTES}
    s = S.Settings()
    lines = ["# Every route in one simulator: the search for net energy", "",
             ("Produced by `tools/simulator_search.py`. Code: `nanofusion/improved/simulator.py` "
             "(`python -m nanofusion sim list | run | search`)."), "",
             "## The account", "",
             f"Each of the {len(S.ROUTES)} routes is wrapped in the same energy account:", "",
             ("electricity → driver (laser light, accelerator beam or pulsed current) → energy in the fuel → "
             "fusion → heat (with a lithium blanket that breeds tritium: +4.78 MeV per D-T fusion) → electricity."),
             "",
             (f"* Driver efficiencies: laser {s.laser_efficiency:.0%} (fusion-laser designs; NIF's own is ~0.5%), "
             f"pulsed power {s.pulsed_power_efficiency:.0%}, accelerators {s.accelerator_efficiency:.0%} (a design "
             f"figure). Light reaching crystal ions or fuel heat: {s.laser_to_ions:.0%}; implosions put "
             f"{s.implosion_efficiency:.0%} of their drive into the fuel. Heat to electricity: "
             f"{s.thermal_efficiency:.0%}."),
             "* **Electric gain** = electricity out / electricity in. Net energy means more than 1.",
             ("* **Stretch**: each design lists what it demands (laser energy and power, field, density, "
             "pressure, current, how long a pinch stays stable, muon supply, ...) against the best ever achieved. "
             "The largest ratio is its stretch. Ideas never shown to work (spin-polarised fuel surviving in a "
             "plasma, a microbubble implosion) count as infinite stretch."),
             (f"* The search samples each route {SAMPLES} times (log-uniform over each range, including the ends, "
             f"with every option on and off), then hill-climbs {REFINE} steps from the best design within each "
             "budget."),
             ("* **The reference route**, `compressed_icf`, is ordinary inertial fusion: an isobaric hot spot "
             "(10 keV, ρR 0.3 g/cm²) inside a dense shell, 1% of the drive in the fuel, burn-up ρR/(ρR + 7 g/cm²). "
             "With NIF-like inputs it gives a target gain of ~2 against NIF's record 4.1, so it is if anything "
             "conservative."), "",
             "## Best electric gain of each route, by how far beyond today's technology it may go", ""]
    lines += leaderboard(fusion_only)
    lines += ["", "Bold: more than 1, net electricity on paper.", "", "## The winning designs", ""]
    for b in S.BUDGETS:
        lines += [f"### {S.budget_label(b).capitalize()}", ""] + design(best_of(fusion_only, b))
    crystal_routes = {"cold_catcher", "crystal_ball", "tamped_ball", "fast_ignition", "pauli_beam", "projectiles",
                      "channelling", "muon_factory", "hollow_shell"}
    lines += ["### The best route that uses exploding crystals, within 100× today", ""]
    lines += design(best_of(fusion_only, 100.0, crystal_routes))
    icf = fusion_only["compressed_icf"]
    pinch = fusion_only["fibre_pinch"]
    dyn = fusion_only["dynamic_pinch"]
    dyn_theory = theory["dynamic_pinch"]
    today = best_of(fusion_only, 1.0)
    crystal100 = best_of(fusion_only, 100.0, crystal_routes)
    spark = next(r for r in idea_rows if r["idea"] == "ignitor" and r["value"] == "crystal spark")
    hot = next(r for r in idea_rows if r["idea"] == "ignitor" and r["value"] == "hot spot")
    drive = {r["value"]: r[S.budget_label(10.0)] for r in idea_rows if r["idea"] == "driver"}
    lines += ["## If the untested theories work", "",
              ("The same search with the untested ideas let in: spin-polarised fuel, alpha channelling, direct "
              "conversion of fusion products, freeing stuck muons with fields, antimatter sparks, light-assisted "
              "tunnelling, the microbubble implosion, and ways to keep a dense pinch stable. Stretch still counts the technology each design needs."), ""]
    lines += leaderboard(theory)
    lines += ["", "Ideas the best designs use:", ""]
    for name in sorted(theory, key=lambda n: -max((o.electric_gain for o in theory[n].values() if o), default=0)):
        used = sorted({u for o in theory[name].values() if o for u in o.untested})
        if used:
            lines.append(f"* {name}: {'; '.join(used)}")
    lines += [""]
    lines += ["## What each idea adds", "",
              ("Each row holds one idea at one setting and lets the search re-optimise everything else. Untested "
              "ideas are let into the search only in their \"on\" rows."), ""]
    lines += idea_lines
    lines += ["", "## What the search found", "",
              (f"* **With today's technology nothing reaches breakeven.** The best is {today.route.replace('_', ' ')} "
              f"at {fmt(today.electric_gain)}."),
              *lead_findings(icf, pinch, dyn, dyn_theory, best_of(fusion_only, 10.0)),
              (f"* **Exploding crystals never help.** The best crystal route within a hundred times today is "
              f"{crystal100.route.replace('_', ' ')} at {fmt(crystal100.electric_gain)}. As the spark for inertial "
              f"fusion, crystals do worse than an ordinary hot spot at every budget (unlimited: "
              f"{fmt(spark['unlimited'] or 0)} against {fmt(hot['unlimited'] or 0)}): the isotropic explosion wastes "
              "the light. As a muon factory they make muons, but not cheaply."),
              (f"* **Lasers beat pulsed-power X-rays** as the driver of a compressed capsule (best within ten times "
              f"today: {fmt(drive['laser'] or 0)} against {fmt(drive['pulsed-power X-rays'] or 0)}). Pulsed power "
              "makes its drive more efficiently, but on Z only ~0.2% of the stored energy reached a capsule."),
              ("* **Small ideas stay small.** Electron screening and the proton–boron avalanche change nothing you "
              "can see; ions heated first help the tamped ball two to five times; deuterium-only fuel costs a factor of "
              "tens; spin-polarised fuel helps but has never been shown to survive in a plasma."),
              *theory_findings(fusion_only, theory, idea_rows),
              ("* **Corrections the search forced.** It found a pinch \"burning\" more than all its fuel, a pinch "
              "at 130,000 Gbar, a pinch too cool for its alphas to make up its X-ray losses, and a compressed "
              "capsule whose hot spot was 600 times the pressure of the shell around it. Burn-up, a pressure "
              "demand, X-ray resupply and an isobaric hot spot are now in the models; those designs are gone."), "",
              f"## With a uranium blanket (k_eff = {URANIUM_K})", "",
              ("A subcritical uranium blanket turns each D-T neutron into fission heat: about "
              f"{fmt(S.ALT.blanket_energy_per_fusion_mev(URANIUM_K) / S.Q_DT_MEV)} times the fusion energy. That "
              "is a fission reactor driven by a neutron source (section 13)."), "",
              "| Route | today's technology | ≤ 10× today |", "|---|---|---|"]
    for name in sorted(uranium, key=lambda n: -max((o.electric_gain for o in uranium[n].values() if o), default=0)):
        by = uranium[name]
        cells = ["–" if by[b] is None else (f"**{fmt(by[b].electric_gain)}**" if by[b].electric_gain >= 1
                                             else fmt(by[b].electric_gain)) for b in (1.0, 10.0)]
        lines.append(f"| {name} | " + " | ".join(cells) + " |")
    lines += ["", ("* **Fission makes several routes \"net energy\"**, muon catalysis with today's accelerators "
              "among them, but over 99% of the heat is fission. That is the conclusion of section 13 again."), ""]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")
    data = {"fusion_only": {n: {S.budget_label(b): (o.summary() if o else None) for b, o in by.items()}
                            for n, by in fusion_only.items()},
            "uranium": {n: {S.budget_label(b): (o.summary() if o else None) for b, o in by.items()}
                        for n, by in uranium.items()},
            "ideas": idea_rows,
            "theory": {n: {S.budget_label(b): (o.summary() if o else None) for b, o in by.items()}
                       for n, by in theory.items()},
            "settings": {"samples": SAMPLES, "refine": REFINE, "seed": SEED, "uranium_k_eff": URANIUM_K}}
    (OUT / "results.json").write_text(json.dumps(data, indent=1, default=float))
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
