"""Making the dense Z-pinch testable, and finding what it takes to make energy.

The fibre pinch (tools/zpinch_experiment.py) stalls: a dense, cold column is too
resistive to take the generator's current, and collapses radiatively. This
script, with the time-dependent model nanofusion/improved/zpinch.py:
* tests levers that don't help (a rising current, a longer pinch, recovering
  the circuit's energy);
* calibrates a different load, a pre-formed gas column (as from a gas puff),
  against Z's deuterium gas-puff record;
* designs a deuterium gas-column experiment Z can do now;
* maps the electric gain over current, stability and tritium share, with
  each design's fuel load and radius optimised.

Writes docs/results/zpinch_experiment/energy.md and energy.json (~3 minutes on 8 cores).

    python3 tools/zpinch_energy.py
"""

from __future__ import annotations

import itertools
import json
import multiprocessing as mp
import sys
import warnings
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import simulator as S
from nanofusion.improved import zpinch as Z

OUT = Path("docs/results/zpinch_experiment")
AMU = 1.66053907e-27
BLANKET = S.ALT.blanket_energy_per_fusion_mev(0.0, uranium=False) / S.Q_DT_MEV
THERMAL = 0.4
KEYS = ("yield_J", "neutrons", "neutrons_dt", "neutrons_dd", "stopped", "t_stop_ns", "min_radius_um", "peak_T_keV",
        "burn_up", "bank_J", "circuit_end_J", "peak_current_MA", "hot_alfven_times", "line_density",
        "ice_burn_bound_J", "burn_width_ns", "T_at_peak_keV")
STABILITY = {"1 Alfvén time": {"growth_limit": 1}, "3": {"growth_limit": 3}, "10": {"growth_limit": 10},
             "30": {"growth_limit": 30}, "sheared flow": {"shear": 0.1}}
GAS_PUFF_RECORD = 3.9e13  # D-D neutrons, Z, ~15-17 MA (Coverdale et al.)


def column(mg_cm: float, radius_mm: float, tritium: float) -> dict:
    ion = (2.014 + 1.002 * tritium) * AMU
    return {"ions_per_m": mg_cm * 1e-4 / ion, "radius_um": radius_mm * 1e3, "wall_mm": max(10.0, 2 * radius_mm),
            "tritium": tritium, "start_kev": 0.005}


def sim(kw: dict) -> dict:
    warnings.simplefilter("ignore")
    r = Z.run(Z.Pinch(**kw))
    return {k: r[k] for k in KEYS} | {"inputs": kw}


def gain(r: dict, pulsed: float = 0.15, ice: bool = False) -> float:
    """Electricity out / in: D-T neutrons into a lithium blanket; each shot pays for the whole bank (and the
    preheat laser at 10%). ice=True adds the upper bound for a liner's frozen fuel."""
    fusion = r["yield_J"] + (r.get("ice_burn_bound_J", 0.0) if ice else 0.0)
    dt = r["neutrons_dt"] + ((r.get("ice_burn_bound_J", 0.0) / (S.Q_DT_MEV * S.MEV)) if ice else 0.0)
    heat = fusion + dt * (BLANKET - 1) * S.Q_DT_MEV * S.MEV
    return heat * THERMAL / (r["bank_J"] / pulsed + r.get("preheat_J", 0.0) / 0.1)


def run_all(jobs: list[dict]) -> list[dict]:
    with ProcessPoolExecutor(8, mp_context=mp.get_context("fork")) as ex:
        return list(ex.map(sim, jobs))


def compute() -> dict:
    out = {}
    # 1. levers that don't help, on the fibre
    lever_jobs = {}
    for rise in (100, 300, 1000):
        lever_jobs[f"rising current, {rise} ns rise"] = {"radius_um": 960, "current_ma": 60, "start_kev": 0.1,
                                                             "rise_stabilised": True, "rise_ns": rise, "hold_ns": rise}
    for length in (10, 30, 100):
        lever_jobs[f"{length} mm long, sheared flow"] = {"radius_um": 960, "current_ma": 60, "start_kev": 0.1, "shear": 0.1,
                                                             "length_mm": length}
    # 2. calibration: Z's deuterium gas puffs (~17 MA)
    cal_jobs = {}
    for mg, rad, lim in itertools.product((0.25, 0.5, 1.0), (10, 20), (1, 2, 3)):
        cal_jobs[f"{mg} mg/cm, {rad} mm, {lim} Alfvén"] = dict(column(mg, rad, 0.0), current_ma=17, length_mm=20,
                                                                growth_limit=lim, hold_ns=150)
    # 3. the experiment on Z now: deuterium gas columns with a 1% tritium trace, 26 MA
    exp_jobs = {}
    for mg, (name, extra) in itertools.product((0.5, 2.0, 7.0), STABILITY.items()):
        exp_jobs[f"{mg} mg/cm | {name}"] = dict(column(mg, 10, 0.01), current_ma=26, length_mm=20, **extra)
    # 4. energy map: D-T gas column and fibre, optimised over fuel load and radius
    map_jobs = []
    for cur, (name, extra), mg, rad in itertools.product((26, 45, 60, 80), STABILITY.items(), (2, 7, 20),
                                                         (3, 10, 20)):
        map_jobs.append((cur, name, "gas column", dict(column(mg, rad, 0.5), current_ma=cur, **extra)))
    for cur, (name, extra), rad in itertools.product((26, 45, 60, 80), STABILITY.items(), (240, 960, 1500)):
        map_jobs.append((cur, name, "fibre", dict(radius_um=rad, current_ma=cur, start_kev=0.1, **extra)))
    # 5. tritium share, at 60 and 80 MA with 10 Alfvén times
    trit_jobs = []
    for cur, f_t, mg, rad in itertools.product((60, 80), (0.5, 0.2, 0.1, 0.05, 0.0), (2, 7, 20), (3, 10, 20)):
        trit_jobs.append((cur, f_t, dict(column(mg, rad, f_t), current_ma=cur, growth_limit=10)))
    flat = list(lever_jobs.values()) + list(cal_jobs.values()) + list(exp_jobs.values()) + \
        [j for *_, j in map_jobs] + [j for *_, j in trit_jobs]
    res = run_all(flat)
    k = 0
    for group, jobs in (("levers", lever_jobs), ("calibration", cal_jobs), ("experiment", exp_jobs)):
        out[group] = {n: res[k + i] for i, n in enumerate(jobs)}
        k += len(jobs)
    best = {}
    for cur, name, load, _ in map_jobs:
        r = res[k]
        k += 1
        key = f"{load}|{cur}|{name}"
        if key not in best or gain(r) > gain(best[key]):
            best[key] = r
    out["map"] = best
    tbest = {}
    for cur, f_t, _ in trit_jobs:
        r = res[k]
        k += 1
        key = f"{cur}|{f_t}"
        if key not in tbest or gain(r) > gain(tbest[key]):
            tbest[key] = r
    out["tritium"] = tbest
    out["liner"] = liner()
    return out


def preheat_kev(energy_j: float, fill_mg_cc: float, radius_um: float, length_mm: float, ion_u: float) -> float:
    """Uniform fuel temperature (keV) from the preheat energy coupled to it."""
    ions = fill_mg_cc * 3.14159265 * (radius_um * 1e-6) ** 2 / (ion_u * AMU) * length_mm * 1e-3
    return energy_j / (3 * ions) / 1.602176634e-16


# Z MagLIF shots (Gomez et al. 2014, 2025; Weis et al. 2026; Lewis et al. 2023): inputs, measured D-D yield, ion
# temperature, burn width, convergence, and the clean 2-D simulated yield.
MAGLIF_SHOTS = {
    "2014 first shots, 19 MA": ({"current_ma": 19, "liner_um": 2790, "liner_aspect": 6, "length_mm": 7.5, "fill_mg_cc": 0.7,
                                     "bz_tesla": 10, "start_kev": preheat_kev(600, 0.7, 2325, 7.5, 2.014)},
                                2.0e12, 2.5, 2.0, 46, 6.1e13),
    "z3289, 20 MA": ({"current_ma": 20, "liner_um": 2790, "liner_aspect": 6, "length_mm": 10, "fill_mg_cc": 1.03, "bz_tesla": 15.9,
                          "start_kev": preheat_kev(1150, 1.03, 2325, 10, 2.014)}, 1.1e13, 3.1, 2.0, 31, 5.2e13),
    "record z3926, 18.5 MA": ({"current_ma": 18.5, "liner_um": 2616, "liner_aspect": 9, "length_mm": 10, "fill_mg_cc": 1.05,
                                   "bz_tesla": 15, "start_kev": preheat_kev(2150, 1.05, 2325, 10, 2.014)},
                              2.3e13, 3.25, 1.75, 40, 5.9e13),
}
# Published high-current designs: inputs and their own simulated yield.
MAGLIF_DESIGNS = {
    "Slutz et al. 2010, D-T gas, 27 MA (1-D: 0.5 MJ)": {"current_ma": 27, "liner_um": 3480, "liner_aspect": 6,
                                                            "length_mm": 5, "fill_mg_cc": 3, "bz_tesla": 30, "start_kev": 0.25},
    "the same at 45 MA (~10 MJ)": {"current_ma": 45, "liner_um": 3480, "liner_aspect": 6, "length_mm": 5, "fill_mg_cc": 3,
                                       "bz_tesla": 30, "start_kev": 0.25},
    "the same at 60 MA (>50 MJ)": {"current_ma": 60, "liner_um": 3480, "liner_aspect": 6, "length_mm": 5, "fill_mg_cc": 3,
                                       "bz_tesla": 30, "start_kev": 0.25},
    "Sefkow et al. 2014, D-T ice, 70 MA (1-D: 7.9 GJ)": {"current_ma": 70, "rise_ns": 110, "liner_um": 4880,
                                                             "liner_aspect": 6, "length_mm": 10, "fill_mg_cc": 5, "bz_tesla": 8,
                                                             "ice_um": 740, "start_kev": preheat_kev(25000, 5, 3330, 10,
                                                                                               2.5)},
    "Pacific Fusion 2025, D-T ice, 56 MA (109 MJ, Be)": {"current_ma": 56, "liner_um": 4430, "liner_aspect": 6,
                                                             "length_mm": 6, "fill_mg_cc": 3, "bz_tesla": 15, "ice_um": 300,
                                                             "start_kev": preheat_kev(18000, 3, 3390, 6, 2.5)},
}


def liner() -> dict:
    """Calibration on Z's MagLIF shots, the published designs, and a map of liner designs."""
    jobs = {}
    for name, (kw, *_) in MAGLIF_SHOTS.items():
        jobs[f"shot|{name}"] = dict(kw, fuel="DD", hold_ns=100)
    for name, kw in MAGLIF_DESIGNS.items():
        jobs[f"design|{name}"] = dict(kw, fuel="DT", hold_ns=150)
    for cur, rad, fill, bz, ice in itertools.product((26, 45, 60, 80), (3000, 4500, 6000), (1, 3), (10, 30),
                                                    (0, 300, 600)):
        jobs[f"map|{cur}|{rad}|{fill}|{bz}|{ice}"] = {"fuel": "DT", "current_ma": cur, "liner_um": rad, "liner_aspect": 6,
                                                          "fill_mg_cc": fill, "bz_tesla": bz, "ice_um": ice, "start_kev": 0.2,
                                                          "hold_ns": 150}
    names = list(jobs)
    res = run_all([jobs[n] for n in names])
    out = {}
    for n, r in zip(names, res):
        kw = jobs[n]
        p = Z.Pinch(**kw)
        r["convergence"] = p.fuel_radius_um / r["min_radius_um"]
        r["preheat_J"] = 3 * r["line_density"] * kw["start_kev"] * 1.602176634e-16 * p.length_mm * 1e-3 / 0.5
        out[n] = r
    return out


def fmt(x: float) -> str:
    if x == 0:
        return "0"
    return f"{x:.2g}" if x >= 0.01 else f"{x:.1e}"


def energy(x: float) -> str:
    if x >= 1e9:
        return f"{x / 1e9:.3g} GJ"
    return f"{x / 1e6:.3g} MJ" if x >= 1e5 else (f"{x / 1e3:.3g} kJ" if x >= 1e2 else f"{x:.2g} J")


def write_liner(lin: dict, mp_: dict) -> list[str]:
    def mev(r: dict) -> str:
        return energy(r["yield_J"])

    lines = ["## With a liner (MagLIF)", "",
             ("A beryllium tube carries the current and implodes onto laser-preheated fuel gas with an axial field, "
              "and optionally a layer of frozen D-T on its inner wall. The liner is incompressible and thick; the "
              "field leaks out as the fuel converges (it grows as convergence^1.57, matching MagLIF's measured B R); "
              "heat and alphas that reach the ice ablate it into the hot fuel. No current flows in the fuel, so "
              "there is no m = 0 mode; the liner's own instability is charged as a demand on convergence."), "",
             "**Calibration on Z's MagLIF shots** (deuterium):", "",
             ("| Shot | D-D neutrons: model / clean 2-D / measured | Ion temperature: model / measured | "
             "Burn width: model / measured | Convergence: model / measured |"), "|---|---|---|---|---|"]
    ratios = []
    for name, (_, measured, t_meas, bw_meas, cr_meas, clean2d) in MAGLIF_SHOTS.items():
        r = lin[f"shot|{name}"]
        ratios.append(r["neutrons_dd"] / clean2d)
        lines.append(f"| {name} | {fmt(r['neutrons_dd'])} / {fmt(clean2d)} / {fmt(measured)} | "
                     f"{r['T_at_peak_keV']:.1f} / {t_meas} keV | {r['burn_width_ns']:.2f} / {bw_meas} ns | "
                     f"{r['convergence']:.0f} / {cr_meas} |")
    lines += ["",
              (f"* **The 0-D liner overshoots at Z's scale**, {min(ratios):.0f}–{max(ratios):.0f}× the clean 2-D "
               "yields: too hot, and its burn far too brief. An incompressible liner bounces off the fuel at once; a "
               "real one decelerates layer by layer and holds the fuel for ~2 ns."), "",
              "**The published high-current designs**, rerun here:", "",
              "| Design (its own simulated yield) | Here, ice by ablation | Here, ice burnt as an ICF shell |",
              "|---|---|---|"]
    for name in MAGLIF_DESIGNS:
        r = lin[f"design|{name}"]
        upper = r["yield_J"] + r["ice_burn_bound_J"]
        lines.append(f"| {name} | {mev(r)} | {energy(upper) if r['ice_burn_bound_J'] else '–'} |")
    lines += ["",
              ("* **At high current the model undershoots**: the gas designs by ~10×, and Sefkow's 70 MA ice design "
               "doesn't ignite at all. The same too-brief stagnation is the likely cause."),
              ("* **Pacific Fusion's 56 MA design** is matched only if the ice burns like an ICF shell. So this "
               "model is for comparing designs and seeing which demand binds, not for absolute yields; the "
               "published 1-D and 2-D designs are the better guide."), "",
              ("**Liner designs for energy** (50/50 D-T, 0.2 keV preheat, best over liner radius 3–6 mm, gas 1–3 "
               "mg/cm³ and field 10–30 T). Electricity out ÷ in, ice burnt by ablation / as an ICF shell:"), "",
              "| Current | Gas only | 300 µm of ice | 600 µm of ice | Gas column, 10 Alfvén times |",
              "|---|---|---|---|---|"]
    for cur in (26, 45, 60, 80):
        cells = []
        for ice in (0, 300, 600):
            runs = [r for k, r in lin.items() if k.startswith(f"map|{cur}|") and k.endswith(f"|{ice}")]
            lo = max(gain(r) for r in runs)
            hi = max(gain(r, ice=True) for r in runs)
            cells.append(fmt(lo) if ice == 0 else f"{fmt(lo)} / {fmt(hi)}")
        lines.append(f"| {cur} MA | {' | '.join(cells)} | {fmt(gain(mp_[f'gas column|{cur}|10']))} |")
    def best_at(cur: int, ice_bound: bool) -> dict:
        return max((r for k, r in lin.items() if k.startswith(f"map|{cur}|")), key=lambda r: gain(r, ice=ice_bound))

    def describe(r: dict, ice_bound: bool) -> str:
        i = r["inputs"]
        field = i["bz_tesla"] * r["convergence"] ** 1.57 / 1e3
        return (f"{fmt(gain(r, ice=ice_bound))}× the electricity in, at {i['current_ma']} MA: a "
                f"{i['liner_um'] / 1e3:g} mm liner, {i['ice_um']} µm of ice, {i['fill_mg_cc']} mg/cm³ of gas and "
                f"{i['bz_tesla']} T; {energy(r['yield_J'] + (r['ice_burn_bound_J'] if ice_bound else 0))} per cm, "
                f"convergence {r['convergence']:.0f}, compressed field {field:.1f} kT")

    b80, b60, b45 = best_at(80, False), best_at(60, False), best_at(45, False)
    gave_up = sum(1 for k, r in lin.items() if k.startswith("map|") and r["stopped"] == "solver gave up")
    lines += ["",
              ("* **Thick ice burns at high current, by the model's own ablation.** Below ~50 MA it barely ablates "
               f"(best {fmt(gain(b45))} at 45 MA): the field that holds the heat in keeps it from the ice. At 60–80 "
               "MA the burning fuel ablates the ice fast enough for the burn to spread."),
              f"* **The best liner design, untested ideas left out:** {describe(b80, False)}.",
              (f"* **Its demands are modest except the current:** convergence {b80['convergence']:.0f} against "
               "MagLIF's measured 30–40, and a compressed field below the 7.1 kT shown on NIF. The bare gas column "
               "needed 160–260× convergence for less."),
              f"* At 60 MA: {describe(b60, False)}.",
              (f"* If the ice burns like an ICF shell (untested), the best is {fmt(gain(best_at(60, True), ice=True))} "
               "at 60 MA."),
              (f"* {gave_up} of the {sum(1 for k in lin if k.startswith('map|'))} liner designs hit the solver's "
               "work cap and count only what they made before it." if gave_up else ""), ""]
    return lines


def write(d: dict) -> list[str]:
    lev, cal, exp, mp_, tr = d["levers"], d["calibration"], d["experiment"], d["map"], d["tritium"]
    lines = ["# Making the dense Z-pinch testable, and what it takes to make energy", "",
             ("Produced by `tools/zpinch_energy.py` with `nanofusion/improved/zpinch.py`. The proposal is "
              "[docs/zpinch_experiment.md](../../zpinch_experiment.md)."), "",
             "## Levers that don't help the fibre", "",
             "A 960 µm D-T fibre at 60 MA (fully flow-stabilised: 1.9× the electricity in):", "",
             "| Lever | Fusion | Current reached | Ended | Electricity out ÷ in |", "|---|---|---|---|---|"]
    for name, r in lev.items():
        lines.append(f"| {name} | {energy(r['yield_J'])} | {r['peak_current_MA']:.1f} MA | {r['stopped']} at "
                     f"{r['t_stop_ns']:.0f} ns | {fmt(gain(r))} |")
    rise = [r for n, r in lev.items() if n.startswith("rising")]
    lines += ["",
              (f"* **Keeping the current rising** (NRL's stabilisation) fails whatever the driver: the column "
              f"radiates, cools and turns resistive, so the current peaks at {min(r['peak_current_MA'] for r in rise):.0f}–"
              f"{max(r['peak_current_MA'] for r in rise):.0f} MA and falls."),
              "* **A longer fibre** has more resistance and inductance, so it draws less current.",
              ("* **Recovering the circuit's energy** after the burn is a different machine: an inductive store that "
              "takes the field's energy back. Krakowski's 0-D fibre-pinch study finds recovering half of it halves "
              "the fusion gain needed (14 → 7.5). Nobody has shown it at this scale; Helion claims >95% without "
              "published numbers."),
              "* **The limit is the load:** a dense, cold column can't take the current.", "",
              "## A different load: a pre-formed gas column", "",
              ("A column of gas a centimetre or two wide (a gas puff, as Z already fires), ionised by the current. It "
              "conducts easily and is heated by its own implosion, so it needs no preheat laser."), "",
              ("**Calibration on Z's deuterium gas puffs** (record 3.9×10¹³ D-D neutrons at 15–17 MA; ion "
              "temperature 9.5 keV; stagnation radius 2.5 mm). Model at 17 MA, 2 cm long:"), "",
              "| Gas column | Held 1 Alfvén time | 2 | 3 |", "|---|---|---|---|"]
    for mg, rad in itertools.product((0.25, 0.5, 1.0), (10, 20)):
        cells = []
        for lim in (1, 2, 3):
            r = cal[f"{mg} mg/cm, {rad} mm, {lim} Alfvén"]
            cells.append(f"{fmt(r['neutrons'])} n, {r['peak_T_keV']:.1f} keV, {r['min_radius_um'] / 1e3:.2g} mm")
        lines.append(f"| {mg} mg/cm, {rad} mm | {' | '.join(cells)} |")
    grid = list(itertools.product((0.25, 0.5, 1.0), (10, 20)))
    one = [cal[f"{mg} mg/cm, {rad} mm, 1 Alfvén"]["neutrons"] for mg, rad in grid]
    two = [cal[f"{mg} mg/cm, {rad} mm, 2 Alfvén"]["neutrons"] for mg, rad in grid]
    lines += ["",
              (f"* **The record sits between one and two Alfvén times** in the model: one gives at most "
              f"{fmt(max(one))} neutrons, two at least {fmt(min(two))}. The pinch heats in its second Alfvén time, "
              "so the yield climbs steeply there."),
              ("* **The model converges too far.** Held two Alfvén times, its columns stagnate at 0.1–0.2 mm, "
              "against the 2.5 mm measured: a perfectly uniform column squeezes harder than a real, unstable "
              "implosion. The simulator charges that squeeze as a demand on technology."), "",
              "## The experiment Z can do now: deuterium gas columns with a 1% tritium trace", "",
              "26 MA, 1 cm wide, 2 cm long. Each cell: D-D neutrons / D-T neutrons from the trace / narrowest radius.",
              "", f"| Fuel load | {' | '.join(STABILITY)} |", "|---|" + "---|" * len(STABILITY)]
    for mg in (0.5, 2.0, 7.0):
        cells = []
        for name in STABILITY:
            r = exp[f"{mg} mg/cm | {name}"]
            cells.append(f"{fmt(r['neutrons_dd'])} / {fmt(r['neutrons_dt'])} / {r['min_radius_um'] / 1e3:.2g} mm")
        lines.append(f"| {mg} mg/cm | {' | '.join(cells)} |")
    lines += ["",
              ("* Gas puffs are Z's routine load, deuterium needs no new licence, and the 1% trace is the next step "
              "after Z's 2016 tritium shot."),
              ("* **What it measures** is the number the energy map below turns on: how many Alfvén times a dense "
              "pinch holds. Neutron yield, the D-T/D-D ratio (ion temperature), X-ray imaging of the stagnated "
              "column and its duration give it directly."), "",
              "## What it takes to make energy", "",
              ("Best electricity out ÷ in over fuel load (2–20 mg/cm) and radius (3–20 mm) for the gas column, and "
               "radius (0.24–1.5 mm) for a preheated fibre. 50/50 D-T, a lithium blanket, 40% of heat to "
               "electricity, 15% of electricity reaching the load, and each shot paying for the whole bank."), "",
              f"| Load | Current | {' | '.join(STABILITY)} |", "|---|---|" + "---|" * len(STABILITY)]
    for load in ("gas column", "fibre"):
        for cur in (26, 45, 60, 80):
            cells = [fmt(gain(mp_[f"{load}|{cur}|{name}"])) for name in STABILITY]
            lines.append(f"| {load} | {cur} MA | {' | '.join(cells)} |")
    best_key = max(mp_, key=lambda k: gain(mp_[k]) if not k.endswith("sheared flow") else 0)
    b = mp_[best_key]
    bi = b["inputs"]
    making = [mp_[f"gas column|{cur}|{name}"] for cur in (60, 80) for name in ("3", "10")]
    squeeze = [r["inputs"]["radius_um"] / r["min_radius_um"] for r in making if gain(r) > 0.5]
    held = [r["hot_alfven_times"] for r in making if gain(r) > 0.5]
    lines += ["",
              (f"* **The gas column needs far less stability than the fibre.** With 3–10 Alfvén times it reaches "
               f"{fmt(gain(mp_['gas column|60|3']))} at 60 MA and {fmt(gain(mp_['gas column|80|10']))} at 80 MA; "
               "the fibre needs about fifty (section 24)."),
              (f"* **But it needs deep convergence.** The columns that make energy burn for "
               f"{min(held):.1f}–{max(held):.1f} Alfvén times after squeezing {min(squeeze):.0f}–{max(squeeze):.0f} "
               "times, to 40–130 µm. Z's gas puffs stagnate at 2.5 mm, and the best implosion of any kind (NIF) "
               "converges about 35 times."),
              f"* **The best tested-physics design**, a {best_key.split('|')[0]} at {best_key.split('|')[1]} MA "
              f"held {best_key.split('|')[2]} Alfvén times: "
              f"{fmt(gain(b))}× the electricity in, "
              + (f"{bi['ions_per_m'] * 2.5 * AMU * 1e4:.3g} mg/cm in a {bi['radius_um'] / 1e3:.3g} mm column"
                 if "ions_per_m" in bi else f"a {bi['radius_um']:.0f} µm fibre")
              + f", {energy(b['yield_J'])} per cm, {b['burn_up']:.0%} burnt, from a {energy(b['bank_J'])} bank.",
              "* **On Z (26 MA) nothing reaches breakeven**, whatever the load or stability.", "",
              "**Tritium share** (gas column, 10 Alfvén times, fuel load and radius optimised):", "",
              "| Current | 50% | 20% | 10% | 5% | deuterium only |", "|---|---|---|---|---|---|"]
    for cur in (60, 80):
        cells = [fmt(gain(tr[f"{cur}|{f_t}"])) for f_t in (0.5, 0.2, 0.1, 0.05, 0.0)]
        lines.append(f"| {cur} MA | {' | '.join(cells)} |")
    lines += ["",
              ("* Lean tritium costs gain steeply: D-T burns ~100 times faster than D-D at these temperatures, and "
              "the tritium D-D breeds in place is too little to make up for it."),
              "",
              ("**A more efficient driver.** Gain scales with the share of electricity that reaches the load: 15% "
              "here (about Z), 17–24% stored-to-load in the Z 800 and Z 300 designs (Stygar et al. 2015), 60% "
              "assumed in LLNL's 2006 Z-pinch plant study."), "",
              "| Design | 15% | 24% | 60% |", "|---|---|---|---|"]
    for cur, name in ((60, "3"), (80, "10")):
        g = gain(mp_[f"gas column|{cur}|{name}"])
        lines.append(f"| gas column, {cur} MA, {name} Alfvén times | {fmt(g)} | {fmt(g * 24 / 15)} | "
                     f"{fmt(g * 60 / 15)} |")
    lines += ["",
              ("**Against MagLIF.** Pacific Fusion's published D-T ice MagLIF designs give 109–380 MJ at 56 MA and "
              f"816 MJ at 65 MA. The gas column gives {energy(mp_['gas column|60|3']['yield_J'])} per cm at 60 MA "
              "here: the same class, without a liner, but only if the bare column holds its convergence and a few "
              "Alfvén times."), ""]
    lines += write_liner(d["liner"], mp_)
    return lines


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / "energy.json"
    if "--reuse" in sys.argv and cache.exists():
        d = json.loads(cache.read_text())
    else:
        d = compute()
        cache.write_text(json.dumps(d, indent=1, default=float))
    lines = write(d)
    (OUT / "energy.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
