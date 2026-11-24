"""Can the ball of exploded crystals be held together? Tamper, magnetic field, self-heating.

The bare ball of docs/energy.md (section 3) needs 11-78 MJ to break even,
because it flies apart at its own sound speed. Magnetised target fusion and
MagLIF hold hot fuel with a heavy wall and insulate it from that wall with a
magnetic field. This script applies the same three tricks to an
uncompressed, solid-density D-T ball: a thick gold tamper, a kilotesla
field, and alpha self-heating. Model: nanofusion/improved/confinement.py.

Writes docs/results/confined/summary.md.

    python3 tools/confined_ball.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import confinement as C  # noqa: E402
from nanofusion.improved import energy as E  # noqa: E402
from nanofusion.improved import physics as P  # noqa: E402

OUT = Path("docs/results/confined")
HOT_ION_TEMPERATURES = (10.0, 14.0, 20.0, 28.0, 40.0, 56.0, 80.0)

# --- inputs from the literature ----------------------------------------------------
# MagLIF on the Z machine: deuterium at stagnation.
MAGLIF = [
    # Gomez et al., PRL 113, 155003 (2014): 2.5 +- 0.8 keV, 0.4 +- 0.2 g/cm3, column radius 50 um, height 6 mm,
    # burn about 2 ns, 2.0e12 D-D neutrons.
    {"shot": "2014", "t_kev": 2.5, "t_err": 0.8, "ions_m3": 0.4e3 / (P.MASS_U["D"] * P.AMU),
     "volume_m3": math.pi * (50e-6) ** 2 * 6e-3, "burn_s": 2.0e-9, "neutrons": 2.0e12},
    # Best shot (z3289): Gomez et al., PRL 125, 155002 (2020), with the stagnation conditions inferred by Knapp et
    # al., Phys. Plasmas 29, 052711 (2022): 2.6 +- 0.3 keV, 1.8 Gbar, 0.05 mm3, 2.0 ns, 1.1e13 neutrons,
    # stated as 2 kJ D-T equivalent.
    {"shot": "2020", "t_kev": 2.6, "t_err": 0.3, "ions_m3": 1.8e14 / (2 * 2.6 * P.J_PER_KEV),
     "volume_m3": 0.05e-9, "burn_s": 2.0e-9, "neutrons": 1.1e13, "dt_equivalent_J": 2.0e3,
     "fuel_heat_J": 1.5 * 1.8e14 * 0.05e-9},
]
# Tamped fuel with no field, as others have estimated it:
# Ruhl & Korn (arXiv:2302.06562): D-T at 0.5 g/cm3 in a 1 mm x 1 mm cylinder inside a high-Z case, 20 keV, about
# 3.4 MJ for fusion to exceed the fuel's energy.
# Lackner et al. (arXiv:2312.13429): perfectly tamped D-T ignites with 250 MJ at 0.5 g/cm3 and 700 MJ at 0.3 g/cm3.
PUBLISHED_TAMPED = {"ruhl_korn_J": 3.4e6, "lackner_J": {0.5: 250e6, 0.3: 700e6}}
# Short-pulse lasers that exist: NIF ARC, 4.4 kJ in four 30 ps beamlets; LFEX, 2 kJ in 1 ps; ELI-NP, 10.2 PW
# (251 J in 22 fs). Strongest laser-driven coil: 800 T claimed (Santos et al., New J. Phys. 17, 083051 (2015)).
BIGGEST_PICOSECOND_LASER_J = 4.4e3
BIGGEST_PEAK_POWER_W = 1.02e16
STRONGEST_COIL_FIELD_T = 800.0
COUPLING = 0.15  # share of the laser light assumed to end up as heat in the fuel


def two_digits(value: float) -> str:
    """Two significant digits without an exponent: 130, 4.1, 0.28."""
    return f"{float(f'{value:.2g}'):,g}"


def scaled(value: float, units: tuple) -> str:
    for unit, scale in units:
        if value >= scale:
            return f"{two_digits(value / scale)} {unit}"
    unit, scale = units[-1]
    return f"{two_digits(value / scale)} {unit}"


def energy_text(joules: float) -> str:
    return scaled(joules, (("TJ", 1e12), ("GJ", 1e9), ("MJ", 1e6), ("kJ", 1e3), ("J", 1.0)))


def power_text(watts: float) -> str:
    return scaled(watts, (("EW", 1e18), ("PW", 1e15), ("TW", 1e12)))


def sci(value: float) -> str:
    """4e18 -> 4×10¹⁸."""
    mantissa, exponent = f"{value:.0e}".split("e")
    return f"{mantissa}×10" + str(int(exponent)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))


def sci2(value: float) -> str:
    """1.1e13 -> 1.1×10¹³."""
    mantissa, exponent = f"{value:.1e}".split("e")
    return f"{mantissa}×10" + str(int(exponent)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))


def time_text(seconds: float) -> str:
    return scaled(seconds, (("µs", 1e-6), ("ns", 1e-9), ("ps", 1e-12)))


def cell(r: dict | None) -> str:
    if r is None:
        return "never"
    start = f"{r['t_i_keV']:g} keV" if r["t_i_keV"] == r["t_e_keV"] else f"ions {r['t_i_keV']:g} keV"
    return f"{energy_text(r['energy_in_J'])} (radius {two_digits(r['radius_m'] * 1e3)} mm, {start})"


def slim(r: dict | None) -> dict | None:
    return None if r is None else {k: v for k, v in r.items() if k != "history"}


def validation() -> tuple[list[str], dict]:
    lines = ["## Check: the bare ball of section 3", "",
             "With no tamper, no field and no self-heating, the new time-dependent model must give the thermal burn "
             "of `energy.array_burn` (fixed temperature, burn time R / 4c_s).", "",
             "| ball radius | temperature | section 3 model: thermal burn / energy in | this model |", "|---|---|---|---|"]
    rows = []
    for radius in (1e-4, 1e-3, 3e-3):
        for t in (11.0, 31.0):
            plasma = E.fuel_plasma("DT", 1.0, t, t)
            volume = 4 / 3 * math.pi * radius**3
            old = (plasma.fusion_power()[0] * radius / (4 * plasma.sound_speed())
                   / (1.5 * (plasma.n_e + sum(plasma.ions.values())) * t * E.KEV))
            new = C.burn(C.Ball(radius_m=radius, t_i_kev=t, alpha_heating=False))["gain"]
            rows.append({"radius_m": radius, "t_keV": t, "array_burn": old, "confinement": new, "volume_m3": volume})
            lines.append(f"| {radius * 1e3:g} mm | {t:g} keV | {old:.3g} | {new:.3g} |")
    lines.append("")
    return lines, {"rows": rows}


def tricks() -> tuple[list[str], dict]:
    cases = [
        ("bare ball, no self-heating (the model of section 3)", dict(tamper="none", alpha_heating=False)),
        ("bare ball with alpha self-heating", dict(tamper="none")),
        ("thick gold tamper, no field", dict(tamper="gold")),
        ("gold tamper + 100 T", dict(tamper="gold", b_tesla=100.0)),
        ("gold tamper + 1 kT", dict(tamper="gold", b_tesla=1e3)),
        ("gold tamper + 3 kT", dict(tamper="gold", b_tesla=3e3)),
        ("gold tamper + 10 kT", dict(tamper="gold", b_tesla=1e4)),
        ("gold tamper + 30 kT", dict(tamper="gold", b_tesla=3e4)),
    ]
    lines = ["## 1. What each trick buys", "",
             "Solid-density D-T, ions and electrons starting at the same temperature. For each kind of ball, the "
             "cheapest one (least energy in the fuel, over starting temperature and radius) that returns 1, 10 or "
             "100 times that energy as fusion. Energy in = heat of the ions and electrons plus the field energy "
             "inside the fuel. It leaves out every loss on the way in.", "",
             "| ball | fusion = energy in | fusion = 10 × energy in | fusion = 100 × energy in |", "|---|---|---|---|"]
    data = {}
    for name, settings in cases:
        row = [C.smallest_ball(g, r_hi=0.5, **settings) for g in (1, 10, 100)]
        data[name] = [slim(r) for r in row]
        lines.append(f"| {name} | " + " | ".join(cell(r) for r in row) + " |")
    bare = data["bare ball with alpha self-heating"]
    best = min((v for k, v in data.items() if "kT" in k), key=lambda v: v[0]["energy_in_J"])
    ratio = two_digits(bare[0]["energy_in_J"] / best[0]["energy_in_J"])
    lines += ["", f"* **The tamper and the field together cut breakeven about {ratio}-fold**, "
              f"from {energy_text(bare[0]['energy_in_J'])} to {energy_text(best[0]['energy_in_J'])}. Neither does it "
              "alone: a tamped ball with no field loses its heat to the gold, and a field cannot hold fuel whose "
              "pressure is a hundred times its own.",
              "* **The best field is 3–10 kT.** Below 1 kT too much heat still leaks and the alphas escape. Above "
              "10 kT the field itself costs as much energy as the heat.",
              "* **High gain stays expensive.** The wall recedes in about a nanosecond, and at solid density an "
              "alpha needs about that long to hand over its energy, so small balls never run away (ignite). A "
              "hundredfold return needs tens of megajoules in a centimetre-sized ball inside centimetres of gold.",
              "* **Where the model stops.** Each burn is followed until the ball's radius has grown fivefold; the "
              "wall model is not trusted further. Balls that ignite are still burning then, so the hundredfold "
              "column understates them (followed ten times further, they return 20–50% more).", ""]
    dd_temps = (20.0, 40.0, 60.0, 100.0, 150.0)
    dd = {"bare": C.smallest_ball(1, fuel="DD", temperatures=dd_temps, r_hi=3.0, tamper="none"),
          "tamped": C.smallest_ball(1, fuel="DD", temperatures=dd_temps, r_hi=3.0, tamper="gold", b_tesla=3e3)}
    data["deuterium only"] = {k: slim(v) for k, v in dd.items()}
    lines += [f"**Without tritium.** Solid deuterium alone, counted generously (every triton and helium-3 it makes "
              f"burns at once): breakeven needs {cell(dd['bare'])} bare, and {cell(dd['tamped'])} with the tamper "
              "and 3 kT. The tricks help just as much, but deuterium burns about a hundred times more slowly, so it "
              "stays out of reach. Lithium deuteride is worse: it radiates more than it fuses (section 3), and "
              "holding it longer does not change that.", ""]
    return lines, data


def history(design: dict) -> tuple[list[str], dict]:
    ball = C.Ball(radius_m=design["radius_m"], t_i_kev=design["t_i_keV"], t_e_kev=design["t_e_keV"],
                  tamper="gold", b_tesla=design["b_tesla"])
    res = C.burn(ball, history=True)
    track = res["history"]
    lines = ["## 2. How the cheapest breakeven ball burns", "",
             f"Gold tamper, {design['b_tesla'] / 1e3:g} kT, radius {design['radius_m'] * 1e3:.2g} mm, "
             f"{energy_text(res['energy_in_J'])} in ({res['mass_g'] * 1e6:.2g} µg of D-T, "
             f"B·R = {res['b_r_tesla_m']:.2g} T m).", "",
             "| time | ion temperature | electron temperature | radius / starting radius | fusion so far / energy in |",
             "|---|---|---|---|---|"]
    times = [p[0] for p in track]
    for t_ns in (0.0, 0.1, 0.3, 1.0, 2.0, 4.0, 8.0, 16.0):
        row = [float(np.interp(t_ns * 1e-9, times, [p[i] for p in track])) for i in range(1, 5)]
        lines.append(f"| {t_ns:g} ns | {row[0]:.3g} keV | {row[1]:.3g} keV | {row[2] / design['radius_m']:.3g} "
                     f"| {row[3] / res['energy_in_J']:.2g} |")
    lines += ["", f"The fuel doubles its volume after {res['time_to_double_volume_s'] * 1e9:.2g} ns. "
              f"Of the energy put in, {res['conduction_J'] / res['energy_in_J']:.0%} is conducted to the wall and "
              f"{res['bremsstrahlung_J'] / res['energy_in_J']:.0%} radiated before the burn ends; the rest pushes "
              f"the gold. The shock in the gold travels {res['tamper_thickness_needed_m'] * 1e3:.2g} mm in that time, "
              f"which is how thick the tamper must be. {res['burn_fraction']:.1%} of the fuel burns.", ""]
    return lines, slim(res)


def less_ideal(base_b: float) -> tuple[list[str], dict]:
    base = dict(tamper="gold", b_tesla=base_b)
    cases = [
        ("as in section 1", dict(base)),
        ("ions heated first, electrons start at 1 keV (what exploding crystals would do)",
         dict(base, t_e_kev=1.0, temperatures=HOT_ION_TEMPERATURES)),
        ("Bohm heat loss across the field instead of classical", dict(base, bohm=True)),
        ("wall recedes 1.5 × faster than the shock formula says", dict(base, wall_speed_factor=1.5)),
        ("gold shell only 0.1 mm thick", dict(base, tamper_mm=0.1)),
        ("gold shell 0.3 mm thick", dict(base, tamper_mm=0.3)),
        ("entrance hole, 1% of the wall", dict(base, hole_fraction=0.01)),
        ("entrance hole, 2% of the wall", dict(base, hole_fraction=0.02)),
        ("entrance hole, 5% of the wall", dict(base, hole_fraction=0.05)),
        ("gold mixed into the fuel, 1 atom per 10,000 fuel ions", dict(base, mix=1e-4)),
        ("gold mixed into the fuel, 1 atom per 1,000", dict(base, mix=1e-3)),
        ("fuel at 20% of solid density (a foam that light can enter)", dict(base, fill=0.2)),
        ("20% foam, ions heated first, 1% hole", dict(base, fill=0.2, hole_fraction=0.01, t_e_kev=1.0,
                                                     temperatures=HOT_ION_TEMPERATURES)),
    ]
    lines = ["## 3. Making it less ideal", "",
             f"The same search with a {base_b / 1e3:g} kT field and one thing changed at a time.", "",
             "| change | fusion = energy in | fusion = 10 × energy in |", "|---|---|---|"]
    data = {}
    for name, settings in cases:
        row = [C.smallest_ball(g, r_hi=0.5, **settings) for g in (1, 10)]
        data[name] = [slim(r) for r in row]
        lines.append(f"| {name} | " + " | ".join(cell(r) for r in row) + " |")
    try:
        E.DT_SPIN_FACTOR = 1.5
        row = [C.smallest_ball(g, r_hi=0.5, **base) for g in (1, 10)]
    finally:
        E.DT_SPIN_FACTOR = 1.0
    data["spin-polarised fuel"] = [slim(r) for r in row]
    lines.append("| spin-polarised fuel (cross section × 1.5, untested) | " + " | ".join(cell(r) for r in row) + " |")
    lines += ["", "* **The entrance hole matters most.** The gold wall recedes 13 times more slowly than the fuel's "
              "sound speed, so fuel streaming out of even a 2% hole is a comparable loss.",
              "* **Heating the ions first helps**, because the electrons take about a nanosecond to catch up and "
              "meanwhile the ions burn at twice the temperature for the same energy.",
              "* **Bohm diffusion and a little gold mix cost almost nothing at breakeven.** The wall's motion sets "
              "the clock there, not heat loss. Both matter more for the larger, longer-burning balls.",
              "* **The wall speed is the number everything hangs on.** The energy needed goes as its cube.", ""]
    return lines, data


def straight_field() -> tuple[list[str], dict]:
    lines = ["## 4. A straight field instead of a ring field", "",
             "Everything above assumes the field lies along the wall everywhere. Only a ring-shaped field around "
             "an axis does that, and making 3 kT at 0.3 mm that way takes about 4.5 MA through the fuel. Laser-driven "
             "coils make a straight field instead. Heat and alphas then run along it to the two end walls at the "
             "unmagnetised rate, so the fuel must be a long column, as in MagLIF. The cheapest breakeven column "
             "for each field, over its radius, its length (2–32 diameters) and its starting temperature:", "",
             "| field | energy in at breakeven | column | lost to the walls | with ions heated first |",
             "|---|---|---|---|---|"]
    data = {}
    for b in (1e3, 3e3, 1e4):
        best = {}
        for label, extra in (("equal", {}), ("ions first", {"t_e_kev": 1.0, "temperatures": HOT_ION_TEMPERATURES})):
            for aspect in (2.0, 4.0, 8.0, 16.0, 32.0):
                r = C.smallest_ball(1, r_hi=0.05, tamper="gold", b_tesla=b, aspect=aspect, field="axial", **extra)
                if r and (label not in best or r["energy_in_J"] < best[label]["energy_in_J"]):
                    best[label] = {**slim(r), "aspect": aspect}
        data[b] = best
        e = best["equal"]
        lines.append(f"| {b / 1e3:g} kT | {energy_text(e['energy_in_J'])} "
                     f"| {two_digits(2 * e['radius_m'] * 1e3)} mm wide, {two_digits(2 * e['aspect'] * e['radius_m'] * 1e3)} mm "
                     f"long, {e['t_i_keV']:g} keV | {e['conduction_J'] / e['energy_in_J']:.0%} "
                     f"| {energy_text(best['ions first']['energy_in_J'])} |")
    ring = C.smallest_ball(1, r_hi=0.5, tamper="gold", b_tesla=3e3)
    lines += ["", f"A straight field costs about {two_digits(data[3e3]['equal']['energy_in_J'] / ring['energy_in_J'])} "
              "times more energy than the ideal ring field at 3 kT, because the ends leak.", ""]
    return lines, data


def crystal_for_ion_energy(mean_ion_kev: float) -> tuple[float, float]:
    """D-T crystal radius (m) whose Coulomb explosion gives that mean ion energy, and its stripping intensity (W/cm^2)."""
    rho = E.charge_density("DT")
    per_r2 = 0.9 * 0.6 * P.E_CHARGE * rho / (3 * P.EPS0) / P.J_PER_KEV
    radius = math.sqrt(mean_ion_kev / per_r2)
    return radius, P.intensity_w_cm2(rho * radius / (3 * P.EPS0))


def density_scan() -> tuple[list[str], dict]:
    lines = ["## 5. Thinner fuel: from solid crystals to cluster jets", "",
             "The same tamped, magnetised ball with the fuel spread thinner, down to the average density of the "
             "gas-jet cluster experiments of section 6. For each density the field (30 T – 10 kT) that gives the "
             "cheapest breakeven. The last column is the laser power needed just to reach the stripping intensity "
             "of the crystals over the ball's cross-section.", "",
             "| fuel density (× solid) | field | energy in at breakeven | ball radius | time to double its volume "
             "| gold needed | laser power to strip the crystals |", "|---|---|---|---|---|---|---|"]
    data = {}
    for fill in (1.0, 0.2, 0.03, 0.01, 0.001):
        best = None
        for b in (30.0, 100.0, 300.0, 1e3, 3e3, 1e4):
            r = C.smallest_ball(1, fill=fill, tamper="gold", b_tesla=b, r_hi=3.0)
            if r and (best is None or r["energy_in_J"] < best["energy_in_J"]):
                best = r
        _, intensity = crystal_for_ion_energy(3 * best["t_i_keV"])
        power = intensity * math.pi * (best["radius_m"] * 100) ** 2
        data[fill] = {**slim(best), "strip_intensity_W_cm2": intensity, "strip_power_W": power}
        lines.append(f"| {fill:g} | {best['b_tesla'] / 1e3:g} kT | {energy_text(best['energy_in_J'])} "
                     f"| {two_digits(best['radius_m'] * 1e3)} mm | {time_text(best['time_to_double_volume_s'])} "
                     f"| {two_digits(best['tamper_thickness_needed_m'] * 1e3)} mm | {power_text(power)} |")
    fills = sorted(data)
    slope = math.log(data[fills[0]]["energy_in_J"] / data[fills[-1]]["energy_in_J"]) / math.log(fills[0] / fills[-1])
    lines += ["", "* **Density hardly changes the energy**, unlike the bare ball (where it goes as 1/density²). "
              f"A tamped, insulated ball needs roughly energy ∝ density^{slope:.2f}. That is the central point of "
              "magnetised target fusion.",
              "* **Thin fuel buys time, not power.** The burn slows to microseconds, but the ball grows, and "
              "exploding crystals need the same intensity everywhere. Solid density needs the least laser power.", ""]
    return lines, data


def what_it_takes(designs: dict) -> tuple[list[str], dict]:
    lines = ["## 6. What it would take", "",
             f"If {COUPLING:.0%} of the laser light ended up as heat in the fuel. That is about what nanowire "
             "arrays manage: 4% of the light reaches fast deuterons, which by their share of the charge means "
             "about 15% reaches all the ions (docs/results/further/summary.md, section 5).", "",
             "| design | crystals | stripping intensity | laser energy | lowest power that strips the crystals, "
             "and the pulse length at that power | fusion / laser energy |", "|---|---|---|---|---|---|"]
    data = {}
    for name, r in designs.items():
        mean_ion = 1.5 * r["t_i_keV"] if r["t_i_keV"] != r["t_e_keV"] else 3 * r["t_i_keV"]
        radius, intensity = crystal_for_ion_energy(mean_ion)
        laser = r["energy_in_J"] / COUPLING
        power = intensity * math.pi * (r["radius_m"] * 100) ** 2
        data[name] = {"crystal_radius_m": radius, "intensity_W_cm2": intensity, "laser_J": laser, "power_W": power,
                      "pulse_s": laser / power, "target_gain": r["gain"] * COUPLING, "energy_in_J": r["energy_in_J"]}
        lines.append(f"| {name} | {2 * radius * 1e9:.0f} nm | {sci(intensity)} W/cm² | {energy_text(laser)} "
                     f"| {power_text(power)}, {time_text(laser / power)} | {r['gain'] * COUPLING:.2g} |")
    lines += ["", f"For comparison, the most energetic picosecond-class laser shot so far is "
              f"{BIGGEST_PICOSECOND_LASER_J / 1e3:.2g} kJ in 30 ps (NIF's ARC; LFEX gives 2 kJ in 1 ps), and the "
              f"highest peak power is {BIGGEST_PEAK_POWER_W / 1e15:.0f} PW, in 0.25 kJ (ELI-NP). The strongest "
              f"laser-driven coil field claimed is about {STRONGEST_COIL_FIELD_T:.0f} T, and a 2022 reanalysis "
              "puts such coils well below that.", ""]
    return lines, data


def published_checks() -> tuple[list[str], dict]:
    """Compare with a real magnetised experiment and with other people's estimates for tamped fuel."""
    lines = ["## 7. Checks against MagLIF and against published estimates", "",
             "**MagLIF** on Sandia's Z machine reaches this kind of plasma by compression: deuterium near solid "
             "density at 2.5–3 keV, in a column about 0.1 mm wide and 6 mm long, held about 2 ns by its liner, with "
             "field × radius ≈ 0.34 T m. D-D neutrons from this project's reactivity, n²⟨σv⟩/4 × volume × burn time, "
             "for the published stagnation conditions:", "",
             "| shot | conditions | neutrons, this model (range over the temperature error) | measured |",
             "|---|---|---|---|"]
    data: dict = {"maglif": []}
    for m in MAGLIF:
        def neutrons(t: float, m: dict = m) -> float:
            return 0.25 * m["ions_m3"] ** 2 * E.reactivity("D", "D", t) * m["volume_m3"] * m["burn_s"]
        mid, lo, hi = neutrons(m["t_kev"]), neutrons(m["t_kev"] - m["t_err"]), neutrons(m["t_kev"] + m["t_err"])
        rho = m["ions_m3"] * P.MASS_U["D"] * P.AMU / 1e3
        data["maglif"].append({"shot": m["shot"], "model": mid, "low": lo, "high": hi, "measured": m["neutrons"]})
        lines.append(f"| {m['shot']} | {m['t_kev']} ± {m['t_err']} keV, {rho:.1f} g/cm³, {m['volume_m3'] * 1e9:.2g} mm³, "
                     f"{m['burn_s'] * 1e9:g} ns | {sci2(mid)} ({sci2(lo)} to {sci2(hi)}) | {sci2(m['neutrons'])} |")
    best = MAGLIF[1]
    ratio = best["dt_equivalent_J"] / best["fuel_heat_J"]
    data["maglif_dt_equivalent_over_heat"] = ratio
    lines += ["", f"The best shot is equivalent to {energy_text(best['dt_equivalent_J'])} of D-T fusion from fuel "
              f"holding {energy_text(best['fuel_heat_J'])} of heat: a ratio of {ratio:.2f}. The designs above are the "
              "same kind of plasma at 20–40 keV instead of 3 keV, which nobody knows how to reach without compressing.",
              "", "**Tamped fuel with no field**, as others have estimated it:", ""]
    # Ruhl & Korn: 0.5 g/cm3, cylinder with radius = length = 1 mm, 20 keV, high-Z case
    fill = lambda rho: rho * 1e3 / (C.SOLID_DT_IONS * C.A_FUEL * P.AMU)  # noqa: E731
    rk = C.burn(C.Ball(radius_m=1e-3, t_i_kev=20.0, fill=fill(0.5), tamper="gold", aspect=0.5))
    data["ruhl_korn"] = slim(rk)
    lines.append(f"* Ruhl & Korn (arXiv:2302.06562) propose fuel with D-T at 0.5 g/cm³ in a 1 mm × 1 mm cylinder "
                 f"inside a high-Z case, heated to 20 keV, and find that about "
                 f"{energy_text(PUBLISHED_TAMPED['ruhl_korn_J'])} makes fusion exceed the fuel's energy. This model: "
                 f"{energy_text(rk['energy_in_J'])} in that fuel returns {rk['gain']:.2f} of it.")
    data["lackner"] = {}
    for rho, published in PUBLISHED_TAMPED["lackner_J"].items():
        r = C.smallest_ball(10, fill=fill(rho), tamper="gold", r_hi=0.5)
        data["lackner"][rho] = slim(r)
        lines.append(f"* Their critics (Lackner et al., arXiv:2312.13429) find that perfectly tamped D-T at {rho} g/cm³ "
                     "ignites with "
                     f"{energy_text(published)}. This model, with a gold wall rather than a perfect one: "
                     f"{energy_text(r['energy_in_J'])} for a tenfold return.")
    lines += ["", "So without a field this model agrees with both the proposal and its critics. The critics exempt "
              "schemes with an applied magnetic field, and that is the case worked out here.", ""]
    return lines, data


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    lines = ["# Holding the crystal ball together: tamper, magnetic field, self-heating", "",
             "Produced by `tools/confined_ball.py`. Model: `nanofusion/improved/confinement.py` (0-D, two "
             "temperatures: fusion with burn-up, delayed alpha heating, electron-ion exchange, bremsstrahlung, "
             "Braginskii heat loss to the wall, and a gold wall that recedes at the shock piston speed).", ""]
    data: dict = {}
    part, data["validation"] = validation()
    lines += part
    part, data["tricks"] = tricks()
    lines += part
    candidates = {k: v[0] for k, v in data["tricks"].items() if "kT" in k and v[0]}
    best_name = min(candidates, key=lambda k: candidates[k]["energy_in_J"])
    design = data["tricks"]["gold tamper + 3 kT"][0]
    part, data["history"] = history(design)
    lines += part
    part, data["less_ideal"] = less_ideal(3e3)
    lines += part
    part, data["straight_field"] = straight_field()
    lines += part
    part, data["density"] = density_scan()
    lines += part
    li = data["less_ideal"]
    straight = data["straight_field"][3e3]["ions first"]
    designs = {
        "solid D-T, ions heated first, breakeven":
            li["ions heated first, electrons start at 1 keV (what exploding crystals would do)"][0],
        "20% foam, ions heated first, 1% hole, breakeven": li["20% foam, ions heated first, 1% hole"][0],
        "20% foam, ions heated first, 1% hole, tenfold return": li["20% foam, ions heated first, 1% hole"][1],
        "solid D-T column in a straight 3 kT field, ions heated first, breakeven": straight,
    }
    part, data["what_it_takes"] = what_it_takes(designs)
    lines += part
    part, data["checks"] = published_checks()
    lines += part
    data["best_field_case"] = best_name
    text = "\n".join(lines) + "\n"
    (OUT / "summary.md").write_text(text)
    (OUT / "results.json").write_text(json.dumps(data, indent=1, default=str))
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
