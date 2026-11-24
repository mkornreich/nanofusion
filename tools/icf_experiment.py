"""A practical experiment with ordinary inertial fusion: a one-change-at-a-time campaign on NIF.

The simulator's compressed_icf model (an isobaric hot spot in a dense D-T shell,
burning rhoR / (rhoR + 7 g/cm2)) has no mix of ablator into the fuel and no
asymmetry, so it gives a *clean-implosion ceiling*. It is calibrated on NIF's
measured shots:
* the fuel mass (0.22 mg) fixes the shell's areal density at the measured
  hot-spot pressure and design adiabat;
* the energy that fuel then holds, divided by the laser energy, is the coupling
  used for every other laser energy;
* each shot's yield divided by its ceiling is its "quality" (how close it came).

The campaign changes one thing per step, in three phases: today's target with
one knob moved; a thicker fuel layer at NIF's full laser energy; NIF's planned
2.6 MJ upgrade. Each step's yield is predicted as quality x ceiling.

Writes docs/results/icf_experiment/summary.md and results.json.

    python3 tools/icf_experiment.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import simulator as S

OUT = Path("docs/results/icf_experiment")

# --- measured NIF shots (sources in docs/icf_experiment.md) ----------------------------------------------------------
# Laser and yield: LLNL "Achieving Fusion Ignition"; hot-spot pressure: Abu-Shawareb et al., PRL 129, 075001 (2022) and
# PRL 132, 065102 (2024); design adiabat: Kritcher et al., PRE 109, 025204 (2024), Table I (midpoints of 2.55-2.75 and
# 2.22-2.68); DT fuel: 64.5-65 um of ice, 220 ug, in every Hybrid-E capsule (PRL 132; Lindl et al. 2026).
SHOTS = {
    "N210808": {"laser_J": 1.917e6, "yield_J": 1.37e6, "pressure_gbar": 569.0, "adiabat": 2.65},
    "N221204": {"laser_J": 2.05e6, "yield_J": 3.15e6, "pressure_gbar": 608.0, "adiabat": 2.45},
}
# Later shots, same fuel, pressure not published.
LATER = {"N230729": (2.05e6, 3.88e6), "N240210": (2.2e6, 5.2e6), "N250406": (2.08e6, 8.6e6)}
RECORD, TYPICAL = "N250406", "N240210"  # the record's quality, and the best before the graded dopant
FUEL_MASS_KG = 0.22e-6
P0, A0 = 608.0, 2.45  # N221204's pressure and adiabat: the starting point
NIF_MAX_J = 2.2e6  # most laser energy on an ignition target so far (Oct 2023)
EYC_J = 2.6e6  # NIF's Enhanced Yield Capability upgrade (CD-1 approved 2026)
PUBLISHED_EYC = ("LLNL expects \"greater than 30 MJ\"; Lindl et al. (2026) very likely more than 20 MJ, possibly "
                 "50–70 MJ fully optimised; Hurricane's hydro-scaling 24 MJ at 3 MJ of laser")
ICE_DENSITY = 0.255  # g/cm3, solid D-T
INNER_RADIUS_MM = 1.05  # inner radius of today's diamond capsule
NIF_LASER_EFFICIENCY = 2.05 / 350.0  # NIF: ~300-400 MJ from the grid for a 2.05 MJ shot
DESIGN_LASER_EFFICIENCY = 0.10  # a fusion-plant laser
THERMAL = 0.4
BLANKET = S.ALT.blanket_energy_per_fusion_mev(0.0, uranium=False) / S.Q_DT_MEV  # lithium blanket: heat / fusion

# Steps: name, phase, what changes in practice, laser, pressure, adiabat, fuel ("same" = 0.22 mg, or "resized").
STEPS = [
    ("A. Repeat the record", 1, "N250406's target and pulse", 2.08e6, P0, A0, "same"),
    ("B. Lower adiabat", 1, "A with a lower foot in the pulse", 2.08e6, P0, 2.25, "same"),
    ("C. Higher pressure", 1, "A with a faster, rounder implosion", 2.08e6, 700.0, A0, "same"),
    ("D. Thicker fuel", 2, "full laser energy, more ice", NIF_MAX_J, P0, A0, "resized"),
    ("E. Thicker fuel, lower adiabat", 2, "D with B's pulse", NIF_MAX_J, P0, 2.25, "resized"),
    ("F. Thicker fuel, both", 2, "D with B's pulse and C's pressure", NIF_MAX_J, 700.0, 2.25, "resized"),
    ("G. Upgraded laser", 3, "D at 2.6 MJ", EYC_J, P0, A0, "resized"),
    ("H. Upgraded laser, both", 3, "F at 2.6 MJ", EYC_J, 700.0, 2.25, "resized"),
    ("I. Upgraded laser, stretched", 3, "H with 800 Gbar and adiabat 2.0", EYC_J, 800.0, 2.0, "resized"),
]


def ice_um(mass_mg: float) -> float:
    """Thickness of a D-T ice layer of this mass inside today's capsule."""
    volume = mass_mg * 1e-3 / ICE_DENSITY * 1e3  # mm3
    inner = (INNER_RADIUS_MM**3 - volume * 3 / (4 * math.pi)) ** (1 / 3)
    return (INNER_RADIUS_MM - inner) * 1e3


def fixed_mass(pressure: float, adiabat: float, mass_kg: float = FUEL_MASS_KG) -> dict:
    """The model's target holding mass_kg of D-T at this pressure and adiabat."""
    rho = S.icf_shell_density(pressure, adiabat)
    rho_r = S.icf_shell_for_mass(pressure, adiabat, mass_kg)
    t = S.icf_target(rho, rho_r, adiabat)
    return {"rho_g_cc": rho, "rho_r": rho_r, "adiabat": adiabat, "pressure_gbar": pressure,
            "fusion_J": t["fusion_J"], "fuel_J": t["hot"]["energy_J"] + t["cold_J"], "hot_spot_J": t["hot"]["energy_J"],
            "mass_mg": t["mass_kg"] * 1e6, "burnt": t["burnt"], "hot_spot_um": t["hot"]["radius_m"] * 1e6}


def calibrate() -> dict:
    out = {}
    for name, shot in SHOTS.items():
        t = fixed_mass(shot["pressure_gbar"], shot["adiabat"])
        out[name] = {**shot, **t, "total_rho_r": t["rho_r"] + S.HOT_SPOT_RHOR, "coupling": t["fuel_J"] / shot["laser_J"],
                     "ceiling_J": t["fusion_J"], "quality": shot["yield_J"] / t["fusion_J"]}
    return out


def campaign(coupling: float, qualities: dict) -> list[dict]:
    rows = []
    for name, phase, practice, laser, pressure, adiabat, fuel in STEPS:
        if fuel == "same":
            d = fixed_mass(pressure, adiabat)
            if d["fuel_J"] > coupling * laser:
                raise ValueError(f"{name}: the model needs {d['fuel_J'] / coupling / 1e6:.2f} MJ of laser")
        else:
            d = S.icf_design(coupling * laser, pressure, adiabat, adiabats=[adiabat], pressures=[pressure])
        row = {"step": name, "phase": phase, "practice": practice, "laser_J": laser, "fuel": fuel, **d,
               "ceiling_J": d["fusion_J"], "laser_needed_J": d["fuel_J"] / coupling, "ice_um": ice_um(d["mass_mg"])}
        for label, q in qualities.items():
            row[f"yield_{label}_J"] = q * d["fusion_J"]
        y = row["yield_record_J"]
        row["target_gain"] = y / laser
        row["electric_nif"] = y * BLANKET * THERMAL * NIF_LASER_EFFICIENCY / laser
        row["electric_plant_laser"] = y * BLANKET * THERMAL * DESIGN_LASER_EFFICIENCY / laser
        rows.append(row)
    return rows


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cal = calibrate()
    ref = cal["N221204"]
    coupling, ceiling = ref["coupling"], ref["ceiling_J"]
    later = {name: {"laser_J": laser, "yield_J": y, "ceiling_J": ceiling, "quality": y / ceiling}
             for name, (laser, y) in LATER.items()}
    qualities = {"record": later[RECORD]["quality"], "typical": later[TYPICAL]["quality"]}
    rows = campaign(coupling, qualities)
    by = {r["step"][0]: r for r in rows}
    q_rec, q_typ = qualities["record"], qualities["typical"]

    lines = ["# A practical experiment with ordinary inertial fusion: testing what limits NIF's yield", "",
             ("Produced by `tools/icf_experiment.py` from the simulator's `compressed_icf` model "
             "(`nanofusion/improved/simulator.py`). The proposal is [docs/icf_experiment.md](../../icf_experiment.md)."),
             "", "## Calibration on NIF's measured shots", "",
             ("The model has no mixing and no asymmetry, so it gives the yield of a *clean* implosion: a ceiling. "
             "It is matched to NIF's fuel (0.22 mg) at each shot's measured hot-spot pressure and design adiabat:"), "",
             ("| Shot | Laser | Hot-spot pressure | Adiabat | Areal density in the model | Energy in the fuel "
             "| Ceiling | Measured | Quality (measured ÷ ceiling) |"), "|---|---|---|---|---|---|---|---|---|"]
    for name, c in cal.items():
        lines.append(f"| {name} | {c['laser_J'] / 1e6:.2f} MJ | {c['pressure_gbar']:.0f} Gbar | {c['adiabat']:.2f} | "
                     f"{c['total_rho_r']:.2f} g/cm² | {c['fuel_J'] / 1e3:.0f} kJ ({c['coupling']:.1%} of the laser) | "
                     f"{c['ceiling_J'] / 1e6:.1f} MJ | {c['yield_J'] / 1e6:.2f} MJ | {c['quality']:.0%} |")
    lines += ["",
              (f"* The model's areal density, {ref['total_rho_r']:.1f} g/cm², is close to NIF's simulated 1.35–1.41 "
              "g/cm² with alpha heating turned off (Kritcher et al. 2024)."),
              (f"* The coupling, {coupling:.1%} of the laser energy into the fuel, is fitted on N221204 and held "
              "fixed. It absorbs the model's simplifications."),
              (f"* The model's hot spot is the smallest that ignites ({ref['hot_spot_um']:.0f} µm, 10 keV). NIF's "
              "burning hot spots are larger and hotter (~68 µm, 13 keV for N221204): alpha heating, which the model "
              "leaves out."), "",
              ("Later shots used the same fuel; their pressures are not published, so they are compared with "
              "N221204's ceiling:"), "",
              "| Shot | Laser | Ceiling | Measured | Quality |", "|---|---|---|---|---|"]
    for name, c in later.items():
        lines.append(f"| {name} | {c['laser_J'] / 1e6:.2f} MJ | {c['ceiling_J'] / 1e6:.1f} MJ | "
                     f"{c['yield_J'] / 1e6:.2f} MJ | {c['quality']:.0%} |")
    qs = [cal["N210808"]["quality"], ref["quality"]] + [c["quality"] for c in later.values()]
    lines += ["",
              (f"* **Quality is the story of NIF since 2021:** {' → '.join(f'{q:.0%}' for q in qs)}. Fewer capsule "
              "defects, then a graded tungsten dopant that cut the mixing of ablator into the fuel. The record "
              "(N250406) came close to the clean ceiling."),
              ("* Part of the later gains may be higher pressure, which this comparison can't separate. Measuring "
              "pressure on every shot of the campaign separates them."), "",
              "## The campaign", "",
              (f"One change per step. Yield predicted at the record's quality ({q_rec:.0%}) and at the best quality "
              f"before the graded dopant ({TYPICAL}, {q_typ:.0%}). Coupling held at {coupling:.1%}."), "",
              ("| Step | What changes in practice | Laser | Hot-spot pressure | Adiabat | D-T fuel (ice) | Ceiling "
              f"| Yield at {q_rec:.0%} | at {q_typ:.0%} | Target gain |"), "|---|---|---|---|---|---|---|---|---|---|"]
    phases = {1: "**Phase 1: today's target**", 2: "**Phase 2: a thicker fuel layer at full energy**",
              3: "**Phase 3: NIF's 2.6 MJ upgrade**"}
    for r in rows:
        if r["step"][0] in "ADG":
            lines.append(f"| {phases[r['phase']]} | | | | | | | | | |")
        lines.append(f"| {r['step']} | {r['practice']} | {r['laser_J'] / 1e6:.2f} MJ | {r['pressure_gbar']:.0f} Gbar "
                     f"| {r['adiabat']:.2f} | {r['mass_mg']:.2f} mg ({r['ice_um']:.0f} µm) | "
                     f"{r['ceiling_J'] / 1e6:.1f} MJ | **{r['yield_record_J'] / 1e6:.1f} MJ** | "
                     f"{r['yield_typical_J'] / 1e6:.1f} MJ | {r['target_gain']:.1f} |")
    ratio = {k: by[k]["ceiling_J"] / by[base]["ceiling_J"] for k, base in
             (("B", "A"), ("C", "A"), ("D", "A"), ("E", "D"), ("F", "D"), ("G", "D"), ("H", "G"), ("I", "G"))}
    lines += ["",
              "What the model predicts, step by step:",
              (f"* **B: the static model says lowering the adiabat alone does nothing** (×{ratio['B']:.2f}), "
              "because it holds the pressure fixed. **NIF has already falsified this:** Hohenberger et al. (2025) "
              "lowered the adiabat from 2.89 to 2.68 on the same target and got up to 80% more yield. At fixed "
              "drive, a lower adiabat raises the pressure and areal density; the 1-D model (dynamics.md) shows it."),
              (f"* **C: pressure is the lever** (×{ratio['C']:.2f} for +15%). The model even needs less energy in the "
              f"fuel ({by['C']['laser_needed_J'] / 1e6:.2f} MJ of laser against {by['A']['laser_needed_J'] / 1e6:.2f}), "
              "because a smaller hot spot ignites at higher pressure. It doesn't price what the faster implosion costs."),
              f"* **D: the last 6% of laser energy into thicker ice** gives ×{ratio['D']:.2f}.",
              (f"* **E: once the fuel is re-sized, the adiabat matters** (×{ratio['E']:.2f} over D): cheaper "
              "compression leaves energy for more fuel."),
              ("* **F, H, I: pressure again dominates.** These rows assume the higher pressure costs no extra laser "
              "energy, so they are an upside."),
              (f"* **The upgrade.** At 2.6 MJ: {by['G']['yield_record_J'] / 1e6:.0f} MJ at today's pressure, "
              f"{by['H']['yield_record_J'] / 1e6:.0f} MJ with 700 Gbar, {by['I']['yield_record_J'] / 1e6:.0f} MJ "
              f"stretched. Published projections: {PUBLISHED_EYC}."), "",
              "## Energy", "",
              "With a lithium blanket and 40% of the heat turned into electricity:", "",
              "| Step | Target gain | Electricity out ÷ in with NIF's laser (~0.6%) | with a 10% plant laser |",
              "|---|---|---|---|"]
    for r in rows:
        lines.append(f"| {r['step']} | {r['target_gain']:.1f} | {r['electric_nif']:.3f} | "
                     f"{r['electric_plant_laser']:.2f} |")
    best = max(rows, key=lambda r: r["electric_plant_laser"])
    lines += ["",
              (f"No step makes net electricity. The best, {best['step'][:1]}, would return "
              f"{best['electric_plant_laser']:.0%} of its electricity with a plant-grade laser. The simulator's "
              "net-electricity designs (section 21 of docs/energy.md) need about ten times NIF's laser energy "
              "and areal density."), ""]
    (OUT / "summary.md").write_text("\n".join(lines) + "\n")
    data = {"calibration": cal, "later_shots": later, "coupling": coupling, "qualities": qualities,
            "campaign": rows, "blanket_heat_per_fusion": BLANKET}
    (OUT / "results.json").write_text(json.dumps(data, indent=1, default=float))
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
