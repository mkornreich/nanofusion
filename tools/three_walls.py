"""Three more walls, three new attacks: trap the X-rays, hold the fuel with a current, steer through the lattice.

1. X-rays. Proton-boron fails because its bremsstrahlung carries off more
   than fusion makes (section 4). Make the fuel thick enough that the X-rays
   can't get out.
2. Holding the fuel without a laser. Run megaamperes through a fibre of solid
   fuel; the current's own field squeezes and holds it (a dense Z-pinch).
3. Electron drag, via the crystal itself. Ions channelling between the rows of
   atoms of a tritide crystal meet fewer electrons and more interstitial
   tritons.

Writes docs/results/walls/summary.md and results.json.

    python3 tools/three_walls.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import energy as E  # noqa: E402
from nanofusion.improved import walls as W  # noqa: E402

OUT = Path("docs/results/walls")
MEGATON_J = 4.184e15
NIF_DENSITY_G_CC = 1000.0  # design density of NIF's compressed fuel (section 17)
# NIF's best shots, N210808 and N221204: hot-spot rho-R 0.44 g/cm2 and 560-610 Gbar (PRL 129, 075001; PRL 132,
# 065102).
NIF_HOTSPOT_RHOR_G_CM2 = 0.44
NIF_PRESSURE_GBAR = 600.0
Z_CURRENT_A = 26e6  # Sandia's Z machine: ~26 MA from ~20 MJ stored (Stygar et al., PRSTAB 18, 110401 (2015))
Z_STORED_J = 20e6
# Fibre pinches: Sethian et al., PRL 59, 892 (1987): deuterium fibres 80-125 um across, up to 640 kA, stable
# through the 130 ns current rise (~100 growth times), unstable once the current peaked. Zap Energy's FuZE
# (Levitt et al., PRL 132, 155101 (2024)): ~250-370 kA, 1-3 keV, ~1e17 cm^-3, pinch radius 3-10 mm, stable for
# 1,000-10,000 growth times.
SETHIAN = {"current_A": 640e3, "radius_m": 50e-6}
FUZE = {"current_A": 300e3, "T_keV": 1.0, "n_m3": 1e23, "radius_m": 6e-3}
STRONGEST_LASER_FIELD_T = 800.0  # claimed laser-driven coil field (section 12)
# Channelling (flux peaking at the tritons, channelled / random stopping). Generous: a parallel beam at the
# channel centre and the ~1/2 stopping of fast ions. Measured: 1.4-2.2 flux peaking (Wampler & Myers 1999;
# simulations in diamond), stopping 5-13% below random for 50-150 keV protons in Si (arXiv:1907.08519).
GENEROUS_AXIAL = (3.0, 0.5)
MEASURED_AXIAL = (2.0, 0.9)
GENEROUS_PLANAR = (2.0, 0.9)
PLANAR_TO_AXIAL = 2.5  # planar critical angles are 2-3 times smaller than axial ones


def sig2(x: float) -> str:
    return f"{float(f'{x:.2g}'):,g}"


def fmt(x: float) -> str:
    if x == 0:
        return "0"
    if 0.01 <= abs(x) < 10000:
        return sig2(x)
    mantissa, exponent = f"{x:.1e}".split("e")
    mantissa = mantissa.rstrip("0").rstrip(".")
    sup = str(int(exponent)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"10{sup}" if mantissa == "1" else f"{mantissa}×10{sup}"


def length(m: float) -> str:
    for unit, scale in (("m", 1.0), ("cm", 1e-2), ("mm", 1e-3), ("µm", 1e-6), ("nm", 1e-9)):
        if m >= scale:
            return f"{sig2(m / scale)} {unit}"
    return f"{fmt(m)} m"


def mass(kg: float) -> str:
    for unit, scale in (("t", 1e3), ("kg", 1.0), ("g", 1e-3), ("mg", 1e-6), ("µg", 1e-9), ("ng", 1e-12)):
        if kg >= scale:
            return f"{sig2(kg / scale)} {unit}"
    return f"{fmt(kg)} kg"


def joules(x: float) -> str:
    for unit, scale in (("EJ", 1e18), ("PJ", 1e15), ("TJ", 1e12), ("GJ", 1e9), ("MJ", 1e6), ("kJ", 1e3), ("J", 1.0)):
        if x >= scale:
            return f"{sig2(x / scale)} {unit}"
    return f"{fmt(x)} J"


def seconds(x: float) -> str:
    for unit, scale in (("ns", 1e-9), ("ps", 1e-12), ("fs", 1e-15), ("as", 1e-18)):
        if x >= scale:
            return f"{sig2(x / scale)} {unit}"
    return f"{fmt(x)} s"


DENSITIES = (100.0, 1000.0, 1e4, 1e5)


def xrays() -> tuple[list[str], dict]:
    lines = ["## 1. Trap the X-rays: an opaque proton–boron ball", "",
             "**Idea.** Proton–boron makes no neutrons, but its electrons radiate X-rays (bremsstrahlung) faster "
             "than fusion heats the fuel (section 4). A big or dense enough ball is opaque to its own X-rays: they "
             "are reabsorbed inside and only leak out slowly from the surface. Does that save proton–boron?", "",
             "**Model.** A free ball of proton–boron fuel (15% boron) at a fixed density and ion temperature. The "
             "electrons settle where the energy they receive (from the hot ions, and their share of the alphas) "
             "equals what they lose: the X-rays that escape. Escape is followed frequency by frequency (free-free "
             "absorption and Thomson scattering, Osterbrock's escape probability for a sphere). The trapped "
             "frequencies fill with blackbody radiation, and photons near m_e c² make electron–positron pairs; "
             "the ball must pay for both. It \"repays its heat\" if its ions keep themselves hot and it releases at "
             "least its own energy before it flies apart (R / 4c_s).", "",
             "### The smallest balls that repay their heat", "",
             "| Density (g/cm³) | Ion temperature | Radius | Mass | Energy it holds | Areal density ρR | Pressure "
             "| X-rays trapped |", "|---|---|---|---|---|---|---|---|"]
    small = {}
    for rho in DENSITIES:
        b = W.smallest_burning_ball(rho)
        if b is None:
            lines.append(f"| {fmt(rho)} | none | | | | | | |")
            continue
        small[rho] = b
        lines.append(f"| {fmt(rho)} | {b['T_i_keV']:.0f} keV | {length(b['radius_m'])} | {mass(b['mass_kg'])} | "
                     f"{joules(b['energy_J'])} | {b['rhoR_g_cm2']:.2g} g/cm² | {fmt(b['pressure_Pa'] / 1e14)} Gbar | "
                     f"{b['trapped']:.0%} |")
    lines += ["", "### The smallest balls that trap at least half their X-rays (300 keV ions)", "",
              "| Density (g/cm³) | Radius | Mass | Energy it holds | Share of it in X-rays and pairs | Positrons "
              "per electron | Fusion released while it holds together ÷ its energy |", "|---|---|---|---|---|---|---|"]
    opaque = {}
    for rho in DENSITIES:
        o = W.smallest_opaque_ball(rho)
        if o is None:
            lines.append(f"| {fmt(rho)} | none below 1 km | | | | | |")
            continue
        opaque[rho] = o
        bomb = f" ({fmt(o['energy_J'] / MEGATON_J)} megatons)" if o["energy_J"] > 1e-3 * MEGATON_J else ""
        lines.append(f"| {fmt(rho)} | {length(o['radius_m'])} | {mass(o['mass_kg'])} | {joules(o['energy_J'])}{bomb} | "
                     f"{o['radiation_share']:.0%} | {fmt(o['positrons_per_electron'])} | {fmt(o['gain'])} |")
    nif = small.get(NIF_DENSITY_G_CC)
    nif_opaque = opaque.get(NIF_DENSITY_G_CC)
    rhor = [b["rhoR_g_cm2"] for b in small.values()]
    pressures = [b["pressure_Pa"] / 1e14 for b in small.values()]
    trapped = [b["trapped"] for b in small.values()]
    lines += ["",
              f"* **The balls that repay their heat are hot-ion balls at {min(rhor):.0f}–{max(rhor):.0f} g/cm².** "
              f"Ions at 200 keV, electrons near 90 keV. They trap {min(trapped):.0%}–{max(trapped):.0%} of their "
              "X-rays, more the denser they are. That areal density is "
              f"{min(rhor) / NIF_HOTSPOT_RHOR_G_CM2:.0f}–{max(rhor) / NIF_HOTSPOT_RHOR_G_CM2:.0f} times the hot "
              f"spot of NIF's best shots ({NIF_HOTSPOT_RHOR_G_CM2:g} g/cm²), at {fmt(min(pressures) / NIF_PRESSURE_GBAR)}"
              f" to {fmt(max(pressures) / NIF_PRESSURE_GBAR)} times their ~{NIF_PRESSURE_GBAR:.0f} Gbar.",
              "* **This matches the published burn studies.** Lavell et al. (Fusion Sci. Technol. 82 (2026), "
              "[OSTI](https://www.osti.gov/pages/servlets/purl/3366104)) find proton–boron burn propagating at "
              "7,000–9,000 g/cm³ with 14–18 g/cm², and failing at 1,000 g/cm³ with 2 g/cm². Ochs et al. "
              "([arXiv:2511.10885](https://arxiv.org/abs/2511.10885)) find breakeven needs areal densities about two "
              "orders of magnitude beyond the state of the art, at pressures three orders higher.",
              f"* **At NIF's design density** the smallest such ball holds {joules(nif['energy_J'])} "
              f"({mass(nif['mass_kg'])} of fuel)." if nif else "* No ball works at NIF's design density.",
              f"* **Trapping more than half the X-rays costs more than it saves at low density.** At NIF's density "
              f"that takes {mass(nif_opaque['mass_kg'])} of fuel holding {joules(nif_opaque['energy_J'])}, about "
              f"{nif_opaque['energy_J'] / MEGATON_J:.0f} megatons. The trapped X-rays and pairs hold "
              f"{nif_opaque['radiation_share']:.1%} of that, and it returns {fmt(nif_opaque['gain'])} of it." if
              nif_opaque else "",
              "* **Trapping helps only where proton–boron already needs white-dwarf-like matter**, 10⁴–10⁵ g/cm³, "
              "where a micrometre-scale ball is opaque without holding much radiation.", ""]
    return lines, {"smallest_burning": small, "smallest_opaque": opaque}


def pinch() -> tuple[list[str], dict]:
    held = W.alpha_confining_current()
    lines = ["## 2. Hold the fuel with a current: a fibre Z-pinch", "",
             "**Idea.** Section 12 held the fuel with a gold tamper and a magnetic field that something else had "
             "to make. A current does both at once: run megaamperes along a fibre of solid fuel, and the current's "
             "own magnetic field squeezes it and holds it in (a Z-pinch). No laser is needed.", "",
             "**Model.** Bennett's equilibrium: the field's pressure balances the plasma's, μ₀I²/8π = (N_i + N_e)kT, "
             "so each current holds a given amount of fuel per metre at 10 keV. The pinch lives some number of "
             "Alfvén times (radius ÷ Alfvén speed) before instabilities tear it. A narrower pinch is denser and "
             "burns faster, so there is a radius at which fusion repays the plasma's heat. Above "
             f"{held / 1e6:.2f} MA the pinch's own field holds the fusion alphas inside it.", ""]
    sethian_t = W.bennett_temperature(SETHIAN["current_A"], "D2", SETHIAN["radius_m"])
    fuze_n = 4e-7 * math.pi * FUZE["current_A"] ** 2 / (8 * math.pi * 2 * FUZE["T_keV"] * 1.602176634e-16)
    fuze_measured = FUZE["n_m3"] * math.pi * FUZE["radius_m"] ** 2
    lines += ["**Checks against real pinches.**",
              f"* A deuterium fibre 100 µm across at 640 kA (Sethian et al., PRL 59, 892 (1987)) is in Bennett "
              f"equilibrium at {sethian_t * 1e3:.0f} eV. Fibre pinches of that size reached a few hundred eV.",
              f"* At {FUZE['current_A'] / 1e3:.0f} kA and {FUZE['T_keV']:g} keV, Bennett holds {fmt(fuze_n)} ions per "
              f"metre. Zap Energy's FuZE pinch (Levitt et al., PRL 132, 155101 (2024)), ~10¹⁷ cm⁻³ over a 6 mm "
              f"radius, holds {fmt(fuze_measured)}.", "",
              "| Current | Solid D-T fibre it holds at 10 keV | Squeeze to repay the heat in one Alfvén time "
              "(= Alfvén times needed with no squeeze) | That pinch lives | Squeeze for the alphas alone to repay it "
              "| Field energy ÷ plasma heat |", "|---|---|---|---|---|---|"]
    rows = []
    for current in (0.64e6, 1.4e6, 5e6, Z_CURRENT_A, 60e6):
        one = W.fibre_pinch(current)
        alpha = W.fibre_pinch(current, fusion_over_heat=5.0)
        squeeze = f"×{sig2(one['squeeze'])} ({fmt(one['rho_g_cc'])} g/cm³)" if one["squeeze"] >= 1 else "none needed"
        lines.append(f"| {current / 1e6:g} MA | {length(2 * one['fibre_radius_m'])} across | {squeeze} | "
                     f"{seconds(one['lifetime_s'])} | ×{sig2(alpha['squeeze'])} ({fmt(alpha['rho_g_cc'])} g/cm³) | "
                     f"{one['field_J_m'] / one['heat_J_m']:.1f} |")
        rows.append({"one": one, "alpha": alpha})
    lid = W.fibre_pinch(1e6, fuel="6LiD")
    z = rows[3]
    field_per_cm = z["one"]["field_J_m"] / 100
    lines += ["",
              f"* **At the Z machine's {Z_CURRENT_A / 1e6:g} MA the numbers are not crazy.** A D-T fibre "
              f"{length(2 * z['one']['fibre_radius_m'])} across repays its heat if squeezed "
              f"{z['one']['squeeze']:.0f}-fold and held for {seconds(z['one']['lifetime_s'])}, one Alfvén time; "
              f"or held unsqueezed for {z['one']['squeeze']:.0f} Alfvén times.",
              "* **Fibre pinches haven't lasted.** Sethian's deuterium fibres stayed stable for ~100 growth times, "
              "but only while the current rose and only at a few hundred eV; at peak current they went unstable. "
              "Other fibre pinches (Los Alamos, Imperial College) expanded and broke up early, and their neutrons "
              "came mostly from beams, not heat. Zap Energy's sheared-flow pinch lasts 1,000–10,000 growth times, "
              "but at a millionth of solid density and 1–3 keV.",
              f"* **The field costs more than the plasma.** At 26 MA the field around the pinch (out to a return "
              f"conductor 1 cm away) holds {z['one']['field_J_m'] / z['one']['heat_J_m']:.1f} times the plasma's "
              f"heat, {joules(field_per_cm)} per cm, out of the ~{Z_STORED_J / 1e6:.0f} MJ that Z stores. So even "
              "ignition-like conditions are far from net electricity unless the burn multiplies the energy many "
              "times over.",
              f"* **Lithium deuteride fibres are hopeless.** D + ⁶Li burns thousands of times more slowly than D-T: "
              f"a 1 MA LiD fibre would need a {fmt(lid['squeeze'])}-fold squeeze. No LiD or D-T fibre pinch seems "
              "to have been tried.", ""]
    return lines, {"alpha_current_A": held, "rows": rows, "LiD_1MA": lid,
                   "checks": {"sethian_T_keV": sethian_t, "fuze_bennett_ions_per_m": fuze_n,
                              "fuze_measured_ions_per_m": fuze_measured}}


def channels() -> tuple[list[str], dict]:
    # Titanium tritide has the fluorite structure of TiH2 (fcc titanium, a = 4.42-4.45 A, hydrogen in tetrahedral
    # sites). The closest-packed rows, <110>, have Ti every a / sqrt 2 = 3.15 A and the widest critical angle;
    # using it for every axis is generous. Planar angles are 2-3 times smaller (Andersen's lecture notes).
    spacing = 4.45e-10 / math.sqrt(2)
    lines = ["## 3. Steer through the lattice: channelling in a tritide crystal", "",
             "**Idea.** A fast ion entering a crystal almost along a row or plane of atoms is steered gently "
             "between them (channelling). In the open channel it meets fewer electrons, so it loses energy more "
             "slowly, and its flux piles up in the channel centre. In titanium tritide the tritons sit in the "
             "channel centres along ⟨100⟩ (but are hidden behind the rows along ⟨111⟩). Does the crystal itself "
             "beat the cold-catcher limit?", "",
             "| Deuteron energy | Axial critical angle (Lindhard) | Share of an isotropic explosion along an axis "
             "(26 axes) | ... along a plane (13 planes) |", "|---|---|---|---|"]
    angles = {}
    for kev in (100.0, 300.0, 1000.0):
        psi = W.critical_angle(kev, 1, 22, spacing)
        axial = W.axial_acceptance(psi, 26)
        planar = W.planar_acceptance(psi / PLANAR_TO_AXIAL, 13)
        angles[kev] = {"psi_deg": math.degrees(psi), "axial": axial, "planar": planar}
        lines.append(f"| {kev:g} keV | {math.degrees(psi):.2g}° | {axial:.2%} | {planar:.0%} |")
    random = float(E.thick_target(300.0, "D", "TiD2", {"T": 1e29})["gain"][0])
    generous = W.channelled_return(random, [(1.0, *GENEROUS_AXIAL, 1.0)])
    measured = W.channelled_return(random, [(1.0, *MEASURED_AXIAL, 1.0)])
    a = angles[300.0]
    isotropic = W.channelled_return(random, [(a["axial"], *GENEROUS_AXIAL, 0.5), (a["planar"], *GENEROUS_PLANAR, 0.5)])
    field = W.magnetised_electron_field(4.5e28, 1.0)
    lines += ["", f"A 300 keV deuteron stopping in random (non-channelling) titanium tritide returns {random:.4f} of "
              "its energy. With channelling:", "",
              "| Case | Flux at the tritons | Stopping ÷ random | Return |", "|---|---|---|---|",
              f"| Perfectly aligned beam, generous | ×{GENEROUS_AXIAL[0]:g} | {GENEROUS_AXIAL[1]:g} | {generous:.3f} |",
              f"| Perfectly aligned beam, measured values | ×{MEASURED_AXIAL[0]:g} | {MEASURED_AXIAL[1]:g} | "
              f"{measured:.4f} |",
              f"| Exploding crystal (isotropic), generous, half the path channelled | axes ×{GENEROUS_AXIAL[0]:g}, "
              f"planes ×{GENEROUS_PLANAR[0]:g} | {GENEROUS_AXIAL[1]:g}, {GENEROUS_PLANAR[1]:g} | {isotropic:.4f} |",
              "",
              "Measured inputs: flux peaking of 1.4–2.2 at interstitial sites (deuterium in GaN; simulated in "
              "diamond), and channelled stopping only 5–13% below random for 50–150 keV protons in silicon, the "
              "speeds of 100–300 keV deuterons. Stopping halves only above ~1 MeV per nucleon.", "",
              f"* **A perfectly aligned beam does {measured / random:.1f}–{generous / random:.0f} times better**, "
              f"{measured:.4f}–{generous:.3f}. Still {1 / generous:.0f}–{1 / measured:.0f} times short of breakeven. Experiments agree on the "
              "size: 20 keV deuterons on palladium gave about twice the neutrons along a crystal plane as at 30° "
              "to it (Dalkarov et al., Phys. Rev. Accel. Beams 22, 034201 (2019)).",
              f"* **An exploding crystal can't aim.** {a['axial']:.1%} of its ions fall inside an axial critical "
              f"angle and about {a['planar']:.0%} along a plane (generous), and they dechannel partway: "
              f"{isotropic:.4f}, a {isotropic / random - 1:.0%} gain even with generous numbers.",
              f"* **The other way to blind the electrons, a magnetic field, needs ~{fmt(field)} T** in solid D-T "
              "(electrons must gyrate inside a Debye length, i.e. their cyclotron frequency must exceed the plasma "
              f"frequency), about {field / STRONGEST_LASER_FIELD_T:.0f} times the strongest laser-made field claimed.",
              ""]
    return lines, {"angles": angles, "random": random, "aligned_generous": generous, "aligned_measured": measured,
                   "isotropic": isotropic, "magnetised_field_T": field}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    lines = ["# Three more walls: trap the X-rays, hold the fuel with a current, steer through the lattice", "",
             "Produced by `tools/three_walls.py`. Code: `nanofusion/improved/walls.py`.", ""]
    data = {}
    for name, fn in (("xrays", xrays), ("pinch", pinch), ("channels", channels)):
        part, d = fn()
        lines += part
        data[name] = d
    text = "\n".join(line for line in lines) + "\n"
    (OUT / "summary.md").write_text(text)
    (OUT / "results.json").write_text(json.dumps(data, indent=1, default=float))
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
