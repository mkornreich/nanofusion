"""Pauli blocking: does squeezing the fuel until its electrons are degenerate switch off their drag?

Every route in docs/energy.md hits the same wall: electrons drag fast ions to
a stop before they fuse. In fuel squeezed until its electrons are degenerate,
an ion slower than the Fermi velocity can only kick electrons near the Fermi
surface, so the drag per metre stops growing with density while fusion keeps
growing. Three uses:

1. Crystal ions fired into cold squeezed fuel (beam-target).
2. Hot ions in cold degenerate electrons (Son & Fisch's proposal for p-11B).
3. Muon catalysis: a muon stuck to an alpha is knocked free more often.

Writes docs/results/pauli/summary.md and results.json.

    python3 tools/pauli_blocking.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import alternatives as A  # noqa: E402
from nanofusion.improved import degenerate as D  # noqa: E402
from nanofusion.improved import energy as E  # noqa: E402
from nanofusion.improved import physics as P  # noqa: E402

OUT = Path("docs/results/pauli")

# --- inputs from the literature and from earlier sections ------------------------------
SOLID_DT_G_CC = 0.25
# NIF's cold D-T shell: designed for ~1,000 g/cm3 (Edwards et al., LLNL-JRNL-433484 (2010)); 500-800 g/cm3
# measured in 2012 (Mackinnon et al., PRL 108, 215005); about 3 times the Fermi pressure (adiabat ~3) in the
# ignition-era designs (Kritcher et al., PRE 106, 025201 (2022)).
NIF_DENSITY_G_CC = 1000.0
ION_COUPLING = 0.15  # laser energy -> crystal explosion energy, best case (section 7)
IMPLOSION_EFFICIENCY = 0.01  # laser energy -> energy of the squeezed fuel, NIF-class implosions (~1%)
LASER_EFFICIENCY = 0.10  # electricity -> light, fusion laser designs
MUON_COST_GEV = 4.7  # optimised accelerator design (Kelly et al. 2021; section 9)
MUON_LIFETIME_S = 2.197e-6
BEST_MUON_SOURCE_PER_S = 1e8  # most intense muon beams (PSI, J-PARC), order of magnitude
ELECTRICITY_FUSIONS_PER_MUON = 1000.0  # section 9
# Slow protons in aluminium, a real degenerate electron gas (3 valence electrons per atom, r_s = 2.07).
AL_G_CC = 2.699
AL_VALENCE_M3 = 3 * AL_G_CC * 1e3 / (26.9815 * P.AMU)
PSTAR_AL = {10.0: 292.1, 25.0: 413.2, 50.0: 473.4, 100.0: 446.8}  # keV -> electronic stopping, MeV cm^2/g (NIST PSTAR)
METAL_DRAG = 1.6  # measured / linear-response friction on slow protons in aluminium (0.25 / 0.154 atomic units)


def sig2(x: float) -> str:
    """Two significant digits, no exponent: 338.2 -> 340, 0.0123 -> 0.012."""
    return f"{float(f'{x:.2g}'):,g}"


def percent(x: float) -> str:
    """A share as a percentage, keeping the digits that matter near 100%."""
    if x > 0.9995:
        return "≈100%"
    return f"{x:.1%}" if x > 0.99 else f"{x:.0%}"


def fmt(x: float) -> str:
    if x == 0:
        return "0"
    if math.isinf(x):
        return "∞"
    if 0.01 <= abs(x) < 1000:
        return f"{x:.2g}" if abs(x) < 100 else f"{x:.0f}"
    mantissa, exponent = f"{x:.1e}".split("e")
    mantissa = mantissa.rstrip("0").rstrip(".")
    sup = str(int(exponent)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"10{sup}" if mantissa == "1" else f"{mantissa}×10{sup}"


def density_text(rho: float) -> str:
    return f"{rho:,.0f}" if rho >= 10 else f"{rho:g}"


def energy_text(kev: float) -> str:
    return f"{sig2(kev / 1e3)} MeV" if kev >= 1000 else (f"{sig2(kev)} keV" if kev >= 1 else f"{sig2(kev * 1e3)} eV")


def joules(x: float) -> str:
    for unit, scale in (("GJ", 1e9), ("MJ", 1e6), ("kJ", 1e3), ("J", 1.0), ("mJ", 1e-3)):
        if x >= scale:
            return f"{sig2(x / scale)} {unit}"
    return f"{fmt(x)} J"


def seconds(x: float) -> str:
    for unit, scale in (("s", 1.0), ("ms", 1e-3), ("µs", 1e-6), ("ns", 1e-9), ("ps", 1e-12), ("fs", 1e-15)):
        if x >= scale:
            return f"{sig2(x / scale)} {unit}"
    return f"{fmt(x)} s"


def checks() -> tuple[list[str], dict]:
    lines = ["## Model checks", "",
             "**Its two textbook limits.** The electron stopping is Lindhard's random-phase approximation for a "
             "free electron gas at zero temperature, integrated numerically. The stopping number L (dE/dx = "
             "4π Z²e⁴ n_e L / m_e v²):", "",
             "| Electron density (× solid D-T) | L at v = 0.1 v_F | Lindhard–Winther friction | L at v = 20 v_F | "
             "Bethe, ln(2mv²/ħω_p) |", "|---|---|---|---|---|"]
    data = {"limits": []}
    for scale in (10, 1000, 1e5):
        n = 4.5e28 * scale
        g = D.fermi_gas(n)
        j = D.low_speed_coefficient(g["chi2"])
        slow = float(D.stopping_number(0.1, g["chi2"]))
        fast = float(D.stopping_number(20.0, g["chi2"]))
        bethe = math.log(2 * P.M_E * (20 * g["v_F"]) ** 2 / (D.HBAR * g["omega_p"]))
        lines.append(f"| {scale:g} | {slow:.3g} | {j * 1e-3:.3g} | {fast:.3g} | {bethe:.3g} |")
        data["limits"].append({"scale": scale, "slow": slow, "lindhard_winther": j * 1e-3, "fast": fast,
                               "bethe": bethe})
    lines += ["", "The low-speed coefficient also matches Lindhard & Winther's own table (Mat. Fys. Medd. Dan. Vid. "
              "Selsk. 34, no. 4 (1964), Table 1) to within 0.5%.", "",
              "**A real degenerate electron gas: slow protons in aluminium.** Aluminium's three valence "
              "electrons are close to a free gas, and their stopping of slow protons is measured, so this tests "
              "the Pauli-blocked regime directly.", "",
              "| Proton energy | v / v_F | This model (valence electrons only) | Measured (NIST PSTAR) | "
              "Model / measured |", "|---|---|---|---|---|"]
    g = D.fermi_gas(AL_VALENCE_M3)
    data["aluminium"] = []
    for kev, measured in PSTAR_AL.items():
        v = math.sqrt(2 * kev * P.J_PER_KEV / (P.MASS_U["p"] * P.AMU))
        model = float(D.electron_stopping(kev, P.MASS_U["p"], 1, AL_VALENCE_M3)) * 1e-5 / AL_G_CC
        lines.append(f"| {kev:g} keV | {v / g['v_F']:.2f} | {model:.0f} MeV cm²/g | {measured:.0f} MeV cm²/g | "
                     f"{model / measured:.2f} |")
        data["aluminium"].append({"keV": kev, "model": model, "measured": measured})
    lines += ["", "* **Slower than the electrons, the model drags too little.** Below the Fermi velocity it gives "
              f"two-thirds of the measured stopping. Measured slow-proton friction in metals is about {METAL_DRAG:g} "
              "times the linear-response value, because the electrons pile up around the proton (Montanari & "
              "Miraglia, [arXiv:1705.01575](https://arxiv.org/abs/1705.01575)). That pile-up weakens as the "
              "electrons get denser, and squeezed fuel is about six times denser in this sense (r_s ≈ 0.3 against "
              f"2.1). The tables below also give the results with the drag raised {METAL_DRAG:g}-fold.",
              "* **A hint from NIF.** Knock-on ions crossing NIF's cold, partly degenerate fuel showed a "
              "degeneracy-induced shift of their Bragg peak, at 90% confidence (Hayes et al., Nat. Phys. 16, 432 "
              "(2020), [OSTI](https://www.osti.gov/pages/servlets/purl/1881627)).", ""]
    solid = D.best_in_flight("D", "DT", SOLID_DT_G_CC)
    e0 = solid["e0_keV"]
    cold = float(E.thick_target(e0, "D", "DT", {"T": 2.25e28, "D": 2.25e28})["gain"][0])
    lines += ["", f"**At solid density** this free-electron model gives a deuteron in D-T a return of "
              f"{solid['gain']:.3f} at {e0:.0f} keV. The project's cold-matter model, with bound electrons, gives "
              f"{cold:.3f}. Below a few times solid density the electrons are bound, so the free-electron model is "
              "only a rough guide there. Above ~100 times solid density the fuel is pressure-ionised and the model "
              "applies.", ""]
    data["solid"] = {"free_electron": solid["gain"], "cold_model": cold, "e0_keV": e0}
    return lines, data


DENSITIES = (2.5, 25.0, 250.0, NIF_DENSITY_G_CC, 2500.0, 25000.0, 2.5e5)


def beams() -> tuple[list[str], dict]:
    energies = (100.0, 300.0, 1000.0, 3000.0)
    lines = ["## 1. Crystal ions fired into cold, squeezed fuel", "",
             "A deuteron slowing down in cold D-T squeezed to each density. The return is fusion energy (all "
             "17.6 MeV) / the deuteron's energy, as in sections 1–2 of docs/energy.md. Solid D-T is 0.25 g/cm³; "
             f"a NIF-class implosion reaches ~{NIF_DENSITY_G_CC:,.0f} g/cm³.", "",
             "| D-T density (g/cm³) | Fermi energy | Least work to squeeze, per electron | 100 keV | 300 keV | 1 MeV "
             "| 3 MeV | Best (deuteron energy) | Its energy to electrons | Its path |",
             "|---|---|---|---|---|---|---|---|---|---|"]
    rows = []
    for rho in DENSITIES:
        res = D.in_flight(energies, "D", "DT", rho)
        best = D.best_in_flight("D", "DT", rho)
        squeeze = D.squeeze_energy_kev(res["n_e"])
        lines.append(f"| {density_text(rho)} | {energy_text(res['fermi']['E_F_keV'])} | {energy_text(squeeze)} | "
                     + " | ".join(f"{g:.2g}" for g in res["gain"])
                     + f" | **{best['gain']:.2g}** ({best['e0_keV']:.0f} keV) | {best['electron_share']:.0%} | "
                     f"{best['range_m'] * 1e6:.2g} µm |")
        rows.append({"rho_g_cc": rho, "E_F_keV": res["fermi"]["E_F_keV"], "squeeze_keV": squeeze,
                     "gain": dict(zip(energies, res["gain"])), "best": best})
    crossing = _crossing(lambda rho: D.best_in_flight("D", "DT", rho)["gain"], 1.0, 250.0, 2500.0)
    D.ELECTRON_DRAG_FACTOR = METAL_DRAG
    try:
        stronger = {rho: D.best_in_flight("D", "DT", rho)["gain"] for rho in (NIF_DENSITY_G_CC, 2500.0, 25000.0)}
        crossing_strong = _crossing(lambda rho: D.best_in_flight("D", "DT", rho)["gain"], 1.0, 250.0, 25000.0)
    finally:
        D.ELECTRON_DRAG_FACTOR = 1.0
    lines += ["", f"* **The cold-catcher wall breaks.** In cold matter at normal density a deuteron returns at most "
              f"~0.4% of its energy (section 1). At the density NIF's fuel is designed for "
              f"({NIF_DENSITY_G_CC:,.0f} g/cm³) it returns {rows[3]['best']['gain']:.1f}. Above about "
              f"{sig2(crossing)} g/cm³ it returns **more than its own energy**, the first time in this project for "
              "an ion in unheated fuel.",
              "* **Why.** Below the Fermi velocity the drag per metre no longer grows with density. Fusion per metre "
              "does. At the highest densities the deuteron loses most of its energy to fuel ions instead, and that "
              f"caps the return near {rows[-1]['best']['gain']:.0f}.",
              f"* **With the drag raised {METAL_DRAG:g}-fold** (as measured in metals): "
              + ", ".join(f"{g:.2g} at {density_text(rho)} g/cm³" for rho, g in stronger.items())
              + f"; the crossing moves to about {sig2(crossing_strong)} g/cm³.",
              "* **Squeezing buys what heating buys.** Degenerate electrons drag on a slow ion about as much as a "
              "plasma at their Fermi temperature does (section 2's hot catcher), at a fifth of the energy and "
              "without radiating. That is the whole trick.", ""]
    lines += ["### Protons into squeezed boron", "",
              "| Boron density (g/cm³) | Fermi energy | Best return (proton energy) | Its energy to electrons |",
              "|---|---|---|---|"]
    pb = []
    for rho in (2.34, 234.0, 2340.0, 23400.0, 2.34e5):
        best = D.best_in_flight("p", "boron", rho, e_lo=100.0, e_hi=8000.0)
        lines.append(f"| {density_text(rho)} | {energy_text(best['E_F_keV'])} | {best['gain']:.2g} "
                     f"({best['e0_keV']:.0f} keV) | {best['electron_share']:.0%} |")
        pb.append({"rho_g_cc": rho, **best})
    lines += ["", f"* **Proton–boron stays below 1** even at {pb[-1]['rho_g_cc']:,.0f} g/cm³ "
              f"({pb[-1]['gain']:.2g}). Boron's charge of 5 makes collisions with the nuclei 25 times stronger, so "
              "they take over from the electrons as the drag.", ""]
    return lines, {"DT": rows, "crossing_g_cc": crossing, "pB": pb,
                   "metal_drag": {"factor": METAL_DRAG, "best": stronger, "crossing_g_cc": crossing_strong}}


def _crossing(f, target: float, lo: float, hi: float) -> float:
    """Density between lo and hi where f(rho) = target (f rising)."""
    for _ in range(30):
        mid = math.sqrt(lo * hi)
        lo, hi = (mid, hi) if f(mid) < target else (lo, mid)
        if hi / lo < 1.02:
            break
    return math.sqrt(lo * hi)


def costs() -> tuple[list[str], dict]:
    lines = ["## What the squeeze costs", "",
             "The trick lasts only while the fuel stays cold. Once the beam has added about as much heat as the "
             "squeeze put in, the fuel is a hot plasma (section 2), and the same ions would be better spent "
             "igniting it. So the fairest test of the trick alone gives the beam as much energy as the squeeze:", "",
             "* perfect drivers: return = best return / 2;",
             f"* real drivers: {ION_COUPLING:.0%} of the light into ions (section 7), and {IMPLOSION_EFFICIENCY:.0%} "
             "of the light into the squeezed fuel (a NIF-class implosion); then a "
             f"{LASER_EFFICIENCY:.0%}-efficient laser for the wall plug.", "",
             "| D-T density (g/cm³) | Best return per ion | Return, perfect drivers | Return on the light | "
             "From the wall plug | Hot spark that ignites the same fuel |", "|---|---|---|---|---|---|"]
    rows = []
    for rho in (250.0, NIF_DENSITY_G_CC, 2500.0, 25000.0, 2.5e5):
        g = D.best_in_flight("D", "DT", rho)["gain"]
        on_light = g / (1 / ION_COUPLING + 1 / IMPLOSION_EFFICIENCY)
        spark = min(A.ignition_ball(rho, t * 1e3)["energy_J"] for t in (8, 10, 12, 16, 20, 25, 32))
        lines.append(f"| {density_text(rho)} | {g:.2g} | {g / 2:.2g} | {on_light:.2g} | "
                     f"{on_light * LASER_EFFICIENCY:.2g} | {joules(spark)} |")
        rows.append({"rho_g_cc": rho, "gain": g, "perfect": g / 2, "on_light": on_light,
                     "wall_plug": on_light * LASER_EFFICIENCY, "spark_J": spark})
    lines += ["", "* **The squeeze eats the gain.** Even with perfect drivers the fuel returns only about its "
              "best return / 2. Real implosions put about 1% of the light into the fuel, so the light comes back "
              f"at {rows[1]['on_light']:.1%}–{rows[-1]['on_light']:.1%}, and the wall plug at a few parts per "
              "thousand at best.",
              "* **At these densities ignition is cheap.** A hot spark of joules to kilojoules ignites the whole "
              "squeezed fuel (section 11's cheapest hot ball). Ions are far better spent as that spark (fast "
              "ignition, section 7) than on in-flight fusion. The real problem is the squeeze itself: "
              "2–250 times NIF's density, without heating.", ""]
    return lines, {"rows": rows}


def hot_ions() -> tuple[list[str], dict]:
    lines = ["## 2. Hot ions in degenerate electrons (Son & Fisch)", "",
             "**Idea** (Son & Fisch, Phys. Lett. A 329, 76 (2004)). Proton–boron fails because its electrons drain "
             "the hot ions' energy and radiate it (section 4). Degenerate electrons drain less and radiate less. "
             "An ion much slower than the electrons loses energy at a fixed rate, whatever the density, while "
             "fusion grows with density. So above some density the ions' own fusion products should keep them "
             "hot. Son & Fisch found a possible self-burning regime above ~10⁵ g/cm³, with ions near 100 keV and "
             "electrons near 30 keV, and a gain of ~18 on the electrons' Fermi energy (40–200 in a 2006 "
             "follow-up, Phys. Lett. A 356, 72). They noted themselves that a power plant needs far more, and "
             "that squeezing fuel to such densities is an unsolved problem "
             "([preprint](https://bp-pub.pppl.gov/pub_report/2004/PPPL-4003.pdf)).", "",
             "**Test.** Steady power balance per unit volume, for proton–boron with 15% boron:",
             "* the drain: each ion loses energy at dE/dt = −(2η/M)(E − 3T_e/2). η is the Lindhard–Winther "
             "friction of degenerate electrons (with the relativistic mass correction), merging into Spitzer's "
             "once the electrons are hotter than their Fermi energy (generous by up to 2 near the join);",
             "* the heating: the share of the alphas' energy that goes to ions, from the stopping model above "
             "(target ions at rest, which is generous);",
             "* the electrons must radiate what they receive. Pauli blocking cuts their bremsstrahlung by about "
             "2.2 (T/T_F)^1.5 (Son & Fisch give the scaling; the prefactor is derived here). That fixes the "
             "electron temperature;",
             "* the ions keep themselves hot when charged fusion power exceeds the radiation.", "",
             "| Ion temperature | Lowest density where the ions keep themselves hot | × NIF's density | Fermi energy | "
             "Electron temperature (÷ Fermi energy) | Alpha energy to ions | Burn time | Fusion in one burn time ÷ "
             "energy to squeeze and heat |", "|---|---|---|---|---|---|---|---|"]
    rows = []
    for t in (50.0, 75.0, 100.0, 125.0, 150.0):
        rho = D.self_heating_density("pB", t)
        if math.isinf(rho):
            lines.append(f"| {t:g} keV | none below 10⁹ g/cm³ | | | | | | |")
            continue
        h = D.hot_ion_balance("pB", rho, t)
        burn = 1 / h["burn_rate_per_s"]
        gain = h["fusion_W_m3"] * burn / (h["squeeze_J_m3"] + h["ion_heat_J_m3"])
        lines.append(f"| {t:g} keV | {fmt(rho)} g/cm³ | {sig2(rho / NIF_DENSITY_G_CC)} | {energy_text(h['E_F_keV'])} | "
                     f"{energy_text(h['T_e_keV'])} ({h['T_e_over_T_F']:.2f}) | {h['alpha_to_ions']:.0%} | "
                     f"{seconds(burn)} | {gain:.2g} |")
        rows.append({**h, "ideal_gain_per_burn_time": gain})
    hot = {t: D.hot_ion_balance("pB", D.self_heating_density("pB", t), t) for t in (200.0, 300.0)}
    lo_rho = min(r["rho_g_cc"] for r in rows)
    hi_rho = max(r["rho_g_cc"] for r in rows)
    at_100 = next(r for r in rows if r["T_i_keV"] == 100.0)
    lines += ["", f"* **This reproduces Son & Fisch.** With 100 keV ions the model closes the balance at "
              f"{fmt(at_100['rho_g_cc'])} g/cm³ with electrons at {energy_text(at_100['T_e_keV'])}; they found "
              "> 10⁵ g/cm³ and ~30 keV.",
              f"* **Degeneracy does open a window, at white-dwarf density.** With 50–150 keV ions, proton–boron "
              f"keeps itself hot at {fmt(lo_rho)}–{fmt(hi_rho)} g/cm³, {sig2(lo_rho / NIF_DENSITY_G_CC)} to "
              f"{sig2(hi_rho / NIF_DENSITY_G_CC)} times NIF's design density. The fusion released in one burn time "
              f"is {min(r['ideal_gain_per_burn_time'] for r in rows):.0f}–"
              f"{max(r['ideal_gain_per_burn_time'] for r in rows):.1f} times the energy needed to squeeze and "
              "heat the fuel, before any driver losses. It is the same verdict as pycnonuclear D-T (section 10): it "
              "works on paper, but only in matter like a white dwarf's, and only after a 50–160-fold radial squeeze "
              "of a solid without heating it.",
              "* **The electrons are only partly degenerate there.** They must radiate what the ions drain into "
              "them, so they settle at 0.1–0.9 of their Fermi energy. Pauli blocking still cuts their radiation "
              "and their drag enough.",
              "* **Hotter ions don't need degeneracy.** With 200–300 keV ions the balance closes with the electrons "
              "at " + " and ".join(f"{energy_text(h['T_e_keV'])}" for h in hot.values()) + ", "
              + " and ".join(f"{h['T_e_over_T_F']:.0f}" for h in hot.values())
              + " times their Fermi energy: ordinary hot electrons, at any density. That is the narrow classical "
              "hot-ion window found for magnetically confined proton–boron (Liu et al. 2025, cited under \"Why "
              "everything hits the same wall\"). It needs confinement, not squeezing.",
              "* **D-T needs none of this.** Its alphas outrun its bremsstrahlung above about 4 keV at any density, "
              "and fuel that dense ignites from a joule-sized spark.", ""]
    return lines, {"rows": rows, "classical": {str(t): h for t, h in hot.items()}}


def muons() -> tuple[list[str], dict]:
    k = D.stripping_calibration()
    lines = ["## 3. Muon catalysis: freeing stuck muons", "",
             "**Idea.** Muon catalysis came closest of all cold routes (section 9). It is capped by sticking: in "
             "about 0.9% of fusions the muon sticks to the alpha. The stuck ion (μα)⁺ leaves with 3.46 MeV, and "
             "collisions with fuel nuclei knock the muon free while it slows down. At liquid density this "
             "happens about 35% of the time (the \"reactivation\"), leaving an effective sticking near the measured "
             "0.45–0.56%. Electrons slow the ion without freeing the muon: an electron moving with the ion carries "
             "about 0.5 keV, too little for the muon's 11 keV binding. If Pauli blocking switches off the "
             "electrons' drag, the ion meets many more nuclei before it stops.", "",
             "**Earlier work.** Hot plasma lowers the drag too: Jändel et al. (Phys. Rev. A 40, 2799 (1989)) found "
             "reactivation rising with temperature and density, but nearly complete only above ~1 keV (10⁷ K), "
             "thousands of times hotter than any fuel in which muon catalysis has been measured. Cripps, Harms & Goel (Hyperfine Interact. 77, 181 (1993)) found an "
             "energetically useful number of fusions per muon only below ~5 eV and above 100 times liquid density: "
             "the cold, dense corner tested here.", "",
             "**Model.** Stripping cross section: ionisation of hydrogen by protons (Rudd et al., Rev. Mod. Phys. "
             "57, 965 (1985)), scaled to the muonic ion (Born scaling, σ(v) = σ_H(v/2)(m_e/m_r)²/16). Slowing: "
             f"this section's stopping. A single factor (×{k:.2f}) is tuned once so that liquid D-T reproduces "
             "35%. The model leaves out the muon's excited states, which that factor absorbs.", "",
             "| D-T density (g/cm³) | Muons freed | Counting only the fast part (v > muon orbit) | Effective "
             "sticking | Fusions per muon allowed by sticking |", "|---|---|---|---|---|",
             f"| liquid, {D.liquid_density_g_cc():.2f} (calibration) | 35% | | "
             f"{D.STICKING_INITIAL * 0.65:.2%} | {1 / (D.STICKING_INITIAL * 0.65):.0f} |"]
    rows = []
    for rho in (2.5, 25.0, 50.0, 100.0, 250.0, NIF_DENSITY_G_CC):
        r = D.muon_reactivation(rho)
        fast = D.muon_reactivation(rho, fast_only=True)
        lines.append(f"| {density_text(rho)} | {percent(r['reactivation'])} | {percent(fast['reactivation'])} | "
                     f"{sig2(r['sticking'] * 100) if r['sticking'] > 1e-5 else fmt(r['sticking'] * 100)}% | "
                     f"{fmt(r['fusions_limit'])} |")
        rows.append({"rho_g_cc": rho, "reactivation": r["reactivation"], "fast_only": fast["reactivation"],
                     "sticking": r["sticking"], "fusions_limit": r["fusions_limit"]})
    enough = _crossing(lambda rho: D.muon_reactivation(rho)["fusions_limit"], ELECTRICITY_FUSIONS_PER_MUON, 2.5,
                       1000.0)
    # how long the squeeze lasts, and what a 1 ns squeeze would take
    rho = 100.0
    need_time = ELECTRICITY_FUSIONS_PER_MUON / D.DTMU_FUSION_RATE
    cc = A.cold_compression(rho)
    radius = 4 * cc["sound_speed_m_s"] * need_time
    mass_kg = 4 / 3 * math.pi * radius**3 * rho * 1e3
    electrons = mass_kg / (0.5 * (P.MASS_U["D"] + P.MASS_U["T"]) * P.AMU)
    squeeze_j = electrons * cc["electron_energy_keV"] * P.J_PER_KEV
    per_muon_j = ELECTRICITY_FUSIONS_PER_MUON * A.DT_FUSION_MEV * 1.602176634e-13
    muons_needed = squeeze_j / per_muon_j
    rate_needed = muons_needed / MUON_LIFETIME_S
    # with a real implosion (about 1% of the light reaches the fuel) the muons must repay 100 times more; their
    # alphas then heat the core (3/2 T per ion and per electron)
    real_muons = muons_needed / IMPLOSION_EFFICIENCY
    alpha_heat_j = real_muons * ELECTRICITY_FUSIONS_PER_MUON * 3.52e3 * P.J_PER_KEV
    heated_kev = alpha_heat_j / (3 * electrons * P.J_PER_KEV)
    nif_core = 50e-6
    lines += ["", f"* **Squeezing frees the stuck muons.** Above about {enough:,.0f} g/cm³ sticking falls below "
              f"0.1%, enough on paper for the ~{ELECTRICITY_FUSIONS_PER_MUON:,.0f} fusions per muon that "
              f"electricity needs (section 9). At {NIF_DENSITY_G_CC:,.0f} g/cm³ almost no muon stays stuck. Most "
              "of the stripping happens while the ion is still faster than the muon's orbit, where the cross "
              "section is most reliable.",
              "* **But the squeeze is too short.** A squeezed core holds together for R / 4c_s: "
              f"{seconds(D.squeeze_lifetime(rho, nif_core))} for a 50 µm core at {rho:g} g/cm³. Even if every "
              f"step of the cycle were as fast as fusion inside the d-t-μ molecule ({fmt(D.DTMU_FUSION_RATE)} per "
              f"second), a muon could catalyse only ~{D.DTMU_FUSION_RATE * D.squeeze_lifetime(rho, nif_core):.0f} "
              "fusions in that time, fewer than the 150 it already manages in liquid fuel. The cycle's other steps "
              "have never been measured in pressure-ionised fuel, and the molecular resonance that drives them at "
              "liquid density is gone there.",
              f"* **A squeeze long enough ({seconds(need_time)})** needs a core of {radius * 1e3:.2g} mm at "
              f"{rho:g} g/cm³: {mass_kg * 1e6:.2g} mg of D-T and {joules(squeeze_j)} even if squeezed perfectly cold. "
              f"To repay even that, the core must hold {fmt(muons_needed)} muons at once. They live "
              f"{MUON_LIFETIME_S * 1e6:.1f} µs, so the source would have to make {fmt(rate_needed)} per second, about "
              f"{fmt(rate_needed / BEST_MUON_SOURCE_PER_S)} times the most intense muon beams (~10⁸ per second). "
              f"A real implosion puts ~{IMPLOSION_EFFICIENCY:.0%} of the light into the fuel, so it needs "
              f"{1 / IMPLOSION_EFFICIENCY:.0f} times more, and their alphas would then heat the core to "
              f"~{energy_text(heated_kev)}: no longer cold fuel, and far outside anything muon catalysis has been "
              "measured in.", ""]
    return lines, {"rows": rows, "calibration": k, "density_for_1000": enough,
                   "long_squeeze": {"rho_g_cc": rho, "time_s": need_time, "radius_m": radius, "mass_kg": mass_kg,
                                    "squeeze_J": squeeze_j, "muons_needed": muons_needed,
                                    "rate_needed_per_s": rate_needed, "real_muons": real_muons,
                                    "alpha_heating_keV": heated_kev},
                   "nif_core_lifetime_s": D.squeeze_lifetime(rho, nif_core)}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    lines = ["# Pauli blocking: switching off the electrons' drag", "",
             "Produced by `tools/pauli_blocking.py`. Code: `nanofusion/improved/degenerate.py`.", "",
             "Every route in docs/energy.md hits the same wall: electrons drag fast ions to a stop before they "
             "fuse. In fuel squeezed until its electrons are degenerate, an ion slower than the Fermi velocity can "
             "only kick the few electrons near the Fermi surface, because every lower state is full. The drag per "
             "metre then stops growing with density (Fermi & Teller 1947), while fusion per metre keeps growing.",
             ""]
    data = {}
    for name, fn in (("checks", checks), ("beams", beams), ("costs", costs), ("hot_ions", hot_ions),
                     ("muons", muons)):
        part, d = fn()
        lines += part
        data[name] = d
    text = "\n".join(lines) + "\n"
    (OUT / "summary.md").write_text(text)
    (OUT / "results.json").write_text(json.dumps(data, indent=1, default=float))
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
