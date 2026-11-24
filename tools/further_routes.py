"""More routes: a fission blanket, a D-T chain reaction, bare ions in a trap, crystals as projectiles.

1. Fission blanket (fusion-fission hybrid). Every D-T fusion gives a neutron,
   and a neutron in a subcritical uranium or thorium blanket starts fission
   chains worth far more than the fusion itself. How good must the fusion
   source be before the blanket turns it into net electricity?
2. Chain reaction in D-T. Section 4 of docs/energy.md found no avalanche in
   proton-boron fuel. D-T is a better case: its neutron carries 14 MeV and
   knocks on deuterons and tritons that fuse 100 times more easily.
3. No electrons at all. Electrons do the damage in every route so far: they
   stop fast ions and radiate. A cloud of bare ions in a magnetic trap has
   neither problem. What limits it?
4. Crystals as projectiles (cluster-impact fusion). Accelerate whole crystals
   with a voltage instead of exploding them with a laser. Accelerators are
   several times more efficient than lasers.
5. A check of the model against nanowire arrays, the best nanostructured
   fusion targets so far.

Writes docs/results/further/summary.md.

    python3 tools/further_routes.py
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import alternatives as A  # noqa: E402
from nanofusion.improved import energy as E  # noqa: E402
from nanofusion.improved import physics as P  # noqa: E402

OUT = Path("docs/results/further")

# --- inputs from the literature and from earlier sections ------------------------------
LASER_EFFICIENCY = 0.10  # electricity -> light, diode-pumped fusion laser designs (LIFE: 10-15%)
ACCELERATOR_EFFICIENCY = 0.40  # electricity -> beam assumed in Rubbia's energy amplifier (best existing: PSI, 18%)
NIF_LASER_EFFICIENCY = 2.08 / 400.0  # NIF: 300-400 MJ of electricity per 2.08 MJ shot
COUPLING = 0.15  # share of the light assumed to reach the fuel in the crystal-ball designs
K_LIMIT = 0.98  # highest k_eff of accelerator-driven reactor designs (Rubbia's energy amplifier; MYRRHA uses 0.95)
SPALLATION_NEUTRONS_PER_GEV = 20.0  # thick lead target at 1 GeV (Bauer, PSI Bericht 98-06 (1998))
MEASURED_PILE_MEV = 290.0  # heat per 14 MeV neutron in a natural-uranium pile: 1.46 fissions (McEnhill & Weale 1962)
HYBRID_SOURCES = [
    # name, fusion energy / driver energy, driver efficiency
    ("30 nm crystals on a tritide catcher, 6 fs laser (experiments.md)", 3e-5, LASER_EFFICIENCY),
    ("the same if all the light became ion energy (section 1)", 0.004, LASER_EFFICIENCY),
    ("hollow shell implosion (section 8)", 3e-4, LASER_EFFICIENCY),
    ("ball of crystals at its breakeven, bare (section 3) or tamped and magnetised (section 12)", COUPLING,
     LASER_EFFICIENCY),
    ("tamped, magnetised ball returning 10 × its fuel energy (section 12)", 10 * COUPLING, LASER_EFFICIENCY),
    ("muon catalysis at the record 150 fusions per muon (section 9)", 150 * 17.59 / 5000.0, ACCELERATOR_EFFICIENCY),
    ("NIF's best shot, with NIF's own laser", 4.1, NIF_LASER_EFFICIENCY),
    ("NIF's best shot, with a 10% efficient laser", 4.1, LASER_EFFICIENCY),
]
STRONGEST_STEADY_FIELD_T = 45.5  # hybrid magnet, National High Magnetic Field Laboratory
RECORD_LUMINOSITY_CM2_S = 5.1e34  # SuperKEKB, 2024
# Beuhler, Friedlander & Friedman, PRL 63, 1292 (1989): (D2O)N+, N = 25-1300, 200-325 keV, on TiD. Their rate
# for N = 100 at 300 keV, as normalised by Carraro et al. (PRA 42, 1379 (1990)), is 1e-12 to 1e-11 per cluster;
# Carraro et al.'s own thick-target estimate is 4e-39.
BEUHLER = {"kev": 300.0, "sizes": (25, 100, 1300), "claimed_per_cluster": (1e-12, 1e-11), "carraro_n100": 4e-39}
# Curtis et al., Phys. Rev. Research 3, 043181 (2021): deuterated nanowire array, 8 J, 45 fs, 3e21 W/cm2.
# Their particle-in-cell run gives 3.8e12 deuterons above 0.2 MeV and 3e11 above 1 MeV; measured 1.2e7 neutrons.
NANOWIRE = {"laser_J": 8.0, "above_200keV": 3.8e12, "above_1MeV": 3.0e11, "max_keV": 13000.0, "neutrons": 1.2e7,
            "record_neutrons_per_J": 2.2e6}  # the record per joule: Curtis et al., Nat. Commun. 9, 1077 (2018)


def fmt(x: float) -> str:
    if x == 0:
        return "0"
    if 0.01 <= abs(x) < 1000:
        return f"{x:.2g}" if abs(x) < 100 else f"{x:.0f}"
    mantissa, exponent = f"{x:.1e}".split("e")
    mantissa = mantissa.rstrip("0").rstrip(".")
    sup = str(int(exponent)).translate(str.maketrans("-0123456789", "⁻⁰¹²³⁴⁵⁶⁷⁸⁹"))
    return f"10{sup}" if mantissa == "1" else f"{mantissa}×10{sup}"


def round_2(x: float) -> str:
    """Two significant digits with thousands separators: 58633 -> 59,000."""
    return f"{float(f'{x:.2g}'):,.0f}"


def k_text(k: float | None) -> str:
    if k is None:
        return "no uranium needed"
    if k == 0:
        return "0 (uranium, no chain)"
    return f"{k:.2f}" if k < 0.99 else (f"{k:.4f}" if k < 0.9999 else "> 0.9999")


def hybrid() -> tuple[list[str], dict]:
    ks = (0, 0.5, 0.9, 0.95, 0.98)
    lithium = A.blanket_energy_per_fusion_mev(0.0, uranium=False)
    lines = ["## 1. A fission blanket around the fusion source", "",
             "Each D-T fusion sends a 14 MeV neutron into a blanket. In lithium alone it breeds the next triton, for "
             f"{lithium:.0f} MeV per fusion in all. In uranium a 14 MeV neutron does much more: it splits "
             f"uranium-238 directly (about {A.FAST_FISSIONS} fissions) and multiplies to about {A.NEUTRON_WORTH:g} "
             "neutrons. If the blanket is kept just below critical (k_eff < 1), each of those starts chains worth "
             "k / (2.5 (1−k)) fissions at 200 MeV. Heat per fusion:", "",
             "| uranium blanket k_eff | 0 (depleted uranium, no chain) | 0.5 | 0.9 | 0.95 | 0.98 |",
             "|---|---|---|---|---|---|",
             "| heat per fusion | " + " | ".join(f"{A.blanket_energy_per_fusion_mev(k):,.0f} MeV" for k in ks) + " |",
             "| × the fusion energy | "
             + " | ".join(f"{A.blanket_energy_per_fusion_mev(k) / A.DT_FUSION_MEV:,.0f}" for k in ks) + " |",
             "",
             f"Check: a natural-uranium pile measured {MEASURED_PILE_MEV:.0f} MeV per 14 MeV neutron (McEnhill & "
             f"Weale 1962); the k_eff = 0 column gives {A.blanket_energy_per_fusion_mev(0.0) - A.DT_FUSION_MEV + 14.07:.0f}. "
             f"Accelerator-driven reactor designs stop at k_eff = 0.95 (MYRRHA) to {K_LIMIT} (Rubbia's energy "
             "amplifier) to stay safely below critical. Net electricity needs (driver efficiency) × (fusion / driver "
             "energy) × (blanket multiplication) × (40% heat to electricity) ≥ 1. Lasers are taken as 10% efficient "
             "and accelerators as 40%. Both are design figures: NIF's laser manages 0.5% and the best proton "
             "accelerator 18%.", "",
             "| fusion source | fusion / driver energy | electricity per neutron | k_eff for net electricity "
             "| k_eff for 4 × the electricity put in | verdict |",
             "|---|---|---|---|---|---|"]
    rows = []
    for name, gain, eff in HYBRID_SOURCES:
        k1 = A.hybrid_k_needed(gain, eff, 1.0)
        k4 = A.hybrid_k_needed(gain, eff, 4.0)
        per_neutron = A.DT_FUSION_MEV / gain / eff
        if k4 is None or k4 <= K_LIMIT:
            verdict = "works on paper"
        elif k1 is None or k1 <= K_LIMIT:
            verdict = "net electricity on paper, but most of it recirculates"
        else:
            verdict = "no: the blanket would have to be a critical reactor"
        rows.append({"source": name, "gain": gain, "efficiency": eff, "k_net": k1, "k_plant4": k4,
                     "MeV_electric_per_neutron": per_neutron,
                     "fission_share": 0.0 if k1 is None else A.fission_share(k1)})
        lines.append(f"| {name} | {fmt(gain)} | {fmt(per_neutron)} MeV | {k_text(k1)} | {k_text(k4)} | {verdict} |")
    spallation = 1000.0 / SPALLATION_NEUTRONS_PER_GEV / ACCELERATOR_EFFICIENCY
    worst = [r for r in rows if r["k_net"] is not None and r["k_net"] > K_LIMIT]
    lines += ["", f"* **The benchmark is spallation.** A 1 GeV proton makes about {SPALLATION_NEUTRONS_PER_GEV:.0f} "
              f"neutrons in lead: about {spallation:.0f} MeV of electricity per neutron. A fusion source is worth "
              "putting inside a blanket only if its neutrons cost about that, or a few times more, since a 14 MeV "
              "neutron is worth about three ordinary ones.",
              f"* **Exploding crystals are {round_2(min(r['MeV_electric_per_neutron'] for r in worst) / spallation)} to "
              f"{round_2(max(r['MeV_electric_per_neutron'] for r in worst) / spallation)} times too expensive a neutron "
              f"source.** They would need k_eff of {min(r['k_net'] for r in worst):.4f} or more. The margin to a "
              "runaway chain is the delayed-neutron fraction, 0.2–0.7%, so that blanket is a critical reactor in all "
              "but name and the crystals add nothing.",
              "* **Any source that returns a tenth or more of its driver energy works on paper.** Muon catalysis "
              "needs only plain uranium. But then at least "
              f"{A.fission_share(0.0):.0%} of the heat is fission, with fission's fuel cycle and waste. It is a "
              "fission reactor with an expensive neutron source.", ""]
    return lines, {"rows": rows, "spallation_MeV_electric_per_neutron": spallation}


def chain() -> tuple[list[str], dict]:
    mfp = 2.515 * 1.66053906660e-24 / (0.5 * (E.N_ELASTIC_14MEV_B["D"] + E.N_ELASTIC_14MEV_B["T"]) * 1e-24)
    lines = ["## 2. A chain reaction in D-T, carried by its own neutrons and alphas", "",
             "k = D-T fusions triggered per D-T fusion. The 14.1 MeV neutron knocks on deuterons and tritons "
             "(up to 12.5 and 10.6 MeV), and so does the 3.5 MeV alpha; each knock-on may fuse before it stops. "
             "Infinite medium: the neutron is followed through all its collisions, which are peaked forwards as "
             "measured. Alpha knock-ons use Rutherford scattering, × 5 above 0.5 MeV to allow for nuclear "
             "scattering. Two versions of the D-T cross section above 550 keV: "
             "held constant (as elsewhere in this project, an overestimate) and the measured fall-off "
             "(Bosch & Hale's high-energy fit).", "",
             "| medium (infinite) | k, cross section held constant | k, measured fall-off | of which from the alpha "
             "| chance that a 5 MeV deuteron fuses |", "|---|---|---|---|---|"]
    data = {}
    media = [("cold solid D-T", None)]
    media += [(f"solid-density D-T plasma at {t:g} keV", t) for t in (1.0, 3.0, 10.0, 30.0, 100.0)]
    try:
        for name, t in media:
            plasma = None if t is None else E.fuel_plasma("DT", 1.0, t, t)
            row = {}
            for label, high in (("constant", False), ("measured", True)):
                E.DT_HIGH_ENERGY = high
                row[label] = E.dt_chain(hot_plasma=plasma, nuclear_boost=5.0)
            data[name] = row
            m = row["measured"]
            lines.append(f"| {name} | {fmt(row['constant']['k'])} | {fmt(m['k'])} | {fmt(m['k_alpha'])} "
                         f"| {fmt(m['yield_at_5MeV']['D'])} |")
        # a body of finite size: only the first neutron collision counts
        E.DT_HIGH_ENERGY = True
        lines += ["", f"A 14 MeV neutron travels about {mfp:.1f} g/cm² between elastic collisions, which is "
                  f"{mfp / 0.19:.0f} cm of solid D-T (4.7 g/cm² counting all reactions). In a body thinner than "
                  f"that, only a share ρR / {mfp:.1f} "
                  "of the neutrons collide at all. The alpha's part counts only if the alphas and their knock-ons "
                  "stay inside, which needs a field (section 12):", "",
                  "| body | ρR | k from the neutron | k from the alpha, if held in |", "|---|---|---|---|"]
        finite = [("tamped, magnetised breakeven ball of section 12 (20 keV)", 0.005, 20.0),
                  ("the 100-fold-return ball of section 12 (30 keV while burning)", 0.07, 30.0),
                  ("a ball of solid D-T one mean free path in radius (30 keV)", mfp, 30.0)]
        for name, rho_r, t in finite:
            plasma = E.fuel_plasma("DT", 1.0, t, t)
            res = E.dt_chain(hot_plasma=plasma, nuclear_boost=5.0,
                             first_collision_probability=1 - math.exp(-rho_r / mfp))
            data[name] = res
            lines.append(f"| {name} | {rho_r:.2g} g/cm² | {fmt(res['k_neutron'])} | {fmt(res['k_alpha'])} |")
    finally:
        E.DT_HIGH_ENERGY = False
    lines += ["", "* **Cold D-T: k ≈ 0.002.** No chain reaction, by a factor of 600. Electrons stop the knock-ons "
              "first, as in every cold route.",
              "* **Hot D-T: k reaches about 1, but only in a plasma at 30–100 keV that is metres across.** That "
              "is section 2's result again: a fast ion in hot D-T returns 1–3 times its energy, and the neutron "
              "hands its 14 MeV to fast ions only if the fuel is many mean free paths thick. Such a plasma is "
              "already burning. The chain adds to a fire; it cannot light one.",
              "* **In any ball this project can afford, k is a few per cent**, almost all from the alpha, and "
              "only if the alphas are held in. That is a few per cent of extra yield.",
              "* **Published work is less generous still.** Peres & Shvarts (Nucl. Fusion 15, 687 (1975)) found "
              "cold D-T critical only above 3.5×10⁴ g/cm³, and Afek et al. (J. Phys. D 11, 2171 (1978)) found "
              "4,200 g/cm³ at 14–15 keV. A 2026 Monte Carlo study (arXiv:2604.06769) finds D-T critical only "
              "without neutron leakage and above about 200 ρ^-0.23 keV, which is 290 keV at solid density.", ""]
    return lines, {"media": data, "neutron_mfp_g_cm2": mfp}


def trap() -> tuple[list[str], dict]:
    lines = ["## 3. No electrons at all: bare ions in a magnetic trap", "",
             "A cloud of bare deuterons and tritons has no electrons to stop fast ions or to radiate, and "
             "collisions between like charges do not move them across a magnetic field. But the cloud's own "
             "charge limits its density to the Brillouin limit, n = ε₀B² / 2m (Brillouin 1945; Dubin & O'Neil, "
             "Rev. Mod. Phys. 71, 87 (1999)). Column 5 cm in radius and 1 m "
             "long, ions at 60 keV (where D-T reacts fastest):", "",
             "| field | ion density | fusion power per m³ | fusion power of the column | voltage across the column "
             "| time the field must last to repay its own energy |", "|---|---|---|---|---|---|"]
    rows = []
    for b, label in ((10.0, "10 T (ordinary superconducting magnet)"), (20.0, "20 T (best superconducting magnets)"),
                     (STRONGEST_STEADY_FIELD_T, f"{STRONGEST_STEADY_FIELD_T:g} T (strongest steady magnet)"),
                     (1e3, "1 kT (laser-driven coil, lasts ~1 ns)"),
                     (1e5, "100 kT (inside a laser plasma, lasts ~1 ps)")):
        r = A.ion_trap(b, 0.05, 1.0)
        rows.append(r)
        payback = r["field_payback_s"]
        payback_text = (f"{payback / 86400:.0f} days" if payback > 86400 else
                        (f"{payback / 60:.0f} minutes" if payback > 60 else f"{fmt(payback)} s"))
        lines.append(f"| {label} | {fmt(r['density_m3'] * 1e-6)} cm⁻³ | {fmt(r['power_density_W_m3'])} W "
                     f"| {fmt(r['power_W'])} W | {fmt(r['voltage_V'] / 1e6)} MV | {payback_text} |")
    watt = 1.0 / (5.0e-24 * 17.59 * 1.602176634e-13)
    lines += ["", "* **A steady trap makes watts, not megawatts.** At 10 T the limit is 10¹¹ ions per cm³, a few "
              "billionths of air's density, and 7 W per m³. A tokamak reactor needs ~10⁶ W/m³. The power "
              "grows as B⁴, but so does the difficulty: the cloud's own voltage is megavolts at 10 T.",
              "* **Laser-made fields are strong enough but far too brief.** A 100 kT field would have to last "
              "~0.1 s to repay its energy. It lasts a picosecond.",
              "* **The same limit capped the imploding shell of section 8.** There the ions' repulsion limited "
              "how many fit in the focus; here it limits how many fit in the magnet.",
              f"* **Colliding beams** are the same idea in a ring. 1 W of D-T fusion needs a luminosity of "
              f"{fmt(watt)} cm⁻² s⁻¹ at the 5-barn peak. The record for any collider is "
              f"{fmt(RECORD_LUMINOSITY_CM2_S)} cm⁻² s⁻¹ (SuperKEKB), which would give "
              f"{RECORD_LUMINOSITY_CM2_S / watt:.1f} W.", ""]
    return lines, {"rows": rows, "luminosity_for_1W": watt}


def projectiles() -> tuple[list[str], dict]:
    lines = ["## 4. Crystals as projectiles (cluster-impact fusion)", "",
             "**The original claim.** Beuhler, Friedlander & Friedman (Phys. Rev. Lett. 63, 1292 (1989)) fired "
             "singly charged heavy-water clusters, (D₂O)ₙ with n = 25–1300, at 200–325 keV into titanium "
             "deuteride, and reported about 10⁻¹² to 10⁻¹¹ D-D fusions per cluster of 100 molecules. A cluster's "
             f"energy is shared by all its atoms, so each deuteron has very little. At {BEUHLER['kev']:.0f} keV:", "",
             "| cluster | energy per deuteron | fusions per cluster, this model |", "|---|---|---|"]
    claim = []
    for n in BEUHLER["sizes"]:
        r = A.cluster_impact_yield(n, BEUHLER["kev"])
        claim.append({"molecules": n, **r})
        lines.append(f"| (D₂O){n} | {r['keV_per_deuteron'] * 1e3:.0f} eV | {fmt(r['fusions_per_cluster'])} |")
    n100 = next(c for c in claim if c["molecules"] == 100)["fusions_per_cluster"]
    orders = math.log10(BEUHLER["claimed_per_cluster"][0] / n100)
    lines += ["", f"For 100 molecules this model gives {fmt(n100)} per cluster, and an independent estimate at the "
              f"time gave {fmt(BEUHLER['carraro_n100'])} (Carraro et al., Phys. Rev. A 42, 1379 (1990)). The claim "
              f"was about {orders:.0f} orders of magnitude higher. The authors withdrew it in 1992 (Phys. Rev. "
              "Lett. 68, 2108): light ions contaminating the beam, each carrying the full energy, had probably "
              "produced the fusion.", "",
              "**Doing it properly.** A grain can hold positive charge up to a surface field of about 10 V/nm "
              "before it sheds ions. A solid-deuterium crystal charged that far (generous: hydrogen ice is weakly "
              "bound) can be accelerated by a voltage to the 100–300 keV per deuteron where D-T fusion peaks. Its "
              "deuterons then stop in a tritide catcher like any beam:", "",
              "| crystal diameter | charge | voltage for 100 keV per deuteron | for 300 keV | best fusion / kinetic energy "
              "| against the wall plug (accelerator taken as 40% efficient) |", "|---|---|---|---|---|---|"]
    energies = np.geomspace(30, 1000, 200)
    y = P.beam_target(energies, "D", "LiT", "Li6")
    gain = sum(v * P.channel_q(k, "Li6") for k, v in y.items()) / (energies / 1e3)
    best = float(gain.max())
    rows = []
    for r_nm in (1.5, 5.0, 15.0):
        a, b = A.crystal_projectile(r_nm, 100.0), A.crystal_projectile(r_nm, 300.0)
        mv = (a["voltage_V"] / 1e6, b["voltage_V"] / 1e6)
        rows.append({"radius_nm": r_nm, "charge_e": a["charge_e"], "MV_100": mv[0], "MV_300": mv[1]})
        lines.append(f"| {2 * r_nm:g} nm | {a['charge_e']:.0f} e | {fmt(mv[0])} MV | {fmt(mv[1])} MV "
                     f"| {best:.2g} | {best * ACCELERATOR_EFFICIENCY:.2g} |")
    laser_route = 5e-4 * LASER_EFFICIENCY
    lines += ["", f"* **It is the cold catcher again.** Whatever throws the deuteron, it returns at most "
              f"{best:.1%} of its energy before the catcher's electrons stop it.",
              "* **But the driver is far better.** Megavolt accelerators exist (dust accelerators throw "
              "micrometre grains at up to 80 km/s with 2 MV; this needs 3,000 km/s). Taking one as 40% efficient, this "
              f"gives {best * ACCELERATOR_EFFICIENCY:.1g} from the wall plug, against ~{fmt(laser_route)} for "
              f"laser-exploded crystals. It is still about {round(1 / (best * ACCELERATOR_EFFICIENCY), -2):.0f} times short, and a "
              "plain deuteron beam does the same.",
              "* **No collective help.** A 30 nm crystal's track is heated to roughly 1 keV as it stops, which "
              "would cut the electron drag, but the heat spreads over about half a micrometre within the stopping "
              "time (~0.3 ps) and the track is cold again. (A rough estimate, not modelled.)", ""]
    return lines, {"beuhler": claim, "best_cold_gain": best, "crystals": rows}


def nanowires() -> tuple[list[str], dict]:
    w = NANOWIRE
    slope = 800.0 / math.log(w["above_200keV"] / w["above_1MeV"])
    spec = A.spectrum_yield(slope, w["max_keV"], e_min_kev=200.0)
    neutrons = 0.5 * w["above_200keV"] * spec["reactions_per_ion"]
    ion_energy = w["above_200keV"] * spec["mean_keV"] * P.J_PER_KEV
    measured_ratio = w["neutrons"] * 2 * P.Q_VALUE_MEV["DD"] * 1.602176634e-13 / w["laser_J"]
    record_ratio = w["record_neutrons_per_J"] * 2 * P.Q_VALUE_MEV["DD"] * 1.602176634e-13
    # the same deuterons into a tritide catcher
    e = np.linspace(200.0, w["max_keV"], 4000)
    weight = np.exp(-(e - 200.0) / slope)
    weight /= weight.sum()
    y = P.beam_target(e, "D", "LiT", "Li6")
    mev = sum(v * P.channel_q(k, "Li6") for k, v in y.items())
    tritide = float(np.sum(weight * mev)) * 1.602176634e-13 * w["above_200keV"] / w["laser_J"]
    lines = ["## 5. Check against the best nanostructure experiments: nanowire arrays", "",
             "Arrays of deuterated-plastic nanowires hold the record for D-D fusion per joule of laser light: "
             f"{fmt(w['record_neutrons_per_J'])} neutrons per joule (Curtis et al., Nat. Commun. 9, 1077 (2018)), "
             "about 500 times a flat target. Light gets several micrometres down between the wires, so a volume at "
             "a sixth of solid density is heated at once. This is the nearest real thing to the crystal foam of "
             "section 12.", "",
             f"A later shot with {w['laser_J']:g} J gave {fmt(w['neutrons'])} neutrons (Curtis et al., Phys. Rev. "
             f"Research 3, 043181 (2021)). The authors' own plasma simulation gives {fmt(w['above_200keV'])} "
             f"deuterons above 0.2 MeV and {fmt(w['above_1MeV'])} above 1 MeV, an exponential spectrum with a slope "
             f"of {slope:.0f} keV.", "",
             f"* **Those deuterons stopping in cold plastic give {fmt(neutrons)} neutrons in this project's model.** "
             f"Measured: {fmt(w['neutrons'])}. So this yield is ordinary beam-target fusion, as in section 1.",
             f"* **The fast deuterons carry {ion_energy:.2f} J, {ion_energy / w['laser_J']:.0%} of the light.** "
             "Nanostructures absorb most of the light (70–80% measured), but only a few per cent reaches fast ions.",
             f"* **Fusion energy / laser energy: {fmt(measured_ratio)}** for that shot and {fmt(record_ratio)} for the "
             "record. Ten times better than gas-jet clusters (section 6), and still a millionth.",
             f"* **With a tritide catcher** the same deuterons would return {fmt(tritide)} of the laser energy. "
             "The cold-catcher limit holds.", ""]
    return lines, {"slope_keV": slope, "neutrons_model": neutrons, "ion_energy_J": ion_energy,
                   "fusion_over_laser": measured_ratio, "record_fusion_over_laser": record_ratio,
                   "tritide_fusion_over_laser": tritide}


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    lines = ["# More routes: fission blanket, D-T chain reaction, bare ions, crystal projectiles", "",
             "Produced by `tools/further_routes.py`. Code: `nanofusion/improved/alternatives.py` and "
             "`energy.dt_chain`.", ""]
    data = {}
    for name, fn in (("hybrid", hybrid), ("chain", chain), ("trap", trap), ("projectiles", projectiles),
                     ("nanowires", nanowires)):
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
