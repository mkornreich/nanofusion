"""What the static inertial-fusion model leaves out, simulated: the 1-D implosion-and-burn model on NIF's shots.

The static model (compressed_icf; tools/icf_experiment.py) has no time in it.
Here nanofusion/improved/burn1d.py follows the shell from peak velocity
through stagnation, ignition, burn and disassembly, and
* isolates each piece of physics the static model leaves out;
* retrodicts Hohenberger et al.'s adiabat experiment, which the static model
  gets wrong;
* fits each shot's degradation, as a mass of ablator mixed into the hot spot;
* re-predicts the shot campaign with the record shot's degradation held fixed.

Writes docs/results/icf_experiment/dynamics.md and dynamics.json (about 3 minutes on 8 cores).

    python3 tools/icf_dynamics.py
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

from nanofusion.improved import burn1d as B
from nanofusion.improved import simulator as S

OUT = Path("docs/results/icf_experiment")
# LLNL post-shot values (Kritcher et al., PRE 106, 025201 (2022) and PRE 109, 025204 (2024)): implosion velocity,
# design adiabat, remaining ablator (share of the diamond capsule's 3.927 / 4.253 mg).
N08 = {"velocity_kms": 401.0, "adiabat": 2.65, "ablator_mg": 0.038 * 3.927}
N22 = {"velocity_kms": 394.0, "adiabat": 2.45, "ablator_mg": 0.057 * 4.253}
MEASURED = {"N210808": (N08, 1.37e6, 1.92e6), "N221204": (N22, 3.15e6, 2.05e6), "N250406": (N22, 8.6e6, 2.08e6)}
LLNL_1D = {"N210808": "4.2 MJ (HYDRA), 6.2 MJ (xRAGE)", "N221204": "~15 MJ (idealised 1-D at ~2.1 MJ); 13.6 MJ in 2-D "
           "without low modes"}
KEYS = ("yield_J", "burn_up", "burn_pressure_gbar", "burn_T_keV", "peak_rhoR", "burn_rhoR", "burn_width_ps",
        "bang_time_ps", "kinetic_J", "energy_error")


def sim(kw: dict) -> dict:
    warnings.simplefilter("ignore")
    r = B.run(B.Capsule(**kw))
    return {k: float(r[k]) for k in KEYS}


def pool():
    return ProcessPoolExecutor(8, mp_context=mp.get_context("fork"))


def fit_mix(args) -> dict:
    """Mixed-in ablator mass (ng) that brings the clean implosion down to the measured yield (bisection in log)."""
    kw, target = args
    lo, hi = 1.0, 30000.0
    for _ in range(14):
        mid = math.sqrt(lo * hi)
        lo, hi = (mid, hi) if sim(dict(kw, mix_ng=mid))["yield_J"] > target else (lo, mid)
    best = math.sqrt(lo * hi)
    return {"mix_ng": best, **sim(dict(kw, mix_ng=best))}


def at_fixed_kinetic(kw: dict, velocity: float) -> dict:
    """A faster (or slower) shell from a thinner (thicker) ablator left, with the same kinetic energy."""
    ke = 0.5 * (kw.get("fuel_mg", 0.22) + kw["ablator_mg"]) * (kw["velocity_kms"]) ** 2
    return dict(kw, velocity_kms=velocity, ablator_mg=2 * ke / velocity**2 - kw.get("fuel_mg", 0.22))


def with_laser(kw: dict, laser_ratio: float, fuel_mg: float | None = None, ablator_mg: float | None = None) -> dict:
    """Kinetic energy in proportion to the laser energy; the fuel and ablator as given; velocity follows."""
    fuel = kw.get("fuel_mg", 0.22)
    ke = 0.5 * (fuel + kw["ablator_mg"]) * kw["velocity_kms"] ** 2 * laser_ratio
    fuel = fuel_mg or fuel
    abl = ablator_mg or kw["ablator_mg"]
    return dict(kw, fuel_mg=fuel, ablator_mg=abl, velocity_kms=math.sqrt(2 * ke / (fuel + abl)))


def scaled(kw: dict, laser_ratio: float) -> dict:
    """A hydrodynamically scaled capsule: every mass x laser_ratio, every length x laser_ratio^(1/3)."""
    s = laser_ratio ** (1 / 3)
    out = dict(kw, fuel_mg=kw.get("fuel_mg", 0.22) * s**3, ablator_mg=kw["ablator_mg"] * s**3,
               vapour_ug=kw.get("vapour_ug", 1.8) * s**3, inner_radius_um=kw.get("inner_radius_um", 300.0) * s)
    if kw.get("mix_ng"):
        out["mix_ng"] = kw["mix_ng"] * s**2  # mixing comes in through the surface
    return out


def campaign_steps(base: dict) -> list[tuple[str, str, dict]]:
    """The shot campaign, each step as 1-D inputs. base: N221204's implosion with a shot's degradation."""
    r22 = 2.2 / 2.05
    r26 = 2.6 / 2.05
    thick = with_laser(base, r22, fuel_mg=0.22 * r22)
    big = scaled(with_laser(base, 1.0), r26)
    big = dict(big, velocity_kms=base["velocity_kms"])
    return [
        ("A. Repeat the record", "today's target and pulse", base),
        ("B. Lower adiabat", "adiabat 2.45 → 2.25", dict(base, adiabat=2.25)),
        ("C. Faster", "+4% velocity from a thinner ablator, same kinetic energy", at_fixed_kinetic(base, 410.0)),
        ("D0. 2.2 MJ, thicker ablator", "what NIF did: 5% more ablator left", with_laser(base, r22, ablator_mg=base[
            "ablator_mg"] * 90 / 85.4)),
        ("D. 2.2 MJ, thicker fuel", "Landen et al.'s proposal: 7% more D-T", thick),
        ("E. D + lower adiabat", "", dict(thick, adiabat=2.25)),
        ("F. D + lower adiabat + faster", "", at_fixed_kinetic(dict(thick, adiabat=2.25), thick["velocity_kms"] * 1.04)),
        ("G. 2.6 MJ, scaled capsule", "every length × 1.08", big),
        ("H. G + lower adiabat + faster", "", at_fixed_kinetic(dict(big, adiabat=2.25), big["velocity_kms"] * 1.04)),
        ("I. G + adiabat 2.0 + 8% faster", "a stretch", at_fixed_kinetic(dict(big, adiabat=2.0),
                                                                          big["velocity_kms"] * 1.08)),
    ]


def compute() -> dict:
    jobs = {}
    for name, b in (("N210808", N08), ("N221204", N22)):
        jobs[f"{name} clean"] = dict(b)
        jobs[f"{name} no alphas"] = dict(b, alpha_heating=False)
    for ng in (30, 100, 300, 1000):
        jobs[f"mix {ng}"] = dict(N22, mix_ng=ng)
    for rke in (0.01, 0.03, 0.1):
        jobs[f"rke {rke}"] = dict(N22, rke=rke)
    for v in (360, 380, 410, 430):
        jobs[f"velocity {v}"] = at_fixed_kinetic(N22, v)
        jobs[f"velocity {v} no alphas"] = dict(at_fixed_kinetic(N22, v), alpha_heating=False)
    for a in (2.89, 2.68):
        jobs[f"adiabat {a} (N210808 target)"] = dict(N08, adiabat=a)
        jobs[f"adiabat {a} (N210808 target) no alphas"] = dict(N08, adiabat=a, alpha_heating=False)
    jobs["adiabat 2.25 no alphas"] = dict(N22, adiabat=2.25, alpha_heating=False)
    jobs["adiabat 2.25"] = dict(N22, adiabat=2.25)
    names = list(jobs)
    with pool() as ex:
        runs = dict(zip(names, ex.map(sim, [jobs[n] for n in names])))
        fits = dict(zip(MEASURED, ex.map(fit_mix, [(kw, y) for kw, y, _ in MEASURED.values()])))
    for n in names:
        runs[n]["inputs"] = jobs[n]
    # Hohenberger et al.: same N210808 target, adiabat 2.89 -> 2.68, with N210808's degradation
    deg = fits["N210808"]["mix_ng"]
    hoh = {a: dict(N08, adiabat=a, mix_ng=deg) for a in (2.89, 2.68)}
    camp = {}
    camp["record"] = campaign_steps(dict(N22, mix_ng=fits["N250406"]["mix_ng"]))
    camp["clean"] = campaign_steps(dict(N22))
    with pool() as ex:
        hoh_runs = dict(zip(hoh, ex.map(sim, hoh.values())))
        flat = [kw for label in camp for _, _, kw in camp[label]]
        results = list(ex.map(sim, flat))
    k = 0
    campaign = {}
    for label, steps in camp.items():
        campaign[label] = []
        for name, practice, kw in steps:
            campaign[label].append({"step": name, "practice": practice, "inputs": kw, **results[k]})
            k += 1
    return {"runs": runs, "fits": fits, "hohenberger": {str(a): {"inputs": hoh[a], **hoh_runs[a]} for a in hoh},
            "campaign": campaign}


def mj(x: float) -> str:
    return f"{x / 1e6:.3g} MJ" if x >= 1e5 else f"{x / 1e3:.2g} kJ"


def write(d: dict) -> list[str]:
    runs, fits, camp = d["runs"], d["fits"], d["campaign"]
    c08, c22 = runs["N210808 clean"], runs["N221204 clean"]
    n22 = runs["N221204 no alphas"]
    static = {"N210808": None, "N221204": None}
    for name, kw in (("N210808", N08), ("N221204", N22)):
        p = {"N210808": 569.0, "N221204": 608.0}[name]
        rho = S.icf_shell_density(p, kw["adiabat"])
        rr = S.icf_shell_for_mass(p, kw["adiabat"], 0.22e-6)
        static[name] = S.icf_target(rho, rr, kw["adiabat"])["fusion_J"]
    lines = ["# What the static model leaves out, simulated", "",
             ("Produced by `tools/icf_dynamics.py` with the 1-D implosion-and-burn model "
              "`nanofusion/improved/burn1d.py`. The proposal is [docs/icf_experiment.md](../../icf_experiment.md); "
              "the static model's numbers are in [summary.md](summary.md)."), "",
             "## The model", "",
             "* The D-T shell and the ablator left on it coast inward from peak velocity (300 µm, 4.5 g/cm³).",
             "* They stagnate on the D-T gas (1.8 µg), heat it, and ablate fuel into it.",
             "* The hot spot ignites, or not, and the burn spreads and disassembles the capsule.",
             ("* Included: alphas deposited by range, flux-limited conduction, bremsstrahlung with free-free "
             "re-absorption, an ideal-plus-degenerate equation of state, burn-up."),
             ("* Checked: the same yield within 1% at 2.6× the zones and at half the time step; energy conserved to "
             "0.5%."), "",
             "## Calibration: the clean implosion", "",
             ("| Shot | Inputs (LLNL post-shot) | 1-D yield here | LLNL's 1-D | Static ceiling | Measured | "
             "Measured ÷ 1-D |"), "|---|---|---|---|---|---|---|"]
    for name, c, kw in (("N210808", c08, N08), ("N221204", c22, N22)):
        measured = MEASURED[name][1]
        lines.append(f"| {name} | {kw['velocity_kms']:.0f} km/s, adiabat {kw['adiabat']}, "
                     f"{kw['ablator_mg'] * 1e3:.0f} µg of ablator | {mj(c['yield_J'])} | {LLNL_1D[name]} | "
                     f"{mj(static[name])} | {mj(measured)} | {measured / c['yield_J']:.0%} |")
    lines += ["",
              (f"* The record, N250406 (8.6 MJ), reached {8.6e6 / c22['yield_J']:.0%} of the clean N221204-class "
              f"yield. LLNL puts N240210 at about 35% of its 1-D yield; their goal is 80%."),
              (f"* **The static ceiling is too low.** The clean 1-D yield of N221204's design is {mj(c22['yield_J'])}, "
              f"about {c22['yield_J'] / static['N221204']:.1f}× the static model's {mj(static['N221204'])}. So the "
              "static model's \"quality\" overstated how close NIF has come."), "",
              "## The physics the static model leaves out", "",
              "Each piece isolated on N221204's clean implosion:", "",
              "| Physics | Simulated | Effect on the yield |", "|---|---|---|",
              (f"| Alpha heating | off: {mj(n22['yield_J'])} at {n22['burn_T_keV']:.1f} keV | ×{c22['yield_J'] / n22['yield_J']:.0f} "
              f"(N210808: ×{c08['yield_J'] / runs['N210808 no alphas']['yield_J']:.0f}) |"),
              (f"| Burn propagation and disassembly | burn-up {c22['burn_up']:.0%} | the static formula ρR/(ρR + 7) "
              f"gives {n22['peak_rhoR'] / (n22['peak_rhoR'] + 7):.0%} at the same areal density |")]
    for ng in (100, 300, 1000):
        r = runs[f"mix {ng}"]
        lines.append(f"| Ablator mixed into the hot spot | {ng} ng of carbon | ×{r['yield_J'] / c22['yield_J']:.2f} |")
    for rke in (0.01, 0.03, 0.1):
        r = runs[f"rke {rke}"]
        lines.append(f"| Asymmetry: kinetic energy that never stagnates | {rke:.0%} | ×{r['yield_J'] / c22['yield_J']:.2f} |")
    vs = [360, 380, 394, 410, 430]
    p_no = {v: (n22 if v == 394 else runs[f"velocity {v} no alphas"])["burn_pressure_gbar"] for v in vs}
    y = {v: (c22 if v == 394 else runs[f"velocity {v}"])["yield_J"] for v in vs}
    slope = math.log(p_no[430] / p_no[360]) / math.log(430 / 360)
    lines.append(f"| The cost of pressure: velocity at fixed kinetic energy | {', '.join(f'{v}' for v in vs)} km/s | "
                 f"{', '.join(mj(y[v]) for v in vs)}; pressure without alphas ∝ v^{slope:.1f} |")
    a25 = runs["adiabat 2.25"]
    a25n = runs["adiabat 2.25 no alphas"]
    lines.append(f"| The adiabat at fixed velocity | 2.45 → 2.25 | ×{a25['yield_J'] / c22['yield_J']:.2f}; pressure "
                 f"+{a25n['burn_pressure_gbar'] / n22['burn_pressure_gbar'] - 1:.0%}, areal density "
                 f"+{a25n['peak_rhoR'] / n22['peak_rhoR'] - 1:.0%} (no alphas) |")
    lines += ["",
              ("* **Alpha heating is most of the yield.** Without it the capsule makes tens of kJ. With it, the hot "
              "spot runs away to 20–30 keV and the burn spreads through the shell."),
              ("* **A clean capsule shrugs off small degradations; one near the cliff doesn't.** 100 ng of mixed-in "
              "carbon costs a clean N221204 2%, but the measured NIF trend is 30–50% per 100 ng, on degraded shots "
              "closer to the cliff."),
              "* **Velocity has a cliff.** Below ~380 km/s the capsule fails to ignite.",
              "* **The adiabat matters through pressure and areal density**, which the static model held fixed.", "",
              "## Retrodicting Hohenberger et al. (2025): the adiabat test the static model fails", "",
              ("NIF lowered the adiabat alone, 2.89 → 2.68, on N210808's target at the same laser energy. Areal density "
              "(DSR) rose 14% and the yield up to 80%."), "",
              "| Model | Yield ratio, adiabat 2.68 ÷ 2.89 | Areal density ratio |", "|---|---|---|",
              "| Static (fixed pressure and fuel mass) | ×1.01 | ×1.00 |"]
    h89, h68 = runs["adiabat 2.89 (N210808 target)"], runs["adiabat 2.68 (N210808 target)"]
    h89n, h68n = runs["adiabat 2.89 (N210808 target) no alphas"], runs["adiabat 2.68 (N210808 target) no alphas"]
    d89, d68 = d["hohenberger"]["2.89"], d["hohenberger"]["2.68"]
    lines += [(f"| 1-D, clean | ×{h68['yield_J'] / h89['yield_J']:.1f} | ×{h68n['peak_rhoR'] / h89n['peak_rhoR']:.2f} "
              "(no alphas) |"),
              (f"| 1-D, with N210808's degradation | ×{d68['yield_J'] / d89['yield_J']:.1f} | "
              f"×{d68['peak_rhoR'] / d89['peak_rhoR']:.2f} |"),
              "| Measured | up to ×1.8 | ×1.14 |", "",
              ("* The 1-D model gets the direction right, where the static model can't. It overshoots the size: "
              "N210808's target sits on the ignition cliff in 1-D, so small changes in areal density move the yield "
              "a lot. Hohenberger et al. find the same steepness (yield ∝ ρR^11)."), "",
              "## Each shot's degradation", "",
              ("Every degradation (mixing, low-mode asymmetry, defects) is folded into one number: the carbon mixed "
              "into the hot spot that brings the clean 1-D yield down to the measured one."), "",
              "| Shot | Clean 1-D | Measured | Mix-equivalent carbon | Burn-up | Burn-averaged T |",
              "|---|---|---|---|---|---|"]
    for name, f in fits.items():
        clean = c08 if name == "N210808" else c22
        lines.append(f"| {name} | {mj(clean['yield_J'])} | {mj(MEASURED[name][1])} | {f['mix_ng']:.0f} ng | "
                     f"{f['burn_up']:.1%} | {f['burn_T_keV']:.1f} keV |")
    f08, f22, f50 = fits["N210808"]["mix_ng"], fits["N221204"]["mix_ng"], fits["N250406"]["mix_ng"]
    share = MEASURED["N250406"][1] / c22["yield_J"]
    lines += ["",
              (f"* N210808's equivalent, {f08:.0f} ng, is within NIF's measured bound (less than 250 ng of ablator in "
              f"its hot spot and its repeats). N221204's, {f22:.0f} ng, is not: there the number stands in for "
              "low-mode asymmetry, which LLNL finds cost N221204 a factor of 4.5."),
              (f"* **The degraded shots sit on the cliff.** The record's equivalent ({f50:.0f} ng) is only "
              f"{1 - f50 / f22:.0%} below N221204's, yet its yield is {8.6 / 3.15:.1f}× higher. A single mixing "
              "knob is far too steep to hold fixed, so the campaign below is bracketed two ways."), "",
              "## The campaign, re-predicted", "",
              (f"* **Fixed share:** each step's clean 1-D yield × {share:.0%}, the record's share of its own clean "
              "yield (the usual way to carry a degradation forward)."),
              (f"* **Fixed mixing:** the record's {f50:.0f} ng held fixed. Each step that helps pushes the capsule "
              "back up the cliff, so this is an upper bound."), "",
              "| Step | What changes | Velocity | Fuel | Clean 1-D | Fixed share | Fixed mixing |",
              "|---|---|---|---|---|---|---|"]
    for rec, clean in zip(camp["record"], camp["clean"]):
        kw = clean["inputs"]
        lines.append(f"| {rec['step']} | {rec['practice']} | {kw['velocity_kms']:.0f} km/s | "
                     f"{kw.get('fuel_mg', 0.22):.3g} mg | {mj(clean['yield_J'])} | **{mj(share * clean['yield_J'])}** | "
                     f"{mj(rec['yield_J'])} |")
    by = {r["step"].split(".")[0]: r for r in camp["clean"]}
    fix = {r["step"].split(".")[0]: r for r in camp["record"]}

    def ratio(k: str, base: str = "A") -> str:
        return f"×{by[k]['yield_J'] / by[base]['yield_J']:.2f} (fixed mixing ×{fix[k]['yield_J'] / fix[base]['yield_J']:.1f})"

    lines += ["",
              (f"* **B, lower adiabat:** {ratio('B')}. Not a null result: the static model's central claim fails, "
              "as Hohenberger et al. found on NIF."),
              f"* **C, faster at the same kinetic energy:** {ratio('C')}, if the thinner ablator doesn't mix more.",
              (f"* **D vs D0, where to spend 2.2 MJ:** more fuel {mj(share * by['D']['yield_J'])}, more ablator (NIF's "
              f"choice) {mj(share * by['D0']['yield_J'])}. Landen et al. (2024) proposed the switch; the 1-D model "
              "agrees it helps, modestly."),
              (f"* **The 2.6 MJ upgrade:** {mj(share * by['G']['yield_J'])}–{mj(fix['G']['yield_J'])} for a scaled "
              f"capsule, up to {mj(share * by['I']['yield_J'])}–{mj(fix['I']['yield_J'])} stretched. At LLNL's goal of "
              f"80% of the 1-D yield: {mj(0.8 * by['G']['yield_J'])}–{mj(0.8 * by['I']['yield_J'])}. LLNL expects "
              "\"greater than 30 MJ\"; Lindl et al. very likely > 20 MJ, 50–70 MJ optimised."), ""]
    return lines


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    cache = OUT / "dynamics.json"
    if "--reuse" in sys.argv and cache.exists():
        d = json.loads(cache.read_text())
    else:
        d = compute()
        cache.write_text(json.dumps(d, indent=1))
    lines = write(d)
    (OUT / "dynamics.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
