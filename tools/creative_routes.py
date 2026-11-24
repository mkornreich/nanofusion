"""Creative routes: squeeze a crystal until it fuses, make muons with exploding crystals, polarise the fuel.

1. Cold squeeze (pycnonuclear fusion). Nuclei in a crystal lattice tunnel
   into each other through zero-point motion alone, as in white dwarfs.
   Compressing cold fuel costs little energy compared with heating it. How
   dense must it be?
2. Exploding crystals as a muon factory. Muon catalysis came closest to
   breakeven (energy.md, section 9), but muons are expensive. Micrometre
   crystals Coulomb-exploding to ~1 GeV per nucleon could drive the pion
   target, and their isotropic emission is no handicap if the target
   surrounds them.
3. Boosters. Spin-polarised D-T (cross section x 1.5) applied to the two
   routes that came closest among plasma schemes; X-ray stripping of stuck
   muons.

Writes docs/results/creative/summary.md.

    python3 tools/creative_routes.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import alternatives as A  # noqa: E402
from nanofusion.improved import energy as E  # noqa: E402

OUT = Path("docs/results/creative")
RECORD_INTENSITY = 1.1e23  # W/cm^2, CoReLS 2021 (Optica 8, 630)
NIF_DENSITY = 1000.0  # g/cm^3, roughly the peak fuel density of NIF-class implosions
NIF_FUEL_RADIUS = 30e-6  # m, order of a compressed NIF fuel shell
J_PER_GEV = 1.602176634e-10


def cold_squeeze() -> tuple[list[str], dict]:
    rows = []
    lines = ["## 1. Cold squeeze: can a crystal fuse by compression alone?", "",
             "Pycnonuclear rate of Salpeter & Van Horn (ApJ 155, 183 (1969), static bcc lattice), with the "
             "fuel kept at zero temperature, so its pressure is that of the degenerate electrons. Ignition: the "
             "pycnonuclear burn must consume 1% of the fuel before the ball flies apart (τ = R / 4c_s). The alpha "
             "particles then heat it to thermonuclear temperatures, as in a white dwarf igniting.", "",
             "| D-T density (g/cm³) | × solid | reactions per nucleus per s | pressure (Gbar) "
             "| compression in radius from solid | ball radius that ignites | its mass | ideal compression energy |",
             "|---|---|---|---|---|---|---|---|"]
    for rho in (0.21, 20, 1e3, 1e4, 1e5, 1e6, 3e6, 1e7, 3e7):
        c = A.cold_compression(rho, "DT")
        rows.append(c)
        big = c["ignition_radius_m"] > 1e3
        radius = "more than a planet" if c["ignition_radius_m"] > 1e7 else (
            f"{c['ignition_radius_m'] / 1e3:.2g} km" if big else f"{c['ignition_radius_m']:.2g} m")
        mass = "–" if c["ignition_radius_m"] > 1e7 else f"{c['mass_kg']:.2g} kg"
        energy = "–" if c["ignition_radius_m"] > 1e7 else f"{c['compression_J']:.2g} J"
        lines.append(f"| {rho:g} | {rho / 0.21:.2g} | {c['rate_per_nucleus_s']:.2g} | {c['pressure_Gbar']:.2g} "
                     f"| {c['convergence_from_solid']:.3g} | {radius} | {mass} | {energy} |")
    nif = A.cold_compression(NIF_DENSITY, "DT")
    needed_rate = 0.01 * 4 * nif["sound_speed_m_s"] / NIF_FUEL_RADIUS
    enhancement = needed_rate / nif["rate_per_nucleus_s"]
    at = {r["rho_g_cc"]: r for r in rows}
    c7, c37, c36 = at[1e7], at[3e7], at[3e6]
    lines += ["", f"* **It takes about 10⁷ g/cm³**, white-dwarf density. Below that, the fuel flies apart long "
              f"before it burns. At 3×10⁶ g/cm³ the ball would have to be {c36['mass_kg']:.1g} kg, with "
              f"{c36['compression_J'] / 1e12:.1g} TJ of compression. At 10⁷ g/cm³, {c7['mass_kg'] * 1e6:.2g} mg "
              f"would ignite after {c7['compression_J'] / 1e6:.1f} MJ of ideal, shock-free compression. At "
              f"3×10⁷ g/cm³ even a nanogram would, with {c37['compression_J']:.0f} J. The ideal gain is about "
              f"{c7['gain_ideal']:.0f}: fusion releases far more than squeezing the electrons costs.",
              f"* **But the conditions are out of reach.** That is ~10⁴ times the density NIF reaches, at "
              f"~{c7['pressure_Gbar']:.0e} Gbar (NIF's hot spot is a few hundred Gbar). It means shrinking the fuel "
              f"{c7['convergence_from_solid']:.0f}-fold in radius, against ~30–40 at the limit of hydrodynamic "
              "stability. And it must happen without heating it, when every real compression heats by shocks.",
              f"* **At NIF's ~10³ g/cm³** a cold 30 µm fuel shell would need the pycnonuclear rate to be "
              f"{enhancement:.0e} times the standard prediction to ignite.",
              "* **One dissenting theory.** Ichimaru & Kitamura (Phys. Plasmas 6, 2649 (1999)) predicted 'prodigiously "
              "high' p-d rates in liquid metallic hydrogen at ≥ 20 g/cm³ and ≥ 1 Gbar near solidification. If D-T "
              "rates there were high enough to ignite, cold fuel at 20 g/cm³ (where squeezing the electrons costs only "
              "~60 eV each) would burn, and that would change everything. The standard rate there is ~10⁻⁷³ per "
              "nucleus per second. Nobody has compressed hydrogen cold to 1 Gbar (static and ramp experiments reach "
              "~0.005–0.05 Gbar, 0.5–5 TPa), so this is untested, not refuted. Light nuclei sit outside the regime where the "
              "pycnonuclear formulas are reliable, and different lattice models disagree by orders of magnitude.", ""]
    return lines, {"rows": rows, "nif_enhancement_needed": enhancement}


def muon_factory() -> tuple[list[str], dict]:
    fusions = (150, 222)
    lwfa_gev = 20.0 / 280 / J_PER_GEV  # ELI-NP: 20 J -> 280 negative muons
    lines = ["## 2. Exploding crystals as a muon factory", "",
             "Muons are made by pions from ~1 GeV-per-nucleon ions hitting a target. The best accelerator designs "
             "need ~4.7 GeV of beam per stopped muon (Kelly et al. 2021); deuteron beams ~8 GeV (Jändel/Petrov "
             "studies). A Coulomb-exploding ⁶LiD crystal gives every ion at radius r the same energy per nucleon, "
             "rising as r² to the edge. Below ~0.3 GeV per nucleon they make almost no pions. Here muons per GeV of "
             "beam rise linearly from that threshold to the best value at 1 GeV per nucleon.", "",
             "| edge energy (GeV per nucleon) | crystal radius | stripping intensity (× today's record) "
             "| beam per muon (4.7 GeV best case) | laser per muon, 15% → ions | 50% → ions | 100% → ions "
             "| muon gain at 150 fusions per muon (15% / 50% / 100%) |", "|---|---|---|---|---|---|---|---|"]
    rows = []
    for e_max in (0.5, 1.0, 2.0, 4.0):
        r = {eta: A.crystal_muon_factory(e_max, 4.7, eta) for eta in (0.15, 0.5, 1.0)}
        g = [A.muon_gain(fusions[0], r[eta]["laser_GeV_per_muon"]) for eta in (0.15, 0.5, 1.0)]
        rows.append({"e_max": e_max, **r[0.15], "gains_150": g})
        lines.append(f"| {e_max:g} | {r[0.15]['crystal_radius_m'] * 1e6:.2g} µm | {r[0.15]['strip_intensity_W_cm2']:.1g} W/cm² "
                     f"({r[0.15]['strip_intensity_W_cm2'] / RECORD_INTENSITY:.0f}×) | {r[0.15]['beam_GeV_per_muon']:.3g} GeV "
                     f"| {r[0.15]['laser_GeV_per_muon']:.3g} GeV | {r[0.5]['laser_GeV_per_muon']:.3g} GeV "
                     f"| {r[1.0]['laser_GeV_per_muon']:.3g} GeV | {g[0]:.2f} / {g[1]:.2f} / {g[2]:.2f} |")
    best = A.crystal_muon_factory(4.0, 4.7, 0.15)["laser_GeV_per_muon"]
    lines += ["", f"* **Today's laser muon sources** (laser wakefields, ELI-NP) cost ~{lwfa_gev:.0e} GeV of laser light per "
              f"muon. Exploding crystals would cut that ~{lwfa_gev / best:.0e}-fold on paper.",
              "* **But it is still no better than an accelerator.** Even if all the laser energy became ion energy, a "
              "muon costs the accelerator's ~4.7 GeV, and lasers turn electricity into light far less efficiently "
              "(≤ 20% for future diode-pumped systems) than accelerators turn it into beam.",
              f"* **Gain is 0.08–0.55 against the light**, at {fusions[0]} fusions per muon. Against the wall plug "
              "that is ~0.02–0.1, before turning heat into electricity.",
              "* **The crystals are also out of reach.** They would need to be 2.5–5 µm across, stripped at "
              "10²³–10²⁴ W/cm². Relativistic transparency might let the light in at those intensities, but that is "
              "2–15 times today's record.", ""]
    return lines, {"rows": rows, "lwfa_GeV_per_muon": lwfa_gev}


def polarised() -> tuple[list[str], dict]:
    data: dict = {"hot": {}, "burn": {}}
    radii = np.geomspace(1e-6, 1.0, 241)
    lines = ["## 3. Boosters", "", "### Spin-polarised D-T", "",
             "Aligning the D and T spins raises the D-T cross section by 1.5 (theory, never yet tested in a plasma; "
             "the DIII-D polarised-fuel programme aims to measure whether polarisation survives).", "",
             "| route | unpolarised | polarised |", "|---|---|---|"]
    try:
        for factor in (1.0, 1.5):
            E.DT_SPIN_FACTOR = factor
            for te in (10.0, 30.0):
                pl = E.fuel_plasma("DT", 1.0, te, te)
                data["hot"][(factor, te)] = max(E.beam_in_plasma(e, "D", pl)["gain"] for e in (100, 300, 1000, 3000))
            for label, r_nm, generous in (("generous", 25, True), ("conservative", 15, False)):
                for L in radii:
                    res = E.array_burn("DT", r_nm, 1.0, L, fast_in_final_plasma=generous)
                    if res["gain_after_radiation"] >= 1:
                        data["burn"][(factor, label)] = (L, res["energy_in_J"])
                        break
    finally:
        E.DT_SPIN_FACTOR = 1.0
    for te in (10.0, 30.0):
        lines.append(f"| deuteron fired into a {te:g} keV D-T plasma: fusion energy / ion energy (best ion energy) "
                     f"| {data['hot'][(1.0, te)]:.2g} | {data['hot'][(1.5, te)]:.2g} |")
    for label in ("generous", "conservative"):
        u, p = data["burn"][(1.0, label)], data["burn"][(1.5, label)]
        lines.append(f"| ball of exploding D-T crystals ({label} bound): breakeven ball radius and energy "
                     f"| {u[0] * 1e3:.2g} mm, {u[1]:.2g} J | {p[0] * 1e3:.2g} mm, {p[1]:.2g} J |")
    lines += ["", "Polarisation is the biggest lever on the fuel itself in this project. It cuts the breakeven energy of a "
              "ball of exploding crystals 3–26-fold, because burn in the ball rises steeply with reaction rate. "
              "But the ball is still uncompressed inertial fusion. Its energy must arrive in tens of picoseconds, "
              "and light cannot reach crystals deeper than ~10 nm, so nothing changes in practice. For compressed "
              "fusion, polarisation has been studied as a way to lower ignition thresholds (e.g. Temporal et al., "
              "Nucl. Fusion 52, 103011 (2012)); there it helps real schemes, not these crystals.", "",
              "### Stripping stuck muons with X-rays", "",
              "A stuck muon sits in a muonic helium ion bound by ~11 keV. A Missouri S&T thesis found that an X-ray "
              "fluence of 10²³ photons/cm² would lift muon catalysis to >600 fusions per muon (about 450 more, worth "
              "~7.9 GeV). Even at the minimum photon energy, 11 keV, that fluence carries 1.8×10⁸ J per cm². Paying "
              "for it needs ~1.4×10¹⁷ muons in each cm² of irradiated fuel at once. The most intense muon beams "
              "deliver ~10⁸ per second, so that is billions of beam-seconds of muons for every square centimetre. "
              "A 2026 rate-network analysis (arXiv:2606.07077) finds external-field stripping adds at most ~40% "
              "(113 → 157 fusions per muon in its best case).", ""]
    fluence_j = 1e23 * 11e3 * 1.602176634e-19
    data["xray"] = {"fluence_J_cm2": fluence_j, "muons_per_cm2_needed": fluence_j / (450 * 17.59e-3 * J_PER_GEV)}
    data["hot"] = {f"spin {f:g}, {te:g} keV": v for (f, te), v in data["hot"].items()}
    data["burn"] = {f"spin {f:g}, {label}": {"radius_m": v[0], "energy_J": v[1]} for (f, label), v in data["burn"].items()}
    return lines, data


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    lines = ["# Creative routes: cold squeeze, crystal muon factory, polarised fuel", "",
             "Produced by `tools/creative_routes.py`. Code: `nanofusion/improved/alternatives.py`.", ""]
    data = {}
    for name, fn in (("cold_squeeze", cold_squeeze), ("muon_factory", muon_factory), ("boosters", polarised)):
        part, d = fn()
        lines += part
        data[name] = d
    text = "\n".join(lines) + "\n"
    (OUT / "summary.md").write_text(text)
    (OUT / "results.json").write_text(json.dumps(data, indent=1, default=str))
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
