"""Cold fusion routes in crystals: muon catalysis in solid D-T, and NASA's lattice confinement fusion.

Neither needs a hot plasma, which is why they are the last candidates for
fusion in a crystal. Writes docs/results/other/summary.md.

    python3 tools/other_routes.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import alternatives as A  # noqa: E402

OUT = Path("docs/results/other")

# --- inputs from the literature -------------------------------------------------
RECORD_FUSIONS = 150  # Los Alamos (Jones et al. 1986), up to 150 +- 20 per muon
STICKING = {"measured, record conditions (Los Alamos)": 0.0045, "measured directly (PSI)": 0.0056,
            "theory": 0.0065, "hoped for": 0.003}
RECORD_DENSITY = 1.3  # relative to liquid-hydrogen density, the high end of the Los Alamos targets
MUON_COST_GEV = {"optimised accelerator design (Kelly et al. 2021)": 4.7,
                 "traditional estimate (Jones)": 5.0,
                 "pessimistic estimate (arXiv:1008.3010)": 8.0}
ELI_NP = {"laser_J": 20.0, "negative_muons": 280}  # arXiv:2503.20904: 20 J, 28 fs laser shot at ELI-NP
KELLY_Q_ELEC_AT_150 = 0.14  # electricity out / in at 150 fusions per muon (J. Phys. Energy 3, 035003)
HEAT_PER_FUSION_MEV = 26.0  # Kelly et al.: heat per catalysed fusion including the blanket

NASA = {"beam_MeV": 2.9, "beam_mA": 15.0, "fusion_neutrons_per_s": 1.55e3, "theory_neutrons_per_s": 1.2e3,
        "d_density_m3": 8.0e28, "mass_g": 480.0}


def muon_section() -> tuple[list[str], dict]:
    data: dict = {}
    # Cycle rate per unit density that reproduces the record with the measured sticking
    lam_record = A.MUON_DECAY_RATE / (1 / RECORD_FUSIONS - STICKING["measured, record conditions (Los Alamos)"])
    lam_per_density = lam_record / RECORD_DENSITY
    lines = ["## 1. Muon-catalysed fusion in solid or compressed D-T", "",
             "Fusions per muon Y = 1 / (λ₀/λ_c + ω_s), with λ₀ = 4.55×10⁵ /s (muon decay), λ_c the cycle rate and "
             f"ω_s the sticking per fusion. The record, {RECORD_FUSIONS} fusions per muon with ω_s = 0.45%, needs "
             f"λ_c = {lam_record:.2g} /s, or {lam_per_density:.2g} /s per unit of liquid-hydrogen density at "
             f"φ = {RECORD_DENSITY}. Denser fuel cycles proportionally faster; the sticking does not change.", "",
             "### Fusions per muon", "",
             "| density φ (× liquid hydrogen) | " + " | ".join(f"ω_s = {w:.2%} ({k})" for k, w in STICKING.items())
             + " |", "|---|" + "---|" * len(STICKING)]
    table = {}
    for phi in (0.5, 1.0, 1.3, 2.0, 5.0, 10.0, 1e9):
        lam = lam_per_density * phi
        row = [A.fusions_per_muon(lam, w) for w in STICKING.values()]
        table[phi] = row
        label = "unlimited" if phi > 1e6 else f"{phi:g}"
        lines.append(f"| {label} | " + " | ".join(f"{y:.0f}" for y in row) + " |")
    data["fusions_per_muon"] = table
    lines += ["", "Denser fuel helps only until sticking takes over: Y can never exceed 1/ω_s, about 220 "
              "with the measured sticking.", "", "### Energy", "",
              "| muon cost | fusions needed to repay the muon (17.6 MeV each) | gain at the record (150) "
              "| best possible gain (Y = 1/ω_s = 222) | with blanket heat (26 MeV per fusion), best |",
              "|---|---|---|---|---|"]
    y_max = 1 / STICKING["measured, record conditions (Los Alamos)"]
    data["gains"] = {}
    for label, cost in MUON_COST_GEV.items():
        row = {"needed": A.fusions_needed(cost), "record": A.muon_gain(RECORD_FUSIONS, cost),
               "best": A.muon_gain(y_max, cost), "best_heat": A.muon_gain(y_max, cost, HEAT_PER_FUSION_MEV)}
        data["gains"][label] = row
        lines.append(f"| {cost:g} GeV, {label} | {row['needed']:.0f} | {row['record']:.2f} | {row['best']:.2f} "
                     f"| {row['best_heat']:.2f} |")
    laser_gev = ELI_NP["laser_J"] / ELI_NP["negative_muons"] / 1.602176634e-10
    data["laser_muon_GeV"] = laser_gev
    lines.append(f"| {laser_gev:.1g} GeV of laser light per muon: laser-made muons (ELI-NP, 20 J → 280 μ⁻) "
                 f"| {A.fusions_needed(laser_gev):.1g} | {A.muon_gain(RECORD_FUSIONS, laser_gev):.1g} "
                 f"| {A.muon_gain(y_max, laser_gev):.1g} | {A.muon_gain(y_max, laser_gev, HEAT_PER_FUSION_MEV):.1g} |")
    needed_elec = RECORD_FUSIONS / KELLY_Q_ELEC_AT_150
    data["fusions_for_electric_breakeven"] = needed_elec
    lines += ["", f"**Electricity.** Kelly et al. (2021) find electricity out / in = {KELLY_Q_ELEC_AT_150:.0%} "
              f"at 150 fusions per muon for their optimised plant. Since that ratio is proportional to the fusions per "
              f"muon, electrical breakeven needs about {needed_elec:.0f} fusions per muon. That needs sticking below "
              f"{1 / needed_elec:.2%}, about {STICKING['measured, record conditions (Los Alamos)'] * needed_elec:.0f} "
              "times lower than measured.", ""]
    return lines, data


def lattice_section() -> tuple[list[str], dict]:
    data: dict = {}
    beam_w = NASA["beam_MeV"] * NASA["beam_mA"] * 1e-3 * 1e6  # MeV * A = MW -> W
    fusions = 2 * NASA["fusion_neutrons_per_s"]  # half the D-D reactions make a neutron
    fusion_w = fusions * A.DD_FUSION_MEV * 1.602176634e-13
    measured = fusion_w / beam_w
    data["measured"] = {"beam_W": beam_w, "fusion_W": fusion_w, "ratio": measured}
    lines = ["## 2. Lattice confinement fusion (NASA Glenn, 2020)", "",
             "**Measured.** NASA aimed a 2.9 MeV, 15 mA electron beam (Dynamitron, tantalum converter) at 480 g of "
             "erbium deuteride (8×10²² D/cm³) in 16 vials. The fusion neutrons, scaled to the samples, came to "
             f"{NASA['fusion_neutrons_per_s']:.2g} per second, against {NASA['theory_neutrons_per_s']:.2g} predicted by "
             "NASA's screening theory (NASA/TP-2020-5001616, Phys. Rev. C 101, 044610).", "",
             f"* Beam power {beam_w / 1e3:.1f} kW; D-D fusion power {fusion_w:.2g} W "
             f"({fusions:.2g} reactions/s × 3.65 MeV).",
             f"* **Fusion energy / beam energy = {measured:.1g}**, before counting the accelerator's own losses.", "",
             "**Model.** Photoneutrons from photons just above the 2.2246 MeV threshold (mean 145 keV) knock "
             "deuterons to ≤ 8/9 of their energy (mean 64 keV, as NASA calculated). The deuterons then fuse with the "
             "cold lattice deuterons as they stop. The cross section is taken at E + U_e, where U_e is the electron "
             "screening. Stopping in the deuteride uses this project's model, which is probably low by a factor of "
             "2–3 in erbium, so these yields are generous.", "",
             "| screening U_e | fusions per photoneutron, 1-cm vials (0.3 deuteron hits) | in an infinite sample "
             "(every hit, ~7 per neutron) | best possible energy ratio |", "|---|---|---|---|"]
    neutrons = A.photoneutron_energies()
    rows = {}
    for u in (0, 300, 800, 3000):
        fin = A.lattice_fusions_per_neutron(neutrons, "ErD3", NASA["d_density_m3"], u, collisions=0.3)
        inf = A.lattice_fusions_per_neutron(neutrons[:5000], "ErD3", NASA["d_density_m3"], u)
        # Best case: every beam electron's energy becomes photons just above threshold, each splitting a deuteron.
        photon_mev = A.DEUTERON_BINDING_MEV + 2 * float(neutrons.mean()) / 1e3
        best = inf["fusions_per_neutron"] * A.DD_FUSION_MEV / photon_mev
        rows[u] = {"finite": fin["fusions_per_neutron"], "infinite": inf["fusions_per_neutron"],
                   "hits": inf["deuteron_hits_per_neutron"], "best_ratio": best}
        lines.append(f"| {u} eV | {fin['fusions_per_neutron']:.2g} | {inf['fusions_per_neutron']:.2g} | {best:.1g} |")
    data["model"] = rows
    lines += ["", "Screening in the metal, even at 3 keV (four times the largest measured value, 800 eV in palladium), "
              "raises the yield by less than 20%. The deuterons have tens of keV, far above the screening energy. "
              "The best possible ratio assumes that every joule of beam becomes a photon just above threshold and "
              "that every photon splits a deuteron. In the real setup, most of the beam energy becomes heat and "
              "lower-energy X-rays. The samples also absorb only about 10⁻⁴ of the photons that are above threshold "
              "(about 1 mb cross section through ~1 cm). That is where the remaining factor of ~10⁷ goes. NASA's "
              "own detailed calculation of the whole chain predicted 1.2×10³ neutrons/s, close to what it "
              "measured.", ""]
    # What if the neutrons came from elsewhere (for example D-T fusion neutrons)?
    lines += ["**Faster neutrons.** If a 2.45 MeV (D-D) or 14.1 MeV (D-T) neutron did the knocking, in an infinite "
              "erbium deuteride sample with 800 eV screening:", "",
              "| neutron | fusions per neutron | fusion energy / neutron energy |", "|---|---|---|"]
    data["fast"] = {}
    for e_mev in (2.45, 14.1):
        res = A.lattice_fusions_per_neutron(np.full(4000, e_mev * 1e3), "ErD3", NASA["d_density_m3"], 800)
        ratio = res["fusions_per_neutron"] * A.DD_FUSION_MEV / e_mev
        data["fast"][e_mev] = {"fusions": res["fusions_per_neutron"], "ratio": ratio}
        lines.append(f"| {e_mev:g} MeV | {res['fusions_per_neutron']:.2g} | {ratio:.1g} |")
    lines += ["", "So even a fusion neutron returns only about a ten-thousandth of its energy this way. It is the "
              "cold-catcher limit again (energy.md, section 1): a fast deuteron in cold matter loses its energy to "
              "electrons long before it fuses.", ""]
    return lines, data


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    lines = ["# Cold fusion routes in crystals: muon catalysis and lattice confinement", "",
             "Produced by `tools/other_routes.py`. Code: `nanofusion/improved/alternatives.py`.", ""]
    mu_lines, mu_data = muon_section()
    lcf_lines, lcf_data = lattice_section()
    lines += mu_lines + lcf_lines
    text = "\n".join(lines) + "\n"
    (OUT / "summary.md").write_text(text)
    (OUT / "results.json").write_text(json.dumps({"muon": mu_data, "lattice": lcf_data}, indent=1, default=str))
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
