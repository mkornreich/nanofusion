"""Can laser-exploded nanocrystals give net energy? Three increasingly generous tests.

    python3 tools/energy_experiments.py [--out docs/results/energy]

1. cold   ions stop in cold material; the laser is assumed perfect (all its
          energy becomes ion energy). Fusion energy out per energy in, by ion energy.
2. hot    the same ions slow down in an already-hot plasma instead.
3. burn   a whole ball packed with crystals is exploded at once and burns as a
          hot plasma until it flies apart (0-D, with X-ray losses).
Writes summary.md and JSON tables to the output directory.
"""

from __future__ import annotations

import argparse
import json
import math
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import energy as E  # noqa: E402
from nanofusion.improved import physics as P  # noqa: E402


def coulomb_energy_J(fuel: str, radius_m: float) -> float:
    """Electrostatic energy of a fully stripped uniform sphere: (3/5) Q^2 / (4 pi eps0 R)."""
    q = E.charge_density(fuel) * 4 / 3 * math.pi * radius_m**3
    return 0.6 * P.K_COULOMB * q * q / radius_m


def cold_catcher(lines: list[str], data: dict) -> None:
    lines += ["## 1. Ions stopping in a cold catcher, with a perfect laser", "",
              "Fusion energy per unit of ion energy, for one ion fired into cold material "
              "(the best a laser could do is turn all its energy into ion energy).", "",
              "| Ion energy | D into LiT | D into ⁶LiD | Li-6 into ⁶LiD |", "|---|---|---|---|"]
    energies = np.geomspace(10, 2e4, 400)
    curves = {}
    for label, proj, catcher in (("D_LiT", "D", "LiT"), ("D_LiD", "D", "LiD"), ("Li_LiD", "Li", "LiD")):
        y = P.beam_target(energies, proj, catcher, "Li6")
        mev = sum(v * P.channel_q(k, "Li6") for k, v in y.items())
        curves[label] = mev / (energies / 1e3)
    for e in (30, 100, 300, 1000, 3000, 10000):
        i = int(np.argmin(abs(energies - e)))
        lines.append(f"| {e:g} keV | {curves['D_LiT'][i]:.2g} | {curves['D_LiD'][i]:.2g} | {curves['Li_LiD'][i]:.2g} |")
    best = {k: (float(v.max()), float(energies[int(v.argmax())])) for k, v in curves.items()}
    lines += ["", "Best ratio: " + "; ".join(f"{k.replace('_', ' into ')} {g:.2g} at {e:.0f} keV"
                                             for k, (g, e) in best.items()) + ".", ""]
    # Laser efficiency actually needed: Coulomb energy vs laser energy through a crystal at the stripping threshold
    rows = []
    for r_nm in (5, 15, 50):
        r = r_nm * 1e-9
        e_strip = E.charge_density("6LiD") * r / (3 * P.EPS0)
        pulse = P.Pulse(e0=e_strip, wavelength=0.8e-6, fwhm=6e-15)
        laser = pulse.fluence() * math.pi * r * r
        rows.append(f"{2 * r_nm} nm crystal: {coulomb_energy_J('6LiD', r) / laser:.0%}")
    lines += ["How much of a 6 fs pulse at the stripping threshold ends up as explosion energy "
              "(Coulomb energy / laser energy through the crystal): " + ", ".join(rows) + ".", ""]
    data["cold"] = {"best": best, "energies_keV": energies.tolist(), **{k: v.tolist() for k, v in curves.items()}}


def hot_catcher(lines: list[str], data: dict) -> None:
    lines += ["## 2. The same ions fired into an already-hot plasma", "",
              "Fusion energy per unit of ion energy, for ions from a 30 nm ⁶LiD crystal (mean D ≈ 160 keV, "
              "Li ≈ 490 keV) slowing down in a solid-density plasma at electron temperature T_e "
              "(the hot plasma itself is free in this test).", "",
              "| T_e | D into D-T plasma | D into ⁶LiD plasma | Li-6 into ⁶LiD plasma | ion slowing distance in D-T |",
              "|---|---|---|---|---|"]
    e_ion = E.explosion_energy_per_ion_kev("6LiD", 15e-9)
    rows = []
    for te in (0.5, 1, 3, 10, 30, 100):
        dt = E.fuel_plasma("DT", 1.0, te, te)
        lid = E.fuel_plasma("6LiD", 1.0, te, te)
        a = E.beam_in_plasma(e_ion["D"], "D", dt)
        b = E.beam_in_plasma(e_ion["D"], "D", lid)
        c = E.beam_in_plasma(e_ion["Li6"], "Li6", lid)
        rows.append({"T_e": te, "D_in_DT": a["gain"], "D_in_LiD": b["gain"], "Li_in_LiD": c["gain"],
                     "path_m": a["path_m"]})
        lines.append(f"| {te:g} keV | {a['gain']:.2g} | {b['gain']:.2g} | {c['gain']:.2g} | {a['path_m'] * 1e6:.3g} µm |")
    lines.append("")
    data["hot"] = rows


def burn(lines: list[str], data: dict) -> None:
    lines += ["## 3. A ball packed with crystals, exploded at once", "",
              "Energy in = explosion energy of every crystal + heating the electrons to 1 keV. "
              "Fusion out = fast ions slowing down + thermal burn until the ball flies apart, "
              "after subtracting X-ray (bremsstrahlung) losses. Best crystal size and fill for each fuel.", "",
              "Two bounds for the fast-ion phase. *Generous*: the plasma is at its final temperature at once. "
              "*Conservative*: the ions slow down on the initial 1 keV electrons.", "",
              "| Fuel | Bound | Crystal diameter | Fill | Plasma temperature | Fusion / X-ray power "
              "| Ball radius for breakeven | Fuel mass | Energy in at breakeven | Intensity to strip the crystals |",
              "|---|---|---|---|---|---|---|---|---|---|"]
    radii = np.geomspace(1e-6, 1.0, 241)
    out = {}
    for fuel in ("6LiD", "7LiD", "6LiDT", "DT", "B10H14"):
        for bound, generous in (("generous", True), ("conservative", False)):
            best = None
            for r_nm in (1.5, 2.5, 4, 6, 10, 15, 25, 50):
                for fill in (0.3, 0.6, 1.0):
                    gains = np.array([E.array_burn(fuel, r_nm, fill, L, fast_in_final_plasma=generous)
                                      ["gain_after_radiation"] for L in radii])
                    if gains.max() < 1:
                        cand = {"breakeven": False, "max_gain": float(gains.max()), "r_nm": r_nm, "fill": fill,
                                "L": float(radii[int(gains.argmax())])}
                    else:
                        i = int(np.nonzero(gains >= 1)[0][0])
                        cand = {"breakeven": True, "r_nm": r_nm, "fill": fill, "L": float(radii[i])}
                    res = E.array_burn(fuel, r_nm, fill, cand["L"], fast_in_final_plasma=generous)
                    cand.update(res)
                    key = (cand["breakeven"], -cand["energy_in_J"] if cand["breakeven"] else cand["max_gain"])
                    if best is None or key > best[0]:
                        best = (key, cand)
            c = best[1]
            out[f"{fuel}-{bound}"] = c
            head = (f"| {fuel} | {bound} | {2 * c['r_nm']:g} nm | {c['fill']:.0%} | {c['T_eq_keV']:.0f} keV "
                    f"| {c['fusion_to_brem']:.2g} ")
            if c["breakeven"]:
                lines.append(head + f"| {c['ball_radius_m'] * 1e3:.3g} mm | {c['mass_g']:.3g} g "
                             f"| {c['energy_in_J']:.2g} J | {c['stripping_intensity_W_cm2']:.1g} W/cm² |")
            else:
                lines.append(head + f"| never (best gain {c['max_gain']:.2g}) | – | – "
                             f"| {c['stripping_intensity_W_cm2']:.1g} W/cm² |")
    lines.append("")
    data["burn"] = out
    for bound, generous in (("generous", True), ("conservative", False)):
        c = out[f"DT-{bound}"]
        lines += [f"D-T, {2 * c['r_nm']:g} nm crystals, {c['fill']:.0%} fill ({bound}): gain against ball radius", "",
                  "| Ball radius | Gain after X-ray losses | from fast ions | from thermal burn | Energy in |",
                  "|---|---|---|---|---|"]
        for L in (1e-5, 1e-4, 1e-3, 3e-3, 1e-2, 3e-2):
            r = E.array_burn("DT", c["r_nm"], c["fill"], L, fast_in_final_plasma=generous)
            lines.append(f"| {L * 1e3:g} mm | {r['gain_after_radiation']:.2g} | {r['fusion_fast_J'] / r['energy_in_J']:.2g} "
                         f"| {r['fusion_thermal_J'] / r['energy_in_J']:.2g} | {r['energy_in_J']:.2g} J |")
        lines.append("")


def spectrum_gain(fuel: str, radius_nm: float, catcher_material: str, catcher_targets: dict) -> dict:
    """Fusion energy / explosion energy when every ion of an exploding crystal stops in a cold catcher.

    Uniform-sphere Coulomb explosion: ion energies E = E_max (r/R)^2, so
    dN/dE ~ sqrt(E) up to E_max (0.9 of the ideal value, as in the MD runs).
    """
    parts, _ = E.FUELS[fuel]
    rho = E.charge_density(fuel)
    fusion = 0.0
    kinetic = 0.0
    for species, count in parts.items():
        z = E.SPECIES[species][1]
        e_max = 0.9 * z * P.E_CHARGE * rho * (radius_nm * 1e-9) ** 2 / (3 * P.EPS0) / P.J_PER_KEV
        e = np.linspace(e_max / 400, e_max, 400)
        w = np.sqrt(e)
        w /= w.sum()
        res = E.thick_target(e, species, catcher_material, catcher_targets)
        fusion += count * float(np.sum(w * res["fusion_MeV"]))
        kinetic += count * float(np.sum(w * e)) / 1e3
    return {"gain": fusion / kinetic, "proton_max_keV": 0.9 * P.E_CHARGE * rho * (radius_nm * 1e-9) ** 2
            / (3 * P.EPS0) / P.J_PER_KEV}


def proton_boron(lines: list[str], data: dict) -> None:
    lines += ["## 4. Proton–boron (p + ¹¹B → 3 α), the neutron-free fuel", "",
              "Cross section: Nevins & Swain form as refitted in arXiv:2601.00241 (peak 1.08 b at 612 keV). "
              "Fuel: decaborane-11 (B₁₀H₁₄) crystals; catcher: cold decaborane.", "",
              "| Crystal diameter | Fastest proton | Fusion energy / explosion energy (perfect laser, cold catcher) |",
              "|---|---|---|"]
    n_atoms = P._material("decaborane")[0]
    targets = {"B11": n_atoms[5], "p": n_atoms[1]}
    rows = []
    for r_nm in (5, 15, 25, 50, 100):
        g = spectrum_gain("B10H14", r_nm, "decaborane", targets)
        rows.append({"radius_nm": r_nm, **g})
        lines.append(f"| {2 * r_nm:g} nm | {g['proton_max_keV']:.3g} keV | {g['gain']:.2g} |")
    lines.append("")
    best_p = E.thick_target(np.geomspace(10, 20000, 300), "p", "decaborane", {"B11": n_atoms[5]})["gain"].max()
    lines += [f"Best single proton into cold decaborane: {best_p:.2g} of its energy back as fusion.", ""]
    lines += ["**Avalanche test.** Can fusion alphas knock protons into further fusions, enough for a chain "
              "reaction (k ≥ 1)? k = secondary fusions per fusion; the nuclear boost multiplies the Rutherford "
              "elastic cross section, as an upper bound for nuclear scattering.", "",
              "| Medium | k (Rutherford) | k (×5 nuclear boost) |", "|---|---|---|"]
    ava = {}
    cold1 = E.pb11_avalanche("decaborane", nuclear_boost=1)["k"]
    cold5 = E.pb11_avalanche("decaborane", nuclear_boost=5)["k"]
    ava["cold"] = (cold1, cold5)
    lines.append(f"| cold solid decaborane | {cold1:.2g} | {cold5:.2g} |")
    for t in (30, 100, 300):
        pl = E.Plasma({"p": n_atoms[1], "B11": n_atoms[5]}, t_e=t, t_i=t)
        k1 = E.pb11_avalanche(hot_plasma=pl, nuclear_boost=1)["k"]
        k5 = E.pb11_avalanche(hot_plasma=pl, nuclear_boost=5)["k"]
        ava[f"hot{t}"] = (k1, k5)
        lines.append(f"| hot solid-density p-B plasma, T = {t} keV | {k1:.2g} | {k5:.2g} |")
    lines.append("")
    lines += ["A whole ball of exploded decaborane crystals is in the burn table of section 3 (fuel B10H14).", ""]
    data["proton_boron"] = {"crystals": rows, "best_single_proton": float(best_p), "avalanche": ava}


def screening(lines: list[str], data: dict) -> None:
    lines += ["## 5. Electron screening (\"lattice-enhanced\" fusion)", "",
              "Bound electrons in a solid lower the Coulomb barrier by a screening potential U_e "
              "(≈ 0.15 keV for atoms; up to ≈ 1 keV claimed in some metals). Deuterons into a cold tritide catcher:", "",
              "| U_e | gain at 3 keV | at 30 keV | at 300 keV | best over all energies |", "|---|---|---|---|---|"]
    n = P.lid_molecule_density()
    e0 = np.array([3.0, 30.0, 300.0])
    rows = {}
    for ue in (0, 150, 1000, 3000):
        g = E.thick_target(e0, "D", "LiD", {"T": n}, screening_ev=ue)["gain"]
        best = E.thick_target(np.geomspace(1, 3000, 300), "D", "LiD", {"T": n}, screening_ev=ue)["gain"].max()
        rows[ue] = [float(x) for x in g] + [float(best)]
        lines.append(f"| {ue} eV | {g[0]:.2g} | {g[1]:.2g} | {g[2]:.2g} | {best:.3g} |")
    lines.append("")
    data["screening"] = rows


def reproduce_cluster_experiment(lines: list[str], data: dict) -> None:
    """Texas Petawatt deuterium-cluster record (Bang et al., PRL 111, 055002; PRE 87, 023106)."""
    lines += ["## 6. Check against a real experiment: the Texas Petawatt cluster-fusion record", "",
              "Measured: 1.9×10⁷ neutrons per shot from 90–180 J, deuterium clusters, ion temperature kT = 18 keV "
              "(mean ion energy ≈ 14 keV), average atom density 10¹⁸–10¹⁹ cm⁻³, hot volume 1–10 mm³. Same model as "
              "the experimenters: Y = ½ n² ⟨σv⟩ V t, burn time t = plasma radius / mean ion speed, half of D+D "
              "reactions give a neutron.", "",
              "| Density (cm⁻³) | Volume | Predicted neutrons | Measured |", "|---|---|---|---|"]
    t_kev = 18.0
    sv = E.reactivity("D", "D", t_kev)
    v_mean = math.sqrt(8 * t_kev * P.J_PER_KEV / (math.pi * P.MASS_U["D"] * P.AMU))
    rows = []
    for n_cm3 in (1e18, 3e18, 1e19):
        for vol_mm3 in (1, 10):
            vol = vol_mm3 * 1e-9
            r0 = (3 * vol / (4 * math.pi)) ** (1 / 3)
            n = n_cm3 * 1e6
            neutrons = 0.5 * 0.5 * n * n * sv * vol * (r0 / v_mean)
            rows.append({"n_cm3": n_cm3, "volume_mm3": vol_mm3, "neutrons": neutrons})
            lines.append(f"| {n_cm3:.0e} | {vol_mm3} mm³ | {neutrons:.2g} | 1.9×10⁷ |")
    # Cluster size that the Coulomb-explosion model needs for a 14 keV mean ion energy
    rho = E.charge_density("D2")
    per_r2 = 0.9 * 0.6 * P.E_CHARGE * rho / (3 * P.EPS0) / P.J_PER_KEV
    r_needed = math.sqrt(14.0 / per_r2)
    reactions = 2 * 1.9e7
    ratio = reactions * P.Q_VALUE_MEV["DD"] * 1.602176634e-13 / 135.0
    lines += ["", f"A 14 keV mean ion energy corresponds to Coulomb-exploding clusters of {2 * r_needed * 1e9:.0f} nm diameter "
              f"in this model (the experimenters inferred clusters of roughly this size).",
              f"Measured fusion energy / laser energy: ≈ {ratio:.1g}, the same order as this project's estimates "
              f"(10⁻⁸–10⁻⁷ for ions stopping in cold matter).", ""]
    data["texas_petawatt"] = {"rows": rows, "cluster_diameter_nm": 2 * r_needed * 1e9, "energy_ratio": ratio}


def ignitor_deuteron_energy(range_g_cm2: float = 1.2) -> tuple[float, dict]:
    """Deuteron energy (keV) whose range in 400 g/cm3, 5 keV D-T is range_g_cm2; also the range table."""
    rho_g_cm3 = 400.0
    n_dt = rho_g_cm3 * 1e3 / (2.5 * P.AMU)  # ions per m^3 (half D, half T)
    plasma = E.Plasma({"D": n_dt / 2, "T": n_dt / 2}, t_e=5.0, t_i=5.0)
    ranges = {e: E.beam_in_plasma(e, "D", plasma)["path_m"] * rho_g_cm3 * 100 for e in (2e3, 5e3, 1e4, 2e4, 4e4)}
    return float(np.interp(range_g_cm2, list(ranges.values()), list(ranges.keys()))), ranges


def spot_fraction_isotropic(distance_m: float, spot_radius_m: float = 10e-6) -> float:
    """Share of an isotropic source's ions that head into a spot: the solid-angle fraction (1 - cos theta) / 2."""
    theta = math.atan(spot_radius_m / distance_m)
    return (1 - math.cos(theta)) / 2


def fast_ignition(lines: list[str], data: dict) -> None:
    """Could exploding crystals be the ion 'spark' for fast ignition of compressed D-T?"""
    lines += ["## 7. Exploding crystals as the ignition spark for compressed fuel (ion fast ignition)", "",
              "The one published route by which laser-made ions could give net energy: compress D-T to ~400 g/cm³ "
              "with a main driver, then deliver ~10 kJ of ions into a ~20 µm spot within ~20 ps "
              "(Fernández et al., Nucl. Fusion 54, 054006). Designs of this type promise target gains ~100 "
              "(e.g. 250 kJ compression + 80–100 kJ ignitor). Can a layer of Coulomb-exploding LiD crystals be "
              "the ion source?", ""]
    # Deuteron energy whose stopping range in 400 g/cm3, 5 keV D-T is about 1.2 g/cm2
    e_needed, ranges = ignitor_deuteron_energy()
    lines += ["Deuteron stopping range in the compressed fuel (400 g/cm³, 5 keV): "
              + ", ".join(f"{e / 1e3:g} MeV → {r:.2g} g/cm²" for e, r in ranges.items())
              + f". A ~1.2 g/cm² range needs deuterons of about {e_needed / 1e3:.0f} MeV.", ""]
    # Crystal radius whose explosion reaches that energy at the top of its spectrum
    rho_q = E.charge_density("6LiD")
    per_r2 = 0.9 * P.E_CHARGE * rho_q / (3 * P.EPS0) / P.J_PER_KEV
    r = math.sqrt(e_needed / per_r2)
    # Fraction of all explosion energy carried by deuterons within +-25% of e_needed
    e = np.linspace(e_needed / 400, e_needed, 400)
    w = np.sqrt(e)
    d_energy = w * e
    window = (e >= 0.75 * e_needed)
    li_energy = 3 * float(np.sum(w * e))  # Li3+ ions: same count, three times the energy at each radius
    frac_window = float(np.sum(d_energy[window]) / (np.sum(d_energy) + li_energy))
    e_strip = rho_q * r / (3 * P.EPS0)
    pulse = P.Pulse(e0=e_strip, wavelength=0.8e-6, fwhm=6e-15)
    eff = coulomb_energy_J("6LiD", r) / (pulse.fluence() * math.pi * r * r)
    v_lo = math.sqrt(2 * 0.75 * e_needed * P.J_PER_KEV / (P.MASS_U["D"] * P.AMU))
    v_hi = math.sqrt(2 * e_needed * P.J_PER_KEV / (P.MASS_U["D"] * P.AMU))
    lines += [f"Crystal needed: about {2 * r * 1e9:.0f} nm diameter (stripping intensity "
              f"{P.intensity_w_cm2(e_strip):.1g} W/cm²). Of its explosion energy, {frac_window:.0%} is carried by "
              f"deuterons in the useful band (within 25% of {e_needed / 1e3:.0f} MeV); laser → explosion energy "
              f"efficiency ≈ {min(eff, 1):.0%}.", "",
              "| Source distance from the 20 µm spot | Fraction of ions aimed at the spot (isotropic explosion) "
              "| Arrival spread of the useful deuterons | Ignitor laser energy for 10 kJ in the spot "
              "| Target gain (30 MJ yield ÷ 250 kJ compression + ignitor) |",
              "|---|---|---|---|---|"]
    rows = []
    for d_um in (30, 100, 1000):
        d = d_um * 1e-6
        geom = spot_fraction_isotropic(d)
        spread_ps = d * (1 / v_lo - 1 / v_hi) * 1e12
        laser = 10e3 / (min(eff, 1) * frac_window * geom)
        gain = 30e6 / (250e3 + laser)
        rows.append({"distance_um": d_um, "geometric_fraction": geom, "spread_ps": spread_ps, "laser_J": laser,
                     "target_gain": gain})
        lines.append(f"| {d_um} µm | {geom:.2g} | {spread_ps:.2g} ps | {laser:.2g} J | {gain:.2g} |")
    lines += ["", "An isotropic explosion sends almost all ions away from the spot: the spot covers only the solid-angle "
              "fraction (1 − cos θ)/2 of the sky seen from the crystals. Even with the crystals only 30 µm from the "
              "compressed core (not buildable: they would be swept into the imploding fuel and sit behind dense "
              "plasma the ignitor laser cannot cross), the ignitor needs ~20 MJ of laser and the target gain is "
              "~1.4. At 0.1–1 mm it needs 0.2–20 GJ and the target loses energy. Fast-ignition designs instead "
              "assume directional, foil-accelerated ion beams with ~10% efficiency (~100 kJ of ignitor laser), and "
              "those are still undemonstrated. Shaped crystals are tested in tools/shaped_crystals.py.", ""]
    data["fast_ignition"] = {"ranges_g_cm2": ranges, "deuteron_MeV": e_needed / 1e3, "crystal_diameter_nm": 2 * r * 1e9,
                             "useful_fraction": frac_window, "laser_efficiency": eff, "rows": rows}


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--out", default="docs/results/energy")
    args = ap.parse_args()
    out = Path(args.out)
    out.mkdir(parents=True, exist_ok=True)
    lines = ["# Energy experiments (generated by tools/energy_experiments.py)", ""]
    data: dict = {}
    cold_catcher(lines, data)
    hot_catcher(lines, data)
    burn(lines, data)
    proton_boron(lines, data)
    screening(lines, data)
    reproduce_cluster_experiment(lines, data)
    fast_ignition(lines, data)
    (out / "summary.md").write_text("\n".join(lines) + "\n")
    (out / "data.json").write_text(json.dumps(data, indent=1, default=float))
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
