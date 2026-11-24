"""Shaped crystals: can a directional Coulomb explosion rescue ion fast ignition?

Fast ignition needs ~10 kJ of ~4 MeV deuterons in a 20 um spot on compressed
D-T (energy_experiments.py, section 7). A spherical crystal explodes
isotropically, so almost all its ions miss. This script runs Coulomb
explosions of differently shaped crystals in the improved model and measures
where their deuterons go:

* a sphere (reference);
* LiD disks from thick to a single atomic layer;
* double-layer disks: a lithium-only disk with a deuterium "dot" on one face,
  the asymmetric design used for collimated laser-ion beams
  (Esirkepov et al., PRL 89, 175003 (2002)).

Coulomb explosions at fixed shape and density are self-similar: a crystal s
times larger gives the same directions with s^2 times the energy, s^5 times
the total energy, and needs s times the stripping field. So each run is scaled
to the size whose deuterons have the ~4 MeV the ignitor needs, and the laser
energy for 10 kJ in the spot is worked out for sources 30 um to 1 mm away.

    python3 tools/shaped_crystals.py [--analyse-only] [--force]
"""

from __future__ import annotations

import argparse
import json
import math
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))
sys.path.insert(0, str(Path(__file__).resolve().parent))

from energy_experiments import ignitor_deuteron_energy, spot_fraction_isotropic  # noqa: E402

from nanofusion.improved import physics as P  # noqa: E402
from nanofusion.improved.native import coulomb_step  # noqa: E402
from nanofusion.improved.simulation import ImprovedProperties, build_crystal, run_improved  # noqa: E402

OUT = Path("docs/results/shaped")
DISTANCES_UM = (30, 100, 300, 1000)
SPOT_M = 10e-6  # spot radius
IGNITOR_J = 10e3
YIELD_J = 30e6
COMPRESSION_J = 250e3
PULSE_FS = 6.0  # few-cycle pulse: the most generous laser-to-explosion efficiency
SKIN_NM = 11.0  # laser penetration depth into stripped LiD
LARGEST_LASER_W = 1e16  # ~10 PW, the most powerful lasers today

RUNS = {
    "sphere-R2.5": dict(shape="sphere", size=2.5),
    "disk-a3.4-t1.7": dict(shape="disk", size=3.4, thickness_nm=1.7),
    "disk-a5-t1": dict(shape="disk", size=5.0, thickness_nm=1.0),
    "disk-a7-t0.6": dict(shape="disk", size=7.0, thickness_nm=0.6),
    "disk-a10-t0.2": dict(shape="disk", size=10.0, thickness_nm=0.2),
    "dot-a6-t1-d1.5": dict(shape="disk", size=6.0, thickness_nm=1.0, dot_nm=1.5, lithium="Li7"),
    "dot-a6-t1-d2.5": dict(shape="disk", size=6.0, thickness_nm=1.0, dot_nm=2.5, lithium="Li7"),
    "dot-a6-t1-d4": dict(shape="disk", size=6.0, thickness_nm=1.0, dot_nm=4.0, lithium="Li7"),
    "dot-a6-t1-d6": dict(shape="disk", size=6.0, thickness_nm=1.0, dot_nm=6.0, lithium="Li7"),
    "dot-a6-t0.6-d2.5": dict(shape="disk", size=6.0, thickness_nm=0.6, dot_nm=2.5, lithium="Li7"),
    "dot-a6-t2-d2.5": dict(shape="disk", size=6.0, thickness_nm=2.0, dot_nm=2.5, lithium="Li7"),
}


def strip_field(props: ImprovedProperties) -> float:
    """Largest Coulomb field on any ion of the fully stripped crystal (V/m): what the laser must beat."""
    pos, is_li = build_crystal(props)
    q = np.where(is_li, 3.0, 1.0) * P.E_CHARGE
    mass = np.where(is_li, P.MASS_U[props.lithium], P.MASS_U["D"]) * P.AMU
    f, _, _, _ = coulomb_step(pos, np.zeros_like(pos), q, mass, is_li.astype(np.int32), (1e-12) ** 2,
                              (1e-10) ** 2, 1.0, 0.0, 1.0)
    return float(np.max(np.linalg.norm(f, axis=1) / q))


def run_one(name: str, spec: dict) -> dict:
    props = ImprovedProperties(name=name, intensity_w_cm2=0.0, duration_fs=200.0, frames=3, **spec)
    t0 = time.monotonic()
    _, rep = run_improved(props)
    fs = rep.final_state
    v, li, m = fs["velocities"], fs["is_lithium"], fs["mass"]
    ke = 0.5 * m * np.einsum("ij,ij->i", v, v) / P.J_PER_KEV
    d = ~li
    angle = np.degrees(np.arccos(np.clip(v[d, 2] / np.linalg.norm(v[d], axis=1), -1, 1)))
    if props.shape == "sphere":
        radius = (3 * len(v) * P.SPACING**3 / (4 * math.pi)) ** (1 / 3)
        thickness = 2 * radius
    else:
        radius = props.size * 1e-9
        pos0, _ = build_crystal(props)
        thickness = float(np.ptp(pos0[:, 2])) + P.SPACING
    return {
        "inputs": spec, "ions": len(v), "deuterons": int(d.sum()), "steps": rep.steps,
        "wall_seconds": time.monotonic() - t0,
        "energy_check_keV": rep.dynamics["energy_check_ion_KE_minus_potential_drop_keV"],
        "total_keV": float(ke.sum()), "radius_m": radius, "thickness_m": thickness,
        "strip_field_V_m": strip_field(props),
        "d_energy_keV": ke[d].round(6).tolist(), "d_angle_deg": angle.round(4).tolist(),
    }


def run_all(force: bool) -> None:
    OUT.mkdir(parents=True, exist_ok=True)
    for name, spec in RUNS.items():
        path = OUT / f"{name}.json"
        if path.exists() and not force:
            continue
        rec = run_one(name, spec)
        path.write_text(json.dumps(rec))
        print(f"{name}: {rec['ions']} ions, {rec['wall_seconds']:.0f} s", flush=True)


# --- analysis -------------------------------------------------------------------

def design(rec: dict, distance_m: float, e_needed: float) -> dict:
    """Best scaling of one crystal shape as an ignitor source at the given distance.

    Useful energy: deuterons within +-25% of e_needed whose direction lies inside
    the cone from the crystal to the spot (the crystal's +z axis points at it).
    Laser energy per crystal: a few-cycle pulse whose peak field just strips the
    crystal, over the crystal's face.
    """
    e = np.array(rec["d_energy_keV"])
    ang = np.array(rec["d_angle_deg"])
    theta = math.degrees(math.atan(SPOT_M / distance_m))
    best = None
    for s in np.geomspace(1, 1e3, 600):
        es = e * s**2
        band = (es >= 0.75 * e_needed) & (es <= 1.25 * e_needed)
        # Narrow cones hold few simulated ions. Below 10, use the narrowest cone that holds 10 and
        # scale by solid angle (the ion density near the axis is taken as even).
        cone, factor = theta, 1.0
        for ref in (theta, 2.0, 5.0, 10.0):  # wider cones would say nothing about the axis
            if ref < theta:
                continue
            if np.count_nonzero(band & (ang <= ref)) >= 10:
                cone = ref
                factor = (1 - math.cos(math.radians(theta))) / (1 - math.cos(math.radians(ref)))
                break
        else:
            continue
        # Each simulated ion stands for s^3 real ions.
        useful = float(np.sum(es[band & (ang <= cone)])) * factor * s**3
        if useful <= 0:
            continue
        total = rec["total_keV"] * s**5
        field = rec["strip_field_V_m"] * s
        pulse = P.Pulse(e0=field, wavelength=0.8e-6, fwhm=PULSE_FS * 1e-15)
        laser = pulse.fluence() * math.pi * (rec["radius_m"] * s) ** 2 / P.J_PER_KEV  # keV
        efficiency = min(1.0, total / laser)  # laser -> explosion energy
        fraction = useful / total  # explosion energy that is useful
        ignitor = IGNITOR_J / (efficiency * fraction)
        cand = {"scale": s, "fraction_useful": fraction, "laser_to_explosion": efficiency,
                "cone_extrapolated_from_deg": cone if cone > theta else None,
                "ignitor_laser_J": ignitor, "peak_power_W": ignitor / (PULSE_FS * 1e-15),
                "target_gain": YIELD_J / (COMPRESSION_J + ignitor),
                "radius_nm": rec["radius_m"] * s * 1e9, "thickness_nm": rec["thickness_m"] * s * 1e9,
                "strip_intensity_W_cm2": P.intensity_w_cm2(field)}
        if best is None or cand["ignitor_laser_J"] < best["ignitor_laser_J"]:
            best = cand
    return best or {}


def cone_table(rec: dict) -> dict:
    e = np.array(rec["d_energy_keV"])
    ang = np.array(rec["d_angle_deg"])
    out = {}
    for th in (5, 10, 20):
        sel = ang <= th
        out[th] = {"count": float(sel.mean()), "energy_of_all": float(e[sel].sum() / rec["total_keV"])}
    return out


def analyse() -> str:
    e_needed, _ = ignitor_deuteron_energy()
    recs = {n: json.loads((OUT / f"{n}.json").read_text()) for n in RUNS if (OUT / f"{n}.json").exists()}
    out = ["# Shaped crystals as a directional ion source for fast ignition", "",
           f"Ignitor: {IGNITOR_J / 1e3:g} kJ of deuterons within ±25% of {e_needed / 1e3:.1f} MeV into a "
           f"{2 * SPOT_M * 1e6:g} µm spot; target gain = {YIELD_J / 1e6:g} MJ ÷ ({COMPRESSION_J / 1e3:g} kJ compression "
           f"+ ignitor laser). Laser: {PULSE_FS:g} fs pulses just strong enough to strip each crystal. "
           "Produced by `tools/shaped_crystals.py`.", "",
           "## Where the deuterons go", "",
           "Shares of deuterons (count) and of the whole explosion energy carried by deuterons, within a cone "
           "around the +z axis (disks: their axis; dot targets: the dot side). An isotropic source puts "
           "(1 − cos θ)/2 of its ions in a cone of half-angle θ: 0.2% at 5°, 0.8% at 10°, 3% at 20°.", "",
           "| crystal | ions | D energy spread (10–90%, keV) | D within 5° | 10° | 20° "
           "| explosion energy in D within 10° | 20° | energy check (keV) |",
           "|---|---|---|---|---|---|---|---|---|"]
    for name, rec in recs.items():
        c = cone_table(rec)
        e = np.array(rec["d_energy_keV"])
        lo, hi = np.percentile(e, [10, 90])
        out.append(f"| {name} | {rec['ions']} | {lo:.3g}–{hi:.3g} | {c[5]['count']:.3f} | {c[10]['count']:.3f} "
                   f"| {c[20]['count']:.3f} | {c[10]['energy_of_all']:.4f} | {c[20]['energy_of_all']:.4f} "
                   f"| {rec['energy_check_keV']:.2g} of {rec['total_keV']:.3g} |")
    out.append("")
    out += ["## As a fast-ignition source", "",
            "Each crystal is scaled to the size that needs the least laser energy for the ignitor at that distance. "
            "Isotropic reference with the textbook solid-angle fraction: "
            + ", ".join(f"{d} µm → {spot_fraction_isotropic(d * 1e-6):.2g}" for d in DISTANCES_UM) + ".", "",
            "| crystal | source distance | useful share of explosion energy | laser → explosion | ignitor laser "
            "| peak power | target gain | crystal size at that scale (radius × thickness) | stripping intensity "
            "| thin enough for the laser to strip? |",
            "|---|---|---|---|---|---|---|---|---|---|"]
    summary = {}
    for name, rec in recs.items():
        for d_um in DISTANCES_UM:
            b = design(rec, d_um * 1e-6, e_needed)
            if not b:
                out.append(f"| {name} | {d_um} µm | none: no deuterons go this way | – | – | – | – | – | – | – |")
                continue
            summary.setdefault(name, {})[d_um] = b
            ok = "yes" if b["thickness_nm"] <= SKIN_NM else "no"
            note = f" (from the {b['cone_extrapolated_from_deg']:g}° cone)" if b["cone_extrapolated_from_deg"] else ""
            out.append(f"| {name} | {d_um} µm | {b['fraction_useful']:.2g}{note} | {b['laser_to_explosion']:.2g} "
                       f"| {b['ignitor_laser_J']:.2g} J | {b['peak_power_W']:.2g} W | {b['target_gain']:.2g} "
                       f"| {b['radius_nm']:.0f} × {b['thickness_nm']:.0f} nm | {b['strip_intensity_W_cm2']:.1g} W/cm² "
                       f"| {ok} |")
    out.append("")
    best = {}
    for d_um in DISTANCES_UM:
        cands = [(v[d_um]["target_gain"], n) for n, v in summary.items() if d_um in v]
        if cands:
            g, n = max(cands)
            best[d_um] = (n, summary[n][d_um])
    if best:
        out += ["## Best design at each distance", "",
                "| distance | crystal | ignitor laser | peak power (today's largest lasers: ~10¹⁶ W) | target gain |",
                "|---|---|---|---|---|"]
        for d_um, (n, b) in best.items():
            out.append(f"| {d_um} µm | {n} | {b['ignitor_laser_J']:.2g} J | {b['peak_power_W']:.2g} W "
                       f"({b['peak_power_W'] / LARGEST_LASER_W:.0g}× today's largest) | {b['target_gain']:.2g} |")
        out.append("")
    return "\n".join(out) + "\n"


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--analyse-only", action="store_true")
    ap.add_argument("--force", action="store_true")
    args = ap.parse_args()
    if not args.analyse_only:
        run_all(args.force)
    text = analyse()
    (OUT / "summary.md").write_text(text)
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
