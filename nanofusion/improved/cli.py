"""``python -m nanofusion improved``: run the improved model and print a report."""

from __future__ import annotations

import argparse
import json
import sys
from dataclasses import asdict

from .simulation import ImprovedProperties, ImprovedReport, run_improved


def build_parser() -> argparse.ArgumentParser:
    d = ImprovedProperties()
    p = argparse.ArgumentParser(
        prog="nanofusion improved",
        description="Laser-stripped LiD nanocrystal with all-pairs Coulomb forces, a real laser pulse "
                    "and cross-section fusion yields.",
    )
    p.add_argument("--shape", choices=["sphere", "cube", "disk"], default=d.shape)
    p.add_argument("--size", type=float, default=d.size, help="sphere radius in nm, or cube edge in ions")
    p.add_argument("--lithium", choices=["Li6", "Li7"], default=d.lithium)
    p.add_argument("--intensity", type=float, default=d.intensity_w_cm2, help="peak intensity, W/cm^2")
    p.add_argument("--wavelength", type=float, default=d.wavelength_um, help="µm")
    p.add_argument("--fwhm", type=float, default=d.fwhm_fs, help="pulse intensity FWHM, fs")
    p.add_argument("--phase", type=float, default=d.phase_deg, help="carrier phase at the pulse centre, degrees")
    p.add_argument("--centre", type=float, default=d.centre_fs, help="pulse centre relative to stripping, fs")
    p.add_argument("--duration", type=float, default=d.duration_fs, help="simulated time, fs")
    p.add_argument("--target-radius", type=float, default=d.target_radius_nm,
                   help="crystal radius (nm) to extrapolate yields to; 15 = the proposal's 30 nm crystal")
    p.add_argument("--profile", choices=["uniform", "expanded"], default=d.profile,
                   help="uniform crystal, or expanded: pre-expanded by a prepulse (forms shock shells)")
    p.add_argument("--core-density", type=float, default=d.core_density,
                   help="expanded profile: central density as a fraction of solid LiD")
    p.add_argument("--expansion", type=float, default=d.expansion,
                   help="expanded profile: edge radius / original edge radius")
    p.add_argument("--geometry", choices=["single", "pair"], default=d.geometry,
                   help="one crystal, or two side by side along the laser polarisation")
    p.add_argument("--gap", type=float, default=d.gap_nm, help="pair: gap between the crystals, nm")
    p.add_argument("--thickness", type=float, default=d.thickness_nm, help="disk: thickness along its axis, nm")
    p.add_argument("--dot", type=float, default=d.dot_nm,
                   help="disk: radius of a deuterium layer on one face of a lithium-only disk, nm (0 = LiD disk)")
    p.add_argument("-o", "--output", help="write the full report as JSON")
    return p


def format_report(r: ImprovedReport) -> str:
    L, V, D, F, E, S = r.laser, r.validity, r.dynamics, r.fusion, r.energy, r.scaled
    yes = lambda b: "yes" if b else "NO"  # noqa: E731
    lines = [
        f"{r.particles:,} ions ({r.deuterons:,} D), equivalent radius {r.radius_nm:.2f} nm; "
        f"{r.simulated_fs:.1f} fs in {r.steps:,} steps, {r.wall_seconds:.1f} s"
        f"{'' if r.native_kernel else ' (numpy fallback)'}",
        "",
        "Laser",
        f"  peak field {L['peak_field_V_per_m']:.3g} V/m = {L['intensity_W_per_cm2']:.3g} W/cm², "
        f"period {L['period_fs']:.1f} fs, {L['cycles_in_fwhm']:.2f} cycles in FWHM, electron a0 = {L['electron_a0']:.2g}",
        "Validity",
        f"  field needed to strip this crystal {V['field_needed_to_strip_V_per_m']:.3g} V/m -> "
        f"full stripping plausible: {yes(V['full_stripping_plausible'])}",
        f"  for the {S['target_radius_nm'] * 2:g} nm target crystal: {V['field_needed_to_strip_target_V_per_m']:.3g} V/m "
        f"(= {V['intensity_needed_to_strip_target_W_per_cm2']:.2g} W/cm²) -> plausible: "
        f"{yes(V['full_stripping_plausible_at_target_size'])}",
        f"  stripped LiD is {V['electron_density_over_critical']:.0f}x critical density; skin depth "
        f"{V['skin_depth_nm']:.1f} nm",
        "Dynamics",
        f"  D: max {D['max_D_keV']:.3g} keV, mean {D['mean_D_keV']:.3g} keV "
        f"(uniform-sphere Coulomb explosion: {D['coulomb_explosion_theory_max_D_keV']:.3g} keV)",
        f"  Li: max {D['max_Li_keV']:.3g} keV, mean {D['mean_Li_keV']:.3g} keV "
        f"(theory {D['coulomb_explosion_theory_max_Li_keV']:.3g} keV)",
        f"  closest D-Li approach {D['closest_D_Li_approach_A']:.3g} Å at {D['closest_D_Li_approach_time_fs']:.2f} fs; "
        f"highest D-Li pair energy {D['max_D_Li_pair_energy_keV']:.3g} keV",
        f"  energy check: ion KE - potential drop = {D['energy_check_ion_KE_minus_potential_drop_keV']:.4g} keV, "
        f"laser work = {D['laser_work_on_ions_keV']:.4g} keV",
        "Expected fusion reactions (this crystal, one pulse)",
        f"  inside the crystal:  D+Li {F['in_crystal_D_Li']:.3g}   D+D {F['in_crystal_D_D']:.3g}",
        f"  ions stopping in surrounding cold LiD:  D+Li {F['beam_target_D_Li']:.3g}   D+D {F['beam_target_D_D']:.3g}"
        f"   Li+D {F['beam_target_Li_D']:.3g}",
        f"  ... in a tritide (LiT) catcher instead: {F['catchers']['LiT']['total']:.3g} "
        f"({', '.join(f'{k} {v:.2g}' for k, v in F['catchers']['LiT']['channels'].items())})",
        f"  total {F['total']:.3g}; fusion energy {F['fusion_energy_MeV']:.3g} MeV vs laser energy through the "
        f"crystal {E['laser_energy_through_crystal_MeV']:.3g} MeV (gain {_g(E['gain'])})",
        f"Extrapolated to a {S['target_radius_nm'] * 2:g} nm crystal (energies x R², ions x R³; assumes full stripping)",
        f"  {S['deuterons']:,} D, max D energy {S['max_D_keV']:.3g} keV; LiD catcher: D+Li {S['beam_target_D_Li']:.3g}, "
        f"D+D {S['beam_target_D_D']:.3g}, Li+D {S['beam_target_Li_D']:.3g}; gain {_g(S['gain'])}",
        "Crystal size table (reactions per crystal per pulse; gain at the lowest intensity that strips it)",
        "  radius   intensity to strip   max D      LiD catcher (gain)        LiT catcher (gain)",
    ] + [
        f"  {row['radius_nm']:5.0f} nm  {row['intensity_to_strip_W_per_cm2']:9.2g} W/cm²  {row['max_D_keV']:8.3g} keV"
        f"  {row['LiD']['total']:9.3g} ({_g(row['LiD']['gain_at_strip_threshold'])})"
        f"  {row['LiT']['total']:12.3g} ({_g(row['LiT']['gain_at_strip_threshold'])})"
        for row in r.size_table
    ]
    return "\n".join(lines)


def _g(x) -> str:
    return "n/a" if x is None else f"{x:.2g}"


def main(argv: list[str] | None = None) -> int:
    a = build_parser().parse_args(argv)
    props = ImprovedProperties(shape=a.shape, size=a.size, lithium=a.lithium, intensity_w_cm2=a.intensity,
                               wavelength_um=a.wavelength, fwhm_fs=a.fwhm, phase_deg=a.phase,
                               centre_fs=a.centre, duration_fs=a.duration, target_radius_nm=a.target_radius,
                               profile=a.profile, core_density=a.core_density, expansion=a.expansion, geometry=a.geometry, gap_nm=a.gap,
                               thickness_nm=a.thickness, dot_nm=a.dot)
    try:
        props.validate()
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2
    interactive = sys.stderr.isatty()

    def progress(done: int, total: int) -> None:
        if interactive:
            print(f"\r{done / total:5.1%}", end="", file=sys.stderr, flush=True)

    _, report = run_improved(props, progress=progress)
    if interactive:
        print(file=sys.stderr)
    print(format_report(report))
    if a.output:
        with open(a.output, "w", encoding="utf-8") as f:
            json.dump(asdict(report), f, indent=2)
        print(f"saved {a.output}")
    return 0
