"""Parameter scan of the improved model: can a laser-stripped LiD crystal fuse?

    python3 tools/scan_improved.py [--quick] [--out docs/results]

Writes one JSON report per case plus summary.md to the output directory.
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from dataclasses import asdict, replace
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import ImprovedProperties, run_improved  # noqa: E402

BASE = ImprovedProperties(shape="sphere", size=2.0, lithium="Li6", wavelength_um=3.6, fwhm_fs=6.0,
                          duration_fs=40.0, frames=10, target_radius_nm=15.0)


def cases(quick: bool) -> list[tuple[str, ImprovedProperties]]:
    out = []
    intensities = [1e17, 1e18, 1e19, 1e20, 1e21, 1e22]
    for i in intensities:  # A: the guided mechanism vs intensity, both stripping phases
        for phase in (0.0, 90.0):
            out.append((f"A-I{i:.0e}-ph{phase:.0f}", replace(BASE, intensity_w_cm2=i, phase_deg=phase)))
    sizes = [1.0, 1.5, 2.0, 3.0] if not quick else [1.0, 2.0]
    for r in sizes:  # B: crystal size at an intensity that can strip a few-nm crystal
        out.append((f"B-R{r:g}-I1e20", replace(BASE, size=r, intensity_w_cm2=1e20)))
    for i in (1e19, 1e20, 1e21, 1e22):  # D: Li-7, the only isotope with a q/m difference to "guide" D into Li
        for phase in (0.0, 90.0):
            out.append((f"D-Li7-I{i:.0e}-ph{phase:.0f}", replace(BASE, lithium="Li7", intensity_w_cm2=i,
                                                                 phase_deg=phase)))
    for name in ("D-Li7-I1e+22-ph90", "A-I1e+20-ph0"):  # E: same cases with 4x smaller steps
        src = dict(out)[name]
        out.append((f"E-fine-{name}", replace(src, max_step_fs=src.max_step_fs / 4)))
    out.append(("C-0.8um-I1e20", replace(BASE, wavelength_um=0.8, intensity_w_cm2=1e20)))  # Ti:sapphire
    out.append(("C-0.8um-I1e22", replace(BASE, wavelength_um=0.8, intensity_w_cm2=1e22)))
    if quick:
        out = [c for c in out if not c[0].startswith("A") or "ph0" in c[0]]
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--quick", action="store_true")
    ap.add_argument("--out", default="docs/results")
    ap.add_argument("--only", help="run only cases whose name contains this text")
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, props in cases(args.quick):
        if args.only and args.only not in name:
            continue
        path = out_dir / f"{name}.json"
        if path.exists():
            rep = json.loads(path.read_text())
        else:
            t0 = time.monotonic()
            _, report = run_improved(replace(props, name=name))
            rep = asdict(report)
            path.write_text(json.dumps(rep, indent=1))
            print(f"{name}: {report.particles} ions, {report.steps} steps, {time.monotonic() - t0:.0f} s", flush=True)
        rows.append((name, rep))
    lines = [
        "| case | ions | I (W/cm²) | λ (µm) | phase | closest D–Li (Å) | max D–Li pair (keV) | max D (keV) "
        "| in-crystal D+Li | in-crystal D+D | beam-target D+Li | beam-target D+D | gain | 30 nm: max D (keV) "
        "| 30 nm: reactions | 30 nm: gain | strip 30 nm? |",
        "|---" * 17 + "|",
    ]
    for name, r in rows:
        i, d, f, s, v = r["inputs"], r["dynamics"], r["fusion"], r["scaled"], r["validity"]
        g = r["energy"]["gain"]
        lines.append(
            f"| {name} | {r['particles']} | {i['intensity_w_cm2']:.0e} | {i['wavelength_um']} | {i['phase_deg']:.0f}° "
            f"| {d['closest_D_Li_approach_A']:.3g} | {d['max_D_Li_pair_energy_keV']:.3g} | {d['max_D_keV']:.3g} "
            f"| {f['in_crystal_D_Li']:.2g} | {f['in_crystal_D_D']:.2g} | {f['beam_target_D_Li']:.2g} "
            f"| {f['beam_target_D_D']:.2g} | {g:.1g} | {s['max_D_keV']:.3g} "
            f"| {s['beam_target_D_Li'] + s['beam_target_D_D']:.2g} | {s['gain']:.1g} "
            f"| {'yes' if v['full_stripping_plausible_at_target_size'] else 'no'} |"
            if g is not None and s["gain"] is not None else f"| {name} | (laser off) |"
        )
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
