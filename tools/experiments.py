"""More experiments with the improved model, looking for ways to get fusion.

    python3 tools/experiments.py [--only TEXT] [--out docs/results/experiments]

Groups:
  shock   uniform vs pre-expanded crystals at several sizes: do shock
          shells make ions collide inside the crystal, and how does that grow with size?
  pair    two crystals exploding into each other, at several gaps
  wave    Li-7 "guided" collisions with a long-wavelength (10 µm CO2) laser
Each case writes a JSON report; summary.md collects the key numbers.
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

BASE = ImprovedProperties(intensity_w_cm2=1e20, phase_deg=0.0, duration_fs=80.0, frames=10, lithium="Li6")


def cases() -> list[tuple[str, ImprovedProperties]]:
    out = []
    for li in ("Li6", "Li7"):
        for r in (1.0, 1.5, 2.0, 3.0):
            out.append((f"shock-{li}-uniform-R{r:g}", replace(BASE, lithium=li, size=r)))
            out.append((f"shock-{li}-expanded3-R{r:g}", replace(BASE, lithium=li, size=r, profile="expanded")))
        for e in (2.0, 5.0):
            out.append((f"shock-{li}-expanded{e:g}-R2", replace(BASE, lithium=li, size=2.0, profile="expanded",
                                                               expansion=e)))
    for gap in (0.5, 2.0, 5.0, 10.0):
        out.append((f"pair-gap{gap:g}", replace(BASE, size=2.0, geometry="pair", gap_nm=gap)))
    for r in (1.0, 1.5):  # size scaling of colliding crystals, gap = 2.5 R (like pair-gap5 at 2 nm)
        out.append((f"pair-R{r:g}-gap{2.5 * r:g}", replace(BASE, size=r, geometry="pair", gap_nm=2.5 * r,
                                                          duration_fs=80.0)))
    for i in (1e20, 1e21):
        for ph in (0.0, 90.0):
            out.append((f"wave-Li7-10um-I{i:.0e}-ph{ph:.0f}",
                        replace(BASE, lithium="Li7", size=2.0, wavelength_um=10.6, intensity_w_cm2=i,
                                phase_deg=ph, duration_fs=60.0)))
    return out


def main() -> int:
    ap = argparse.ArgumentParser()
    ap.add_argument("--only")
    ap.add_argument("--out", default="docs/results/experiments")
    args = ap.parse_args()
    out_dir = Path(args.out)
    out_dir.mkdir(parents=True, exist_ok=True)
    rows = []
    for name, props in cases():
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
    lines = ["| case | ions | max D (keV) | highest D–Li pair (keV) | closest D–Li (Å) | in-crystal D+D "
             "| in-crystal D+Li | beam-target LiD | beam-target LiT |", "|---" * 9 + "|"]
    for name, r in rows:
        d, f = r["dynamics"], r["fusion"]
        lines.append(f"| {name} | {r['particles']} | {d['max_D_keV']:.3g} | {d['max_D_Li_pair_energy_keV']:.3g} "
                     f"| {d['closest_D_Li_approach_A']:.3g} | {f['in_crystal_D_D']:.2g} | {f['in_crystal_D_Li']:.2g} "
                     f"| {f['catchers']['LiD']['total']:.2g} | {f['catchers']['LiT']['total']:.2g} |")
    (out_dir / "summary.md").write_text("\n".join(lines) + "\n")
    print("\n".join(lines))
    return 0


if __name__ == "__main__":
    sys.exit(main())
