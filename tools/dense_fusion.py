"""Fusion in dense matter across all regimes: is there a cheaper path than hot fuel?

Uses the all-regime rate of Gasques et al. (Phys. Rev. C 72, 025806 (2005)):
classical thermonuclear, strongly screened, thermo-pycnonuclear, thermally
enhanced pycnonuclear and zero-temperature pycnonuclear burning, with their
optimal, maximum- and minimum-rate model variants. For every density and
temperature it finds the smallest ball of D-T that burns 1% of itself before
flying apart, and what it costs to compress and heat it. The question is
whether strong screening or zero-point motion in dense, cold fuel (a crystal
or a quantum liquid) ever beats ordinary hot fuel. It also tests Ichimaru's
prediction of very fast reactions in liquid metallic hydrogen.

Writes docs/results/dense/summary.md.

    python3 tools/dense_fusion.py
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1]))

from nanofusion.improved import alternatives as A  # noqa: E402
from nanofusion.improved import energy as E  # noqa: E402

OUT = Path("docs/results/dense")
ICHIMARU_T_EV = 950 * 8.617333e-5  # 950 K
OPTIMAL = "optimal (static bcc lattice)"


def validation() -> list[str]:
    lines = ["## Check: the hot-plasma limit", "",
             "At low density and high temperature the formula must give the ordinary Maxwellian rate "
             "(n⟨σv⟩/4 per nucleus for D-T, n⟨σv⟩/2 for D-D).", "",
             "| T | D-T: formula / Maxwellian | D-D: formula / Maxwellian |", "|---|---|---|"]
    for t in (2e3, 5e3, 1e4, 3e4):
        n_dt = 1e3 / (2.515 * 1.66053906660e-27)
        n_dd = 1e3 / (2.0141 * 1.66053906660e-27)
        dt = A.dense_fusion_rate(1.0, t, "DT") / (n_dt * E.reactivity("D", "T", t / 1e3) / 4)
        dd = A.dense_fusion_rate(1.0, t, "DD") / (n_dd * E.reactivity("D", "D", t / 1e3) / 2)
        lines.append(f"| {t / 1e3:g} keV | {dt:.2f} | {dd:.2f} |")
    lines += ["", f"At zero temperature it gives {A.dense_fusion_rate(1e7, 0.0):.2g} /s per nucleus at 10⁷ g/cm³, "
              f"against {A.pycnonuclear_rate_per_nucleus(1e7, 'DT'):.2g} /s from Salpeter & Van Horn directly.", ""]
    return lines


def rate_map() -> list[str]:
    temps = [0.0, 0.1, 1.0, 10.0, 100.0, 1e3, 1e4]
    lines = ["## D-T reactions per nucleus per second (optimal model), log₁₀", "",
             "| density (g/cm³) | " + " | ".join("0" if t == 0 else (f"{t:g} eV" if t < 1e3 else f"{t / 1e3:g} keV")
                                                for t in temps) + " |",
             "|---|" + "---|" * len(temps)]
    for rho in (0.21, 20.0, 1e3, 1e5, 1e7):
        cells = []
        for t in temps:
            r = A.dense_fusion_rate(rho, t)
            cells.append(f"{np.log10(r):.0f}" if r > 0 else "−∞")
        lines.append(f"| {rho:g} | " + " | ".join(cells) + " |")
    lines += ["", "A rate of ~10⁹–10¹¹ per second is needed to burn before a small ball flies apart. Hot fuel gets "
              "there at ~10 keV once compressed to ~10³ g/cm³. Cold fuel gets there only above ~10⁷ g/cm³.", ""]
    return lines


def ichimaru_corner() -> tuple[list[str], dict]:
    # density at which cold D-T reaches 1 Gbar
    rhos = np.geomspace(1, 1e3, 400)
    pressures = [A.cold_compression(r)["pressure_Gbar"] for r in rhos]
    rho_1gbar = float(np.interp(1.0, pressures, rhos))
    lines = ["## Ichimaru's corner: liquid metallic hydrogen near freezing", "",
             "Ichimaru & Kitamura (Phys. Plasmas 6, 2649 (1999)) predicted very high p-d rates at ≥ 20 g/cm³, ≥ 1 Gbar, "
             f"~950 K. Cold D-T reaches 1 Gbar at {rho_1gbar:.0f} g/cm³. Reactions per nucleus per second, at 950 K:",
             "", "| model | D-T, 20 g/cm³ | D-T, " + f"{rho_1gbar:.0f} g/cm³ | D-D, 20 g/cm³ | time to burn 1% of D-T at "
             f"{rho_1gbar:.0f} g/cm³ |", "|---|---|---|---|---|"]
    data = {"rho_1gbar": rho_1gbar}
    for model in A.DENSE_MODELS:
        r20 = A.dense_fusion_rate(20.0, ICHIMARU_T_EV, "DT", model)
        r1g = A.dense_fusion_rate(rho_1gbar, ICHIMARU_T_EV, "DT", model)
        rdd = A.dense_fusion_rate(20.0, ICHIMARU_T_EV, "DD", model)
        years = 0.01 / r1g / 3.156e7 if r1g > 0 else np.inf
        data[model] = {"DT_20": r20, "DT_1Gbar": r1g, "DD_20": rdd}
        lines.append(f"| {model} | {r20:.1e} | {r1g:.1e} | {rdd:.1e} | {years:.0e} years |")
    lines += ["", "Every model, including the one built on Ichimaru's own Monte Carlo screening and the 'maximum rate' "
              "envelope, puts these rates near 10⁻⁶⁰ or below. The age of the Universe is ~10¹⁰ years. Ichimaru's "
              "p-d rates may be 'prodigiously high' by astrophysical standards, but not for energy. The reason is "
              "that 950 K sits far below the ion plasma temperature here (~1 eV). Ions then sit in their quantum "
              "ground states, and only zero-point motion helps them tunnel.", ""]
    return lines, data


def cheapest_ignition() -> tuple[list[str], dict]:
    temps = np.concatenate([[0.0], np.geomspace(0.01, 1e5, 141)])
    data = {}
    lines = ["## Hot against cold fuel at the same density", "",
             "The smallest D-T ball that ignites: it burns 1% of itself before flying apart "
             "(rate × R / 4c_s ≥ 0.01) and is thick enough to stop its own alpha particles (ρR ≥ 0.3 g/cm²). "
             "Energy = what the ball holds after compressing and heating it from cold solid, without the "
             "driver's losses. At 100 g/cm³ the cheapest hot ball holds ~20 kJ, comparable with NIF's hot spot "
             "(~10 kJ), whose laser needed ~2 MJ. Tiny balls at extreme densities are only a scaling exercise.", "",
             "| density | model | cheapest hot ball (≥ 1 keV) | cheapest warm ball (10 eV – 1 keV) "
             "| cheapest cold, crystal-like ball (< 10 eV) |", "|---|---|---|---|---|"]

    def fmt(e, t):
        if not np.isfinite(e) or e > 1e40:
            return "never (bigger than a planet)"
        t_s = "0 K" if t == 0 else (f"{t:.2g} eV" if t < 1e3 else f"{t / 1e3:.2g} keV")
        return f"{e:.1g} J at {t_s}"

    for rho in (1e2, 1e3, 1e5, 1e7):
        for model in A.DENSE_MODELS:
            es = np.array([A.ignition_ball(rho, t, "DT", model)["energy_J"] for t in temps])
            row = {}
            for key, mask in (("hot", temps >= 1e3), ("warm", (temps >= 10) & (temps < 1e3)), ("cold", temps < 10)):
                sub = np.where(mask, es, np.inf)
                j = int(np.argmin(sub))
                row[key] = (float(sub[j]), float(temps[j]))
            data[f"{rho:g} {model}"] = row
            lines.append(f"| {rho:g} g/cm³ | {model} | {fmt(*row['hot'])} | {fmt(*row['warm'])} | {fmt(*row['cold'])} |")
    lines += ["", "* **Hot fuel wins at every density and in every model.** The cheapest balls burn at 16–25 keV, "
              "the ordinary thermonuclear regime, where screening hardly matters (Γ ≪ 1).",
              "* **Warm fuel (10 eV – 1 keV) always costs far more**, and its cheapest point is always the hottest edge "
              "of the range. There is no sweet spot where screening pays: strong screening only grows when the ions "
              "are cold (Γ ≫ 1), and by then bare tunnelling has collapsed by far more orders than screening "
              "restores.",
              "* **Cold, crystal-like fuel** never ignites below ~10⁷ g/cm³ in any model (the maximum-rate model "
              "allows a ~7 km ball at 10⁵ g/cm³, holding 10³² J). Even at 10⁷ g/cm³ it needs a ~µm ball, 7 J to 10¹¹ J depending "
              "on the model, against ~10⁻⁵ J for a hot ball at the same density.",
              ""]
    return lines, data


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    lines = ["# Fusion in dense matter: is there a cheaper path than hot fuel?", "",
             "Produced by `tools/dense_fusion.py`. Rate model: Gasques et al., Phys. Rev. C 72, 025806 (2005), "
             "Eq. (28) with their Table I; code in `nanofusion/improved/alternatives.py`.", ""]
    lines += validation() + rate_map()
    ich_lines, ich = ichimaru_corner()
    ign_lines, ign = cheapest_ignition()
    lines += ich_lines + ign_lines
    text = "\n".join(lines) + "\n"
    (OUT / "summary.md").write_text(text)
    (OUT / "results.json").write_text(json.dumps({"ichimaru": ich, "ignition": ign}, indent=1, default=str))
    print(text)
    return 0


if __name__ == "__main__":
    sys.exit(main())
