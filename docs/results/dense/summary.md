# Fusion in dense matter: is there a cheaper path than hot fuel?

Produced by `tools/dense_fusion.py`. Rate model: Gasques et al., Phys. Rev. C 72, 025806 (2005), Eq. (28) with their Table I; code in `nanofusion/improved/alternatives.py`.

## Check: the hot-plasma limit

At low density and high temperature the formula must give the ordinary Maxwellian rate (n⟨σv⟩/4 per nucleus for D-T, n⟨σv⟩/2 for D-D).

| T | D-T: formula / Maxwellian | D-D: formula / Maxwellian |
|---|---|---|
| 2 keV | 0.94 | 0.97 |
| 5 keV | 0.88 | 0.95 |
| 10 keV | 0.95 | 0.92 |
| 30 keV | 1.40 | 0.86 |

At zero temperature it gives 7.7e+10 /s per nucleus at 10⁷ g/cm³, against 7.6e+10 /s from Salpeter & Van Horn directly.

## D-T reactions per nucleus per second (optimal model), log₁₀

| density (g/cm³) | 0 | 0.1 eV | 1 eV | 10 eV | 100 eV | 1 keV | 10 keV |
|---|---|---|---|---|---|---|---|
| 0.21 | -177 | -139 | -70 | -28 | -7 | 2 | 6 |
| 20 | -73 | -73 | -56 | -25 | -5 | 4 | 8 |
| 1000 | -28 | -28 | -28 | -18 | -3 | 6 | 10 |
| 100000 | -2 | -2 | -2 | -2 | 1 | 8 | 12 |
| 1e+07 | 11 | 11 | 11 | 11 | 11 | 11 | 14 |

A rate of ~10⁹–10¹¹ per second is needed to burn before a small ball flies apart. Hot fuel gets there at ~10 keV once compressed to ~10³ g/cm³. Cold fuel gets there only above ~10⁷ g/cm³.

## Ichimaru's corner: liquid metallic hydrogen near freezing

Ichimaru & Kitamura (Phys. Plasmas 6, 2649 (1999)) predicted very high p-d rates at ≥ 20 g/cm³, ≥ 1 Gbar, ~950 K. Cold D-T reaches 1 Gbar at 40 g/cm³. Reactions per nucleus per second, at 950 K:

| model | D-T, 20 g/cm³ | D-T, 40 g/cm³ | D-D, 20 g/cm³ | time to burn 1% of D-T at 40 g/cm³ |
|---|---|---|---|---|
| optimal (static bcc lattice) | 2.2e-73 | 2.0e-63 | 3.3e-64 | 2e+53 years |
| maximum rate | 7.9e-66 | 1.4e-56 | 1.9e-57 | 2e+46 years |
| minimum rate | 1.1e-74 | 1.1e-64 | 1.8e-65 | 3e+54 years |
| Ogata-Iyetomi-Ichimaru (Monte Carlo) | 2.0e-68 | 3.0e-59 | 4.1e-60 | 1e+49 years |

Every model, including the one built on Ichimaru's own Monte Carlo screening and the 'maximum rate' envelope, puts these rates near 10⁻⁶⁰ or below. The age of the Universe is ~10¹⁰ years. Ichimaru's p-d rates may be 'prodigiously high' by astrophysical standards, but not for energy. The reason is that 950 K sits far below the ion plasma temperature here (~1 eV). Ions then sit in their quantum ground states, and only zero-point motion helps them tunnel.

## Hot against cold fuel at the same density

The smallest D-T ball that ignites: it burns 1% of itself before flying apart (rate × R / 4c_s ≥ 0.01) and is thick enough to stop its own alpha particles (ρR ≥ 0.3 g/cm²). Energy = what the ball holds after compressing and heating it from cold solid, without the driver's losses. At 100 g/cm³ the cheapest hot ball holds ~20 kJ, comparable with NIF's hot spot (~10 kJ), whose laser needed ~2 MJ. Tiny balls at extreme densities are only a scaling exercise.

| density | model | cheapest hot ball (≥ 1 keV) | cheapest warm ball (10 eV – 1 keV) | cheapest cold, crystal-like ball (< 10 eV) |
|---|---|---|---|---|
| 100 g/cm³ | optimal (static bcc lattice) | 2e+04 J at 16 keV | 1e+16 J at 8.9e+02 eV | never (bigger than a planet) |
| 100 g/cm³ | maximum rate | 2e+04 J at 16 keV | 1e+16 J at 8.9e+02 eV | never (bigger than a planet) |
| 100 g/cm³ | minimum rate | 2e+04 J at 16 keV | 1e+16 J at 8.9e+02 eV | never (bigger than a planet) |
| 100 g/cm³ | Ogata-Iyetomi-Ichimaru (Monte Carlo) | 2e+04 J at 16 keV | 1e+16 J at 8.9e+02 eV | never (bigger than a planet) |
| 1000 g/cm³ | optimal (static bcc lattice) | 2e+02 J at 16 keV | 1e+14 J at 8.9e+02 eV | never (bigger than a planet) |
| 1000 g/cm³ | maximum rate | 2e+02 J at 16 keV | 1e+14 J at 8.9e+02 eV | never (bigger than a planet) |
| 1000 g/cm³ | minimum rate | 2e+02 J at 16 keV | 1e+14 J at 8.9e+02 eV | never (bigger than a planet) |
| 1000 g/cm³ | Ogata-Iyetomi-Ichimaru (Monte Carlo) | 2e+02 J at 16 keV | 1e+14 J at 8.9e+02 eV | never (bigger than a planet) |
| 100000 g/cm³ | optimal (static bcc lattice) | 0.02 J at 16 keV | 2e+11 J at 8.9e+02 eV | never (bigger than a planet) |
| 100000 g/cm³ | maximum rate | 0.02 J at 16 keV | 2e+11 J at 8.9e+02 eV | 1e+32 J at 0 K |
| 100000 g/cm³ | minimum rate | 0.02 J at 16 keV | 2e+11 J at 8.9e+02 eV | never (bigger than a planet) |
| 100000 g/cm³ | Ogata-Iyetomi-Ichimaru (Monte Carlo) | 0.02 J at 16 keV | 2e+11 J at 8.9e+02 eV | never (bigger than a planet) |
| 1e+07 g/cm³ | optimal (static bcc lattice) | 1e-05 J at 25 keV | 2e+05 J at 8.9e+02 eV | 2e+06 J at 0 K |
| 1e+07 g/cm³ | maximum rate | 1e-05 J at 25 keV | 6 J at 8.9e+02 eV | 7 J at 0 K |
| 1e+07 g/cm³ | minimum rate | 1e-05 J at 25 keV | 4e+06 J at 8.9e+02 eV | 1e+09 J at 0 K |
| 1e+07 g/cm³ | Ogata-Iyetomi-Ichimaru (Monte Carlo) | 1e-05 J at 25 keV | 6e+04 J at 8.9e+02 eV | 1e+11 J at 0 K |

* **Hot fuel wins at every density and in every model.** The cheapest balls burn at 16–25 keV, the ordinary thermonuclear regime, where screening hardly matters (Γ ≪ 1).
* **Warm fuel (10 eV – 1 keV) always costs far more**, and its cheapest point is always the hottest edge of the range. There is no sweet spot where screening pays: strong screening only grows when the ions are cold (Γ ≫ 1), and by then bare tunnelling has collapsed by far more orders than screening restores.
* **Cold, crystal-like fuel** never ignites below ~10⁷ g/cm³ in any model (the maximum-rate model allows a ~7 km ball at 10⁵ g/cm³, holding 10³² J). Even at 10⁷ g/cm³ it needs a ~µm ball, 7 J to 10¹¹ J depending on the model, against ~10⁻⁵ J for a hot ball at the same density.

