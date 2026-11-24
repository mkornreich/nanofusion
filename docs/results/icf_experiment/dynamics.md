# What the static model leaves out, simulated

Produced by `tools/icf_dynamics.py` with the 1-D implosion-and-burn model `nanofusion/improved/burn1d.py`. The proposal is [docs/icf_experiment.md](../../icf_experiment.md); the static model's numbers are in [summary.md](summary.md).

## The model

* The D-T shell and the ablator left on it coast inward from peak velocity (300 µm, 4.5 g/cm³).
* They stagnate on the D-T gas (1.8 µg), heat it, and ablate fuel into it.
* The hot spot ignites, or not, and the burn spreads and disassembles the capsule.
* Included: alphas deposited by range, flux-limited conduction, bremsstrahlung with free-free re-absorption, an ideal-plus-degenerate equation of state, burn-up.
* Checked: the same yield within 1% at 2.6× the zones and at half the time step; energy conserved to 0.5%.

## Calibration: the clean implosion

| Shot | Inputs (LLNL post-shot) | 1-D yield here | LLNL's 1-D | Static ceiling | Measured | Measured ÷ 1-D |
|---|---|---|---|---|---|---|
| N210808 | 401 km/s, adiabat 2.65, 149 µg of ablator | 5.06 MJ | 4.2 MJ (HYDRA), 6.2 MJ (xRAGE) | 10.1 MJ | 1.37 MJ | 27% |
| N221204 | 394 km/s, adiabat 2.45, 242 µg of ablator | 20 MJ | ~15 MJ (idealised 1-D at ~2.1 MJ); 13.6 MJ in 2-D without low modes | 10.9 MJ | 3.15 MJ | 16% |

* The record, N250406 (8.6 MJ), reached 43% of the clean N221204-class yield. LLNL puts N240210 at about 35% of its 1-D yield; their goal is 80%.
* **The static ceiling is too low.** The clean 1-D yield of N221204's design is 20 MJ, about 1.8× the static model's 10.9 MJ. So the static model's "quality" overstated how close NIF has come.

## The physics the static model leaves out

Each piece isolated on N221204's clean implosion:

| Physics | Simulated | Effect on the yield |
|---|---|---|
| Alpha heating | off: 38 kJ at 4.0 keV | ×532 (N210808: ×149) |
| Burn propagation and disassembly | burn-up 27% | the static formula ρR/(ρR + 7) gives 21% at the same areal density |
| Ablator mixed into the hot spot | 100 ng of carbon | ×0.98 |
| Ablator mixed into the hot spot | 300 ng of carbon | ×0.94 |
| Ablator mixed into the hot spot | 1000 ng of carbon | ×0.26 |
| Asymmetry: kinetic energy that never stagnates | 1% | ×0.97 |
| Asymmetry: kinetic energy that never stagnates | 3% | ×0.88 |
| Asymmetry: kinetic energy that never stagnates | 10% | ×0.02 |
| The cost of pressure: velocity at fixed kinetic energy | 360, 380, 394, 410, 430 km/s | 0.144 MJ, 14 MJ, 20 MJ, 21.4 MJ, 21.7 MJ; pressure without alphas ∝ v^2.4 |
| The adiabat at fixed velocity | 2.45 → 2.25 | ×1.11; pressure +6%, areal density +5% (no alphas) |

* **Alpha heating is most of the yield.** Without it the capsule makes tens of kJ. With it, the hot spot runs away to 20–30 keV and the burn spreads through the shell.
* **A clean capsule shrugs off small degradations; one near the cliff doesn't.** 100 ng of mixed-in carbon costs a clean N221204 2%, but the measured NIF trend is 30–50% per 100 ng, on degraded shots closer to the cliff.
* **Velocity has a cliff.** Below ~380 km/s the capsule fails to ignite.
* **The adiabat matters through pressure and areal density**, which the static model held fixed.

## Retrodicting Hohenberger et al. (2025): the adiabat test the static model fails

NIF lowered the adiabat alone, 2.89 → 2.68, on N210808's target at the same laser energy. Areal density (DSR) rose 14% and the yield up to 80%.

| Model | Yield ratio, adiabat 2.68 ÷ 2.89 | Areal density ratio |
|---|---|---|
| Static (fixed pressure and fuel mass) | ×1.01 | ×1.00 |
| 1-D, clean | ×6.1 | ×1.04 (no alphas) |
| 1-D, with N210808's degradation | ×2.8 | ×1.03 |
| Measured | up to ×1.8 | ×1.14 |

* The 1-D model gets the direction right, where the static model can't. It overshoots the size: N210808's target sits on the ignition cliff in 1-D, so small changes in areal density move the yield a lot. Hohenberger et al. find the same steepness (yield ∝ ρR^11).

## Each shot's degradation

Every degradation (mixing, low-mode asymmetry, defects) is folded into one number: the carbon mixed into the hot spot that brings the clean 1-D yield down to the measured one.

| Shot | Clean 1-D | Measured | Mix-equivalent carbon | Burn-up | Burn-averaged T |
|---|---|---|---|---|---|
| N210808 | 5.06 MJ | 1.37 MJ | 220 ng | 1.8% | 8.8 keV |
| N221204 | 20 MJ | 3.15 MJ | 1035 ng | 4.2% | 10.1 keV |
| N250406 | 20 MJ | 8.6 MJ | 946 ng | 11.5% | 15.2 keV |

* N210808's equivalent, 220 ng, is within NIF's measured bound (less than 250 ng of ablator in its hot spot and its repeats). N221204's, 1035 ng, is not: there the number stands in for low-mode asymmetry, which LLNL finds cost N221204 a factor of 4.5.
* **The degraded shots sit on the cliff.** The record's equivalent (946 ng) is only 9% below N221204's, yet its yield is 2.7× higher. A single mixing knob is far too steep to hold fixed, so the campaign below is bracketed two ways.

## The campaign, re-predicted

* **Fixed share:** each step's clean 1-D yield × 43%, the record's share of its own clean yield (the usual way to carry a degradation forward).
* **Fixed mixing:** the record's 946 ng held fixed. Each step that helps pushes the capsule back up the cliff, so this is an upper bound.

| Step | What changes | Velocity | Fuel | Clean 1-D | Fixed share | Fixed mixing |
|---|---|---|---|---|---|---|
| A. Repeat the record | today's target and pulse | 394 km/s | 0.22 mg | 20 MJ | **8.6 MJ** | 8.59 MJ |
| B. Lower adiabat | adiabat 2.45 → 2.25 | 394 km/s | 0.22 mg | 22.2 MJ | **9.54 MJ** | 18.5 MJ |
| C. Faster | +4% velocity from a thinner ablator, same kinetic energy | 410 km/s | 0.22 mg | 21.4 MJ | **9.22 MJ** | 19.4 MJ |
| D0. 2.2 MJ, thicker ablator | what NIF did: 5% more ablator left | 403 km/s | 0.22 mg | 22.3 MJ | **9.6 MJ** | 20 MJ |
| D. 2.2 MJ, thicker fuel | Landen et al.'s proposal: 7% more D-T | 401 km/s | 0.236 mg | 24.2 MJ | **10.4 MJ** | 21.9 MJ |
| E. D + lower adiabat |  | 401 km/s | 0.236 mg | 25.7 MJ | **11.1 MJ** | 24.6 MJ |
| F. D + lower adiabat + faster |  | 417 km/s | 0.236 mg | 25.8 MJ | **11.1 MJ** | 25.5 MJ |
| G. 2.6 MJ, scaled capsule | every length × 1.08 | 394 km/s | 0.279 mg | 29.6 MJ | **12.8 MJ** | 27.5 MJ |
| H. G + lower adiabat + faster |  | 410 km/s | 0.279 mg | 31.3 MJ | **13.5 MJ** | 31 MJ |
| I. G + adiabat 2.0 + 8% faster | a stretch | 426 km/s | 0.279 mg | 32.4 MJ | **14 MJ** | 32.5 MJ |

* **B, lower adiabat:** ×1.11 (fixed mixing ×2.2). Not a null result: the static model's central claim fails, as Hohenberger et al. found on NIF.
* **C, faster at the same kinetic energy:** ×1.07 (fixed mixing ×2.3), if the thinner ablator doesn't mix more.
* **D vs D0, where to spend 2.2 MJ:** more fuel 10.4 MJ, more ablator (NIF's choice) 9.6 MJ. Landen et al. (2024) proposed the switch; the 1-D model agrees it helps, modestly.
* **The 2.6 MJ upgrade:** 12.8 MJ–27.5 MJ for a scaled capsule, up to 14 MJ–32.5 MJ stretched. At LLNL's goal of 80% of the 1-D yield: 23.7 MJ–25.9 MJ. LLNL expects "greater than 30 MJ"; Lindl et al. very likely > 20 MJ, 50–70 MJ optimised.

