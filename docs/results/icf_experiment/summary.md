# A practical experiment with ordinary inertial fusion: testing what limits NIF's yield

Produced by `tools/icf_experiment.py` from the simulator's `compressed_icf` model (`nanofusion/improved/simulator.py`). The proposal is [docs/icf_experiment.md](../../icf_experiment.md).

## Calibration on NIF's measured shots

The model has no mixing and no asymmetry, so it gives the yield of a *clean* implosion: a ceiling. It is matched to NIF's fuel (0.22 mg) at each shot's measured hot-spot pressure and design adiabat:

| Shot | Laser | Hot-spot pressure | Adiabat | Areal density in the model | Energy in the fuel | Ceiling | Measured | Quality (measured ÷ ceiling) |
|---|---|---|---|---|---|---|---|---|
| N210808 | 1.92 MJ | 569 Gbar | 2.65 | 1.10 g/cm² | 41 kJ (2.1% of the laser) | 10.1 MJ | 1.37 MJ | 14% |
| N221204 | 2.05 MJ | 608 Gbar | 2.45 | 1.21 g/cm² | 38 kJ (1.8% of the laser) | 10.9 MJ | 3.15 MJ | 29% |

* The model's areal density, 1.2 g/cm², is close to NIF's simulated 1.35–1.41 g/cm² with alpha heating turned off (Kritcher et al. 2024).
* The coupling, 1.8% of the laser energy into the fuel, is fitted on N221204 and held fixed. It absorbs the model's simplifications.
* The model's hot spot is the smallest that ignites (38 µm, 10 keV). NIF's burning hot spots are larger and hotter (~68 µm, 13 keV for N221204): alpha heating, which the model leaves out.

Later shots used the same fuel; their pressures are not published, so they are compared with N221204's ceiling:

| Shot | Laser | Ceiling | Measured | Quality |
|---|---|---|---|---|
| N230729 | 2.05 MJ | 10.9 MJ | 3.88 MJ | 36% |
| N240210 | 2.20 MJ | 10.9 MJ | 5.20 MJ | 48% |
| N250406 | 2.08 MJ | 10.9 MJ | 8.60 MJ | 79% |

* **Quality is the story of NIF since 2021:** 14% → 29% → 36% → 48% → 79%. Fewer capsule defects, then a graded tungsten dopant that cut the mixing of ablator into the fuel. The record (N250406) came close to the clean ceiling.
* Part of the later gains may be higher pressure, which this comparison can't separate. Measuring pressure on every shot of the campaign separates them.

## The campaign

One change per step. Yield predicted at the record's quality (79%) and at the best quality before the graded dopant (N240210, 48%). Coupling held at 1.8%.

| Step | What changes in practice | Laser | Hot-spot pressure | Adiabat | D-T fuel (ice) | Ceiling | Yield at 79% | at 48% | Target gain |
|---|---|---|---|---|---|---|---|---|---|
| **Phase 1: today's target** | | | | | | | | | |
| A. Repeat the record | N250406's target and pulse | 2.08 MJ | 608 Gbar | 2.45 | 0.22 mg (66 µm) | 10.9 MJ | **8.6 MJ** | 5.2 MJ | 4.1 |
| B. Lower adiabat | A with a lower foot in the pulse | 2.08 MJ | 608 Gbar | 2.25 | 0.22 mg (66 µm) | 11.0 MJ | **8.6 MJ** | 5.2 MJ | 4.2 |
| C. Higher pressure | A with a faster, rounder implosion | 2.08 MJ | 700 Gbar | 2.45 | 0.22 mg (66 µm) | 12.7 MJ | **10.0 MJ** | 6.1 MJ | 4.8 |
| **Phase 2: a thicker fuel layer at full energy** | | | | | | | | | |
| D. Thicker fuel | full laser energy, more ice | 2.20 MJ | 608 Gbar | 2.45 | 0.25 mg (77 µm) | 13.6 MJ | **10.7 MJ** | 6.5 MJ | 4.9 |
| E. Thicker fuel, lower adiabat | D with B's pulse | 2.20 MJ | 608 Gbar | 2.25 | 0.27 mg (81 µm) | 14.7 MJ | **11.6 MJ** | 7.0 MJ | 5.3 |
| F. Thicker fuel, both | D with B's pulse and C's pressure | 2.20 MJ | 700 Gbar | 2.25 | 0.31 mg (95 µm) | 21.5 MJ | **16.9 MJ** | 10.2 MJ | 7.7 |
| **Phase 3: NIF's 2.6 MJ upgrade** | | | | | | | | | |
| G. Upgraded laser | D at 2.6 MJ | 2.60 MJ | 608 Gbar | 2.45 | 0.34 mg (107 µm) | 21.6 MJ | **17.0 MJ** | 10.3 MJ | 6.5 |
| H. Upgraded laser, both | F at 2.6 MJ | 2.60 MJ | 700 Gbar | 2.25 | 0.39 mg (126 µm) | 31.3 MJ | **24.6 MJ** | 14.9 MJ | 9.5 |
| I. Upgraded laser, stretched | H with 800 Gbar and adiabat 2.0 | 2.60 MJ | 800 Gbar | 2.00 | 0.44 mg (144 µm) | 42.0 MJ | **33.1 MJ** | 20.0 MJ | 12.7 |

What the model predicts, step by step:
* **B: the static model says lowering the adiabat alone does nothing** (×1.01), because it holds the pressure fixed. **NIF has already falsified this:** Hohenberger et al. (2025) lowered the adiabat from 2.89 to 2.68 on the same target and got up to 80% more yield. At fixed drive, a lower adiabat raises the pressure and areal density; the 1-D model (dynamics.md) shows it.
* **C: pressure is the lever** (×1.17 for +15%). The model even needs less energy in the fuel (1.85 MJ of laser against 2.05), because a smaller hot spot ignites at higher pressure. It doesn't price what the faster implosion costs.
* **D: the last 6% of laser energy into thicker ice** gives ×1.25.
* **E: once the fuel is re-sized, the adiabat matters** (×1.08 over D): cheaper compression leaves energy for more fuel.
* **F, H, I: pressure again dominates.** These rows assume the higher pressure costs no extra laser energy, so they are an upside.
* **The upgrade.** At 2.6 MJ: 17 MJ at today's pressure, 25 MJ with 700 Gbar, 33 MJ stretched. Published projections: LLNL expects "greater than 30 MJ"; Lindl et al. (2026) very likely more than 20 MJ, possibly 50–70 MJ fully optimised; Hurricane's hydro-scaling 24 MJ at 3 MJ of laser.

## Energy

With a lithium blanket and 40% of the heat turned into electricity:

| Step | Target gain | Electricity out ÷ in with NIF's laser (~0.6%) | with a 10% plant laser |
|---|---|---|---|
| A. Repeat the record | 4.1 | 0.012 | 0.21 |
| B. Lower adiabat | 4.2 | 0.012 | 0.21 |
| C. Higher pressure | 4.8 | 0.014 | 0.25 |
| D. Thicker fuel | 4.9 | 0.015 | 0.25 |
| E. Thicker fuel, lower adiabat | 5.3 | 0.016 | 0.27 |
| F. Thicker fuel, both | 7.7 | 0.023 | 0.39 |
| G. Upgraded laser | 6.5 | 0.020 | 0.33 |
| H. Upgraded laser, both | 9.5 | 0.028 | 0.48 |
| I. Upgraded laser, stretched | 12.7 | 0.038 | 0.65 |

No step makes net electricity. The best, I, would return 65% of its electricity with a plant-grade laser. The simulator's net-electricity designs (section 21 of docs/energy.md) need about ten times NIF's laser energy and areal density.

