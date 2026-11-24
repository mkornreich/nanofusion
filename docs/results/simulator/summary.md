# Every route in one simulator: the search for net energy

Produced by `tools/simulator_search.py`. Code: `nanofusion/improved/simulator.py` (`python -m nanofusion sim list | run | search`).

## The account

Each of the 22 routes is wrapped in the same energy account:

electricity → driver (laser light, accelerator beam or pulsed current) → energy in the fuel → fusion → heat (with a lithium blanket that breeds tritium: +4.78 MeV per D-T fusion) → electricity.

* Driver efficiencies: laser 10% (fusion-laser designs; NIF's own is ~0.5%), pulsed power 15%, accelerators 40% (a design figure). Light reaching crystal ions or fuel heat: 15%; implosions put 1% of their drive into the fuel. Heat to electricity: 40%.
* **Electric gain** = electricity out / electricity in. Net energy means more than 1.
* **Stretch**: each design lists what it demands (laser energy and power, field, density, pressure, current, how long a pinch stays stable, muon supply, ...) against the best ever achieved. The largest ratio is its stretch. Ideas never shown to work (spin-polarised fuel surviving in a plasma, a microbubble implosion) count as infinite stretch.
* The search samples each route 300 times (log-uniform over each range, including the ends, with every option on and off), then hill-climbs 100 steps from the best design within each budget.
* **The reference route**, `compressed_icf`, is ordinary inertial fusion: an isobaric hot spot (10 keV, ρR 0.3 g/cm²) inside a dense shell, 1% of the drive in the fuel, burn-up ρR/(ρR + 7 g/cm²). With NIF-like inputs it gives a target gain of ~2 against NIF's record 4.1, so it is if anything conservative.

## Best electric gain of each route, by how far beyond today's technology it may go

| Route (section) | today's technology | ≤ 10× today | ≤ 100× today | unlimited |
|---|---|---|---|---|
| dynamic_pinch (24, 25) | 3.1×10⁻³ | **22** | **22** | **22** |
| compressed_icf (7, 11) | – | **1.5** | **3** | **7.1** |
| fibre_pinch (19) | 0.01 | 0.89 | **2.8** | **5.7** |
| pycnonuclear (10) | – | – | – | **5.3** |
| tamped_ball (12) | 8.6×10⁻⁴ | 4.5×10⁻³ | 0.012 | **1.3** |
| muon_squeezed (17) | – | – | – | **1.1** |
| fast_ignition (7) | – | – | – | 0.41 |
| crystal_ball (3) | 7×10⁻⁴ | 1.4×10⁻³ | 3.1×10⁻³ | 0.19 |
| muon_catalysis (9) | 0.12 | 0.16 | 0.16 | 0.16 |
| magnetised_liner (12, 19) | 1.2×10⁻⁴ | 2.1×10⁻⁴ | 4.9×10⁻⁴ | 0.046 |
| son_fisch_pB (17) | – | 6.3×10⁻³ | 0.01 | 0.015 |
| muon_factory (9, 10) | – | 5.8×10⁻³ | 5.9×10⁻³ | 5.9×10⁻³ |
| channelling (20) | 5.4×10⁻³ | 5.4×10⁻³ | 5.4×10⁻³ | 5.4×10⁻³ |
| opaque_pB (18) | – | – | 1.4×10⁻³ | 4.3×10⁻³ |
| pauli_beam (17) | 4.1×10⁻⁴ | 1.4×10⁻³ | 2.1×10⁻³ | 2.3×10⁻³ |
| projectiles (16) | – | – | 9×10⁻⁴ | 9×10⁻⁴ |
| cold_catcher (1, 5) | 8.8×10⁻⁴ | 8.8×10⁻⁴ | 8.8×10⁻⁴ | 8.8×10⁻⁴ |
| pB_catcher (4, 5) | 3×10⁻⁴ | 3×10⁻⁴ | 3×10⁻⁴ | 3×10⁻⁴ |
| lattice_confinement (9) | 4.8×10⁻⁸ | 4.8×10⁻⁸ | 4.8×10⁻⁸ | 4.8×10⁻⁸ |
| hollow_shell (8) | – | – | – | – |
| antimatter_spark (22) | – | – | – | – |
| laser_assisted (22) | – | – | – | – |

Bold: more than 1, net electricity on paper.

## The winning designs

### Today's technology

**muon_catalysis** (Muon-catalysed fusion in liquid or solid D-T): phi 1.4, sticking PSI 1987 (0.45%), rf_stripping none.

* Fusion 4.4×10⁻¹⁰ J from 7.5×10⁻¹⁰ J of accelerator; electricity out / in **0.12**. 155 fusions per muon; section 9.
* liquid or solid D-T density (× liquid hydrogen): 1.4 against 1.4 (solid D-T), ×1

### ≤ 10× today

**dynamic_pinch** (D-T pinch followed in time, driven by a generator circuit): current_ma 34, load liner, radius_um 43, column_mm 3.4, line_mg_cm 0.05, liner_mm 4.5, liner_aspect 12, fill_mg_cc 0.3, bz_tesla 4.3, ice_um 580, ice_burn ablation, tritium 0.5, preheat_kev 0.017, alfven_times 1.2, stability breaks.

* Fusion 2.9×10⁹ J from 9.9×10⁶ J of pulsed power, 470 J of laser; electricity out / in **22**. Blown apart at 126 ns; 35.0 MA reached; 24.6% burnt; hottest 81 keV; narrowest 20 µm; zpinch.py.
* pinch current (A): 3.4×10⁷ against 2.6×10⁷ (Sandia's Z), ×1.3
* energy into the pinch (J): 2.6×10⁶ against 3×10⁶ (Z (2.7 MJ of X-rays from ~20 MJ stored)), ×1
* radial squeeze of the fuel: 170 against 35 (NIF-class implosions, 30–40), ×5 ← the hardest
* laser energy per shot (J): 470 against 2.1×10⁶ (NIF), ×1
* magnetic field on the target (T): 4.3 against 800 (laser-driven coil, claimed), ×1
* compressed magnetic field (T): 1.4×10⁴ against 7,100 (NIF magnetised implosion (Sio et al. 2023)), ×2
* liner aspect ratio (radius / wall): 12 against 9 (Z's record MagLIF shots (coated liners, mass-matched AR 9)), ×1.3

### ≤ 100× today

**dynamic_pinch** (D-T pinch followed in time, driven by a generator circuit): current_ma 34, load liner, radius_um 43, column_mm 3.4, line_mg_cm 0.05, liner_mm 4.5, liner_aspect 12, fill_mg_cc 0.3, bz_tesla 4.3, ice_um 580, ice_burn ablation, tritium 0.5, preheat_kev 0.017, alfven_times 1.2, stability breaks.

* Fusion 2.9×10⁹ J from 9.9×10⁶ J of pulsed power, 470 J of laser; electricity out / in **22**. Blown apart at 126 ns; 35.0 MA reached; 24.6% burnt; hottest 81 keV; narrowest 20 µm; zpinch.py.
* pinch current (A): 3.4×10⁷ against 2.6×10⁷ (Sandia's Z), ×1.3
* energy into the pinch (J): 2.6×10⁶ against 3×10⁶ (Z (2.7 MJ of X-rays from ~20 MJ stored)), ×1
* radial squeeze of the fuel: 170 against 35 (NIF-class implosions, 30–40), ×5 ← the hardest
* laser energy per shot (J): 470 against 2.1×10⁶ (NIF), ×1
* magnetic field on the target (T): 4.3 against 800 (laser-driven coil, claimed), ×1
* compressed magnetic field (T): 1.4×10⁴ against 7,100 (NIF magnetised implosion (Sio et al. 2023)), ×2
* liner aspect ratio (radius / wall): 12 against 9 (Z's record MagLIF shots (coated liners, mass-matched AR 9)), ×1.3

### Unlimited

**dynamic_pinch** (D-T pinch followed in time, driven by a generator circuit): current_ma 34, load liner, radius_um 43, column_mm 3.4, line_mg_cm 0.05, liner_mm 4.5, liner_aspect 12, fill_mg_cc 0.3, bz_tesla 4.3, ice_um 580, ice_burn ablation, tritium 0.5, preheat_kev 0.017, alfven_times 1.2, stability breaks.

* Fusion 2.9×10⁹ J from 9.9×10⁶ J of pulsed power, 470 J of laser; electricity out / in **22**. Blown apart at 126 ns; 35.0 MA reached; 24.6% burnt; hottest 81 keV; narrowest 20 µm; zpinch.py.
* pinch current (A): 3.4×10⁷ against 2.6×10⁷ (Sandia's Z), ×1.3
* energy into the pinch (J): 2.6×10⁶ against 3×10⁶ (Z (2.7 MJ of X-rays from ~20 MJ stored)), ×1
* radial squeeze of the fuel: 170 against 35 (NIF-class implosions, 30–40), ×5 ← the hardest
* laser energy per shot (J): 470 against 2.1×10⁶ (NIF), ×1
* magnetic field on the target (T): 4.3 against 800 (laser-driven coil, claimed), ×1
* compressed magnetic field (T): 1.4×10⁴ against 7,100 (NIF magnetised implosion (Sio et al. 2023)), ×2
* liner aspect ratio (radius / wall): 12 against 9 (Z's record MagLIF shots (coated liners, mass-matched AR 9)), ×1.3

### The best route that uses exploding crystals, within 100× today

**tamped_ball** (Ball held by a gold tamper and a magnetic field): radius_mm 0.35, t_i_kev 30, b_tesla 3,900, fill 1, ions first, fuel DT.

* Fusion 10⁵ J from 4.1×10⁵ J of laser; electricity out / in **0.012**. Gold tamper, 1% entrance hole; burn died out; section 12.
* laser energy per shot (J): 4.1×10⁵ against 2.1×10⁶ (NIF), ×1
* picosecond laser energy (J): 4.1×10⁵ against 4,400 (NIF ARC, four 30 ps beamlets), ×93 ← the hardest
* laser peak power (W): 1.2×10¹⁶ against 10¹⁶ (ELI-NP, 10 PW), ×1.1
* laser intensity (W/cm²): 3×10¹⁸ against 1.1×10²³ (2021 record), ×1
* magnetic field on the target (T): 3,900 against 800 (laser-driven coil, claimed), ×4.9

## If the untested theories work

The same search with the untested ideas let in: spin-polarised fuel, alpha channelling, direct conversion of fusion products, freeing stuck muons with fields, antimatter sparks, light-assisted tunnelling, the microbubble implosion, and ways to keep a dense pinch stable. Stretch still counts the technology each design needs.

| Route (section) | today's technology | ≤ 10× today | ≤ 100× today | unlimited |
|---|---|---|---|---|
| dynamic_pinch (24, 25) | 0.22 | **37** | **37** | **37** |
| compressed_icf (7, 11) | – | **1.8** | **3.5** | **8.2** |
| fibre_pinch (19) | 0.018 | **1.2** | **3.4** | **7.1** |
| pycnonuclear (10) | – | – | – | **5.3** |
| tamped_ball (12) | 1.3×10⁻³ | 6.8×10⁻³ | 0.019 | **1.9** |
| muon_squeezed (17) | – | – | – | **1.1** |
| laser_assisted (22) | – | 0.032 | 0.091 | 0.85 |
| fast_ignition (7) | – | – | – | 0.41 |
| antimatter_spark (22) | – | 0.16 | 0.39 | 0.39 |
| muon_catalysis (9) | 0.14 | 0.2 | 0.2 | 0.2 |
| crystal_ball (3) | 7.9×10⁻⁴ | 1.6×10⁻³ | 3.1×10⁻³ | 0.19 |
| magnetised_liner (12, 19) | 1.8×10⁻⁴ | 3.3×10⁻⁴ | 7.4×10⁻⁴ | 0.077 |
| son_fisch_pB (17) | – | 9.4×10⁻³ | 0.015 | 0.023 |
| muon_factory (9, 10) | – | 0.02 | 0.021 | 0.021 |
| opaque_pB (18) | – | – | 2×10⁻³ | 6.9×10⁻³ |
| channelling (20) | 5.4×10⁻³ | 5.4×10⁻³ | 5.4×10⁻³ | 5.4×10⁻³ |
| pauli_beam (17) | 7.3×10⁻⁴ | 2.1×10⁻³ | 3.2×10⁻³ | 3.5×10⁻³ |
| projectiles (16) | – | – | 9×10⁻⁴ | 9×10⁻⁴ |
| cold_catcher (1, 5) | 8.8×10⁻⁴ | 8.8×10⁻⁴ | 8.8×10⁻⁴ | 8.8×10⁻⁴ |
| hollow_shell (8) | 4.3×10⁻⁴ | 4.3×10⁻⁴ | 4.3×10⁻⁴ | 4.3×10⁻⁴ |
| pB_catcher (4, 5) | 3×10⁻⁴ | 3×10⁻⁴ | 3×10⁻⁴ | 3×10⁻⁴ |
| lattice_confinement (9) | 4.8×10⁻⁸ | 4.8×10⁻⁸ | 4.8×10⁻⁸ | 4.8×10⁻⁸ |

Ideas the best designs use:

* dynamic_pinch: burn spreading into a liner's frozen fuel (high-gain MagLIF; simulated only); sheared-flow stabilisation of a pinch at fibre density (shown only below 0.6 MA, 5e24 m⁻³)
* compressed_icf: spin-polarised D-T that survives in the fuel
* fibre_pinch: spin-polarised D-T that survives in the fuel
* tamped_ball: spin-polarised D-T that survives in the fuel
* laser_assisted: light-assisted tunnelling (Queisser & Schützhold)
* antimatter_spark: an antiproton factory 10,000 times more efficient than CERN's; antiprotons light a hot spot by fission (ICAN, AIMStar)
* muon_catalysis: stuck muons re-accelerated and freed by fields (Kulsrud, TEK)
* magnetised_liner: spin-polarised D-T that survives in the fuel
* son_fisch_pB: alpha channelling (waves hand alpha energy to the ions); direct conversion of fusion products at 60%
* muon_factory: stuck muons re-accelerated and freed by fields (Kulsrud, TEK)
* opaque_pB: alpha channelling (waves hand alpha energy to the ions); direct conversion of fusion products at 60%
* pauli_beam: spin-polarised D-T that survives in the fuel
* hollow_shell: a microbubble implosion (simulated only)

## What each idea adds

Each row holds one idea at one setting and lets the search re-optimise everything else. Untested ideas are let into the search only in their "on" rows.

| Route | Idea | Setting | today's technology | ≤ 10× today | ≤ 100× today | unlimited |
|---|---|---|---|---|---|---|
| compressed_icf | ignitor | hot spot | – | 1.5 | 3 | 7.1 |
| compressed_icf | ignitor | crystal spark | – | – | – | 0.64 |
| compressed_icf | driver | laser | – | 1.5 | 3 | 7.1 |
| compressed_icf | driver | pulsed-power X-rays | – | 0.61 | 1.4 | 3.7 |
| compressed_icf | spin (untested) | False | – | 1.5 | 3 | 7.1 |
| compressed_icf | spin (untested) | True | – | 1.8 | 3.5 | 8.2 |
| tamped_ball | ions_first | False | 3.4×10⁻⁴ | 7.1×10⁻⁴ | 6.9×10⁻³ | 1.2 |
| tamped_ball | ions_first | True | 8.6×10⁻⁴ | 3.3×10⁻³ | 0.011 | 1.3 |
| tamped_ball | fuel | DT | 8.6×10⁻⁴ | 4.5×10⁻³ | 0.012 | 1.3 |
| tamped_ball | fuel | DD | 2×10⁻⁵ | 1.8×10⁻⁴ | 3×10⁻⁴ | 0.015 |
| tamped_ball | spin (untested) | False | 8.6×10⁻⁴ | 4.5×10⁻³ | 0.012 | 1.3 |
| tamped_ball | spin (untested) | True | 1.3×10⁻³ | 6.8×10⁻³ | 0.019 | 1.8 |
| fibre_pinch | spin (untested) | False | 0.01 | 0.89 | 2.8 | 5.7 |
| fibre_pinch | spin (untested) | True | 0.017 | 1.2 | 3.3 | 6.9 |
| magnetised_liner | spin (untested) | False | 1.2×10⁻⁴ | 2.1×10⁻⁴ | 4.9×10⁻⁴ | 0.046 |
| magnetised_liner | spin (untested) | True | 1.8×10⁻⁴ | 3.2×10⁻⁴ | 7.4×10⁻⁴ | 0.077 |
| fast_ignition | shaped | False | – | – | – | 0.073 |
| fast_ignition | shaped | True | – | – | – | 0.41 |
| cold_catcher | source | crystals | 3.3×10⁻⁵ | 3.3×10⁻⁵ | 3.3×10⁻⁵ | 3.3×10⁻⁵ |
| cold_catcher | source | accelerator | 8.8×10⁻⁴ | 8.8×10⁻⁴ | 8.8×10⁻⁴ | 8.8×10⁻⁴ |
| cold_catcher | screening_ev | 0.0 | 8.8×10⁻⁴ | 8.8×10⁻⁴ | 8.8×10⁻⁴ | 8.8×10⁻⁴ |
| cold_catcher | screening_ev | 1500.0 | 8.8×10⁻⁴ | 8.8×10⁻⁴ | 8.8×10⁻⁴ | 8.8×10⁻⁴ |
| pB_catcher | avalanche | False | 3×10⁻⁴ | 3×10⁻⁴ | 3×10⁻⁴ | 3×10⁻⁴ |
| pB_catcher | avalanche | True | 3×10⁻⁴ | 3×10⁻⁴ | 3×10⁻⁴ | 3×10⁻⁴ |
| pB_catcher | screening_ev | 0.0 | 3×10⁻⁴ | 3×10⁻⁴ | 3×10⁻⁴ | 3×10⁻⁴ |
| pB_catcher | screening_ev | 1500.0 | 3×10⁻⁴ | 3×10⁻⁴ | 3×10⁻⁴ | 3×10⁻⁴ |
| muon_catalysis | sticking | PSI 1987 (0.45%) | 0.12 | 0.16 | 0.16 | 0.16 |
| muon_catalysis | sticking | theory (0.65%) | 0.09 | 0.11 | 0.11 | 0.11 |
| muon_catalysis | rf_stripping (untested) | none | 0.12 | 0.16 | 0.16 | 0.16 |
| muon_catalysis | rf_stripping (untested) | Kulsrud/TEK | 0.14 | 0.2 | 0.2 | 0.2 |
| opaque_pB | alpha_channelling (untested) | 0.0 | – | – | 1.4×10⁻³ | 4.3×10⁻³ |
| opaque_pB | alpha_channelling (untested) | 0.5 | – | – | 1.4×10⁻³ | 4.5×10⁻³ |
| opaque_pB | alpha_channelling (untested) | 0.9 | – | – | 1.4×10⁻³ | 4.6×10⁻³ |
| opaque_pB | direct_conversion (untested) | False | – | – | 1.4×10⁻³ | 4.3×10⁻³ |
| opaque_pB | direct_conversion (untested) | True | – | – | 2×10⁻³ | 6.5×10⁻³ |
| son_fisch_pB | alpha_channelling (untested) | 0.0 | – | 6.3×10⁻³ | 0.01 | 0.015 |
| son_fisch_pB | alpha_channelling (untested) | 0.5 | – | 6.6×10⁻³ | 0.011 | 0.015 |
| son_fisch_pB | alpha_channelling (untested) | 0.9 | – | 6.6×10⁻³ | 0.011 | 0.015 |
| son_fisch_pB | direct_conversion (untested) | False | – | 6.3×10⁻³ | 0.01 | 0.015 |
| son_fisch_pB | direct_conversion (untested) | True | – | 9.4×10⁻³ | 0.015 | 0.022 |
| dynamic_pinch | load | fibre | 2.3×10⁻⁷ | 2.3 | 5.2 | 5.2 |
| dynamic_pinch | load | gas column | 2.6×10⁻⁷ | 1.9 | 4.2 | 4.2 |
| dynamic_pinch | load | liner | 2.5×10⁻³ | 12 | 15 | 18 |
| dynamic_pinch | stability (untested) | breaks | 3.1×10⁻³ | 22 | 22 | 22 |
| dynamic_pinch | stability (untested) | rising current | 3.1×10⁻³ | 22 | 22 | 22 |
| dynamic_pinch | stability (untested) | sheared flow | 0.22 | 22 | 22 | 22 |
| dynamic_pinch | ice_burn (untested) | ablation | 3.1×10⁻³ | 22 | 22 | 22 |
| dynamic_pinch | ice_burn (untested) | as an ICF shell | 3.1×10⁻³ | 34 | 35 | 36 |

## What the search found

* **With today's technology nothing reaches breakeven.** The best is muon catalysis at 0.12.
* **Ordinary inertial fusion is the most proven route beyond that.** Within ten times today's best on every demand it reaches 1.5 (target gain 29); within a hundred, 3. The ten-times design needs a drive of 19 MJ, a cold shell of 5.8 g/cm² (9.6 times NIF's) at an adiabat of 1. That matches the standard view that inertial fusion power needs target gains of tens to a hundred and an efficient driver.
* **The dense Z-pinch, static**: 0.89 within ten times today, 2.8 within a hundred, if it stays stable for many Alfvén times.
* **Followed in time, the pinch leads within ten times today** (dynamic_pinch: 22, against the static 0.89; 22 within a hundred). Its best design is a 4.5 mm liner with a wall 1/12 of its radius, 580 µm of D-T ice and 4.3 T at 34 MA, converging 170× (NIF's best: 35). That is in effect cylindrical inertial fusion, in a 0-D model with no liner instability, so its convergence is the number to doubt. A cold fibre can't take the current; a gas column or a liner can. With untested ideas it reaches 37 within ten times today.
* **Exploding crystals never help.** The best crystal route within a hundred times today is tamped ball at 0.012. As the spark for inertial fusion, crystals do worse than an ordinary hot spot at every budget (unlimited: 0.64 against 7.1): the isotropic explosion wastes the light. As a muon factory they make muons, but not cheaply.
* **Lasers beat pulsed-power X-rays** as the driver of a compressed capsule (best within ten times today: 1.5 against 0.61). Pulsed power makes its drive more efficiently, but on Z only ~0.2% of the stored energy reached a capsule.
* **Small ideas stay small.** Electron screening and the proton–boron avalanche change nothing you can see; ions heated first help the tamped ball two to five times; deuterium-only fuel costs a factor of tens; spin-polarised fuel helps but has never been shown to survive in a plasma.
* **Untested theory moves the numbers, not the verdict.** With it, the best within ten times today is dynamic pinch at 37 (with burn spreading into a liner's frozen fuel (high-gain MagLIF; simulated only)); within a hundred, dynamic pinch at 37 (with burn spreading into a liner's frozen fuel (high-gain MagLIF; simulated only)).
* **Proton–boron with alpha channelling and direct conversion** reaches 0.023 with no limits on technology, against 0.015 without them. Channelling feeds the ions, but the X-rays still win.
* **Antimatter sparks** reach at most 0.39, with 1.3×10¹⁴ antiprotons per shot, 22 times the most ever stored, from a factory 10,000 times more efficient than CERN's. Making antiprotons costs more than the fusion they spark.
* **Light-assisted tunnelling** reaches at most 0.85, with no limits on technology: holding 5×10²¹ W/cm² over the fuel costs more light than the extra fusion returns, even with the most optimistic tunnelling boost.
* **Freeing stuck muons with fields** takes muon catalysis from 0.12 to 0.14 with today's technology: twice the fusions per muon, but keeping the fuel below 25 K costs ~30 MeV of electricity per 17.6 MeV fusion.
* **Corrections the search forced.** It found a pinch "burning" more than all its fuel, a pinch at 130,000 Gbar, a pinch too cool for its alphas to make up its X-ray losses, and a compressed capsule whose hot spot was 600 times the pressure of the shell around it. Burn-up, a pressure demand, X-ray resupply and an isobaric hot spot are now in the models; those designs are gone.

## With a uranium blanket (k_eff = 0.95)

A subcritical uranium blanket turns each D-T neutron into fission heat: about 270 times the fusion energy. That is a fission reactor driven by a neutron source (section 13).

| Route | today's technology | ≤ 10× today |
|---|---|---|
| dynamic_pinch | 4.3×10⁻⁵ | **940** |
| compressed_icf | – | **310** |
| fibre_pinch | **2.1** | **75** |
| muon_catalysis | **25** | **34** |
| muon_factory | – | **1.2** |
| channelling | **1.2** | **1.2** |
| pauli_beam | 0.075 | 0.28 |
| tamped_ball | 0.16 | 0.26 |
| cold_catcher | 0.18 | 0.18 |
| crystal_ball | 0.043 | 0.14 |
| magnetised_liner | 9.4×10⁻³ | 0.024 |
| son_fisch_pB | – | 5.7×10⁻³ |
| pB_catcher | 3×10⁻⁴ | 3×10⁻⁴ |
| lattice_confinement | 4.8×10⁻⁸ | 4.8×10⁻⁸ |
| fast_ignition | – | – |
| muon_squeezed | – | – |
| pycnonuclear | – | – |
| projectiles | – | – |
| opaque_pB | – | – |
| hollow_shell | – | – |
| antimatter_spark | – | – |
| laser_assisted | – | – |

* **Fission makes several routes "net energy"**, muon catalysis with today's accelerators among them, but over 99% of the heat is fission. That is the conclusion of section 13 again.

