# An experiment for the thin, ice-lined liner

Produced by `tools/liner_experiment.py` with `nanofusion/improved/zpinch.py`. The proposal is [docs/liner_experiment.md](../../liner_experiment.md).

## The design, as the simulator found it

| | |
|---|---|
| Liner | beryllium, outer radius 4.50 mm, wall 375 µm (1/12 of the radius), 1 cm long |
| Fuel | 583 µm of 50/50 D-T ice on the inner wall; 0.3 mg/cm³ of D-T gas inside |
| Field, preheat | 4.3 T axial; 17 eV (almost none) |
| Driver | 34 MA in 100 ns into a fixed load; bank 9.92 MJ |
| Fusion | 2.91 GJ (24% of the fuel burnt), hottest 83 keV |
| Convergence | 177×: the fuel from 3.54 mm to 20 µm |
| Electricity out ÷ in | 22.3 (15% of electricity to the load, 40% heat to electricity) |

## Its weak point: convergence

If the liner breaks up (magneto-Rayleigh-Taylor) before the fuel converges this far, compression stops there:

| Liner breaks up at | Fusion | Hottest | Electricity out ÷ in |
|---|---|---|---|
| 30× | 2.5 J | 1.1 keV | 1.9e-08 |
| 40× | 25 J | 1.6 keV | 1.9e-07 |
| 60× | 4.7e+02 J | 2.5 keV | 3.6e-06 |
| 100× | 12.2 kJ | 4.6 keV | 9.4e-05 |
| 130× | 62.2 kJ | 6.2 keV | 4.8e-04 |
| 150× | 166 kJ | 7.3 keV | 1.3e-03 |
| 170× | 533 kJ | 8.6 keV | 4.1e-03 |
| never (177×) | 2.91 GJ | 83 keV | 22.3 |

* **All or nothing.** The fuel only gets hot past ~100× and ignites near the full ~170×. At MagLIF's demonstrated 30–40×, the design makes tens of joules.

## Other sensitivities, one at a time

| Change | Fusion | Convergence | Electricity out ÷ in |
|---|---|---|---|
| 20 MA | 5.4e+02 J | 553× | 1.2e-05 |
| 26 MA | 18.4 kJ | 863× | 2.4e-04 |
| 30 MA | 172 kJ | 376× | 1.7e-03 |
| 34.06 MA (design) | 2.91 GJ | 177× | 22.3 |
| 40 MA | 2.68 GJ | 236× | 15 |
| 50 MA | 2.58 GJ | 286× | 9.19 |
| 6 wall ratio | 24.7 kJ | 895× | 1.9e-04 |
| 9 wall ratio | 358 kJ | 244× | 2.8e-03 |
| 12 wall ratio (design) | 2.91 GJ | 177× | 22.3 |
| 0 µm of ice | 1.5 MJ | 452× | 0.0115 |
| 200 µm of ice | 825 MJ | 237× | 6.35 |
| 400 µm of ice | 1.94 GJ | 197× | 14.9 |
| 583.3 µm of ice (design) | 2.91 GJ | 177× | 22.3 |
| 0.5 T | 9.12 kJ | 392× | 7.0e-05 |
| 4.283 T (design) | 2.91 GJ | 177× | 22.3 |
| 10 T | 2.38 GJ | 167× | 18.3 |
| 20 T | 31.6 MJ | 156× | 0.243 |
| 30 T | 7.36 MJ | 144× | 0.0566 |
| 0.1 mg/cm³ of gas | 2.58 GJ | 265× | 19.9 |
| 0.3 mg/cm³ of gas (design) | 2.91 GJ | 177× | 22.3 |
| 1 mg/cm³ of gas | 26.7 kJ | 247× | 2.1e-04 |
| 3 mg/cm³ of gas | 8.08 kJ | 202× | 6.2e-05 |

* **A knife-edge.** Small changes kill it: 30 MA cuts the fusion 17,000-fold, a wall ratio of 9 cuts the fusion 8,100-fold, 1 mg/cm³ of gas cuts the fusion 110,000-fold, 20 T cuts the fusion 92-fold. An optimiser's narrow window like this is more often a model artefact than robust physics.
* **Below ~34 MA the fuel over-converges cold** (hundreds of times) without igniting: the model has nothing to stop it, a real liner would break up first.

## Stage 1 on Z now: the same target in deuterium

Deuterium ice and gas, the same liner, at the current Z can drive. D-D fuel doesn't ignite, so this tests the implosion, not the burn. Each cell: D-D neutrons / stagnation radius / convergence.

| Current, wall ratio | Liner intact | Breaks up at 40× | Breaks up at 100× |
|---|---|---|---|
| 20 MA, 1/6 | 1.0e+10 / 16 µm / 201× | 1.2e+09 / 79 µm / 40× | 1.2e+09 / 32 µm / 100× |
| 20 MA, 1/9 | 2.5e+11 / 11 µm / 303× | 2.9e+10 / 85 µm / 40× | 3.7e+10 / 34 µm / 100× |
| 20 MA, 1/12 | 1.5e+12 / 9 µm / 389× | 8.2e+10 / 88 µm / 40× | 3.1e+11 / 35 µm / 100× |
| 26 MA, 1/6 | 1.6e+13 / 5 µm / 603× | 6.8e+10 / 79 µm / 40× | 4.0e+11 / 32 µm / 100× |
| 26 MA, 1/9 | 3.7e+13 / 5 µm / 677× | 1.5e+11 / 85 µm / 40× | 5.5e+12 / 34 µm / 100× |
| 26 MA, 1/12 | 5.3e+13 / 6 µm / 623× | 2.1e+11 / 88 µm / 40× | 3.2e+13 / 35 µm / 100× |

* **The stagnation radius is the measurement.** With the liner intact the model squeezes the fuel to 5–16 µm (201–677×, more than any implosion has reached: an upper bound). Breaking up at MagLIF's ~40×, it stops at 79–88 µm. X-ray imaging of the stagnated column resolves the difference.
* **The neutrons confirm it.** At 26 MA with the design's thin wall: 5.3e+13 D-D neutrons if it holds, 2.1e+11 if it breaks at 40×. Secondary D-T neutrons (from D-D tritons) give the fuel's areal density.

