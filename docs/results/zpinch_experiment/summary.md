# A practical dense Z-pinch experiment

Produced by `tools/zpinch_experiment.py` with the time-dependent pinch model `nanofusion/improved/zpinch.py`. The proposal is [docs/zpinch_experiment.md](../../zpinch_experiment.md).

## Checks

* **Pease-Braginskii current** (X-rays beat ohmic heating): 1.32–1.58 MA for a Coulomb logarithm of 7–10, against Haines' 0.433 √lnΛ = 1.15–1.37 MA.
* **Bennett temperature of FuZE (2018):** 200 kA, ~1e23 per m³, radius ~3 mm give 2.2 keV; measured 1–2 keV.
* **The 1980s–90s deuterium fibres**, cold start. Thermonuclear D-D neutrons from the model, against the measured neutrons, which were mostly beam-target from instabilities:

| Fibre experiment | Breaks in ~1 Alfvén time | Stable while the current rises | Sheared flow | Measured (mostly beam-target) |
|---|---|---|---|---|
| NRL Poseidon (1987) | < 1 | < 1 | 3.7e+05 | 8e9 (1e12 claimed) |
| LANL HDZP-II | < 1 | 2.9e+04 | 2.5e+12 | 5e9 |
| MAGPIE | < 1 | 3.6e+03 | 5e+10 | 1e8–5e9 |

The model's thermonuclear yields are far below the measured ones in every scenario, which agrees with the finding that those neutrons came from beams made by the instabilities.

## Three barriers the static model left out

A 240 µm D-T fibre, 1 cm long, on Z (26 MA into a fixed load in 100 ns; 15 nH of machine inductance):

| Case | Current reached | Narrowest | Hottest | Ended | Fusion | X-rays | From the bank |
|---|---|---|---|---|---|---|---|
| Cold fibre, stable while the current rises | 4.2 MA | 43.1 µm | 0.31 keV | broke up at 12 ns | 0.038 J | 8.48 kJ | 252 kJ of 7.76 MJ |
| Preheated to 100 eV, breaks in ~1 Alfvén time (1980s fibres) | 5.5 MA | 240 µm | 0.16 keV | broke up at 13 ns | 0 J | 1.02 kJ | 368 kJ of 7.76 MJ |
| Preheated to 100 eV, stable while the current rises (NRL 1987) | 5.5 MA | 83 µm | 0.62 keV | broke up at 14 ns | 0.23 J | 3.4 kJ | 390 kJ of 7.76 MJ |
| Preheated to 100 eV, sheared flow, fully stable (untested at this density) | 8.5 MA | 5.46 µm | 1.33 keV | end of pulse at 150 ns | 1.33 MJ | 6.4 MJ | 6.89 MJ of 7.76 MJ |

1. **A cold fibre can't take the current.** At ~10 eV its resistance is far above the generator's impedance, so the current stalls at 4.2 MA. Real fibres grew a hot, thin corona that carried the current instead, and broke up. Preheating the fibre (100 eV costs ~5 kJ per cm, within Z-Beamlet's reach, as in MagLIF) makes the dense column conduct.
2. **The collapse caps its own current.** Above the Pease-Braginskii current the column radiates faster than it heats and collapses. The collapse's back-voltage holds Z to 8.5 MA, and most of the energy leaves as X-rays.
3. **Stability decides everything, and even full stability isn't enough on Z.** Fully flow-stabilised, the column makes 1.33 MJ (burn-up 0.9%, 1.3 keV at most): no ignition.
* The static route (26 MA, 10 keV, one Alfvén time, squeeze 6) assumed a Bennett pinch already hot: 504 kJ of fusion for an electric gain of 0.0089. The dynamic model shows the current never gets there.

## The experiment: preheated deuterium fibres on Z, with a 1% tritium trace

Each row is one fibre radius; each cell is what the model predicts in one stability scenario: thermonuclear D-D neutrons / D-T neutrons from the tritium trace / narrowest radius.

| Fibre radius | breaks in ~1 Alfvén time (1980s fibres) | stable while the current rises (NRL 1987) | sheared flow, fully stable (untested at this density) |
|---|---|---|---|
| 60 µm | < 1 / < 1 / 60 µm | 1.3e+09 / 1.8e+09 / 27 µm | 1.2e+15 / 2.3e+15 / 0.8 µm |
| 120 µm | < 1 / < 1 / 1.2e+02 µm | 2.2e+09 / 2.5e+09 / 48 µm | 3e+15 / 3.7e+15 / 3.4 µm |
| 240 µm | < 1 / < 1 / 2.4e+02 µm | 2e+09 / 2e+09 / 90 µm | 1.4e+16 / 1.7e+16 / 5.5 µm |

* Without preheat (120 µm, stable during the rise): 2.6e+08 D-D neutrons, at 2.3 MA.
* The scenarios differ by orders of magnitude in thermonuclear neutrons, and the D-T/D-D ratio measures the ion temperature. Beam-target neutrons, the 1980s background, are anisotropic and spectrally shifted; thermonuclear ones are not.

## Outlook: D-T, fully flow-stabilised

| Machine | Fibre radius | Current reached | Fusion (1 cm) | Burn-up | Hottest | Bank | Electricity out ÷ in |
|---|---|---|---|---|---|---|---|
| 26 MA | 120 µm | 4.1 MA | 268 kJ | 0.7% | 2.2 keV | 8.23 MJ | 0.0025 |
| 26 MA | 240 µm | 8.5 MA | 1.33 MJ | 0.9% | 1.3 keV | 7.76 MJ | 0.013 |
| 26 MA | 480 µm | 15.3 MA | 814 kJ | 0.1% | 1.0 keV | 7.29 MJ | 0.0085 |
| 26 MA | 960 µm | 20.1 MA | 6.97 kJ | 0.0% | 0.8 keV | 6.82 MJ | 7.8e-05 |
| 26 MA | 1920 µm | 28.4 MA | 0.94 J | 0.0% | 0.3 keV | 6.35 MJ | 1.1e-08 |
| 60 MA | 120 µm | 5.6 MA | 13.5 MJ | 35.1% | 1.8 keV | 43.8 MJ | 0.023 |
| 60 MA | 240 µm | 29.0 MA | 29.6 MJ | 19.3% | 33.4 keV | 41.3 MJ | 0.055 |
| 60 MA | 480 µm | 43.0 MA | 185 MJ | 30.1% | 84.9 keV | 38.8 MJ | 0.36 |
| 60 MA | 960 µm | 69.8 MA | 924 MJ | 37.6% | 118.1 keV | 36.3 MJ | 1.9 |
| 60 MA | 1920 µm | 46.7 MA | 312 kJ | 0.0% | 0.9 keV | 33.8 MJ | 0.0007 |

* Each shot pays for the whole bank: once the column has burnt and blown apart, the energy left in it is wasted. With a lithium blanket, 40% of heat to electricity and 15% of electricity reaching the load.
* **On Z, even fully flow-stabilised, the best is 0.013 of breakeven.**
* **On a 60 MA machine, wide fibres ignite** and burn up to a third of their fuel. The best, 960 µm at 60 MA, returns 1.9× the electricity: net energy, but only if the pinch stays stable for its whole burn, which no dense pinch has done.

