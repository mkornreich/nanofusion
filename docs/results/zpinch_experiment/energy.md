# Making the dense Z-pinch testable, and what it takes to make energy

Produced by `tools/zpinch_energy.py` with `nanofusion/improved/zpinch.py`. The proposal is [docs/zpinch_experiment.md](../../zpinch_experiment.md).

## Levers that don't help the fibre

A 960 µm D-T fibre at 60 MA (fully flow-stabilised: 1.9× the electricity in):

| Lever | Fusion | Current reached | Ended | Electricity out ÷ in |
|---|---|---|---|---|
| rising current, 100 ns rise | 19 J | 22.6 MA | broke up at 27 ns | 4.0e-08 |
| rising current, 300 ns rise | 0.05 J | 21.1 MA | broke up at 68 ns | 1.0e-10 |
| rising current, 1000 ns rise | 0 J | 10.0 MA | blown apart at 82 ns | 0 |
| 10 mm long, sheared flow | 858 MJ | 69.0 MA | blown apart at 73 ns | 1.8 |
| 30 mm long, sheared flow | 36.2 MJ | 32.9 MA | end of pulse at 150 ns | 0.05 |
| 100 mm long, sheared flow | 13.2 MJ | 28.8 MA | end of pulse at 150 ns | 8.4e-03 |

* **Keeping the current rising** (NRL's stabilisation) fails whatever the driver: the column radiates, cools and turns resistive, so the current peaks at 10–23 MA and falls.
* **A longer fibre** has more resistance and inductance, so it draws less current.
* **Recovering the circuit's energy** after the burn is a different machine: an inductive store that takes the field's energy back. Krakowski's 0-D fibre-pinch study finds recovering half of it halves the fusion gain needed (14 → 7.5). Nobody has shown it at this scale; Helion claims >95% without published numbers.
* **The limit is the load:** a dense, cold column can't take the current.

## A different load: a pre-formed gas column

A column of gas a centimetre or two wide (a gas puff, as Z already fires), ionised by the current. It conducts easily and is heated by its own implosion, so it needs no preheat laser.

**Calibration on Z's deuterium gas puffs** (record 3.9×10¹³ D-D neutrons at 15–17 MA; ion temperature 9.5 keV; stagnation radius 2.5 mm). Model at 17 MA, 2 cm long:

| Gas column | Held 1 Alfvén time | 2 | 3 |
|---|---|---|---|
| 0.25 mg/cm, 10 mm | 8.5e+07 n, 0.6 keV, 0.65 mm | 2.7e+14 n, 6.4 keV, 0.11 mm | 3.4e+14 n, 6.4 keV, 0.11 mm |
| 0.25 mg/cm, 20 mm | 1.3e+08 n, 0.7 keV, 1.1 mm | 8.5e+14 n, 12.0 keV, 0.12 mm | 1.1e+15 n, 12.0 keV, 0.12 mm |
| 0.5 mg/cm, 10 mm | 8.3e+08 n, 0.6 keV, 0.56 mm | 3.2e+14 n, 4.1 keV, 0.13 mm | 3.4e+14 n, 4.1 keV, 0.13 mm |
| 0.5 mg/cm, 20 mm | 1.7e+09 n, 0.8 keV, 0.89 mm | 1.2e+15 n, 7.3 keV, 0.16 mm | 1.2e+15 n, 7.3 keV, 0.16 mm |
| 1.0 mg/cm, 10 mm | 6.6e+09 n, 0.7 keV, 0.48 mm | 2.3e+14 n, 2.5 keV, 0.16 mm | 2.3e+14 n, 2.5 keV, 0.16 mm |
| 1.0 mg/cm, 20 mm | 2.4e+10 n, 0.9 keV, 0.73 mm | 8e+14 n, 4.1 keV, 0.22 mm | 8.1e+14 n, 4.1 keV, 0.22 mm |

* **The record sits between one and two Alfvén times** in the model: one gives at most 2.4e+10 neutrons, two at least 2.3e+14. The pinch heats in its second Alfvén time, so the yield climbs steeply there.
* **The model converges too far.** Held two Alfvén times, its columns stagnate at 0.1–0.2 mm, against the 2.5 mm measured: a perfectly uniform column squeezes harder than a real, unstable implosion. The simulator charges that squeeze as a demand on technology.

## The experiment Z can do now: deuterium gas columns with a 1% tritium trace

26 MA, 1 cm wide, 2 cm long. Each cell: D-D neutrons / D-T neutrons from the trace / narrowest radius.

| Fuel load | 1 Alfvén time | 3 | 10 | 30 | sheared flow |
|---|---|---|---|---|---|
| 0.5 mg/cm | 3.8e+08 / 3.8e+08 / 0.63 mm | 1.8e+15 / 6.1e+15 / 0.096 mm | 8e+15 / 3e+16 / 0.076 mm | 5.8e+16 / 2.2e+17 / 0.062 mm | 1.4e+17 / 5.2e+17 / 0.061 mm |
| 2.0 mg/cm | 2.5e+10 / 2.6e+10 / 0.44 mm | 1.3e+15 / 2.7e+15 / 0.13 mm | 3.1e+15 / 6.3e+15 / 0.08 mm | 3.4e+15 / 6.9e+15 / 0.032 mm | 3.5e+15 / 6.9e+15 / 0.031 mm |
| 7.0 mg/cm | 8.8e+10 / 7.8e+10 / 0.38 mm | 2.6e+13 / 2.9e+13 / 0.16 mm | 2.8e+13 / 3.1e+13 / 0.15 mm | 2.8e+13 / 3.1e+13 / 0.15 mm | 2.8e+13 / 3.1e+13 / 0.15 mm |

* Gas puffs are Z's routine load, deuterium needs no new licence, and the 1% trace is the next step after Z's 2016 tritium shot.
* **What it measures** is the number the energy map below turns on: how many Alfvén times a dense pinch holds. Neutron yield, the D-T/D-D ratio (ion temperature), X-ray imaging of the stagnated column and its duration give it directly.

## What it takes to make energy

Best electricity out ÷ in over fuel load (2–20 mg/cm) and radius (3–20 mm) for the gas column, and radius (0.24–1.5 mm) for a preheated fibre. 50/50 D-T, a lithium blanket, 40% of heat to electricity, 15% of electricity reaching the load, and each shot paying for the whole bank.

| Load | Current | 1 Alfvén time | 3 | 10 | 30 | sheared flow |
|---|---|---|---|---|---|---|
| gas column | 26 MA | 3.5e-07 | 0.044 | 0.084 | 0.12 | 0.12 |
| gas column | 45 MA | 1.9e-07 | 0.48 | 0.48 | 0.48 | 0.86 |
| gas column | 60 MA | 7.2e-08 | 0.85 | 0.85 | 0.85 | 2 |
| gas column | 80 MA | 2.4e-08 | 2.6 | 3.3 | 3.3 | 3.8 |
| fibre | 26 MA | 5.5e-11 | 4.2e-05 | 5.1e-05 | 6.3e-05 | 0.013 |
| fibre | 45 MA | 0 | 8.5e-05 | 1.2e-04 | 2.0e-04 | 0.12 |
| fibre | 60 MA | 0 | 5.1e-05 | 1.6e-04 | 3.2e-04 | 1.8 |
| fibre | 80 MA | 0 | 1.2e-07 | 2.0e-04 | 4.2e-04 | 3.2 |

* **The gas column needs far less stability than the fibre.** With 3–10 Alfvén times it reaches 0.85 at 60 MA and 3.3 at 80 MA; the fibre needs about fifty (section 24).
* **But it needs deep convergence.** The columns that make energy burn for 2.5–7.5 Alfvén times after squeezing 162–263 times, to 40–130 µm. Z's gas puffs stagnate at 2.5 mm, and the best implosion of any kind (NIF) converges about 35 times.
* **The best tested-physics design**, a gas column at 80 MA held 10 Alfvén times: 3.3× the electricity in, 19.9 mg/cm in a 10 mm column, 2.37 GJ per cm, 35% burnt, from a 54 MJ bank.
* **On Z (26 MA) nothing reaches breakeven**, whatever the load or stability.

**Tritium share** (gas column, 10 Alfvén times, fuel load and radius optimised):

| Current | 50% | 20% | 10% | 5% | deuterium only |
|---|---|---|---|---|---|
| 60 MA | 0.85 | 0.41 | 0.13 | 0.092 | 1.5e-03 |
| 80 MA | 3.3 | 0.43 | 0.34 | 0.084 | 2.8e-03 |

* Lean tritium costs gain steeply: D-T burns ~100 times faster than D-D at these temperatures, and the tritium D-D breeds in place is too little to make up for it.

**A more efficient driver.** Gain scales with the share of electricity that reaches the load: 15% here (about Z), 17–24% stored-to-load in the Z 800 and Z 300 designs (Stygar et al. 2015), 60% assumed in LLNL's 2006 Z-pinch plant study.

| Design | 15% | 24% | 60% |
|---|---|---|---|
| gas column, 60 MA, 3 Alfvén times | 0.85 | 1.4 | 3.4 |
| gas column, 80 MA, 10 Alfvén times | 3.3 | 5.4 | 13 |

**Against MagLIF.** Pacific Fusion's published D-T ice MagLIF designs give 109–380 MJ at 56 MA and 816 MJ at 65 MA. The gas column gives 341 MJ per cm at 60 MA here: the same class, without a liner, but only if the bare column holds its convergence and a few Alfvén times.

## With a liner (MagLIF)

A beryllium tube carries the current and implodes onto laser-preheated fuel gas with an axial field, and optionally a layer of frozen D-T on its inner wall. The liner is incompressible and thick; the field leaks out as the fuel converges (it grows as convergence^1.57, matching MagLIF's measured B R); heat and alphas that reach the ice ablate it into the hot fuel. No current flows in the fuel, so there is no m = 0 mode; the liner's own instability is charged as a demand on convergence.

**Calibration on Z's MagLIF shots** (deuterium):

| Shot | D-D neutrons: model / clean 2-D / measured | Ion temperature: model / measured | Burn width: model / measured | Convergence: model / measured |
|---|---|---|---|---|
| 2014 first shots, 19 MA | 5.5e+14 / 6.1e+13 / 2e+12 | 10.0 / 2.5 keV | 0.04 / 2.0 ns | 81 / 46 |
| z3289, 20 MA | 6.1e+14 / 5.2e+13 / 1.1e+13 | 7.1 / 3.1 keV | 0.06 / 2.0 ns | 66 / 31 |
| record z3926, 18.5 MA | 6.2e+14 / 5.9e+13 / 2.3e+13 | 7.4 / 3.25 keV | 0.12 / 1.75 ns | 42 / 40 |

* **The 0-D liner overshoots at Z's scale**, 9–12× the clean 2-D yields: too hot, and its burn far too brief. An incompressible liner bounces off the fuel at once; a real one decelerates layer by layer and holds the fuel for ~2 ns.

**The published high-current designs**, rerun here:

| Design (its own simulated yield) | Here, ice by ablation | Here, ice burnt as an ICF shell |
|---|---|---|
| Slutz et al. 2010, D-T gas, 27 MA (1-D: 0.5 MJ) | 47.1 kJ | – |
| the same at 45 MA (~10 MJ) | 1.04 MJ | – |
| the same at 60 MA (>50 MJ) | 4.82 MJ | – |
| Sefkow et al. 2014, D-T ice, 70 MA (1-D: 7.9 GJ) | 1.4 MJ | – |
| Pacific Fusion 2025, D-T ice, 56 MA (109 MJ, Be) | 3.33 MJ | 97.4 MJ |

* **At high current the model undershoots**: the gas designs by ~10×, and Sefkow's 70 MA ice design doesn't ignite at all. The same too-brief stagnation is the likely cause.
* **Pacific Fusion's 56 MA design** is matched only if the ice burns like an ICF shell. So this model is for comparing designs and seeing which demand binds, not for absolute yields; the published 1-D and 2-D designs are the better guide.

**Liner designs for energy** (50/50 D-T, 0.2 keV preheat, best over liner radius 3–6 mm, gas 1–3 mg/cm³ and field 10–30 T). Electricity out ÷ in, ice burnt by ablation / as an ICF shell:

| Current | Gas only | 300 µm of ice | 600 µm of ice | Gas column, 10 Alfvén times |
|---|---|---|---|---|
| 26 MA | 6.0e-03 | 7.1e-03 / 7.1e-03 | 7.2e-03 / 7.2e-03 | 0.084 |
| 45 MA | 0.014 | 0.045 / 0.9 | 0.1 / 2.8 | 0.48 |
| 60 MA | 0.029 | 1.3 / 1.5 | 3.5 / 6 | 0.85 |
| 80 MA | 0.039 | 1.2 / 1.7 | 4.2 / 5.9 | 3.3 |

* **Thick ice burns at high current, by the model's own ablation.** Below ~50 MA it barely ablates (best 0.1 at 45 MA): the field that holds the heat in keeps it from the ice. At 60–80 MA the burning fuel ablates the ice fast enough for the burn to spread.
* **The best liner design, untested ideas left out:** 4.2× the electricity in, at 80 MA: a 6 mm liner, 600 µm of ice, 1 mg/cm³ of gas and 10 T; 2.89 GJ per cm, convergence 43, compressed field 3.7 kT.
* **Its demands are modest except the current:** convergence 43 against MagLIF's measured 30–40, and a compressed field below the 7.1 kT shown on NIF. The bare gas column needed 160–260× convergence for less.
* At 60 MA: 3.5× the electricity in, at 60 MA: a 3 mm liner, 600 µm of ice, 3 mg/cm³ of gas and 30 T; 1.48 GJ per cm, convergence 131, compressed field 62.9 kT.
* If the ice burns like an ICF shell (untested), the best is 6 at 60 MA.
* 1 of the 144 liner designs hit the solver's work cap and count only what they made before it.

