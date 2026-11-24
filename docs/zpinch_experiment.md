# A practical dense Z-pinch experiment

A first experiment on the dense D-T Z-pinch, the route that came second in
the simulator's search (sections 19 and 21 of [energy.md](energy.md)). It is
designed with a new time-dependent pinch model, and checked against the
literature for whether it is novel and possible.

```bash
python3 tools/zpinch_experiment.py
```

That writes the numbers below to
[results/zpinch_experiment/summary.md](results/zpinch_experiment/summary.md).

## Is it novel? Is it possible?

**Nobody has fielded D-T in a dense equilibrium Z-pinch.** What has been done:
* **Frozen deuterium fibres, 1985–1999,** at 0.35–1.4 MA: NRL's Poseidon,
  Los Alamos' HDZP-I/II, Imperial College's MAGPIE, Kalif and S-300 (Klír's
  thesis has the tables).
  * They all failed the same way. A hot, thin corona carried the current, and
    m = 0 necks broke the pinch.
  * Their neutrons (at most ~10¹⁰) came mostly from beams made by the
    instabilities, not from thermonuclear burn.
  * The line of work was dropped around 2000.
* **Sheared-flow-stabilised pinches** (FuZE, Zap Energy):
  * stable for about 10³ Alfvén times;
  * but only at ≤ 0.6 MA and ≤ 5×10²⁴ per m³, and only with deuterium;
  * their "Q" is a D-T-equivalent figure, and no D-T is planned.
* **Tritium on pulsed power.** Sandia's Z fired contained trace tritium
  (~0.1%) once, in 2016. Going to 50/50 D-T was "at least three years" away
  then. Papers from 2024 say Z is "beginning to use D-T fuel", but no D-T
  yield has been published. Dense plasma focus facilities run D-T (Los
  Alamos: 2 MA, 1.8×10¹² neutrons; Nevada's "Sodium").

So a D-T dense pinch is **novel, but not possible today**: tritium is the
blocker. What is open and possible:
* deuterium fibres revisited on a modern driver, with modern diagnostics;
* a 1–3% tritium trace, the next step after Z's 2016 shot;
* flow stabilisation at mega-ampere currents and fibre densities.

## The model

`nanofusion/improved/zpinch.py` follows a fibre or gas column in time:
* **The driver** is a capacitor bank behind the machine's 15 nH. It is sized
  to give the stated peak current in 100 ns into a fixed load. The pinch's own
  inductance and resistance change what it actually draws.
* **The column** is squeezed by its field and held out by its pressure, with
  degenerate electrons beyond solid density.
* **Heating:** ohmic heating, compression, and the alphas its self-field
  holds in.
* **Losses:** X-rays with the column's own opacity (capped at a black body),
  and heat conducted along the field to the electrodes.
* **Instabilities:** the m = 0 mode grows at v_A/a. Three scenarios:
  * it breaks after about one Alfvén time, as the 1980s fibres did;
  * it is stable while the current rises, as NRL reported in 1987 (2-D MHD
    later blamed a cold core under an unstable corona);
  * a sheared axial flow stabilises it fully (untested at this density).

Checks:
* **Pease–Braginskii current** (X-rays beat ohmic heating): 1.32–1.58 MA,
  against Haines' 1.15–1.37 MA.
* **Bennett temperature of FuZE** (2018): 2.2 keV, against 1–2 keV measured.
* **The 1980s fibres:** the model's thermonuclear D-D yields are far below the
  measured neutrons in every scenario. That agrees with the finding that
  those neutrons were beam-target.

## What the static model missed

The simulator's `fibre_pinch` route assumed a Bennett pinch already hot,
held for some number of Alfvén times. Followed in time, a 240 µm D-T fibre on
Z (26 MA) hits three barriers:
1. **A cold fibre can't take the current.** Its resistance far exceeds the
   generator's impedance, so the current stalls at ~4 MA. That is why real
   fibres grew a corona. Preheating to ~100 eV (5 kJ of heat per cm, as
   MagLIF already does with a laser) makes the dense column conduct.
2. **The collapse caps its own current.** Above the Pease–Braginskii current
   the column radiates faster than it heats and collapses. Its back-voltage
   holds Z to ~8.5 MA, and most of the energy leaves as X-rays.
3. **Stability decides everything**:
   * broken in an Alfvén time: no fusion;
   * stable while the current rises: well under 1 J, because the collapse
     stops the current rising after ~14 ns;
   * fully flow-stabilised: 1.3 MJ from a 7.8 MJ bank, without ignition.

## The experiment

**Where:** Z, at full current, with the preheat laser MagLIF already uses.

**Targets:**
* frozen deuterium fibres with 1% tritium;
* radius 60, 120 or 240 µm, 1 cm long;
* preheated to ~100 eV. That costs ~0.65, 2.6 and 10 kJ of laser at 50%
  coupling.

**Shots:**
* each radius with and without preheat, two shots each: 12 shots;
* without preheat it repeats the 1980s result on a modern driver;
* with preheat it tests the barriers above.

**Measure:**
* the current and load voltage: the collapse's predicted back-voltage;
* X-ray power and imaging: the radiative collapse, and the narrowest radius;
* neutron yield, isotropy and spectra on several lines of sight:
  thermonuclear neutrons are isotropic and centred, beam-target ones are not;
* the D-T/D-D neutron ratio from the tritium trace: the ion temperature.

**Predictions** (thermonuclear D-D neutrons / D-T neutrons from the trace /
narrowest radius):

| Fibre radius | Breaks in ~1 Alfvén time | Stable while the current rises | Sheared flow, fully stable |
|---|---|---|---|
| 60 µm | < 1 / < 1 / 60 µm | 1.3×10⁹ / 1.8×10⁹ / 27 µm | 1.2×10¹⁵ / 2.3×10¹⁵ / 0.8 µm |
| 120 µm | < 1 / < 1 / 120 µm | 2.2×10⁹ / 2.5×10⁹ / 48 µm | 3×10¹⁵ / 3.7×10¹⁵ / 3.4 µm |
| 240 µm | < 1 / < 1 / 240 µm | 2×10⁹ / 2×10⁹ / 90 µm | 1.4×10¹⁶ / 1.7×10¹⁶ / 5.5 µm |

**How to read it:**
* **Thermonuclear neutrons near 10⁹, with a D-T/D-D ratio that matches the
  spectra:** the column held together while the current rose. That is new: no
  fibre has shown it.
* **Anisotropic, shifted neutrons with no thermonuclear signature:** the 1980s
  behaviour; the preheat didn't help.
* **Anything near 10¹⁴ or more:** stability far beyond any dense pinch so far.
  Z's best D-D gas puff gave 3.9×10¹³.
* **The collapse's back-voltage in the current trace**, at the predicted
  time, tests the radiative collapse itself, whatever the neutrons do.

## What it means for energy

The simulator now has the time-dependent pinch as a route, `dynamic_pinch`.
Untested ways to keep it stable are kept apart, like other untested theory.

| D-T, fully flow-stabilised | Best fibre | Fusion (1 cm) | Electricity out ÷ in |
|---|---|---|---|
| Z (26 MA) | 240 µm | 1.3 MJ | 0.013 |
| A 60 MA machine (studied, not built) | 960 µm | 0.92 GJ, 38% burnt | 1.9 |

* Each shot pays for the whole bank, with a lithium blanket, 40% of heat to
  electricity and 15% of electricity reaching the load.
* **The simulator's search** over the route (section 21 of
  [energy.md](energy.md)):

  | Load, best within ten times today | Electricity out ÷ in | The design |
  |---|---|---|
  | A fibre | ~2 | needs tens of Alfvén times of stability |
  | A gas column | ~2 | needs 160–260× convergence |
  | A liner (section below) | 22 | 34 MA, a thin liner (wall 1/12 of its radius), 580 µm of D-T ice, 4.3 T; converges ~170× |
  | The same, ice burnt like an ICF shell (untested) | 37 | |

  The static route gave 0.89. The liner's lead rests on a convergence five
  times NIF's best, in a model with no liner instability.
* **Net energy needs ~70 MA and a pinch that stays stable through its whole
  burn.** No dense pinch has come close.

## Making it testable, and making energy

```bash
python3 tools/zpinch_energy.py
```

Numbers: [results/zpinch_experiment/energy.md](results/zpinch_experiment/energy.md).

**Why the fibre stalls.** The limit is the load, not the machine. A dense,
cold column is too resistive to take the current, and it collapses
radiatively. Three levers don't help:
* **Keeping the current rising** (NRL's stabilisation) fails whatever the
  driver. The column cools, turns resistive, and the current peaks at
  10–23 MA of 60 and falls.
* **Longer fibres** draw even less current: a 3 cm fibre at 60 MA gives
  0.05 of breakeven, against 1.8 for 1 cm.
* **Recovering the field's energy** after the burn needs a different
  machine, an inductive store. Krakowski's study finds that recovering half
  of it halves the fusion gain needed. Nobody has shown it at this scale.

**A different load: a pre-formed gas column.** A centimetre-wide column of
gas, as in the gas puffs Z already fires:
* it conducts easily and is heated by its own implosion, so it needs no
  preheat laser;
* calibrated on Z's deuterium gas-puff record (3.9×10¹³ neutrons at
  15–17 MA), the model puts that record between one and two Alfvén times of
  stability;
* but it converges too far: 0.1–0.2 mm against the 2.5 mm measured.

**Stage 1, possible now: deuterium gas columns on Z with a 1% tritium
trace.** Gas puffs are Z's routine load, and deuterium needs no new licence.
At 26 MA, a 1 cm column, 2 cm long:

| Fuel load | Held ~1 Alfvén time | 3 | 10 or more |
|---|---|---|---|
| 0.5 mg/cm | 4×10⁸ D-D neutrons, stagnates at 0.6 mm | 2×10¹⁵, 0.1 mm | up to 10¹⁷, 0.06 mm |
| 2 mg/cm | 3×10¹⁰, 0.4 mm | 1×10¹⁵, 0.13 mm | 3×10¹⁵, 0.03–0.08 mm |
| 7 mg/cm | 9×10¹⁰, 0.4 mm | 3×10¹³, 0.16 mm | 3×10¹³, 0.15 mm |

It measures the two numbers the energy path turns on:
* **how many Alfvén times** a dense pinch holds: from neutron yield, X-ray
  imaging and the duration of the stagnated column;
* **how far it converges:** from the stagnation radius.

The D-T/D-D ratio from the tritium trace gives the ion temperature.

**The decision it makes.**
* If Z's columns still stagnate near 2.5 mm with at most ~10¹⁴ neutrons, a
  bare pinch can't reach the convergence energy needs. Then the fuel belongs
  inside a liner, which is MagLIF. Sandia is already modelling a frozen
  deuterium fibre inside a MagLIF liner.
* If they reach ~0.1 mm and 10¹⁵ neutrons or more, the bare-pinch path is
  open.

**Stage 2: D-T.**
* Z's 2024 diagnostic papers say it is beginning to use D-T fuel, though no
  D-T yield has been published. Nevada's "Sodium" plasma focus runs D-T.
* On Z, a D-T gas column that holds 3–30 Alfvén times would release
  3–9 MJ per cm. That is more than the 5.7 MJ its bank holds at the load,
  though still only 4–12% of the electricity at the wall.

**Stage 3: energy, on a driver above 60 MA.**
* The only one being built is Pacific Fusion's Demonstration System (>60 MA,
  80 MJ stored, target 2030). NNSA's successor to Z is stalled.
* D-T gas columns at 60–80 MA, holding 3–10 Alfvén times:

| Design | Fusion per cm | Burnt | Electricity out ÷ in, 15% of electricity to the load | at 24% (Z 300 design) |
|---|---|---|---|---|
| 60 MA, 3 Alfvén times, 7 mg/cm in 1 cm | 341 MJ | 14% | 0.85 | 1.4 |
| 80 MA, 10 Alfvén times, 20 mg/cm in 1 cm | 2.4 GJ | 35% | 3.3 | 5.4 |

* **The hidden demand is convergence.** These columns squeeze 160–260
  times, to 40–130 µm. Z's gas puffs stop at 2.5 mm, and the best implosion of
  any kind (NIF) converges about 35 times.
* **Fifty-fifty D-T is needed.** At 20% tritium the gain falls about
  eightfold. Deuterium alone gives about 0.002.
* **Against MagLIF:** Pacific Fusion's published D-T ice MagLIF designs give
  109–380 MJ at 56 MA and 816 MJ at 65 MA. The bare gas column is in the same
  class only if it holds its convergence. Earlier simulations of D-T Z-pinches
  also put breakeven at 40–80 MA.

**Stage 4: a power plant.**
* It needs a driver that fires repeatedly: Sandia's LTD cavities have fired
  tens of thousands of shots at ~0.1 Hz, but a plant needs ~10⁹.
* It needs a transmission line rebuilt every shot. Sandia's Z-IFE study
  assumed 3 GJ shots at 0.1 Hz, with recyclable lines costing a few dollars.
* It needs a tritium-breeding blanket, and ideally energy recovery.

### With a liner

The research pointed to a liner, so the model now has one (`liner_um` in
`zpinch.py`; numbers in
[energy.md](results/zpinch_experiment/energy.md#with-a-liner-maglif)):
* **The setup:** a beryllium tube carries the current and implodes onto
  laser-preheated fuel gas with an axial field, as in MagLIF. Optionally a
  layer of frozen D-T lines its inner wall.
* **The field:** it leaks out as the fuel converges. Its growth is calibrated
  to MagLIF's measured field × radius, 0.25–0.34 T·m.
* **The ice:** heat and alphas that reach it ablate it into the hot fuel.
* **Instabilities:** no current flows in the fuel, so there is no m = 0 mode.
  The liner's own instability is charged as a demand on convergence.

**Calibration.** On Z's MagLIF shots (19–20 MA, deuterium), the model's
clean yields are 9–12× the clean 2-D simulations:
* it runs too hot, 7–10 keV against the 2.5–3.3 keV measured;
* its burn is far too brief, 0.04–0.12 ns against ~2 ns.

The likely cause: an incompressible liner bounces off the fuel at once, while
a real one decelerates layer by layer. At high current it errs the other
way. It gives about a tenth of the published 1-D yields for Slutz's gas
designs, and doesn't ignite Sefkow's 70 MA ice design. It matches Pacific
Fusion's 56 MA design (109 MJ) only if the ice burns like an ICF shell. So
treat its numbers as a comparison between designs, not as forecasts.

**What it finds:**

| Liner design, 50/50 D-T | Electricity out ÷ in | Fusion per cm | Convergence | Compressed field |
|---|---|---|---|---|
| 45 MA, best with ice | 0.1 | 25 MJ | 48 | 4.4 kT |
| 60 MA, 3 mm liner, 600 µm of ice | 3.5 | 1.5 GJ | 131 | 63 kT |
| 80 MA, 6 mm liner, 600 µm of ice, 10 T | 4.2 | 2.9 GJ | 43 | 3.7 kT |

* **Thick ice burns at high current, by the model's own ablation.** Below
  ~50 MA the field that holds the heat in also keeps it from the ice. At
  60–80 MA the burning fuel ablates the ice fast enough for the burn to
  spread.
* **At 80 MA, the liner's demands are modest except the current.** Its
  convergence (43) is close to MagLIF's measured 30–40, and its compressed
  field is below the 7.1 kT shown on NIF. The bare gas column needed
  160–260× convergence for less.
* **If the ice burns like an ICF shell** (untested in a liner), the best is
  6 at 60 MA.

**Can it be made to generate energy?** In the model, yes. The most credible
form is a liner with thick D-T ice on a driver of ~80 MA:
* it converges no further than MagLIF already does;
* it needs no stabilisation of the fuel, because no current flows in it;
* the driver doesn't exist yet. Pacific Fusion's >60 MA machine is being
  built for ~2030, and 80 MA is beyond it;
* fifty-fifty D-T is needed.

The simulator's search goes further: 22× the electricity in at only 34 MA.
But that needs a liner thinner than any fielded (wall 1/12 of its radius) and
~170× convergence, exactly where this model is blindest. The 80 MA design,
at MagLIF's own convergence, is the one to believe.

The bare gas column is the alternative: it needs no liner, but it must
converge 160–260× and hold 3–10 Alfvén times. Stage 1 on Z can test the gas
column with deuterium now. Z's existing MagLIF programme tests the liner.

## Limits

* **One zone.** The column is uniform. Real radiative collapse is uneven, and
  Haines finds that ion heating by m = 0 modes stops it at high current.
  The sub-micron collapses here are the ideal.
* **No corona and no beam-target neutrons**: the processes that dominated the
  1980s experiments.
* **The driver** is a fixed 15 nH and a 100 ns rise; real machines differ.
* **One centimetre.** Longer fibres change the inductance and the
  end losses.
* **Gas columns converge too far.** They stagnate at 0.1–0.2 mm, where Z's
  gas puffs stop at 2.5 mm. The energy designs need 160–260× convergence, so
  the simulator charges it as stretch.
* **The liner is incompressible and 0-D.** Its stagnation is far too brief
  (0.05–0.15 ns against ~2 ns measured). It overshoots Z's MagLIF yields 9–12×
  and undershoots published high-current designs. The ice burns only by a
  simple ablation model.
* **Sheared flow** is now doubtful even in theory. Gyrokinetic simulations of
  FuZE find no full stabilisation, and the kink needs trans-Alfvénic flow.
