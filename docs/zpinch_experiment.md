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
  then, and no D-T MagLIF shot was found since. Only dense plasma focus
  facilities run D-T routinely (Los Alamos: 2 MA, 1.8×10¹² neutrons).

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

  | Budget | Electricity out ÷ in | The design |
  |---|---|---|
  | ≤ 10× today | 0.002 | 3.6 MA, 63 µm fibre, ~9 Alfvén times: it breaks before it burns |
  | ≤ 100× today | 3.7 | 71 MA, 1.5 mm fibre, 46 Alfvén times; burns 40% of its fuel |
  | ≤ 10× today, with sheared flow (untested) | 4.2 | 68 MA, 1.5 mm fibre; the best of any route at that budget |

  The static route gave 0.89 and 2.8 for the first two budgets. Followed in
  time, the pinch needs more stability than the static model asked for, but
  burns better once it has it.
* **Net energy needs ~70 MA and a pinch that stays stable through its whole
  burn.** No dense pinch has come close.

## Limits

* **One zone.** The column is uniform. Real radiative collapse is uneven, and
  Haines finds that ion heating by m = 0 modes stops it at high current.
  The sub-micron collapses here are the ideal.
* **No corona and no beam-target neutrons**: the processes that dominated the
  1980s experiments.
* **The driver** is a fixed 15 nH and a 100 ns rise; real machines differ.
* **One centimetre.** Longer fibres change the inductance and the
  end losses.
