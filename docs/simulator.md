# One simulator for every energy route

Every experiment in [energy.md](energy.md) asks the same question: does more
energy come out than goes in? `nanofusion/improved/simulator.py` puts all of
them into one energy account, so that any design can be run, compared with the
others, and searched for the best result.

```bash
python3 -m nanofusion sim list
python3 -m nanofusion sim run fibre_pinch current_ma=26 alfven_times=1 squeeze=6
python3 -m nanofusion sim search --samples 300
python3 tools/simulator_search.py
```

The last command writes [results/simulator/summary.md](results/simulator/summary.md).

## The account

Every route turns its parameters into the same chain:

1. Electricity becomes driver energy: laser light, an accelerator beam, or a
   pulsed current.
2. Part of that reaches the fuel.
3. The route's own physics model turns it into fusion energy. These are the
   models of the sections in energy.md; the simulator only wires them together.
4. A blanket turns fusion into heat:
   * none;
   * lithium, which breeds tritium (+4.78 MeV per D-T fusion, the default);
   * uranium, a subcritical fission blanket (`--blanket uranium --k-eff 0.95`).
5. Heat becomes electricity.

**Electric gain** = electricity out ÷ electricity in. Net energy means more
than 1.

| Efficiency | Default | Note |
|---|---|---|
| Electricity → laser light | 10% | fusion-laser designs; NIF's own is ~0.5% (`--laser-efficiency`) |
| Electricity → pulsed current | 15% | about what Sandia's Z delivers to its load |
| Electricity → accelerator beam | 40% | a design figure; the best today is 18% |
| Light → crystal ions or fuel heat | 15% | best case (sections 7 and 12) |
| Light → energy of squeezed fuel | 1% | NIF-class implosions |
| Heat → electricity | 40% | (`--thermal-efficiency`) |

## Stretch: how far beyond today's technology

Each design also lists what it demands of technology, against the best ever
achieved:

| Demand | Best so far |
|---|---|
| Laser energy per shot | 2.08 MJ (NIF) |
| Picosecond laser energy | 4.4 kJ (NIF ARC) |
| Laser peak power | 10 PW (ELI-NP) |
| Laser intensity | 1.1×10²³ W/cm² |
| Magnetic field on the target | 800 T (laser-driven coil, claimed) |
| Cold squeezed density | 1,000 g/cm³ (NIF design) |
| Radial squeeze of the fuel | 35 (NIF-class implosions) |
| Hot-spot areal density and pressure | 0.44 g/cm², 600 Gbar (NIF's best shots) |
| Pinch current, energy into the pinch | 26 MA, ~3 MJ (Sandia's Z) |
| Stable Alfvén times of a dense hot pinch | 1 (fibre pinches broke up at peak current) |
| Muons per second | 10⁸ (PSI, J-PARC) |
| Liquid or solid D-T density | 1.45 × liquid hydrogen (solid D-T) |
| Projectile speed | 80 km/s (dust accelerators) |
| Crystals' distance from a compressed core | ≥ 0.1 mm |
| Uranium blanket k_eff | 0.98 |

A design's **stretch** is its largest ratio of demand to best. 1 means today's
best would do; 10 means some demand is ten times beyond it. Spin-polarised
fuel, never shown to survive in a plasma, counts as infinite stretch.

## The routes

| Route | Section | What it varies |
|---|---|---|
| `cold_catcher` | 1 | ion energy, tritide or LiD catcher |
| `crystal_ball` | 3 | fuel, crystal size, fill, ball radius |
| `fast_ignition` | 7 | crystals' distance from the spark, shaped crystals |
| `muon_catalysis` | 9 | fuel density, sticking |
| `pycnonuclear` | 10 | cold density |
| `tamped_ball` | 12 | radius, temperature, field, fill, ions heated first, D-T or D-D, spin |
| `projectiles` | 16 | crystal size, energy per deuteron |
| `pauli_beam` | 17 | squeezed density, fuel size, electron drag, spin |
| `son_fisch_pB` | 17 | density, ion temperature, fuel size |
| `muon_squeezed` | 17 | density, core size, muons in the core |
| `opaque_pB` | 18 | density, ion temperature, ball radius |
| `fibre_pinch` | 19 | current, temperature, Alfvén times stable, squeeze, spin |
| `channelling` | 20 | ion energy, aligned beam or exploding crystal, measured or generous values |

Each run also notes which demand limits it.

## What the search finds

See [results/simulator/summary.md](results/simulator/summary.md) for the full
tables.
* **With today's technology nothing comes close.** The best is muon
  catalysis, at about a tenth of the electricity it costs.
* **One route crosses net electricity on paper: the dense D-T Z-pinch.**
  * Within ten times today's best on every demand it reaches ~0.9.
  * Allowed a hundred times, it reaches ~3. That means holding a 100 MA pinch
    stable for about 100 Alfvén times, which no dense, hot pinch has done.
* **The crystal routes stay below a percent of breakeven** within a hundred
  times today's technology. Picosecond laser energy holds them back.
* **A uranium blanket makes several routes "net energy"**, but nearly all of
  the heat is then fission.

## Limits

* The simulator is only as good as the models it wires together. Each has the
  limits listed in its section of energy.md.
* The search found three flaws in the pinch model, now fixed:
  * fuel burning past 100%;
  * pressures far beyond any experiment;
  * a pinch too cool for its alphas to make up its X-ray losses.

  Other models may still have flaws the search hasn't exposed.
* Stretch compares each demand separately. Meeting several at once is harder
  than any one alone.
