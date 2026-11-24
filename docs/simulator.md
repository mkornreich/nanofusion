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
| Pulsed power → fuel of an X-ray-driven capsule | 0.35% | on Z ~0.2% of the stored energy reached a capsule |
| Pulsed power → fuel of a magnetised liner (MagLIF-like) | 0.5% | 13.5 kJ of hot-spot heat from ~20 MJ stored |
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
| Compressed magnetic field | 7.1 kT (NIF magnetised implosion, 2023) |
| Cold squeezed density | 1,000 g/cm³ (NIF design) |
| Radial squeeze of the fuel | 35 (NIF-class implosions) |
| Hot-spot areal density and pressure | 0.44 g/cm², 600 Gbar (NIF's best shots) |
| Areal density of the cold fuel | 0.6 g/cm² (NIF's best shots) |
| Fuel adiabat (lower is harder) | 1.5 (NIF's 2012 low-foot designs) |
| Electron screening potential | 1.5 keV (⁶Li(d,α) in palladium) |
| Pinch current, energy into the pinch | 26 MA, ~3 MJ (Sandia's Z) |
| Stable Alfvén times of a dense hot pinch | 1 (fibre pinches broke up at peak current) |
| Muons per second | 10⁸ (PSI, J-PARC) |
| Liquid or solid D-T density | 1.45 × liquid hydrogen (solid D-T) |
| Projectile speed | 80 km/s (dust accelerators) |
| Crystals' distance from a compressed core | ≥ 0.1 mm |
| Uranium blanket k_eff | 0.98 |

## Untested theory

Some ideas have never been shown to work at all, whatever the technology.
They are kept apart from stretch:
* a design that relies on one is flagged with the idea's name;
* the main search leaves such designs out;
* `search --theory` (or `search(..., theory=True)`) lets them in, so the
  leaderboard shows what would happen if the theory worked.

| Untested idea | Where | Source |
|---|---|---|
| Spin-polarised D-T that survives in the fuel (cross section ×1.5) | `spin` on the D-T routes | Kulsrud et al. 1982; never tested in a plasma |
| Alpha channelling: waves hand the alphas' energy to the ions | `alpha_channelling` on `opaque_pB`, `son_fisch_pB` | Fisch & Rax 1992; claims "more than half" |
| Direct conversion of fusion products at 60% | `direct_conversion` on the same | designs claim 60–70%; ~48–65% shown on low-energy ions and beams |
| Stuck muons freed by field re-acceleration (sticking 0.1%) | `rf_stripping` on the muon routes | Kulsrud; TEK concept. Needs fuel below 25 K: ~30 MeV of cooling per fusion, charged here |
| Antiprotons light a hot spot by fission | `antimatter_spark` | ICAN-II, AIMStar (Penn State) |
| Light-assisted tunnelling | `laser_assisted` | Queisser & Schützhold 2019; Wang 2020 (optimistic), Thomson et al. 2024 (pessimistic) |
| A microbubble implosion | `hollow_shell` | simulated only |
| A dense pinch stable while its current rises | `stability` on `dynamic_pinch` | NRL 1987; 2-D MHD blamed a cold core under an unstable corona |
| Sheared-flow stabilisation at fibre density | `stability` on `dynamic_pinch` | shown for ~10³ Alfvén times only below 0.6 MA and 5×10²⁴ per m³ (FuZE) |
| Burn spreading into a liner's frozen fuel like an ICF shell | `ice_burn` on `dynamic_pinch` | high-gain MagLIF designs (Slutz & Vesey 2012); simulated only |

A design's **stretch** is its largest ratio of demand to best. 1 means today's
best would do; 10 means some demand is ten times beyond it.

## The routes

| Route | Section | What it varies |
|---|---|---|
| `compressed_icf` | 7, 11 | **the reference**: ordinary inertial fusion. Shell density, areal density, adiabat; hot-spot or crystal-spark ignitor; laser or pulsed-power X-ray drive; spin |
| `magnetised_liner` | 12, 19 | MagLIF-like column: radius, length, temperature, density, compressed field, spin |
| `muon_factory` | 9, 10 | exploding crystals make the muons: explosion energy, fuel density, sticking |
| `hollow_shell` | 8 | ion energy, size of the focus |
| `lattice_confinement` | 9 | NASA's measurement, or a perfect photon source |
| `pB_catcher` | 4, 5 | proton energy, boron or decaborane, avalanche, screening, crystals or accelerator |
| `antimatter_spark` | 22 | compressed D-T lit by antiprotons: shell density, areal density, adiabat, CERN-like or theoretical antiproton factory |
| `laser_assisted` | 22 | cool D-T bathed in 5×10²¹ W/cm² light: density, temperature, size, optimistic or pessimistic tunnelling boost |
| `cold_catcher` | 1, 5 | ion energy, tritide or LiD catcher, crystals or accelerator, screening |
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
| `dynamic_pinch` | 24, 25 | a pinch followed in time on a generator circuit: current; a fibre (radius), a gas column (radius, fuel per cm) or a MagLIF-like liner (radius, wall, gas, axial field, ice); tritium share; laser preheat; Alfvén times stable; stability idea |
| `channelling` | 20 | ion energy, aligned beam or exploding crystal, measured or generous values |

Each run also notes which demand limits it.

**The reference route.** `compressed_icf` is ordinary inertial fusion, so that
every other idea can be measured against it:
* an isobaric hot spot (10 keV, ρR 0.3 g/cm²) inside a dense shell;
* the shell's density and adiabat set the pressure;
* the fuel burns the share ρR/(ρR + 7 g/cm²).

With NIF-like inputs it gives 0.15 mg of fuel, a 39 µm hot spot at 590 Gbar
and a target gain of 1.7, against NIF's record 4.1. So it is, if anything,
conservative.

Calibrated on NIF's own shots, the same model designs a practical shot campaign
for NIF: [icf_experiment.md](icf_experiment.md). It uses `icf_design`,
`icf_shell_for_mass` and `icf_shell_density`, which size the target for a given
fuel energy, pressure and adiabat.

**Two time-dependent models** check the static ones. They are too slow to search
over thousands of designs, except the pinch, which the search samples at a
quarter of the usual rate (`Route.effort`).
* `burn1d.py`: a 1-D implosion and burn of a NIF capsule, matching LLNL's 1-D
  yields.
  * It shows that the static `compressed_icf` misses the adiabat's effect on
    pressure.
  * Its ceiling is about half the clean yield of NIF's best design.
* `zpinch.py`: a pinch on a generator circuit, used by the `dynamic_pinch`
  route.
  * A cold fibre can't take the current.
  * A preheated one collapses radiatively.
  * On Z even a fully stable pinch reaches ~1% of breakeven; the static
    `fibre_pinch` assumed a hot Bennett pinch from the start.
  * It also takes a pre-formed gas column, and a MagLIF-like liner with
    frozen D-T ice. Its field leak is calibrated to MagLIF's measured
    field × radius, but as a 0-D model it is a comparison between designs,
    not a forecast.

**Combinations.** Options can be switched on and off: the ignitor and driver
of `compressed_icf`, spin polarisation on every D-T route, screening and an
accelerator instead of crystals on the catchers, the avalanche, ions heated
first in the tamped ball. `search(..., fixed={...})` holds an idea at one
setting and re-optimises everything else around it;
`tools/simulator_search.py` prints that comparison for every idea.

## What the search finds

See [results/simulator/summary.md](results/simulator/summary.md) for the full
tables, including what each idea adds.

| Best electricity out / in | Today | ≤ 10× today | ≤ 100× today | Unlimited |
|---|---|---|---|---|
| Ordinary inertial fusion | – | **1.5** | **3.0** | **7.1** |
| Dense D-T Z-pinch (static) | 0.01 | 0.89 | **2.8** | **5.7** |
| Dense D-T Z-pinch, followed in time (best: a liner) | 0.003 | **22** | **22** | **22** |
| Muon catalysis | 0.12 | 0.16 | 0.16 | 0.16 |
| Best route using exploding crystals | 0.005 | 0.005 | 0.012 | 1.3 |

* **With today's technology nothing reaches breakeven.** The best is muon
  catalysis, at about a tenth of the electricity it costs.
* **Ordinary inertial fusion is the most proven route beyond that.** Within
  ten times today on every demand it reaches 1.5, with a 19 MJ laser and a
  shell about ten times NIF's areal density.
* **The time-dependent pinch leads on paper**, 22 within ten times today.
  * Its best design is a 4.5 mm beryllium liner, its wall 1/12 of its radius,
    lined with 580 µm of D-T ice, with a 4.3 T field, at 34 MA. It burns 25%
    of its fuel.
  * But it converges ~170×, against NIF's best 35, in a 0-D model with no
    liner instability. That is in effect cylindrical inertial fusion, and its
    convergence is the number to doubt.
  * Bare fibres and gas columns reach about 2 within ten times today.
  * The static pinch gives 0.89.
* **Combining the crystals with other ideas doesn't rescue them.**
  * As a fusion spark they are worse than a hot spot.
  * As a muon factory they reach 0.006 of breakeven.
  * Their best route stays at about 1%.
* **Lasers beat pulsed-power X-rays** for driving a capsule.
* **Screening and the avalanche change nothing visible.** Heating the ions
  first helps two to three times.
* **A uranium blanket makes several routes "net energy"**, but nearly all of
  the heat is then fission.
* **Untested theory moves the numbers, not the verdict.**
  * If the liner's ice burns like an ICF shell, the time-dependent pinch
    reaches 37 within ten times today, against 22 without.
  * Spin-polarised fuel lifts inertial fusion and the pinch by 15–35%.
  * Proton–boron with alpha channelling and direct conversion reaches 2% of
    breakeven.
  * Antimatter sparks reach 0.39, needing 22 times the most antiprotons ever
    stored, even from a factory 10,000 times better than CERN's.
  * Light-assisted tunnelling reaches 0.85 only with no limits on technology.
  * Freeing stuck muons with fields: 0.14 against 0.12.

## Limits

* The simulator is only as good as the models it wires together. Each has the
  limits listed in its section of energy.md.
* The search found four flaws in the models, now fixed:
  * a pinch's fuel burning past 100%;
  * pinch pressures far beyond any experiment;
  * a pinch too cool for its alphas to make up its X-ray losses;
  * a capsule's hot spot far above the pressure of the shell around it.

  Other models may still have flaws the search hasn't exposed.
* Stretch compares each demand separately. Meeting several at once is harder
  than any one alone.
