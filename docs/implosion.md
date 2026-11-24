# Experiment: imploding a hollow shell onto its centre

The original proposal hoped the laser would steer ions inside the crystal into
head-on collisions. [experiments.md](experiments.md) found that this doesn't
happen: ⁶Li and D have the same charge-to-mass ratio, and ⁷Li collisions stay
below about 90 keV. Two crystals exploding into each other gave real head-on
collisions, but too few.

The strongest form of "guided head-on collisions" is **spherical focusing**.
Put the fuel ions on the inside wall of a spherical cavity and pull them all
towards its centre at once, so that they meet in a tiny focus. Two published
schemes do exactly this:

* **Microbubble implosion** ([Murakami, Arefiev & Zosa, Sci. Rep. 8, 7537
  (2018)](https://www.nature.com/articles/s41598-018-25594-3)). A
  10²⁰–10²² W/cm² laser fills a sub-micron bubble in solid hydrogen with hot
  electrons. Their charge pulls the bubble-wall ions inwards, and the ions pile
  up at the centre at up to ~10⁵ times solid density. Follow-ups
  ([Phys. Plasmas 26, 043112 (2019)](https://doi.org/10.1063/1.5093043);
  [Zosa & Murakami](https://arxiv.org/abs/2512.03379)) aim at ultrahigh fields
  and proton beams, not fusion.
* **The periodically oscillating plasma sphere (POPS)**
  ([Barnes & Nebel, Phys. Plasmas 5, 2498 (1998)](https://pubs.aip.org/aip/pop/article-abstract/5/7/2498/444378);
  [Park et al., PRL 95, 015003 (2005)](https://www.osti.gov/etdeweb/biblio/20696429)).
  Ions oscillate in the harmonic well of a uniform electron cloud and all
  reach the centre together every half period. It is a refinement of the
  inertial-electrostatic "fusor", which [Rider (Phys. Plasmas 2, 1853
  (1995))](https://www.osti.gov/biblio/165283) argued cannot reach net energy.

Code: `nanofusion/improved/implosion.py` (also `python -m nanofusion
implosion --help`). Runs and analysis: `tools/implosion_experiments.py`.
Numbers: [results/implosion/summary.md](results/implosion/summary.md).

```bash
python3 tools/implosion_experiments.py
```

## The model

* **Hot electrons.** A fixed, uniform negative charge −Q_e fills the cavity
  of radius R.
  * Inside it the field is linear, so every ion that starts at rest reaches
    the centre after the same time, π/(2ω) with ω² = zeQ_e/(4πε₀ m R³),
    wherever it starts. That makes it a perfect lens.
  * The density is 10²⁸ electrons/m³, typical of hot electrons from a
    10²⁰–10²¹ W/cm² laser and a sixth of the solid's electron density.
* **Ions.** Fully stripped fuel ions sit in a layer on the cavity wall, taken
  at random from the solid's lattice sites.
  * Their total charge is Q_e/κ. The **charge ratio κ** sets how many ions
    each electron pulls.
  * Every ion feels every other ion (the same C kernel as the crystal model),
    so the runs show what limits the focus.
* **Fusion in flight** uses the crystal model's binary-collision estimator,
  with a reaction volume smaller than the focus (checked in section F).
* **Energy in** is counted two ways.
  * The ions' peak kinetic energy: the most optimistic "input".
  * The electrostatic energy of the electron cloud with its balancing positive
    charge at the wall, kQ_e²/(10R). The laser has to supply at least this,
    and the gain against it is the fairer measure.
* **Not modelled.** Electrons don't move. Real hot electrons partly follow the
  ions inwards, and they leak away and cool within about a picosecond.

## Short answer

**The focusing works, but the ions' own repulsion caps the fusion, and the
best case is no better than a crystal firing into a cold catcher.**

* **The lens is good.** The ions of each species arrive at the centre
  together, at the predicted moment, and the focus shrinks to between 36% and
  0.5% of the cavity radius.
* **Space charge sets the limit.** Only N ≈ 4πε₀ r_c E/e² ions of energy E fit
  into a focus of radius r_c before they stop each other. So an ion crossing
  the focus meets at most 4ε₀E/(e² r_c) partners per m². In the runs it meets
  40–100% of that.
* **That bounds the gain.** The fusion energy from the focus, per unit of ion
  energy, is at most **0.5 pm / r_c for D-T** and 8 fm / r_c for D-D.
  * Breaking even on the ions' own energy needs a D-T focus under 0.5 pm
    across, about a hundredth of an atom.
  * The sharpest focus here is about 0.5% of the cavity radius.
* **Packing in fewer ions sharpens the focus but doesn't help.** The number of
  reactions in flight stayed at 2–5×10⁻¹³ (8 nm cavity) while the number of
  ions changed 150-fold.
* **The best gain against the energy the laser must supply** is about
  **3×10⁻⁴**, even with D-T. 97% of it comes from ions springing back out into
  the cold wall, which is the same beam-target fusion as a crystal firing into
  a catcher. The focus itself contributes 5×10⁻⁵ of the ions' energy.

## Results

### A. What stops the focus (D₂, 8 nm cavity, varying the electron-to-ion charge ratio κ)

| κ | ions | focus radius / R | smooth shell would stop at | partners crossed / space-charge bound | D+D in flight | D+D in the wall afterwards |
|---|---|---|---|---|---|---|
| 2 | 10,724 | 0.36 | 0.37 | 1.1 | 2×10⁻¹³ | 7×10⁻¹¹ |
| 5 | 4,290 | 0.17 | 0.17 | 0.83 | 2×10⁻¹³ | 9×10⁻¹² |
| 10 | 2,145 | 0.097 | 0.092 | 0.74 | 2×10⁻¹³ | 9×10⁻¹³ |
| 30 | 715 | 0.042 | 0.032 | 0.68 | 3×10⁻¹³ | 10⁻¹⁴ |
| 100 | 215 | 0.017 | 0.010 | 0.46 | 5×10⁻¹³ | 10⁻¹⁶ |
| 300 | 72 | 0.008 | 0.003 | 0.41 | 2×10⁻¹³ | 2×10⁻¹⁸ |
| 1000 | 22 | 0.005 | 0.001 | too few ions to sample | – | – |

* **Up to κ ≈ 10** the shell stops exactly where a smooth charged shell
  would: (√(1+4/κ) − 1)/2 of the radius. That is pure space charge.
* **At larger κ** the layer's graininess (the ions are discrete) blurs the
  focus.
* **The wall column.** After the bounce, the compressed shell springs back out
  (the "proton flash" of the microbubble papers). Ions reach the wall with up
  to 2–5 times the energy the electrons gave them. With many ions (small κ)
  most of the fusion happens there, not in the focus.

### B. Scaling to bigger cavities

At fixed κ and electron density, a cavity s times larger has s³ ions with s²
times the energy, so the pair exposures move to s² times the energy with s⁴
times the weight. The runs check this:
* κ = 10, R = 4–13 nm: the focus stays at 0.095–0.104 R. Each run predicts the
  next larger one to within a factor of 2.
* κ = 100, R = 6–18 nm: the focus tightens slightly with size (0.020 → 0.015 R),
  and the prediction holds within a factor of 3.

### C. Initial ion temperature (8 nm, κ = 100)

| temperature | 0 | 0.01 eV | 0.1 eV | 1 eV | 10 eV |
|---|---|---|---|---|---|
| focus radius / R | 0.017 | 0.017 | 0.020 | 0.032 | 0.084 |

A cold solid (meV zero-point motion) costs nothing. Any heating from a
prepulse quickly blurs the focus.

### D. Fuel: deuterium and tritium don't arrive together

* **Arrival time depends on mass.** In the harmonic well it is π/(2ω), and
  ω ∝ √(z/m). So tritium arrives √(3/2) ≈ 1.22 times later than deuterium
  (29.8 against 37.5 fs in the 8 nm run).
  * Each deuteron then crosses 55 times less tritium than deuterium.
  * Moving the layers doesn't fix this, because the arrival time doesn't
    depend on starting radius.
* **⁶Li³⁺ and D do arrive together**, since they have the same charge-to-mass
  ratio, but D+⁶Li fusion is too weak at these energies.
* **Best case.** Assuming D and T could somehow be made to arrive together, I
  applied D–T cross sections to the deuterium dynamics.

Best gain for cavities small enough that no ion exceeds 4.9 MeV (the end of
the measured cross sections), with an electron density of 10²⁸ m⁻³:

| fuel | best cavity radius | gain vs field energy | gain vs ion energy | of which from the focus | from ions springing back into the wall |
|---|---|---|---|---|---|
| D₂ | 270 nm | 2×10⁻⁵ | 1.4×10⁻⁴ | 4×10⁻⁷ | 100% |
| D-T as it is | 105 nm | 2×10⁻⁴ | 1.4×10⁻³ | 1.3×10⁻⁶ | 100% |
| ⁶LiD | 180 nm | 9×10⁻⁶ | 6×10⁻⁵ | 6×10⁻⁷ | 99% |
| D-T arriving together (best case) | 80 nm | 2.6×10⁻⁴ | 1.8×10⁻³ | 4.6×10⁻⁵ | 97% |

The focus-only gain for the best case is a third of the space-charge bound
for its focus size (0.042 R = 3.4 nm, giving 0.5 pm / 3.4 nm = 1.4×10⁻⁴).

### E. Repeated passes

Hot electrons last about a picosecond, roughly eight oscillations at this
density (one period is 117 fs for deuterium). Over four passes (8 nm, κ = 100) the focus blurs from 0.017 R to
0.034, 0.044 and 0.062 R, as close encounters at each focus heat the ions.
So repeated passes add a few times the first pass's yield, not orders of
magnitude.

### F. Numerical check

Halving or doubling the reaction volume (0.01–0.05 nm radius) changes the
in-flight D+D estimate by less than a factor of 1.5. The two estimators,
summing the cross section directly and using the exposure histogram, agree to
within 30% here. At the lowest energies, where the cross section climbs
steeply within one histogram bin, they differ by up to 2×. Energy is conserved
to 10⁻⁵ in every run.

## Why it can't be fixed

* **More ions or fewer?**
  * Fewer ions per electron (larger κ) gives a sharper focus, but the ions
    then carry only about 5/κ of the field energy.
  * More ions (smaller κ) stall at a fraction of the radius.
  * Either way, the gain against the field energy stays around 10⁻⁴ or less.
* **The only way past the space-charge bound is to neutralise the focus.**
  Electrons would have to be compressed along with the ions.
  * The hot electrons can't do it: to fill the cavity they must be hotter than
    the potential that would pull them into the focus, so they gather there
    only a little.
  * Cold electrons compressed with the ions are a neutral plasma being
    compressed. That is ordinary inertial-confinement fusion. It needs
    ρR ≈ 0.3 g/cm², about 10 µg of hot fuel, which is trillions of times more
    than the imploding layer of an 80 nm bubble.
* **The same as the fusor.** This is the nanoscale version of the critique
  Rider made of electrostatic-confinement fusors in 1995. Ions scatter and
  space charge builds up long before they fuse often enough.

## Caveats

* **Fixed electrons.**
  * Real hot electrons would partly follow the ions in. That helps the focus,
    but only by a Boltzmann factor of order one, as above.
  * They also leave or cool within about a picosecond, which ends the
    oscillation.
* **Field energy is a lower bound on the laser energy.** The hot electrons'
  kinetic energy and the laser-to-electron conversion (10–30% at best) both
  come on top, so real gains are lower still.
* **Extrapolation.** Large cavities are extrapolated from 4–18 nm runs. The
  scaling was checked to within a factor of 2–3 across that range.
* **Non-relativistic ions; cross sections held constant above their fitted
  ranges.** The best-gain search stops where any ion would exceed 4.9 MeV.
* **Electron density.** At a fixed ion energy:
  * the focus part of the gain grows as √(electron density);
  * the dominant wall part doesn't change.

  Even the solid's own electron density (6×10²⁸ m⁻³) would raise only the
  focus part, about 2.4 times.
