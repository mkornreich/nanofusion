# Experiments: trying to make it fuse

This follows [improved-model.md](improved-model.md). There, a laser-stripped
30 nm LiD crystal gave about 0.5 fusion reactions per pulse, and only when its
exploding ions hit surrounding LiD. Here I tried six ways to get more:

1. pre-expanded crystals, to cause shock shells;
2. two crystals exploding into each other;
3. a long-wavelength laser for the guided collisions;
4. counting lithium ions as projectiles too;
5. a tritium-bearing catcher;
6. larger crystals.

All runs use the improved model: every ion interacting, a real laser pulse,
fusion from measured cross sections. The raw reports are in
[results/experiments/](results/experiments/) and the size analysis in
[results/experiments/size_scaling.md](results/experiments/size_scaling.md).
To reproduce:

```bash
python3 tools/experiments.py
python3 tools/analyze_experiments.py
```

## Summary

| Idea | Does it help? | Best result |
|---|---|---|
| Pre-expanded crystal (shock shells) | **No.** Its hardest collisions don't grow with size, and it stores less Coulomb energy | Fusion inside the crystal about 10⁻¹⁷ at most, even scaled to 100 nm |
| Lithium-7 guided collisions, 10 µm laser | **No.** Collisions reach about 9 keV (55 keV in the unphysical 90° case) | Fusion inside the crystal ≤ 10⁻¹⁴ |
| Two crystals exploding into each other | **Yes, but small.** Collision energies grow as R² | 10⁻⁴ reactions per pair of 30 nm crystals; about 5 per pair of 200 nm crystals |
| Spherical focusing: a hollow shell imploding onto its centre | **Focuses, but the ions' own repulsion caps the yield** | See [implosion.md](implosion.md) |
| Count lithium ions hitting the catcher | **Somewhat.** Adds 10–25% | — |
| Tritium catcher (LiT instead of LiD) | **Yes, about 60×** for 20–30 nm crystals | 33–43 reactions per 30 nm crystal per pulse |
| Larger crystals | **Yes, strongly.** Explosion energy grows as R² | 10⁴ reactions at 100 nm diameter, 10⁶ at 200 nm |

**Net energy:** see [energy.md](energy.md). None of these routes, nor a hot
catcher or a whole ball of exploding crystals, gives net energy short of
ordinary inertial-confinement fusion.

**Bottom line:** this setup can be made to produce plenty of fusion, through
the crystal's Coulomb explosion driving ions into a fuel-bearing catcher,
especially tritium. Nothing makes the crystal fuse *by itself*, and nothing
comes near energy breakeven.
* The best gain found, fusion energy over the laser energy passing through the
  crystal, is about 10⁻³, for 200 nm crystals.
* That case already needs 3×10²¹ W/cm² and a fully stripped 200 nm crystal,
  both beyond what is physically plausible.

## 1. Shock shells: pre-expanded crystals

**Idea.** Kaplan, Dubetsky & Shkolnikov
([PRL 91, 143401](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.91.143401))
showed that if a cluster's density falls off towards its edge, its Coulomb
explosion forms "shock shells". Inner ions overtake outer ones and collide,
and a [follow-up](https://arxiv.org/pdf/physics/0604136) proposed fusion
inside each cluster this way. A weak prepulse that lets the crystal expand
before the main pulse would create such a profile.

**Model.** The crystal is stretched radially, r' = r·(1 + (E−1)(r/R)²), so
that:
* the centre stays at solid density;
* the edge is 63 times thinner for E = 3;
* no ions start closer together than in the solid.

The shell-crossing condition (local density below a third of the mean inside)
is met in the outer 26–40% of the shells.

**Result.** It doesn't help at these sizes:

| Lithium | Crystal | Hardest D–Li collision | Fusion inside the crystal (D+D), 3 nm |
|---|---|---|---|
| Li-6 | solid | 0.045 keV | 3×10⁻⁸⁹ |
| Li-6 | pre-expanded (E = 3) | 0.035 keV | 8×10⁻¹⁰⁵ |
| Li-7 | solid | 1.3 keV | 1×10⁻²⁴ |
| Li-7 | pre-expanded (E = 3) | 1.8 keV | 5×10⁻²² |

* The overtaking does happen, but the overtaking ions differ in speed only
  slightly.
* Stretching the crystal lowers its stored Coulomb energy, so all ions end up
  slower.
* **The key test:** if these collisions came from the explosion, their energy
  would grow as R². From 1 to 3 nm it grows as R^0 to R^0.5. The R² rule also
  fails to predict the 3 nm runs from the 2 nm runs, by up to a factor of 10⁶.
  So fusion inside a single crystal only grows with the number of ions, R³:
  about 10⁻¹⁷ at 100 nm radius even in the best (Li-7) case.

## 2. Two crystals exploding into each other

**Idea.** Ions from two neighbouring crystals fly head-on into each other at
up to twice the explosion speed.

**Result.** This is the one in-flight mechanism that grows with size.

| Gap (2 nm crystals) | 0.5 nm | 2 nm | 5 nm | 10 nm |
|---|---|---|---|---|
| Hardest D–Li collision | 0.75 keV | 2.3 keV | 6.2 keV | 8.6 keV |
| D+D reactions in the collision | 10⁻²⁴ | 2×10⁻¹⁷ | 3×10⁻¹⁴ | 5×10⁻¹³ |

* With the gap scaled to 2.5 R, collision energies grow as R^1.7 (D–D) and
  R^2.1 (D–Li), as expected.
* The R² rule predicts the 2 nm pair from the 1.5 nm pair to within a factor
  of 1.5, so the extrapolation can be trusted:

| Crystal radius (gap 2.5 R) | 5 nm | 15 nm | 30 nm | 50 nm | 100 nm |
|---|---|---|---|---|---|
| Collision reactions per pair | 10⁻⁸ | 9×10⁻⁵ | 0.01 | 0.2 | 4.5 |

This is real fusion between crystals, but it stays thousands to hundreds of
thousands of times below what the same ions produce when they stop in solid
LiD.

## 3. Long-wavelength laser for Li-7 guiding

**Idea.** At a fixed intensity, the speed difference the laser drives between
D and Li-7 grows with wavelength. So I tried a 10.6 µm CO₂ laser.

**Result.**
* At 10²¹ W/cm² (crest), D–Li collisions reach 8.6 keV, against 5 keV at
  3.6 µm. Fusion inside the crystal is 2×10⁻²².
* Only the unphysical 90° stripping phase (ions "born" at a field zero) gets
  to 55 keV, and then just 10⁻¹⁴.
* The crystal explodes in about 2 fs, before a long-wavelength cycle (35 fs)
  can do much.

## 4–6. Catcher material, lithium projectiles, crystal size

These all use the 3 nm crystal (13,301 ions), extrapolated by radius. The
explosion part of the ion velocities is scaled with R; this rule was validated
earlier from 1 to 3 nm. The laser is set at the lowest intensity that strips
each crystal size (6 fs pulse). Li-6 results; Li-7 is similar.

| Crystal radius | Intensity to strip | Fastest D | Reactions, LiD catcher | Reactions, LiT catcher | Fusion / laser energy |
|---|---|---|---|---|---|
| 5 nm | 7×10¹⁸ W/cm² | 33 keV | 2×10⁻⁵ | 0.003 | 2×10⁻⁷ (LiT) |
| 10 nm | 3×10¹⁹ | 132 keV | 0.02 | 2.9 | 1.5×10⁻⁵ (LiT) |
| **15 nm** (the proposal's 30 nm crystal) | 6×10¹⁹ | 298 keV | **0.57** | **33** | 3×10⁻⁵ (LiT) |
| 25 nm | 2×10²⁰ | 0.8 MeV | 46 | 310 | 4×10⁻⁵ |
| 50 nm | 7×10²⁰ | 3.3 MeV | 1.2×10⁴ | 1.2×10⁴ | 1×10⁻⁴ |
| 100 nm | 3×10²¹ | 13 MeV | 1.3×10⁶ | 1.1×10⁶ | 8×10⁻⁴ |

* **Lithium projectiles.** Li³⁺ ions get three times the explosion energy.
  Hitting the deuterium in a LiD catcher, they add 10% at 15 nm and about 25%
  at 50–100 nm. (They are not counted for LiT; no Li+T data was used.)
* **Tritium catcher.** D+T is the easiest fusion reaction (5 barn at 64 keV),
  so a LiT catcher gives about 60 times more reactions for 15 nm crystals.
  Above about 50 nm the D+⁶Li reaction catches up and the two catchers are
  similar.
* **Size.** Between 15 and 100 nm the yield per crystal rises about as R⁷,
  because the cross sections climb steeply with energy, and the gain about as
  R³–R⁴.
  * On paper, breakeven would come at roughly 1 µm crystals. But those would
    need about 3×10²³ W/cm² (the stripping intensity rises as R²), with ions
    at GeV energies far outside any measured cross section.
  * The laser also only reaches about 11 nm into stripped LiD, so crystals
    much larger than the proposal's would not actually be fully stripped.
  * Treat the 50–100 nm rows as an upper bound, and anything larger as
    meaningless.

## What a real experiment along these lines could give

The most practical version:
* 30 nm ⁶LiD crystals spread on a lithium-tritide (or deuterated) catcher;
* a 10 µm focal spot at about 10²⁰ W/cm², which needs a near-infrared laser
  (Ti:sapphire), not a mid-infrared one.

**Yield.** The spot covers about 9×10⁴ crystals, so about 3×10⁶ fusion
reactions per shot with tritium, or about 5×10⁴ with deuterium. That's easily
measurable, similar to the laser-cluster fusion experiments that followed the
1999 Nature paper.

**Energy.** A 6 fs pulse needs about 0.5 J; longer, more common pulses need
proportionally more. Against that, the fusion energy is about 10⁻⁵ J with
tritium.

**So:** fusion yes, energy no, by about five orders of magnitude even in this
best case. The idea of guiding ions in the crystal into head-on collisions
plays no part in it.

## Caveats

* **Electrons** are not simulated. Real electrons would shield part of the
  explosion and absorb laser energy, which lowers yields.
* **Cross sections** are held constant above their fitted ranges (4.9 MeV for
  D+D and D+Li, 550 keV for D+T). That matters only for the 50–100 nm rows.
* **Stopping powers** use a simple model (factor ~2 accuracy). A target
  without fuel, such as the glass of the proposed fibre, would give almost no
  beam–target fusion.
* **Gains** count only the laser energy that passes through the crystals
  themselves, which is optimistic.
