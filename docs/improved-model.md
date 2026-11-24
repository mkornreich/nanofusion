# Can it produce fusion? Results from the improved model

The question: could a laser-stripped lithium-deuteride nanocrystal, as in
Kornreich and Hartmann's "guided nuclear laser fusion" proposal (see
[what-it-simulates.md](what-it-simulates.md)), produce fusion if the physics is
done properly?

## Short answer

**The proposed "guided" mechanism does not produce fusion, but a different
mechanism in the same setup does, though far too weakly to be an energy
source.**

1. **Guiding the deuterons into lithium doesn't work.**
   * With **lithium-6**, the isotope that makes LiD a fusion fuel, there is
     nothing to guide: Li-6³⁺ and D⁺ have almost the same charge-to-mass ratio
     (0.4988 vs 0.4965 e/u), so the laser pushes both the same way.
   * Only **lithium-7** (0.428 e/u), the isotope the original code used, lets
     the laser drive deuterons into their lithium neighbours. Even at
     10²² W/cm², far beyond any mid-infrared laser, colliding pairs reach only
     about 90 keV. That gives about 10⁻¹¹ fusion reactions per crystal per
     pulse.
2. **The crystal blows itself apart first.** Once the electrons are gone, the
   nuclei repel each other so hard that the crystal explodes in about 2 fs,
   faster than the laser can guide anything. This "Coulomb explosion" is what
   the original code's neighbour-only forces left out.
3. **The explosion is the real fusion mechanism.** A fully stripped 30 nm
   crystal launches its deuterons at up to about 300 keV.
   * If they fly into surrounding LiD, that gives about **0.5 fusion reactions
     per crystal per laser pulse**, mostly D+D, half of which emit a neutron.
   * This is the mechanism of the real 1999 deuterium-cluster experiments,
     which detected fusion neutrons.
   * It needs at least about 6×10¹⁹ W/cm² just to strip a 30 nm crystal.
4. **Energy gain is about 10⁻⁷.** The laser energy passing through each
   crystal is about 10 million times the fusion energy it yields. At higher
   intensity the gain only gets worse.

So the setup can make *detectable* fusion, like table-top cluster-fusion
experiments, but not *useful* fusion energy.

## Results

The full scan is in [results/summary.md](results/summary.md), with one JSON
report per case. Unless noted, every run uses:
* a 2 nm-radius sphere (3,911 ions);
* a 6 fs pulse at 3.6 µm, as in the proposal;
* 40 fs of simulated time.

Phase 0° means the electrons are stripped at a field crest, which is where
field ionisation actually happens. 90° means stripping at a field zero, which
gives the ions the largest net push from the laser (the original code's
situation), so 90° results are an upper bound.

### 1. Fusion inside the crystal ("guided" collisions)

| Lithium | Intensity (W/cm²) | Phase | Closest D–Li approach | Highest D–Li collision energy | Reactions per crystal |
|---|---|---|---|---|---|
| Li-6 | 10²⁰ | 0° | 2.04 Å (never closer than at the start) | 0.04 keV | 10⁻²⁴⁹ |
| Li-6 | 10²² | 0° | 0.66 Å | 0.17 keV | 10⁻¹²⁷ |
| Li-6 | 10²² | 90° | 0.34 Å | 0.21 keV | 10⁻¹¹⁷ |
| Li-7 | 10¹⁹ | 0° | 0.40 Å | 0.22 keV | 10⁻¹¹⁷ |
| Li-7 | 10²⁰ | 0° | 0.088 Å | 1.2 keV | 10⁻⁵⁴ |
| Li-7 | 10²¹ | 0° | 0.010 Å | 5 keV | 10⁻²⁷ |
| Li-7 | 10²¹ | 90° | 85 fm | 52 keV | 3×10⁻¹⁹ |
| Li-7 | 10²² | 0° | 77 fm | 88 keV | 6×10⁻¹² |
| Li-7 | 10²² | 90° | 620 fm | 86 keV | 3×10⁻¹¹ |

* The guiding effect is real for Li-7. At 10²² W/cm² deuterons are driven
  to within 77 fm of a lithium nucleus.
* But a D+Li pair needs about 2.2 MeV to touch at 2 fm, and tunnelling at
  90 keV is still unlikely: the D+⁷Li cross section there is about 0.2 mb.
* Even scaled to a 30 nm crystal, with about 400 times more deuterons, this
  is only about 10⁻⁸ reactions per pulse.

### 2. Coulomb explosion, then ions stopping in surrounding LiD

Extrapolated to the proposal's 30 nm crystal:

| Case | Max D energy | D+Li | D+D | Total reactions | Fusion / laser energy |
|---|---|---|---|---|---|
| Li-6, 10²⁰, 0° | 281 keV | 0.05 | 0.41 | **0.46** | 9×10⁻⁸ |
| Li-7, 10²⁰, 0° | 482 keV | 0.04 | 0.82 | 0.86 | 1×10⁻⁷ |
| Li-6, 10²⁰, 90° | 601 keV | 0.40 | 0.95 | 1.4 | 4×10⁻⁷ |
| Li-6, 10²¹, 0° | 284 keV | 0.05 | 0.41 | 0.46 | 9×10⁻⁹ |
| Li-6, 10²², 0° | 299 keV | 0.05 | 0.42 | 0.47 | 1×10⁻⁹ |
| Li-6, 10²², 90° (upper bound) | 9.7 MeV | 1,800 | 410 | 2,200 | 2×10⁻⁵ |
| Li-6, 10²⁰, 0°, 0.8 µm laser | 281 keV | 0.05 | 0.41 | 0.46 | 9×10⁻⁸ |

* **The size extrapolation is robust.** Simulating 1, 1.5, 2 or 3 nm crystals
  (485 to 13,301 ions) and scaling each to 30 nm gives the same answer,
  0.46 ± 0.01 reactions. The extrapolated top deuteron energy (281 keV) is 87%
  of the uniform-sphere value for 15 nm (322 keV).
* **Intensity hardly matters once the crystal is stripped**, because the
  energy comes from the crystal's own charge. A stronger laser only wastes
  more energy, so the gain falls.
* **The 90° rows are laser acceleration, not fusion physics.** They assume the
  electrons vanish at a field zero, which leaves every ion with a net
  MeV-scale push. Field ionisation happens near crests, so treat these rows as
  an upper bound.
* **Wavelength doesn't matter.** 0.8 µm gives the same result as 3.6 µm.

### Conditions the extrapolation assumes (and the proposal would have to meet)

* **Stripping a 30 nm crystal** needs a field of 2.1×10¹³ V/m, about
  6×10¹⁹ W/cm².
  * Below that the crystal keeps electrons and the explosion is much weaker.
    The proposal's implied 10¹⁷–10¹⁸ W/cm² can't even fully strip a 2 nm
    crystal.
  * The most powerful mid-infrared (3.9 µm) lasers of the period gave about
    20 mJ in under 100 fs
    ([Photonics 8, 290](https://dx.doi.org/10.3390/photonics8080290)).
    Focused to a few µm, that is roughly 10¹⁸ W/cm² *(inference)*.
* **Laser penetration.** Fully stripped LiD is about 2,700 times denser than
  the density 3.6 µm light can penetrate; the light reaches only about 11 nm
  in. A 30 nm crystal may not be stripped all the way through.
* **A target for the ions.** The beam–target numbers assume every deuteron
  stops in LiD. In the proposed hollow-core fibre they would mostly hit glass,
  and the yield would drop.
* **Electrons are not simulated.** Real electrons would partly shield the
  explosion and absorb laser energy, making things worse for fusion.

### What would be needed for useful energy

The fusion energy per crystal stays about 10⁷ times below the laser energy
passing through it. The highest ratio found, 2×10⁻⁵, came from the unphysical
90° case at 10²² W/cm².

* At about 0.5 reactions per crystal, a 10 µm laser spot covering about
  9×10⁴ crystals would give about 4×10⁴ reactions per pulse. That is easily
  detectable, and similar to the 1999 cluster experiments.
* But that pulse would need about 0.5 J of laser energy to deliver about
  4×10⁻⁸ J of fusion energy.
* No setting in the scan comes within several orders of magnitude of
  breakeven.


## More experiments

[experiments.md](experiments.md) tries six ways to get more fusion:
* pre-expanded crystals (shock shells);
* colliding crystals;
* a 10 µm laser;
* lithium projectiles;
* a tritium catcher;
* larger crystals.

The crystal never fuses by itself. But driving its exploding ions into a
tritium catcher gives about 33 reactions per 30 nm crystal per pulse, and
larger crystals give far more. The best fusion-to-laser energy ratio found is
about 10⁻³.

## What the improved model changes

| Original (2013) model | Improved model |
|---|---|
| Coulomb force only from neighbours within 2 lattice steps | Every ion feels every other ion (exact direct sum, C/OpenMP kernel), so the stripped crystal can Coulomb-explode |
| Field `E0·sin(ωt)` switched on abruptly at t = 0 and never off, which gives every ion a permanent push | Gaussian pulse built from a vector potential, with zero net area (as any real light pulse). The carrier phase at the moment of stripping is a parameter |
| "Lorentz" term missing the field amplitude | Magnetic force q·v×B with B = E/c |
| Explicit Euler, every value rounded to 6 significant digits | Velocity-Verlet in double precision, with a time step that adapts to approaching pairs |
| "Fusion" = a pair closer than 2 fm, a classical criterion that ignores tunnelling | Measured fusion cross sections, which include tunnelling, applied (a) to every D–Li and D–D pair inside the crystal during the run and (b) to every deuteron flying out and stopping in surrounding cold LiD |
| Li-7 only | Li-6 (enriched, the fusion fuel in LiD) or Li-7 |
| No validity checks | Checks whether the laser can strip the crystal at all, and reports the fusion energy against the laser energy |

Kept from the proposal: every electron is removed at t = 0, and the nuclei
start at rest on LiD lattice sites (2.04 Å apart). Electrons are **not**
simulated, so this is the most favourable reading of the proposal.

### Data and methods

* D+D: Bosch & Hale, *Nucl. Fusion* 32, 611 (1992). D+⁶Li → 2α: bare S(0) =
  21.7 MeV·b ([Phys. Rev. C 91, 014601](https://arxiv.org/abs/1406.5735)).
  D+⁷Li: S(0) = 5.4 MeV·b, the n0 + n1 channels only, so a lower bound
  ([Phys. Rev. C 73, 015801](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.73.015801)).
  Cross sections are held constant above 4.9 MeV, where the formulas stop being
  valid.
* **Inside the crystal**, each D–Li or D–D pair closer than 4.08 Å adds a
  reaction probability σ(E)·v_rel·dt / V, with V the volume of that 4.08 Å
  sphere. E is the pair's centre-of-mass energy, kinetic plus Coulomb.
* **Beam–target**, after the explosion, assumes every deuteron stops in cold
  LiD. That is a best case: in the proposed hollow-core fibre they would
  mostly hit glass.
  * Stopping power: Lindhard–Scharff at low energy, Bethe at high energy, with
    I = 36.5 eV ([PDG, LiH](https://pdg.lbl.gov/2024/AtomicNuclearProperties/HTML/lithium_hydride_LiH.html)).
    About 75 keV/µm at its peak, expected accurate to about a factor of 2.
  * Li ions hitting the target are not counted.
* **Extrapolating** from the simulated 1–3 nm crystals to the proposal's 30 nm
  crystal:
  * Each species' mean velocity, the laser drift, is kept as it is.
  * The rest of each velocity, the Coulomb-explosion part, is scaled with
    crystal radius R, so explosion energies grow as R².
  * Ion numbers are scaled as R³.

### Checks

* **Energy conservation:** with the laser off, ion kinetic energy equals the
  drop in Coulomb potential energy to about 10⁻⁵, and to about 10⁻⁸ for small
  crystals.
* **Laser work:** at a stripping phase of 90°, where the laser leaves a net
  drift, the energy the ions gain equals the laser's work to 0.1%.
* **Coulomb explosion against theory:** the fastest deuteron reaches 86%, 92%
  and 94% of the uniform-sphere result q·ρ·R²/(3ε₀) for 1, 1.5 and 2 nm
  crystals (run 200 fs). The gap shrinks with size, which supports the R²
  extrapolation.
* **Real pulse:** the pulse's net field area is zero to 10⁻⁶.
* **Kernel:** the numpy fallback reproduces the C kernel to 10⁻⁹.
* **Convergence:** two key cases rerun with 4× smaller steps give the same
  numbers (see the table).

All of these are in `tests/test_improved.py` or the scan below.
