# Experiments: trying to get net energy

[experiments.md](experiments.md) showed that laser-exploded LiD crystals can
make plenty of *fusion reactions*, but the fusion energy stayed below the laser
energy by a factor of at least 1,000. This page tries to close that gap:
* three tests of increasing generosity (sections 1–3);
* every published idea I could find for beating the limits (sections 4–11):
  proton–boron fuel and its "avalanche", electron screening, a check against a
  real record experiment, fast ignition of compressed fuel, spherical
  focusing of ions in a hollow shell, muon-catalysed fusion, NASA's lattice
  confinement fusion, three more speculative routes (cold pycnonuclear
  compression, crystals as a muon factory, spin-polarised fuel), and a map of
  fusion in dense matter across every burning regime.

Full numbers are in [results/energy/summary.md](results/energy/summary.md).
Code: `nanofusion/improved/energy.py`. To reproduce:

```bash
python3 tools/energy_experiments.py
```

## Short answer

**No version of this scheme gives net energy unless it turns into ordinary
inertial-confinement fusion.** That means heating or compressing a few mm of
deuterium-tritium fuel, with megajoules delivered in under a nanosecond. At
that point the nanocrystals are at best an inefficient way of heating it. The
published literature reaches the same conclusion for every related idea.

| Test | Best fusion energy out / energy in | Verdict |
|---|---|---|
| Ions stopping in a cold catcher, perfect laser | 0.004 (D into tritide); 0.0014 (p into boron) | Never: an ion in cold matter returns at most ~0.4% of its energy |
| ... with a real laser (1–15% of the pulse becomes explosion energy) | ≤ 5×10⁻⁴ | Never |
| Ions into an already-hot D-T plasma (3–30 keV) | 1–2.8 | Only with tritium, and only inside an existing fusion plasma |
| Protons into hot p-B plasma (100 keV) | up to 0.42 at 4 MeV (agrees with a 2026 study's ≤ 0.40) | Never above 1 |
| Proton–boron "avalanche" (alphas knock protons into more fusions) | k ≤ 0.02 secondary fusions per fusion (literature: ~10⁻²) | No chain reaction (needs k ≥ 1) |
| Electron screening ("lattice-enhanced" fusion) | 0.0041 → 0.0043 even with an extreme 3 keV screening | No effect where fusion happens |
| Whole ball of crystals burning as a plasma | Breakeven at 1.6–4 mm of D-T, 11–78 MJ in | Possible on paper, but it is uncompressed ICF |
| Same, with proton–boron or ⁷LiD fuel | never: X-ray losses exceed fusion power | Never |
| Crystals as the ignition spark for compressed D-T (fast ignition) | target gain ~1.4 if the crystals sat 30 µm from the core; ≤ 0.14 at buildable distances | Not buildable: an isotropic explosion wastes > 97% of the ions |
| ... with shaped (double-layer) crystals that make a beam | target gain ~1.3 at 100 µm, but it needs a 23 MJ, 6 fs ignitor laser (4×10²¹ W) | Not buildable: 400,000 times today's most powerful laser |
| Hollow shell imploding onto its centre (microbubble implosion, POPS) | 3×10⁻⁴ of the field energy, even for D-T; the focus itself gives 5×10⁻⁵ of the ion energy | Never: the ions' own repulsion caps how many fit in the focus |
| Muon-catalysed fusion in solid or compressed D-T | 0.5 of the muon's cost at the record (150 fusions per muon); at most 0.8 (1.2 counting blanket heat); laser-made muons 10⁻⁸ | Closest of all, but sticking caps it below electrical breakeven, which needs ~1,000 fusions per muon |
| NASA's lattice confinement fusion (gamma rays on erbium deuteride) | 4×10⁻¹⁴ measured; ≤ 3×10⁻⁷ even with a perfect photon source | Never: knocked-on deuterons stop long before they fuse; screening adds < 20% |
| Squeezing a cold D-T crystal until it fuses (pycnonuclear) | ideal gain ~30, but only at ~10⁷ g/cm³ | Needs 10⁴ × NIF's density, ~6×10⁸ Gbar and a 360-fold radial squeeze, all without heating |
| Exploding µm crystals as a muon factory for muon catalysis | 0.08–0.55 of the laser light | No better than an accelerator, and needs 2–15 × today's record laser intensity |
| Cold or warm dense fuel, all burning regimes (screened, pycnonuclear), four rate models | hot fuel ignites more cheaply at every density in every model; liquid metallic hydrogen burns 1% in > 10⁴⁶ years | Never: no regime where screening or zero-point motion beats heat |
| Spin-polarised D-T (cross section × 1.5) | breakeven energy of a crystal ball 11 → 0.4 MJ (generous), 78 → 23 MJ (conservative); hot-plasma ion gain 3.1 → 4.7 | Biggest single lever found, but still uncompressed fusion delivered in picoseconds |

### Where this sits among real results (fusion energy / laser energy)

| System | Ratio | Source |
|---|---|---|
| NASA lattice confinement fusion, 2020 (fusion / electron-beam energy) | ~4×10⁻¹⁴ (measured) | [NASA/TP-2020-5001616](https://ntrs.nasa.gov/citations/20205001616) |
| This project: 30 nm LiD crystals, ions stopping in LiD | ~10⁻⁷ | [experiments.md](experiments.md) |
| Texas Petawatt deuterium-cluster record (1.9×10⁷ n, 90–180 J) | ~2×10⁻⁷ (measured); reproduced by this model, section 6 | [PRL 111, 055002](https://arxiv.org/abs/1308.3798) |
| This project: same crystals, tritide catcher | ~3×10⁻⁵ | [experiments.md](experiments.md) |
| Best laser proton–boron shots (PALS, LFEX), ~10⁸ α per joule | ~10⁻⁴ | [PRE 101, 013204](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.101.013204) |
| Muon-catalysed fusion record (Los Alamos, 150 fusions per muon; fusion / muon cost) | ~0.5; electricity out / in ~0.14 for an optimised plant | [Kelly et al. 2021](https://doi.org/10.1088/2515-7655/abfb4b) |
| National Ignition Facility, 7 Apr 2025: 8.6 MJ from 2.08 MJ | 4.1 (target gain); ~0.02 from the wall plug | [LLNL](https://lmf.llnl.gov/science/achieving-fusion-ignition), [AIP Adv. 15, 095013](https://pubs.aip.org/aip/adv/article/15/9/095013/3361996) |
| Needed for a power plant | target gain ~60–140 (laser efficiency × gain ≳ 10) | [National Academies 2013](https://www.nationalacademies.org/read/18289/chapter/4) |

## 1. Cold catcher, with a perfect laser

Assume the laser turns *all* its energy into ion energy, and the ions stop in
cold material.

| Ion energy | D into LiT (tritide) | D into ⁶LiD | Li-6 into ⁶LiD |
|---|---|---|---|
| 100 keV | 0.0027 | 4×10⁻⁶ | 10⁻¹¹ |
| 300 keV | 0.0037 | 4×10⁻⁵ | 7×10⁻⁸ |
| 1 MeV | 0.0024 | 6×10⁻⁴ | 10⁻⁵ |
| 3 MeV | 0.0044 | 0.003 | 10⁻⁴ |

* **Why it can't work.** A fast ion in cold matter spends almost all its
  energy kicking electrons, and only rarely fuses. The best it can return is
  about 0.4% of its energy. The ~1–2% values at 10–20 MeV are outside the
  measured cross sections.
* **And the laser isn't perfect.** For a 6 fs pulse at the stripping
  threshold, only 1% (10 nm crystal), 4% (30 nm) or 13% (100 nm) of the laser
  energy passing through the crystal becomes explosion energy.

## 2. Fire the ions into a plasma that is already hot

Hot electrons drag on ions far less: the slowing time grows as T_e^1.5. So an
ion gets many more chances to fuse before it stops. This is the
"neutral-beam-driven" fusion idea from tokamak reactor studies. Ions from a
30 nm crystal (mean deuteron about 160 keV) into a solid-density plasma:

| Electron temperature | D into D-T plasma | D into ⁶LiD plasma | Stopping distance (D-T) |
|---|---|---|---|
| 1 keV | 0.3 | 4×10⁻⁴ | 0.5 mm |
| 3 keV | 1.0 | 0.0013 | 1.3 mm |
| 10 keV | 2.3 | 0.003 | 2.7 mm |
| 30 keV | 2.8 | 0.004 | 3.1 mm |
| 100 keV | 0.6 | 0.0015 | 0.8 mm |

* **Energy multiplication is real, but only with tritium.** In D-T at
  10–30 keV the ions return 2–3 times their energy. In a lithium deuteride
  plasma they return under 0.5%, because D+Li reactions are far too slow at
  these energies.
* **It doesn't solve the problem.** It needs a millimetre-sized, solid-density
  plasma that is already at fusion temperature. Making and holding that plasma
  costs far more than the ions return. That is exactly the problem fusion
  energy has to solve in the first place.

## 3. A whole ball of crystals, exploded at once

**The model.** A ball is packed with crystals and the laser strips all of
them.
* **Energy in:** the crystals' explosion energy (the minimum the laser must
  supply) plus heating the electrons to 1 keV.
* **Fusion out:** fast ions slowing down, then a thermal burn until the ball
  flies apart (time ≈ radius / 4 sound speeds), minus X-ray (bremsstrahlung)
  losses.
* **Two bounds for the fast-ion phase:** the plasma is at its final
  temperature at once (generous), or the ions slow down on the initial cold
  electrons (conservative).
* **Search:** crystal size and fill fraction are chosen to minimise the
  energy needed.

| Fuel | Bound | Best crystals | Plasma temperature | Breakeven ball radius | Fuel mass | Energy in at breakeven |
|---|---|---|---|---|---|---|
| D-T (solid D-T crystals) | generous | 50 nm | 31 keV | 1.6 mm | 3 mg | 11 MJ |
| D-T | conservative | 30 nm | 11 keV | 4.2 mm | 59 mg | 78 MJ |
| ⁶Li(D,T) (tritium in LiD) | both | 12 nm | 13 keV | 18 mm | 20 g | 26 GJ |
| ⁶LiD (no tritium) | generous / conservative | 100 nm | 850 keV | 16 / 32 cm | 13 / 100 kg | 1 / 10 PJ |
| ⁷LiD (the original's isotope) | both | – | – | never: X-ray losses exceed fusion power | – | – |
| Decaborane (proton–boron) | both | 100–200 nm | 290–1150 keV | never: fusion is only ~0.4 of X-ray losses | – | – |

D-T, 30 nm crystals, conservative bound:

| Ball radius | Energy in | Fusion out / energy in | of which fast ions |
|---|---|---|---|
| 0.1 mm | 1 kJ | 0.03 | 0.007 |
| 1 mm | 1 MJ | 0.25 | 0.007 |
| 3 mm | 28 MJ | 0.73 | 0.007 |
| 1 cm | 1 GJ | 2.4 | 0.007 |
| 3 cm | 28 GJ | 7.2 | 0.007 |

**What this shows**

* **The gain comes from ordinary thermal burn, not from the crystals.** It
  grows with the ball's size (its areal density ρR), just as the Lawson
  criterion says. The fast ions from the crystals add under 1% in the
  conservative case. The crystals are only a way to heat the fuel.
* **Without tritium there is no realistic path.** ⁶LiD breaks even only at
  about 850 keV and petajoules. That is beyond the model's validity: the cross
  sections there are extrapolated, and relativistic X-ray losses are
  underestimated. ⁷LiD never breaks even.
* **This is uncompressed ICF.** The cheapest case is 11–78 MJ for a few mm of
  D-T at normal solid density. Real laser fusion compresses the fuel about a
  thousandfold instead. The energy needed for a given ρR falls as 1/density²,
  which is how the National Ignition Facility reached net fusion gain in
  December 2022 with about 2 MJ of laser. Exploding nanocrystals only push
  matter outward; they cannot compress anything.

**Why even the paper breakeven can't be built**

* **Penetration.** Every crystal in a millimetre-scale ball would have to be
  stripped at about 2–6×10¹⁸ W/cm². But light only gets about 10 nm into
  stripped LiD or D-T, so there is no known way to reach crystals more than a
  few tens of nanometres deep.
* **Timing.** The energy has to arrive faster than the ball flies apart (about
  0.3–0.7 ns), which means more than 10¹⁷ W sustained for nanoseconds.
  * The highest-energy laser, the National Ignition Facility, delivers about
    2 MJ at about 5×10¹⁴ W.
  * Petawatt lasers reach higher power, but only for femtoseconds to
    picoseconds and with kilojoules, not megajoules.
* **Laser efficiency.** Only 1–13% of the laser energy becomes explosion
  energy, so the laser would need 10–100 times the energies in the table.

## 4. Proton–boron fuel (p + ¹¹B → 3 α), and the "avalanche"

Proton–boron fusion makes no neutrons, which is why several laser-fusion
companies (HB11 Energy, Marvel Fusion) and the NAPLIFE project pursue it with
nanostructured targets.

**Setup.** I added the p–¹¹B cross section (Nevins & Swain form as refitted in
[arXiv:2601.00241](https://arxiv.org/html/2601.00241v1): 0.1 b at 148 keV,
1.08 b at 612 keV), and decaborane (B₁₀H₁₄) crystals as the fuel.

**Cold material.** A proton stopping in cold decaborane returns at most 0.14%
of its energy as fusion. Exploding decaborane crystals return at most 0.02%,
even for 200 nm crystals.

**Hot plasma.** Protons into a solid-density p–B plasma do better:

| Electron temperature | 0.6 MeV proton | 2 MeV | 4 MeV |
|---|---|---|---|
| 1 keV | 0.002 | 0.002 | 0.002 |
| 10 keV | 0.036 | 0.035 | 0.030 |
| 30 keV | 0.11 | 0.13 | 0.12 |
| 100 keV | 0.22 | 0.35 | 0.42 |

This matches an independent 2026 Monte Carlo study: "additional energy from
fast fusions is unlikely to exceed 40% of the initial proton beam energy",
best at 4 MeV ([arXiv:2604.06769](https://arxiv.org/abs/2604.06769)).

**The avalanche claim.** Proposals ([Eliezer et al. 2016](https://ui.adsabs.harvard.edu/abs/2016PhPl...23e0704E/abstract))
argued that each fusion's alphas knock protons into the 600 keV resonance and
trigger about 3 more fusions, a chain reaction. I computed the knock-on chain
(Rutherford scattering, plus a ×5 allowance for nuclear scattering):

| Medium | Secondary fusions per fusion, k |
|---|---|
| cold solid decaborane | 4×10⁻⁸ |
| hot solid-density p-B plasma, 30 keV | 0.004 |
| hot, 100–300 keV | 0.016–0.017 |

A chain reaction needs k ≥ 1. Published analyses agree:
* [Belloni 2021](https://arxiv.org/abs/2012.14533) finds k ~ 10⁻².
* The 2026 study finds no realistic regime supports a chain reaction.

**Whole ball of exploded decaborane crystals.** Its fusion power is only about
0.4 of its X-ray power, so it never breaks even.

**Nanostructures don't help.** The main critique of nanostructured p-B
schemes ([Lackner et al.](https://arxiv.org/abs/2302.01628)) is exactly
section 1's limit: in solid or sub-solid material, electrons stop the protons
long before they fuse. Only strong pre-compression helps, and nanostructures
can't provide it. The NAPLIFE claims of fusion in gold-nanorod-doped polymers
rest on spectroscopy without detected nuclear products, and haven't been
replicated ([review](https://arxiv.org/abs/2609.01366)).

## 5. Electron screening ("lattice-enhanced" fusion)

**What it is.** Electrons bound in a solid partly cancel the Coulomb barrier,
by a screening potential U_e. Measured values ([review](https://www.mdpi.com/2571-712X/7/3/50)):
* about 0.3 keV in many metals;
* about 1.5 keV for ⁶Li(d,α) in palladium;
* only ~25 eV in D₂ gas.

NASA's "lattice confinement" experiments in erbium and titanium deuterides
reported neutrons but no energy balance
([PRC 101, 044610](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.101.044610)).

**Result.** Deuterons into a cold tritide catcher:

| U_e | 3 keV ion | 30 keV | 300 keV | best over all energies |
|---|---|---|---|---|
| none | 1.7×10⁻¹¹ | 1.1×10⁻⁴ | 0.0037 | 0.0041 |
| 1 keV | 3.6×10⁻⁹ | 1.4×10⁻⁴ | 0.0037 | 0.0042 |
| 3 keV (beyond any measurement) | 7×10⁻⁷ | 2.3×10⁻⁴ | 0.0037 | 0.0043 |

Screening multiplies slow-ion fusion enormously, but slow ions barely fuse
anyway. Where the energy balance is decided (100–300 keV), it changes nothing.
The laser-stripped crystals also have no bound electrons to screen with.

## 6. Check against a real experiment: the Texas Petawatt record

To make sure this framework isn't simply pessimistic, I reproduced the
deuterium-cluster fusion record: 1.9×10⁷ neutrons from 90–180 J, ion
temperature 18 keV ([Bang et al.](https://arxiv.org/abs/1308.3798),
[PRE 87, 023106](https://journals.aps.org/pre/abstract/10.1103/PhysRevE.87.023106)).
I used the experimenters' own burn model with this project's D–D reactivity:

| Atom density | 1 mm³ hot volume | 10 mm³ |
|---|---|---|
| 10¹⁸ cm⁻³ | 4×10⁵ neutrons | 9×10⁶ |
| 3×10¹⁸ cm⁻³ | 4×10⁶ | 9×10⁷ |
| 10¹⁹ cm⁻³ | 4×10⁷ | 9×10⁸ |

* The measured 1.9×10⁷ sits in the middle of the published density and
  volume ranges.
* The Coulomb-explosion model needs 17 nm clusters to give the measured
  14 keV mean ion energy.
* The measured fusion-to-laser energy ratio, about 2×10⁻⁷, is the same as
  this project's estimates for exploding crystals.

## 7. Crystals as the ignition spark for compressed fuel (fast ignition)

**Idea.** The one published route where laser-made ions could deliver net
energy is ion fast ignition:
1. Compress D-T to about 400 g/cm³ with a main laser.
2. Then put about 10 kJ of ions into a ~20 µm spot within ~20 ps
   ([Fernández et al. 2014](https://iopscience.iop.org/article/10.1088/0029-5515/54/5/054006);
   [Atzeni et al. 2002](https://www.osti.gov/etdeweb/biblio/20243458)).

Designs of this kind aim for target gains around 100, for example 250 kJ of
compression plus 80–100 kJ of ignitor
([HiPER](https://www.osti.gov/etdeweb/biblio/21120388)).

**Could exploding crystals be the ignitor?**
* Deuterons need about 4 MeV to stop within ~1.2 g/cm² of the compressed fuel
  (this model's hot-plasma stopping at 5 keV).
* That calls for ~120 nm LiD crystals, stripped at about 10²¹ W/cm².
* Of their explosion energy, 13% is in deuterons near 4 MeV, and the laser
  turns about 15% of its energy into explosion energy.
* But the explosion is **isotropic**, so only a small fraction of the ions
  head for the spot:

| Crystals' distance from the spot | Ions aimed at the spot | Ignitor laser for 10 kJ | Target gain (30 MJ ÷ all laser energy) |
|---|---|---|---|
| 30 µm (not buildable) | 2.6% | ~20 MJ | ~1.4 |
| 100 µm | 0.25% | ~200 MJ | 0.14 |
| 1 mm | 0.0025% | ~20 GJ | 0.0014 |

*Correction:* an earlier version of this table gave 11%, 1% and 0.01%, and
gains of ~6, 0.6 and 0.006. It used (r/d)² for the share of ions heading into
the spot. The correct solid-angle fraction is (1 − cos θ)/2 ≈ (r/d)²/4, so
those earlier numbers were four times too generous.

**Why it isn't buildable.** Crystals 30 µm from a core imploding at hundreds
of km/s would be swept into it, behind dense plasma the ignitor laser cannot
cross. Real fast-ignition designs assume *directional* ion beams accelerated
from foils at ~10% efficiency. Those are themselves undemonstrated: today's
foil sources typically reach 1–5%
([review](https://arxiv.org/abs/2609.19571)). An isotropic Coulomb explosion
is a far worse source.

### Can shaped crystals make a beam?

I ran Coulomb explosions of 11 crystal shapes in the improved model and
measured where the deuterons go. Code: `tools/shaped_crystals.py`. Numbers:
[results/shaped/summary.md](results/shaped/summary.md).
* **A LiD disk does not beam.** Thick or thin, its deuterons spread almost
  evenly: 1.1–1.8 times the isotropic share within 10° of its axis.
  * The disk's mid-plane is a symmetry plane, so ions near it get no push
    along the axis.
  * A disk one atomic layer thick explodes entirely sideways, as a flat fan.
* **A double-layer target does beam.** This is the asymmetric design used for
  collimated laser-ion beams ([Esirkepov et al., PRL 89, 175003
  (2002)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.89.175003)):
  a lithium-only disk with a deuterium "dot" on one face.
  * Its deuterons come out within about ±10% of one energy.
  * There is a trade-off in the dot's size:
    * a small dot sends 40% of its deuterons within 10° of the axis, but they
      carry only 1.5% of the explosion energy;
    * a dot covering the whole face carries 18%, but only 5% lands within 10°.

    The heavy disk takes the rest of the energy.

Scaled to the ~4.4 MeV deuterons the ignitor needs, with every crystal aimed
at the spot and 6 fs pulses (the most efficient case):

| Source distance | Best shaped crystal | Sphere | Ignitor laser (best) | Its peak power | Target gain (best) |
|---|---|---|---|---|---|
| 30 µm (not buildable) | double-layer | 2.4 | 2.3 MJ | 4×10²⁰ W | 12 |
| 100 µm | double-layer | 0.2 | 23 MJ | 4×10²¹ W | 1.3 |
| 300 µm | double-layer | 0.02 | 0.2 GJ | 3×10²² W | 0.15 |
| 1 mm | double-layer | 0.002 | 2 GJ | 4×10²³ W | 0.013 |

The "Sphere" column is the sphere's target gain, for comparison.

* **Shaping helps 5–10×, but not enough.** At any buildable distance the
  target still barely breaks even.
* **The ignitor laser would have to deliver megajoules in femtoseconds.**
  That is 10²⁰–10²³ W, 40,000 to 40 million times the most powerful lasers
  today (~10¹⁶ W). Longer pulses lower the power but waste proportionally more
  energy.
* **The laser can't fully strip them.** At full scale the best crystals are
  20–45 nm thick, beyond the ~11 nm the laser reaches into stripped LiD.

## 8. Spherical focusing: a hollow shell imploding onto its centre

**Idea.** This is the original "guided head-on collisions" idea in its
strongest form. A laser fills a nanoscale cavity in the fuel with hot
electrons, and their charge pulls the wall's ions to the centre all at once.
That is [microbubble implosion](https://www.nature.com/articles/s41598-018-25594-3),
and with repeated oscillations the
[periodically oscillating plasma sphere](https://pubs.aip.org/aip/pop/article-abstract/5/7/2498/444378).
Full details are in [implosion.md](implosion.md).

**Result (24 runs, 22 to 10,700 ions).**
* The lens works: the focus shrinks to 0.5–36% of the cavity radius.
* But only N ≈ 4πε₀ r_c E/e² ions fit into a focus of radius r_c. So the fusion
  energy from the focus is at most 0.5 pm / r_c times the ion energy for D-T.
  Breaking even would need a focus a hundredth of an atom across.
* Deuterium and tritium don't even reach the centre together: their
  charge-to-mass ratios differ, so tritium arrives 22% later.
* The best case, assuming D and T could be made to arrive together, gives
  3×10⁻⁴ of the field energy the laser must supply. 97% of that comes from ions
  springing back into the cold wall, which is ordinary beam-target fusion.
* Getting past the bound means neutralising the focus with compressed
  electrons. That is inertial-confinement fusion again.

## 9. Cold fusion in a crystal: muon catalysis and lattice confinement

These two routes need no hot plasma, which makes them the last candidates
for fusion inside a crystal. The models are literature-based, not this
project's ion simulation. Code: `nanofusion/improved/alternatives.py` and
`tools/other_routes.py`. Numbers:
[results/other/summary.md](results/other/summary.md).

### Muon-catalysed fusion in solid D-T

**How it works.** A negative muon replaces an electron and binds a deuteron
and a triton 207 times closer than an ordinary molecule. They fuse within a
picosecond, and the muon is usually freed to do it again. It stops when the
muon decays (2.2 µs) or sticks to the alpha particle, about 0.45% of the time
in the record experiment.

**Measured inputs.**
* The record is about 150 fusions per muon, at Los Alamos
  ([Jones et al. 1986](https://www.osti.gov/biblio/6978999)).
* Sticking was measured at 0.45% (Los Alamos), 0.56% (PSI, direct), and 0.65%
  in theory.
* A muon costs 4.7 GeV in an optimised accelerator design
  ([Kelly et al. 2021](https://doi.org/10.1088/2515-7655/abfb4b)), 5 GeV in
  the traditional estimate, and 8 GeV in a pessimistic one.
* The best laser-driven source made about 280 negative muons from a 20 J shot
  ([ELI-NP, 2025](https://arxiv.org/abs/2503.20904)).

**Results.**

| | Fusions per muon | Fusion energy / muon cost (4.7–5 GeV) |
|---|---|---|
| Record conditions | 150 | 0.53–0.56 |
| Fuel 10 times denser than liquid hydrogen, e.g. in a diamond anvil cell | 209 | 0.74–0.78 |
| Unlimited density (sticking is the only loss) | 222 | 0.78–0.83, or 1.2 counting blanket heat |
| Laser-made muons (~4×10⁸ GeV of laser light each) | 150 | 6×10⁻⁹ |

* **This is the closest any cold route comes.** A muon can return about half
  of what it cost.
* **Sticking caps it.** Squeezing the fuel speeds up the cycle, but the muon
  still sticks once per ~220 fusions.
* **Electricity is much further off.** It needs about 1,000 fusions per muon,
  scaling Kelly et al.'s 14% electricity out/in at 150. That means sticking
  five times lower than anything measured.
* **Laser-made muons are hopeless.** They cost about 10⁸ times more than
  accelerator muons.
* **Fission blankets** (μCF neutrons breeding plutonium) can show a net gain,
  but that is fission energy.

### NASA's lattice confinement fusion

**How it works** ([NASA/TP-2020-5001616](https://ntrs.nasa.gov/citations/20205001616);
[Phys. Rev. C 101, 044610](https://journals.aps.org/prc/abstract/10.1103/PhysRevC.101.044610)).
Gamma rays split deuterons in erbium or titanium deuteride. The resulting
neutrons, about 145 keV, knock other deuterons to tens of keV. Those
deuterons fuse with cold lattice deuterons, helped by electron screening.

**What NASA measured.** A 2.9 MeV, 15 mA electron beam (43.5 kW) on 480 g of
ErD₃ gave about 1.6×10³ fusion neutrons per second. That is about 2 nW of
fusion power: **4×10⁻¹⁴ of the beam energy.**

**Model.** I followed each photoneutron through its deuteron collisions and
took each deuteron's thick-target yield, with screening:
* 3×10⁻⁸ fusions per photoneutron in NASA's ~1 cm vials;
* 1.6×10⁻⁷ in an infinite sample;
* screening adds under 20%, even at 3 keV, because the deuterons have tens of
  keV.

Even if every joule of beam became a photon just above threshold, and every
photon split a deuteron, the return would be at most **3×10⁻⁷**. A 14 MeV
fusion neutron doing the knocking would return only 10⁻⁴ of its energy.

This is the cold-catcher limit again: fast deuterons in cold matter lose
their energy to electrons before they fuse.

## 10. Three more speculative routes

Code: `nanofusion/improved/alternatives.py` and `tools/creative_routes.py`.
Numbers: [results/creative/summary.md](results/creative/summary.md).

### Squeeze a cold crystal until it fuses (pycnonuclear fusion)

**Idea.** In a crystal lattice, nuclei tunnel into each other through their
zero-point motion alone. That is how dense white dwarfs burn. Compressing
cold fuel costs only the energy of squeezing its electrons, far less than
heating it.

**Model.**
* The Salpeter & Van Horn (1969) pycnonuclear rate for a static lattice.
* Fuel at zero temperature, held up only by degenerate electrons.
* A ball ignites if it burns 1% of itself before flying apart. Its alpha
  particles then heat it to thermonuclear temperatures.

| D-T density | Burn rate per nucleus | Pressure | Radial squeeze from solid | Smallest ball that ignites |
|---|---|---|---|---|
| 10³ g/cm³ (NIF-class) | 4×10⁻²⁹ /s | 220 Gbar | 17× | bigger than a planet |
| 10⁵ g/cm³ | 0.01 /s | 4.5×10⁵ Gbar | 78× | 3,400 km |
| 10⁶ g/cm³ | 4×10⁵ /s | 2×10⁷ Gbar | 170× | 0.19 m (29,000 tonnes) |
| 10⁷ g/cm³ | 8×10¹⁰ /s | 6×10⁸ Gbar | 360× | 1.6 µm (0.16 mg), 1.8 MJ to compress |

* **On paper, an ideal gain of ~30.** Fusion releases far more than squeezing
  the electrons costs. But only at white-dwarf density, ~10⁴ times what NIF
  reaches.
* **The squeeze is impossible.** It means shrinking the fuel 360-fold in
  radius without heating it. Hydrodynamic instabilities limit real implosions
  to ~30–40-fold, and every real compression heats the fuel by shocks.
* **At NIF's density** the standard rate would have to be 10³⁷ times higher.
* **One dissenting theory.** Ichimaru & Kitamura (Phys. Plasmas 6, 2649 (1999))
  predicted "prodigiously high" p-d rates in liquid metallic hydrogen at
  ≥ 20 g/cm³ and ≥ 1 Gbar. Nobody has compressed hydrogen cold beyond
  ~0.05 Gbar, so this is untested rather than refuted. It is the one place in
  this project where an unexplored corner could still matter.

### Exploding crystals as a muon factory

**Idea.** Muon catalysis came closest of everything here (section 9), but
muons are expensive. They are made by pions from ~1 GeV-per-nucleon ions
hitting a target. A 2–5 µm ⁶LiD crystal would explode to that energy, and
every ion at the same radius gets the same energy per nucleon. The
isotropic explosion is no handicap here: the pion target can surround the
crystals.

**Result.**
* **About 34 GeV of laser light per muon**, if 15% of the light becomes ion
  energy. Even at 100% it is 4.8 GeV, no better than an accelerator.
* **Muon gain 0.08–0.55 against the light**, at 150 fusions per muon. Lasers
  convert electricity to light at ≤ 20%, so against the wall plug it is
  ~0.02–0.1.
* **Still ten million times cheaper** than today's laser-made muons (~4×10⁸ GeV
  each).
* **The laser doesn't exist yet.** The crystals would need stripping at
  4×10²³–2×10²⁴ W/cm², 4–15 times the 2021 record.

### Boosters

**Spin-polarised D-T.**
* Aligning the spins raises the D-T cross section by 1.5. This is theory; the
  DIII-D programme aims to test whether polarisation survives in a plasma.
* It is the biggest single lever found here. A ball of exploding crystals
  breaks even at 0.4 MJ instead of 11 MJ (generous bound), or 23 MJ instead of
  78 MJ (conservative).
* A fast deuteron in a 30 keV D-T plasma returns 4.7 times its energy instead
  of 3.1.
* It doesn't change the verdict: the ball is still uncompressed fusion, with
  its energy needed in picoseconds deep inside, where light can't reach.

**Freeing stuck muons with X-rays.**
* A thesis found that 10²³ photons/cm² would raise muon catalysis to >600
  fusions per muon.
* That fluence carries ≥ 1.8×10⁸ J per cm². Paying for it needs ~10¹⁷ muons
  in each cm² of fuel at once, billions of seconds of the most intense muon
  beam.
* A 2026 analysis ([arXiv:2606.07077](https://arxiv.org/abs/2606.07077))
  finds external-field stripping adds at most ~40%.

## 11. Fusion in dense matter, every regime at once

**The open question.** Section 10 left one corner untested. Ichimaru &
Kitamura (Phys. Plasmas 6, 2649 (1999)) predicted "prodigiously high"
reaction rates in liquid metallic hydrogen at ~20 g/cm³ and ~1 Gbar. More
generally: could strong screening in dense, cool fuel, or zero-point motion
in a crystal, ever ignite more cheaply than ordinary hot fuel?

**Model.** The single formula of [Gasques et al., Phys. Rev. C 72, 025806
(2005)](https://arxiv.org/abs/astro-ph/0506386) covers all five burning
regimes:
* classical thermonuclear;
* strongly screened (Monte Carlo screening);
* thermo-pycnonuclear;
* thermally enhanced pycnonuclear;
* zero-temperature pycnonuclear.

It comes with four model variants: optimal, maximum-rate, minimum-rate, and
one built on Ogata, Iyetomi & Ichimaru's own Monte Carlo screening. In hot
fuel it reproduces the Maxwellian D-T and D-D rates to within 5–14% (40% at
30 keV). At zero temperature it matches Salpeter & Van Horn to within 2%.
Code: `nanofusion/improved/alternatives.py` and `tools/dense_fusion.py`.
Numbers: [results/dense/summary.md](results/dense/summary.md).

**Ichimaru's corner.** D-T at 950 K, at 20 g/cm³ and at 40 g/cm³ (where cold
D-T reaches 1 Gbar), burns 10⁻⁵⁶ to 10⁻⁷⁴ times per nucleus per second,
depending on the model. Burning 1% would take more than 10⁴⁶ years. At 950 K
the ions are far below their quantum "plasma temperature" (~1 eV) and sit in
their ground states. Only zero-point motion helps them tunnel, and at these
densities that is hopeless.

**Hot against cold at the same density.** I looked for the smallest ball that
burns 1% of itself before flying apart and is thick enough to stop its own
alphas (ρR ≥ 0.3 g/cm²). The energy is what the ball holds after
compression and heating; the cheapest hot ball at 100 g/cm³ holds ~20 kJ,
comparable with NIF's ~10 kJ hot spot.

| Density | Cheapest hot ball (≥ 1 keV) | Cheapest warm ball (10 eV – 1 keV) | Cheapest cold ball (< 10 eV) |
|---|---|---|---|
| 100 g/cm³ | 2×10⁴ J at 16 keV | 10¹⁶ J | never |
| 1,000 g/cm³ | 200 J at 16 keV | 10¹⁴ J | never |
| 10⁵ g/cm³ | 0.02 J at 16 keV | 2×10¹¹ J | never (maximum-rate model: a 7 km ball) |
| 10⁷ g/cm³ | 10⁻⁵ J at 25 keV | 6 J – 4 MJ | 7 J – 10¹¹ J, by model |

* **Hot fuel wins at every density, in every model.** Warm fuel's cheapest
  point is always the hottest edge of its range, so there is no sweet spot
  where screening pays.
* **Why.** Strong screening only grows once the ions are cold (Γ ≫ 1). By then
  the bare tunnelling rate has collapsed by far more orders of magnitude than
  screening can restore.
* **Caveat.** These formulas were tuned for carbon in white dwarfs. For
  hydrogen isotopes the true uncertainty is larger than the spread between the
  four models. But the gaps here are tens of orders of magnitude.

## Why everything hits the same wall

**Rider's argument.** Rider's general analysis of fusion out of thermal
equilibrium ([Phys. Plasmas 4, 1039 (1997)](https://www.osti.gov/biblio/527873))
explains all of this at once. Fast ions scatter and slow down many times more
often than they fuse:
* about 100 hard Coulomb scatterings per fusion even at the best energy;
* in cold matter, far more electron collisions on top of that.

So any scheme that relies on fast ions in cold or unconfined fuel must keep
re-supplying more energy than fusion returns.

**The only escape.** Fuel that is hot, dense and large enough, a Lawson-type
condition, burns in thermal equilibrium. That is what the National Ignition
Facility achieves by compression: 8.6 MJ from 2.08 MJ of laser in 2025. Even
there, a power plant needs about 15–35 times more target gain and a far more
efficient laser.

**Later challenges to Rider.**
* [Rostoker et al. 2003](https://ui.adsabs.harvard.edu/abs/2003JFuE...22...83R)
  disputed his conclusion for colliding-beam reactors.
* [Liu et al. 2025](https://arxiv.org/abs/2405.13260) find a narrow net-gain
  window for very hot magnetically confined p-B plasmas.

Neither applies to laser-exploded nanocrystals.

## Conclusion

The nanocrystal approach can make detectable fusion: about 30 reactions per
30 nm crystal per pulse with a tritium catcher. This model reproduces the real
laser-cluster record within its measured range. But nothing tested here gives
net energy:
* **Fast ions in cold matter** return at most ~0.4% of their energy as fusion.
  Proton–boron, electron screening and avalanches don't change that.
* **Hot-plasma tricks** work only with tritium inside an existing fusion
  plasma.
* **Burning a whole volume of crystals** reduces to uncompressed
  inertial-confinement fusion: megajoules to gigajoules delivered faster and
  deeper than any laser can.
* **As a fast-ignition spark**, an isotropic explosion wastes almost all its
  ions. Shaped crystals can beam, but they need lasers far beyond any that
  exist.
* **Spherical focusing**, the best possible version of the original "guided
  head-on collisions", is capped by the ions' own repulsion. It ends up no
  better than a crystal firing into a cold catcher.
* **Cold fusion in crystals without a plasma.**
  * Muon catalysis comes closest: a muon returns about half its cost, and at
    most 0.8 (1.2 counting blanket heat). Sticking keeps it well short of
    producing electricity, and laser-made muons are 10⁸ times too expensive.
  * NASA's lattice confinement fusion returns 4×10⁻¹⁴ of its beam energy.
* **Speculative routes.**
  * Squeezing a cold crystal until it fuses would work on paper, but only at
    white-dwarf density.
  * Crystals could make muons ten million times cheaper than today's laser
    sources, yet no cheaper than an accelerator.
  * Spin-polarised fuel cuts breakeven energies 3–26-fold, without changing
    the verdict.
* **Dense, cold or warm fuel.** In every burning regime and every rate model,
  hot fuel ignites more cheaply than warm or cold fuel at the same density.
  The last open corner, liquid metallic hydrogen, would take > 10⁴⁶ years to
  burn 1%.

## Caveats (generous or conservative?)

* **Generous:**
  * perfect laser-to-ion conversion in tests 1–2;
  * instant, uniform stripping of a whole ball in test 3;
  * ignoring electron energy beyond 1 keV;
  * the generous fast-ion bound.
* **Conservative:**
  * no alpha-particle self-heating in test 3 (it matters only for balls with
    ρR ≳ 0.3 g/cm², which are larger than the breakeven sizes here);
  * D+⁶Li counted for the (d,α) channel only, maybe a factor of 2 low.
* **Model limits:**
  * 0-D burn with fixed temperature and a fixed burn time;
  * cross sections held constant above their measured ranges;
  * non-relativistic bremsstrahlung, which underestimates losses at hundreds
    of keV.

None of these change the conclusion by the orders of magnitude involved.
