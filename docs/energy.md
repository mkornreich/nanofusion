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
  fusion in dense matter across every burning regime;
* five more attempts (sections 12–16): holding the ball of crystals together
  with a tamper and a magnetic field, a fission blanket, a chain reaction in
  D-T, bare ions in a magnetic trap, and crystals fired as projectiles;
* one new lever (section 17): squeezing the fuel until its electrons are
  degenerate, so that Pauli blocking weakens their drag. It is tested on
  crystal ions, on hot proton–boron ions and on muon catalysis;
* three walls no earlier section touched (sections 18–20): trapping
  proton–boron's X-rays inside an opaque ball, holding the fuel with a
  current instead of a laser (a fibre Z-pinch), and steering ions through a
  tritide crystal's channels;
* every route in one simulator, with ordinary inertial fusion as a
  yardstick and the ideas combined, searched for the best electricity out per
  electricity in (section 21, [simulator.md](simulator.md));
* untested theory in the same simulator, searched as if it worked
  (section 22).

Full numbers are in [results/energy/summary.md](results/energy/summary.md).
Code: `nanofusion/improved/energy.py`. To reproduce:

```bash
python3 tools/energy_experiments.py
```

## Short answer

**No version of this scheme gives net energy.** The best it can do is turn
into a known kind of fusion target, with the nanocrystals as an inefficient
way of heating it:
* ordinary inertial-confinement fusion: a few mm of deuterium-tritium fuel,
  with megajoules delivered in under a nanosecond (section 3);
* magnetised target fusion: a sub-millimetre of fuel held by a heavy wall
  and a kilotesla magnetic field, with tens to hundreds of kilojoules
  delivered in picoseconds (section 12). This is the biggest lever found,
  about a thousandfold. On paper it still returns only 0.15 of the laser
  energy at breakeven, and 1.5 at best, for 180 MJ of picosecond light.

Squeezing the fuel until its electrons are degenerate (section 17) is the
only lever that weakens the electrons' drag itself. A fast ion then returns
its own energy at the density NIF's fuel is designed for. But the squeeze
costs more than it returns, and fuel that dense ignites the ordinary way from
a small hot spark.

The published literature reaches the same conclusion for every related idea.

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
| Spin-polarised D-T (cross section × 1.5) | breakeven energy of a crystal ball 11 → 0.4 MJ (generous), 78 → 23 MJ (conservative); hot-plasma ion gain 3.1 → 4.7 | Biggest lever on the fuel itself, but still uncompressed fusion delivered in picoseconds |
| Ball of crystals held by a gold tamper and a 3 kT magnetic field (magnetised target fusion without compression) | fusion equals the fuel's energy at 50 kJ instead of 60 MJ (10 kJ if the crystals heat the ions first; 0.1–0.7 MJ with a coil's straight field); 0.15 of the laser energy at breakeven, 1.5 for 180 MJ of laser | Biggest lever found, about 1,200-fold, and still no net energy: it needs a picosecond laser 15–170 times more energetic than any, a field 4 times the strongest laser-made one claimed, and an entrance hole under 1% of the wall |
| ... the same without tritium | deuterium alone: 620 MJ | Never |
| Fission blanket around the fusion source | crystals need k_eff ≥ 0.998; a plasma at breakeven needs 0.92 and muon catalysis only plain uranium, with over 90% of the heat from fission | Never for crystals; otherwise a fission reactor with a costly neutron source |
| Chain reaction in D-T through its own neutrons and alphas | 0.002 new fusions per fusion in cold D-T; about 1 only in a 30–100 keV plasma metres across | No chain reaction in anything that can be built |
| Bare ions in a magnetic trap (no electrons at all) | 7 W per m³ at 10 T; the record collider luminosity would give 0.7 W | Never: the ions' own charge caps the density |
| Whole crystals as projectiles (cluster-impact fusion) | 0.4% of the kinetic energy; 0.0016 from the wall plug | Never: the cold-catcher limit again, though 30 times better than the laser route |
| Crystal ions fired into cold fuel squeezed until its electrons are degenerate (Pauli blocking) | the ion's own energy at 1,000 g/cm³ (250 times the cold-catcher limit), up to ~5 when denser; 0.004–0.035 of the light once the squeeze is paid for | Never: the squeeze eats the gain, and fuel that dense ignites from a joule-to-kilojoule spark anyway |
| Hot proton–boron ions in degenerate electrons (Son & Fisch 2004) | the ions keep themselves hot at 10⁵–4×10⁶ g/cm³; fusion per burn time 3–5.5 times the squeeze and heating | Only on paper, at white-dwarf density, like pycnonuclear D-T |
| Muon catalysis in squeezed fuel (stuck muons knocked free) | sticking below 0.1% above ~60 g/cm³; but ≤ ~140 fusions per muon in a 140 ps squeeze | Never: a squeeze is far too short, and a long one needs 10¹¹ times today's muon beams |
| Proton–boron ball made opaque to its own X-rays | the balls that repay their heat need ρR 17–25 g/cm² (40–60 × NIF's hot spot) and trap 10–65%; trapping half at 1,000 g/cm³ takes 1.8 kg holding ~8 megatons, returning 0.3% | Only at white-dwarf-like densities; trapped X-rays and pairs soak up the energy |
| Fibre of D-T fuel pinched by its own current (Z-pinch, no laser) | at Z's 26 MA, a sixfold squeeze held one Alfvén time (37 ps) repays the plasma's heat; the field costs 7.5 times that heat | Not with today's pinches; the field energy demands a large burn multiplication |
| Every route in one simulator, with the ideas combined (section 21) | today's technology: 0.12 electricity out / in (muon catalysis); within 10× of it on every demand: 1.5 (ordinary inertial fusion), 0.89 (dense Z-pinch); crystal routes ≤ 0.012 within 100× | Only ordinary inertial fusion and the Z-pinch cross 1, and only well beyond today; no combination with exploding crystals gets near |
| Untested theory in the simulator: alpha channelling, direct conversion, freeing stuck muons, antimatter sparks, light-assisted tunnelling, spin-polarised fuel (section 22) | proton–boron 0.023 at best; antimatter sparks 0.39 with 22 × the most antiprotons ever stored; light-assisted tunnelling 0.85 with no limits; muons 0.14 | Moves the numbers, not the verdict |
| Ions channelled through a tritide crystal | 0.008–0.023 for a perfectly aligned beam (2–6 × random); +9% at most for an exploding crystal | Never: the cold-catcher limit, slightly dented |

### Where this sits among real results (fusion energy / laser energy)

| System | Ratio | Source |
|---|---|---|
| NASA lattice confinement fusion, 2020 (fusion / electron-beam energy) | ~4×10⁻¹⁴ (measured) | [NASA/TP-2020-5001616](https://ntrs.nasa.gov/citations/20205001616) |
| This project: 30 nm LiD crystals, ions stopping in LiD | ~10⁻⁷ | [experiments.md](experiments.md) |
| Texas Petawatt deuterium-cluster record (1.9×10⁷ n, 90–180 J) | ~2×10⁻⁷ (measured); reproduced by this model, section 6 | [PRL 111, 055002](https://arxiv.org/abs/1308.3798) |
| Nanowire-array record (2.2×10⁶ neutrons per joule) | ~3×10⁻⁶ (measured); an 8 J shot reproduced by this model, section 6 | [Nat. Commun. 9, 1077](https://www.nature.com/articles/s41467-018-03445-z) |
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

### A second check: nanowire arrays

Arrays of deuterated-plastic nanowires hold the record for D-D fusion per
joule of laser light: 2.2×10⁶ neutrons per joule, about 500 times a flat
target ([Curtis et al. 2018](https://www.nature.com/articles/s41467-018-03445-z)).
Light gets several micrometres down between the wires
([Bargsten et al. 2017](https://www.science.org/doi/10.1126/sciadv.1601558)),
so a volume at a sixth of solid density is heated at once. This is the best
nanostructured fusion target so far.

A later 8 J shot gave 1.2×10⁷ neutrons
([Curtis et al. 2021](https://journals.aps.org/prresearch/abstract/10.1103/PhysRevResearch.3.043181)).
The authors' own plasma simulation gives 3.8×10¹² deuterons above 0.2 MeV and
3×10¹¹ above 1 MeV. Code: `tools/further_routes.py`.
* **Those deuterons stopping in cold plastic give 1.2×10⁷ neutrons in this
  model.** That matches the measurement. The yield is ordinary beam-target
  fusion, as in section 1.
* **The fast deuterons carry 4% of the light.** By their share of the
  charge, all the ions together get about 15%. Nanostructures absorb 70–80%
  of the light ([Park et al. 2021](https://doi.org/10.1063/5.0035174)), but
  most of it stays in electrons.
* **Fusion energy / laser energy is 2–3×10⁻⁶.** That is ten times better than
  gas-jet clusters, and still a few millionths. With a tritide catcher the
  same deuterons would return 1.2×10⁻⁴.

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
muon decays (2.2 µs) or sticks to the alpha particle, about 0.5% of the time.

**Measured inputs.**
* The record is about 150 fusions per muon, at Los Alamos
  ([Jones et al. 1986](https://www.osti.gov/biblio/6978999)).
* Sticking was measured at 0.45% (PSI,
  [Breunlich et al. 1987](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.58.329)),
  0.56% (PSI, direct), and 0.50–0.53% in later runs; theory gives 0.56–0.65%.
  Los Alamos reported lower values at high density, which PSI did not
  reproduce. *Correction:* an earlier version credited the 0.45% to Los
  Alamos.
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
* It is the biggest lever on the fuel itself. A ball of exploding crystals
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

## 12. Holding the ball together: tamper, magnetic field, self-heating

**Idea.** Section 3's ball needed 11–78 MJ because nothing holds it. It flies
apart at its own sound speed, within tens of picoseconds for a
sub-millimetre ball. Fusion research has three known ways to keep hot fuel
together longer without compressing it. Used together they are magnetised
target fusion, and Sandia's
[MagLIF](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.155002)
experiment is the working example.
* **A heavy tamper.** Fuel inside a thick gold shell can expand only as fast
  as it shoves the gold, about 13 times slower than its own sound speed.
* **A magnetic field.** Hot fuel touching a cold wall loses its heat to it in
  picoseconds. A field of a few kilotesla lying along the wall slows that
  heat flow by orders of magnitude. It also holds in the fusion alphas once
  the ball is several of their gyro-radii across (field × radius ≳ 0.5 T m:
  [Basko et al. 2000](https://www.osti.gov/etdeweb/servlets/purl/21328998)).
* **Self-heating.** Trapped alphas reheat the fuel, so a large enough ball
  runs away (ignites).

**Model.** A time-dependent version of section 3's ball:
* separate ion and electron temperatures, and burn-up of the fuel;
* alpha heating, with the nanosecond delay alphas need to slow down at solid
  density;
* X-ray losses and Braginskii heat conduction to the wall;
* a gold wall that recedes at the speed of the shock that the fuel's
  pressure drives into it.

With all three tricks off it reproduces section 3 to within 1%. Code:
`nanofusion/improved/confinement.py` and `tools/confined_ball.py`. For one
ball: `python3 -m nanofusion ball --help`. Numbers:
[results/confined/summary.md](results/confined/summary.md).

**Result.** The cheapest solid-density D-T ball whose fusion output is 1, 10
or 100 times the energy put into the fuel (its heat plus the field inside
it):

| Ball | Fusion = energy in | 10 × | 100 × |
|---|---|---|---|
| Bare, no self-heating (section 3's model) | 65 MJ | 130 GJ | never |
| Bare, with self-heating | 60 MJ | 16 GJ | 2.7 TJ |
| Thick gold tamper, no field | 23 MJ | 2.5 GJ | 8.2 GJ |
| Gold tamper + 1 kT | 160 kJ | 8.4 MJ | 47 MJ |
| Gold tamper + 3 kT | 52 kJ | 4.4 MJ | 44 MJ |
| Gold tamper + 10 kT | 49 kJ | 5.3 MJ | 91 MJ |
| Gold tamper + 3 kT, ions heated first | 10 kJ | 3 MJ | – |
| Gold tamper + 3 kT, straight field | 680 kJ | – | – |
| ... straight field, ions heated first | 110 kJ | – | – |

* **Breakeven falls about 1,200-fold**, from 60 MJ to about 50 kJ, in a ball
  0.6 mm across holding 22 µg of D-T. It doubles its volume in 1 ns instead
  of flying apart in 60 ps.
* **Both tricks are needed.** With a tamper alone the fuel's heat leaks into
  the gold. A field alone cannot hold fuel whose pressure is a hundred times
  its own. The best field is 3–10 kT; above that the field costs as much
  energy as the heat.
* **Heating the ions first cuts it to 10 kJ.** Exploding crystals do this:
  the Texas Petawatt plasma of section 6 had 18 keV ions. Electrons take
  about a nanosecond to catch up, and meanwhile the ions burn at twice the
  temperature for the same energy.
* **A straight field costs about 13 times more.** The numbers above assume
  the field lies along the wall everywhere. Only a ring-shaped field around
  an axis does that, and it takes about 4.5 MA through the fuel. Coils make
  a straight field. Heat then runs along it to the two ends, so the fuel
  must be a column 8 times longer than wide, as in MagLIF.
* **High gain stays expensive.** A tenfold return needs 3–4 MJ and a
  hundredfold one 44 MJ, in a centimetre-sized ball inside centimetres of
  gold. At solid density an alpha needs about a nanosecond to hand over its
  energy, which is as long as the wall holds, so small balls never ignite.
* **Density hardly matters**: energy ∝ density^-0.25, against density^-2 for
  the bare ball. Fuel a thousand times thinner, like a cluster jet, needs
  290 kJ in a 1 cm ball for a microsecond, inside tens of centimetres of
  gold. That is the central point of magnetised target fusion. But exploding crystals need the same intensity
  everywhere, so the thinner the fuel, the more laser power: 12 PW at solid
  density, 4,000 PW at a thousandth of it.
* **Without tritium it stays out of reach.** Deuterium alone, counted
  generously, needs 620 MJ even with the tamper and 3 kT (1.4 TJ bare).
  Lithium deuteride radiates more than it fuses, and holding it longer does
  not change that.

**How fragile is it?** One change at a time, 3 kT ring field
([full table](results/confined/summary.md)):

| Change | Energy in at breakeven |
|---|---|
| none | 52 kJ |
| entrance hole of 1% / 2% / 5% of the wall | 80 / 120 / 340 kJ |
| gold shell 0.3 mm / 0.1 mm thick instead of 1 mm | 55 / 110 kJ |
| wall recedes 1.5 × faster than the shock formula says | 140 kJ |
| fuel as a foam at 20% of solid density, so that light can enter | 85 kJ |
| Bohm heat loss instead of classical; gold mixed in at 1 per 1,000 ions | 52; 54 kJ |
| spin-polarised fuel (section 10) | 18 kJ |

**Check against MagLIF.** MagLIF reaches this kind of plasma by squeezing
it with a 20 MA current: deuterium at 0.4–0.7 g/cm³ and 2.5–3 keV, in a
column about 0.1 mm wide and 6 mm long, held about 2 ns, with
field × radius ≈ 0.34 T m
([Gomez et al. 2014](https://link.aps.org/accepted/10.1103/PhysRevLett.113.155003),
[2020](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.125.155002);
[Knapp et al. 2022](https://www.osti.gov/pages/servlets/purl/1871989)).
* For the published conditions of the 2014 shot, this project's D-D
  reactivity gives 5×10¹² neutrons (10¹² to 1.5×10¹³ over the temperature
  error). Measured: 2×10¹².
* For the best shot (2020): 2×10¹³ (10¹³ to 3×10¹³). Measured: 1.1×10¹³.
* That best shot is equivalent to 2 kJ of D-T fusion from fuel holding
  14 kJ of heat, a ratio of 0.15. The designs above are the same kind of
  plasma at 20–40 keV instead of 3 keV.

**Check against published estimates without a field.** Section 12 is close
to a published proposal. [Ruhl & Korn](https://arxiv.org/abs/2302.06562)
(Marvel Fusion) would heat uncompressed fuel inside a high-Z case with
laser-driven nanostructures, with no magnetic field.
* They find about 3.4 MJ for fusion to exceed the fuel's energy, in a
  1 mm × 1 mm cylinder with D-T at 0.5 g/cm³ and 20 keV. This model gives
  3.6 MJ returning 0.77 for that design.
* Their critics ([Lackner et al.](https://arxiv.org/abs/2302.01628), with a
  [follow-up](https://arxiv.org/abs/2312.13429)) find that perfectly tamped
  D-T needs 250–700 MJ to ignite at 0.5–0.3 g/cm³. This model, with a real
  gold wall, gives 0.5–1.2 GJ for a tenfold return.
* The critics exempt schemes with an applied magnetic field. That is the
  case worked out here.

**What it would take.** If 15% of the light became heat in the fuel, as in
the nanowire experiments of section 6, and the crystals heat the ions first.
The power is the least that reaches the crystals' stripping intensity
(4×10¹⁸ W/cm² for 40 nm D-T crystals) across the fuel:

| Design | Laser energy | Least power, and pulse length | Fusion / laser energy |
|---|---|---|---|
| Solid D-T, ring field, breakeven | 66 kJ | 3.8 PW, 17 ps | 0.15 |
| 20% foam with a 1% hole, ring field, breakeven | 340 kJ | 32 PW, 10 ps | 0.15 |
| Solid D-T column, straight 3 kT field, breakeven | 750 kJ | 3.7 PW, 200 ps | 0.15 |
| 20% foam with a 1% hole, ring field, tenfold return | 180 MJ | 1,900 PW, 97 ps | 1.5 |

**Why it still isn't net energy**
* **The gain is far too low.** The breakeven designs return 0.15 of the
  laser energy, and the tenfold one 1.5 for 180 MJ of picosecond light. A
  power plant needs 60–140.
  Inside a fission blanket it could give net electricity on paper, but that
  is section 13's story.
* **The laser doesn't exist.** It must deliver 66–750 kJ in 10–200 ps.
  * The most energetic picosecond-class shot so far is 4.4 kJ in 30 ps
    ([NIF's ARC](https://www.osti.gov/servlets/purl/1343818)). Osaka's LFEX
    gives 2 kJ in 1 ps ([Kawanaka 2016](https://inis.iaea.org/records/b3pe5-97221)).
  * The highest peak power, 10 PW, comes with only 0.25 kJ
    ([ELI-NP](https://doi.org/10.1017/hpl.2022.11)).
  * The 50–70 kJ ignitor lasers designed for fast ignition
    ([HiPER](https://www.osti.gov/etdeweb/biblio/21153113)) were never built.
* **The field doesn't exist.** 3 kT over a cubic millimetre is about four
  times the strongest laser-driven coil field claimed, 0.8 kT
  ([Santos et al. 2015](https://arxiv.org/abs/1503.00247)), and a later
  reanalysis puts such coils far lower (Peebles et al., Phys. Plasmas 29,
  080501 (2022)). About 5 kT has been reached only by squeezing a seed field
  inside an implosion ([Moody et al. 2022](https://www.osti.gov/servlets/purl/1901997)).
* **The light must get in and the fuel must not get out.** The entrance hole
  has to stay under about 1% of the wall: a 0.2 mm hole in a 1 mm cavity.
  The fuel has to be a foam of 40 nm D-T crystals at a fifth of solid
  density, which nobody makes. Nanowire arrays show such heating a few
  micrometres deep
  ([Bargsten et al. 2017](https://www.science.org/doi/10.1126/sciadv.1601558)),
  not a millimetre.
* **Left out of the model**, all of it harmful: gold blown into the fuel by
  X-rays and alphas, the Nernst effect carrying field to the wall,
  instabilities, and the quarter of the ions that start within one
  gyro-radius of the wall in the smallest balls.
* **It is magnetised target fusion.** As in section 3, the crystals are only
  a way of heating.

This is still the largest lever found in this project: three to four orders
of magnitude in energy, where everything else moved the answer by factors of
a few.

## 13. A fission blanket around the fusion source

**Idea.** Every D-T fusion gives a 14 MeV neutron. In a uranium blanket that
neutron releases far more than the fusion itself, and more again if the
blanket is kept just below critical. This is the fusion–fission hybrid
([Bethe 1979](http://large.stanford.edu/publications/coal/references/docs/bethe.pdf);
Livermore's [LIFE](https://web.mit.edu/fusion-fission/HybridsPubli/LIFE_TOFE/Moses.pdf)
design). Accelerator-driven reactors such as Rubbia's
[energy amplifier](https://www.osti.gov/etdeweb/servlets/purl/602571) do the
same with spallation neutrons. Could a blanket turn even a poor fusion
source into net electricity?

**Model.** Code: `nanofusion/improved/alternatives.py` and
`tools/further_routes.py`. Numbers:
[results/further/summary.md](results/further/summary.md).
* Lithium alone breeds the next triton: 22 MeV per fusion in all.
* In uranium a 14 MeV neutron splits uranium-238 directly (about 1.2
  fissions) and multiplies to about 3 neutrons. With no chain at all that is
  262 MeV per fusion, 15 times the fusion energy. A natural-uranium pile
  [measured](https://ncsp.llnl.gov/sites/ncsp/files/2021-11/la10860_ref_139.pdf)
  290 MeV per neutron; this formula gives 259.
* In a blanket of multiplication k_eff, each of those neutrons starts chains
  worth k / (2.5 (1−k)) fissions: 4,800 MeV per fusion at k_eff = 0.95 and
  12,000 MeV at 0.98. Designs stop at 0.95
  ([MYRRHA](https://www.gen-4.org/sites/default/files/2024-09/GIF%20SCK%20CEN%20%20MYRRHA-Prof.%20Hamid_3-21-18.pdf))
  to 0.98 (Rubbia) to stay safely below critical.
* Net electricity needs (driver efficiency) × (fusion / driver energy) ×
  (blanket multiplication) × (40% heat to electricity) ≥ 1. Lasers are taken
  as 10% efficient and accelerators as 40%. Both are design figures: NIF's
  laser manages 0.5% and the best proton accelerator 18%.

| Fusion source | Fusion / driver energy | Electricity per neutron | k_eff for net electricity | k_eff for 4 × the electricity put in |
|---|---|---|---|---|
| 30 nm crystals on a tritide catcher ([experiments.md](experiments.md)) | 3×10⁻⁵ | 5.9 TeV | > 0.9999 | > 0.9999 |
| The same if all the light became ion energy (section 1) | 0.004 | 44 GeV | 0.9978 | 0.9995 |
| Hollow shell implosion (section 8) | 3×10⁻⁴ | 590 GeV | 0.9998 | > 0.9999 |
| Ball of crystals at its breakeven, bare or tamped (sections 3 and 12) | 0.15 | 1.2 GeV | 0.92 | 0.98 |
| Tamped, magnetised ball returning 10 × its fuel energy (section 12) | 1.5 | 117 MeV | 0.11 | 0.79 |
| Muon catalysis at the record 150 fusions per muon (section 9) | 0.53 | 83 MeV | 0 (uranium, no chain) | 0.70 |
| NIF's best shot, with NIF's own laser | 4.1 | 825 MeV | 0.88 | 0.97 |
| NIF's best shot, with a 10% efficient laser | 4.1 | 43 MeV | 0 (uranium, no chain) | 0.41 |
| For comparison: spallation neutrons from a 1 GeV proton beam | – | 125 MeV | – | – |

* **The benchmark is spallation**: a 1 GeV proton makes about 20 neutrons in
  lead ([Bauer 1998](https://www.osti.gov/etdeweb/servlets/purl/20020531)),
  or 125 MeV of electricity per neutron. A fusion source belongs inside a
  blanket only if its neutrons cost about that, or a few times more, since a
  14 MeV neutron is worth about three ordinary ones.
* **Exploding crystals are 350 to 47,000 times too expensive a neutron
  source.** They would need k_eff of 0.998 or more. The margin to a runaway
  chain is the delayed-neutron fraction, 0.2–0.7%. So that blanket is a
  critical reactor in all but name, and the crystals add nothing to it.
* **Any source that returns a tenth or more of its driver energy works on
  paper.** Muon catalysis needs only plain uranium, which is the known
  result for muon-catalysed hybrids (Petrov, Nature 285, 466 (1980); see
  [Tajima & Eliezer](https://www.osti.gov/servlets/purl/6964059)). So would
  section 12's ball, if it could be built.
* **But it is not fusion energy.** At least 91% of the heat is fission, with
  fission's fuel cycle and waste. It is a fission reactor with an expensive
  neutron source.

## 14. A chain reaction in D-T

**Idea.** Section 4 found no avalanche in proton–boron fuel. D-T is a better
case. Its neutron carries 14 MeV and knocks deuterons and tritons on to as
much as 12.5 MeV, and they fuse a hundred times more easily than protons and
boron. The 3.5 MeV alpha knocks ions on too. Does each fusion trigger
another? The question is old (Peres & Shvarts, Nucl. Fusion 15, 687 (1975);
reviewed in [Belloni 2021](https://arxiv.org/abs/2012.14533)).

**Model.** k = fusions triggered per fusion (`energy.dt_chain`).
* The neutron is followed through all its collisions in an infinite medium.
  They are peaked forwards, as measured.
* Alpha knock-ons use Rutherford scattering, × 5 above 0.5 MeV to allow for
  nuclear scattering.
* Each knock-on fuses with the thick-target probability (cold fuel) or the
  hot-plasma slowing-down probability of section 2.
* The D-T cross section above 550 keV follows the measured fall-off (Bosch &
  Hale's high-energy fit). Elsewhere in this project it is held constant,
  which overstates it; that version is shown for comparison.

| Infinite medium | k (measured fall-off) | k (cross section held constant) | Chance that a 5 MeV deuteron fuses |
|---|---|---|---|
| cold solid D-T | 0.0017 | 0.0023 | 0.04% |
| solid-density D-T plasma, 1 keV | 0.05 | 0.05 | 0.6% |
| 3 keV | 0.17 | 0.17 | 2% |
| 10 keV | 0.49 | 0.52 | 7% |
| 30 keV | 0.93 | 1.0 | 15% |
| 100 keV | 1.3 | 1.7 | 29% |

* **Cold D-T: no chain reaction, by a factor of about 600.** Electrons stop
  the knock-ons first, as in every cold route.
* **Hot D-T: k reaches 1, but only in an infinite plasma at 30–100 keV.**
  This is section 2 again: a fast ion in hot D-T returns 1–3 times its
  energy. And the neutron hands over its energy only in fuel many mean free
  paths thick. One mean free path is 5.3 g/cm², or 28 cm of solid D-T.
* **In the balls of section 12, k is about 0.03–0.05**, nearly all from
  alphas, and only if the field holds them in. The neutron adds 10⁻⁴ to
  2×10⁻³. The chain feeds a fire that is already burning. It cannot light
  one.
* **Published work is less generous still.** Peres & Shvarts found cold D-T
  critical only above 3.5×10⁴ g/cm³, and Afek et al. (J. Phys. D 11, 2171
  (1978)) found 4,200 g/cm³ at 14–15 keV. A 2026 Monte Carlo study
  ([arXiv:2604.06769](https://arxiv.org/abs/2604.06769)) finds D-T critical
  only without neutron leakage and above about 200 ρ^-0.23 keV, which is
  290 keV at solid density.

## 15. No electrons at all: bare ions in a magnetic trap

**Idea.** Electrons do the damage in every route so far: they stop fast ions
and they radiate. A cloud of bare deuterons and tritons has neither problem,
and collisions between like charges do not carry ions across a magnetic
field. On paper its energy balance is unlimited. Penning-trap fusion was
proposed on this basis
([Barnes, Nebel & Turner 1993](https://doi.org/10.1063/1.860837)).

**The catch.** The cloud's own charge limits its density to the Brillouin
limit, n = ε₀B² / 2m
([Dubin & O'Neil 1999](https://nnp.ucsd.edu/pdf_files/rmp71_87_99.pdf)).
A column 5 cm in radius and 1 m long, with ions at 60 keV:

| Field | Ion density | Fusion power per m³ | Power of the column | Voltage across it | Time the field must last to repay its own energy |
|---|---|---|---|---|---|
| 10 T (ordinary superconducting magnet) | 10¹¹ cm⁻³ | 7 W | 0.06 W | 1.2 MV | 65 days |
| 45.5 T (strongest steady magnet) | 2×10¹² cm⁻³ | 3 kW | 24 W | 25 MV | 3 days |
| 1 kT (laser-driven coil, lasts ~1 ns) | 10¹⁵ cm⁻³ | 0.7 GW | 5.5 MW | 12 GV | 9 minutes |
| 100 kT (inside a laser plasma, lasts ~1 ps) | 10¹⁹ cm⁻³ | 7×10¹⁶ W | 5.5×10¹⁴ W | 10¹⁴ V | 0.06 s |

* **A steady trap makes watts.** At 10 T the limit is a few billionths of
  air's density. A tokamak reactor makes about 10⁶ W/m³.
* **Laser-made fields are strong enough but far too brief**, by about
  eleven orders of magnitude.
* **It is section 8's limit again.** There the ions' repulsion capped how
  many fit in the focus. Here it caps how many fit in the magnet.
* **Colliding beams** are the same idea in a ring. One watt of D-T fusion
  needs a luminosity of 7×10³⁴ cm⁻² s⁻¹. The record for any collider is
  5.1×10³⁴ ([SuperKEKB, 2024](https://arxiv.org/abs/2507.13617)), which
  would give 0.7 W.
* **Published assessments agree.** Trap schemes that squeeze the ions into
  a dense focus reach high gain only under optimistic assumptions
  ([Chacón et al. 2000](https://doi.org/10.1063/1.1310199)), and holding
  ions out of equilibrium costs more power than they return
  ([Nevins 1995](https://doi.org/10.1063/1.871080); Rider, below).

## 16. Crystals as projectiles (cluster-impact fusion)

**Idea.** Don't explode the crystal: charge it and accelerate it whole with
a voltage. Accelerators turn electricity into beam far more efficiently
than lasers turn it into ion energy.

**The original claim.** [Beuhler, Friedlander & Friedman
(1989)](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.63.1292)
fired singly charged heavy-water clusters of 25–1,300 molecules at
200–325 keV into titanium deuteride. They reported about 10⁻¹² to 10⁻¹¹
fusions per cluster of 100 molecules.
* A cluster shares its energy among all its atoms. At 300 keV each deuteron
  of a 100-molecule cluster has only 0.3 keV.
* This model then gives 8×10⁻³⁹ fusions per cluster. An independent estimate
  at the time gave 4×10⁻³⁹
  ([Carraro et al. 1990](https://authors.library.caltech.edu/records/j00c8-pp470)).
  The claim was 26 orders of magnitude higher.
* The authors withdrew it in 1992
  ([erratum](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.68.2108)):
  light ions contaminating the beam had probably caused the fusion.

**Doing it properly.** A grain holds positive charge up to a surface field
of about 10 V/nm before it sheds ions
([Dzhanoev et al. 2016](https://arxiv.org/abs/1603.08565)). A deuterium
crystal charged that far could be accelerated to the 100–300 keV per deuteron
where D-T fusion peaks. That is generous, since hydrogen ice is weakly bound.

| Crystal diameter | Charge | Voltage for 100 keV per deuteron | Best fusion / kinetic energy | From the wall plug |
|---|---|---|---|---|
| 3 nm | 16 e | 5 MV | 0.004 | 0.0016 |
| 10 nm | 170 e | 18 MV | 0.004 | 0.0016 |
| 30 nm | 1,600 e | 54 MV | 0.004 | 0.0016 |

* **It is the cold catcher again.** Whatever throws the deuteron, the
  catcher's electrons stop it after it has returned at most 0.4% of its
  energy.
* **The driver is far better.** With a 40% efficient accelerator that is
  0.0016 from the wall plug, against about 5×10⁻⁵ for laser-exploded
  crystals. It is still 600 times short, and a plain deuteron beam does the
  same.
* **The speed is far beyond dust accelerators.** They throw micrometre
  grains at up to 80 km/s with 2 MV
  ([Mocker et al. 2011](https://doi.org/10.1063/1.3637461)). This needs
  3,000 km/s.
* **No collective help.** A 30 nm crystal heats its track to roughly 1 keV,
  which would cut the electron drag. But that heat spreads over about half a
  micrometre within the 0.3 ps the crystal takes to stop (a rough estimate,
  not modelled).

## 17. Switching off the electrons' drag: Pauli blocking in squeezed fuel

**Idea.** Every route above fails the same way: electrons drag a fast ion to
a stop long before it fuses (see "Why everything hits the same wall" below).
Squeezing the fuel changes the electrons themselves:
* squeezed hard enough, the electrons are degenerate, with every low-energy
  state full;
* an ion slower than their Fermi velocity can only kick the few electrons
  near the top;
* so the drag per metre stops growing with density
  ([Fermi & Teller 1947](https://journals.aps.org/pr/abstract/10.1103/PhysRev.72.399)),
  while fusion per metre keeps growing.

I tested three uses:
1. crystal ions fired into cold, squeezed fuel;
2. hot ions in degenerate electrons, proposed for proton–boron by Son & Fisch;
3. muon catalysis, where a stuck muon gets knocked free more often.

Code: `nanofusion/improved/degenerate.py` and `tools/pauli_blocking.py`.
Numbers: [results/pauli/summary.md](results/pauli/summary.md).

**Model.**
* **Electron drag:** Lindhard's random-phase dielectric function of a free
  electron gas at zero temperature, integrated numerically.
  * It reproduces Lindhard & Winther's slow-ion friction (their table, to
    0.5%) and Bethe's formula for fast ions.
  * Checked on a real degenerate electron gas, slow protons in aluminium, it
    gives 0.67–1.06 of the measured stopping
    ([NIST PSTAR](https://physics.nist.gov/PhysRefData/Star/Text/PSTAR.html)).
  * Below the Fermi velocity this kind of model drags too little: measured
    friction in metals is about 1.6 times higher
    ([Montanari & Miraglia](https://arxiv.org/abs/1705.01575)). Squeezed fuel
    is more weakly coupled, so the gap should shrink, but every result is
    also run with the drag raised 1.6-fold.
* **Ion drag:** Coulomb collisions with the fuel's nuclei, which Pauli
  blocking does not touch.
* **NIF may already have seen the effect.** Knock-on ions crossing its cold,
  partly degenerate fuel showed a degeneracy-induced shift of their Bragg
  peak, at 90% confidence
  ([Hayes et al., Nat. Phys. 16, 432 (2020)](https://www.osti.gov/pages/servlets/purl/1881627)).

### Crystal ions into squeezed cold fuel

A deuteron slowing down in cold D-T (solid D-T is 0.25 g/cm³; NIF's fuel is
designed for about 1,000 g/cm³ and has reached 500–800):

| D-T density (g/cm³) | Fermi energy | Best return (deuteron energy) | Share of its energy lost to electrons |
|---|---|---|---|
| 0.25 (cold solid, section 1's model) | – | 0.006 (236 keV) | ~100% |
| 25 | 120 eV | 0.064 (160 keV) | 91% |
| 250 | 560 eV | 0.39 (167 keV) | 75% |
| 1,000 | 1.4 keV | **1.0** (175 keV) | 56% |
| 2,500 | 2.6 keV | 1.7 (191 keV) | 43% |
| 25,000 | 12 keV | 3.8 (236 keV) | 15% |
| 250,000 | 53 keV | 4.8 (269 keV) | 3% |

* **The cold-catcher wall breaks.** At the density NIF's fuel is designed
  for, a deuteron returns its own energy, 250 times the 0.4% limit of
  section 1. That is the first time in this project for an ion in unheated
  fuel. With the drag raised 1.6-fold it returns 0.73, and the crossing moves
  from ~940 to ~1,600 g/cm³.
* **It saturates near 5.** In the densest fuel the deuteron loses most of
  its energy to nuclei instead of electrons.
* **Protons in squeezed boron stay below 1** (0.74 at 234,000 g/cm³).
  Boron's charge of 5 makes collisions with the nuclei the drag.
* **Why it works.** Degenerate electrons drag on a slow ion about as much as
  a plasma at their Fermi temperature does (section 2's hot catcher). The
  squeeze gets that at a fifth of the energy, and without radiation.

**What the squeeze costs.** The trick lasts only while the fuel stays cold,
so the fair test gives the beam as much energy as the squeeze. Real drivers
here means 15% of the light into ions (section 7), and 1% of the light into
the squeezed fuel, as in NIF-class implosions:

| D-T density (g/cm³) | Return, perfect drivers | Return on the light, real drivers | Hot spark that ignites the same fuel |
|---|---|---|---|
| 250 | 0.2 | 0.004 | 3.3 kJ |
| 1,000 | 0.52 | 0.010 | 210 J |
| 2,500 | 0.87 | 0.016 | 33 J |
| 25,000 | 1.9 | 0.035 | 0.34 J |

* **The squeeze eats the gain.** The light comes back at 0.4–3.5%, and the
  wall plug at a few parts per thousand at best.
* **Fuel that dense ignites from a small spark.** A hot spot of joules to
  kilojoules ignites it (section 11's cheapest hot ball). The ions are far
  better spent as that spark (fast ignition, section 7). The real obstacle is
  the squeeze: up to 25 times NIF's design density, without heating.

### Hot ions in degenerate electrons (Son & Fisch's proton–boron)

**Idea**
([Son & Fisch, Phys. Lett. A 329, 76 (2004)](https://bp-pub.pppl.gov/pub_report/2004/PPPL-4003.pdf)).
Proton–boron fails because its electrons drain the hot ions' energy and
radiate it (section 4). Degenerate electrons drain at a fixed rate whatever
the density, and radiate less, while fusion grows with density. So above
some density the ions' own fusion products should keep them hot. Son &
Fisch found such a regime above ~10⁵ g/cm³, with ions near 100 keV and
electrons near 30 keV.

**Test.** A steady power balance:
* the ions lose energy to the electrons at the degenerate (Lindhard–Winther)
  rate;
* they get back the alphas' share that goes to ions;
* the electrons must radiate what they receive. Pauli blocking cuts their
  bremsstrahlung by about 2.2 (T/T_F)^1.5 (Son & Fisch give the scaling; the
  prefactor is derived here), which fixes the electron temperature.

| Ion temperature | Lowest density where the ions keep themselves hot | × NIF's design density | Electron temperature (÷ Fermi energy) | Fusion in one burn time ÷ energy to squeeze and heat |
|---|---|---|---|---|
| 50 keV | 3.9×10⁶ g/cm³ | 3,900 | 27 keV (0.08) | 3 |
| 75 keV | 8×10⁵ g/cm³ | 800 | 29 keV (0.20) | 5.1 |
| 100 keV | 3.5×10⁵ g/cm³ | 350 | 31 keV (0.35) | 5.5 |
| 150 keV | 1.1×10⁵ g/cm³ | 110 | 40 keV (0.94) | 4.9 |

* **This reproduces Son & Fisch.** With 100 keV ions, 3.5×10⁵ g/cm³ and
  electrons at 31 keV, against their > 10⁵ g/cm³ and ~30 keV.
* **The window is real, but at white-dwarf density.** It needs 110–3,900
  times NIF's design density, and a 50–160-fold radial squeeze of a solid
  without heating it. That is the same verdict as pycnonuclear D-T
  (section 10).
* **Son & Fisch said so themselves.** They found a gain of ~18 on the
  electrons' Fermi energy (40–200 in a 2006 follow-up). A power plant needs
  far more, and they listed the squeeze as unsolved. A 2025 analysis finds
  that proton–boron inertial fusion needs about 100 times today's
  compressed areal density
  ([Ochs et al., arXiv:2511.10885](https://arxiv.org/abs/2511.10885)).
* **Hotter ions don't need degeneracy.** With 200–300 keV ions the balance
  closes with electrons at 77–120 keV, 5–17 times their Fermi energy:
  ordinary hot electrons. That is the narrow classical window of Liu et al.
  (below), which needs confinement, not squeezing.

### Muon catalysis: freeing stuck muons

**Idea.** Muon catalysis came closest of the cold routes (section 9), and
sticking caps it:
* in about 0.9% of fusions the muon sticks to the alpha;
* the stuck ion leaves with 3.46 MeV, and collisions with fuel nuclei knock
  the muon free while it slows down, about 35% of the time at liquid
  density;
* electrons slow the stuck ion without freeing the muon, because an electron
  moving with it carries only ~0.5 keV against the muon's 11 keV binding.

With the electrons' drag switched off, the stuck ion meets many more nuclei
before it stops.

**Model.**
* Stripping cross section: proton-on-hydrogen ionisation (Rudd et al. 1985),
  scaled to the muonic ion.
* Slowing: the stopping above.
* One factor, tuned so that liquid D-T reproduces 35%.

| D-T density (g/cm³) | Stuck muons freed | Effective sticking | Fusions per muon allowed by sticking |
|---|---|---|---|
| liquid, 0.18 | 35% | 0.58% | 171 |
| 25 | 68% | 0.28% | 351 |
| 50 | 85% | 0.14% | 718 |
| 100 | 96% | 0.035% | 2,900 |
| 250 | 99.9% | 0.0009% | 1.1×10⁵ |

* **Squeezing frees the stuck muons.** Above ~60 g/cm³ sticking falls below
  0.1%, enough on paper for the ~1,000 fusions per muon that electricity
  needs (section 9). Most of the stripping happens while the ion is faster
  than the muon's orbit, where the cross section is most reliable. This
  matches an earlier finding: Cripps, Harms & Goel (Hyperfine Interact. 77,
  181 (1993)) found useful fusion numbers only below ~5 eV and above 100
  times liquid density.
* **But a squeeze is far too short.**
  * A 50 µm core at 100 g/cm³ holds together for 140 ps.
  * Even at the d-t-μ molecule's own fusion rate (10¹² per second), a muon
    could catalyse only ~140 fusions in that time, fewer than in liquid
    fuel.
  * The cycle's other steps have never been measured in pressure-ionised
    fuel, where the molecular resonance that drives them is gone.
* **A squeeze long enough for 1,000 fusions (1 ns) is out of reach.**
  * It needs a 0.35 mm core: 18 mg of D-T and 130 kJ even if squeezed
    perfectly cold.
  * Repaying that takes 4.5×10¹³ muons in the core at once. They live
    2.2 µs, so the source would have to make 2×10¹⁹ per second, about 10¹¹
    times the most intense muon beams.

**Verdict.** Pauli blocking is the first lever in this project that weakens
the electrons' drag itself instead of going around it, and the lever is
large: 250-fold for a fast ion, and sticking cut tenfold or more for muons. But
every use needs cold fuel squeezed to NIF's design density or far beyond,
and fuel that dense is better used the ordinary way: a small hot spark
ignites it.

## 18. Trap the X-rays: an opaque proton–boron ball

Sections 18–20 each attack a wall no earlier section touched. Code:
`nanofusion/improved/walls.py` and `tools/three_walls.py`. Numbers:
[results/walls/summary.md](results/walls/summary.md).

**Idea.** Proton–boron makes no neutrons, but its electrons radiate X-rays
faster than fusion heats the fuel (section 4). A big or dense enough ball is
opaque to its own X-rays: they are reabsorbed inside and only leak out slowly
through the surface. Does that save proton–boron?

**Model.**
* A free ball of proton–boron (15% boron) at a fixed density and ion
  temperature.
* The electrons settle where what they receive (the hot ions' drain and
  part of the alpha energy) equals the X-rays that escape.
* Escape is followed frequency by frequency: free-free absorption, Thomson
  scattering, and the escape probability from a sphere. Low-frequency X-rays
  are trapped first.
* The trapped frequencies fill with blackbody radiation, and photons near
  m_e c² make electron–positron pairs. The ball must pay for both.
* The ball "repays its heat" if its ions keep themselves hot and it
  releases at least its own energy before flying apart.

**The smallest balls that repay their heat** (all with 200 keV ions):

| Density (g/cm³) | Radius | Mass | Energy it holds | Areal density ρR | X-rays trapped |
|---|---|---|---|---|---|
| 100 | 2.5 mm | 6.5 g | 180 GJ | 25 g/cm² | 10% |
| 1,000 | 190 µm | 30 mg | 720 MJ | 19 g/cm² | 22% |
| 10⁴ | 19 µm | 0.29 mg | 6.5 MJ | 19 g/cm² | 43% |
| 10⁵ | 1.7 µm | 1.9 µg | 41 kJ | 17 g/cm² | 65% |

**The smallest balls that trap at least half their X-rays** (300 keV ions):

| Density (g/cm³) | Radius | Mass | Energy it holds | Share of it in X-rays and pairs | Fusion before it flies apart ÷ its energy |
|---|---|---|---|---|---|
| 100 | 26 cm | 7.4 t | 1.8×10²¹ J (440,000 megatons) | 100% | 3×10⁻⁴ |
| 1,000 | 7.5 mm | 1.8 kg | 3.6×10¹⁶ J (8.5 megatons) | 99.8% | 0.003 |
| 10⁴ | 46 µm | 4.2 mg | 190 MJ | 36% | 1.5 |
| 10⁵ | 1 µm | 0.42 µg | 12 kJ | 3% | 0.6 |

* **The balls that repay their heat are hot-ion balls at 17–25 g/cm².**
  * Their ions are at 200 keV and their electrons near 90 keV.
  * They trap 10–65% of their X-rays, more the denser they are.
  * That is 40–60 times the areal density of the hot spot in NIF's best
    shots (0.44 g/cm²), at 26 to 23,000 times its ~600 Gbar.
* **This matches the published studies.**
  * Lavell et al. find proton–boron burn propagating at 7,000–9,000 g/cm³
    with 14–18 g/cm², and failing at 1,000 g/cm³ with 2 g/cm²
    ([Fusion Sci. Technol. 82 (2026)](https://www.osti.gov/pages/servlets/purl/3366104)).
  * Ochs et al. find breakeven needs areal densities about 100 times the
    state of the art, at about 1,000 times the pressure
    ([arXiv:2511.10885](https://arxiv.org/abs/2511.10885)).
* **Trapping most X-rays costs more than it saves, except in
  white-dwarf-like matter.**
  * At NIF's design density, trapping half the X-rays takes 1.8 kg of fuel
    holding about 8 megatons.
  * The trapped X-rays and pairs hold 99.8% of that energy, and the ball
    returns 0.3% of it.
  * Only at 10⁴–10⁵ g/cm³ is a ball small enough to be opaque without
    holding much radiation.

## 19. Hold the fuel with a current: a fibre Z-pinch

**Idea.** Section 12 held the fuel with a gold tamper and a field that
something else had to make. A current does both at once. Run megaamperes
along a fibre of solid fuel, and the current's own magnetic field squeezes
the fibre and holds it (a Z-pinch). No laser is needed.

**Model.**
* Bennett's equilibrium, μ₀I²/8π = (N_i + N_e)kT: each current holds a given
  amount of fuel per metre at 10 keV.
* The pinch lives some number of Alfvén times (radius ÷ Alfvén speed)
  before instabilities tear it.
* A narrower pinch is denser and burns faster, so there is a radius at which
  fusion repays the plasma's heat.
* Above 1.35 MA the pinch's own field holds the fusion alphas inside it.

**Checks.** Bennett's relation puts a 100 µm deuterium fibre at 640 kA at
140 eV; fibre pinches of that kind reached a few hundred eV
([Sethian et al. 1987](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.59.892)).
At 300 kA and 1 keV it holds 1.4×10¹⁹ ions per metre; Zap Energy's FuZE
pinch holds about 1.1×10¹⁹
([Levitt et al. 2024](https://www.osti.gov/pages/servlets/purl/2337575)).

| Current | Solid D-T fibre it holds at 10 keV | Squeeze to repay the heat in one Alfvén time (= Alfvén times needed unsqueezed) | That pinch lives | Squeeze for the alphas alone to repay it | Field energy ÷ plasma heat |
|---|---|---|---|---|---|
| 0.64 MA | 13 µm across | ×240 (1.1×10⁴ g/cm³) | 22 fs | ×1,200 | 17 |
| 1.4 MA | 29 µm across | ×110 (2,300 g/cm³) | 110 fs | ×560 | 15 |
| 5 MA | 110 µm across | ×31 (180 g/cm³) | 1.4 ps | ×160 | 12 |
| 26 MA (Sandia's Z) | 550 µm across | ×6 (6.8 g/cm³) | 37 ps | ×30 | 7.5 |
| 60 MA | 1.3 mm across | ×2.6 (1.3 g/cm³) | 200 ps | ×13 | 5.3 |

* **At the Z machine's current the numbers are not crazy.** A D-T fibre
  half a millimetre across repays its heat if squeezed sixfold and held for
  one Alfvén time, or held unsqueezed for six.
* **But fibre pinches haven't lasted.**
  * Sethian's deuterium fibres stayed stable for ~100 growth times, but only
    while the current rose and only at a few hundred eV. At peak current
    they went unstable.
  * Other fibre pinches (Los Alamos, Imperial College) expanded and broke up
    early, and their neutrons came mostly from beams, not heat.
  * Zap Energy's sheared-flow pinch lasts 1,000–10,000 growth times, but at
    a millionth of solid density and 1–3 keV.
* **The field costs more than the plasma.** At 26 MA the field around the
  pinch holds 7.5 times the plasma's heat, 3.8 MJ per cm, out of the ~20 MJ
  that Z stores. Even ignition-like conditions are far from net electricity
  unless the burn multiplies the energy many times over.
* **Lithium deuteride fibres are hopeless.** D + ⁶Li burns thousands of times
  more slowly than D-T: a 1 MA LiD fibre would need a 160,000-fold squeeze.

## 20. Steer through the lattice: channelling in a tritide crystal

**Idea.** A fast ion entering a crystal almost along a row or plane of atoms
is steered gently between them (channelling). In the open channel it meets
fewer electrons, so it loses energy more slowly, and its flux piles up in the
channel centre. In titanium tritide the tritons sit in the channel centres
along ⟨100⟩. Can the crystal itself beat the cold-catcher limit?

| Deuteron energy | Axial critical angle (Lindhard) | Share of an isotropic explosion along an axis | ... along a plane (generous) |
|---|---|---|---|
| 100 keV | 2.6° | 1.3% | 23% |
| 300 keV | 1.5° | 0.44% | 13% |
| 1 MeV | 0.81° | 0.13% | 7% |

A 300 keV deuteron stopping in non-channelling titanium tritide returns 0.0038
of its energy. With channelling:

| Case | Flux at the tritons | Stopping ÷ random | Return |
|---|---|---|---|
| Perfectly aligned beam, generous | ×3 | 0.5 | 0.023 |
| Perfectly aligned beam, measured values | ×2 | 0.9 | 0.0084 |
| Exploding crystal, generous, half the path channelled | axes ×3, planes ×2 | 0.5, 0.9 | 0.0041 |

The measured values: flux peaking of 1.4–2.2 at interstitial sites, and
channelled stopping only 5–13% below random at the speeds of 100–300 keV
deuterons. Stopping halves only above ~1 MeV per nucleon.

* **A perfectly aligned beam does 2–6 times better** (0.008–0.023), still
  40–120 times short of breakeven. Experiments agree on the size: 20 keV
  deuterons on palladium gave about twice the neutrons along a crystal plane
  as at 30° to it
  ([Dalkarov et al. 2019](https://journals.aps.org/prab/abstract/10.1103/PhysRevAccelBeams.22.034201)).
* **An exploding crystal can't aim.** Under 0.5% of its ions fall inside an
  axial critical angle and ~13% along a plane, and they dechannel partway: a
  9% gain even with generous numbers.
* **Blinding the electrons with a magnetic field instead** needs their
  cyclotron frequency to beat the plasma frequency: about 68,000 T in solid
  D-T, 85 times the strongest laser-made field claimed.

## 21. Every route in one simulator, and the ideas combined

**Idea.** Put every experiment above into one energy account, add ordinary
inertial fusion as a yardstick, combine the ideas, and let a search try
thousands of designs for the best electricity out per electricity in. Code:
`nanofusion/improved/simulator.py` (`python -m nanofusion sim`), guide:
[simulator.md](simulator.md), numbers:
[results/simulator/summary.md](results/simulator/summary.md).

**The account.**
* electricity → driver (laser, accelerator or pulsed current) → fuel →
  fusion → heat (with a lithium blanket) → electricity;
* realistic-to-generous efficiencies: 10% lasers, 15% pulsed power, 40%
  accelerators, 40% heat to electricity.

**Stretch.** Each design's demands (laser energy and power, field, density,
pressure, current, pinch stability, muon supply, ...) are compared with the
best ever achieved. The largest ratio is its stretch. Designs are ranked
within stretch budgets. Ideas never shown to work are flagged and left out
(section 22).

**The yardstick.** The `compressed_icf` route is ordinary inertial fusion: an
isobaric hot spot (10 keV, ρR 0.3 g/cm²) in a dense shell, burning the share
ρR/(ρR + 7 g/cm²). NIF-like inputs give a target gain of 1.7, against NIF's
record 4.1, so it is if anything conservative.

**New routes and combinations** beyond sections 1–20:
* a MagLIF-like magnetised liner (section 12's fuel, squeezed by pulsed power);
* exploding crystals as the muon factory for muon catalysis;
* inertial fusion lit by a crystal spark instead of a hot spot, or driven by
  pulsed-power X-rays instead of a laser;
* the hollow shell, NASA's lattice confinement, and a proton–boron catcher
  with screening and the avalanche.

| Best electricity out / in | Today's technology | ≤ 10× today | ≤ 100× today | Unlimited |
|---|---|---|---|---|
| Ordinary inertial fusion (yardstick) | – | **1.5** | **3.0** | **7.1** |
| Dense D-T Z-pinch | 0.01 | 0.89 | **2.8** | **5.7** |
| Dense D-T Z-pinch followed in time, best with a liner (sections 24–25) | 0.003 | **22** | **22** | **22** |
| Muon catalysis | 0.12 | 0.16 | 0.16 | 0.16 |
| Best route using exploding crystals | 0.005 | 0.005 | 0.012 | 1.3 (tamped ball) |

Untested theory is left out of this table; section 22 adds it.

* **Nothing reaches breakeven with today's technology.** The best is muon
  catalysis, which returns about a tenth of the electricity it costs.
* **Ordinary inertial fusion is the most proven route beyond that.** Within
  ten times today's best on every demand it reaches 1.5. That needs a 19 MJ
  laser, a cold shell of 6 g/cm² (about ten times NIF's) and a low adiabat, as
  power-plant studies of inertial fusion also find.
* **On paper the time-dependent pinch leads** (section 25): 22 within ten
  times today.
  * The design is a thin beryllium liner lined with D-T ice, at 34 MA.
  * But it converges ~170×, five times NIF's best, in a 0-D model with no
    liner instability. Its convergence is the number to doubt.
* **The static Z-pinch** gives 0.89, if it stays stable for many Alfvén
  times, which no dense, hot pinch has done.
* **Exploding crystals never help.**
  * As the spark for inertial fusion they do worse than an ordinary hot spot
    at every budget (0.64 against 7.1 even with no limits).
  * As a muon factory they make muons, but at 0.006 of breakeven.
  * Their best route, the tamped ball, stays at about 1% of breakeven within
    a hundred times today's technology.
* **Lasers beat pulsed-power X-rays** for driving a capsule (1.5 against 0.61
  within ten times today). On Z only ~0.2% of the stored energy reached a
  capsule.
* **Small ideas stay small.**
  * Electron screening and the proton–boron avalanche change nothing visible.
  * Heating the ions first helps the tamped ball two to five times.
  * Deuterium-only fuel costs a factor of 30–100.
* **The search exposed four flaws in the models**, now fixed:
  * a pinch "burning" more than all its fuel;
  * a pinch at 130,000 Gbar;
  * a pinch too cool for its alphas to make up its X-ray losses;
  * a capsule whose hot spot was 600 times the pressure of the shell around
    it.
* **A uranium blanket** (k_eff 0.95) makes several routes "net energy",
  muon catalysis included, but over 99% of the heat is then fission.

## 22. Untested theory, as if it worked

**Idea.** Some proposals have never been shown to work at all. The simulator
keeps them apart from technology: a design that relies on one is flagged,
left out of section 21's search, and searched separately as if the theory
held (`python -m nanofusion sim search --theory`). Numbers:
[results/simulator/summary.md](results/simulator/summary.md).

| Untested idea | What it promises | Model and source |
|---|---|---|
| Alpha channelling | waves hand the alphas' energy to the fuel ions instead of the electrons, the key hope for proton–boron | that share (up to 90%) of the alphas' electron heating goes to ions ([Fisch & Rax 1992](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.69.612); "more than half" claimed in later work) |
| Direct conversion | charged fusion products become electricity directly | 60%, as designs claim (60–70%); ~48–65% has been shown only on low-energy ions and beams |
| Freeing stuck muons with fields | stuck (μα)⁺ ions re-accelerated until stripped; sticking 0.1% | Kulsrud; the TEK concept. Both papers find the fuel must stay below 25 K, which costs ~30 MeV of cooling per 17.6 MeV fusion; charged here |
| Antimatter sparks | antiprotons annihilate on uranium in compressed D-T; fission fragments light the hot spot | ~0.17 GeV deposited per antiproton; CERN's ~2 mJ of beam per antiproton, or a factory 10,000 times better (ICAN-II, AIMStar) |
| Light-assisted tunnelling | intense light raises the D-T tunnelling probability in cool fuel | 5×10²¹ W/cm² held over the fuel; boost up to 10⁹ at 1 keV ([Queisser & Schützhold 2019](https://arxiv.org/abs/1902.04905); Wang 2020), or +70% (Thomson et al. 2024) |
| Spin-polarised fuel | D-T cross section ×1.5 | never tested in a plasma |

| Best electricity out / in, theory allowed | Today's technology | ≤ 10× today | ≤ 100× today | Unlimited |
|---|---|---|---|---|
| Ordinary inertial fusion (with spin-polarised fuel) | – | **1.8** | **3.5** | **8.2** |
| Dense D-T Z-pinch (with spin-polarised fuel) | 0.018 | **1.2** | **3.4** | **7.1** |
| Dense D-T Z-pinch followed in time (liner ice burnt as an ICF shell) | 0.22 | **37** | **37** | **37** |
| Antimatter spark (theoretical factory) | – | 0.16 | 0.39 | 0.39 |
| Light-assisted tunnelling | – | 0.032 | 0.091 | 0.85 |
| Muon catalysis (stuck muons freed) | 0.14 | 0.2 | 0.2 | 0.2 |
| Proton–boron (channelling, direct conversion) | – | 0.009 | 0.015 | 0.023 |

* **Theory moves the numbers, not the verdict.**
  * Spin-polarised fuel lifts inertial fusion and the pinch by 15–35%.
  * If a liner's ice burns like an ICF shell, the time-dependent pinch goes
    from 22 to 37 within ten times today.
  * Every other untested idea stays far below breakeven.
* **Proton–boron gets about 50% better and stays at 2% of breakeven.**
  Channelling feeds the ions, but the X-rays still win. Direct conversion
  helps more than channelling here.
* **Antimatter sparks need more antiprotons than have ever existed at once.**
  A NIF-like target needs ~10¹⁵ antiprotons, and the best design 1.3×10¹⁴,
  22 times the most ever stored. Even from a factory 10,000 times better than
  CERN's, making them costs more than the fusion they spark. Livermore's
  independent estimate also needed at least ~3×10¹⁵ antiprotons.
* **Light-assisted tunnelling costs more light than it returns.** Holding
  5×10²¹ W/cm² over the fuel for its whole burn costs gigajoules. The optimistic
  boost does help: 0.091 against 0.037 with the pessimistic one for the same
  ball. Both are far below breakeven.
* **Freeing stuck muons doubles the fusions per muon**, but the cooling it
  needs leaves 0.14 against 0.12. A 2026 analysis of X-ray stripping finds at
  most +40% ([arXiv:2606.07077](https://arxiv.org/abs/2606.07077)).

## 23. A practical experiment: what limits NIF's yield

**Idea.** Ordinary inertial fusion leads every search above, so design a real
experiment with it: a campaign of shots on NIF, one change at a time. Then
check whether it is novel, and simulate what the simple model leaves out.
Proposal: [icf_experiment.md](icf_experiment.md). Code: `tools/icf_experiment.py`
(static model) and `tools/icf_dynamics.py` (1-D model). Numbers:
[summary.md](results/icf_experiment/summary.md) and
[dynamics.md](results/icf_experiment/dynamics.md).

**Two models.**
* **Static:** the simulator's `compressed_icf`, an isobaric hot spot in a
  cold shell, calibrated on NIF's fuel and measured pressures.
* **1-D:** `nanofusion/improved/burn1d.py` follows the shell from peak
  velocity through stagnation, alpha heating, burn and disassembly. It uses
  NIF's published velocities, leftover ablator and gas fill. It matches
  LLNL's own 1-D yields: 5.1 MJ for N210808 against their 4.2–6.2 MJ.

**Is it novel? Mostly not.**
* Measured ÷ clean yield is NIF's standard yield-over-clean.
* NIF has scanned every knob.
* Landen et al. (2024) proposed spending extra laser energy on fuel instead of
  ablator.
* **The static model's central claim is falsified.** It said a lower adiabat
  alone does nothing (×1.01). Hohenberger et al. (2025) did exactly that and
  got up to 80% more yield. The 1-D model gets the direction right: at fixed
  drive, a lower adiabat raises pressure and areal density.

**What the static model leaves out**, simulated on N221204's design:

| Physics | Effect |
|---|---|
| Alpha heating | ×530 |
| Burn propagation | 27% burnt, against 21% from the static formula |
| Velocity at fixed kinetic energy | pressure ∝ v^2.4; below ~380 km/s, no ignition |
| Mixed-in ablator | 100 ng: −2% on a clean capsule; 1000 ng: −74% |
| Asymmetry (unstagnated kinetic energy) | 3%: −12%; 10%: −98% |

* **The static ceiling misled.** The clean 1-D yield of N221204's design is
  20 MJ, not 10.9. The record reached 43% of it, not 79%.

**The campaign, re-predicted.** Each step's clean 1-D yield is taken at the
record's 43%, with the record's degradation held fixed as an upper bound:

| Step | Predicted yield |
|---|---|
| Repeat the record | 8.6 MJ |
| Lower adiabat (2.45 → 2.25) | 9.5–18.5 MJ |
| 2.2 MJ spent on more ablator (NIF's choice) | 9.6–20 MJ |
| 2.2 MJ spent on more fuel (Landen's proposal) | 10.4–21.9 MJ |
| 2.6 MJ upgrade, scaled capsule | 12.8–27.5 MJ |

* **The open test is fuel against ablator at 2.2 MJ.** It hasn't been shot,
  as far as the search found. The 1-D model gives fuel an 8% edge.
* **No net electricity.** The plant designs of section 21 need about ten times
  NIF's laser energy and areal density.

## 24. A practical dense Z-pinch experiment

**Idea.** The dense D-T Z-pinch came second in section 21's search. Follow
it in time, design a first experiment, and check whether it is novel and
possible. Proposal: [zpinch_experiment.md](zpinch_experiment.md). Code:
`nanofusion/improved/zpinch.py`, `tools/zpinch_experiment.py`. Numbers:
[results/zpinch_experiment/summary.md](results/zpinch_experiment/summary.md).

**The model.**
* A generator circuit drives a fibre that heats ohmically and by compression.
* It radiates with its own opacity, and loses heat to the electrodes.
* It burns, with alphas held by its self-field.
* The m = 0 mode grows in Alfvén times.
* It reproduces the Pease–Braginskii current (1.3–1.6 MA, against Haines'
  1.15–1.37) and FuZE's Bennett temperature.

**Three barriers the static `fibre_pinch` route missed**, on Z (26 MA, a
240 µm D-T fibre):
1. **A cold fibre can't take the current.** Its resistance stalls the current
   at ~4 MA, which is why the 1980s fibres grew a current-carrying corona.
   Laser preheat to ~100 eV fixes this.
2. **The radiative collapse caps its own current** at ~8.5 MA. Most of the
   energy leaves as X-rays.
3. **Stability decides everything:**
   * broken in an Alfvén time: nothing;
   * stable only while the current rises: under 1 J;
   * fully flow-stabilised: 1.3 MJ from a 7.8 MJ bank, with no ignition.

**Novel? Possible?**
* No one has fielded D-T in a dense equilibrium pinch.
* Fibre pinches stopped around 2000.
* Flow-stabilised pinches have reached only 0.6 MA, in deuterium.
* Tritium on Z has been a trace once (2016).
* So **D-T is novel but not possible today**. Preheated deuterium fibres on Z
  with a 1% tritium trace are. The model's three scenarios predict anything
  from no thermonuclear neutrons to ~10⁹ to ~10¹⁵. Neutron isotropy, the
  D-T/D-D ratio and the collapse's back-voltage tell them apart.

**Energy.**
* On Z, even fully flow-stabilised, the best is 1% of breakeven.
* On a 60 MA machine (studied, not built), a ~1 mm fibre ignites and
  returns 1.9× its electricity, but only if it stays stable through its whole
  burn.
* As the simulator route `dynamic_pinch`, a bare fibre reaches about 2
  within ten times today, against the static route's 0.89. With a liner
  (section 25) the route reaches 22.

## 25. Making the dense pinch testable, and what energy takes

**Idea.** The fibre pinch of section 24 stalls on Z. Find a version Z can
test now, and what it would take to make energy. Code:
`tools/zpinch_energy.py`; numbers:
[results/zpinch_experiment/energy.md](results/zpinch_experiment/energy.md);
plan: [zpinch_experiment.md](zpinch_experiment.md).

**Why the fibre stalls.** A dense, cold column is too resistive to take the
current. Keeping the current rising, a longer fibre, or a stiffer driver
doesn't help: the column cools, its resistance climbs, and the current
falls.

**A gas column instead.** A centimetre-wide column of gas, as in Z's gas
puffs, conducts easily and heats by its own implosion.
* The model puts Z's deuterium gas-puff record (3.9×10¹³ neutrons) between
  one and two Alfvén times of stability.
* It converges too far: 0.1–0.2 mm against 2.5 mm measured.

| D-T gas column | Fusion per cm | Electricity out ÷ in |
|---|---|---|
| Z, 26 MA, held 3–30 Alfvén times | 3–9 MJ | 0.04–0.12 |
| 60 MA, held 3 Alfvén times | 341 MJ | 0.85 |
| 80 MA, held 10 Alfvén times | 2.4 GJ | 3.3 |

* **The energy designs need deep convergence**, 160–260×, where the best
  implosion of any kind converges about 35×. They also need a few Alfvén
  times of stability and fifty-fifty D-T: 20% tritium costs a factor of ~8.
* **A more efficient driver** scales gain directly: at the 24% of the Z 300
  design, 60 MA gives 1.4.
* **Possible now:** deuterium gas columns on Z with a 1% tritium trace. They
  measure how many Alfvén times the pinch holds and how far it converges,
  the two numbers energy turns on.
* **The literature points the same way:**
  * the only driver above 60 MA being built is Pacific Fusion's (target
    2030);
  * sheared-flow stabilisation is doubtful at these densities;
  * a fuel column inside a MagLIF liner is the likelier winner. Pacific
    Fusion's own designs give 109–816 MJ at 56–65 MA, the same class as the
    gas column here.

**With a liner.** The model now has one: a MagLIF-like beryllium tube
imploding onto preheated, magnetised gas, optionally lined with frozen D-T
ice.
* The field's leak is calibrated to MagLIF's measured field × radius.
* As a 0-D model it overshoots Z's MagLIF yields 9–12× (too hot, burn too
  brief) and undershoots published high-current designs. So read its
  numbers as a comparison between designs.

| Liner, 50/50 D-T | Electricity out ÷ in | Convergence | Compressed field |
|---|---|---|---|
| 60 MA, 600 µm of ice | 3.5 | 131 | 63 kT |
| 80 MA, 600 µm of ice, 10 T | 4.2 | 43 | 3.7 kT |

* **Thick ice burns at 60–80 MA**, by the model's own ablation; below ~50 MA
  the field keeps the heat from it.
* **At 80 MA the liner needs no more convergence than MagLIF already shows**,
  and no stabilisation of the fuel. It is the most credible energy design
  here, on a driver that doesn't exist yet.
* **The simulator's search goes further on paper**: 22 at 34 MA, with a thinner
  liner and ~170× convergence. That is where a 0-D model with no liner
  instability is least trustworthy.

## 26. An experiment for the thin, ice-lined liner

**Idea.** Test the simulator's best design within ten times today: a 4.5 mm
beryllium liner with a wall 1/12 of its radius, lined with 580 µm of D-T ice,
4.3 T, at 34 MA. In the model it gives 2.9 GJ, 22× the electricity in.
Proposal: [liner_experiment.md](liner_experiment.md); code:
`tools/liner_experiment.py`; numbers:
[results/liner_experiment/summary.md](results/liner_experiment/summary.md).

**Has it been done? No, and what's new is what's least credible.**
* The closest published design is Sefkow et al.'s 70 MA ice burner (2014):
  the same geometry, twice the current, converging only 24–35×.
* This design ignites at 34 MA, below every published threshold (42–56 MA).
* It converges ~170×, where designs stay below 40× and the best radiographed
  liner reached 20 ± 3.
* No ice-lined liner has been imploded on Z.

**The model's own verdict.**
* **All or nothing:** if the liner breaks up at MagLIF's ~40×, the shot makes
  25 J instead of 2.9 GJ.
* **A knife-edge:** 30 MA, a 1/9 wall, or MagLIF's usual 1 mg/cm³ of gas each
  cut the yield 8,000–110,000-fold.

**The experiment.**
* **Stage 1, on Z now, in deuterium (~$20–50M):** image the stagnated
  column. The model predicts ~10 µm if the liner holds, ~90 µm if it breaks
  like MagLIF's. Anything near 0.1 mm ends the design.
* **Stage 2:** D-T on a ~34 MA driver, which doesn't exist yet. Pacific
  Fusion's >60 MA machine, run below full current, could serve.

## 27. The thin liner in 2-D and 3-D MHD

**Idea.** Section 26's design needs its liner intact to ~170× convergence, and
the 0-D model has no liner instability. So implode it in Athena++, Princeton's
open MHD code, with a ripple on the liner to seed the magneto-Rayleigh-Taylor
instability. Code: `tools/athena_liner.py` and `tools/athena/`. Numbers:
[results/liner_athena/summary.md](results/liner_athena/summary.md).

**Result: the liner breaks up (2-D).**
* With no ripple, the fuel converges 69×; the 0-D model says 177×.
* A 0.1 µm ripple, about a diamond-turned finish, cuts that to 25×. A 1 µm
  ripple cuts it to 9×.
* Beryllium spikes reach the gas 15–20 ns before stagnation. They mix about a
  fifth of the gas's mass of beryllium into it.
* At 25×, the 0-D model's 2.9 GJ falls to 0.5 J.

**3-D (a quick look, half resolution) is no better.**
* The same 1 µm ripple, seeded also with patterns running around the liner,
  left 7× convergence at 0.1 keV.
* Those patterns grew fastest, because the axial field can't stiffen them.
  They pinched the gas off the axis, and by 135 ns beryllium filled it.

**Caveats.**
* Ideal MHD.
* The beryllium is a soft ideal gas, with no strength or resistivity.
* A prescribed current.
* 0.5 mm of liner.

These could exaggerate the instability by a factor of a few, not by the 7–20×
in convergence the design is missing. Stage 1 of section 26 on Z would most
likely see the column stop near 0.1 mm, ending the design.

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

**What moved the wall, and what didn't (sections 12–16).** Only confinement
did. A tamper and a magnetic field hold the fuel about 15 times longer, and
that cuts the breakeven energy a thousandfold. It is a Lawson-type fix: the
fuel is hot and held, so its ions thermalise instead of being lost. The
other four routes change the driver (an accelerator instead of a laser),
the product (neutrons into a blanket), or the plasma (no electrons). Each
meets a limit already found: electrons stop fast ions, and like charges keep
bare ions apart.

**Pauli blocking dents the wall itself (section 17).** In fuel squeezed until
its electrons are degenerate, the drag on a slow ion stops growing with
density. A fast ion then returns its own energy at NIF's design density, and stuck
muons are knocked free. But it takes cold fuel at or far beyond NIF's
density, and fuel that dense is better ignited the ordinary way.

**The other walls (sections 18–20).** Trapping proton–boron's X-rays makes
the trapped radiation itself the heat sink, except in white-dwarf-like matter. A current can hold fuel as well
as a tamper and field can, at the price of a field that holds several times
the plasma's heat. A crystal's channels raise the return two- to sixfold at
best, for ions that could hardly be aimed into them.

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
* **Holding the ball together** with a gold tamper and a kilotesla field is
  the biggest lever found. On paper, breakeven on the fuel's energy falls
  from 60 MJ to about 50 kJ, or 10 kJ if the crystals heat the ions first,
  and 0.1–0.7 MJ with the straight field a real coil makes. It is still not
  net energy: 0.15 of the laser energy comes back at breakeven and 1.5 at
  best, from a laser, a field and a target that do not exist. And it is
  magnetised target fusion, with the crystals only as the heater.
* **Four more routes fail outright.**
  * A fission blanket fed by crystal-made neutrons would have to be a
    critical reactor in all but name.
  * D-T supports no chain reaction through its own neutrons and alphas.
  * Bare ions in a magnetic trap are capped by their own charge at watts per
    cubic metre.
  * Crystals fired as projectiles hit the same 0.4% cold-catcher limit.
* **Pauli blocking** in squeezed, degenerate fuel weakens the electrons' drag
  itself, the first lever here that does.
  * A fast ion returns its own energy at 1,000 g/cm³, but paying for the
    squeeze leaves 0.4–3.5% of the light.
  * Son & Fisch's proton–boron keeps itself hot only at white-dwarf density.
  * Stuck muons are freed above ~60 g/cm³, but a squeeze lasts ~0.1 ns, far
    too short for the 1,000 fusions per muon that electricity needs.
* **Untested theory** (section 22) mostly moves the numbers, not the
  verdict: proton–boron stays at 2% of breakeven, antimatter sparks need more
  antiprotons than have ever existed at once, and light-assisted tunnelling
  costs more light than it returns.
* **Every route in one simulator** (section 21): with today's technology
  the best is a tenth of breakeven. Within ten times today's technology
  ordinary inertial fusion crosses net electricity (1.5). The time-dependent
  pinch with a liner goes further on paper (22), but only by converging ~170×
  in a model with no liner instability. Combining exploding crystals with any other idea stays below
  ~1% of breakeven.
* **A practical experiment** (section 23): a one-change-at-a-time shot
  campaign on NIF. It is mostly not novel, and its static model's claim that
  the adiabat doesn't matter was already falsified on NIF. A 1-D
  implosion-and-burn model matches LLNL's 1-D yields, puts the 2025 record at
  43% of a clean implosion, and predicts 13–28 MJ after the 2.6 MJ upgrade.
  The open test is spending extra laser energy on fuel rather than ablator.
* **A dense Z-pinch experiment** (section 24): followed in time, a D-T fibre
  on Z can't take the current cold, collapses radiatively when preheated, and
  reaches 1% of breakeven even if fully stabilised. D-T pinches are untried
  and, for want of tritium on Z, not possible yet. Preheated deuterium fibres
  with a tritium trace are.
* **Making the pinch testable and energy-producing** (section 25): a gas
  column instead of a fibre takes the current and heats itself. Z can test it
  now with deuterium. In the model it makes net electricity at 60–80 MA, but
  only with 160–260× convergence and a few Alfvén times of stability. A
  MagLIF-like liner with thick D-T ice does better: 4.2× the electricity in
  at 80 MA, converging no further than MagLIF already does, on a driver
  beyond the ~60 MA being built for ~2030.
* **An experiment for the search's best liner** (section 26): 34 MA, a thin
  wall and thick ice give 22× in the model. But it needs ~170× convergence
  (no liner has shown more than ~40×) and sits on a knife-edge. It hasn't
  been tried. A deuterium test on Z, imaging the stagnated column, would
  settle it for ~$20–50M.
* **Three more walls attacked, none breached.**
  * An opaque proton–boron ball traps its X-rays, but the trapped radiation
    and pairs soak up its energy. The balls that do repay their heat need
    40–60 times the areal density of NIF's hot spot, as other studies find.
  * A fibre Z-pinch holds fuel with its own current. At the Z machine's
    current the squeeze is modest, but the pinch must survive its
    instabilities, and its field holds several times the plasma's heat.
  * Channelling in a tritide crystal returns at most 2.3% of an aligned
    ion's energy (0.8% with measured values), and exploding crystals can't
    aim.

## Caveats (generous or conservative?)

* **Generous:**
  * perfect laser-to-ion conversion in tests 1–2;
  * instant, uniform stripping of a whole ball in test 3;
  * ignoring electron energy beyond 1 keV;
  * the generous fast-ion bound;
  * in section 12: a field lying along the wall everywhere, no gold mixed
    into the fuel, perfect symmetry, and 15% of the light reaching the fuel.
* **Conservative:**
  * no alpha-particle self-heating in test 3 (it matters only for balls with
    ρR ≳ 0.3 g/cm², which are larger than the breakeven sizes here);
  * D+⁶Li counted for the (d,α) channel only, maybe a factor of 2 low;
  * section 12 stops following a ball once its radius has grown fivefold,
    because the wall model is not trusted further. Balls that ignite are
    still burning then, and would return 20–50% more.
* **Section 17:**
  * generous: electrons at zero temperature; target ions at rest when sharing
    the alphas' energy; ions and electrons merging their drag smoothly near
    the Fermi temperature (up to 2 times generous there);
  * uncertain either way: linear response, which drags 1.6 times too little
    on slow protons in metals (every result is also given with that much
    more drag); the bremsstrahlung prefactor for degenerate electrons,
    derived here; and a single-state model of the stuck muon, tuned once at
    liquid density.
* **Sections 18–20:**
  * generous: the opaque ball ignores Compton scattering's energy exchange
    and the annihilation of its pairs; the pinch is uniform, at its Bennett
    temperature from the start;
    channelling uses the widest critical angle for every axis;
  * rough: the opaque ball's frequency-resolved escape (scattering folded
    into an effective depth) and the pinch's lifetime in Alfvén times are
    order-of-magnitude models.
* **Model limits:**
  * 0-D burn with fixed temperature and a fixed burn time;
  * cross sections held constant above their measured ranges;
  * non-relativistic bremsstrahlung, which underestimates losses at hundreds
    of keV;
  * section 12 is also 0-D, and its answer goes as the cube of the speed at
    which the gold wall recedes (taken from gold's shock data, good to
    perhaps 20%).

None of these change the conclusion by the orders of magnitude involved.
