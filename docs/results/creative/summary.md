# Creative routes: cold squeeze, crystal muon factory, polarised fuel

Produced by `tools/creative_routes.py`. Code: `nanofusion/improved/alternatives.py`.

## 1. Cold squeeze: can a crystal fuse by compression alone?

Pycnonuclear rate of Salpeter & Van Horn (ApJ 155, 183 (1969), static bcc lattice), with the fuel kept at zero temperature, so its pressure is that of the degenerate electrons. Ignition: the pycnonuclear burn must consume 1% of the fuel before the ball flies apart (τ = R / 4c_s). The alpha particles then heat it to thermonuclear temperatures, as in a white dwarf igniting.

| D-T density (g/cm³) | × solid | reactions per nucleus per s | pressure (Gbar) | compression in radius from solid | ball radius that ignites | its mass | ideal compression energy |
|---|---|---|---|---|---|---|---|
| 0.21 | 1 | 3e-177 | 0.00016 | 1 | more than a planet | – | – |
| 20 | 95 | 2.3e-73 | 0.32 | 4.57 | more than a planet | – | – |
| 1000 | 4.8e+03 | 4.2e-29 | 2.2e+02 | 16.8 | more than a planet | – | – |
| 10000 | 4.8e+04 | 1.6e-13 | 9.9e+03 | 36.2 | more than a planet | – | – |
| 100000 | 4.8e+05 | 0.01 | 4.5e+05 | 78.1 | 3.4e+03 km | 1.6e+28 kg | 1.1e+40 J |
| 1e+06 | 4.8e+06 | 3.6e+05 | 1.8e+07 | 168 | 0.19 m | 2.9e+07 kg | 8.5e+19 J |
| 3e+06 | 1.4e+07 | 2.1e+08 | 1e+08 | 243 | 0.00043 m | 1 kg | 5.7e+12 J |
| 1e+07 | 4.8e+07 | 7.6e+10 | 6.1e+08 | 362 | 1.6e-06 m | 1.6e-07 kg | 1.8e+06 J |
| 3e+07 | 1.4e+08 | 7.2e+12 | 2.9e+09 | 523 | 2.1e-08 m | 1.1e-12 kg | 22 J |

* **It takes about 10⁷ g/cm³**, white-dwarf density. Below that, the fuel flies apart long before it burns. At 3×10⁶ g/cm³ the ball would have to be 1 kg, with 6 TJ of compression. At 10⁷ g/cm³, 0.16 mg would ignite after 1.8 MJ of ideal, shock-free compression. At 3×10⁷ g/cm³ even a nanogram would, with 22 J. The ideal gain is about 30: fusion releases far more than squeezing the electrons costs.
* **But the conditions are out of reach.** That is ~10⁴ times the density NIF reaches, at ~6e+08 Gbar (NIF's hot spot is a few hundred Gbar). It means shrinking the fuel 362-fold in radius, against ~30–40 at the limit of hydrodynamic stability. And it must happen without heating it, when every real compression heats by shocks.
* **At NIF's ~10³ g/cm³** a cold 30 µm fuel shell would need the pycnonuclear rate to be 6e+36 times the standard prediction to ignite.
* **One dissenting theory.** Ichimaru & Kitamura (Phys. Plasmas 6, 2649 (1999)) predicted 'prodigiously high' p-d rates in liquid metallic hydrogen at ≥ 20 g/cm³ and ≥ 1 Gbar near solidification. If D-T rates there were high enough to ignite, cold fuel at 20 g/cm³ (where squeezing the electrons costs only ~60 eV each) would burn, and that would change everything. The standard rate there is ~10⁻⁷³ per nucleus per second. Nobody has compressed hydrogen cold to 1 Gbar (static and ramp experiments reach ~0.005–0.05 Gbar, 0.5–5 TPa), so this is untested, not refuted. Light nuclei sit outside the regime where the pycnonuclear formulas are reliable, and different lattice models disagree by orders of magnitude.

## 2. Exploding crystals as a muon factory

Muons are made by pions from ~1 GeV-per-nucleon ions hitting a target. The best accelerator designs need ~4.7 GeV of beam per stopped muon (Kelly et al. 2021); deuteron beams ~8 GeV (Jändel/Petrov studies). A Coulomb-exploding ⁶LiD crystal gives every ion at radius r the same energy per nucleon, rising as r² to the edge. Below ~0.3 GeV per nucleon they make almost no pions. Here muons per GeV of beam rise linearly from that threshold to the best value at 1 GeV per nucleon.

| edge energy (GeV per nucleon) | crystal radius | stripping intensity (× today's record) | beam per muon (4.7 GeV best case) | laser per muon, 15% → ions | 50% → ions | 100% → ions | muon gain at 150 fusions per muon (15% / 50% / 100%) |
|---|---|---|---|---|---|---|---|
| 0.5 | 0.89 µm | 2e+23 W/cm² (2×) | 40.6 GeV | 270 GeV | 81.1 GeV | 40.6 GeV | 0.01 / 0.03 / 0.07 |
| 1 | 1.3 µm | 4e+23 W/cm² (4×) | 7.86 GeV | 52.4 GeV | 15.7 GeV | 7.86 GeV | 0.05 / 0.17 / 0.34 |
| 2 | 1.8 µm | 8e+23 W/cm² (8×) | 5.06 GeV | 33.7 GeV | 10.1 GeV | 5.06 GeV | 0.08 / 0.26 / 0.52 |
| 4 | 2.5 µm | 2e+24 W/cm² (15×) | 4.76 GeV | 31.7 GeV | 9.52 GeV | 4.76 GeV | 0.08 / 0.28 / 0.55 |

* **Today's laser muon sources** (laser wakefields, ELI-NP) cost ~4e+08 GeV of laser light per muon. Exploding crystals would cut that ~1e+07-fold on paper.
* **But it is still no better than an accelerator.** Even if all the laser energy became ion energy, a muon costs the accelerator's ~4.7 GeV, and lasers turn electricity into light far less efficiently (≤ 20% for future diode-pumped systems) than accelerators turn it into beam.
* **Gain is 0.08–0.55 against the light**, at 150 fusions per muon. Against the wall plug that is ~0.02–0.1, before turning heat into electricity.
* **The crystals are also out of reach.** They would need to be 2.5–5 µm across, stripped at 10²³–10²⁴ W/cm². Relativistic transparency might let the light in at those intensities, but that is 2–15 times today's record.

## 3. Boosters

### Spin-polarised D-T

Aligning the D and T spins raises the D-T cross section by 1.5 (theory, never yet tested in a plasma; the DIII-D polarised-fuel programme aims to measure whether polarisation survives).

| route | unpolarised | polarised |
|---|---|---|
| deuteron fired into a 10 keV D-T plasma: fusion energy / ion energy (best ion energy) | 2.2 | 3.3 |
| deuteron fired into a 30 keV D-T plasma: fusion energy / ion energy (best ion energy) | 3.1 | 4.7 |
| ball of exploding D-T crystals (generous bound): breakeven ball radius and energy | 1.6 mm, 1.1e+07 J | 0.53 mm, 4.2e+05 J |
| ball of exploding D-T crystals (conservative bound): breakeven ball radius and energy | 4.2 mm, 7.8e+07 J | 2.8 mm, 2.3e+07 J |

Polarisation is the biggest lever on the fuel itself in this project. It cuts the breakeven energy of a ball of exploding crystals 3–26-fold, because burn in the ball rises steeply with reaction rate. But the ball is still uncompressed inertial fusion. Its energy must arrive in tens of picoseconds, and light cannot reach crystals deeper than ~10 nm, so nothing changes in practice. For compressed fusion, polarisation has been studied as a way to lower ignition thresholds (e.g. Temporal et al., Nucl. Fusion 52, 103011 (2012)); there it helps real schemes, not these crystals.

### Stripping stuck muons with X-rays

A stuck muon sits in a muonic helium ion bound by ~11 keV. A Missouri S&T thesis found that an X-ray fluence of 10²³ photons/cm² would lift muon catalysis to >600 fusions per muon (about 450 more, worth ~7.9 GeV). Even at the minimum photon energy, 11 keV, that fluence carries 1.8×10⁸ J per cm². Paying for it needs ~1.4×10¹⁷ muons in each cm² of irradiated fuel at once. The most intense muon beams deliver ~10⁸ per second, so that is billions of beam-seconds of muons for every square centimetre. A 2026 rate-network analysis (arXiv:2606.07077) finds external-field stripping adds at most ~40% (113 → 157 fusions per muon in its best case).

