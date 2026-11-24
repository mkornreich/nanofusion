# Cold fusion routes in crystals: muon catalysis and lattice confinement

Produced by `tools/other_routes.py`. Code: `nanofusion/improved/alternatives.py`.

## 1. Muon-catalysed fusion in solid or compressed D-T

Fusions per muon Y = 1 / (λ₀/λ_c + ω_s), with λ₀ = 4.55×10⁵ /s (muon decay), λ_c the cycle rate and ω_s the sticking per fusion. The record, 150 fusions per muon with ω_s = 0.45%, needs λ_c = 2.1e+08 /s, or 1.6e+08 /s per unit of liquid-hydrogen density at φ = 1.3. Denser fuel cycles proportionally faster; the sticking does not change.

### Fusions per muon

| density φ (× liquid hydrogen) | ω_s = 0.45% (measured (PSI, Breunlich et al. 1987)) | ω_s = 0.56% (measured directly (PSI)) | ω_s = 0.65% (theory) | ω_s = 0.30% (hoped for) |
|---|---|---|---|---|
| 0.5 | 99 | 89 | 82 | 116 |
| 1 | 137 | 119 | 107 | 172 |
| 1.3 | 150 | 129 | 115 | 194 |
| 2 | 169 | 143 | 126 | 227 |
| 5 | 197 | 162 | 142 | 281 |
| 10 | 209 | 170 | 147 | 305 |
| unlimited | 222 | 179 | 154 | 333 |

Denser fuel helps only until sticking takes over: Y can never exceed 1/ω_s, about 220 with the measured sticking.

### Energy

| muon cost | fusions needed to repay the muon (17.6 MeV each) | gain at the record (150) | best possible gain (Y = 1/ω_s = 222) | with blanket heat (26 MeV per fusion), best |
|---|---|---|---|---|
| 4.7 GeV, optimised accelerator design (Kelly et al. 2021) | 267 | 0.56 | 0.83 | 1.23 |
| 5 GeV, traditional estimate (Jones) | 284 | 0.53 | 0.78 | 1.16 |
| 8 GeV, pessimistic estimate (arXiv:1008.3010) | 455 | 0.33 | 0.49 | 0.72 |
| 4e+08 GeV of laser light per muon: laser-made muons (ELI-NP, 20 J → 280 μ⁻) | 3e+10 | 6e-09 | 9e-09 | 1e-08 |

**Electricity.** Kelly et al. (2021) find electricity out / in = 14% at 150 fusions per muon for their optimised plant. Since that ratio is proportional to the fusions per muon, electrical breakeven needs about 1071 fusions per muon. That needs sticking below 0.09%, about 5 times lower than measured.

## 2. Lattice confinement fusion (NASA Glenn, 2020)

**Measured.** NASA aimed a 2.9 MeV, 15 mA electron beam (Dynamitron, tantalum converter) at 480 g of erbium deuteride (8×10²² D/cm³) in 16 vials. The fusion neutrons, scaled to the samples, came to 1.6e+03 per second, against 1.2e+03 predicted by NASA's screening theory (NASA/TP-2020-5001616, Phys. Rev. C 101, 044610).

* Beam power 43.5 kW; D-D fusion power 1.8e-09 W (3.1e+03 reactions/s × 3.65 MeV).
* **Fusion energy / beam energy = 4e-14**, before counting the accelerator's own losses.

**Model.** Photoneutrons from photons just above the 2.2246 MeV threshold (mean 145 keV) knock deuterons to ≤ 8/9 of their energy (mean 64 keV, as NASA calculated). The deuterons then fuse with the cold lattice deuterons as they stop. The cross section is taken at E + U_e, where U_e is the electron screening. Stopping in the deuteride uses this project's model, which is probably low by a factor of 2–3 in erbium, so these yields are generous.

| screening U_e | fusions per photoneutron, 1-cm vials (0.3 deuteron hits) | in an infinite sample (every hit, ~7 per neutron) | best possible energy ratio |
|---|---|---|---|
| 0 eV | 3.1e-08 | 1.5e-07 | 2e-07 |
| 300 eV | 3.2e-08 | 1.6e-07 | 2e-07 |
| 800 eV | 3.2e-08 | 1.6e-07 | 2e-07 |
| 3000 eV | 3.6e-08 | 1.8e-07 | 3e-07 |

Screening in the metal, even at 3 keV (four times the largest measured value, 800 eV in palladium), raises the yield by less than 20%. The deuterons have tens of keV, far above the screening energy. The best possible ratio assumes that every joule of beam becomes a photon just above threshold and that every photon splits a deuteron. In the real setup, most of the beam energy becomes heat and lower-energy X-rays. The samples also absorb only about 10⁻⁴ of the photons that are above threshold (about 1 mb cross section through ~1 cm). That is where the remaining factor of ~10⁷ goes. NASA's own detailed calculation of the whole chain predicted 1.2×10³ neutrons/s, close to what it measured.

**Faster neutrons.** If a 2.45 MeV (D-D) or 14.1 MeV (D-T) neutron did the knocking, in an infinite erbium deuteride sample with 800 eV screening:

| neutron | fusions per neutron | fusion energy / neutron energy |
|---|---|---|
| 2.45 MeV | 2.6e-05 | 4e-05 |
| 14.1 MeV | 0.00036 | 9e-05 |

So even a fusion neutron returns only about a ten-thousandth of its energy this way. It is the cold-catcher limit again (energy.md, section 1): a fast deuteron in cold matter loses its energy to electrons long before it fuses.

