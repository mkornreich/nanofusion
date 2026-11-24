# More routes: fission blanket, D-T chain reaction, bare ions, crystal projectiles

Produced by `tools/further_routes.py`. Code: `nanofusion/improved/alternatives.py` and `energy.dt_chain`.

## 1. A fission blanket around the fusion source

Each D-T fusion sends a 14 MeV neutron into a blanket. In lithium alone it breeds the next triton, for 22 MeV per fusion in all. In uranium a 14 MeV neutron does much more: it splits uranium-238 directly (about 1.2 fissions) and multiplies to about 3 neutrons. If the blanket is kept just below critical (k_eff < 1), each of those starts chains worth k / (2.5 (1−k)) fissions at 200 MeV. Heat per fusion:

| uranium blanket k_eff | 0 (depleted uranium, no chain) | 0.5 | 0.9 | 0.95 | 0.98 |
|---|---|---|---|---|---|
| heat per fusion | 262 MeV | 502 MeV | 2,422 MeV | 4,822 MeV | 12,022 MeV |
| × the fusion energy | 15 | 29 | 138 | 274 | 683 |

Check: a natural-uranium pile measured 290 MeV per 14 MeV neutron (McEnhill & Weale 1962); the k_eff = 0 column gives 259. Accelerator-driven reactor designs stop at k_eff = 0.95 (MYRRHA) to 0.98 (Rubbia's energy amplifier) to stay safely below critical. Net electricity needs (driver efficiency) × (fusion / driver energy) × (blanket multiplication) × (40% heat to electricity) ≥ 1. Lasers are taken as 10% efficient and accelerators as 40%. Both are design figures: NIF's laser manages 0.5% and the best proton accelerator 18%.

| fusion source | fusion / driver energy | electricity per neutron | k_eff for net electricity | k_eff for 4 × the electricity put in | verdict |
|---|---|---|---|---|---|
| 30 nm crystals on a tritide catcher, 6 fs laser (experiments.md) | 3×10⁻⁵ | 5.9×10⁶ MeV | > 0.9999 | > 0.9999 | no: the blanket would have to be a critical reactor |
| the same if all the light became ion energy (section 1) | 4×10⁻³ | 4.4×10⁴ MeV | 0.9978 | 0.9995 | no: the blanket would have to be a critical reactor |
| hollow shell implosion (section 8) | 3×10⁻⁴ | 5.9×10⁵ MeV | 0.9998 | > 0.9999 | no: the blanket would have to be a critical reactor |
| ball of crystals at its breakeven, bare (section 3) or tamped and magnetised (section 12) | 0.15 | 1.2×10³ MeV | 0.92 | 0.98 | works on paper |
| tamped, magnetised ball returning 10 × its fuel energy (section 12) | 1.5 | 117 MeV | 0.11 | 0.79 | works on paper |
| muon catalysis at the record 150 fusions per muon (section 9) | 0.53 | 83 MeV | 0 (uranium, no chain) | 0.70 | works on paper |
| NIF's best shot, with NIF's own laser | 4.1 | 825 MeV | 0.88 | 0.97 | works on paper |
| NIF's best shot, with a 10% efficient laser | 4.1 | 43 MeV | 0 (uranium, no chain) | 0.41 | works on paper |

* **The benchmark is spallation.** A 1 GeV proton makes about 20 neutrons in lead: about 125 MeV of electricity per neutron. A fusion source is worth putting inside a blanket only if its neutrons cost about that, or a few times more, since a 14 MeV neutron is worth about three ordinary ones.
* **Exploding crystals are 350 to 47,000 times too expensive a neutron source.** They would need k_eff of 0.9978 or more. The margin to a runaway chain is the delayed-neutron fraction, 0.2–0.7%, so that blanket is a critical reactor in all but name and the crystals add nothing.
* **Any source that returns a tenth or more of its driver energy works on paper.** Muon catalysis needs only plain uranium. But then at least 91% of the heat is fission, with fission's fuel cycle and waste. It is a fission reactor with an expensive neutron source.

## 2. A chain reaction in D-T, carried by its own neutrons and alphas

k = D-T fusions triggered per D-T fusion. The 14.1 MeV neutron knocks on deuterons and tritons (up to 12.5 and 10.6 MeV), and so does the 3.5 MeV alpha; each knock-on may fuse before it stops. Infinite medium: the neutron is followed through all its collisions, which are peaked forwards as measured. Alpha knock-ons use Rutherford scattering, × 5 above 0.5 MeV to allow for nuclear scattering. Two versions of the D-T cross section above 550 keV: held constant (as elsewhere in this project, an overestimate) and the measured fall-off (Bosch & Hale's high-energy fit).

| medium (infinite) | k, cross section held constant | k, measured fall-off | of which from the alpha | chance that a 5 MeV deuteron fuses |
|---|---|---|---|---|
| cold solid D-T | 2.3×10⁻³ | 1.7×10⁻³ | 2.2×10⁻⁷ | 4×10⁻⁴ |
| solid-density D-T plasma at 1 keV | 0.048 | 0.047 | 8.4×10⁻⁵ | 5.7×10⁻³ |
| solid-density D-T plasma at 3 keV | 0.17 | 0.17 | 1.1×10⁻³ | 0.021 |
| solid-density D-T plasma at 10 keV | 0.52 | 0.49 | 0.011 | 0.068 |
| solid-density D-T plasma at 30 keV | 1 | 0.93 | 0.048 | 0.15 |
| solid-density D-T plasma at 100 keV | 1.7 | 1.3 | 0.088 | 0.29 |

A 14 MeV neutron travels about 5.3 g/cm² between elastic collisions, which is 28 cm of solid D-T (4.7 g/cm² counting all reactions). In a body thinner than that, only a share ρR / 5.3 of the neutrons collide at all. The alpha's part counts only if the alphas and their knock-ons stay inside, which needs a field (section 12):

| body | ρR | k from the neutron | k from the alpha, if held in |
|---|---|---|---|
| tamped, magnetised breakeven ball of section 12 (20 keV) | 0.005 g/cm² | 1.1×10⁻⁴ | 0.03 |
| the 100-fold-return ball of section 12 (30 keV while burning) | 0.07 g/cm² | 1.9×10⁻³ | 0.048 |
| a ball of solid D-T one mean free path in radius (30 keV) | 5.3 g/cm² | 0.092 | 0.048 |

* **Cold D-T: k ≈ 0.002.** No chain reaction, by a factor of 600. Electrons stop the knock-ons first, as in every cold route.
* **Hot D-T: k reaches about 1, but only in a plasma at 30–100 keV that is metres across.** That is section 2's result again: a fast ion in hot D-T returns 1–3 times its energy, and the neutron hands its 14 MeV to fast ions only if the fuel is many mean free paths thick. Such a plasma is already burning. The chain adds to a fire; it cannot light one.
* **In any ball this project can afford, k is a few per cent**, almost all from the alpha, and only if the alphas are held in. That is a few per cent of extra yield.
* **Published work is less generous still.** Peres & Shvarts (Nucl. Fusion 15, 687 (1975)) found cold D-T critical only above 3.5×10⁴ g/cm³, and Afek et al. (J. Phys. D 11, 2171 (1978)) found 4,200 g/cm³ at 14–15 keV. A 2026 Monte Carlo study (arXiv:2604.06769) finds D-T critical only without neutron leakage and above about 200 ρ^-0.23 keV, which is 290 keV at solid density.

## 3. No electrons at all: bare ions in a magnetic trap

A cloud of bare deuterons and tritons has no electrons to stop fast ions or to radiate, and collisions between like charges do not move them across a magnetic field. But the cloud's own charge limits its density to the Brillouin limit, n = ε₀B² / 2m (Brillouin 1945; Dubin & O'Neil, Rev. Mod. Phys. 71, 87 (1999)). Column 5 cm in radius and 1 m long, ions at 60 keV (where D-T reacts fastest):

| field | ion density | fusion power per m³ | fusion power of the column | voltage across the column | time the field must last to repay its own energy |
|---|---|---|---|---|---|
| 10 T (ordinary superconducting magnet) | 1.1×10¹¹ cm⁻³ | 7.1 W | 0.055 W | 1.2 MV | 65 days |
| 20 T (best superconducting magnets) | 4.2×10¹¹ cm⁻³ | 113 W | 0.89 W | 4.8 MV | 16 days |
| 45.5 T (strongest steady magnet) | 2.2×10¹² cm⁻³ | 3×10³ W | 24 W | 25 MV | 3 days |
| 1 kT (laser-driven coil, lasts ~1 ns) | 1.1×10¹⁵ cm⁻³ | 7.1×10⁸ W | 5.5×10⁶ W | 1.2×10⁴ MV | 9 minutes |
| 100 kT (inside a laser plasma, lasts ~1 ps) | 1.1×10¹⁹ cm⁻³ | 7.1×10¹⁶ W | 5.5×10¹⁴ W | 1.2×10⁸ MV | 0.056 s |

* **A steady trap makes watts, not megawatts.** At 10 T the limit is 10¹¹ ions per cm³, a few billionths of air's density, and 7 W per m³. A tokamak reactor needs ~10⁶ W/m³. The power grows as B⁴, but so does the difficulty: the cloud's own voltage is megavolts at 10 T.
* **Laser-made fields are strong enough but far too brief.** A 100 kT field would have to last ~0.1 s to repay its energy. It lasts a picosecond.
* **The same limit capped the imploding shell of section 8.** There the ions' repulsion limited how many fit in the focus; here it limits how many fit in the magnet.
* **Colliding beams** are the same idea in a ring. 1 W of D-T fusion needs a luminosity of 7.1×10³⁴ cm⁻² s⁻¹ at the 5-barn peak. The record for any collider is 5.1×10³⁴ cm⁻² s⁻¹ (SuperKEKB), which would give 0.7 W.

## 4. Crystals as projectiles (cluster-impact fusion)

**The original claim.** Beuhler, Friedlander & Friedman (Phys. Rev. Lett. 63, 1292 (1989)) fired singly charged heavy-water clusters, (D₂O)ₙ with n = 25–1300, at 200–325 keV into titanium deuteride, and reported about 10⁻¹² to 10⁻¹¹ D-D fusions per cluster of 100 molecules. A cluster's energy is shared by all its atoms, so each deuteron has very little. At 300 keV:

| cluster | energy per deuteron | fusions per cluster, this model |
|---|---|---|
| (D₂O)25 | 1207 eV | 7×10⁻²² |
| (D₂O)100 | 302 eV | 7.6×10⁻³⁹ |
| (D₂O)1300 | 23 eV | 1.7×10⁻⁸⁷ |

For 100 molecules this model gives 7.6×10⁻³⁹ per cluster, and an independent estimate at the time gave 4×10⁻³⁹ (Carraro et al., Phys. Rev. A 42, 1379 (1990)). The claim was about 26 orders of magnitude higher. The authors withdrew it in 1992 (Phys. Rev. Lett. 68, 2108): light ions contaminating the beam, each carrying the full energy, had probably produced the fusion.

**Doing it properly.** A grain can hold positive charge up to a surface field of about 10 V/nm before it sheds ions. A solid-deuterium crystal charged that far (generous: hydrogen ice is weakly bound) can be accelerated by a voltage to the 100–300 keV per deuteron where D-T fusion peaks. Its deuterons then stop in a tritide catcher like any beam:

| crystal diameter | charge | voltage for 100 keV per deuteron | for 300 keV | best fusion / kinetic energy | against the wall plug (accelerator taken as 40% efficient) |
|---|---|---|---|---|---|
| 3 nm | 16 e | 5.4 MV | 16 MV | 0.0041 | 0.0016 |
| 10 nm | 174 e | 18 MV | 54 MV | 0.0041 | 0.0016 |
| 30 nm | 1563 e | 54 MV | 163 MV | 0.0041 | 0.0016 |

* **It is the cold catcher again.** Whatever throws the deuteron, it returns at most 0.4% of its energy before the catcher's electrons stop it.
* **But the driver is far better.** Megavolt accelerators exist (dust accelerators throw micrometre grains at up to 80 km/s with 2 MV; this needs 3,000 km/s). Taking one as 40% efficient, this gives 0.002 from the wall plug, against ~5×10⁻⁵ for laser-exploded crystals. It is still about 600 times short, and a plain deuteron beam does the same.
* **No collective help.** A 30 nm crystal's track is heated to roughly 1 keV as it stops, which would cut the electron drag, but the heat spreads over about half a micrometre within the stopping time (~0.3 ps) and the track is cold again. (A rough estimate, not modelled.)

## 5. Check against the best nanostructure experiments: nanowire arrays

Arrays of deuterated-plastic nanowires hold the record for D-D fusion per joule of laser light: 2.2×10⁶ neutrons per joule (Curtis et al., Nat. Commun. 9, 1077 (2018)), about 500 times a flat target. Light gets several micrometres down between the wires, so a volume at a sixth of solid density is heated at once. This is the nearest real thing to the crystal foam of section 12.

A later shot with 8 J gave 1.2×10⁷ neutrons (Curtis et al., Phys. Rev. Research 3, 043181 (2021)). The authors' own plasma simulation gives 3.8×10¹² deuterons above 0.2 MeV and 3×10¹¹ above 1 MeV, an exponential spectrum with a slope of 315 keV.

* **Those deuterons stopping in cold plastic give 1.2×10⁷ neutrons in this project's model.** Measured: 1.2×10⁷. So this yield is ordinary beam-target fusion, as in section 1.
* **The fast deuterons carry 0.31 J, 4% of the light.** Nanostructures absorb most of the light (70–80% measured), but only a few per cent reaches fast ions.
* **Fusion energy / laser energy: 1.8×10⁻⁶** for that shot and 2.6×10⁻⁶ for the record. Ten times better than gas-jet clusters (section 6), and still a millionth.
* **With a tritide catcher** the same deuterons would return 1.2×10⁻⁴ of the laser energy. The cold-catcher limit holds.

