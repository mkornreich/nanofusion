# Pauli blocking: switching off the electrons' drag

Produced by `tools/pauli_blocking.py`. Code: `nanofusion/improved/degenerate.py`.

Every route in docs/energy.md hits the same wall: electrons drag fast ions to a stop before they fuse. In fuel squeezed until its electrons are degenerate, an ion slower than the Fermi velocity can only kick the few electrons near the Fermi surface, because every lower state is full. The drag per metre then stops growing with density (Fermi & Teller 1947), while fusion per metre keeps growing.

## Model checks

**Its two textbook limits.** The electron stopping is Lindhard's random-phase approximation for a free electron gas at zero temperature, integrated numerically. The stopping number L (dE/dx = 4π Z²e⁴ n_e L / m_e v²):

| Electron density (× solid D-T) | L at v = 0.1 v_F | Lindhard–Winther friction | L at v = 20 v_F | Bethe, ln(2mv²/ħω_p) |
|---|---|---|---|---|
| 10 | 0.000455 | 0.000453 | 7.23 | 7.23 |
| 1000 | 0.00104 | 0.00104 | 8 | 7.99 |
| 100000 | 0.00175 | 0.00174 | 8.76 | 8.76 |

The low-speed coefficient also matches Lindhard & Winther's own table (Mat. Fys. Medd. Dan. Vid. Selsk. 34, no. 4 (1964), Table 1) to within 0.5%.

**A real degenerate electron gas: slow protons in aluminium.** Aluminium's three valence electrons are close to a free gas, and their stopping of slow protons is measured, so this tests the Pauli-blocked regime directly.

| Proton energy | v / v_F | This model (valence electrons only) | Measured (NIST PSTAR) | Model / measured |
|---|---|---|---|---|
| 10 keV | 0.68 | 196 MeV cm²/g | 292 MeV cm²/g | 0.67 |
| 25 keV | 1.08 | 321 MeV cm²/g | 413 MeV cm²/g | 0.78 |
| 50 keV | 1.53 | 500 MeV cm²/g | 473 MeV cm²/g | 1.06 |
| 100 keV | 2.16 | 397 MeV cm²/g | 447 MeV cm²/g | 0.89 |

* **Slower than the electrons, the model drags too little.** Below the Fermi velocity it gives two-thirds of the measured stopping. Measured slow-proton friction in metals is about 1.6 times the linear-response value, because the electrons pile up around the proton (Montanari & Miraglia, [arXiv:1705.01575](https://arxiv.org/abs/1705.01575)). That pile-up weakens as the electrons get denser, and squeezed fuel is about six times denser in this sense (r_s ≈ 0.3 against 2.1). The tables below also give the results with the drag raised 1.6-fold.
* **A hint from NIF.** Knock-on ions crossing NIF's cold, partly degenerate fuel showed a degeneracy-induced shift of their Bragg peak, at 90% confidence (Hayes et al., Nat. Phys. 16, 432 (2020), [OSTI](https://www.osti.gov/pages/servlets/purl/1881627)).


**At solid density** this free-electron model gives a deuteron in D-T a return of 0.003 at 236 keV. The project's cold-matter model, with bound electrons, gives 0.006. Below a few times solid density the electrons are bound, so the free-electron model is only a rough guide there. Above ~100 times solid density the fuel is pressure-ionised and the model applies.

## 1. Crystal ions fired into cold, squeezed fuel

A deuteron slowing down in cold D-T squeezed to each density. The return is fusion energy (all 17.6 MeV) / the deuteron's energy, as in sections 1–2 of docs/energy.md. Solid D-T is 0.25 g/cm³; a NIF-class implosion reaches ~1,000 g/cm³.

| D-T density (g/cm³) | Fermi energy | Least work to squeeze, per electron | 100 keV | 300 keV | 1 MeV | 3 MeV | Best (deuteron energy) | Its energy to electrons | Its path |
|---|---|---|---|---|---|---|---|---|---|
| 2.5 | 26 eV | 16 eV | 0.0075 | 0.0077 | 0.0037 | 0.0021 | **0.0099** (160 keV) | 97% | 1.4 µm |
| 25 | 120 eV | 72 eV | 0.048 | 0.05 | 0.019 | 0.0077 | **0.064** (160 keV) | 91% | 0.83 µm |
| 250 | 560 eV | 330 eV | 0.28 | 0.32 | 0.12 | 0.045 | **0.39** (167 keV) | 75% | 0.44 µm |
| 1,000 | 1.4 keV | 840 eV | 0.68 | 0.87 | 0.35 | 0.13 | **1** (175 keV) | 56% | 0.27 µm |
| 2,500 | 2.6 keV | 1.6 keV | 1 | 1.5 | 0.64 | 0.25 | **1.7** (191 keV) | 43% | 0.19 µm |
| 25,000 | 12 keV | 7.2 keV | 1.7 | 3.7 | 2 | 0.89 | **3.8** (236 keV) | 15% | 0.053 µm |
| 250,000 | 53 keV | 32 keV | 2 | 4.8 | 3.4 | 2.2 | **4.8** (269 keV) | 3% | 0.0085 µm |

* **The cold-catcher wall breaks.** In cold matter at normal density a deuteron returns at most ~0.4% of its energy (section 1). At the density NIF's fuel is designed for (1,000 g/cm³) it returns 1.0. Above about 940 g/cm³ it returns **more than its own energy**, the first time in this project for an ion in unheated fuel.
* **Why.** Below the Fermi velocity the drag per metre no longer grows with density. Fusion per metre does. At the highest densities the deuteron loses most of its energy to fuel ions instead, and that caps the return near 5.
* **With the drag raised 1.6-fold** (as measured in metals): 0.73 at 1,000 g/cm³, 1.3 at 2,500 g/cm³, 3.4 at 25,000 g/cm³; the crossing moves to about 1,600 g/cm³.
* **Squeezing buys what heating buys.** Degenerate electrons drag on a slow ion about as much as a plasma at their Fermi temperature does (section 2's hot catcher), at a fifth of the energy and without radiating. That is the whole trick.

### Protons into squeezed boron

| Boron density (g/cm³) | Fermi energy | Best return (proton energy) | Its energy to electrons |
|---|---|---|---|
| 2.34 | 27 eV | 0.0015 (5335 keV) | 100% |
| 234 | 580 eV | 0.0045 (816 keV) | 96% |
| 2,340 | 2.7 keV | 0.032 (816 keV) | 87% |
| 23,400 | 12 keV | 0.21 (816 keV) | 62% |
| 234,000 | 55 keV | 0.74 (4278 keV) | 66% |

* **Proton–boron stays below 1** even at 234,000 g/cm³ (0.74). Boron's charge of 5 makes collisions with the nuclei 25 times stronger, so they take over from the electrons as the drag.

## What the squeeze costs

The trick lasts only while the fuel stays cold. Once the beam has added about as much heat as the squeeze put in, the fuel is a hot plasma (section 2), and the same ions would be better spent igniting it. So the fairest test of the trick alone gives the beam as much energy as the squeeze:

* perfect drivers: return = best return / 2;
* real drivers: 15% of the light into ions (section 7), and 1% of the light into the squeezed fuel (a NIF-class implosion); then a 10%-efficient laser for the wall plug.

| D-T density (g/cm³) | Best return per ion | Return, perfect drivers | Return on the light | From the wall plug | Hot spark that ignites the same fuel |
|---|---|---|---|---|---|
| 250 | 0.39 | 0.2 | 0.0037 | 0.00037 | 3.3 kJ |
| 1,000 | 1 | 0.52 | 0.0097 | 0.00097 | 210 J |
| 2,500 | 1.7 | 0.87 | 0.016 | 0.0016 | 33 J |
| 25,000 | 3.8 | 1.9 | 0.035 | 0.0035 | 340 mJ |
| 250,000 | 4.8 | 2.4 | 0.045 | 0.0045 | 4.5 mJ |

* **The squeeze eats the gain.** Even with perfect drivers the fuel returns only about its best return / 2. Real implosions put about 1% of the light into the fuel, so the light comes back at 1.0%–4.5%, and the wall plug at a few parts per thousand at best.
* **At these densities ignition is cheap.** A hot spark of joules to kilojoules ignites the whole squeezed fuel (section 11's cheapest hot ball). Ions are far better spent as that spark (fast ignition, section 7) than on in-flight fusion. The real problem is the squeeze itself: 2–250 times NIF's density, without heating.

## 2. Hot ions in degenerate electrons (Son & Fisch)

**Idea** (Son & Fisch, Phys. Lett. A 329, 76 (2004)). Proton–boron fails because its electrons drain the hot ions' energy and radiate it (section 4). Degenerate electrons drain less and radiate less. An ion much slower than the electrons loses energy at a fixed rate, whatever the density, while fusion grows with density. So above some density the ions' own fusion products should keep them hot. Son & Fisch found a possible self-burning regime above ~10⁵ g/cm³, with ions near 100 keV and electrons near 30 keV, and a gain of ~18 on the electrons' Fermi energy (40–200 in a 2006 follow-up, Phys. Lett. A 356, 72). They noted themselves that a power plant needs far more, and that squeezing fuel to such densities is an unsolved problem ([preprint](https://bp-pub.pppl.gov/pub_report/2004/PPPL-4003.pdf)).

**Test.** Steady power balance per unit volume, for proton–boron with 15% boron:
* the drain: each ion loses energy at dE/dt = −(2η/M)(E − 3T_e/2). η is the Lindhard–Winther friction of degenerate electrons (with the relativistic mass correction), merging into Spitzer's once the electrons are hotter than their Fermi energy (generous by up to 2 near the join);
* the heating: the share of the alphas' energy that goes to ions, from the stopping model above (target ions at rest, which is generous);
* the electrons must radiate what they receive. Pauli blocking cuts their bremsstrahlung by about 2.2 (T/T_F)^1.5 (Son & Fisch give the scaling; the prefactor is derived here). That fixes the electron temperature;
* the ions keep themselves hot when charged fusion power exceeds the radiation.

| Ion temperature | Lowest density where the ions keep themselves hot | × NIF's density | Fermi energy | Electron temperature (÷ Fermi energy) | Alpha energy to ions | Burn time | Fusion in one burn time ÷ energy to squeeze and heat |
|---|---|---|---|---|---|---|---|
| 50 keV | 3.9×10⁶ g/cm³ | 3,900 | 350 keV | 27 keV (0.08) | 99% | 220 fs | 3 |
| 75 keV | 8×10⁵ g/cm³ | 800 | 150 keV | 29 keV (0.20) | 96% | 240 fs | 5.1 |
| 100 keV | 3.5×10⁵ g/cm³ | 350 | 88 keV | 31 keV (0.35) | 92% | 220 fs | 5.5 |
| 125 keV | 2.1×10⁵ g/cm³ | 210 | 65 keV | 33 keV (0.51) | 89% | 210 fs | 5.2 |
| 150 keV | 1.1×10⁵ g/cm³ | 110 | 43 keV | 40 keV (0.94) | 82% | 280 fs | 4.9 |

* **This reproduces Son & Fisch.** With 100 keV ions the model closes the balance at 3.5×10⁵ g/cm³ with electrons at 31 keV; they found > 10⁵ g/cm³ and ~30 keV.
* **Degeneracy does open a window, at white-dwarf density.** With 50–150 keV ions, proton–boron keeps itself hot at 1.1×10⁵–3.9×10⁶ g/cm³, 110 to 3,900 times NIF's design density. The fusion released in one burn time is 3–5.5 times the energy needed to squeeze and heat the fuel, before any driver losses. It is the same verdict as pycnonuclear D-T (section 10): it works on paper, but only in matter like a white dwarf's, and only after a 50–160-fold radial squeeze of a solid without heating it.
* **The electrons are only partly degenerate there.** They must radiate what the ions drain into them, so they settle at 0.1–0.9 of their Fermi energy. Pauli blocking still cuts their radiation and their drag enough.
* **Hotter ions don't need degeneracy.** With 200–300 keV ions the balance closes with the electrons at 77 keV and 120 keV, 5 and 17 times their Fermi energy: ordinary hot electrons, at any density. That is the narrow classical hot-ion window found for magnetically confined proton–boron (Liu et al. 2025, cited under "Why everything hits the same wall"). It needs confinement, not squeezing.
* **D-T needs none of this.** Its alphas outrun its bremsstrahlung above about 4 keV at any density, and fuel that dense ignites from a joule-sized spark.

## 3. Muon catalysis: freeing stuck muons

**Idea.** Muon catalysis came closest of all cold routes (section 9). It is capped by sticking: in about 0.9% of fusions the muon sticks to the alpha. The stuck ion (μα)⁺ leaves with 3.46 MeV, and collisions with fuel nuclei knock the muon free while it slows down. At liquid density this happens about 35% of the time (the "reactivation"), leaving an effective sticking near the measured 0.45–0.56%. Electrons slow the ion without freeing the muon: an electron moving with the ion carries about 0.5 keV, too little for the muon's 11 keV binding. If Pauli blocking switches off the electrons' drag, the ion meets many more nuclei before it stops.

**Earlier work.** Hot plasma lowers the drag too: Jändel et al. (Phys. Rev. A 40, 2799 (1989)) found reactivation rising with temperature and density, but nearly complete only above ~1 keV (10⁷ K), thousands of times hotter than any fuel in which muon catalysis has been measured. Cripps, Harms & Goel (Hyperfine Interact. 77, 181 (1993)) found an energetically useful number of fusions per muon only below ~5 eV and above 100 times liquid density: the cold, dense corner tested here.

**Model.** Stripping cross section: ionisation of hydrogen by protons (Rudd et al., Rev. Mod. Phys. 57, 965 (1985)), scaled to the muonic ion (Born scaling, σ(v) = σ_H(v/2)(m_e/m_r)²/16). Slowing: this section's stopping. A single factor (×1.73) is tuned once so that liquid D-T reproduces 35%. The model leaves out the muon's excited states, which that factor absorbs.

| D-T density (g/cm³) | Muons freed | Counting only the fast part (v > muon orbit) | Effective sticking | Fusions per muon allowed by sticking |
|---|---|---|---|---|
| liquid, 0.18 (calibration) | 35% | | 0.58% | 171 |
| 2.5 | 37% | 36% | 0.56% | 177 |
| 25 | 68% | 62% | 0.28% | 351 |
| 50 | 85% | 79% | 0.14% | 718 |
| 100 | 96% | 93% | 0.035% | 2.9×10³ |
| 250 | 99.9% | 99.7% | 8.8×10⁻⁴% | 1.1×10⁵ |
| 1,000 | ≈100% | ≈100% | 3.8×10⁻¹⁰% | 2.7×10¹¹ |

* **Squeezing frees the stuck muons.** Above about 61 g/cm³ sticking falls below 0.1%, enough on paper for the ~1,000 fusions per muon that electricity needs (section 9). At 1,000 g/cm³ almost no muon stays stuck. Most of the stripping happens while the ion is still faster than the muon's orbit, where the cross section is most reliable.
* **But the squeeze is too short.** A squeezed core holds together for R / 4c_s: 140 ps for a 50 µm core at 100 g/cm³. Even if every step of the cycle were as fast as fusion inside the d-t-μ molecule (10¹² per second), a muon could catalyse only ~142 fusions in that time, fewer than the 150 it already manages in liquid fuel. The cycle's other steps have never been measured in pressure-ionised fuel, and the molecular resonance that drives them at liquid density is gone there.
* **A squeeze long enough (1 ns)** needs a core of 0.35 mm at 100 g/cm³: 18 mg of D-T and 130 kJ even if squeezed perfectly cold. To repay even that, the core must hold 4.5×10¹³ muons at once. They live 2.2 µs, so the source would have to make 2.1×10¹⁹ per second, about 2.1×10¹¹ times the most intense muon beams (~10⁸ per second). A real implosion puts ~1% of the light into the fuel, so it needs 100 times more, and their alphas would then heat the core to ~1.2 keV: no longer cold fuel, and far outside anything muon catalysis has been measured in.

