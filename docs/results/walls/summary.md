# Three more walls: trap the X-rays, hold the fuel with a current, steer through the lattice

Produced by `tools/three_walls.py`. Code: `nanofusion/improved/walls.py`.

## 1. Trap the X-rays: an opaque proton–boron ball

**Idea.** Proton–boron makes no neutrons, but its electrons radiate X-rays (bremsstrahlung) faster than fusion heats the fuel (section 4). A big or dense enough ball is opaque to its own X-rays: they are reabsorbed inside and only leak out slowly from the surface. Does that save proton–boron?

**Model.** A free ball of proton–boron fuel (15% boron) at a fixed density and ion temperature. The electrons settle where the energy they receive (from the hot ions, and their share of the alphas) equals what they lose: the X-rays that escape. Escape is followed frequency by frequency (free-free absorption and Thomson scattering, Osterbrock's escape probability for a sphere). The trapped frequencies fill with blackbody radiation, and photons near m_e c² make electron–positron pairs; the ball must pay for both. It "repays its heat" if its ions keep themselves hot and it releases at least its own energy before it flies apart (R / 4c_s).

### The smallest balls that repay their heat

| Density (g/cm³) | Ion temperature | Radius | Mass | Energy it holds | Areal density ρR | Pressure | X-rays trapped |
|---|---|---|---|---|---|---|---|
| 100 | 200 keV | 2.5 mm | 6.5 g | 180 GJ | 25 g/cm² | 1.6×10⁴ Gbar | 10% |
| 1,000 | 200 keV | 190 µm | 30 mg | 720 MJ | 19 g/cm² | 1.4×10⁵ Gbar | 22% |
| 10⁴ | 200 keV | 19 µm | 290 µg | 6.5 MJ | 19 g/cm² | 1.4×10⁶ Gbar | 43% |
| 10⁵ | 200 keV | 1.7 µm | 1.9 µg | 41 kJ | 17 g/cm² | 1.4×10⁷ Gbar | 65% |

### The smallest balls that trap at least half their X-rays (300 keV ions)

| Density (g/cm³) | Radius | Mass | Energy it holds | Share of it in X-rays and pairs | Positrons per electron | Fusion released while it holds together ÷ its energy |
|---|---|---|---|---|---|---|
| 100 | 26 cm | 7.4 t | 1,800 EJ (4.4×10⁵ megatons) | 100% | 1,500 | 2.7×10⁻⁴ |
| 1,000 | 7.5 mm | 1.8 kg | 36 PJ (8.5 megatons) | 100% | 120 | 3.4×10⁻³ |
| 10⁴ | 46 µm | 4.2 mg | 190 MJ | 36% | 1.8×10⁻⁴ | 1.5 |
| 10⁵ | 1 µm | 420 ng | 12 kJ | 3% | 1.2×10⁻⁷ | 0.6 |

* **The balls that repay their heat are hot-ion balls at 17–25 g/cm².** Ions at 200 keV, electrons near 90 keV. They trap 10%–65% of their X-rays, more the denser they are. That areal density is 38–57 times the hot spot of NIF's best shots (0.44 g/cm²), at 26 to 2.3×10⁴ times their ~600 Gbar.
* **This matches the published burn studies.** Lavell et al. (Fusion Sci. Technol. 82 (2026), [OSTI](https://www.osti.gov/pages/servlets/purl/3366104)) find proton–boron burn propagating at 7,000–9,000 g/cm³ with 14–18 g/cm², and failing at 1,000 g/cm³ with 2 g/cm². Ochs et al. ([arXiv:2511.10885](https://arxiv.org/abs/2511.10885)) find breakeven needs areal densities about two orders of magnitude beyond the state of the art, at pressures three orders higher.
* **At NIF's design density** the smallest such ball holds 720 MJ (30 mg of fuel).
* **Trapping more than half the X-rays costs more than it saves at low density.** At NIF's density that takes 1.8 kg of fuel holding 36 PJ, about 8 megatons. The trapped X-rays and pairs hold 99.8% of that, and it returns 3.4×10⁻³ of it.
* **Trapping helps only where proton–boron already needs white-dwarf-like matter**, 10⁴–10⁵ g/cm³, where a micrometre-scale ball is opaque without holding much radiation.

## 2. Hold the fuel with a current: a fibre Z-pinch

**Idea.** Section 12 held the fuel with a gold tamper and a magnetic field that something else had to make. A current does both at once: run megaamperes along a fibre of solid fuel, and the current's own magnetic field squeezes it and holds it in (a Z-pinch). No laser is needed.

**Model.** Bennett's equilibrium: the field's pressure balances the plasma's, μ₀I²/8π = (N_i + N_e)kT, so each current holds a given amount of fuel per metre at 10 keV. The pinch lives some number of Alfvén times (radius ÷ Alfvén speed) before instabilities tear it. A narrower pinch is denser and burns faster, so there is a radius at which fusion repays the plasma's heat. Above 1.35 MA the pinch's own field holds the fusion alphas inside it.

**Checks against real pinches.**
* A deuterium fibre 100 µm across at 640 kA (Sethian et al., PRL 59, 892 (1987)) is in Bennett equilibrium at 136 eV. Fibre pinches of that size reached a few hundred eV.
* At 300 kA and 1 keV, Bennett holds 1.4×10¹⁹ ions per metre. Zap Energy's FuZE pinch (Levitt et al., PRL 132, 155101 (2024)), ~10¹⁷ cm⁻³ over a 6 mm radius, holds 1.1×10¹⁹.

| Current | Solid D-T fibre it holds at 10 keV | Squeeze to repay the heat in one Alfvén time (= Alfvén times needed with no squeeze) | That pinch lives | Squeeze for the alphas alone to repay it | Field energy ÷ plasma heat |
|---|---|---|---|---|---|
| 0.64 MA | 13 µm across | ×240 (1.1×10⁴ g/cm³) | 22 fs | ×1,200 (2.8×10⁵ g/cm³) | 17.4 |
| 1.4 MA | 29 µm across | ×110 (2,300 g/cm³) | 110 fs | ×560 (5.9×10⁴ g/cm³) | 15.3 |
| 5 MA | 110 µm across | ×31 (180 g/cm³) | 1.4 ps | ×160 (4,600 g/cm³) | 11.9 |
| 26 MA | 550 µm across | ×6 (6.8 g/cm³) | 37 ps | ×30 (170 g/cm³) | 7.5 |
| 60 MA | 1.3 mm across | ×2.6 (1.3 g/cm³) | 200 ps | ×13 (32 g/cm³) | 5.3 |

* **At the Z machine's 26 MA the numbers are not crazy.** A D-T fibre 550 µm across repays its heat if squeezed 6-fold and held for 37 ps, one Alfvén time; or held unsqueezed for 6 Alfvén times.
* **Fibre pinches haven't lasted.** Sethian's deuterium fibres stayed stable for ~100 growth times, but only while the current rose and only at a few hundred eV; at peak current they went unstable. Other fibre pinches (Los Alamos, Imperial College) expanded and broke up early, and their neutrons came mostly from beams, not heat. Zap Energy's sheared-flow pinch lasts 1,000–10,000 growth times, but at a millionth of solid density and 1–3 keV.
* **The field costs more than the plasma.** At 26 MA the field around the pinch (out to a return conductor 1 cm away) holds 7.5 times the plasma's heat, 3.8 MJ per cm, out of the ~20 MJ that Z stores. So even ignition-like conditions are far from net electricity unless the burn multiplies the energy many times over.
* **Lithium deuteride fibres are hopeless.** D + ⁶Li burns thousands of times more slowly than D-T: a 1 MA LiD fibre would need a 1.6×10⁵-fold squeeze. No LiD or D-T fibre pinch seems to have been tried.

## 3. Steer through the lattice: channelling in a tritide crystal

**Idea.** A fast ion entering a crystal almost along a row or plane of atoms is steered gently between them (channelling). In the open channel it meets fewer electrons, so it loses energy more slowly, and its flux piles up in the channel centre. In titanium tritide the tritons sit in the channel centres along ⟨100⟩ (but are hidden behind the rows along ⟨111⟩). Does the crystal itself beat the cold-catcher limit?

| Deuteron energy | Axial critical angle (Lindhard) | Share of an isotropic explosion along an axis (26 axes) | ... along a plane (13 planes) |
|---|---|---|---|
| 100 keV | 2.6° | 1.31% | 23% |
| 300 keV | 1.5° | 0.44% | 13% |
| 1000 keV | 0.81° | 0.13% | 7% |

A 300 keV deuteron stopping in random (non-channelling) titanium tritide returns 0.0038 of its energy. With channelling:

| Case | Flux at the tritons | Stopping ÷ random | Return |
|---|---|---|---|
| Perfectly aligned beam, generous | ×3 | 0.5 | 0.023 |
| Perfectly aligned beam, measured values | ×2 | 0.9 | 0.0084 |
| Exploding crystal (isotropic), generous, half the path channelled | axes ×3, planes ×2 | 0.5, 0.9 | 0.0041 |

Measured inputs: flux peaking of 1.4–2.2 at interstitial sites (deuterium in GaN; simulated in diamond), and channelled stopping only 5–13% below random for 50–150 keV protons in silicon, the speeds of 100–300 keV deuterons. Stopping halves only above ~1 MeV per nucleon.

* **A perfectly aligned beam does 2.2–6 times better**, 0.0084–0.023. Still 44–119 times short of breakeven. Experiments agree on the size: 20 keV deuterons on palladium gave about twice the neutrons along a crystal plane as at 30° to it (Dalkarov et al., Phys. Rev. Accel. Beams 22, 034201 (2019)).
* **An exploding crystal can't aim.** 0.4% of its ions fall inside an axial critical angle and about 13% along a plane (generous), and they dechannel partway: 0.0041, a 9% gain even with generous numbers.
* **The other way to blind the electrons, a magnetic field, needs ~6.8×10⁴ T** in solid D-T (electrons must gyrate inside a Debye length, i.e. their cyclotron frequency must exceed the plasma frequency), about 85 times the strongest laser-made field claimed.

