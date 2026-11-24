# Every route in one simulator: the search for net energy

Produced by `tools/simulator_search.py`. Code: `nanofusion/improved/simulator.py` (`python -m nanofusion sim list | run | search`).

## The account

Each route of docs/energy.md is wrapped in the same energy account:

electricity → driver (laser light, accelerator beam or pulsed current) → energy in the fuel → fusion → heat (with a lithium blanket that breeds tritium: +4.78 MeV per D-T fusion) → electricity.

* Driver efficiencies: laser 10% (fusion-laser designs; NIF's own is ~0.5%), pulsed power 15%, accelerators 40% (a design figure). Light reaching crystal ions or fuel heat: 15%; light reaching squeezed fuel: 1%. Heat to electricity: 40%.
* **Electric gain** = electricity out / electricity in. Net energy means more than 1.
* **Stretch**: each design lists what it demands (laser energy and power, field, density, pressure, current, how long a pinch stays stable, muon supply, ...) against the best ever achieved. The largest ratio is its stretch. Spin-polarised fuel, never shown to survive in a plasma, counts as infinite stretch.
* The search samples each route 300 times (log-uniform over each range, including the ends, with and without the boosters), then hill-climbs 100 steps from the best design within each budget.

## Best electric gain of each route, by how far beyond today's technology it may go

| Route (section) | today's technology | ≤ 10× today | ≤ 100× today | unlimited |
|---|---|---|---|---|
| fibre_pinch (19) | 0.016 | 0.86 | **2.8** | **5.9** |
| pycnonuclear (10) | – | – | – | **5.3** |
| tamped_ball (12) | 8.3×10⁻⁴ | 5.8×10⁻³ | 0.01 | **1.7** |
| muon_squeezed (17) | – | – | – | **1.1** |
| fast_ignition (7) | – | – | – | 0.41 |
| crystal_ball (3) | 6.4×10⁻⁴ | 1.6×10⁻³ | 2.6×10⁻³ | 0.19 |
| muon_catalysis (9) | 0.12 | 0.16 | 0.16 | 0.16 |
| son_fisch_pB (17) | – | 5.9×10⁻³ | 0.01 | 0.015 |
| channelling (20) | 5.4×10⁻³ | 5.4×10⁻³ | 5.4×10⁻³ | 5.4×10⁻³ |
| opaque_pB (18) | – | – | 1.4×10⁻³ | 4.3×10⁻³ |
| pauli_beam (17) | 4.8×10⁻⁴ | 1.4×10⁻³ | 2.1×10⁻³ | 3.5×10⁻³ |
| projectiles (16) | – | – | 9×10⁻⁴ | 9×10⁻⁴ |
| cold_catcher (1) | 3.3×10⁻⁵ | 3.3×10⁻⁵ | 3.3×10⁻⁵ | 3.3×10⁻⁵ |

Bold: more than 1, net electricity on paper.

## The winning designs

### Today's technology

**muon_catalysis** (Muon-catalysed fusion in liquid or solid D-T): phi 1.4, sticking PSI 1987 (0.45%).

* Fusion 4.4×10⁻¹⁰ J from 7.5×10⁻¹⁰ J of accelerator; electricity out / in **0.12**. 155 fusions per muon; section 9.
* liquid or solid D-T density (× liquid hydrogen): 1.4 against 1.4 (solid D-T), ×1

### ≤ 10× today

**fibre_pinch** (D-T fibre pinched by its own current): current_ma 62, t_kev 9.7, alfven_times 9.8, squeeze 48.

* Fusion 3.3×10⁸ J from 2.9×10⁷ J of pulsed power; electricity out / in **0.86**. 1 cm of D-T fibre; fusion / plasma heat 110; X-rays to resupply 0 x the heat; section 19.
* pinch current (A): 6.2×10⁷ against 2.6×10⁷ (Sandia's Z), ×2.4
* energy into the pinch (J): 2.9×10⁷ against 3×10⁶ (Z (2.7 MJ of X-rays from ~20 MJ stored)), ×9.8 ← the hardest
* stable Alfvén times of a dense hot pinch: 9.8 against 1 (fibre pinches broke up at peak current), ×9.8
* radial squeeze of the fuel: 48 against 35 (NIF-class implosions, 30–40), ×1.4
* hot-spot pressure (Gbar): 3,200 against 600 (NIF's best shots), ×5.3

### ≤ 100× today

**fibre_pinch** (D-T fibre pinched by its own current): current_ma 100, t_kev 7.1, alfven_times 100, squeeze 26.

* Fusion 2.3×10⁹ J from 6.3×10⁷ J of pulsed power; electricity out / in **2.8**. 1 cm of D-T fibre; fusion / plasma heat 310; X-rays to resupply 0 x the heat; section 19.
* pinch current (A): 10⁸ against 2.6×10⁷ (Sandia's Z), ×3.8
* energy into the pinch (J): 6.3×10⁷ against 3×10⁶ (Z (2.7 MJ of X-rays from ~20 MJ stored)), ×21
* stable Alfvén times of a dense hot pinch: 100 against 1 (fibre pinches broke up at peak current), ×100 ← the hardest
* radial squeeze of the fuel: 26 against 35 (NIF-class implosions, 30–40), ×1
* hot-spot pressure (Gbar): 720 against 600 (NIF's best shots), ×1.2

### Unlimited

**fibre_pinch** (D-T fibre pinched by its own current): current_ma 100, t_kev 6.2, alfven_times 2,200, squeeze 1.

* Fusion 2.3×10⁹ J from 3×10⁷ J of pulsed power; electricity out / in **5.9**. 1 cm of D-T fibre; fusion / plasma heat 310; X-rays to resupply 0 x the heat; section 19.
* pinch current (A): 10⁸ against 2.6×10⁷ (Sandia's Z), ×3.8
* energy into the pinch (J): 3×10⁷ against 3×10⁶ (Z (2.7 MJ of X-rays from ~20 MJ stored)), ×10
* stable Alfvén times of a dense hot pinch: 2,200 against 1 (fibre pinches broke up at peak current), ×2,200 ← the hardest
* radial squeeze of the fuel: 1 against 35 (NIF-class implosions, 30–40), ×1
* hot-spot pressure (Gbar): 0.89 against 600 (NIF's best shots), ×1

## What the search found

* **With today's technology nothing comes close.** The best is muon catalysis at 0.12: a muon returns about a tenth of the electricity it cost.
* **One route crosses 1 on paper: the dense D-T Z-pinch** (section 19). Within ten times today's best on every demand at once it reaches 0.86. Allowed a hundred times, it reaches 2.8, limited by stable Alfvén times of a dense hot pinch. That needs the pinch held stable about 100 Alfvén times at 100 MA. No dense, hot pinch has stayed stable past its peak current; Zap Energy's sheared-flow pinches last 1,000–10,000 growth times, but at a millionth of solid density. The 1980 Los Alamos dense Z-pinch reactor study also found net power on paper, assuming stability.
* **The crystal routes never get near.** Exploding crystals stay below a percent of breakeven within 100 times today's technology, held back by picosecond laser energy. The tamped, magnetised ball (section 12), the best of them, needs a picosecond laser far larger than any.
* **Unlimited technology** lets several routes pass 1 on paper: pycnonuclear burning at white-dwarf density, a magnetised ball with spin-polarised fuel, muons in squeezed fuel. Each needs something many orders of magnitude beyond anything built.
* **Corrections the search forced.** Its first runs found a pinch "burning" more than all its fuel, a pinch at 130,000 Gbar, and one at 3.6 keV whose X-rays would have cost 17 times its drive. Fuel burn-up, a pressure demand and X-ray resupply are now in the model; those designs are gone.

## With a uranium blanket (k_eff = 0.95)

A subcritical uranium blanket turns each D-T neutron into fission heat: about 270 times the fusion energy. That is a fission reactor driven by a neutron source (section 13).

| Route | today's technology | ≤ 10× today |
|---|---|---|
| fibre_pinch | 0.44 | **130** |
| muon_catalysis | **25** | **34** |
| tamped_ball | 0.17 | **1.3** |
| channelling | **1.2** | **1.2** |
| pauli_beam | 0.083 | 0.29 |
| crystal_ball | 0.1 | 0.16 |
| cold_catcher | 6.8×10⁻³ | 6.8×10⁻³ |
| son_fisch_pB | – | 5.9×10⁻³ |
| fast_ignition | – | – |
| muon_squeezed | – | – |
| pycnonuclear | – | – |
| projectiles | – | – |
| opaque_pB | – | – |

* **Fission makes several routes "net energy"**, muon catalysis with today's accelerators among them, but over 99% of the heat is fission. That is the conclusion of section 13 again.

