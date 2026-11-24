# A practical experiment with ordinary inertial fusion

A shot campaign for the National Ignition Facility (NIF), designed with the
simulator's models of ordinary inertial fusion and calibrated on NIF's own
shots. It was checked against the literature for novelty, and against a 1-D
implosion-and-burn simulation for what the simple model leaves out.

```bash
python3 tools/icf_experiment.py    # the static model: results/icf_experiment/summary.md
python3 tools/icf_dynamics.py      # the 1-D model: results/icf_experiment/dynamics.md (~3 min)
```

## Is it novel?

**Mostly not**, and one of its first claims turned out to be wrong:
* **The metric is standard.** "Quality" (measured yield ÷ a clean
  implosion's yield) is the yield-over-clean ratio NIF has used since Betti
  et al. (2010). LLNL puts N240210 at ~35% of its 1-D yield, and aims for 80%.
* **The knobs have been scanned.** NIF has run one-knob scans of adiabat,
  velocity, coast time and ablator thickness (Haan & Lindl 2026; Hurricane
  2024).
* **The lower-adiabat step has been done, and it falsified the static model.**
  Hohenberger et al. (PRL 2025) lowered the adiabat alone, 2.89 → 2.68, on the
  same target at the same laser energy. They got up to 80% more yield. The
  static model predicted ×1.01.
* **Spending extra laser energy on fuel instead of ablator** was proposed by
  Landen et al. (2024). That proposal is the campaign's open question.

What this project adds is two cheap, open models, calibrated on published
shots, and a bracketed prediction for each step.

## Two models

**The static model** (`compressed_icf` in `nanofusion/improved/simulator.py`):
* an isobaric hot spot at 10 keV inside a cold shell;
* the shell burns the share ρR / (ρR + 7 g/cm²);
* it gives a "ceiling" for the yield.

**The 1-D model** (`nanofusion/improved/burn1d.py`):
* The shell and the ablator left on it coast in from peak velocity, using
  NIF's published values: 383–401 km/s, 4–6% of the ablator left, 0.45 mg/cm³
  of D-T gas.
* They stagnate, ignite or don't, burn and fly apart.
* It includes alpha heating, conduction, bremsstrahlung with re-absorption,
  and burn-up.
* It gives the same yield to 1% at 2.6× the zones, and conserves energy to
  0.5%.

| Shot | Measured | Static ceiling | 1-D clean | LLNL's own 1-D | Measured ÷ 1-D |
|---|---|---|---|---|---|
| N210808 (2021) | 1.37 MJ | 10.1 MJ | 5.1 MJ | 4.2–6.2 MJ | 27% |
| N221204 (first ignition) | 3.15 MJ | 10.9 MJ | 20 MJ | ~15 MJ | 16% |
| N250406 (record) | 8.6 MJ | 10.9 MJ | 20 MJ | – | 43% |

* **The 1-D model matches LLNL's 1-D codes**: within their range for N210808,
  and a third high for N221204's design.
* **The static ceiling is wrong both ways.** It is twice too high for
  N210808, a capsule on the edge of ignition, and about half the clean yield
  for N221204's. So the static "quality" (79% for the record) overstated how
  close NIF has come. Against the 1-D model, the record reached 43%.

## What the static model leaves out

Each piece isolated in the 1-D model on N221204's clean implosion
([dynamics.md](results/icf_experiment/dynamics.md)):

| Physics | Effect |
|---|---|
| Alpha heating | ×530 (38 kJ without it) |
| Burn propagation | 27% of the fuel burnt, against 21% from the static formula |
| The cost of pressure | at fixed kinetic energy, pressure grows as v^2.4; below ~380 km/s the capsule fails to ignite |
| The adiabat, at fixed velocity | 2.45 → 2.25: pressure +6%, areal density +5%, yield ×1.11 |
| Ablator mixed into the hot spot | 100 ng: −2%; 1000 ng: −74% |
| Asymmetry (kinetic energy that never stagnates) | 3%: −12%; 10%: −98% |

**Retrodicting Hohenberger et al.** (adiabat 2.89 → 2.68 on N210808's target):

| | Yield | Areal density |
|---|---|---|
| Static model | ×1.01 | ×1.00 |
| 1-D, clean | ×6.1 | ×1.04 |
| 1-D, with N210808's degradation | ×2.8 | ×1.03 |
| Measured | up to ×1.8 | ×1.14 |

The 1-D model gets the direction right and overshoots the size. A clean 1-D
capsule at N210808's conditions sits on the ignition cliff, where small
changes move the yield a lot; Hohenberger et al. find yield ∝ ρR¹¹ there.

## The shots

Each step's yield is bracketed two ways:
* **Fixed share:** the step's clean 1-D yield × 43%, the record's share of
  its own clean yield. This is the usual way to carry a degradation forward.
* **Fixed mixing:** the record's degradation (946 ng of mix-equivalent
  carbon) held fixed. Every step that helps pushes the capsule up the cliff,
  so this is an upper bound.

| Step | What changes | Fixed share | Fixed mixing | Status |
|---|---|---|---|---|
| A | repeat the record | **8.6 MJ** | 8.6 MJ | routine |
| B | adiabat 2.45 → 2.25, same target | **9.5 MJ** | 18.5 MJ | done in kind (Hohenberger 2025): repeat on the graded-dopant target |
| C | +4% velocity from a thinner ablator, same kinetic energy | **9.2 MJ** | 19.4 MJ | scanned in the CH era: yield peaked, then fell as more ablator burnt off |
| D0 | 2.2 MJ on 5% more ablator (what NIF did) | **9.6 MJ** | 20 MJ | done (N240210: 5.2 MJ) |
| **D** | **2.2 MJ on 7% more D-T** | **10.4 MJ** | 21.9 MJ | **proposed (Landen 2024); not found shot** |
| E | D + adiabat 2.25 | **11.1 MJ** | 24.6 MJ | – |
| G | 2.6 MJ upgrade, capsule scaled ×1.08 | **12.8 MJ** | 27.5 MJ | after the upgrade |
| I | G + adiabat 2.0 + 8% faster | **14 MJ** | 32.5 MJ | a stretch |

* **The open test is D against D0.** It asks where the last 7% of laser
  energy should go: into more fuel, or into more ablator, which is NIF's choice
  so far. The 1-D model gives more fuel an 8% edge. Landen et al. suggest a
  dudded (THD) shot first. A pair of each settles it.
* **At the upgrade:** 13–28 MJ across the bracket, and 24–26 MJ if quality
  reaches LLNL's goal of 80% of the 1-D yield. LLNL expects "greater than
  30 MJ"; Lindl et al. very likely more than 20 MJ, and 50–70 MJ optimised.

**Shot budget.** NIF fires about a dozen ignition shots a year.
* Steps A, B, D0 and D twice each are 8 D-T shots, plus THD shots for D.
* The 2.6 MJ steps wait for the upgrade (approved 2026, no date published).

## What to measure on every shot

| Measurement | NIF diagnostic | What it tests |
|---|---|---|
| Neutron yield | activation detectors, magnetic recoil spectrometer, neutron time-of-flight | yield ÷ the 1-D yield |
| Hot-spot pressure | 0-D inference from yield, ion temperature, burn width and size | the cost-of-pressure curve |
| Fuel areal density | down-scattered neutron ratio | the 1-D model's areal-density changes (B: +3–5%) |
| Ion temperature | neutron peak widths | the burn |
| Bang time, burn width | gamma reaction history, gas Cherenkov | the 1-D model's timing |
| Hot-spot shape, hot-spot velocity | neutron and X-ray imaging, neutron spectra | the asymmetry share of the degradation |
| Ablator mixed in | carbon gamma rays; X-ray brightness against neutron yield | the mixing share |

## Success and falsification

* **Repeatability:** step A's two shots set the spread.
* **D against D0:**
  * if more fuel beats more ablator by more than the spread, Landen's
    proposal holds;
  * if it loses, extra ablator buys more than its 1-D value, through
    stability.
* **The 1-D model's areal-density prediction** for B (+3% with alpha
  heating, +5% without) is tested directly by the down-scattered ratio,
  whatever the yield does. D and D0 differ by under 1% in areal density, so
  only the yield can separate them.
* **If the measured share of the 1-D yield moves with a knob,** that knob's
  real cost is in mixing or asymmetry. The 1-D model then needs that
  degradation tied to the knob.

## What it means for energy

With a lithium blanket and 40% of the heat turned into electricity, even the
stretch step returns well under the electricity a plant-grade 10% laser
would use. The simulator's net-electricity designs (section 21 of
[energy.md](energy.md)) need about ten times NIF's laser energy and areal
density. The campaign's value is a tested answer to which knobs scale yield.

## What the models still leave out

* **2-D and 3-D asymmetry.** The 1-D model folds every degradation into one
  mixing mass. For N221204 that stands in for low-mode asymmetry, which LLNL
  finds cost a factor of 4.5. That one knob is far steeper than real
  degradation, hence the bracket.
* **Two temperatures, radiation transport beyond straight-line absorption,
  the ablator's partial ionisation and tungsten dopant.**
* **The drive.** Runs start at peak velocity. The starting radius and density
  in flight are not published, so they were fitted to LLNL's 1-D yields.
* **Larger capsules** at 2.6 MJ need larger hohlraums. The kinetic energy is
  assumed to scale with the laser energy.

## Sources

Shots and design:
* Abu-Shawareb et al., PRL 129, 075001 (2022) and PRL 132, 065102 (2024):
  [manuscript](https://www.osti.gov/servlets/purl/1880954),
  [manuscript](https://escholarship.org/content/qt2px009hq/qt2px009hq.pdf)
* Kritcher et al., PRE 106, 025201 (2022) and PRE 109, 025204 (2024):
  velocities, remaining ablator, degradations,
  [article](https://journals.aps.org/pre/pdf/10.1103/PhysRevE.106.025201),
  [article](https://physics.aps.org/featured-article-pdf/10.1103/PhysRevE.109.025204)
* Lindl et al., Phys. Plasmas 33, 050503 (2026), key metrics and upgrade
  projections: [manuscript](https://www.osti.gov/pages/servlets/purl/3367337)
* Clark, Hohenberger et al., PRL 137, 225101 (2026), graded tungsten dopant:
  [abstract](https://journals.aps.org/prl/abstract/10.1103/2jhg-pclg)
* Casey et al., hot-spot velocity and mode-1 degradation:
  [manuscript](https://www.osti.gov/pages/servlets/purl/2283276)
* Ruof et al., measured implosion velocities (2DConA):
  [manuscript](https://www.osti.gov/pages/servlets/purl/3367476)
* Weber, Kritcher et al., arXiv:2609.19752, N221204 energy balance:
  [arXiv](https://arxiv.org/html/2609.19752v1)
* Regan et al., PRL 111, 045001 (2013) and Ma et al., Phys. Plasmas 24,
  056311 (2017), mix and yield:
  [manuscript](https://www.osti.gov/pages/servlets/purl/1088864),
  [manuscript](https://www.osti.gov/servlets/purl/1369464)
* LLNL, [Achieving Fusion Ignition](https://lmf.llnl.gov/science/achieving-fusion-ignition);
  [FY2025 report](https://annual.llnl.gov/fy-2025/national-ignition-facility-2025)

Prior work on the campaign's ideas:
* Betti et al., PRL 104, 135002 (2010), yield over clean:
  [abstract](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.104.135002)
* Christopherson et al., PRE 99, 021201 (2019):
  [manuscript](https://link.aps.org/accepted/10.1103/PhysRevE.99.021201)
* Hohenberger et al., PRL 135, 035101 (2025), the adiabat test:
  [manuscript](https://www.osti.gov/servlets/purl/2589558)
* Landen et al., Phys. Plasmas 31, 062712 (2024), more fuel instead of
  ablator: [article](https://pubs.aip.org/aip/pop/article/31/6/062712/3299549/What-next-Further-implosion-space-exploration-on)
* Haan & Lindl, Phys. Plasmas (2026), one-knob scans:
  [manuscript](https://www.osti.gov/pages/servlets/purl/3367335)
* Betti & Hurricane, Nature Physics 12, 435 (2016):
  [manuscript](https://www.osti.gov/servlets/purl/1255527)
* Hurricane et al., PPCF (2019), [article](https://iopscience.iop.org/article/10.1088/1361-6587/aaed71);
  Hurricane, hydro-scaling to 3 MJ: [manuscript](https://www.osti.gov/pages/servlets/purl/2497918)
* Thomas et al., Phys. Plasmas 27, 112708 (2020), pulse shaping:
  [abstract](https://www.osti.gov/pages/biblio/1711433-experiments-explore-influence-pulse-shaping-national-ignition-facility)
* Baker et al., Phys. Plasmas 30, 032702 (2023), thicker ice:
  [abstract](https://www.osti.gov/pages/biblio/2280485)
* Spears et al., Science 389, 727 (2025), Bayesian degradation inference:
  [manuscript](https://www.osti.gov/pages/servlets/purl/2589559)
