# An experiment for the thin, ice-lined liner

The simulator's search (section 21 of [energy.md](energy.md)) found its best
design within ten times today's technology in a time-dependent pinch with a
liner (`dynamic_pinch`, `nanofusion/improved/zpinch.py`). This document asks
whether it has been tried, how fragile it is, and what experiment would test
it.

```bash
python3 tools/liner_experiment.py
```

Numbers: [results/liner_experiment/summary.md](results/liner_experiment/summary.md).

## The design

| | |
|---|---|
| Liner | beryllium, outer radius 4.5 mm, wall 0.37 mm (1/12 of the radius), 1 cm long |
| Fuel | 580 µm of 50/50 D-T ice on the inner wall; D-T vapour inside at 0.3 mg/cm³ |
| Field, preheat | 4.3 T axial; almost no preheat (17 eV) |
| Driver | 34 MA in ~100 ns; a ~10 MJ bank |
| In the model | 2.9 GJ of fusion, 25% of the fuel burnt, fuel converging 177× (3.5 mm → 20 µm) |
| Electricity out ÷ in | 22 (15% of electricity to the load, 40% of heat to electricity) |

In effect it is cylindrical inertial fusion with a little magnetisation.

## Has it been done?

**No, neither as an experiment nor as a published design.** The closest is
the 70 MA "ice burner" of Sefkow et al. (2014):
* Sefkow's design had almost the same geometry: beryllium, 4.9 mm radius,
  1 cm long, ~0.7 mm of ice.
* But it used twice the current, a thicker wall (1/6), an 8 T field and
  25 kJ of preheat, and converged only 24–35×.

The optimiser has traded current for convergence. Three things about this
design are new, and they are also the least credible:

| | This design | The literature |
|---|---|---|
| Current at which the ice ignites | 34 MA | 42–56 MA for ice burners (Slutz et al. 2018); Pacific Fusion: ≥50 MA |
| Fuel convergence | ~170× | designs keep it below 40; measured 30–50 at stagnation; best radiographed liner convergence 20 ± 3 |
| Liner wall | 1/12 of the radius | designs use 1/6; Z's record shots used coated liners, ~1/8.4 overall |

Other prior work:
* Slutz & Vesey (2012), Slutz (2016, 2022) and Pacific Fusion (2025, 2026)
  designed ice-layered liners at 56–70 MA.
* Lindemuth (2024) studied an unmagnetised, unpreheated liner at 0.3 mg/cm³,
  converging 32×. He calls convergence of 80–120× "highly unattainable".
* Slutz (2022) notes that D-T vapour above the ice (0.3 mg/cm³, this
  design's fill) is too thin for usual MagLIF designs.

**On the experimental side:**
* No ice-lined liner has been imploded on Z.
* Sandia has grown deuterium ice in bench cells and fielded cryogenic *gas*
  fills at ~72 K. Its tool for laying ice inside a liner is not ready.
* Thin liners on Z need coatings (the record shots used parylene-coated
  1/10.6 beryllium). Stagnation grows more unstable as the wall thins.

## What the model says

**All or nothing on convergence.** If the liner breaks up before the fuel
converges this far, compression stops there:

| Liner breaks up at | Fusion |
|---|---|
| 40× (about what MagLIF has shown) | 25 J |
| 60× | 0.5 kJ |
| 100× | 12 kJ |
| 150× | 0.17 MJ |
| never (177×) | 2.9 GJ |

**A knife-edge.** Small changes kill it:

| Change | Fusion drops by |
|---|---|
| 30 MA instead of 34 | ~17,000× |
| a wall 1/9 of the radius instead of 1/12 | ~8,100× |
| 1 mg/cm³ of gas (MagLIF's usual) instead of 0.3 | ~110,000× |
| 20 T instead of 4.3 | ~92× |

An optimiser's window this narrow is more often a model artefact than robust
physics. The model has no liner instability. So it lets thin, fast liners
squeeze far further than real ones, and it ignites wherever the squeeze is
deepest.

## The experiment

Each stage costs more than the last, and each can stop the programme.

**Stage 1, on Z now: the same target in deuterium (~$20–50M, ~12 shots).**
* **Targets:** coated beryllium liners, wall 1/10–1/12 of the radius, lined
  with deuterium ice, at the ~20 MA Z drives into MagLIF loads. Add a 1/6
  wall control.
* **Development needed:** deuterium ice laid inside a liner. Bench growth is
  shown; the tool isn't ready.
* **Measure:**
  * the stagnated column's radius (X-ray and neutron imaging);
  * liner integrity in flight (radiography);
  * D-D neutrons;
  * secondary D-T neutrons from D-D tritons, which give the fuel's areal
    density and field.
* **What the model predicts at 20 MA**, thin wall:

| Outcome | Stagnation radius | D-D neutrons |
|---|---|---|
| The liner holds (the model alone) | ~10 µm (~390×) | ~1.5×10¹² |
| It breaks up at 100× | ~35 µm | ~3×10¹¹ |
| It breaks up at 40× (MagLIF-like) | ~90 µm | ~8×10¹⁰ |

* **The decision:**
  * If the column stops near 0.1 mm, as every liner so far has, the design is
    dead: it needs ~170×.
  * Only a column below ~30 µm (100× or more) justifies going on.

**Stage 2: D-T, on a driver of ~34 MA.**
* Z delivers ~18–22 MA to such loads; even its unfunded upgrade ("ZX")
  wouldn't reach 34 MA.
* The machines that could are Pacific Fusion's >60 MA Demonstration System,
  run below full current (ground broken 2026, target ~2030), or NNSA's
  unbuilt next generation.
* It also needs uncontained 50/50 D-T on a pulsed-power machine. Z has fired
  only trace tritium.

**Stage 3: an ignition shot at 34 MA**, only if both earlier stages pass.

## What it would mean for energy

If it worked as modelled, each shot would release ~2.9 GJ (0.7 t TNT
equivalent) from a ~10 MJ bank. That is 22× the electricity in: the best of
any design here, at a current only 1.3× Z's. The literature and the model's
own fragility say the convergence it needs is the most likely thing to fail.
The 80 MA liner (section 25), at MagLIF's own ~43× convergence, is the design
to believe. Stage 1 is the cheap way to find out which.

## Limits

* **0-D, incompressible liner, no instability.** Convergence is the model's
  free variable, which is why the knife-edge appears.
* **The convergence limit** is a sharp cut-off. Real break-up degrades
  compression gradually and mixes liner into the fuel.
* **Calibration:** at Z's scale the liner model overshoots clean 2-D MagLIF
  yields ~10×. At high current it undershoots published designs (section
  25).

## Sources

* Sefkow et al., Phys. Plasmas 21, 072711 (2014): [article](https://pubs.aip.org/aip/pop/article/21/7/072711/915357)
* Slutz & Vesey, PRL 108, 025003 (2012): [abstract](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.108.025003)
* Slutz et al., Phys. Plasmas 23, 022702 (2016): [article](https://pubs.aip.org/aip/pop/article/23/2/022702/318337)
* Slutz et al. (2018), ice-burning thresholds: [manuscript](https://www.osti.gov/pages/servlets/purl/1467390)
* Slutz et al. (2022), a 60 MA ice design: [manuscript](https://www.osti.gov/pages/servlets/purl/1855793)
* Slutz et al., Phys. Plasmas 17, 056303 (2010): [article](https://pubs.aip.org/aip/pop/article/17/5/056303/922150)
* Pacific Fusion (2025), AMPS: [arXiv:2504.10680](https://arxiv.org/html/2504.10680); (2026): [arXiv:2604.21149](https://arxiv.org/html/2604.21149)
* Lindemuth, Phys. Plasmas 31, 072707 (2024): [article](https://pubs.aip.org/aip/pop/article/31/7/072707/3303929)
* Weis et al., Phys. Plasmas 33, 062711 (2026), record shots and ice layers in simulation: [article](https://pubs.aip.org/aip/pop/article/33/6/062711/3395915)
* Liner convergence 20 ± 3 radiographed: [manuscript](https://link.aps.org/accepted/10.1103/PhysRevLett.116.065001)
* Stagnation convergence from X-ray imaging: [article](https://pubs.aip.org/aip/pop/article/31/2/022703/3266623)
* MagLIF overview (convergence below 40): [manuscript](https://www.osti.gov/pages/servlets/purl/1834113)
* Aspect-ratio comparison: [manuscript](https://www.osti.gov/servlets/purl/1476822)
* Cryogenic MagLIF fills: [manuscript](https://www.osti.gov/pages/servlets/purl/2311505); deuterium ice growth, Awe et al., RSI 92, 073505 (2021): [article](https://pubs.aip.org/aip/rsi/article/92/7/073505/1031015)
* Z current and plans: [OSTI](https://www.osti.gov/servlets/purl/1429410); next-generation pulsed power: [OSTI](https://www.osti.gov/servlets/purl/2003885/)
* Tritium on Z: [Sandia Lab News, 2016](https://www.sandia.gov/labnews/2016/10/27/tritium/)
