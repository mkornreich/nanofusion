# What NanoFusion was built to simulate

The original project has no documentation: no comments about its purpose and
no README. This write-up draws on four sources:

* the constants and model in the code;
* the saved runs and the file dates in the archive;
* public statements by the author, found by a web search;
* the physics literature, checked with the ported engine.

Statements marked *(inference)* are conclusions, not documented facts.

## In one paragraph

NanoFusion is a numerical test of **"guided nuclear laser fusion"**. That is an
idea by Philipp Kornreich, Professor Emeritus of Electrical Engineering at
Syracuse University, and his colleague Carlos R. P. Hartmann. The code's
author tag and Windows user name, `pkornrei`, match.

The idea:
1. Fire an extremely short, strong laser pulse at a lithium-deuteride (LiD)
   nanocrystal.
2. The pulse strips the electrons in under a femtosecond. The heavy nuclei are
   briefly left on their lattice sites.
3. The laser field, aligned with a crystal axis, pushes the light deuterium
   nuclei faster than the lithium nuclei, so each deuteron is "guided" into
   the lithium nucleus ahead of it.

The program puts bare D and Li nuclei on the LiD lattice, shakes them with the
laser field for about 20 fs, and reports "fusion" if any pair comes within
2 fm, roughly the range of the nuclear force.

## Where the idea comes from

* **Paper:** "Guided Nuclear Laser Fusion" by Philipp Kornreich and Carlos
  Hartmann, presented at Photonics North 2011 (Ottawa, 16–18 May 2011). His
  Syracuse faculty page lists it, and describes his research as "the theory
  associated with nuclear laser fusion using femto second scale deterministic
  ion motion in a LiD crystal rather than random motion in a plasma"
  ([Syracuse ECS profile](https://eng-cs.syr.edu/directory/?peopleid=3123)).
  I could not find the paper's text or abstract online, or whether it appeared
  in the conference proceedings
  ([SPIE vol. 8007](http://spie.org/Publications/Proceedings/Volume/8007)).
* **Interview** (EEWeb, 15 April 2011,
  [link](https://www.edn.com/eeweb-content/interview-with-philip-kornreich/)).
  Kornreich describes "Laser Guided Nuclear Fusion in an Optical Fiber"
  (a hollow-core fibre), worked on with Professor Hartmann:
  * A pulse of about **6 fs** hits a LiD nanocrystal with sides of about
    **30 nm**.
  * The electrons are stripped "in the first half femto second". The ions need
    "about 30 femto seconds to move significantly", so for a moment they stay
    on their crystal sites.
  * The field is aligned with the crystal axis and accelerates the deuterium
    faster than the lithium, so "the Deuterium ions are GUIDED to collide with
    the Lithium ions and fuse".
  * They claim it could run at about 2,000 pulses per second, essentially at
    room temperature.
  * They were "looking for some lab with the appropriated laser to try the
    experiment", needing a long-wavelength, very strong femtosecond laser.
  * "We have some detailed calculations." The work was never patented.
* **Collaborator:** Carlos R. P. Hartmann was a professor in the same Syracuse
  college. He died in 2015
  ([SU News](https://news.syr.edu/blog/2015/04/20/professor-carlos-r-p-hartmann-dies-service-saturday-28959/)).
* **The code came after the paper.** The saved runs date from October 2013,
  more than two years after the 2011 talk. NanoFusion looks like the numerical
  check of that proposal *(inference)*. It matches it point for point: stripped
  nuclei on fixed lattice sites, a field along the x (crystal) axis, and a
  "long" 3.6 µm wavelength. One difference: the code runs a continuous field
  for 20 fs, not a 6 fs pulse.

## How the code maps to physics

| In the code | Physical meaning |
|---|---|
| Alternating D and Li on a cubic grid, 2.04 Å apart | **Lithium deuteride**, which has the rock-salt (NaCl) structure with lattice constant ≈ 4.07 Å. The Li–D nearest-neighbour distance is a/2 ≈ 2.04 Å ([OSTI](https://www.osti.gov/servlets/purl/1502300), [OSTI review](https://www.osti.gov/servlets/purl/4062196)). |
| Charges +1 (D) and +3 (Li), no electrons | Fully ionised nuclei: the model assumes the laser has already stripped every electron. |
| Lithium mass 7.016 u | Lithium-7, the main natural isotope. |
| Force q·E₀·sin(2πct/λ), λ = 3.6 µm | A laser's electric field. 3.6 µm is the short-wavelength edge of the **deuterium-fluoride (DF) chemical laser** band, 3.6–4.0 µm ([MIRACL](https://en.wikipedia.org/wiki/MIRACL)). The period is 12 fs. |
| E₀ = 1×10¹² and 3×10¹² V/m in the saved runs | Peak intensities of 1.3×10¹⁷ and 1.2×10¹⁸ W/cm² (I = ½ε₀cE₀²). These are ultra-intense, femtosecond-pulse levels. |
| "Include Lorentz force" option | Meant to be the laser's magnetic force, q·v×B with B = E/c for a wave travelling along z. |
| Time step 2×10⁻¹⁷ s, 1000 steps | 20 fs in total, about 1.7 laser cycles. |
| Fusion threshold 2×10⁻¹⁵ m | Nuclear-force range; pairs closer than 2 fm count as "fused". |
| Reactions this would enable | D + ⁷Li → 2α + n (Q ≈ 15.1 MeV) and D + ⁶Li → 2α (Q ≈ 22.4 MeV), the lithium-deuteride reactions ([ref](https://handwiki.org/wiki/Physics:Aneutronic_fusion)). |

**History from the archive.** (Background: the idea was presented in 2011, see above.)
* The saved runs date from 2–7 October 2013: 4×4×4 to 10×10×10 lattices, fields of 1×10¹² and 3×10¹² V/m, and one run with the Lorentz option.
* The engine was reworked until July 2014 (separate X/Y/Z sizes, log-time stepping). The NetBeans application was last run on 9 July 2014.
* The project was archived in June 2022.

## The real experiments it resembles

The proposal differs from the closest established physics, **laser-driven
Coulomb explosion of atomic clusters**:

* In 1999, 35 fs laser pulses on deuterium clusters about 5 nm across produced
  D–D fusion neutrons ([Nature 398, 489](https://www.nature.com/articles/19037)).
* There, the laser strips the electrons and the remaining ball of positive
  ions blows itself apart. Its ions reach keV energies, and fusion happens when
  ions from different clusters collide.
* Molecular-dynamics studies of laser-stripped clusters like (D₂O)ₙ model this
  at 10¹⁶–10¹⁸ W/cm²
  ([PRL 87, 033401](https://journals.aps.org/prl/abstract/10.1103/PhysRevLett.87.033401),
  [PNAS](https://www.pnas.org/doi/pdf/10.1073/pnas.0508622103)).

The guided-fusion idea instead wants an ordered, head-on D→Li collision inside
the crystal, driven by the difference in how hard the field pushes D and Li,
before the crystal turns into a random plasma. I found no publication testing
it experimentally. Related lines of work include
[laser-driven ⁷Li(d,n) neutron sources](https://www.nature.com/articles/srep27363)
and [deuterated nanowire targets](https://www.nature.com/articles/s41467-018-03445-z).

## Could it have found fusion? The guiding works, but far too weakly

I checked the ported engine, which matches the Java engine exactly:

* **The field force is right.** A free deuteron at 1×10¹² V/m reaches
  1.37×10⁵ m/s after 20 fs; the engine gives 1.368×10⁵ m/s against an analytic
  1.372×10⁵ m/s.
* **The guiding effect is real.** Deuterium's charge-to-mass ratio is 14%
  higher than lithium's. Without Coulomb repulsion, a deuteron would gain
  2.75 Å on its lithium neighbour in 20 fs at 1×10¹² V/m, more than one
  lattice spacing.
* **How close they actually get.** In the full simulation (6×4×4 lattice,
  1000 steps), the closest D–Li approach over the whole run was **1.31 Å** at
  1×10¹² V/m and **0.56 Å** at 3×10¹² V/m. Fusion needs **2 fm**, about 28,000
  times closer.
* **Energy budget.** The relative D–Li energy is only about 5–80 eV.
  Bringing D and Li to 2 fm means climbing a Coulomb barrier of
  1.44 MeV·fm × 3 / 2 fm ≈ **2.2 MeV**. Quantum tunnelling doesn't close the
  gap: the probability at about 50 eV is of order e⁻⁵⁰⁰.
* **Field needed to make it work.** For the field's differential push on D
  versus Li to supply 2.2 MeV over one lattice spacing, it would have to be
  about **10¹⁷ V/m**, roughly **10²⁷ W/cm²**. That is about 10,000 times the
  most intense laser ever demonstrated, 1.1×10²³ W/cm² in 2021
  ([Optica](https://www.optica.org/about/newsroom/news_releases/2021/researchers_produce_laser_pulses_with_record-break/)).
  At such fields the "ions stay on their sites" assumption would also fail
  *(inference)*.

So, within its own assumptions, the model shows the guiding effect but cannot
produce fusion at realistic laser strengths. At the default 2 fm threshold, a
"YES" could only come from a numerical artifact.

## What the model leaves out, and why that matters

* **Electrons.** Real LiD is an ionic crystal (Li⁺D⁻). In real experiments the
  electrons absorb the laser energy. Fully ionised LiD would be about 2,800 times
  the density the 3.6 µm light can penetrate (critical density), so the light
  would reach only about 10 nm into it.
* **Long-range Coulomb forces.** The engine only lets each nucleus feel
  neighbours within two lattice steps, about 4 Å. But the energy of a Coulomb
  explosion comes from the field of the *whole* charged cluster. For a fully
  stripped LiD sphere, a surface deuteron gains about 14 eV × R², with R in Å:
  about **1.4 keV for a 1 nm cluster and 140 keV for a 10 nm cluster**. This is
  how real experiments reach fusion-relevant energies, and the neighbour cut-off
  removes it.
* **Pulse length and timescale.** The run covers 1.7 cycles with the field
  switched on abruptly. Real explosions take 100 fs to picoseconds, and the
  fusion happens over nanoseconds, between clusters.

## Numerical issues in the original code (kept in the port, to match it)

* **6-significant-digit rounding.** Positions are rounded every step. Between
  1 and 10 nm from the origin, a position can only change in steps of 10 fm,
  which is coarser than the 2 fm threshold. A particle slower than about
  250 m/s there never moves at all.
* **Explicit Euler integration.** It does not conserve energy, and distances
  are checked only once per step. At 10⁵ m/s a particle moves about 2 pm per
  step, so a brief close approach could be skipped.
* **The "Lorentz" term is missing the field.** The code adds (q·e/c)·(v_z, 0,
  −v_x). The intended force is q·e·(E₀/c)·sin(ωt)·(−v_z, 0, v_x), so the code
  omits E₀·sin(ωt). The term as written has units of charge, not force, and is
  about 10⁻¹² of the intended size. Even the correct term would be only about
  v/c ≈ 10⁻⁴ of the electric force on these ions.
* **The threshold box was ignored.** The Java form never read the field. The
  port reads it.

## Tested with an improved model

[improved-model.md](improved-model.md) redoes the simulation with every ion
interacting, a real laser pulse and measured fusion cross sections. The
guided mechanism fails:
* Li-6 can't be guided at all (same charge-to-mass ratio as D).
* Li-7 reaches only about 90 keV collisions even at 10²² W/cm².

But the crystal's Coulomb explosion could give about 0.5 fusion reactions per
30 nm crystal per pulse, if its deuterons hit surrounding LiD. That is
detectable, yet about 10⁷ times short of energy breakeven.

## What it would take to test the idea properly

* All-pairs Coulomb forces, or the cluster's mean field.
* Electrons, or at least a model of ionisation and electron expulsion.
* A pulse envelope.
* Double-precision floats.
* A symplectic integrator.
* Collision checks between steps.
* Collisions between ions from neighbouring clusters, which is where cluster
  fusion actually happens.
