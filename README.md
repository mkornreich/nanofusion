# NanoFusion Simulation (Python, in the browser)

Simulates a tiny **lithium-deuteride crystal hit by an intense mid-infrared
laser**. Bare deuterium and lithium nuclei sit on LiD's lattice sites, 2.04 Å
apart, and repel each other through Coulomb forces. The laser's oscillating
electric field (wavelength 3.6 µm) drives them, with an optional Lorentz
force. The run reports **fusion** if any two neighbouring nuclei come closer
than a threshold distance (default 2×10⁻¹⁵ m).

It is a numerical test of "guided nuclear laser fusion", an idea Philipp
Kornreich and Carlos Hartmann (Syracuse University) presented at Photonics
North 2011. [docs/what-it-simulates.md](docs/what-it-simulates.md) covers:
* where the idea comes from;
* how each setting maps to physics;
* why the model shows the "guiding" effect but cannot reach fusion energies at
  realistic laser strengths.

This is a Python port of the original NetBeans Platform (Java 6) application.
That application is kept unchanged in [`java/`](java/), extracted from
`nanofusion.adding`.

## Running

Requirements: Python 3.10+ and numpy. No other packages and no internet
connection are needed.

```bash
python3 -m nanofusion
```

This starts a local web server on http://localhost:8765/ and opens it in your
browser; press Ctrl+C to stop it. You can also open a saved simulation at
start-up, including files written by the old Java app:

```bash
python3 -m nanofusion java/FusionSimulation/sim6x6x6lor.json
```

Options: `--port N`, `--no-browser`, and `--host ADDRESS`. By default the
server listens only on this computer (127.0.0.1).

To run a simulation without the browser:

```bash
python3 -m nanofusion run --size 6 6 6 -n 1000 --efield 3e12 0 0 --lorentz -o run.json
```

Run `python3 -m nanofusion run --help` for all options. A 6×6×6 lattice takes
about 8 s per 1000 steps.

### The improved model

The original model can't answer "can this produce fusion?". It leaves out the
Coulomb explosion, pushes the ions with an unphysical always-on field, and
calls a pair "fused" when it comes within a fixed distance. The **Improved
physics** model (a switch at the top of the form, or the command below)
changes that:

* Every ion feels every other ion. A small C/OpenMP kernel compiles on first
  use; without a compiler it falls back to numpy.
* The laser is a real Gaussian pulse with the correct magnetic force.
* Integration is velocity-Verlet in double precision.
* Fusion yields come from measured cross sections, both inside the crystal and
  for ions that fly out into surrounding LiD.
* It checks whether the laser can really strip the crystal, and reports the
  energy balance.

```bash
python3 -m nanofusion improved --size 2 --intensity 1e20 --lithium Li7 --phase 90
```

The results of a parameter scan (`python3 tools/scan_improved.py`) are in
[docs/improved-model.md](docs/improved-model.md). Further experiments
(`python3 tools/experiments.py`, then `python3 tools/analyze_experiments.py`)
are in [docs/experiments.md](docs/experiments.md): shock shells, colliding
crystals, a tritium catcher, and crystal size. Whether any of it can give net
energy (`python3 tools/energy_experiments.py`) is in
[docs/energy.md](docs/energy.md). It can't, short of ordinary
inertial-confinement fusion. That page also tests proton–boron fuel and its
"avalanche", electron screening, a reproduction of the Texas Petawatt
cluster-fusion record, and ion fast ignition, with literature sources.
[docs/implosion.md](docs/implosion.md) pushes the original "guided head-on
collisions" idea to its limit: a hollow fuel shell imploding onto its centre
(`python3 tools/implosion_experiments.py`, or
`python3 -m nanofusion implosion --help` for single runs). The ions' own
repulsion caps it at about 3×10⁻⁴ of the energy put in.
`python3 tools/shaped_crystals.py` tests whether disk-shaped and
double-layer crystals can beam their ions into a fast-ignition spot
(results in the fast-ignition section of docs/energy.md). They can, but only
with lasers far beyond any that exist. `python3 -m nanofusion improved
--shape disk --thickness 0.6 --dot 2` runs one.
`python3 tools/other_routes.py` models the two "cold" routes, muon-catalysed
fusion in solid D-T and NASA's lattice confinement fusion (section 9 of
docs/energy.md). `python3 tools/creative_routes.py` tries three speculative
routes: squeezing a cold crystal until it fuses, exploding crystals as a muon
factory, and spin-polarised fuel (section 10). `python3 tools/dense_fusion.py` maps fusion in
dense matter across all burning regimes (section 11).
`python3 tools/confined_ball.py` asks whether a ball of exploded crystals can
be held together, with a gold tamper, a kilotesla magnetic field and heating
by its own alphas (section 12). That cuts the breakeven energy about a
thousandfold on paper, the biggest lever found, and still gives no net
energy. `python3 -m nanofusion ball --help` burns a single ball.
`python3 tools/further_routes.py` tests four more routes: a fission blanket
around the fusion source, a chain reaction in D-T, bare ions in a magnetic
trap, and crystals fired as projectiles (sections 13–16). It also checks the
model against nanowire-array experiments (section 6).
`python3 tools/pauli_blocking.py` tries one new lever: squeezing the fuel
until its electrons are degenerate, so that Pauli blocking weakens their drag
on fast ions (section 17). A crystal ion then returns its own energy at the
density NIF's fuel is designed for, 250 times better than in cold matter. It
also tests Son & Fisch's degenerate proton–boron, and freeing stuck muons in
muon catalysis. None gives net energy: the squeeze costs more than it returns.
`python3 tools/three_walls.py` attacks three walls nothing else touched
(sections 18–20): trapping proton–boron's X-rays in an opaque ball, holding a
D-T fibre with its own current (a Z-pinch, no laser), and channelling ions
through a tritide crystal. None breaks through.
`python3 -m nanofusion sim` puts every route into one simulator with a common
energy account (list, run, search; see [docs/simulator.md](docs/simulator.md)),
and `python3 tools/simulator_search.py` searches them all for the best
electricity out per electricity in (section 21). It also runs ordinary
inertial fusion as a yardstick and combines the ideas (crystal sparks, muon
factories, pulsed-power drivers, magnetised liners, spin-polarised fuel,
blankets). With today's technology the best is a tenth of breakeven. Within
ten times today's technology only ordinary inertial fusion (1.5) crosses net
electricity, with the dense Z-pinch next; no combination with exploding
crystals gets above ~1% of breakeven. Untested theory (alpha channelling,
direct conversion, freeing stuck muons with fields, antimatter sparks,
light-assisted tunnelling, spin-polarised fuel) is kept apart and searched
separately (`python3 -m nanofusion sim search --theory`, section 22). It mostly
moves the numbers, not the verdict. The exception is a dense pinch kept stable
by sheared flow, which would lead within ten times today's technology.
`python3 tools/icf_experiment.py` and `python3 tools/icf_dynamics.py` design a
practical experiment with ordinary inertial fusion: a one-change-at-a-time shot
campaign on NIF (section 23; [docs/icf_experiment.md](docs/icf_experiment.md)).
A 1-D implosion-and-burn model (`nanofusion/improved/burn1d.py`) simulates what
the static model leaves out and matches LLNL's own 1-D yields. The campaign is
mostly not novel: NIF has already falsified the static model's claim that the
adiabat doesn't matter. The open test is spending extra laser energy on fuel
rather than ablator.
`python3 tools/zpinch_experiment.py` does the same for a dense D-T Z-pinch,
with a time-dependent pinch model (`nanofusion/improved/zpinch.py`; section 24;
[docs/zpinch_experiment.md](docs/zpinch_experiment.md)). Followed in time, a
fibre on Z can't take the current cold and collapses radiatively when
preheated. On Z it reaches 1% of breakeven even if fully stabilised; at
~70 MA, with tens of Alfvén times of stability, a millimetre-thick fibre
ignites (`dynamic_pinch` route: 3.7 within a hundred times today). D-T pinches
have never been tried and need tritium Z doesn't yet handle; preheated
deuterium fibres with a tritium trace are possible today.

### The page

* **Parameters** (left): name, iterations, time step, precision (significant
  digits), lattice size, laser field (V/m), fusion threshold, Lorentz force, and
  log-time stepping. Press **Run**: a progress bar with **Cancel** appears, then
  **Fusion?** shows YES or NO with details of the first close approach.
* **Visualizer** (right): a 3D view (drag to rotate, wheel to zoom,
  double-click to reset) or a 2D slice through one lattice layer (pick it with
  *Z slice*). Step through frames with *Restart / Previous / Next / Last*, the
  frame slider, or ←/→. *Animate* (or Space) plays the frames at the chosen
  rate. Positions are in ångström.
* **Open… / Save**: upload a saved `.json` file, or download the current
  results. The tab title shows `•` while results are unsaved, and the browser
  warns before you leave the page with unsaved results.
* **About the model**: a short version of
  [docs/what-it-simulates.md](docs/what-it-simulates.md).

The simulation runs in the Python server, which keeps the exact numerics
described below. The page only edits parameters and draws the frames. It uses
no frameworks and loads nothing from the internet, and it follows the system
light/dark theme.

## Faithful numerics

The Java engine computes everything in `BigDecimal`. Some operations are exact
and others round to *precision* significant digits. At the default precision
of 6, that rounding is about as large as one step's movement, so it shapes the
results. `nanofusion/engine.py` repeats the same operations in the same order
with Python's `decimal` module. Its output matches the original Java engine
**character for character**, which `tests/test_engine.py` checks against
reference runs of the Java code in `tests/fixtures/`. That includes a
1000-step run checked during development.

The Java engine had a few bugs that are reproduced on purpose, because they
change the numbers:

* `1/c` is rounded to 10 decimal places (`3.3E-9`).
* Positions and velocities are rounded each step.
* Charge products are rounded in actor-then-particle order.

These Java bugs are fixed, because the original either crashed or ignored the
input:

* The **fusion threshold** field is now used. Java always used the default.
* The **Lorentz force** works without an electric field. Java threw
  `NullPointerException`.
* Clearing the electric-field boxes removes the field. Java kept the previous
  run's value.
* Missing or zero **damping factor** with log time gets an error message
  instead of a crash.

## Files

Saved files are standard JSON with one frame per line (`properties`, then
`results`). Number formatting matches the Java output. The loader also accepts
the original app's non-standard layout (`properties: {...}` with a trailing
comma). If a file is truncated, the complete frames are loaded and a warning is
shown. `java/FusionSimulation/sim10x10x10.json` is such a file: it stops after
569 of 1001 frames.

## Layout

```
nanofusion/            Python package
  constants.py         physical constants (exact Decimals, as in Java)
  model.py             SimulatorProperties, Frame, SimulationResult
  engine.py            the simulation (port of Simulator/Particle/Lattice)
  fileio.py            save/load JSON (new and original formats)
  forms.py             parameter form <-> SimulatorProperties, with validation
  cli.py               `python3 -m nanofusion run`
  improved/            improved physics model: kernel.c (all-pairs forces), physics.py
                       (cross sections, stopping, laser pulse), simulation.py, cli.py,
                       energy.py (reactivities, hot-plasma slowing, burn energy balance),
                       confinement.py (tamped, magnetised ball), alternatives.py (other routes)
  web/server.py        local web server and JSON API (standard library only)
  web/static/          the page: index.html, app.js, viewer.js (canvas plots), style.css
docs/                  what-it-simulates.md (purpose, sources), improved-model.md (scan results)
tools/                 scan_improved.py, experiments.py, analyze_experiments.py, energy_experiments.py
tests/                 pytest suite (python3 -m pytest)
  fixtures/            reference outputs produced by the original Java engine
  java_reference/      Java runner + script that regenerates the fixtures
java/                  the original NetBeans project, extracted, for reference
```

To regenerate the Java reference outputs (needs a JDK):

```bash
python3 tests/java_reference/generate_fixtures.py
```
