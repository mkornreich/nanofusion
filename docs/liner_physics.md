# The liner in Athena++ with the physics a yield needs

The ideal-MHD runs ([liner_experiment.md](liner_experiment.md),
`tools/athena_liner.py`) test whether the thin, ice-lined liner holds
together. They can't say what it yields: they have no burn, and no way for
the fuel to lose heat or field. This build adds every piece a yield depends
on, each one switchable, so the yield can be measured and so can what each
piece is worth.

The quantity to maximise is the one the 0-D simulator uses
([simulator.md](simulator.md)): **electricity out ÷ electricity in**. For a
liner it is

    electric gain = 0.40 × 1.27 × Y / (E_bank / 0.15)

where:

* Y is the fusion yield;
* 1.27 is the lithium blanket (+4.78 MeV per D-T fusion);
* 40% of heat becomes electricity;
* the bank is charged at 15%.

The 0-D model gives 22 for this design (2.9 GJ from a 9.9 MJ bank).

Code:

* `tools/athena/liner.cpp`: the problem generator, with the physics compiled
  in the physics build;
* `liner_eos.cpp`: the equation of state;
* `liner_physics.hpp`: the shared formulas;
* `athena-physics.patch`: a 29-line change to Athena++ itself;
* `athinput.liner_physics`: the input, with every piece switched on.

Checks: `tools/athena/test_liner_physics.cpp`. Results:
[results/liner_athena/summary.md](results/liner_athena/summary.md)
(section "With the physics a yield needs").

```bash
git -C ../athena worktree add ../athena-physics      # a second Athena++ tree
python3 tools/athena_liner.py build-physics
python3 tools/athena_liner.py run "physics flat"     # all on
python3 tools/athena_liner.py run "physics flat, no alpha heating"   # one piece off
python3 tools/athena_liner.py report
```

## First result: the flat liner with everything on

Run `physics flat`: 1/8 of the length, no ripple, 3.2 h on 2 threads, 2026-11-25.

| | 0-D model (`zpinch`) | Athena++, all physics |
|---|---|---|
| Energy into the load | 2.6 MJ | 3.0 MJ (peak) |
| Fuel convergence | 177× (20 µm) | 84× (30 µm rms) |
| Bang time | 121 ns | 128 ns, 2.6 ns wide |
| Temperature | ignites, 80 keV peak | 4.3 keV T_i, 4.1 keV T_e (burn-weighted) |
| Field × radius | — | B_z 12 kT, BR 0.37 T m |
| Yield (1 cm) | 2.9 GJ, 25% burnt | **36 kJ**, 3×10⁻⁶ burnt |
| Electric gain | 22 | **0.0003** |

**Stagnation (130 ns).**

| Region | Radius | State |
|---|---|---|
| Hot gas | out to 28 µm | up to 4.4 Gbar and 5 keV; ρR 0.007 g/cm² |
| D-T ice | 28–127 µm | dense; ρR 0.30 g/cm² |
| Beryllium tamper | 127–520 µm | ρR 1.5 g/cm²; up to 90 g/cm³ |

* **The gas is too thin to hold its alphas.** The field holds them only
  marginally: BR is 0.37 T m, against the 0.27 T m at which an alpha's
  orbit is as wide as the column.
* **So it never runs away.** Most alpha energy goes into the ice: 7.3 kJ
  of alpha energy over the run.
* **Radiation took 82 kJ.**
* **The 0-D model ignites** because it squeezes twice as far, with fewer
  losses.

The runs with one piece switched off will show which loss costs the most.

## What is added

| Piece | Model | Switch |
|---|---|---|
| Equation of state | Beryllium and D-T. Each has a cold curve with no tension, ideal ions, and electrons that go from degenerate to ideal, with Thomas–Fermi ionisation (More's fit) and its ionisation energy. Built into tables at start-up. In mixed cells the cold parts share one pressure, side by side; the thermal parts add at one temperature. | `<hydro> eos_model` |
| Burn | D-T at the ion temperature (Bosch & Hale). The fuel becomes helium ash, and the neutrons leave with their mass and momentum. | `burn` |
| Alpha heating | The alphas lose energy on the electrons and ions as in the 0-D model: range, critical energy, a degeneracy floor, and the share to the ions. They spread across the field by radial diffusion, held in by their Larmor radius and averaged over z. So escaping alphas heat the cold ice, which is how the burn spreads. | `alpha` = 1 (0: they stop where they are made) |
| Radiation | Bremsstrahlung, recombination, and the beryllium K-shell lines (van Regemorter, coronal), carried along r (below). | `radiation` |
| Resistivity | Spitzer–Braginskii (α⊥), with the Lee–More minimum mean free path in cold dense matter. The vacuum gets its own resistivity, so the drive field can fill it. | `resistivity` |
| Electron heat conduction | Braginskii κ∥ and κ⊥ along and across the local field (Sharma–Hammett limiters), flux-limited to 0.1 of free streaming. Super-time-stepped (below). | `conduction` |
| Two temperatures | The electrons' energy is a passive scalar. Shocks heat the ions. Ohmic heat, conduction, radiation and the electrons' share of the alphas act on the electrons. Ions and electrons exchange energy at the Spitzer rate (as burn1d). | `two_temp` |
| Nernst effect | The field is carried down the electron temperature gradient at the Braginskii β∧ velocity, as an extra EMF through Athena++'s resistive EMF (the patch). Upwinded. | `nernst` |
| Circuit | A capacitor bank behind 15 nH, sized as the 0-D model sizes it (`zpinch.circuit`). The domain's voltage is the power crossing its outer boundary over the current, so the bank, the circuit and the plasma share one energy account. | `circuit` |

**Conduction is super-time-stepped.**

* **Why.** Done explicitly, it made the flat run stall at 96 ns: once the first shock reflected
  off the axis, the hot, thin gas in the 3 µm cells there cut the step 30×.
* **How.** Once per cycle, after the hydro step, the whole mesh is gathered onto one grid. The
  heat flow is advanced over the next step in RKL2 stages (Meyer, Balsara & Aslam 2014), with
  conductivities, field direction and flux limits held at their start-of-step values.
* **Hand-back.** The result goes back to each cell as a heating rate over the next step.
* **Cost.** Before stagnation it takes 2–3 stages.

**Radiation is carried along r.**

* **How.** Once per cycle, each radial column is swept outward and inward in 17 photon groups.
  * Each cell emits its optically thin power with a bremsstrahlung spectrum.
  * Its opacity in each group follows from Kirchhoff's law, so it falls as E⁻³ above the
    cell's temperature.
  * A cell of optical depth d lets out (½ − E₃(d))/2d of its own emission through each face,
    and passes 2E₃(d) of a beam crossing it. The rest it absorbs.
  * The inward beam crosses the axis and goes out through the column again.
* **What it gives.**
  * Opaque matter radiates as a black body from its surface only.
  * The compressed shell absorbs the hot spot's soft X-rays.
  * Energy is conserved exactly. The rays' spreading in a cylinder is left out.
* **Why.** An earlier version let every cell lose up to a 0.1 mm-deep black body's power. The
  opaque liner then radiated from its whole volume, ~100× too much at stagnation, which was
  half of all the energy delivered.

**Restarts.** Runs write restart files every 5 ns.

* `run` resumes a case from its newest one.
* `python3 tools/athena_liner.py job "<case>"` prints a stand-alone script that runs or resumes
  a case on another machine.

History outputs, for the whole 1 cm liner:

* fusion power, and its burn-weighted T_i, T_e, B_z and radius;
* the burn in cells where beryllium at least doubles the bremsstrahlung;
* yield, alpha energy made and deposited, energy radiated;
* current, load voltage, energy into the domain, and energy left in the bank.

## Checks

* **Beryllium Hugoniot.** It matches the measured one to 2–6% from 0.4 to
  11 Mbar. The sound speed at rest is 8.4 km/s.
* **D-T ice Hugoniot.** It compresses 3× at 0.5 Mbar and 3.5× at 1.7 Mbar.
  Liquid deuterium reaches ~4.3× at 1 Mbar. The model is stiffer, because it
  has no molecular dissociation.
* **Ideal MHD unchanged.**
  * The physics build, with everything off and a γ = 5/3 gas, reproduces the
    ideal-MHD flat run to 6 digits through 60 ns.
  * The ideal-MHD build of the new `liner.cpp` reproduces it bit for bit.
* **Energy.** With every piece on, the change in the domain's energy plus
  what it radiated equals the energy that crossed its boundary, to 4–5
  digits (2×10⁻⁵ at 40 ns).
* **Radiation transport.**
  * A thin column lets out 99.4% of what it emits.
  * An opaque slab radiates σT⁴ from its surface, with no net loss inside.
  * Over the cells of a column, the net losses add up exactly to what escapes.
  * At the 130 ns state of an earlier run, the hot fuel loses 3.4 TW and the
    ice absorbs it. The liner's net loss is 0.3 TW, against 32 TW with the
    old black-body cap.
* **Conduction.** To 40 ns the super-time-stepped solve agrees with the
  explicit one:
  * the gas's thermal energy to 10⁻⁵;
  * the current to 10⁻⁴;
  * the total radiated energy to 1%;
  * the electrons' energy to 7%, from splitting the step.
* **Restarts.** A run resumed from its 5 ns restart file matches the
  uninterrupted run at 12 ns to 10⁻⁶ (the history's printed digits).
* **Rates.**
  * Reactivity: 1.14×10⁻²² m³/s at 10 keV.
  * Bremsstrahlung: as burn1d's 5.35×10⁻³⁷ n_e Σ n Z² √T.
  * Electron collision time: as the NRL formulary.
* **Circuit.** The current rises as the 0-D model's: 12.6 MA at 23 ns,
  against 12.1 MA for its quarter sine.

## Limits

* **Resistivity is explicit.**
  * Where diffusion would cross a cell in less than a time step, it is
    slowed to that.
  * This mostly affects the vacuum: the drive field fills it in nanoseconds
    rather than at once.
  * Early in the rise it can also cap the liner's resistivity by up to ~2×.
* **The vacuum is a thin plasma** (0.1 mg/cm³), as in the ideal-MHD runs. The
  field pushes it inward, and during the rise the field at the liner runs
  10–20% below the vacuum value. The ideal-MHD runs share this.
* **Alpha transport is averaged** over z and φ and treated as diffusion.
  That is fine for a flat liner, but it smears the alphas of a rippled one.
* **Transport coefficients** are Braginskii's for Z = 1. The Righi–Leduc,
  Ettingshausen and Biermann terms are left out, as is the thermal force
  along the field.
* **Radiation moves only along r,** in a slab-like way, and is split from the
  rest of the step.
  * Opacities come from the emission (Kirchhoff): free–free, bound–free and
    the K-shell lines are lumped together.
  * Ionisation is LTE.
  * The line model covers only beryllium's K shell.
* **The equation of state is analytic,** not a SESAME table. The beryllium
  has no strength, helium ash counts as D-T, and D:T stays 50:50.
* **The geometry is periodic along z,** so there are no end losses of heat,
  field or alphas through the electrodes. The 0-D model includes end losses.
* **Conduction is split from the rest** and held fixed over each step. It
  is first order in time.
* **Cost.** A physics run costs ~4× an ideal-MHD run per cell and cycle.
