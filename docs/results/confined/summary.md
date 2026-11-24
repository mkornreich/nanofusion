# Holding the crystal ball together: tamper, magnetic field, self-heating

Produced by `tools/confined_ball.py`. Model: `nanofusion/improved/confinement.py` (0-D, two temperatures: fusion with burn-up, delayed alpha heating, electron-ion exchange, bremsstrahlung, Braginskii heat loss to the wall, and a gold wall that recedes at the shock piston speed).

## Check: the bare ball of section 3

With no tamper, no field and no self-heating, the new time-dependent model must give the thermal burn of `energy.array_burn` (fixed temperature, burn time R / 4c_s).

| ball radius | temperature | section 3 model: thermal burn / energy in | this model |
|---|---|---|---|
| 0.1 mm | 11 keV | 0.0235 | 0.0235 |
| 0.1 mm | 31 keV | 0.0238 | 0.0238 |
| 1 mm | 11 keV | 0.235 | 0.235 |
| 1 mm | 31 keV | 0.238 | 0.237 |
| 3 mm | 11 keV | 0.706 | 0.7 |
| 3 mm | 31 keV | 0.715 | 0.707 |

## 1. What each trick buys

Solid-density D-T, ions and electrons starting at the same temperature. For each kind of ball, the cheapest one (least energy in the fuel, over starting temperature and radius) that returns 1, 10 or 100 times that energy as fusion. Energy in = heat of the ions and electrons plus the field energy inside the fuel. It leaves out every loss on the way in.

| ball | fusion = energy in | fusion = 10 × energy in | fusion = 100 × energy in |
|---|---|---|---|
| bare ball, no self-heating (the model of section 3) | 65 MJ (radius 3.7 mm, 14 keV) | 130 GJ (radius 41 mm, 20 keV) | never |
| bare ball with alpha self-heating | 60 MJ (radius 3.6 mm, 14 keV) | 16 GJ (radius 29 mm, 7 keV) | 2.7 TJ (radius 180 mm, 5 keV) |
| thick gold tamper, no field | 23 MJ (radius 2.1 mm, 28 keV) | 2.5 GJ (radius 14 mm, 10 keV) | 8.2 GJ (radius 23 mm, 7 keV) |
| gold tamper + 100 T | 2.6 MJ (radius 1.1 mm, 20 keV) | 120 MJ (radius 5.1 mm, 10 keV) | 480 MJ (radius 9.1 mm, 7 keV) |
| gold tamper + 1 kT | 160 kJ (radius 0.4 mm, 28 keV) | 8.4 MJ (radius 1.9 mm, 14 keV) | 47 MJ (radius 4.2 mm, 7 keV) |
| gold tamper + 3 kT | 52 kJ (radius 0.3 mm, 20 keV) | 4.4 MJ (radius 1.5 mm, 14 keV) | 44 MJ (radius 4.1 mm, 7 keV) |
| gold tamper + 10 kT | 49 kJ (radius 0.26 mm, 28 keV) | 5.3 MJ (radius 1.5 mm, 14 keV) | 91 MJ (radius 4.8 mm, 7 keV) |
| gold tamper + 30 kT | 280 kJ (radius 0.38 mm, 40 keV) | 27 MJ (radius 2 mm, 20 keV) | 65 GJ (radius 32 mm, 5 keV) |

* **The tamper and the field together cut breakeven about 1,200-fold**, from 60 MJ to 49 kJ. Neither does it alone: a tamped ball with no field loses its heat to the gold, and a field cannot hold fuel whose pressure is a hundred times its own.
* **The best field is 3–10 kT.** Below 1 kT too much heat still leaks and the alphas escape. Above 10 kT the field itself costs as much energy as the heat.
* **High gain stays expensive.** The wall recedes in about a nanosecond, and at solid density an alpha needs about that long to hand over its energy, so small balls never run away (ignite). A hundredfold return needs tens of megajoules in a centimetre-sized ball inside centimetres of gold.
* **Where the model stops.** Each burn is followed until the ball's radius has grown fivefold; the wall model is not trusted further. Balls that ignite are still burning then, so the hundredfold column understates them (followed ten times further, they return 20–50% more).

**Without tritium.** Solid deuterium alone, counted generously (every triton and helium-3 it makes burns at once): breakeven needs 1.4 TJ (radius 83 mm, 20 keV) bare, and 620 MJ (radius 5 mm, 40 keV) with the tamper and 3 kT. The tricks help just as much, but deuterium burns about a hundred times more slowly, so it stays out of reach. Lithium deuteride is worse: it radiates more than it fuses (section 3), and holding it longer does not change that.

## 2. How the cheapest breakeven ball burns

Gold tamper, 3 kT, radius 0.3 mm, 52 kJ in (22 µg of D-T, B·R = 0.91 T m).

| time | ion temperature | electron temperature | radius / starting radius | fusion so far / energy in |
|---|---|---|---|---|
| 0 ns | 20 keV | 20 keV | 1 | 0 |
| 0.1 ns | 18.6 keV | 18.7 keV | 1.03 | 0.13 |
| 0.3 ns | 16.5 keV | 16.7 keV | 1.09 | 0.32 |
| 1 ns | 12.3 keV | 12.9 keV | 1.26 | 0.67 |
| 2 ns | 9.32 keV | 10.2 keV | 1.43 | 0.85 |
| 4 ns | 6.54 keV | 7.43 keV | 1.66 | 0.96 |
| 8 ns | 4.25 keV | 4.84 keV | 1.95 | 1 |
| 16 ns | 2.53 keV | 2.78 keV | 2.29 | 1 |

The fuel doubles its volume after 1 ns. Of the energy put in, 16% is conducted to the wall and 3% radiated before the burn ends; the rest pushes the gold. The shock in the gold travels 0.95 mm in that time, which is how thick the tamper must be. 0.7% of the fuel burns.

## 3. Making it less ideal

The same search with a 3 kT field and one thing changed at a time.

| change | fusion = energy in | fusion = 10 × energy in |
|---|---|---|
| as in section 1 | 52 kJ (radius 0.3 mm, 20 keV) | 4.4 MJ (radius 1.5 mm, 14 keV) |
| ions heated first, electrons start at 1 keV (what exploding crystals would do) | 9.9 kJ (radius 0.17 mm, ions 40 keV) | 3 MJ (radius 1.3 mm, ions 28 keV) |
| Bohm heat loss across the field instead of classical | 52 kJ (radius 0.3 mm, 20 keV) | 4.8 MJ (radius 1.5 mm, 14 keV) |
| wall recedes 1.5 × faster than the shock formula says | 140 kJ (radius 0.43 mm, 20 keV) | 14 MJ (radius 2.2 mm, 14 keV) |
| gold shell only 0.1 mm thick | 110 kJ (radius 0.4 mm, 20 keV) | 3.1 GJ (radius 13 mm, 14 keV) |
| gold shell 0.3 mm thick | 55 kJ (radius 0.31 mm, 20 keV) | 150 MJ (radius 4.8 mm, 14 keV) |
| entrance hole, 1% of the wall | 80 kJ (radius 0.35 mm, 20 keV) | 9.2 MJ (radius 1.9 mm, 14 keV) |
| entrance hole, 2% of the wall | 120 kJ (radius 0.41 mm, 20 keV) | 18 MJ (radius 2.4 mm, 14 keV) |
| entrance hole, 5% of the wall | 340 kJ (radius 0.57 mm, 20 keV) | 54 MJ (radius 3.9 mm, 10 keV) |
| gold mixed into the fuel, 1 atom per 10,000 fuel ions | 52 kJ (radius 0.3 mm, 20 keV) | 4.8 MJ (radius 1.5 mm, 14 keV) |
| gold mixed into the fuel, 1 atom per 1,000 | 54 kJ (radius 0.28 mm, 28 keV) | 6.8 MJ (radius 1.5 mm, 20 keV) |
| fuel at 20% of solid density (a foam that light can enter) | 85 kJ (radius 0.55 mm, 28 keV) | 8.2 MJ (radius 3.1 mm, 14 keV) |
| 20% foam, ions heated first, 1% hole | 50 kJ (radius 0.51 mm, ions 40 keV) | 27 MJ (radius 4.6 mm, ions 28 keV) |
| spin-polarised fuel (cross section × 1.5, untested) | 18 kJ (radius 0.19 mm, 28 keV) | 1.8 MJ (radius 1.1 mm, 14 keV) |

* **The entrance hole matters most.** The gold wall recedes 13 times more slowly than the fuel's sound speed, so fuel streaming out of even a 2% hole is a comparable loss.
* **Heating the ions first helps**, because the electrons take about a nanosecond to catch up and meanwhile the ions burn at twice the temperature for the same energy.
* **Bohm diffusion and a little gold mix cost almost nothing at breakeven.** The wall's motion sets the clock there, not heat loss. Both matter more for the larger, longer-burning balls.
* **The wall speed is the number everything hangs on.** The energy needed goes as its cube.

## 4. A straight field instead of a ring field

Everything above assumes the field lies along the wall everywhere. Only a ring-shaped field around an axis does that, and making 3 kT at 0.3 mm that way takes about 4.5 MA through the fuel. Laser-driven coils make a straight field instead. Heat and alphas then run along it to the two end walls at the unmagnetised rate, so the fuel must be a long column, as in MagLIF. The cheapest breakeven column for each field, over its radius, its length (2–32 diameters) and its starting temperature:

| field | energy in at breakeven | column | lost to the walls | with ions heated first |
|---|---|---|---|---|
| 1 kT | 1.3 MJ | 0.89 mm wide, 3.5 mm long, 28 keV | 73% | 260 kJ |
| 3 kT | 680 kJ | 0.56 mm wide, 4.5 mm long, 28 keV | 60% | 110 kJ |
| 10 kT | 810 kJ | 0.58 mm wide, 4.7 mm long, 28 keV | 55% | 130 kJ |

A straight field costs about 13 times more energy than the ideal ring field at 3 kT, because the ends leak.

## 5. Thinner fuel: from solid crystals to cluster jets

The same tamped, magnetised ball with the fuel spread thinner, down to the average density of the gas-jet cluster experiments of section 6. For each density the field (30 T – 10 kT) that gives the cheapest breakeven. The last column is the laser power needed just to reach the stripping intensity of the crystals over the ball's cross-section.

| fuel density (× solid) | field | energy in at breakeven | ball radius | time to double its volume | gold needed | laser power to strip the crystals |
|---|---|---|---|---|---|---|
| 1 | 10 kT | 51 kJ | 0.27 mm | 820 ps | 2 mm | 12 PW |
| 0.2 | 3 kT | 85 kJ | 0.55 mm | 3.6 ns | 4 mm | 53 PW |
| 0.03 | 1 kT | 170 kJ | 1.3 mm | 28 ns | 12 mm | 290 PW |
| 0.01 | 1 kT | 240 kJ | 2.1 mm | 69 ns | 62 mm | 760 PW |
| 0.001 | 0.3 kT | 290 kJ | 4.8 mm | 730 ns | 310 mm | 4 EW |

* **Density hardly changes the energy**, unlike the bare ball (where it goes as 1/density²). A tamped, insulated ball needs roughly energy ∝ density^-0.25. That is the central point of magnetised target fusion.
* **Thin fuel buys time, not power.** The burn slows to microseconds, but the ball grows, and exploding crystals need the same intensity everywhere. Solid density needs the least laser power.

## 6. What it would take

If 15% of the laser light ended up as heat in the fuel. That is about what nanowire arrays manage: 4% of the light reaches fast deuterons, which by their share of the charge means about 15% reaches all the ions (docs/results/further/summary.md, section 5).

| design | crystals | stripping intensity | laser energy | lowest power that strips the crystals, and the pulse length at that power | fusion / laser energy |
|---|---|---|---|---|---|
| solid D-T, ions heated first, breakeven | 40 nm | 4×10¹⁸ W/cm² | 66 kJ | 3.8 PW, 17 ps | 0.15 |
| 20% foam, ions heated first, 1% hole, breakeven | 40 nm | 4×10¹⁸ W/cm² | 340 kJ | 32 PW, 10 ps | 0.15 |
| 20% foam, ions heated first, 1% hole, tenfold return | 34 nm | 3×10¹⁸ W/cm² | 180 MJ | 1.9 EW, 97 ps | 1.5 |
| solid D-T column in a straight 3 kT field, ions heated first, breakeven | 40 nm | 4×10¹⁸ W/cm² | 750 kJ | 3.7 PW, 200 ps | 0.15 |

For comparison, the most energetic picosecond-class laser shot so far is 4.4 kJ in 30 ps (NIF's ARC; LFEX gives 2 kJ in 1 ps), and the highest peak power is 10 PW, in 0.25 kJ (ELI-NP). The strongest laser-driven coil field claimed is about 800 T, and a 2022 reanalysis puts such coils well below that.

## 7. Checks against MagLIF and against published estimates

**MagLIF** on Sandia's Z machine reaches this kind of plasma by compression: deuterium near solid density at 2.5–3 keV, in a column about 0.1 mm wide and 6 mm long, held about 2 ns by its liner, with field × radius ≈ 0.34 T m. D-D neutrons from this project's reactivity, n²⟨σv⟩/4 × volume × burn time, for the published stagnation conditions:

| shot | conditions | neutrons, this model (range over the temperature error) | measured |
|---|---|---|---|
| 2014 | 2.5 ± 0.8 keV, 0.4 g/cm³, 0.047 mm³, 2 ns | 5.3×10¹² (1.0×10¹² to 1.5×10¹³) | 2.0×10¹² |
| 2020 | 2.6 ± 0.3 keV, 0.7 g/cm³, 0.05 mm³, 2 ns | 2.2×10¹³ (1.3×10¹³ to 3.3×10¹³) | 1.1×10¹³ |

The best shot is equivalent to 2 kJ of D-T fusion from fuel holding 14 kJ of heat: a ratio of 0.15. The designs above are the same kind of plasma at 20–40 keV instead of 3 keV, which nobody knows how to reach without compressing.

**Tamped fuel with no field**, as others have estimated it:

* Ruhl & Korn (arXiv:2302.06562) propose fuel with D-T at 0.5 g/cm³ in a 1 mm × 1 mm cylinder inside a high-Z case, heated to 20 keV, and find that about 3.4 MJ makes fusion exceed the fuel's energy. This model: 3.6 MJ in that fuel returns 0.77 of it.
* Their critics (Lackner et al., arXiv:2312.13429) find that perfectly tamped D-T at 0.5 g/cm³ ignites with 250 MJ. This model, with a gold wall rather than a perfect one: 490 MJ for a tenfold return.
* Their critics (Lackner et al., arXiv:2312.13429) find that perfectly tamped D-T at 0.3 g/cm³ ignites with 700 MJ. This model, with a gold wall rather than a perfect one: 1.2 GJ for a tenfold return.

So without a field this model agrees with both the proposal and its critics. The critics exempt schemes with an applied magnetic field, and that is the case worked out here.

