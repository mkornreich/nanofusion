# Hollow-shell Coulomb implosions

Hot-electron density 1e+28 m⁻³ throughout. Produced by `tools/implosion_experiments.py`.

## A. Electron-to-ion charge ratio κ (D₂, cavity radius 8 nm)

| κ | ions | energy at centre, no space charge (keV) | mean ion energy reached (keV) | focus radius / R (run) | smooth-shell bounce / R | column per ion (m⁻²) | column / space-charge bound | D+D in flight | D+D when the ions spring back into the wall | gain vs ion energy | gain vs field energy |
|---|---|---|---|---|---|---|---|---|---|---|---|
| 2 | 10724 | 1.93 | 0.717 | 0.359 | 0.366 | 6e+19 | 1.09 | 2e-13 | 6.8e-11 | 3.2e-11 | 3e-11 |
| 5 | 4290 | 1.93 | 1.16 | 0.171 | 0.171 | 1.6e+20 | 0.83 | 2.3e-13 | 9.2e-12 | 7e-12 | 4.2e-12 |
| 10 | 2145 | 1.93 | 1.41 | 0.097 | 0.0916 | 3e+20 | 0.74 | 2.2e-13 | 9.2e-13 | 1.5e-12 | 5.5e-13 |
| 30 | 715 | 1.93 | 1.68 | 0.0419 | 0.0323 | 7.6e+20 | 0.68 | 3.3e-13 | 1.3e-14 | 1.4e-12 | 2.1e-13 |
| 100 | 215 | 1.93 | 1.84 | 0.0165 | 0.00993 | 1.4e+21 | 0.46 | 4.7e-13 | 1.2e-16 | 3.5e-12 | 1.7e-13 |
| 300 | 72 | 1.93 | 1.91 | 0.00822 | 0.00335 | 2.7e+21 | 0.41 | 2e-13 | 2e-18 | 1.7e-12 | 2.9e-14 |
| 1000 | 22 | 1.93 | 1.96 | 0.00484 | 0.00102 | too few ions to sample | | | | | |

## B. Size scaling, κ = 10 (D₂)

| R (nm) | ions | mean ion energy (keV) | focus radius / R | D+D in flight (run) | D+D predicted from the next smaller run |
|---|---|---|---|---|---|
| 4 | 269 | 0.371 | 0.104 | 1e-21 | – |
| 6 | 905 | 0.808 | 0.0998 | 1.8e-16 | 3.3e-16 |
| 9 | 3054 | 1.78 | 0.0951 | 3.9e-12 | 1.9e-12 |
| 13 | 9203 | 3.67 | 0.0945 | 2.4e-09 | 1.7e-09 |

Extrapolated from the R = 13 nm run (same κ and hot-electron density):

| R (nm) | ions | energy at centre (keV) | D+D in flight | fusion energy from ions springing back into the wall | gain vs ion energy | gain vs field energy | field energy (J) |
|---|---|---|---|---|---|---|---|
| 20 | 3.4e+04 | 12.1 | 4.8e-07 | 95% | 1.3e-07 | 4.6e-08 | 1.3e-10 |
| 30 | 1.1e+05 | 27.1 | 2.2e-05 | 97% | 1.1e-06 | 3.9e-07 | 9.8e-10 |
| 50 | 5.2e+05 | 75.4 | 0.00099 | 98% | 6.4e-06 | 2.3e-06 | 1.3e-08 |
| 75 | 1.8e+06 | 170 | 0.012 | 99% | 1.8e-05 | 6.6e-06 | 9.6e-08 |
| 100 | 4.2e+06 | 302 | 0.059 | 99% | 3.5e-05 | 1.3e-05 | 4e-07 |
| 150 | 1.4e+07 | 679 | 0.47 | 100% | 8.1e-05 | 2.9e-05 | 3.1e-06 |
| 200 | 3.4e+07 | 1.21e+03 | 1.8 | 100% | 0.00014 | 5e-05 | 1.3e-05 |
| 300 | 1.1e+08 | 2.71e+03 | 11 | 100% | 0.00027 | 9.6e-05 | 9.8e-05 |

## B. Size scaling, κ = 100 (D₂)

| R (nm) | ions | mean ion energy (keV) | focus radius / R | D+D in flight (run) | D+D predicted from the next smaller run |
|---|---|---|---|---|---|
| 6 | 91 | 1.05 | 0.0198 | 9.6e-17 | – |
| 9 | 306 | 2.32 | 0.0166 | 3.4e-12 | 1e-12 |
| 13 | 921 | 4.81 | 0.0155 | 1.5e-09 | 1.2e-09 |
| 18 | 2443 | 9.16 | 0.0147 | 1.1e-07 | 1e-07 |

Extrapolated from the R = 18 nm run (same κ and hot-electron density):

| R (nm) | ions | energy at centre (keV) | D+D in flight | fusion energy from ions springing back into the wall | gain vs ion energy | gain vs field energy | field energy (J) |
|---|---|---|---|---|---|---|---|
| 20 | 3.4e+03 | 12.1 | 3.4e-07 | 18% | 4e-08 | 1.9e-09 | 1.3e-10 |
| 30 | 1.1e+04 | 27.1 | 1.6e-05 | 44% | 3.5e-07 | 1.7e-08 | 9.8e-10 |
| 50 | 5.2e+04 | 75.4 | 0.00071 | 72% | 2.5e-06 | 1.2e-07 | 1.3e-08 |
| 75 | 1.8e+05 | 170 | 0.0088 | 85% | 7.6e-06 | 3.6e-07 | 9.6e-08 |
| 100 | 4.2e+05 | 302 | 0.043 | 91% | 1.5e-05 | 7e-07 | 4e-07 |
| 150 | 1.4e+06 | 679 | 0.34 | 96% | 3.5e-05 | 1.7e-06 | 3.1e-06 |
| 200 | 3.4e+06 | 1.21e+03 | 1.3 | 98% | 6.3e-05 | 2.9e-06 | 1.3e-05 |
| 300 | 1.1e+07 | 2.71e+03 | 7.7 | 99% | 0.00013 | 6.2e-06 | 9.8e-05 |

## C. Initial ion temperature (D₂, R = 8 nm, κ = 100)

| temperature (eV) | focus radius / R | D+D in flight | gain vs ion energy |
|---|---|---|---|
| 0 | 0.0165 | 4.7e-13 | 3.5e-12 |
| 0.01 | 0.0171 | 4.1e-13 | 3.7e-12 |
| 0.1 | 0.0196 | 2.7e-13 | 2.8e-12 |
| 1 | 0.0322 | 6.2e-13 | 7.9e-12 |
| 10 | 0.0842 | 1.9e-13 | 1.2e-12 |

## D. Fuel (R = 8 nm, κ = 30)

| fuel | ions | D focus time (fs) | partner focus time (fs) | column of D seen by D (m⁻²) | column of partner seen by D (m⁻²) | reactions in flight | gain vs ion energy |
|---|---|---|---|---|---|---|---|
| D2 | 715 | 30.0 | 30.0 | 7.6e+20 | 7.6e+20 | 4.5e-13 | 1.4e-12 |
| DT | 715 | 29.8 | 37.5 | 1.1e+21 | 1.9e+19 | 2.5e-11 | 3.7e-10 |
| 6LiD | 354 | 30.3 | 30.0 | 1.1e+20 | 1.3e+20 | 2.3e-15 | 8e-15 |

Fuels extrapolated from their R = 8 nm, κ = 30 runs:

| fuel | R (nm) | energy at centre, D (keV) | reactions in flight | share from ions springing back into the wall | gain vs ion energy | gain vs field energy |
|---|---|---|---|---|---|---|
| D2 | 20 | 12.1 | 3.2e-07 | 63% | 2.7e-08 | 3.9e-09 |
| D2 | 30 | 27.1 | 1.6e-05 | 81% | 3.4e-07 | 5e-08 |
| D2 | 50 | 75.4 | 0.00074 | 92% | 2.9e-06 | 4.2e-07 |
| D2 | 100 | 302 | 0.047 | 98% | 1.9e-05 | 2.7e-06 |
| D2 | 200 | 1.21e+03 | 1.5 | 99% | 8.1e-05 | 1.2e-05 |
| DT | 20 | 12.1 | 5.9e-06 | 61% | 2.2e-06 | 3.2e-07 |
| DT | 30 | 27.1 | 0.00023 | 88% | 3.6e-05 | 5.2e-06 |
| DT | 50 | 75.4 | 0.0045 | 99% | 0.00047 | 6.8e-05 |
| DT | 100 | 302 | 0.062 | 100% | 0.0014 | 0.0002 |
| DT | 200 | 1.21e+03 | 1.3 | 100% | 0.00099 | 0.00014 |
| 6LiD | 20 | 12.1 | 5.4e-09 | 82% | 9.3e-10 | 1.3e-10 |
| 6LiD | 30 | 27.1 | 3.3e-07 | 91% | 1.6e-08 | 2.3e-09 |
| 6LiD | 50 | 75.4 | 2.2e-05 | 93% | 1.7e-07 | 2.5e-08 |
| 6LiD | 100 | 302 | 0.0043 | 97% | 5.4e-06 | 7.8e-07 |
| 6LiD | 200 | 1.21e+03 | 0.36 | 99% | 9e-05 | 1.3e-05 |
| DT, arriving together | 20 | 12.1 | 9.9e-06 | 75% | 5.8e-06 | 8.4e-07 |
| DT, arriving together | 30 | 27.1 | 0.00074 | 85% | 9.9e-05 | 1.4e-05 |
| DT, arriving together | 50 | 75.4 | 0.036 | 95% | 0.001 | 0.00015 |
| DT, arriving together | 100 | 302 | 0.72 | 98% | 0.0017 | 0.00025 |
| DT, arriving together | 200 | 1.21e+03 | 5.6 | 99% | 0.001 | 0.00015 |

Best gain against the field energy, for cavities small enough that no ion exceeds 4.9 MeV (the end of the measured cross sections):

| fuel | best R (nm) | gain vs field energy | gain vs ion energy | of which from the focus | share from ions springing back into the wall | fusion reactions per implosion |
|---|---|---|---|---|---|---|
| D2 | 273 | 2.1e-05 | 0.00014 | 3.9e-07 | 100% | 5.9 in flight |
| DT | 105 | 0.0002 | 0.0014 | 1.3e-06 | 100% | 0.075 in flight |
| 6LiD | 177 | 8.7e-06 | 6e-05 | 5.5e-07 | 99% | 0.18 in flight |
| DT, arriving together | 82 | 0.00026 | 0.0018 | 4.6e-05 | 97% | 0.35 in flight |

## Space-charge bound

Only so many ions fit in a focus of radius r_c before their own repulsion stops the rest: N e²/(4πε₀ r_c) ≤ E. So the matter an ion crosses is at most 4ε₀E/(e² r_c), and the fusion energy per unit of ion energy is at most

* D₂: 2ε₀σQ/(e² r_c) = 7.84 fm / r_c (σ at its maximum, 0.19 b);
* D–T (50/50, arriving together): ε₀σQ/(e² r_c) = 0.492 pm / r_c (σ = 5.1 b at 64 keV).

Breaking even on the ions' own energy needs a focus smaller than these lengths; an atom is about 50,000 fm (50 pm) across.

## E. Repeated passes (D₂, R = 8 nm, κ = 100, two oscillation periods)

| pass | time (fs) | focus radius / R | D+D in this pass |
|---|---|---|---|
| 1 | 29.7 | 0.0165 | 3.8e-13 |
| 2 | 88.9 | 0.0335 | 2.9e-12 |
| 3 | 148.5 | 0.0443 | 3.3e-14 |
| 4 | 208.3 | 0.0616 | 4.9e-13 |

## F. Numerical check: reaction-volume radius (D₂, R = 8 nm, κ = 100)

| radius (nm) | D+D in flight (exact estimator) | D+D from the exposure histogram |
|---|---|---|
| 0.01 | 5.9e-13 | 5e-13 |
| 0.0198 | 4.7e-13 | 3.8e-13 |
| 0.05 | 4e-13 | 5.2e-13 |

