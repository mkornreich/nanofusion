# Shaped crystals as a directional ion source for fast ignition

Ignitor: 10 kJ of deuterons within ±25% of 4.4 MeV into a 20 µm spot; target gain = 30 MJ ÷ (250 kJ compression + ignitor laser). Laser: 6 fs pulses just strong enough to strip each crystal. Produced by `tools/shaped_crystals.py`.

## Where the deuterons go

Shares of deuterons (count) and of the whole explosion energy carried by deuterons, within a cone around the +z axis (disks: their axis; dot targets: the dot side). An isotropic source puts (1 − cos θ)/2 of its ions in a cone of half-angle θ: 0.2% at 5°, 0.8% at 10°, 3% at 20°.

| crystal | ions | D energy spread (10–90%, keV) | D within 5° | 10° | 20° | explosion energy in D within 10° | 20° | energy check (keV) |
|---|---|---|---|---|---|---|---|---|
| sphere-R2.5 | 7809 | 1.97–8.29 | 0.002 | 0.007 | 0.032 | 0.0016 | 0.0079 | 0.01 of 8.2e+04 |
| disk-a3.4-t1.7 | 7893 | 1.63–6.9 | 0.003 | 0.010 | 0.040 | 0.0026 | 0.0108 | 0.054 of 7.34e+04 |
| disk-a5-t1 | 9425 | 1.27–6.2 | 0.002 | 0.009 | 0.036 | 0.0027 | 0.0107 | 0.13 of 7.84e+04 |
| disk-a7-t0.6 | 11115 | 1.24–5.59 | 0.003 | 0.014 | 0.050 | 0.0045 | 0.0163 | 0.2 of 8.02e+04 |
| disk-a10-t0.2 | 7533 | 0.275–3.28 | 0.000 | 0.000 | 0.000 | 0.0000 | 0.0000 | 0.045 of 2.47e+04 |
| dot-a6-t1-d1.5 | 6999 | 6.23–7.03 | 0.119 | 0.390 | 1.000 | 0.0057 | 0.0152 | 0.19 of 7.86e+04 |
| dot-a6-t1-d2.5 | 7303 | 5.98–7 | 0.060 | 0.210 | 0.742 | 0.0079 | 0.0284 | 0.19 of 8.12e+04 |
| dot-a6-t1-d4 | 8023 | 5.87–6.97 | 0.024 | 0.107 | 0.391 | 0.0093 | 0.0341 | 0.2 of 8.72e+04 |
| dot-a6-t1-d6 | 9547 | 6.02–6.91 | 0.014 | 0.053 | 0.190 | 0.0095 | 0.0340 | 0.22 of 9.84e+04 |
| dot-a6-t0.6-d2.5 | 4559 | 3.19–4.31 | 0.044 | 0.160 | 0.601 | 0.0092 | 0.0349 | 0.082 of 3.05e+04 |
| dot-a6-t2-d2.5 | 12753 | 11.3–12 | 0.077 | 0.285 | 0.859 | 0.0064 | 0.0194 | 0.45 of 2.48e+05 |

## As a fast-ignition source

Each crystal is scaled to the size that needs the least laser energy for the ignitor at that distance. Isotropic reference with the textbook solid-angle fraction: 30 µm → 0.026, 100 µm → 0.0025, 300 µm → 0.00028, 1000 µm → 2.5e-05.

| crystal | source distance | useful share of explosion energy | laser → explosion | ignitor laser | peak power | target gain | crystal size at that scale (radius × thickness) | stripping intensity | thin enough for the laser to strip? |
|---|---|---|---|---|---|---|---|---|---|
| sphere-R2.5 | 30 µm | 0.0048 | 0.17 | 1.2e+07 J | 2.1e+21 W | 2.4 | 65 × 130 nm | 1e+21 W/cm² | no |
| sphere-R2.5 | 100 µm | 0.00039 (from the 10° cone) | 0.17 | 1.5e+08 J | 2.6e+22 W | 0.2 | 65 × 130 nm | 1e+21 W/cm² | no |
| sphere-R2.5 | 300 µm | 4.3e-05 (from the 10° cone) | 0.17 | 1.4e+09 J | 2.3e+23 W | 0.022 | 65 × 130 nm | 1e+21 W/cm² | no |
| sphere-R2.5 | 1000 µm | 3.9e-06 (from the 10° cone) | 0.17 | 1.5e+10 J | 2.5e+24 W | 0.002 | 65 × 130 nm | 1e+21 W/cm² | no |
| disk-a3.4-t1.7 | 30 µm | 0.0063 | 0.13 | 1.2e+07 J | 2e+21 W | 2.5 | 96 × 52 nm | 9e+20 W/cm² | no |
| disk-a3.4-t1.7 | 100 µm | 0.00063 (from the 10° cone) | 0.13 | 1.2e+08 J | 2e+22 W | 0.25 | 96 × 52 nm | 9e+20 W/cm² | no |
| disk-a3.4-t1.7 | 300 µm | 7.1e-05 (from the 10° cone) | 0.13 | 1.1e+09 J | 1.8e+23 W | 0.028 | 96 × 52 nm | 9e+20 W/cm² | no |
| disk-a3.4-t1.7 | 1000 µm | 6.4e-06 (from the 10° cone) | 0.13 | 1.2e+10 J | 2e+24 W | 0.0025 | 96 × 52 nm | 9e+20 W/cm² | no |
| disk-a5-t1 | 30 µm | 0.0099 | 0.12 | 8.4e+06 J | 1.4e+21 W | 3.5 | 154 × 31 nm | 6e+20 W/cm² | no |
| disk-a5-t1 | 100 µm | 0.00093 | 0.12 | 8.8e+07 J | 1.5e+22 W | 0.34 | 155 × 32 nm | 6e+20 W/cm² | no |
| disk-a5-t1 | 300 µm | 9.7e-05 (from the 10° cone) | 0.12 | 8.4e+08 J | 1.4e+23 W | 0.036 | 155 × 32 nm | 6e+20 W/cm² | no |
| disk-a5-t1 | 1000 µm | 8.8e-06 (from the 10° cone) | 0.12 | 9.4e+09 J | 1.6e+24 W | 0.0032 | 155 × 32 nm | 6e+20 W/cm² | no |
| disk-a7-t0.6 | 30 µm | 0.014 | 0.13 | 5.5e+06 J | 9.2e+20 W | 5.2 | 236 × 21 nm | 4e+20 W/cm² | no |
| disk-a7-t0.6 | 100 µm | 0.0014 | 0.13 | 5.5e+07 J | 9.2e+21 W | 0.54 | 239 × 21 nm | 4e+20 W/cm² | no |
| disk-a7-t0.6 | 300 µm | 0.00014 (from the 5° cone) | 0.13 | 5.7e+08 J | 9.5e+22 W | 0.053 | 239 × 21 nm | 4e+20 W/cm² | no |
| disk-a7-t0.6 | 1000 µm | 1.2e-05 (from the 5° cone) | 0.13 | 6.3e+09 J | 1.1e+24 W | 0.0048 | 239 × 21 nm | 4e+20 W/cm² | no |
| disk-a10-t0.2 | 30 µm | none: no deuterons go this way | – | – | – | – | – | – | – |
| disk-a10-t0.2 | 100 µm | none: no deuterons go this way | – | – | – | – | – | – | – |
| disk-a10-t0.2 | 300 µm | none: no deuterons go this way | – | – | – | – | – | – | – |
| disk-a10-t0.2 | 1000 µm | none: no deuterons go this way | – | – | – | – | – | – | – |
| dot-a6-t1-d1.5 | 30 µm | 0.015 | 0.12 | 5.3e+06 J | 8.9e+20 W | 5.4 | 166 × 34 nm | 3e+20 W/cm² | no |
| dot-a6-t1-d1.5 | 100 µm | 0.0021 | 0.13 | 3.8e+07 J | 6.4e+21 W | 0.78 | 170 × 35 nm | 3e+20 W/cm² | no |
| dot-a6-t1-d1.5 | 300 µm | 0.00025 (from the 5° cone) | 0.13 | 3.1e+08 J | 5.2e+22 W | 0.097 | 170 × 35 nm | 3e+20 W/cm² | no |
| dot-a6-t1-d1.5 | 1000 µm | 2.3e-05 (from the 5° cone) | 0.13 | 3.4e+09 J | 5.7e+23 W | 0.0087 | 170 × 35 nm | 3e+20 W/cm² | no |
| dot-a6-t1-d2.5 | 30 µm | 0.025 | 0.13 | 3.2e+06 J | 5.3e+20 W | 8.7 | 168 × 34 nm | 3e+20 W/cm² | no |
| dot-a6-t1-d2.5 | 100 µm | 0.0028 | 0.13 | 2.7e+07 J | 4.5e+21 W | 1.1 | 170 × 35 nm | 3e+20 W/cm² | no |
| dot-a6-t1-d2.5 | 300 µm | 0.00033 (from the 5° cone) | 0.13 | 2.3e+08 J | 3.9e+22 W | 0.13 | 170 × 35 nm | 3e+20 W/cm² | no |
| dot-a6-t1-d2.5 | 1000 µm | 3e-05 (from the 5° cone) | 0.13 | 2.6e+09 J | 4.3e+23 W | 0.012 | 170 × 35 nm | 3e+20 W/cm² | no |
| dot-a6-t1-d4 | 30 µm | 0.029 | 0.13 | 2.6e+06 J | 4.3e+20 W | 11 | 170 × 35 nm | 3e+20 W/cm² | no |
| dot-a6-t1-d4 | 100 µm | 0.0033 | 0.13 | 2.3e+07 J | 3.8e+21 W | 1.3 | 170 × 35 nm | 3e+20 W/cm² | no |
| dot-a6-t1-d4 | 300 µm | 0.00031 (from the 5° cone) | 0.13 | 2.5e+08 J | 4.1e+22 W | 0.12 | 170 × 35 nm | 3e+20 W/cm² | no |
| dot-a6-t1-d4 | 1000 µm | 2.8e-05 (from the 5° cone) | 0.13 | 2.7e+09 J | 4.5e+23 W | 0.011 | 170 × 35 nm | 3e+20 W/cm² | no |
| dot-a6-t1-d6 | 30 µm | 0.03 | 0.13 | 2.6e+06 J | 4.3e+20 W | 11 | 168 × 34 nm | 4e+20 W/cm² | no |
| dot-a6-t1-d6 | 100 µm | 0.003 | 0.13 | 2.6e+07 J | 4.3e+21 W | 1.2 | 168 × 34 nm | 4e+20 W/cm² | no |
| dot-a6-t1-d6 | 300 µm | 0.00035 (from the 5° cone) | 0.13 | 2.2e+08 J | 3.6e+22 W | 0.14 | 168 × 34 nm | 4e+20 W/cm² | no |
| dot-a6-t1-d6 | 1000 µm | 3.1e-05 (from the 5° cone) | 0.13 | 2.4e+09 J | 4e+23 W | 0.012 | 168 × 34 nm | 4e+20 W/cm² | no |
| dot-a6-t0.6-d2.5 | 30 µm | 0.031 | 0.14 | 2.3e+06 J | 3.9e+20 W | 12 | 217 × 29 nm | 2e+20 W/cm² | no |
| dot-a6-t0.6-d2.5 | 100 µm | 0.0029 | 0.14 | 2.4e+07 J | 4e+21 W | 1.2 | 222 × 30 nm | 2e+20 W/cm² | no |
| dot-a6-t0.6-d2.5 | 300 µm | 0.00035 (from the 5° cone) | 0.14 | 2e+08 J | 3.3e+22 W | 0.15 | 222 × 30 nm | 2e+20 W/cm² | no |
| dot-a6-t0.6-d2.5 | 1000 µm | 3.1e-05 (from the 5° cone) | 0.14 | 2.2e+09 J | 3.7e+23 W | 0.013 | 222 × 30 nm | 2e+20 W/cm² | no |
| dot-a6-t2-d2.5 | 30 µm | 0.017 | 0.13 | 4.7e+06 J | 7.8e+20 W | 6.1 | 127 × 43 nm | 4e+20 W/cm² | no |
| dot-a6-t2-d2.5 | 100 µm | 0.0021 | 0.13 | 3.8e+07 J | 6.3e+21 W | 0.79 | 129 × 44 nm | 4e+20 W/cm² | no |
| dot-a6-t2-d2.5 | 300 µm | 0.00025 (from the 5° cone) | 0.13 | 3.2e+08 J | 5.3e+22 W | 0.095 | 129 × 44 nm | 4e+20 W/cm² | no |
| dot-a6-t2-d2.5 | 1000 µm | 2.3e-05 (from the 5° cone) | 0.13 | 3.5e+09 J | 5.8e+23 W | 0.0086 | 129 × 44 nm | 4e+20 W/cm² | no |

## Best design at each distance

| distance | crystal | ignitor laser | peak power (today's largest lasers: ~10¹⁶ W) | target gain |
|---|---|---|---|---|
| 30 µm | dot-a6-t0.6-d2.5 | 2.3e+06 J | 3.9e+20 W (4e+04× today's largest) | 12 |
| 100 µm | dot-a6-t1-d4 | 2.3e+07 J | 3.8e+21 W (4e+05× today's largest) | 1.3 |
| 300 µm | dot-a6-t0.6-d2.5 | 2e+08 J | 3.3e+22 W (3e+06× today's largest) | 0.15 |
| 1000 µm | dot-a6-t0.6-d2.5 | 2.2e+09 J | 3.7e+23 W (4e+07× today's largest) | 0.013 |

