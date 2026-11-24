/*
 * All-pairs Coulomb forces for the improved NanoFusion model, plus in-crystal
 * fusion-rate bookkeeping. Compiled on first use (see native.py); OpenMP is
 * optional. SI units throughout.
 *
 * For every ordered pair the force on i is k*q_i*q_j*(x_i-x_j)/(r^2+s^2)^(3/2)
 * (Plummer softening s, tiny compared with the lattice spacing).
 *
 * For unordered pairs closer than r_c the kernel also adds the binary-collision
 * reaction rate sigma(E)*|v_rel|/V_c, i.e. the probability per second that the
 * pair fuses if the partner were spread uniformly over the volume V_c around
 * the projectile. E is the centre-of-mass energy, kinetic + Coulomb potential,
 * in keV. Cross sections: Bosch & Hale (1992) for D+D, sigma = S/E *
 * exp(-B_G/sqrt(E)) with constant S for D+Li.
 */
#include <math.h>

#define K_COULOMB 8.9875517873681764e9
#define J_PER_KEV 1.602176634e-16
#define MB_TO_M2 1e-31
#define E_CAP_KEV 4900.0  /* cross sections are held constant above this CM energy */
#define NBINS 64          /* pair-exposure histogram: 8 bins per decade, 1e-4 .. 1e4 keV */

/* Bosch-Hale S-factors (keV mb), D(d,p)T and D(d,n)3He; B coefficients are 0. */
static double sigma_dd_mb(double e_kev)
{
    if (e_kev < 1e-3) return 0.0;
    const double bg = 31.3970;
    double e = e_kev < E_CAP_KEV ? e_kev : E_CAP_KEV;
    double s_p = 5.5576e4 + e * (2.1054e2 + e * (-3.2638e-2 + e * (1.4987e-6 + e * 1.8181e-10)));
    double s_n = 5.3701e4 + e * (3.3027e2 + e * (-1.2706e-1 + e * (2.9327e-5 + e * -2.5151e-9)));
    return (s_p + s_n) / (e * exp(bg / sqrt(e)));
}

static double sigma_dli_mb(double e_kev, double s_kev_mb, double bg)
{
    if (e_kev < 1e-3) return 0.0;
    double e = e_kev < E_CAP_KEV ? e_kev : E_CAP_KEV;
    return s_kev_mb / (e * exp(bg / sqrt(e)));
}

/*
 * n          number of ions
 * x, v       positions (m) and velocities (m/s), n*3 each
 * q, mass    charge (C) and mass (kg) per ion
 * kind       0 = deuterium, 1 = lithium
 * soft2      softening length squared (m^2)
 * rc2        reaction-volume radius squared (m^2); inv_vc = 1/V_c (m^-3)
 * s_dli, bg_dli   D+Li S-factor (keV mb) and Gamow constant (keV^0.5)
 * force      out: n*3 Coulomb forces (N)
 * stats      out[0] D+Li rate (1/s), [1] D+D rate (1/s), [2] min D-Li distance (m),
 *            [3] max D-Li CM energy within r_c (keV), [4] min r/|v_rel| over approaching
 *            pairs within r_c (s), [5] Coulomb potential energy (J), [6] min distance over all pairs (m)
 */
void coulomb_step(int n, const double *x, const double *v, const double *q, const double *mass,
                  const int *kind, double soft2, double rc2, double inv_vc,
                  double s_dli, double bg_dli, double *force, double *stats,
                  double *hist_dd, double *hist_dli)
{
    /* hist_*: for close pairs, sum of |v_rel| / V_c binned by relative kinetic
       energy (1/2 mu v_rel^2). Times dt and a cross section, this gives the
       expected reactions in that energy range. */
    for (int b = 0; b < NBINS; b++) { hist_dd[b] = 0.0; hist_dli[b] = 0.0; }
    double rate_dli = 0.0, rate_dd = 0.0, min_dli = 1e300, max_e_dli = 0.0, min_tau = 1e300, pot = 0.0;
    double min_r2 = 1e300;

#pragma omp parallel for schedule(dynamic, 16) reduction(+:rate_dli,rate_dd,pot) \
    reduction(min:min_dli,min_tau,min_r2) reduction(max:max_e_dli) \
    reduction(+:hist_dd[0:NBINS], hist_dli[0:NBINS])
    for (int i = 0; i < n; i++) {
        const double xi = x[3*i], yi = x[3*i+1], zi = x[3*i+2];
        const double qi = q[i];
        double fx = 0.0, fy = 0.0, fz = 0.0;
        for (int j = 0; j < n; j++) {
            if (j == i) continue;
            const double dx = xi - x[3*j], dy = yi - x[3*j+1], dz = zi - x[3*j+2];
            const double r2 = dx*dx + dy*dy + dz*dz;
            const double rs2 = r2 + soft2;
            const double inv = 1.0 / sqrt(rs2);
            const double c = K_COULOMB * qi * q[j] * inv * inv * inv;
            fx += c * dx; fy += c * dy; fz += c * dz;
            if (j < i) continue;  /* pair quantities once per unordered pair */
            pot += K_COULOMB * qi * q[j] * inv;
            if (r2 < min_r2) min_r2 = r2;
            if (r2 > rc2) continue;
            const int dli = kind[i] != kind[j];
            const double r = sqrt(r2);
            if (dli && r < min_dli) min_dli = r;
            const double dvx = v[3*i] - v[3*j], dvy = v[3*i+1] - v[3*j+1], dvz = v[3*i+2] - v[3*j+2];
            const double v2 = dvx*dvx + dvy*dvy + dvz*dvz;
            const double vrel = sqrt(v2);
            if (vrel > 0.0 && dx*dvx + dy*dvy + dz*dvz < 0.0) {
                const double tau = r / vrel;
                if (tau < min_tau) min_tau = tau;
            }
            if (kind[i] == 1 && kind[j] == 1) continue;  /* Li+Li: negligible (Z=3 x 3) */
            const double mu = mass[i] * mass[j] / (mass[i] + mass[j]);
            const double ek_kev = 0.5 * mu * v2 / J_PER_KEV;
            const double e_kev = ek_kev + K_COULOMB * qi * q[j] / (r > 0 ? r : 1e-30) / J_PER_KEV;
            if (ek_kev > 1e-4) {
                int b = (int)((log10(ek_kev) + 4.0) * 8.0);
                if (b >= NBINS) b = NBINS - 1;
                if (dli) hist_dli[b] += vrel * inv_vc; else hist_dd[b] += vrel * inv_vc;
            }
            if (dli) {
                if (e_kev > max_e_dli) max_e_dli = e_kev;
                rate_dli += sigma_dli_mb(e_kev, s_dli, bg_dli) * MB_TO_M2 * vrel * inv_vc;
            } else {
                rate_dd += sigma_dd_mb(e_kev) * MB_TO_M2 * vrel * inv_vc;
            }
        }
        force[3*i] = fx; force[3*i+1] = fy; force[3*i+2] = fz;
    }
    stats[0] = rate_dli; stats[1] = rate_dd; stats[2] = min_dli;
    stats[3] = max_e_dli; stats[4] = min_tau; stats[5] = pot; stats[6] = sqrt(min_r2);
}
