//========================================================================================
// Athena++ problem generator: a magnetised liner implosion in r-z, for nanofusion.
//
// A beryllium liner (outer radius r_out, inner r_in) lined with D-T ice down to r_gas and
// filled with D-T gas, in an axial field bz. The current enters as B_phi at the outer
// radial boundary, rising as sin^2 to i_peak at t_rise and then held. The liner's outer
// surface carries a multi-mode ripple (rms amplitude `ripple`), the seed for the
// magneto-Rayleigh-Taylor instability: axial modes n = 1..modes, and for 3-D runs also
// azimuthal and helical modes m = 1..mmax, n = -modes..modes, all of equal amplitude
// with random phases. mseed >= 0 keeps only the modes with m <= mseed (mseed = 0: the
// axisymmetric part of the same seed, for a 2-D run to compare with the 3-D one).
//
// With two cells in phi over the full circle and no azimuthal modes the solution stays
// axisymmetric, i.e. an r-z simulation (Athena++ has no x1-x3 2-D). With more cells in
// phi it is a full 3-D simulation.
//
// Units: length mm, time ns, density g/cm3. Then velocity 1e6 m/s, pressure 1e15 Pa, and
// B_code = B_SI / 35450 T (Athena's B^2/2 magnetic pressure).
//
// Passive scalars tag the material: r0 liner, r1 ice, r2 gas. History outputs give the
// gas and ice masses and their mass-weighted <R^2>, for the convergence.
//
// Two builds:
// * ideal MHD (the original): --prob=liner --coord=cylindrical -b --nscalars=3 --flux=hlld -omp
// * with the physics a yield needs: --prob=liner --coord=cylindrical -b --nscalars=5
//   --flux=hlld --eos=general/liner_eos -omp, on Athena++ with tools/athena/athena-physics.patch
//   and liner_eos.cpp / liner_physics.hpp copied in (see tools/athena_liner.py). It adds
//   two more scalars, r3 helium ash and r4 the electrons' energy (per unit mass), and,
//   each switched on in <problem>:
//     burn        D-T burn (Bosch-Hale, ion temperature), fuel depletion into ash; the
//                 neutrons leave
//     alpha       alpha heating: alphas slow down where they are made (0) or spread by a
//                 radial diffusion model across the field (1, averaged over z and phi),
//                 shared between ions and electrons
//     radiation   bremsstrahlung, recombination and beryllium K-shell lines, carried
//                 along r in 17 photon groups with Kirchhoff opacities (lp::ColumnTransport):
//                 opaque matter radiates from its surface, the shell absorbs the hot spot's
//                 soft X-rays
//     resistivity Spitzer-Braginskii (alpha_perp), with the Lee-More minimum mean free
//                 path in cold dense matter, capped at eta_max; the vacuum gets eta_vac
//                 so the drive field fills it within ~ns. Explicit: where diffusion would
//                 cross a cell in less than a time step it is slowed to that (needs eta_ohm = 1)
//     conduction  electron heat flux, Braginskii parallel and perpendicular, flux-limited,
//                 super-time-stepped (RKL2) once per cycle on one grid
//     two_temp    separate ion and electron temperatures: shocks heat ions; ohmic heat,
//                 conduction, radiation and most alpha heat act on the electrons
//     nernst      the field carried down the electron temperature gradient (needs eta_ohm)
//     circuit     the current from a capacitor bank behind driver_nh, sized as the 0-D
//                 model sizes it (zpinch.circuit), instead of the prescribed sin^2 rise
//   The equation of state (liner_eos.cpp) has cold curves, ionisation and degenerate
//   electrons for beryllium and D-T. History outputs add the fusion power and its
//   burn-weighted temperatures, field and radius, and the energy budget.
//========================================================================================

#include <algorithm>
#include <cmath>
#include <limits>
#include <random>
#include <sstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "../athena.hpp"
#include "../athena_arrays.hpp"
#include "../bvals/bvals.hpp"
#include "../coordinates/coordinates.hpp"
#include "../eos/eos.hpp"
#include "../field/field.hpp"
#include "../hydro/hydro.hpp"
#include "../mesh/mesh.hpp"
#include "../parameter_input.hpp"
#include "../scalars/scalars.hpp"

#if !MAGNETIC_FIELDS_ENABLED
#error "liner.cpp requires magnetic fields (-b)"
#endif
#if NSCALARS < 3
#error "liner.cpp needs --nscalars=3 (liner, ice, gas)"
#endif

#define LINER_PHYSICS (GENERAL_EOS && NSCALARS >= 5)
#if LINER_PHYSICS
#include "../eos/general/liner_physics.hpp"
#include "../field/field_diffusion/field_diffusion.hpp"
#endif

namespace {
Real r_out, r_in, r_gas, rho_liner, rho_ice, rho_gas, rho_vac, p0, bz0;
Real i_peak, t_rise;            // MA, ns
struct Mode { int m, n; Real amp, phase; };
std::vector<Mode> seed_modes;
Real z_len;

// The prescribed current's field at radius r (mm) and time t (ns), in code units: mu0 I / (2 pi r).
Real BPhi(Real r, Real t) {
  Real s = (t < t_rise) ? std::sin(0.5*M_PI*t/t_rise) : 1.0;
  return i_peak*1.0e6*s*s*5.642e-9/r;
}

// The liner's outer surface at (z, phi): r_out plus the ripple.
Real Surface(Real z, Real phi) {
  Real d = 0.0;
  for (const Mode &md : seed_modes)
    d += md.amp*std::cos(2.0*M_PI*md.n*z/z_len + md.m*phi + md.phase);
  return r_out + d;
}

#if LINER_PHYSICS
// ---- physics switches and parameters ---------------------------------------------------
bool burn_on, radiation_on, resistivity_on, conduction_on, two_temp_on, nernst_on, circuit_on;
int alpha_mode;          // -1 off, 0 local, 1 radial diffusion
Real eta_vac, eta_max, flux_limit, cond_safety;
Real eps_e_vac;          // the vacuum's electron energy (code units) at p0
// the circuit (SI): bank, machine inductance, inductance between the domain and the return wall
Real cap_c, l_machine, l_outside, length_m, r_domain_m;
// Mesh-level state kept in ruser_mesh_data[0] (saved in restart files)
enum MeshState {kCurrent, kVcap, kFluxPrev, kCurrentPrev, kYield, kAlphaMade, kAlphaDep, kRadiated,
                kDelivered, kBank0, kNState};
// per-cell plasma state cached in ruser_meshblock_data[0]
enum Plasma {kTe, kTi, kNe, kNi, kZeff, kTauE, kXe, kEta, kKN, kKpar, kKperp, kFvac, kZbe, kZl,
             kNfuel, kNbe, kNhe, kNlight, kAlphaTau, kAlphaFi, kAlphaD, kRate, kRad, kEpsE, kEpsI, kEpsTh, kCondRate, kRadNet,
             kNPlasma};
// the alphas' energy density (J/m3) per global radial cell, from the last cycle
std::vector<Real> alpha_eps;
// the domain's real length (m) per simulated length: the history and the circuit use the real liner
Real length_scale;
Real load_voltage = 0;   // V, over the last cycle
bool physics_source;     // some physics needs the source function (and keeps the plasma state fresh)

Real Clamp01(Real x) { return std::min(static_cast<Real>(1.0), std::max(static_cast<Real>(0.0), x)); }

// Plasma state of cells [kl..ku][jl..ju][il..iu] from the primitives; t_guess reuses the
// last electron temperature as the starting point.
void PlasmaState(MeshBlock *pmb, const AthenaArray<Real> &w, const AthenaArray<Real> &r,
                 const AthenaArray<Real> &bcc, int il, int iu, int jl, int ju, int kl, int ku) {
  AthenaArray<Real> &c = pmb->ruser_meshblock_data[0];
  bool valid = pmb->ruser_meshblock_data[1](0) > 0;
  // one block around the whole periodic circle (the r-z runs): its phi ghosts are copies
  Mesh *pm = pmb->pmy_mesh;
  bool full_phi = pm->f2 && pmb->block_size.nx2 == pm->mesh_size.nx2
                  && pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic;
  int jl0 = jl, ju0 = ju;
  if (full_phi) { jl = std::max(jl, pmb->js); ju = std::min(ju, pmb->je); }
  for (int k = kl; k <= ku; ++k)
    for (int j = jl; j <= ju; ++j)
      for (int i = il; i <= iu; ++i) {
        Real rho = w(IDN, k, j, i), p = w(IPR, k, j, i);
        Real f0 = Clamp01(r(0, k, j, i)), f1 = Clamp01(r(1, k, j, i)), f2 = Clamp01(r(2, k, j, i));
        Real f3 = Clamp01(r(3, k, j, i));
        Real fvac = std::max(static_cast<Real>(0.0), 1.0 - f0 - f1 - f2 - f3);
        double rho_si = rho * lp::U_RHO, rho_be = f0 * rho, rho_l = rho - rho_be;
        double n_be = f0 * rho_si / (lp::A_BE * lp::AMU);
        double n_fuel = (f1 + f2) * rho_si / (lp::A_DT * lp::AMU);
        double n_vac = fvac * rho_si / (lp::A_DT * lp::AMU);
        double n_he = f3 * rho_si / (lp::A_HE * lp::AMU);
        double n_i_eos = n_be + rho_l * lp::U_RHO / (lp::A_DT * lp::AMU);  // the EOS's ion count
        double e;
        lp::FromRhoP(rho, f0, p * lp::U_P, &e, nullptr);
        double eps_th = std::max(rho_si * (e - lp::ColdEnergy(rho, f0)), 0.0);
        double te, ti, zb, zl, eps_e;
        if (two_temp_on) {
          eps_e = std::max(r(4, k, j, i) * rho * lp::U_P, 1e-30);
          double guess = valid ? std::log10(std::max(c(kTe, k, j, i), 1e-3))
                               : std::log10(lp::FromRhoE(rho, f0, e).t_ev);
          te = lp::ElectronTemperature(rho_be, rho_l, eps_e, guess, &zb, &zl);
          double eps_i = eps_th - eps_e;
          ti = (eps_i > 0) ? eps_i / (1.5 * n_i_eos * lp::QE) : te;
          c(kEpsI, k, j, i) = std::max(eps_i, 0.0);
          c(kEpsTh, k, j, i) = eps_th;
        } else {
          te = ti = lp::FromRhoE(rho, f0, e).t_ev;
          eps_e = lp::ElectronEnergy(rho_be, rho_l, te);
          lp::Fwd fb = lp::Forward(lp::BE, std::max(rho_be, 1e-30), std::log10(te));
          lp::Fwd fl = lp::Forward(lp::LIGHT, std::max(rho_l, 1e-30), std::log10(te));
          zb = fb.z; zl = fl.z;
          c(kEpsI, k, j, i) = std::max(eps_th - eps_e, 0.0);
          c(kEpsTh, k, j, i) = eps_th;
        }
        double ne = std::max(zb * n_be + zl * (n_fuel + n_vac + 2 * n_he), 1.0);
        double nz2 = zb * zb * n_be + zl * zl * (n_fuel + n_vac + 4 * n_he);
        double ni = n_be + n_fuel + n_vac + n_he;
        double zeff = nz2 / ne;
        double sza = (zb * zb * n_be / lp::A_BE + zl * zl * (n_fuel + n_vac) * lp::DT_Z2_OVER_A
                      + 4 * zl * zl * n_he / lp::A_HE) / ne;
        double tau = lp::TauE(ne, ni, te, zeff);
        double b_t = std::sqrt(SQR(bcc(IB1, k, j, i)) + SQR(bcc(IB2, k, j, i)) + SQR(bcc(IB3, k, j, i)))
                     * lp::U_B;
        double x = lp::QE * b_t / lp::ME * tau;
        double eta = lp::AlphaPerp(x) * lp::ME / (ne * lp::QE * lp::QE * tau);
        eta = std::min((1 - fvac) * eta + fvac * eta_vac, static_cast<double>(eta_max));
        double plasma = 1 - fvac;   // conduction and the Nernst effect only in real material
        double kbase = ne * te * lp::QE * tau / lp::ME;
        lp::Alpha al{0.0, 1.0, 0.0};
        double d_alpha = 0;
        if (burn_on && alpha_mode >= 0) {
          al = lp::AlphaSlowing(ne, te, sza);
          double r_larmor = lp::ALPHA_LARMOR_TM / std::max(b_t, 1e-6);
          d_alpha = al.range * al.range / (3 * al.tau_e) / (1 + SQR(al.range / r_larmor));
        }
        double rate = 0.25 * n_fuel * n_fuel * lp::SigmaVDT(ti * 1e-3);
        double rad = 0;
        if (radiation_on) {
          rad = lp::ThinRadiation(ne, te, n_be, zb, n_fuel + n_vac, zl, n_he) * plasma;
        }
        c(kTe, k, j, i) = te; c(kTi, k, j, i) = ti; c(kNe, k, j, i) = ne; c(kNi, k, j, i) = ni;
        c(kZeff, k, j, i) = zeff; c(kTauE, k, j, i) = tau; c(kXe, k, j, i) = x;
        c(kEta, k, j, i) = eta; c(kKN, k, j, i) = tau / lp::ME * lp::BetaWedgeOverX(x) * plasma;
        c(kKpar, k, j, i) = lp::KAPPA_PAR * kbase * plasma;
        c(kKperp, k, j, i) = lp::KappaPerp(x) * kbase * plasma;
        c(kFvac, k, j, i) = fvac; c(kZbe, k, j, i) = zb; c(kZl, k, j, i) = zl;
        c(kNfuel, k, j, i) = n_fuel; c(kNbe, k, j, i) = n_be; c(kNhe, k, j, i) = n_he;
        c(kNlight, k, j, i) = n_fuel + n_vac;
        c(kAlphaTau, k, j, i) = al.tau_e; c(kAlphaFi, k, j, i) = al.f_ions; c(kAlphaD, k, j, i) = d_alpha;
        c(kRate, k, j, i) = burn_on ? rate : 0.0; c(kRad, k, j, i) = rad; c(kEpsE, k, j, i) = eps_e;
      }
  if (full_phi && (jl0 < jl || ju0 > ju)) {
    for (int n = 0; n < kNPlasma; ++n)
      for (int k = kl; k <= ku; ++k)
        for (int i = il; i <= iu; ++i) {
          if (jl0 < jl) c(n, k, jl0, i) = c(n, k, pmb->je, i);
          if (ju0 > ju) c(n, k, ju0, i) = c(n, k, pmb->js, i);
        }
  }
  pmb->ruser_meshblock_data[1](0) = 1.0;
}

// The cells whose state the fluxes and edges need: the block and one ghost layer.
void PlasmaRange(MeshBlock *pmb, int *il, int *iu, int *jl, int *ju, int *kl, int *ku) {
  *il = pmb->is - 1; *iu = pmb->ie + 1;
  *jl = pmb->js; *ju = pmb->je; *kl = pmb->ks; *ku = pmb->ke;
  if (pmb->pmy_mesh->f2) { *jl -= 1; *ju += 1; }
  if (pmb->pmy_mesh->f3) { *kl -= 1; *ku += 1; }
}
// The callbacks use the state the source function left (half a step old); without the
// source function, or before its first call, compute it now.
void EnsurePlasma(MeshBlock *pmb) {
  if (physics_source && pmb->ruser_meshblock_data[1](0) > 0) return;
  int il, iu, jl, ju, kl, ku;
  PlasmaRange(pmb, &il, &iu, &jl, &ju, &kl, &ku);
  PlasmaState(pmb, pmb->phydro->w, pmb->pscalars->r, pmb->pfield->bcc, il, iu, jl, ju, kl, ku);
}

// cell-to-cell distances (code units) along each direction
inline Real Dx1(Coordinates *pc, int i) { return pc->x1v(i) - pc->x1v(i - 1); }
inline Real Dx2(Coordinates *pc, int i, int j) { return pc->x1v(i) * (pc->x2v(j) - pc->x2v(j - 1)); }
inline Real Dx3(Coordinates *pc, int k) { return pc->x3v(k) - pc->x3v(k - 1); }

// monotonised-central average of four one-sided differences (Sharma & Hammett 2007)
inline Real MC4(Real a, Real b, Real c, Real d) {
  auto mc = [](Real x, Real y) {
    if (x * y <= 0) return static_cast<Real>(0.0);
    Real s = (x > 0) ? 1.0 : -1.0;
    return s * std::min(std::min(2 * std::abs(x), 2 * std::abs(y)), 0.5 * std::abs(x + y));
  };
  return mc(mc(a, b), mc(c, d));
}

// Electron temperature gradient (eV per code length) at cell (k,j,i), central.
void GradTe(MeshBlock *pmb, int k, int j, int i, Real g[3]) {
  AthenaArray<Real> &c = pmb->ruser_meshblock_data[0];
  Coordinates *pc = pmb->pcoord;
  g[0] = (c(kTe, k, j, i + 1) - c(kTe, k, j, i - 1)) / (pc->x1v(i + 1) - pc->x1v(i - 1));
  g[1] = pmb->pmy_mesh->f2 ? (c(kTe, k, j + 1, i) - c(kTe, k, j - 1, i))
                             / (pc->x1v(i) * (pc->x2v(j + 1) - pc->x2v(j - 1))) : 0.0;
  g[2] = pmb->pmy_mesh->f3 ? (c(kTe, k + 1, j, i) - c(kTe, k - 1, j, i))
                             / (pc->x3v(k + 1) - pc->x3v(k - 1)) : 0.0;
}

// ---- resistivity and the Nernst effect -------------------------------------------------
void LinerDiffusivity(FieldDiffusion *pfdif, MeshBlock *pmb, const AthenaArray<Real> &w,
                      const AthenaArray<Real> &bmag, int is, int ie, int js, int je, int ks, int ke) {
  EnsurePlasma(pmb);
  AthenaArray<Real> &c = pmb->ruser_meshblock_data[0];
  int il, iu, jl, ju, kl, ku;
  PlasmaRange(pmb, &il, &iu, &jl, &ju, &kl, &ku);
  Coordinates *pc = pmb->pcoord;
  Mesh *pm = pmb->pmy_mesh;
  Real fac = pm->f3 ? 1.0 / 6.0 : (pm->f2 ? 0.25 : 0.5);   // as FieldDiffusion::NewDiffusionDt
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie; ++i) {
        int kk = std::min(std::max(k, kl), ku), jj = std::min(std::max(j, jl), ju);
        int ii = std::min(std::max(i, il), iu);
        Real eta = resistivity_on ? c(kEta, kk, jj, ii) / (lp::MU0 * lp::U_DIFF) : 0.0;
        // cap: diffusion faster than ~a cell per time step would set the (explicit) step;
        // it only slows diffusion over a few cells (the vacuum still fills in ~ns)
        Real len = pc->dx1f(ii);
        if (pm->f2) len = std::min(len, pc->x1v(ii) * pc->dx2f(jj));
        if (pm->f3) len = std::min(len, pc->dx3f(kk));
        if (pm->dt > 0 && pm->dt < 1e10) eta = std::min(eta, 0.9 * fac * pm->cfl_number * len * len / pm->dt);
        if (nernst_on && ii > il && ii < iu && (!pmb->pmy_mesh->f2 || (jj > jl && jj < ju))
            && (!pmb->pmy_mesh->f3 || (kk > kl && kk < ku))) {
          // upwinding the Nernst advection: numerical diffusivity |v_N| dx / 2
          Real g[3]; GradTe(pmb, kk, jj, ii, g);
          Real gmag = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) * lp::QE / lp::U_LEN;  // J/m
          Real vn = c(kKN, kk, jj, ii) * gmag / lp::U_VEL;                                      // code
          Real dx = std::min(pc->dx1f(ii), pmb->pmy_mesh->f3 ? pc->dx3f(kk) : pc->dx1f(ii));
          eta += 0.5 * vn * dx;
        }
        pfdif->etaB(FieldDiffusion::ohmic, k, j, i) = eta;
      }
}

// Nernst EMF on the cell edges: E = K_N grad(T_e) x B (code units: K_N grad T in m/s / U_VEL).
void NernstEMF(FieldDiffusion *pfdif, MeshBlock *pmb, const AthenaArray<Real> &w,
               const AthenaArray<Real> &bcc, EdgeField &e) {
  if (!nernst_on) return;
  EnsurePlasma(pmb);
  AthenaArray<Real> &c = pmb->ruser_meshblock_data[0];
  Coordinates *pc = pmb->pcoord;
  const Real s = lp::QE / lp::U_LEN / lp::U_VEL;  // eV per code length -> J/m, then per U_VEL
  int is = pmb->is, ie = pmb->ie, js = pmb->js, je = pmb->je, ks = pmb->ks, ke = pmb->ke;
  // e1 at (x1v(i), x2f(j), x3f(k)): K (dT2 B3 - dT3 B2)
  for (int k = ks; k <= ke + 1; ++k)
    for (int j = js; j <= je + 1; ++j)
      for (int i = is; i <= ie; ++i) {
        Real kn = 0.25 * (c(kKN, k, j, i) + c(kKN, k, j - 1, i) + c(kKN, k - 1, j, i) + c(kKN, k - 1, j - 1, i));
        Real t2 = 0.5 * ((c(kTe, k, j, i) - c(kTe, k, j - 1, i)) + (c(kTe, k - 1, j, i) - c(kTe, k - 1, j - 1, i)))
                  / Dx2(pc, i, j);
        Real t3 = 0.5 * ((c(kTe, k, j, i) - c(kTe, k - 1, j, i)) + (c(kTe, k, j - 1, i) - c(kTe, k - 1, j - 1, i)))
                  / Dx3(pc, k);
        Real b2 = 0.25 * (bcc(IB2, k, j, i) + bcc(IB2, k, j - 1, i) + bcc(IB2, k - 1, j, i) + bcc(IB2, k - 1, j - 1, i));
        Real b3 = 0.25 * (bcc(IB3, k, j, i) + bcc(IB3, k, j - 1, i) + bcc(IB3, k - 1, j, i) + bcc(IB3, k - 1, j - 1, i));
        e.x1e(k, j, i) += kn * s * (t2 * b3 - t3 * b2);
      }
  // e2 at (x1f(i), x2v(j), x3f(k)): K (dT3 B1 - dT1 B3)
  for (int k = ks; k <= ke + 1; ++k)
    for (int j = js; j <= je; ++j)
      for (int i = is; i <= ie + 1; ++i) {
        Real kn = 0.25 * (c(kKN, k, j, i) + c(kKN, k, j, i - 1) + c(kKN, k - 1, j, i) + c(kKN, k - 1, j, i - 1));
        Real t3 = 0.5 * ((c(kTe, k, j, i) - c(kTe, k - 1, j, i)) + (c(kTe, k, j, i - 1) - c(kTe, k - 1, j, i - 1)))
                  / Dx3(pc, k);
        Real t1 = 0.5 * ((c(kTe, k, j, i) - c(kTe, k, j, i - 1)) + (c(kTe, k - 1, j, i) - c(kTe, k - 1, j, i - 1)))
                  / Dx1(pc, i);
        Real b1 = 0.25 * (bcc(IB1, k, j, i) + bcc(IB1, k, j, i - 1) + bcc(IB1, k - 1, j, i) + bcc(IB1, k - 1, j, i - 1));
        Real b3 = 0.25 * (bcc(IB3, k, j, i) + bcc(IB3, k, j, i - 1) + bcc(IB3, k - 1, j, i) + bcc(IB3, k - 1, j, i - 1));
        e.x2e(k, j, i) += kn * s * (t3 * b1 - t1 * b3);
      }
  // e3 at (x1f(i), x2f(j), x3v(k)): K (dT1 B2 - dT2 B1)
  for (int k = ks; k <= ke; ++k)
    for (int j = js; j <= je + 1; ++j)
      for (int i = is; i <= ie + 1; ++i) {
        Real kn = 0.25 * (c(kKN, k, j, i) + c(kKN, k, j, i - 1) + c(kKN, k, j - 1, i) + c(kKN, k, j - 1, i - 1));
        Real t1 = 0.5 * ((c(kTe, k, j, i) - c(kTe, k, j, i - 1)) + (c(kTe, k, j - 1, i) - c(kTe, k, j - 1, i - 1)))
                  / Dx1(pc, i);
        Real rf = pc->x1f(i);
        Real t2 = 0.5 * ((c(kTe, k, j, i) - c(kTe, k, j - 1, i)) + (c(kTe, k, j, i - 1) - c(kTe, k, j - 1, i - 1)))
                  / (rf * (pc->x2v(j) - pc->x2v(j - 1)));
        Real b1 = 0.25 * (bcc(IB1, k, j, i) + bcc(IB1, k, j, i - 1) + bcc(IB1, k, j - 1, i) + bcc(IB1, k, j - 1, i - 1));
        Real b2 = 0.25 * (bcc(IB2, k, j, i) + bcc(IB2, k, j, i - 1) + bcc(IB2, k, j - 1, i) + bcc(IB2, k, j - 1, i - 1));
        e.x3e(k, j, i) += kn * s * (t1 * b2 - t2 * b1);
      }
}

// ---- electron heat conduction: super-time-stepped on one grid --------------------------
// Explicit conduction is stiff where thin gas is hot (near the axis, after the first shock
// reflects there): its step can fall 100x below the hydro step. So once per cycle, after the
// hydro step, the whole mesh is gathered onto one grid (one process, uniform blocks) and the
// heat flow q = -kperp grad T - (kpar - kperp) b (b . grad T), flux-limited to flux_limit
// n_e T_e v_te, is advanced over the next step by RKL2 super-time-stepping (Meyer, Balsara &
// Aslam 2014): s stages, s set by the explicit limit. Conductivities, field direction, flux
// limits and heat capacity are frozen over the step. The result goes back as a heating rate
// per cell, which the source term applies over the next step (to the total energy, and to
// the electrons in the two-temperature model): conservative, and the blocks stay consistent.
struct CondGrid {
  int n1 = 0, n2 = 0, n3 = 0;
  bool per2 = true, per3 = true;
  Real dphi = 0, dz = 0;                     // rad, m
  std::vector<Real> rv, rf;                  // m
  std::vector<Real> t0, cv, kpar, kperp, qm, b1, b2, b3;      // per cell, frozen
  std::vector<Real> y0, y1, y2, ynew, m0, d, f1, f2, f3;      // per cell / face, the stages
};
CondGrid cg;
int cond_stages = 0;   // the last cycle's (history)

inline int Cell(int i, int j, int k) { return (k * cg.n2 + j) * cg.n1 + i; }
inline int Wrap(int j, int n, bool per) { return per ? (j + n) % n : std::min(std::max(j, 0), n - 1); }

// d = dT/dt (eV/s) from the heat flow, for temperatures y (eV)
void CondRate(const std::vector<Real> &y, int nthreads) {
  const int n1 = cg.n1, n2 = cg.n2, n3 = cg.n3;
  auto Y = [&](int i, int j, int k) { return y[Cell(i, j, k)]; };
  // limited heat flux (W/m2) along `dir` through a face between cells a and b (b after a)
  auto flux = [&](int a, int b, int dir, Real g1, Real g2, Real g3) {
    Real bx = cg.b1[a] + cg.b1[b], by = cg.b2[a] + cg.b2[b], bz = cg.b3[a] + cg.b3[b];
    Real bm = std::sqrt(bx * bx + by * by + bz * bz);
    Real bg = 0, bd = 0;
    if (bm > 0) {
      bx /= bm; by /= bm; bz /= bm;
      bg = bx * g1 + by * g2 + bz * g3;
      bd = (dir == 0) ? bx : (dir == 1) ? by : bz;
    }
    Real kpar = 0.5 * (cg.kpar[a] + cg.kpar[b]), kperp = 0.5 * (cg.kperp[a] + cg.kperp[b]);
    Real gd = (dir == 0) ? g1 : (dir == 1) ? g2 : g3;
    Real q = -(kperp * gd + (kpar - kperp) * bd * bg) * lp::QE;
    Real qm = std::max(cg.qm[a], cg.qm[b]);
    return q / (1 + std::abs(q) / std::max(qm, static_cast<Real>(1e-30)));
  };
  auto diff1 = [&](int i, int j, int k) {   // one-sided r differences, 0 off the grid
    Real up = (i + 1 < n1) ? Y(i + 1, j, k) - Y(i, j, k) : 0.0;
    Real dn = (i > 0) ? Y(i, j, k) - Y(i - 1, j, k) : 0.0;
    return std::make_pair(up, dn);
  };
#pragma omp parallel for num_threads(nthreads) schedule(static)
  for (int k = 0; k < n3; ++k) {
    int kp = Wrap(k + 1, n3, cg.per3), km = Wrap(k - 1, n3, cg.per3);
    for (int j = 0; j < n2; ++j) {
      int jp = Wrap(j + 1, n2, cg.per2), jm = Wrap(j - 1, n2, cg.per2);
      // r faces: i = 1..n1-1 (the axis wall and the outer edge carry no heat)
      for (int i = 0; i <= n1; ++i) {
        Real q = 0;
        if (i > 0 && i < n1) {
          Real g1 = (Y(i, j, k) - Y(i - 1, j, k)) / (cg.rv[i] - cg.rv[i - 1]);
          Real g2 = (n2 > 1) ? MC4(Y(i, jp, k) - Y(i, j, k), Y(i, j, k) - Y(i, jm, k),
                                   Y(i - 1, jp, k) - Y(i - 1, j, k), Y(i - 1, j, k) - Y(i - 1, jm, k))
                               / (cg.rf[i] * cg.dphi) : 0.0;
          Real g3 = (n3 > 1) ? MC4(Y(i, j, kp) - Y(i, j, k), Y(i, j, k) - Y(i, j, km),
                                   Y(i - 1, j, kp) - Y(i - 1, j, k), Y(i - 1, j, k) - Y(i - 1, j, km))
                               / cg.dz : 0.0;
          q = flux(Cell(i - 1, j, k), Cell(i, j, k), 0, g1, g2, g3) * cg.rf[i] * cg.dphi * cg.dz;
        }
        cg.f1[(k * n2 + j) * (n1 + 1) + i] = q;
      }
      for (int i = 0; i < n1; ++i) {
        Real dr = cg.rf[i + 1] - cg.rf[i];
        // phi face before cell j, and z face before cell k (periodic: face 0 is face n)
        Real q2 = 0, q3 = 0;
        if (n2 > 1 && (cg.per2 || j > 0)) {
          int jb = Wrap(j - 1, n2, cg.per2);
          Real g2 = (Y(i, j, k) - Y(i, jb, k)) / (cg.rv[i] * cg.dphi);
          auto a = diff1(i, j, k), bb = diff1(i, jb, k);
          Real g1 = MC4(a.first, a.second, bb.first, bb.second) / dr;
          Real g3 = (n3 > 1) ? MC4(Y(i, j, kp) - Y(i, j, k), Y(i, j, k) - Y(i, j, km),
                                   Y(i, jb, kp) - Y(i, jb, k), Y(i, jb, k) - Y(i, jb, km)) / cg.dz : 0.0;
          q2 = flux(Cell(i, jb, k), Cell(i, j, k), 1, g1, g2, g3) * dr * cg.dz;
        }
        if (n3 > 1 && (cg.per3 || k > 0)) {
          int kb = Wrap(k - 1, n3, cg.per3);
          Real g3 = (Y(i, j, k) - Y(i, j, kb)) / cg.dz;
          auto a = diff1(i, j, k), bb = diff1(i, j, kb);
          Real g1 = MC4(a.first, a.second, bb.first, bb.second) / dr;
          Real g2 = (n2 > 1) ? MC4(Y(i, jp, k) - Y(i, j, k), Y(i, j, k) - Y(i, jm, k),
                                   Y(i, jp, kb) - Y(i, j, kb), Y(i, j, kb) - Y(i, jm, kb))
                               / (cg.rv[i] * cg.dphi) : 0.0;
          q3 = flux(Cell(i, j, kb), Cell(i, j, k), 2, g1, g2, g3)
               * 0.5 * (cg.rf[i + 1] * cg.rf[i + 1] - cg.rf[i] * cg.rf[i]) * cg.dphi;
        }
        cg.f2[Cell(i, j, k)] = q2;   // through the face before cell j
        cg.f3[Cell(i, j, k)] = q3;   // through the face before cell k
      }
    }
  }
#pragma omp parallel for num_threads(nthreads) schedule(static)
  for (int k = 0; k < n3; ++k) {
    int kn = Wrap(k + 1, n3, cg.per3);
    for (int j = 0; j < n2; ++j) {
      int jn = Wrap(j + 1, n2, cg.per2);
      for (int i = 0; i < n1; ++i) {
        int c = Cell(i, j, k);
        Real out = cg.f1[(k * n2 + j) * (n1 + 1) + i + 1] - cg.f1[(k * n2 + j) * (n1 + 1) + i];
        if (n2 > 1) out += ((cg.per2 || j + 1 < n2) ? cg.f2[Cell(i, jn, k)] : 0.0) - cg.f2[c];
        if (n3 > 1) out += ((cg.per3 || k + 1 < n3) ? cg.f3[Cell(i, j, kn)] : 0.0) - cg.f3[c];
        Real vol = 0.5 * (cg.rf[i + 1] * cg.rf[i + 1] - cg.rf[i] * cg.rf[i]) * cg.dphi * cg.dz;
        cg.d[c] = -out / (vol * cg.cv[c]);
      }
    }
  }
}

// ---- radiation: carried along r, column by column (lp::ColumnTransport) ------------------
// Once per cycle, after the hydro step: each cell's net radiative loss (W/m3), which the
// source term applies over the next step. Returns the power that leaves the domain (W, for
// the simulated slice).
Real RadiationTransport(Mesh *pm) {
  const int n1 = pm->mesh_size.nx1, n2 = pm->mesh_size.nx2, n3 = pm->mesh_size.nx3;
  static std::vector<Real> vol, dr, te, eta, loss;
  const size_t n = static_cast<size_t>(n1) * n2 * n3;
  if (vol.size() != n) for (auto *v : {&vol, &dr, &te, &eta, &loss}) v->assign(n, 0.0);
  auto at = [&](int i, int j, int k) { return (static_cast<size_t>(k) * n2 + j) * n1 + i; };
  for (int b = 0; b < pm->nblocal; ++b) {
    MeshBlock *pmb = pm->my_blocks(b);
    AthenaArray<Real> &c = pmb->ruser_meshblock_data[0];
    Coordinates *pc = pmb->pcoord;
    int i0 = static_cast<int>(pmb->loc.lx1) * pmb->block_size.nx1 - pmb->is;
    int j0 = static_cast<int>(pmb->loc.lx2) * pmb->block_size.nx2 - pmb->js;
    int k0 = static_cast<int>(pmb->loc.lx3) * pmb->block_size.nx3 - pmb->ks;
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j)
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          size_t g = at(i + i0, j + j0, k + k0);
          vol[g] = pc->GetCellVolume(k, j, i) * lp::U_VOL;
          dr[g] = pc->dx1f(i) * lp::U_LEN;
          te[g] = c(kTe, k, j, i);
          eta[g] = c(kRad, k, j, i);
        }
  }
  Real escaped = 0;
  const int ncol = n2 * n3;
#pragma omp parallel num_threads(pm->GetNumMeshThreads()) reduction(+:escaped)
  {
    std::vector<double> work;
#pragma omp for schedule(dynamic)
    for (int col = 0; col < ncol; ++col) {
      size_t c0 = static_cast<size_t>(col) * n1;
      escaped += lp::ColumnTransport(n1, &vol[c0], &dr[c0], &te[c0], &eta[c0], &loss[c0], work);
    }
  }
  for (int b = 0; b < pm->nblocal; ++b) {
    MeshBlock *pmb = pm->my_blocks(b);
    AthenaArray<Real> &c = pmb->ruser_meshblock_data[0];
    int i0 = static_cast<int>(pmb->loc.lx1) * pmb->block_size.nx1 - pmb->is;
    int j0 = static_cast<int>(pmb->loc.lx2) * pmb->block_size.nx2 - pmb->js;
    int k0 = static_cast<int>(pmb->loc.lx3) * pmb->block_size.nx3 - pmb->ks;
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j)
        for (int i = pmb->is; i <= pmb->ie; ++i)
          c(kRadNet, k, j, i) = loss[at(i + i0, j + j0, k + k0)];
  }
  return escaped;
}

// Gather, super-time-step over tau (s), and leave the heating rate (W/m3) in each block's cache.
void ConductSTS(Mesh *pm, Real tau) {
  const int n1 = pm->mesh_size.nx1, n2 = pm->mesh_size.nx2, n3 = pm->mesh_size.nx3;
  const int nthreads = pm->GetNumMeshThreads();
  if (cg.n1 != n1) {   // first call: the grid
    cg.n1 = n1; cg.n2 = n2; cg.n3 = n3;
    cg.per2 = pm->mesh_bcs[BoundaryFace::inner_x2] == BoundaryFlag::periodic;
    cg.per3 = pm->mesh_bcs[BoundaryFace::inner_x3] == BoundaryFlag::periodic;
    cg.dphi = (pm->mesh_size.x2max - pm->mesh_size.x2min) / n2;
    cg.dz = (pm->mesh_size.x3max - pm->mesh_size.x3min) / n3 * lp::U_LEN;
    cg.rv.assign(n1, 0.0); cg.rf.assign(n1 + 1, 0.0);
    int n = n1 * n2 * n3;
    for (auto *v : {&cg.t0, &cg.cv, &cg.kpar, &cg.kperp, &cg.qm, &cg.b1, &cg.b2, &cg.b3, &cg.y0, &cg.y1,
                    &cg.y2, &cg.ynew, &cg.m0, &cg.d, &cg.f2, &cg.f3})
      v->assign(n, 0.0);
    cg.f1.assign((n1 + 1) * n2 * n3, 0.0);
  }
  // gather (the plasma state was just computed for the end of the step)
  Real dt_expl = std::numeric_limits<Real>::max();
  for (int b = 0; b < pm->nblocal; ++b) {
    MeshBlock *pmb = pm->my_blocks(b);
    AthenaArray<Real> &c = pmb->ruser_meshblock_data[0];
    AthenaArray<Real> &bcc = pmb->pfield->bcc;
    Coordinates *pc = pmb->pcoord;
    int i0 = static_cast<int>(pmb->loc.lx1) * pmb->block_size.nx1 - pmb->is;
    int j0 = static_cast<int>(pmb->loc.lx2) * pmb->block_size.nx2 - pmb->js;
    int k0 = static_cast<int>(pmb->loc.lx3) * pmb->block_size.nx3 - pmb->ks;
    for (int i = pmb->is; i <= pmb->ie; ++i) {
      cg.rv[i + i0] = pc->x1v(i) * lp::U_LEN;
      cg.rf[i + i0] = pc->x1f(i) * lp::U_LEN; cg.rf[i + i0 + 1] = pc->x1f(i + 1) * lp::U_LEN;
    }
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j)
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          int g = Cell(i + i0, j + j0, k + k0);
          Real ne = c(kNe, k, j, i), te = c(kTe, k, j, i), plasma = 1 - c(kFvac, k, j, i);
          cg.t0[g] = te;
          cg.cv[g] = 1.5 * (two_temp_on ? ne : ne + c(kNi, k, j, i)) * lp::QE;   // J/m3/eV
          cg.kpar[g] = c(kKpar, k, j, i); cg.kperp[g] = c(kKperp, k, j, i);
          Real tj = te * lp::QE;
          cg.qm[g] = flux_limit * ne * tj * std::sqrt(tj / lp::ME) * plasma;
          Real bx = bcc(IB1, k, j, i), by = bcc(IB2, k, j, i), bz = bcc(IB3, k, j, i);
          Real bm = std::sqrt(bx * bx + by * by + bz * bz);
          cg.b1[g] = bm > 0 ? bx / bm : 0; cg.b2[g] = bm > 0 ? by / bm : 0; cg.b3[g] = bm > 0 ? bz / bm : 0;
          // explicit limit: the conductivity the flux limiter allows, over the heat capacity
          if (plasma > 0.5) {
            Real gr[3]; GradTe(pmb, k, j, i, gr);
            Real gm = std::sqrt(gr[0] * gr[0] + gr[1] * gr[1] + gr[2] * gr[2]) * lp::QE / lp::U_LEN;
            Real kmax = std::max(cg.kpar[g], cg.kperp[g]);
            if (gm > 0) kmax = std::min(kmax, cg.qm[g] / gm);
            Real chi = kmax * lp::QE / cg.cv[g];   // m2/s
            Real inv = 1 / SQR(pc->dx1f(i) * lp::U_LEN);
            if (n2 > 1) inv += 1 / SQR(pc->x1v(i) * pc->dx2f(j) * lp::U_LEN);
            if (n3 > 1) inv += 1 / SQR(pc->dx3f(k) * lp::U_LEN);
            if (chi > 0) dt_expl = std::min(dt_expl, static_cast<Real>(1 / (2 * chi * inv)));
          }
        }
  }
  dt_expl *= cond_safety;
  // RKL2 stages: (s^2 + s - 2) / 4 >= tau / dt_expl
  int s = static_cast<int>(std::ceil(0.5 * (-1 + std::sqrt(9 + 16 * tau / dt_expl))));
  s = std::max(2, std::min(s, 2000));
  cond_stages = s;
  const Real w1 = 4.0 / (s * s + s - 2.0);
  auto bj = [](int j) { return (j <= 2) ? 1.0 / 3.0 : (j * j + j - 2.0) / (2.0 * j * (j + 1.0)); };
  const int n = n1 * n2 * n3;
  cg.y0 = cg.t0;
  CondRate(cg.y0, nthreads);
  cg.m0 = cg.d;
  for (int c = 0; c < n; ++c) cg.y1[c] = cg.y0[c] + bj(1) * w1 * tau * cg.m0[c];   // Y1
  cg.y2 = cg.y0;                                                                    // Y(j-2)
  for (int j = 2; j <= s; ++j) {
    Real mu = (2.0 * j - 1) / j * bj(j) / bj(j - 1);
    Real nu = -(j - 1.0) / j * bj(j) / bj(j - 2);
    Real mt = mu * w1, gt = -(1 - bj(j - 1)) * mt;
    CondRate(cg.y1, nthreads);
#pragma omp parallel for num_threads(nthreads) schedule(static)
    for (int c = 0; c < n; ++c)
      cg.ynew[c] = mu * cg.y1[c] + nu * cg.y2[c] + (1 - mu - nu) * cg.y0[c] + mt * tau * cg.d[c]
                   + gt * tau * cg.m0[c];
    std::swap(cg.y2, cg.y1);   // Y(j-2) <- Y(j-1)
    std::swap(cg.y1, cg.ynew); // Y(j-1) <- Y(j)
  }
  // the heating rate, back into each block's cache
  for (int b = 0; b < pm->nblocal; ++b) {
    MeshBlock *pmb = pm->my_blocks(b);
    AthenaArray<Real> &c = pmb->ruser_meshblock_data[0];
    int i0 = static_cast<int>(pmb->loc.lx1) * pmb->block_size.nx1 - pmb->is;
    int j0 = static_cast<int>(pmb->loc.lx2) * pmb->block_size.nx2 - pmb->js;
    int k0 = static_cast<int>(pmb->loc.lx3) * pmb->block_size.nx3 - pmb->ks;
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j)
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          int g = Cell(i + i0, j + j0, k + k0);
          c(kCondRate, k, j, i) = cg.cv[g] * (cg.y1[g] - cg.t0[g]) / tau;
        }
  }
}

// ---- the source terms --------------------------------------------------------------
void LinerSource(MeshBlock *pmb, const Real time, const Real dt, const AthenaArray<Real> &prim,
                 const AthenaArray<Real> &prim_scalar, const AthenaArray<Real> &bcc,
                 AthenaArray<Real> &cons, AthenaArray<Real> &cons_scalar) {
  int il, iu, jl, ju, kl, ku;
  PlasmaRange(pmb, &il, &iu, &jl, &ju, &kl, &ku);
  PlasmaState(pmb, prim, prim_scalar, bcc, il, iu, jl, ju, kl, ku);
  AthenaArray<Real> &c = pmb->ruser_meshblock_data[0];
  Coordinates *pc = pmb->pcoord;
  const bool f2 = pmb->pmy_mesh->f2, f3 = pmb->pmy_mesh->f3;
  const Real dt_s = dt * lp::U_TIME;
  int nx1g = pmb->pmy_mesh->mesh_size.nx1;
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        Real vol = pc->GetCellVolume(k, j, i);
        double ne = c(kNe, k, j, i), te = c(kTe, k, j, i), ti = c(kTi, k, j, i);
        double eps_e = c(kEpsE, k, j, i);
        double d_tot = 0, d_e = 0;   // energy changes (J/m3) of the gas and of its electrons
        // electron heat conduction: the rate the super-time-step left (ConductSTS)
        if (conduction_on) {
          double q = c(kCondRate, k, j, i) * dt_s;   // J/m3 over this stage
          q = std::max(q, -0.9 * eps_e);
          d_tot += q; d_e += q;
        }
        // radiation: the net loss the transport left (RadiationTransport); the electrons emit
        // and absorb, and lose no more than they have
        if (radiation_on) {
          double p_rad = c(kRadNet, k, j, i) * dt_s;
          double lost = (p_rad > 0) ? p_rad / (1 + p_rad / std::max(eps_e, 1e-30)) : p_rad;
          d_tot -= lost; d_e -= lost;
        }
        // alpha heating: where they are made, or from the radial diffusion model
        if (burn_on && alpha_mode >= 0) {
          double dep;
          if (alpha_mode == 0 || alpha_eps.empty()) {
            dep = c(kRate, k, j, i) * lp::E_ALPHA;
          } else {
            int ig = static_cast<int>(pmb->loc.lx1) * pmb->block_size.nx1 + (i - pmb->is);
            dep = (ig < nx1g) ? alpha_eps[ig] / c(kAlphaTau, k, j, i) : 0.0;
          }
          d_tot += dep * dt_s;
          d_e += dep * dt_s * (1 - c(kAlphaFi, k, j, i));
        }
        // burn: fuel -> ash; the neutron's mass, momentum and kinetic energy leave
        if (burn_on) {
          double n_fuel = c(kNfuel, k, j, i);
          double reactions = std::min(c(kRate, k, j, i) * dt_s, 0.45 * n_fuel);
          if (reactions > 0) {
            Real dfuel = reactions * 2 * lp::A_DT * lp::AMU / lp::U_RHO;   // code density
            Real dash = reactions * lp::A_HE * lp::AMU / lp::U_RHO;
            Real dneut = dfuel - dash;
            Real s1 = cons_scalar(1, k, j, i), s2 = cons_scalar(2, k, j, i);
            Real share = (s1 + s2 > 0) ? s1 / (s1 + s2) : 0.5;
            cons_scalar(1, k, j, i) = std::max(static_cast<Real>(0.0), s1 - dfuel * share);
            cons_scalar(2, k, j, i) = std::max(static_cast<Real>(0.0), s2 - dfuel * (1 - share));
            cons_scalar(3, k, j, i) += dash;
            Real v1 = prim(IVX, k, j, i), v2 = prim(IVY, k, j, i), v3 = prim(IVZ, k, j, i);
            cons(IDN, k, j, i) -= dneut;
            cons(IM1, k, j, i) -= dneut * v1;
            cons(IM2, k, j, i) -= dneut * v2;
            cons(IM3, k, j, i) -= dneut * v3;
            cons(IEN, k, j, i) -= 0.5 * dneut * (v1 * v1 + v2 * v2 + v3 * v3);
          }
        }
        if (two_temp_on) {
          // compression of the electrons: p_e dV (the scalar carries only advection)
          Real dv = (pc->GetFace1Area(k, j, i + 1) * 0.5 * (prim(IVX, k, j, i + 1) + prim(IVX, k, j, i))
                     - pc->GetFace1Area(k, j, i) * 0.5 * (prim(IVX, k, j, i) + prim(IVX, k, j, i - 1)));
          if (f2) dv += pc->GetFace2Area(k, j + 1, i) * 0.5 * (prim(IVY, k, j + 1, i) + prim(IVY, k, j, i))
                        - pc->GetFace2Area(k, j, i) * 0.5 * (prim(IVY, k, j, i) + prim(IVY, k, j - 1, i));
          if (f3) dv += pc->GetFace3Area(k + 1, j, i) * 0.5 * (prim(IVZ, k + 1, j, i) + prim(IVZ, k, j, i))
                        - pc->GetFace3Area(k, j, i) * 0.5 * (prim(IVZ, k, j, i) + prim(IVZ, k - 1, j, i));
          double divv = dv / vol / lp::U_TIME;   // 1/s
          double pe = std::min(ne * te * lp::QE, 2.0 / 3.0 * eps_e);
          d_e += eps_e * std::expm1(-(pe / std::max(eps_e, 1e-30)) * divv * dt_s);
          // ohmic heating (already in the total energy through the Poynting flux) and the
          // work of the Nernst field on the current go to the electrons
          double ohmic = 0;   // J/m3 this step
          if (resistivity_on || nernst_on) {
            Real jr = 0, jp = 0, jz = 0;   // curl B, code units
            Real dbz_dr = (bcc(IB3, k, j, i + 1) - bcc(IB3, k, j, i - 1)) / (pc->x1v(i + 1) - pc->x1v(i - 1));
            Real drbp_dr = (pc->x1v(i + 1) * bcc(IB2, k, j, i + 1) - pc->x1v(i - 1) * bcc(IB2, k, j, i - 1))
                           / (pc->x1v(i + 1) - pc->x1v(i - 1)) / pc->x1v(i);
            Real dbp_dz = 0, dbr_dz = 0, dbz_dp = 0, dbr_dp = 0;
            if (f3) {
              Real hz = pc->x3v(k + 1) - pc->x3v(k - 1);
              dbp_dz = (bcc(IB2, k + 1, j, i) - bcc(IB2, k - 1, j, i)) / hz;
              dbr_dz = (bcc(IB1, k + 1, j, i) - bcc(IB1, k - 1, j, i)) / hz;
            }
            if (f2) {
              Real hp = pc->x1v(i) * (pc->x2v(j + 1) - pc->x2v(j - 1));
              dbz_dp = (bcc(IB3, k, j + 1, i) - bcc(IB3, k, j - 1, i)) / hp;
              dbr_dp = (bcc(IB1, k, j + 1, i) - bcc(IB1, k, j - 1, i)) / hp;
            }
            jr = dbz_dp - dbp_dz; jp = dbr_dz - dbz_dr; jz = drbp_dr - dbr_dp;
            double heat = 0;   // code power density
            // the resistivity Athena++ applied this step (capped; with the Nernst upwinding), whose
            // dissipation is already in the total energy through the Poynting flux
            Real eta_used = pmb->pfield->fdif.etaB(FieldDiffusion::ohmic, k, j, i);
            ohmic = eta_used * (jr * jr + jp * jp + jz * jz) * lp::U_POWER * dt_s;
            heat += eta_used * (jr * jr + jp * jp + jz * jz);
            if (nernst_on) {
              Real g[3]; GradTe(pmb, k, j, i, g);
              Real kn = c(kKN, k, j, i) * lp::QE / lp::U_LEN / lp::U_VEL;
              Real e1 = kn * (g[1] * bcc(IB3, k, j, i) - g[2] * bcc(IB2, k, j, i));
              Real e2 = kn * (g[2] * bcc(IB1, k, j, i) - g[0] * bcc(IB3, k, j, i));
              Real e3 = kn * (g[0] * bcc(IB2, k, j, i) - g[1] * bcc(IB1, k, j, i));
              heat += jr * e1 + jp * e2 + jz * e3;
            }
            d_e += heat * lp::U_POWER * dt_s;
          }
          // ion-electron equilibration, exact for fixed heat capacities
          double ce = 1.5 * ne * lp::QE, ci = 1.5 * c(kNi, k, j, i) * lp::QE;   // J/m3/eV
          double s_za = (c(kZbe, k, j, i) * c(kZbe, k, j, i) * c(kNbe, k, j, i) / lp::A_BE
                         + c(kZl, k, j, i) * c(kZl, k, j, i) * c(kNlight, k, j, i) * lp::DT_Z2_OVER_A
                         + 4 * c(kZl, k, j, i) * c(kZl, k, j, i) * c(kNhe, k, j, i) / lp::A_HE);
          double tau1 = lp::TauE(ne, c(kNi, k, j, i), te, 1.0);
          double nu = 2 * lp::ME / lp::AMU * s_za / (ne * tau1);   // electron energy relaxation rate
          double t_eq = (ce * te + ci * ti) / (ce + ci);
          double dtemp = (te - ti) * std::exp(-nu * (1 + ce / ci) * dt_s);
          double te_new = t_eq + ci / (ce + ci) * dtemp;
          d_e += ce * (te_new - te);
          // keep the electrons between a cold floor and all of the thermal energy
          double eps_new = eps_e + d_e;
          double eps_max = c(kEpsTh, k, j, i) + d_tot + ohmic - 1.5 * c(kNi, k, j, i) * lp::QE * 1e-3;
          double eps_min = 1.5 * ne * lp::QE * 1e-3;
          eps_new = std::min(std::max(eps_new, eps_min), std::max(eps_max, eps_min));
          cons_scalar(4, k, j, i) = eps_new / lp::U_P;
        }
        cons(IEN, k, j, i) += d_tot / lp::U_P;
      }
}

// ---- time step: the Nernst advection (conduction is super-time-stepped) -------------
Real LinerTimeStep(MeshBlock *pmb) {
  if (!nernst_on) return std::numeric_limits<Real>::max();
  EnsurePlasma(pmb);
  AthenaArray<Real> &c = pmb->ruser_meshblock_data[0];
  Coordinates *pc = pmb->pcoord;
  Real dt = std::numeric_limits<Real>::max();
  for (int k = pmb->ks; k <= pmb->ke; ++k)
    for (int j = pmb->js; j <= pmb->je; ++j)
      for (int i = pmb->is; i <= pmb->ie; ++i) {
        if (c(kFvac, k, j, i) > 0.5) continue;
        Real g[3]; GradTe(pmb, k, j, i, g);
        double gm = std::sqrt(g[0] * g[0] + g[1] * g[1] + g[2] * g[2]) * lp::QE / lp::U_LEN;
        Real dxmin = std::min(pc->dx1f(i), pmb->pmy_mesh->f3 ? pc->dx3f(k) : pc->dx1f(i));
        if (nernst_on) {   // the field moves at v_N = K_N |grad T|: a Courant limit
          Real vn = c(kKN, k, j, i) * gm / lp::U_VEL;
          // (0.05: the upwinding diffusivity v_N dx / 2 then stays below the Ohmic step limit)
          if (vn > 0) dt = std::min(dt, static_cast<Real>(0.05 * dxmin / vn));
        }
      }
  return dt;
}
#endif  // LINER_PHYSICS

void DriveBoundary(MeshBlock *pmb, Coordinates *pco, AthenaArray<Real> &prim,
                   FaceField &b, Real time, Real dt,
                   int il, int iu, int jl, int ju, int kl, int ku, int ngh) {
  for (int k=kl; k<=ku; ++k) {
    for (int j=jl; j<=ju; ++j) {
      for (int i=1; i<=ngh; ++i) {
        prim(IDN,k,j,iu+i) = rho_vac;
        prim(IVX,k,j,iu+i) = 0.0;
        prim(IVY,k,j,iu+i) = 0.0;
        prim(IVZ,k,j,iu+i) = 0.0;
        prim(IPR,k,j,iu+i) = p0;
#if LINER_PHYSICS
        for (int n = 0; n < NSCALARS; ++n) pmb->pscalars->r(n, k, j, iu + i) = 0.0;
        pmb->pscalars->r(4, k, j, iu + i) = eps_e_vac / rho_vac;
#endif
      }
    }
  }
  Real current_ma = -1;
#if LINER_PHYSICS
  if (circuit_on) current_ma = pmb->pmy_mesh->ruser_mesh_data[0](kCurrent) * 1e-6;
#endif
  for (int k=kl; k<=ku; ++k)
    for (int j=jl; j<=ju; ++j)
      for (int i=2; i<=ngh+1; ++i) b.x1f(k,j,iu+i) = 0.0;
  for (int k=kl; k<=ku; ++k)
    for (int j=jl; j<=ju+1; ++j)
      for (int i=1; i<=ngh; ++i) {
        Real r = pco->x1v(iu+i);
        b.x2f(k,j,iu+i) = (current_ma >= 0) ? current_ma*1.0e6*5.642e-9/r : BPhi(r, time);
      }
  for (int k=kl; k<=ku+1; ++k)
    for (int j=jl; j<=ju; ++j)
      for (int i=1; i<=ngh; ++i) b.x3f(k,j,iu+i) = bz0;
}

Real Sum(MeshBlock *pmb, int n, bool r2) {
  Real s = 0.0;
  for (int k=pmb->ks; k<=pmb->ke; ++k)
    for (int j=pmb->js; j<=pmb->je; ++j)
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real w = pmb->pscalars->s(n,k,j,i)*pmb->pcoord->GetCellVolume(k,j,i);
        Real r = pmb->pcoord->x1v(i);
        s += r2 ? w*r*r : w;
      }
  return s;
}
Real GasMass(MeshBlock *pmb, int iout) { return Sum(pmb, 2, false); }
Real GasR2(MeshBlock *pmb, int iout) { return Sum(pmb, 2, true); }
Real IceMass(MeshBlock *pmb, int iout) { return Sum(pmb, 1, false); }
Real IceR2(MeshBlock *pmb, int iout) { return Sum(pmb, 1, true); }
Real GasThermal(MeshBlock *pmb, int iout) {  // sum of (gas fraction) * p * dV: the gas's thermal energy x 2/3
  Real s = 0.0;
  for (int k=pmb->ks; k<=pmb->ke; ++k)
    for (int j=pmb->js; j<=pmb->je; ++j)
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real f = pmb->pscalars->s(2,k,j,i)/pmb->phydro->u(IDN,k,j,i);
        s += f*pmb->phydro->w(IPR,k,j,i)*pmb->pcoord->GetCellVolume(k,j,i);
      }
  return s;
}
Real MaxGasDensity(MeshBlock *pmb, int iout) {
  Real m = 0.0;
  for (int k=pmb->ks; k<=pmb->ke; ++k)
    for (int j=pmb->js; j<=pmb->je; ++j)
      for (int i=pmb->is; i<=pmb->ie; ++i) m = std::max(m, pmb->pscalars->s(2,k,j,i));
  return m;
}

#if LINER_PHYSICS
// ---- history: the fusion power and what it depends on, for the whole liner (SI) -------
// 0 fusion power (W), 1 its share in cells where beryllium at least doubles the
// bremsstrahlung, 2-5 fusion power times T_i (keV), T_e (keV), |B_z| (T) and r (m), for
// burn-weighted means, 6 ash mass (kg); global: 7 current (MA), 8 load voltage (kV),
// 9 yield (J), 10 alpha energy made (J), 11 deposited (J), 12 radiated (J), 13 energy
// delivered into the domain (J), 14 energy left in the bank (J). Everything is for the
// whole liner (length_mm), not just the simulated slice.
Real BurnSum(MeshBlock *pmb, int what) {
  EnsurePlasma(pmb);
  AthenaArray<Real> &c = pmb->ruser_meshblock_data[0];
  Real s = 0;
  for (int k=pmb->ks; k<=pmb->ke; ++k)
    for (int j=pmb->js; j<=pmb->je; ++j)
      for (int i=pmb->is; i<=pmb->ie; ++i) {
        Real v = pmb->pcoord->GetCellVolume(k,j,i) * lp::U_VOL * length_scale;   // m3 of the real liner
        Real pf = c(kRate, k, j, i) * lp::E_FUSION * v;
        switch (what) {
          case 0: s += pf; break;
          case 1: if (c(kNbe, k, j, i) * lp::A_BE >= 0.225 * c(kNfuel, k, j, i) * lp::A_DT) s += pf; break;
          case 2: s += pf * c(kTi, k, j, i) * 1e-3; break;
          case 3: s += pf * c(kTe, k, j, i) * 1e-3; break;
          case 4: s += pf * std::abs(pmb->pfield->bcc(IB3, k, j, i)) * lp::U_B; break;
          case 5: s += pf * pmb->pcoord->x1v(i) * lp::U_LEN; break;
          case 6: s += pmb->pscalars->s(3, k, j, i) * lp::U_RHO * v; break;
        }
      }
  return s;
}
Real HistBurn(MeshBlock *pmb, int iout) { return BurnSum(pmb, iout - 6); }
Real HistGlobal(MeshBlock *pmb, int iout) {
  if (pmb->gid != 0) return 0.0;
  AthenaArray<Real> &m = pmb->pmy_mesh->ruser_mesh_data[0];
  switch (iout - 13) {
    case 0: return m(kCurrent) * 1e-6;
    case 1: return load_voltage * 1e-3;
    case 2: return m(kYield);
    case 3: return m(kAlphaMade);
    case 4: return m(kAlphaDep);
    case 5: return m(kRadiated);
    case 6: return m(kDelivered);
    case 7: return 0.5 * cap_c * m(kVcap) * m(kVcap);
    case 8: return cond_stages;
  }
  return 0.0;
}
#endif
}  // namespace

void Mesh::InitUserMeshData(ParameterInput *pin) {
  r_out = pin->GetReal("problem", "r_out");
  Real aspect = pin->GetReal("problem", "aspect");
  r_in = r_out*(1.0 - 1.0/aspect);
  r_gas = r_in - pin->GetReal("problem", "ice");
  rho_liner = pin->GetReal("problem", "rho_liner");
  rho_ice = pin->GetReal("problem", "rho_ice");
  rho_gas = pin->GetReal("problem", "rho_gas");
  rho_vac = pin->GetReal("problem", "rho_vac");
  p0 = pin->GetReal("problem", "p0");
  bz0 = pin->GetReal("problem", "bz")/35450.0;
  i_peak = pin->GetReal("problem", "i_peak");
  t_rise = pin->GetReal("problem", "t_rise");
  z_len = mesh_size.x3max - mesh_size.x3min;
  int modes = pin->GetOrAddInteger("problem", "modes", 32);
  int mmax = pin->GetOrAddInteger("problem", "mmax", 0);
  int mseed = pin->GetOrAddInteger("problem", "mseed", -1);
  if (mseed < 0) mseed = mmax;
  Real ripple = pin->GetOrAddReal("problem", "ripple", 0.0);  // rms, mm
  std::mt19937 gen(pin->GetOrAddInteger("problem", "seed", 1));
  std::uniform_real_distribution<Real> phase(0.0, 2.0*M_PI);
  int count = modes + mmax*(2*modes + 1);
  Real amp = count > 0 ? ripple*std::sqrt(2.0/count) : 0.0;  // equal amplitudes, rms = ripple over all modes
  seed_modes.clear();
  for (int n = 1; n <= modes; ++n) seed_modes.push_back({0, n, amp, phase(gen)});  // drawn first, as before
  for (int m = 1; m <= mmax; ++m)
    for (int n = -modes; n <= modes; ++n) {
      Real ph = phase(gen);
      if (m <= mseed) seed_modes.push_back({m, n, amp, ph});
    }
  EnrollUserBoundaryFunction(BoundaryFace::outer_x1, DriveBoundary);
#if LINER_PHYSICS
  burn_on = pin->GetOrAddBoolean("problem", "burn", false);
  alpha_mode = pin->GetOrAddInteger("problem", "alpha", -1);
  radiation_on = pin->GetOrAddBoolean("problem", "radiation", false);
  resistivity_on = pin->GetOrAddBoolean("problem", "resistivity", false);
  conduction_on = pin->GetOrAddBoolean("problem", "conduction", false);
  two_temp_on = pin->GetOrAddBoolean("problem", "two_temp", false);
  nernst_on = pin->GetOrAddBoolean("problem", "nernst", false);
  circuit_on = pin->GetOrAddBoolean("problem", "circuit", false);
  eta_vac = pin->GetOrAddReal("problem", "eta_vac", 1.26e-3);   // Ohm m: 1000 m2/s
  eta_max = pin->GetOrAddReal("problem", "eta_max", 1.26e-3);   // Ohm m
  flux_limit = pin->GetOrAddReal("problem", "flux_limit", 0.1);
  cond_safety = pin->GetOrAddReal("problem", "cond_safety", 0.4);
  if ((resistivity_on || nernst_on) && pin->GetOrAddReal("problem", "eta_ohm", 0.0) == 0.0)
    throw std::runtime_error("liner: resistivity or nernst needs <problem> eta_ohm = 1 (switches on field diffusion)");
  length_m = pin->GetOrAddReal("problem", "length_mm", 10.0) * 1e-3;
  length_scale = length_m / (z_len * lp::U_LEN);
  r_domain_m = mesh_size.x1max * lp::U_LEN;
  // the bank, as zpinch.circuit sizes it: a quarter sine to i_peak in t_rise into the
  // starting load, behind driver_nh; the return wall at wall_mm
  Real wall_m = pin->GetOrAddReal("problem", "wall_mm", 10.0) * 1e-3;
  l_machine = pin->GetOrAddReal("problem", "driver_nh", 15.0) * 1e-9;
  Real l_load = lp::MU0 * length_m / (2 * lp::kPi) * (std::log(wall_m / (r_out * lp::U_LEN)) + 0.25);
  Real l_sc = l_machine + l_load;
  cap_c = SQR(2 * t_rise * 1e-9 / lp::kPi) / l_sc;
  Real v0 = i_peak * 1e6 * std::sqrt(l_sc / cap_c);
  l_outside = lp::MU0 * length_m / (2 * lp::kPi) * std::log(std::max(wall_m / r_domain_m, 1.0));
  // the vacuum's electron energy at p0 (for the boundary)
  {
    double e; lp::FromRhoP(rho_vac, 0.0, p0 * lp::U_P, &e, nullptr);
    double t = lp::FromRhoE(rho_vac, 0.0, e).t_ev;
    eps_e_vac = lp::ElectronEnergy(0.0, rho_vac, t) / lp::U_P;
  }
  AllocateRealUserMeshDataField(1);
  ruser_mesh_data[0].NewAthenaArray(kNState);
  ruser_mesh_data[0](kVcap) = v0;
  ruser_mesh_data[0](kBank0) = 0.5 * cap_c * v0 * v0;
  alpha_eps.assign(mesh_size.nx1, 0.0);
  physics_source = burn_on || radiation_on || resistivity_on || conduction_on || two_temp_on || nernst_on;
  if (physics_source) EnrollUserExplicitSourceFunction(LinerSource);
  if (nernst_on) EnrollUserTimeStepFunction(LinerTimeStep);
  // the liner's diffusivity (zero with resistivity and nernst off) whenever field
  // diffusion is on (eta_ohm != 0; its value is not used)
  if (pin->GetOrAddReal("problem", "eta_ohm", 0.0) != 0.0) EnrollFieldDiffusivity(LinerDiffusivity);
  user_diffusion_emf = NernstEMF;
  AllocateUserHistoryOutput(22);
#else
  AllocateUserHistoryOutput(6);
#endif
  EnrollUserHistoryOutput(0, GasMass, "gas_mass");
  EnrollUserHistoryOutput(1, GasR2, "gas_r2");
  EnrollUserHistoryOutput(2, IceMass, "ice_mass");
  EnrollUserHistoryOutput(3, IceR2, "ice_r2");
  EnrollUserHistoryOutput(4, GasThermal, "gas_pv");
  EnrollUserHistoryOutput(5, MaxGasDensity, "gas_rho_max", UserHistoryOperation::max);
#if LINER_PHYSICS
  const char *burn_names[7] = {"fus_W", "fus_mixed_W", "fus_Ti", "fus_Te", "fus_Bz", "fus_r", "ash_kg"};
  for (int n = 0; n < 7; ++n) EnrollUserHistoryOutput(6 + n, HistBurn, burn_names[n]);
  const char *global_names[8] = {"current_MA", "vload_kV", "yield_J", "alpha_made_J", "alpha_dep_J",
                                 "radiated_J", "delivered_J", "bank_J"};
  for (int n = 0; n < 8; ++n) EnrollUserHistoryOutput(13 + n, HistGlobal, global_names[n]);
  EnrollUserHistoryOutput(21, HistGlobal, "cond_stages");   // RKL2 stages of the last conduction step
#endif
}

#if LINER_PHYSICS
void MeshBlock::InitUserMeshBlockData(ParameterInput *pin) {
  AllocateRealUserMeshBlockDataField(2);
  ruser_meshblock_data[0].NewAthenaArray(kNPlasma, ncells3, ncells2, ncells1);
  ruser_meshblock_data[1].NewAthenaArray(1);
  ruser_meshblock_data[1](0) = 0.0;   // the plasma state is not computed yet
}

// After each cycle: the circuit, the alphas' radial transport, and the energy totals.
void Mesh::UserWorkInLoop() {
  AthenaArray<Real> &m = ruser_mesh_data[0];
  const Real dt_s = dt * lp::U_TIME;
  int nb = nblocal;
  int nx1g = mesh_size.nx1;
  // the plasma state at the end of the step (the primitives are up to date, ghosts included):
  // for the conduction, the history and the next step's first stage
  if (physics_source) {
#pragma omp parallel for num_threads(GetNumMeshThreads()) schedule(dynamic)
    for (int b = 0; b < nb; ++b) {
      MeshBlock *pmb = my_blocks(b);
      int il, iu, jl, ju, kl, ku;
      PlasmaRange(pmb, &il, &iu, &jl, &ju, &kl, &ku);
      PlasmaState(pmb, pmb->phydro->w, pmb->pscalars->r, pmb->pfield->bcc, il, iu, jl, ju, kl, ku);
    }
  }
  Real escaped = (radiation_on && physics_source) ? RadiationTransport(this) * length_scale : 0.0;
  if (conduction_on && dt_s > 0) ConductSTS(this, dt_s);
  // energy totals from the cached rates (the last stage's)
  Real fus = 0, made = 0, dep = 0;
  std::vector<Real> src(nx1g, 0.0), loss(nx1g, 0.0), dif(nx1g, 0.0), vsum(nx1g, 0.0), rf(nx1g + 1, 0.0);
  Real p_in = 0;   // W into the domain through its outer boundary (the Poynting flux and any gas)
  for (int b = 0; b < nb; ++b) {
    MeshBlock *pmb = my_blocks(b);
    AthenaArray<Real> &c = pmb->ruser_meshblock_data[0];
    Coordinates *pc = pmb->pcoord;
    for (int k = pmb->ks; k <= pmb->ke; ++k)
      for (int j = pmb->js; j <= pmb->je; ++j)
        for (int i = pmb->is; i <= pmb->ie; ++i) {
          Real v = pc->GetCellVolume(k, j, i) * lp::U_VOL * length_scale;
          Real rate = c(kRate, k, j, i);
          fus += rate * lp::E_FUSION * v;
          made += rate * lp::E_ALPHA * v;
          int ig = static_cast<int>(pmb->loc.lx1) * pmb->block_size.nx1 + (i - pmb->is);
          if (alpha_mode == 1) {
            dep += alpha_eps[ig] / c(kAlphaTau, k, j, i) * v;
            src[ig] += rate * lp::E_ALPHA * v;
            loss[ig] += v / c(kAlphaTau, k, j, i);
            dif[ig] += c(kAlphaD, k, j, i) * v;
            vsum[ig] += v;
            rf[ig] = pc->x1f(i) * lp::U_LEN; rf[ig + 1] = pc->x1f(i + 1) * lp::U_LEN;
          } else if (alpha_mode == 0) {
            dep += rate * lp::E_ALPHA * v;
          }
        }
    if (pmb->pbval->block_bcs[BoundaryFace::outer_x1] == BoundaryFlag::user) {
      // the last stage's fluxes are the ones the update used (VL2)
      AthenaArray<Real> &f1 = pmb->phydro->flux[X1DIR];
      for (int k = pmb->ks; k <= pmb->ke; ++k)
        for (int j = pmb->js; j <= pmb->je; ++j)
          p_in -= f1(IEN, k, j, pmb->ie + 1) * pc->GetFace1Area(k, j, pmb->ie + 1) * lp::U_FLUX * 1e-6
                  * length_scale;
    }
  }
  if (burn_on) {
    m(kYield) += fus * dt_s; m(kAlphaMade) += made * dt_s; m(kAlphaDep) += dep * dt_s;
  }
  if (radiation_on) m(kRadiated) += escaped * dt_s;
  // alphas: steady radial diffusion with loss, -1/r d/dr (r D de/dr) + e L = S
  if (burn_on && alpha_mode == 1) {
    int n = nx1g;
    std::vector<Real> a(n, 0.0), bdiag(n, 0.0), cc(n, 0.0), d(n, 0.0);
    for (int i = 0; i < n; ++i) {
      Real vol = 0.5 * (rf[i + 1] * rf[i + 1] - rf[i] * rf[i]);   // per radian and length
      Real L = vsum[i] > 0 ? loss[i] / vsum[i] : 0.0;
      Real S = vsum[i] > 0 ? src[i] / vsum[i] : 0.0;
      bdiag[i] += L * vol; d[i] = S * vol;   // += : the coupling to cell i-1 is already in
      if (i + 1 < n) {
        Real dm = 0.5 * ((vsum[i] > 0 ? dif[i] / vsum[i] : 0) + (vsum[i + 1] > 0 ? dif[i + 1] / vsum[i + 1] : 0));
        Real rc0 = 0.5 * (rf[i] + rf[i + 1]), rc1 = 0.5 * (rf[i + 1] + rf[i + 2]);
        Real g = rf[i + 1] * dm / std::max(rc1 - rc0, 1e-12);
        bdiag[i] += g; bdiag[i + 1] += g; cc[i] = -g; a[i + 1] = -g;
      }
    }
    // the outer edge: the alphas leave (energy density zero just outside)
    for (int i = 1; i < n; ++i) {   // Thomas algorithm
      Real w = a[i] / bdiag[i - 1];
      bdiag[i] -= w * cc[i - 1]; d[i] -= w * d[i - 1];
    }
    alpha_eps[n - 1] = d[n - 1] / bdiag[n - 1];
    for (int i = n - 2; i >= 0; --i) alpha_eps[i] = (d[i] - cc[i] * alpha_eps[i + 1]) / bdiag[i];
    for (int i = 0; i < n; ++i) alpha_eps[i] = std::max(alpha_eps[i], static_cast<Real>(0.0));
  }
  // the circuit: bank -> machine inductance -> the inductance outside the domain -> the domain.
  // The domain's voltage is the power that crossed its boundary over the current, so the
  // bank, the circuit and the plasma share one energy account.
  Real cur = m(kCurrent), vcap = m(kVcap);
  m(kDelivered) += p_in * dt_s;
  Real v_dom = (std::abs(cur) > 1e4) ? p_in / cur : 0.0;
  if (circuit_on) {
    // until there is current, the domain is a fixed inductance (the vacuum gap)
    Real l_start = (std::abs(cur) > 1e4) ? 0.0
                   : lp::MU0 * length_m / (2 * lp::kPi) * std::log(r_domain_m / (r_out * lp::U_LEN));
    Real didt = (vcap - v_dom) / (l_machine + l_outside + l_start);
    load_voltage = v_dom + l_outside * didt;
    Real cur_new = cur + dt_s * didt;
    m(kCurrentPrev) = cur;
    m(kCurrent) = cur_new;
    m(kVcap) = vcap - cur * dt_s / cap_c;
  } else {
    m(kCurrentPrev) = cur;
    Real t1 = (time + dt), s = (t1 < t_rise) ? std::sin(0.5 * M_PI * t1 / t_rise) : 1.0;
    m(kCurrent) = i_peak * 1e6 * s * s;
    load_voltage = v_dom + l_outside * ((dt_s > 0) ? (m(kCurrent) - cur) / dt_s : 0.0);
  }
}
#endif

void MeshBlock::ProblemGenerator(ParameterInput *pin) {
#if LINER_PHYSICS
  Real gm1 = 0.0;
#else
  Real gm1 = peos->GetGamma() - 1.0;
#endif
  for (int k=ks; k<=ke; ++k) {
    Real z = pcoord->x3v(k);
    for (int j=js; j<=je; ++j) {
      Real surf = Surface(z, pcoord->x2v(j));
      for (int i=is; i<=ie; ++i) {
        // the share of each cell's radial width in each material (the ripple is far below a cell)
        Real lo = pcoord->x1f(i), hi = pcoord->x1f(i+1), w = hi - lo;
        auto share = [&](Real a, Real b2) {
          return std::max(0.0, std::min(hi, b2) - std::max(lo, a))/w; };
        Real f_gas = share(0.0, r_gas), f_ice = share(r_gas, r_in), f_lin = share(r_in, surf);
        Real f_vac = std::max(0.0, 1.0 - f_gas - f_ice - f_lin);
        Real m_gas = f_gas*rho_gas, m_ice = f_ice*rho_ice, m_lin = f_lin*rho_liner;
        Real rho = m_gas + m_ice + m_lin + f_vac*rho_vac;
        phydro->u(IDN,k,j,i) = rho;
        phydro->u(IM1,k,j,i) = 0.0;
        phydro->u(IM2,k,j,i) = 0.0;
        phydro->u(IM3,k,j,i) = 0.0;
        pscalars->s(0,k,j,i) = m_lin;
        pscalars->s(1,k,j,i) = m_ice;
        pscalars->s(2,k,j,i) = m_gas;
#if LINER_PHYSICS
        // everything at p0, through the equation of state; electrons at the same temperature
        pscalars->s(3,k,j,i) = 0.0;
        for (int n = 0; n < NSCALARS; ++n) pscalars->r(n,k,j,i) = pscalars->s(n,k,j,i)/rho;
        peos->SetCell(k, j, i);
        phydro->u(IEN,k,j,i) = peos->EgasFromRhoP(rho, p0) + 0.5*bz0*bz0;
        double e; lp::FromRhoP(rho, m_lin/rho, p0*lp::U_P, &e, nullptr);
        double t = lp::FromRhoE(rho, m_lin/rho, e).t_ev;
        pscalars->s(4,k,j,i) = lp::ElectronEnergy(m_lin, rho - m_lin, t)/lp::U_P;
        (void)gm1;
#else
        phydro->u(IEN,k,j,i) = p0/gm1 + 0.5*bz0*bz0;
#endif
      }
    }
  }
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie+1; ++i) pfield->b.x1f(k,j,i) = 0.0;
  for (int k=ks; k<=ke; ++k)
    for (int j=js; j<=je+1; ++j)
      for (int i=is; i<=ie; ++i) pfield->b.x2f(k,j,i) = 0.0;
  for (int k=ks; k<=ke+1; ++k)
    for (int j=js; j<=je; ++j)
      for (int i=is; i<=ie; ++i) pfield->b.x3f(k,j,i) = bz0;
}
