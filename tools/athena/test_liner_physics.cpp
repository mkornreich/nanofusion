// Checks of liner_physics.hpp: g++ -O2 -std=c++11 test_liner_physics.cpp -o /tmp/tlp && /tmp/tlp
#include <chrono>
#include <cstdio>
#include "liner_physics.hpp"

using namespace lp;

// Hugoniot from (rho0, e0, p0): the p at density rho with e - e0 = (p + p0)(1/rho0 - 1/rho)/2.
static double Hugoniot(double rho0, double fbe, double rho, double *t_ev) {
  double e0; FromRhoP(rho0, fbe, 1e5, &e0, nullptr);
  double v0 = 1 / (rho0 * U_RHO), v = 1 / (rho * U_RHO);
  double lo = FromRhoE(rho, fbe, 0.0).p_cold + 1e5, hi = 1e16;
  for (int it = 0; it < 200; ++it) {
    double p = std::sqrt(lo * hi), e; FromRhoP(rho, fbe, p, &e, nullptr);
    double f = e - e0 - 0.5 * (p + 1e5) * (v0 - v);
    if (f > 0) hi = p; else lo = p;
  }
  double p = std::sqrt(lo * hi), e; FromRhoP(rho, fbe, p, &e, nullptr);
  *t_ev = FromRhoE(rho, fbe, e).t_ev;
  return p;
}

int main() {
  auto t0 = std::chrono::steady_clock::now();
  Tables(0); Tables(1);
  auto t1 = std::chrono::steady_clock::now();
  std::printf("tables built in %.2f s\n", std::chrono::duration<double>(t1 - t0).count());

  std::printf("\nThomas-Fermi zbar\n");
  for (double t : {0.0, 1.0, 10.0, 100.0, 1000.0})
    std::printf("  Be 1.85 g/cc %7.0f eV: %.2f   D-T 0.25 g/cc: %.2f   D-T 3e-4 g/cc: %.2f\n", t,
                TFZbar(4, A_BE, 1.85, t), TFZbar(1, A_DT, 0.25, t), TFZbar(1, A_DT, 3e-4, t));

  std::printf("\nStarting states at 4e8 Pa\n");
  struct { const char *name; double rho, fbe; } st[] = {{"beryllium", 1.85, 1.0}, {"ice", 0.25, 0.0},
                                                        {"gas", 3e-4, 0.0}, {"vacuum", 1e-4, 0.0}};
  for (auto &s : st) {
    double e, c2; FromRhoP(s.rho, s.fbe, 4e8, &e, &c2);
    State x = FromRhoE(s.rho, s.fbe, e);
    std::printf("  %-9s e = %.3e J/kg  T = %.3g eV  p back = %.4e Pa  c = %.3g km/s\n", s.name, e, x.t_ev,
                x.p, std::sqrt(x.csq) / 1e3);
  }

  std::printf("\nBeryllium Hugoniot (LASL fit Us = 8.00 + 1.124 up km/s)\n");
  for (double up : {2.0, 4.0, 6.0, 10.0, 20.0}) {
    double us = 7.998 + 1.124 * up, eta = us / (us - up), p_data = 1850 * us * up * 1e6;
    double t, p = Hugoniot(1.85, 1.0, 1.85 * eta, &t);
    std::printf("  up %4.0f km/s: rho/rho0 %.2f  p data %.3g Mbar  model %.3g Mbar  (T %.3g eV)\n", up, eta,
                p_data / 1e11, p / 1e11, t);
  }
  std::printf("\nD-T ice Hugoniot from 0.25 g/cc\n");
  for (double eta : {2.0, 3.0, 3.5, 3.9}) {
    double t, p = Hugoniot(0.25, 0.0, 0.25 * eta, &t);
    std::printf("  rho/rho0 %.1f: %.3g Mbar (T %.3g eV)\n", eta, p / 1e11, t);
  }

  std::printf("\nRound trips p -> e -> p (relative error)\n");
  double worst = 0;
  for (double fbe : {0.0, 0.001, 0.003, 0.1, 0.5, 0.9, 0.997, 1.0})
    for (double rho : {1e-6, 1e-3, 0.3, 2.0, 30.0})
      for (double p : {1e3, 1e8, 1e11, 1e14}) {
        if (p < 1.01 * FromRhoE(rho, fbe, 0.0).p_cold) continue;  // below the cold curve: not a state
        double e, c2; FromRhoP(rho, fbe, p, &e, &c2);
        double p2 = FromRhoE(rho, fbe, e).p;
        double err = std::abs(p2 / p - 1);
        if (err > worst) worst = err;
        if (err > 1e-2) std::printf("  fbe %.3f rho %.0e p %.0e: back %.4e (%.2g)\n", fbe, rho, p, p2, err);
      }
  std::printf("  worst %.2g\n", worst);

  std::printf("\nMixed cells: across the pure threshold (rho 2 g/cc, e 1e6 J/kg)\n");
  for (double fbe : {0.0, 0.0019, 0.0021, 0.01, 0.5, 0.99, 0.9979, 0.9981, 1.0}) {
    State s = FromRhoE(2.0, fbe, 1e6);
    std::printf("  fbe %.4f: p %.4e Pa  T %.4g eV  c %.3g km/s\n", fbe, s.p, s.t_ev, std::sqrt(s.csq) / 1e3);
  }

  std::printf("\nElectron temperature round trip\n");
  for (double t : {0.05, 3.0, 100.0, 5000.0}) {
    double eps = ElectronEnergy(0.5, 1.0, t), zb, zl;
    double back = ElectronTemperature(0.5, 1.0, eps, 0.0, &zb, &zl);
    std::printf("  T %8.2f eV -> eps %.3e J/m3 -> %.4f eV (zbar Be %.2f, light %.2f)\n", t, eps, back, zb, zl);
  }

  std::printf("\nTransport and rates\n");
  double ne = 3e29, te = 5000;
  double tau = TauE(ne, ne, te, 1.0);
  std::printf("  tau_e(3e29, 5 keV) = %.3g s  (expect ~5.7e-14)\n", tau);
  std::printf("  <sigma v> D-T at 10 keV = %.3g m3/s (expect 1.1e-22)\n", SigmaVDT(10.0));
  Alpha a = AlphaSlowing(ne, te, DT_Z2_OVER_A);
  std::printf("  alpha at 3e29, 5 keV: range %.3g um (rho R %.3g g/cm2), tau_E %.3g ps, to ions %.2f\n",
              a.range * 1e6, a.range * 100 * ne * A_DT * AMU / 1e3 * 1e-3 * 1e3, a.tau_e * 1e12, a.f_ions);
  double br = ThinRadiation(ne, te, 0, 0, ne, 1, 0);
  std::printf("  D-T bremsstrahlung at 3e29, 5 keV: %.3g W/m3 (5.35e-37 ne^2 sqrt(T) = %.3g)\n", br,
              5.35e-37 * ne * ne * std::sqrt(5.0));
  std::printf("  Be at 1 g/cc, 100 eV: thin %.3g W/m3, limited (100 um) %.3g\n",
              ThinRadiation(2.7e29, 100, 6.7e28, TFZbar(4, A_BE, 1.0, 100), 0, 1, 0),
              LimitedRadiation(ThinRadiation(2.7e29, 100, 6.7e28, TFZbar(4, A_BE, 1.0, 100), 0, 1, 0), 100, 1e-4));

  // speed
  auto s0 = std::chrono::steady_clock::now();
  double acc = 0; int n = 0;
  for (int it = 0; it < 200000; ++it) {
    double rho = 1e-3 * (1 + it % 1000), e = 1e7 * (1 + it % 777);
    acc += FromRhoE(rho, 0.0, e).p; ++n;
  }
  auto s1 = std::chrono::steady_clock::now();
  for (int it = 0; it < 20000; ++it) {
    double rho = 1e-3 * (1 + it % 1000), e = 1e7 * (1 + it % 777);
    acc += FromRhoE(rho, 0.3, e).p;
  }
  auto s2 = std::chrono::steady_clock::now();
  std::printf("\nspeed: pure %.0f ns/call, mixed %.0f ns/call (%g)\n",
              std::chrono::duration<double>(s1 - s0).count() / n * 1e9,
              std::chrono::duration<double>(s2 - s1).count() / 20000 * 1e9, acc > 0 ? 1.0 : 0.0);
  return 0;
}
