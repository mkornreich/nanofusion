//========================================================================================
// Plasma physics for the nanofusion liner runs (tools/athena/liner.cpp and liner_eos.cpp).
//
// Header-only and SI inside; code units (mm, ns, g/cm3) are converted at the edges.
//
// * Equation of state: two materials, beryllium and "light" (D-T; the vacuum and helium
//   ash are treated as D-T here). Each is
//     cold curve (no tension) + ideal ions + electrons that go from degenerate to ideal,
//     with Thomas-Fermi ionisation (More's fit) and the ionisation energy that goes with it.
//   Built into tables at start-up. Mixed cells: the cold parts share one pressure (the
//   materials sit side by side), the thermal parts add (Dalton) at one temperature.
// * Collisions, Braginskii transport (Z = 1 coefficients), the Lee-More minimum mean free
//   path for cold dense matter, D-T reactivity (Bosch & Hale), bremsstrahlung,
//   recombination and beryllium K-shell line radiation, alpha slowing down: the same
//   formulas as nanofusion/improved (burn1d.py, confinement.py) where those exist.
//========================================================================================
#ifndef NANOFUSION_LINER_PHYSICS_HPP_
#define NANOFUSION_LINER_PHYSICS_HPP_

#include <algorithm>
#include <cmath>
#include <vector>

namespace lp {

// --- constants (SI) ------------------------------------------------------------------
constexpr double AMU = 1.66053907e-27;
constexpr double QE = 1.602176634e-19;     // also J per eV
constexpr double ME = 9.1093837e-31;
constexpr double HBAR = 1.054571817e-34;
constexpr double MU0 = 1.25663706e-6;
constexpr double SIGMA_SB = 5.670374e-8;
constexpr double kPi = 3.14159265358979323846;
constexpr double KEV = 1.602176634e-16;

// --- code units: mm, ns, g/cm3 ----------------------------------------------------------
constexpr double U_RHO = 1e3;        // kg/m3
constexpr double U_LEN = 1e-3;       // m
constexpr double U_TIME = 1e-9;      // s
constexpr double U_VEL = 1e6;        // m/s
constexpr double U_P = 1e15;         // Pa = J/m3
constexpr double U_E = 1e12;         // J/kg
constexpr double U_POWER = 1e24;     // W/m3
constexpr double U_FLUX = 1e21;      // W/m2
constexpr double U_DIFF = 1e3;       // m2/s
constexpr double U_B = 35449.077;    // T per code unit: sqrt(mu0 * 1e15 Pa)
constexpr double U_VOL = 1e-9;       // m3 per mm3

// --- species ------------------------------------------------------------------------------
constexpr double A_BE = 9.0122, A_DT = 2.5, A_HE = 4.0026;
constexpr double E_FUSION = 17.59e3 * KEV;   // J per D-T reaction
constexpr double E_ALPHA_KEV = 3520.0;
constexpr double E_ALPHA = E_ALPHA_KEV * KEV;
constexpr double DT_Z2_OVER_A = 0.5 / 2.014 + 0.5 / 3.016;  // per D-T ion

struct Material {
  double Z, A, rho0;      // nuclear charge, mass (u), reference density (g/cm3)
  double a_cold, n_cold;  // cold curve p = a (x^5/3 - x^n), x = rho / rho0 > 1 (Pa)
  int nion;
  double ion_ev[4];       // successive ionisation energies (eV)
};

// Beryllium: bulk modulus 130 GPa at rho0, Fermi-gas stiffness (all 4 electrons free) far above.
// Light (D-T): solid D-T density; the cold curve is the electrons' degeneracy pressure above it.
inline double FermiPressure(double ne) {  // 2/5 n_e E_F (Pa)
  double ef = HBAR * HBAR * std::pow(3 * kPi * kPi * ne, 2.0 / 3.0) / (2 * ME);
  return 0.4 * ne * ef;
}
struct Materials { Material m[2]; };
inline Materials MakeMaterials() {
  Materials s = {{{4.0, A_BE, 1.85, 0.0, 0.0, 4, {9.3227, 18.2111, 153.8961, 217.7186}},
                  {1.0, A_DT, 0.25, 0.0, 0.0, 1, {13.598, 0, 0, 0}}}};
  double pf_be = FermiPressure(4.0 * 1.85e3 / (A_BE * AMU));
  s.m[0].a_cold = pf_be;
  s.m[0].n_cold = 5.0 / 3.0 - 130e9 / pf_be;
  s.m[1].a_cold = FermiPressure(0.25e3 / (A_DT * AMU));
  s.m[1].n_cold = 0.0;
  return s;
}
inline const Material &Mat(int m) {
  static const Materials mats = MakeMaterials();  // initialised once, thread-safe (C++11)
  return mats.m[m];
}
enum {BE = 0, LIGHT = 1};

// Thomas-Fermi mean ionisation, More (1985) fit as given by Salzmann; rho g/cm3, T eV.
inline double TFZbar(double Z, double A, double rho, double t_ev) {
  const double a1 = 0.003323, a2 = 0.9718, a3 = 9.26148e-5, a4 = 3.10165;
  const double b0 = -1.7630, b1 = 1.43175, b2 = 0.31546, c1 = -0.366667, c2 = 0.983333;
  double r = std::max(rho, 1e-30) / (Z * A);
  double t0 = std::max(t_ev, 0.0) / std::pow(Z, 4.0 / 3.0);
  double tf = t0 / (1 + t0);
  double aa = a1 * std::pow(t0, a2) + a3 * std::pow(t0, a4);
  double b = -std::exp(b0 + b1 * tf + b2 * std::pow(tf, 7));
  double c = c1 * tf + c2;
  double q1 = aa * std::pow(r, b);
  double q = std::pow(std::pow(r, c) + std::pow(q1, c), 1.0 / c);
  double x = 14.3139 * std::pow(q, 0.6624);
  return Z * x / (1 + x + std::sqrt(1 + 2 * x));
}

// Energy (eV) to strip an atom down to mean charge zbar (linear between stages).
inline double IonisationEnergy(const Material &m, double zbar) {
  double e = 0.0;
  for (int s = 0; s < m.nion; ++s) {
    if (zbar >= s + 1) e += m.ion_ev[s];
    else { e += (zbar - s) * m.ion_ev[s]; break; }
  }
  return std::max(e, 0.0);
}

// Cold curves (zero-temperature pressure above the reference density; no tension below it):
// * beryllium: p = a (x^5/3 - x^n), a = its Fermi pressure with all 4 electrons free, n set by
//   the bulk modulus (130 GPa). Matches the measured Hugoniot to a few % up to 11 Mbar.
// * light (D-T): Vinet (K0 = 0.3 GPa, K0' = 7, solid hydrogen-like) at low compression,
//   blending into the electrons' Fermi pressure, p_F(rho) - p_F(rho0), at high compression.
inline double ColdPressureFormula(int m, double x) {
  const Material &mat = Mat(m);
  if (x <= 1.0) return 0.0;
  double x53 = std::pow(x, 5.0 / 3.0);
  if (m == BE) return mat.a_cold * (x53 - std::pow(x, mat.n_cold));
  double eta = std::cbrt(1.0 / x), k0 = 0.3e9, k1 = 7.0;
  double vinet = 3 * k0 * (1 - eta) / (eta * eta) * std::exp(1.5 * (k1 - 1) * (1 - eta));
  double fermi = mat.a_cold * (x53 - 1);
  double s = (x - 1) * (x - 1) / ((x - 1) * (x - 1) + 25.0);
  return vinet + (fermi - vinet) * s;
}
struct ColdTable { double dl; int n; std::vector<double> p, e, k; };  // over log10 x
inline ColdTable MakeColdTable(int m) {
  ColdTable t; t.dl = 2e-4; t.n = 20001;   // x from 1 to 1e4
  t.p.resize(t.n); t.e.resize(t.n); t.k.resize(t.n);
  double rho0 = Mat(m).rho0 * U_RHO;
  for (int i = 0; i < t.n; ++i) t.p[i] = ColdPressureFormula(m, std::pow(10.0, i * t.dl));
  t.e[0] = 0;
  for (int i = 1; i < t.n; ++i) {  // e = integral p / rho^2 drho = integral p / (rho0 x^2) dx
    double xa = std::pow(10.0, (i - 1) * t.dl), xb = std::pow(10.0, i * t.dl);
    t.e[i] = t.e[i - 1] + 0.5 * (t.p[i - 1] / (xa * xa) + t.p[i] / (xb * xb)) * (xb - xa) / rho0;
  }
  for (int i = 0; i < t.n; ++i) {  // K = rho dp/drho = dp / dln x
    int a = std::max(i - 1, 0), b = std::min(i + 1, t.n - 1);
    t.k[i] = (t.p[b] - t.p[a]) / ((b - a) * t.dl * std::log(10.0));
  }
  return t;
}
inline const ColdTable &ColdTab(int m) {
  static const ColdTable tabs[2] = {MakeColdTable(BE), MakeColdTable(LIGHT)};
  return tabs[m];
}
// Cold curve: pressure (Pa), energy (J/kg) and dp/drho (m2/s2) at density rho (g/cm3).
inline void Cold(const Material &m, double rho, double *p, double *e, double *dpdr) {
  int id = (m.Z > 2) ? BE : LIGHT;
  const ColdTable &t = ColdTab(id);
  double rho0 = m.rho0 * U_RHO;
  double x = rho / m.rho0;
  if (x <= 1.0) {
    // no tension; the sound speed keeps the solid's stiffness just below rho0 (ramped to zero
    // at 0.9 rho0) so that it does not jump where the material starts to be compressed
    *p = 0; *e = 0;
    *dpdr = t.k[0] / rho0 * std::max(0.0, (x - 0.9) / 0.1);
    return;
  }
  double s = std::log10(x) / t.dl;
  int i = std::min(static_cast<int>(s), t.n - 2);
  double w = s - i;
  *p = (1 - w) * t.p[i] + w * t.p[i + 1];
  *e = (1 - w) * t.e[i] + w * t.e[i + 1];
  *dpdr = ((1 - w) * t.k[i] + w * t.k[i + 1]) / (rho * U_RHO);
}
// Inverse cold curve: compression x (>= 1) at cold pressure p (Pa).
inline double ColdCompression(int m, double p) {
  const ColdTable &t = ColdTab(m);
  if (p <= 0) return 1.0;
  if (p >= t.p[t.n - 1]) return std::pow(10.0, (t.n - 1) * t.dl);
  int lo = 0, hi = t.n - 1;
  while (hi - lo > 1) { int mid = (lo + hi) / 2; if (t.p[mid] <= p) lo = mid; else hi = mid; }
  double w = (p - t.p[lo]) / std::max(t.p[hi] - t.p[lo], 1e-300);
  return std::pow(10.0, (lo + w) * t.dl);
}

// Thermal part at (rho g/cm3, T eV): pressure (Pa), energy (J/kg), the electrons' share of
// that energy (J/kg, with the ionisation energy) and the mean ionisation.
inline void Thermal(const Material &m, double rho, double t_ev, double *p, double *e,
                    double *e_el, double *zbar) {
  double ni = rho * U_RHO / (m.A * AMU);
  double z = TFZbar(m.Z, m.A, rho, t_ev);
  double z0 = TFZbar(m.Z, m.A, rho, 0.0);
  double ne = z * ni;
  double kt = t_ev * QE;
  double c = 0.4 * HBAR * HBAR * std::pow(3 * kPi * kPi * ne, 2.0 / 3.0) / (2 * ME);
  double pe = ne * kt * kt / (std::sqrt(kt * kt + c * c) + c);  // n_e (sqrt(T^2 + c^2) - c)
  double pi = ni * kt;
  double eion = ni * std::max(0.0, IonisationEnergy(m, z) - IonisationEnergy(m, z0)) * QE;
  double rho_si = rho * U_RHO;
  *p = pi + pe;
  *e_el = (1.5 * pe + eion) / rho_si;
  *e = *e_el + 1.5 * pi / rho_si;
  *zbar = z;
}

// --- tables ---------------------------------------------------------------------------
struct Grid {
  double x0, dx; int n;
  double at(int i) const { return x0 + i * dx; }
  // index and weight for bilinear interpolation, clamped to the grid
  void locate(double x, int *i, double *w) const {
    double s = (x - x0) / dx;
    if (s <= 0) { *i = 0; *w = 0; return; }
    if (s >= n - 1) { *i = n - 2; *w = 1; return; }
    *i = static_cast<int>(s); *w = s - *i;
  }
};

struct Table2 {
  Grid gx, gy; std::vector<double> v;
  void init(Grid x, Grid y) { gx = x; gy = y; v.assign(static_cast<size_t>(x.n) * y.n, 0.0); }
  double &at(int i, int j) { return v[static_cast<size_t>(i) * gy.n + j]; }
  double at(int i, int j) const { return v[static_cast<size_t>(i) * gy.n + j]; }
  double interp(int i, double wi, int j, double wj) const {
    return (1 - wi) * ((1 - wj) * at(i, j) + wj * at(i, j + 1))
        + wi * ((1 - wj) * at(i + 1, j) + wj * at(i + 1, j + 1));
  }
};

// Per material: forward tables over (log10 rho, log10 T) and inverse tables over
// (log10 rho, log10 e_th) and (log10 rho, log10 p_th) for single-material cells.
struct MaterialTables {
  Table2 le, lp, lee, zb;      // log10 e_th (J/kg), log10 p_th (Pa), log10 e_el (J/kg), zbar
  Table2 inv_e_lp, inv_e_lt, inv_e_gam;   // over (lrho, le): log10 p_th, log10 T, Gamma_th
  Table2 inv_p_le, inv_p_gam;             // over (lrho, lp): log10 e_th, Gamma_th
};

inline const Grid &RhoGrid() { static Grid g{-10.0, 0.05, 281}; return g; }   // 1e-10..1e4 g/cm3
inline const Grid &TGrid() { static Grid g{-3.0, 0.025, 361}; return g; }     // 1 meV..1 MeV
inline const Grid &EGrid() { static Grid g{1.0, 0.025, 641}; return g; }      // 10..1e17 J/kg
inline const Grid &PGrid() { static Grid g{-6.0, 0.025, 1281}; return g; }    // 1e-6..1e26 Pa

// Invert a monotone column y(t) (both log10) for target y; returns the fractional
// index into the column. Outside the column, extrapolate with slope 1 (ideal gas).
inline double InvertColumn(const double *col, int n, double y, double *extra) {
  *extra = 0;
  if (y <= col[0]) { *extra = y - col[0]; return 0.0; }
  if (y >= col[n - 1]) { *extra = y - col[n - 1]; return n - 1.0; }
  int lo = 0, hi = n - 1;
  while (hi - lo > 1) { int mid = (lo + hi) / 2; if (col[mid] <= y) lo = mid; else hi = mid; }
  double d = col[hi] - col[lo];
  return lo + (d > 0 ? (y - col[lo]) / d : 0.0);
}

inline void BuildTables(int m, MaterialTables *t) {
  const Material &mat = Mat(m);
  const Grid &gr = RhoGrid(), &gt = TGrid(), &ge = EGrid(), &gp = PGrid();
  t->le.init(gr, gt); t->lp.init(gr, gt); t->lee.init(gr, gt); t->zb.init(gr, gt);
  for (int i = 0; i < gr.n; ++i) {
    double rho = std::pow(10.0, gr.at(i));
    for (int j = 0; j < gt.n; ++j) {
      double p, e, ee, z;
      Thermal(mat, rho, std::pow(10.0, gt.at(j)), &p, &e, &ee, &z);
      t->le.at(i, j) = std::log10(std::max(e, 1e-300));
      t->lp.at(i, j) = std::log10(std::max(p, 1e-300));
      t->lee.at(i, j) = std::log10(std::max(ee, 1e-300));
      t->zb.at(i, j) = z;
    }
    // keep the columns strictly increasing (round-off at the cold end)
    for (int j = 1; j < gt.n; ++j) {
      t->le.at(i, j) = std::max(t->le.at(i, j), t->le.at(i, j - 1) + 1e-9);
      t->lp.at(i, j) = std::max(t->lp.at(i, j), t->lp.at(i, j - 1) + 1e-9);
    }
  }
  t->inv_e_lp.init(gr, ge); t->inv_e_lt.init(gr, ge); t->inv_e_gam.init(gr, ge);
  t->inv_p_le.init(gr, gp); t->inv_p_gam.init(gr, gp);
  std::vector<double> ce(gt.n), cp(gt.n);
  for (int i = 0; i < gr.n; ++i) {
    for (int j = 0; j < gt.n; ++j) { ce[j] = t->le.at(i, j); cp[j] = t->lp.at(i, j); }
    for (int j = 0; j < ge.n; ++j) {
      double extra, s = InvertColumn(ce.data(), gt.n, ge.at(j), &extra);
      int k = std::min(static_cast<int>(s), gt.n - 2); double w = s - k;
      t->inv_e_lt.at(i, j) = (1 - w) * gt.at(k) + w * gt.at(k + 1) + extra;
      t->inv_e_lp.at(i, j) = (1 - w) * cp[k] + w * cp[k + 1] + extra;
    }
    for (int j = 0; j < gp.n; ++j) {
      double extra, s = InvertColumn(cp.data(), gt.n, gp.at(j), &extra);
      int k = std::min(static_cast<int>(s), gt.n - 2); double w = s - k;
      t->inv_p_le.at(i, j) = (1 - w) * ce[k] + w * ce[k + 1] + extra;
    }
  }
  // Gamma_th = dln p/dln rho|e + (p / rho e) dln p/dln e|rho, from the inverse table
  for (int i = 0; i < gr.n; ++i) {
    int i0 = std::max(i - 1, 0), i1 = std::min(i + 1, gr.n - 1);
    for (int j = 0; j < ge.n; ++j) {
      int j0 = std::max(j - 1, 0), j1 = std::min(j + 1, ge.n - 1);
      double dr = (t->inv_e_lp.at(i1, j) - t->inv_e_lp.at(i0, j)) / ((i1 - i0) * gr.dx);
      double de = (t->inv_e_lp.at(i, j1) - t->inv_e_lp.at(i, j0)) / ((j1 - j0) * ge.dx);
      double ratio = std::pow(10.0, t->inv_e_lp.at(i, j) - ge.at(j) - gr.at(i) - 3.0);  // p/(rho e)
      t->inv_e_gam.at(i, j) = std::max(1.0, std::min(3.0, dr + ratio * de));
    }
    for (int j = 0; j < gp.n; ++j) {
      // Gamma at the e that gives this p
      double le = t->inv_p_le.at(i, j);
      int jj; double w; ge.locate(le, &jj, &w);
      t->inv_p_gam.at(i, j) = (1 - w) * t->inv_e_gam.at(i, jj) + w * t->inv_e_gam.at(i, jj + 1);
    }
  }
}

struct AllTables { MaterialTables t[2]; };
inline AllTables *MakeTables() {
  AllTables *a = new AllTables;
  BuildTables(BE, &a->t[0]);
  BuildTables(LIGHT, &a->t[1]);
  return a;
}
inline const MaterialTables &Tables(int m) {
  static const AllTables *tabs = MakeTables();  // built once, thread-safe (C++11)
  return tabs->t[m];
}

// --- single-material EOS (fast path); rho g/cm3, specific energies J/kg, p Pa ------------
// Thermal pressure from thermal energy, with Gamma_th; below the table's density the
// intensive values of its lowest density are used.
inline void PureFromE(int m, double rho, double eth, double *pth, double *gam, double *t_ev) {
  const MaterialTables &t = Tables(m);
  const Grid &gr = RhoGrid(), &ge = EGrid();
  double lr = std::log10(std::max(rho, 1e-300));
  double lrc = std::max(lr, gr.x0);
  double le = std::log10(std::max(eth, 1e-300));
  int i, j; double wi, wj;
  gr.locate(lrc, &i, &wi);
  double lec = std::min(std::max(le, ge.x0), ge.at(ge.n - 1));
  ge.locate(lec, &j, &wj);
  double extra = le - lec;  // outside the energy grid: ideal-gas slope
  double lpv = t.inv_e_lp.interp(i, wi, j, wj) + extra + (lr - lrc);
  *pth = std::pow(10.0, lpv);
  if (gam) *gam = t.inv_e_gam.interp(i, wi, j, wj);
  if (t_ev) *t_ev = std::pow(10.0, t.inv_e_lt.interp(i, wi, j, wj) + extra);
}
inline void PureFromP(int m, double rho, double pth, double *eth, double *gam) {
  const MaterialTables &t = Tables(m);
  const Grid &gr = RhoGrid(), &gp = PGrid();
  double lr = std::log10(std::max(rho, 1e-300));
  double lrc = std::max(lr, gr.x0);
  double lpv = std::log10(std::max(pth, 1e-300)) - (lr - lrc);
  int i, j; double wi, wj;
  gr.locate(lrc, &i, &wi);
  double lpc = std::min(std::max(lpv, gp.x0), gp.at(gp.n - 1));
  gp.locate(lpc, &j, &wj);
  *eth = std::pow(10.0, t.inv_p_le.interp(i, wi, j, wj) + (lpv - lpc));
  if (gam) *gam = t.inv_p_gam.interp(i, wi, j, wj);
}

// Forward table values at (rho, log10 T): log10 e_th, log10 p_th, log10 e_el, zbar and
// d(log10 e)/d(log10 T), d(log10 p)/d(log10 T), d(log10 e_el)/d(log10 T).
struct Fwd { double le, lp, lee, z, dle, dlp, dlee; };
inline Fwd Forward(int m, double rho, double lt) {
  const MaterialTables &t = Tables(m);
  const Grid &gr = RhoGrid(), &gt = TGrid();
  double lr = std::log10(std::max(rho, 1e-300));
  double lrc = std::max(lr, gr.x0);
  double ltc = std::min(std::max(lt, gt.x0), gt.at(gt.n - 1));
  int i, j; double wi, wj;
  gr.locate(lrc, &i, &wi); gt.locate(ltc, &j, &wj);
  Fwd f;
  double extra = lt - ltc;
  f.le = t.le.interp(i, wi, j, wj) + extra;
  f.lp = t.lp.interp(i, wi, j, wj) + extra + (lr - lrc);
  f.lee = t.lee.interp(i, wi, j, wj) + extra;
  f.z = t.zb.interp(i, wi, j, wj);
  auto slope = [&](const Table2 &tb) {
    if (extra != 0) return 1.0;
    double a = (1 - wi) * tb.at(i, j) + wi * tb.at(i + 1, j);
    double b = (1 - wi) * tb.at(i, j + 1) + wi * tb.at(i + 1, j + 1);
    return (b - a) / gt.dx;
  };
  f.dle = slope(t.le); f.dlp = slope(t.lp); f.dlee = slope(t.lee);
  return f;
}

// --- mixtures of beryllium and light material --------------------------------------
// Cold part: both materials at one pressure, side by side (volume fractions add to one).
struct ColdMix { double p, e, k; };   // Pa, J/kg, bulk modulus (Pa)
inline ColdMix ColdOfMix(double rho, double fbe) {
  const Material &mb = Mat(BE), &ml = Mat(LIGHT);
  double fl = 1 - fbe;
  ColdMix c{0, 0, 0};
  if (fbe / mb.rho0 + fl / ml.rho0 <= 1.0 / rho) return c;  // room to spare: no cold pressure
  if (fbe < 1e-12 || fl < 1e-12) {
    int m = fbe > 0.5 ? BE : LIGHT;
    double dp; Cold(Mat(m), rho, &c.p, &c.e, &dp);
    c.k = rho * U_RHO * dp; return c;
  }
  // solve fbe / rho_be(p) + fl / rho_l(p) = 1 / rho (Newton on p, bracketed)
  auto inv = [](const Material &m, double p, double *rho_m, double *k_m, double *e_m) {
    *rho_m = m.rho0 * ColdCompression(m.Z > 2 ? BE : LIGHT, p);
    double pp, dp; Cold(m, *rho_m, &pp, e_m, &dp);
    *k_m = std::max(*rho_m * U_RHO * dp, 1e3);
  };
  double lo = 0.0, hi = 1e20, p = 1e9;
  for (int it = 0; it < 80; ++it) {
    double rb, kb, eb, rl, kl, el;
    inv(mb, p, &rb, &kb, &eb); inv(ml, p, &rl, &kl, &el);
    double g = fbe / rb + fl / rl - 1.0 / rho;     // > 0: needs more pressure
    if (g > 0) lo = p; else hi = p;
    double dg = -(fbe / (rb * kb) + fl / (rl * kl));  // d(1/rho_m)/dp = -1 / (rho_m K_m)
    double pn = p - g / dg;
    if (!(pn > lo && pn < hi)) pn = (hi < 1e19) ? 0.5 * (lo + hi) : 2 * p + 1e5;
    if (std::abs(pn - p) < 1e-10 * p) { p = pn; break; }
    p = pn;
  }
  double rb, kb, eb, rl, kl, el;
  inv(mb, p, &rb, &kb, &eb); inv(ml, p, &rl, &kl, &el);
  c.p = p; c.e = fbe * eb + fl * el;
  double ab = fbe * rho / rb, al = fl * rho / rl;   // volume fractions
  c.k = 1.0 / (ab / kb + al / kl);
  return c;
}

// Mixture thermal part at one temperature (Dalton): solve for log10 T from the thermal
// energy (by_energy) or the thermal pressure.
inline double MixTemperature(double rho, double fbe, double target, bool by_energy,
                             double *eth, double *pth) {
  double fl = 1 - fbe, rb = fbe * rho, rl = fl * rho;
  double lt = 1.0;
  double ltarget = std::log10(std::max(target, 1e-300));
  double rho_si = rho * U_RHO;
  for (int it = 0; it < 40; ++it) {
    Fwd b = Forward(BE, std::max(rb, 1e-30), lt), l = Forward(LIGHT, std::max(rl, 1e-30), lt);
    double vb, vl, db, dl;
    if (by_energy) {
      vb = fbe * std::pow(10.0, b.le); vl = fl * std::pow(10.0, l.le); db = b.dle; dl = l.dle;
    } else {   // partial pressures at the partial densities
      vb = (rb > 1e-30) ? std::pow(10.0, b.lp) : 0; vl = (rl > 1e-30) ? std::pow(10.0, l.lp) : 0;
      db = b.dlp; dl = l.dlp;
    }
    double v = vb + vl;
    double lv = std::log10(std::max(v, 1e-300));
    double slope = (vb * db + vl * dl) / std::max(v, 1e-300);
    double step = (ltarget - lv) / std::max(slope, 0.2);
    step = std::max(-2.0, std::min(2.0, step));
    lt += step;
    if (std::abs(step) < 1e-9) break;
  }
  Fwd b = Forward(BE, std::max(rb, 1e-30), lt), l = Forward(LIGHT, std::max(rl, 1e-30), lt);
  *eth = fbe * std::pow(10.0, b.le) + fl * std::pow(10.0, l.le);
  *pth = ((rb > 1e-30) ? std::pow(10.0, b.lp) : 0) + ((rl > 1e-30) ? std::pow(10.0, l.lp) : 0);
  (void)rho_si;
  return lt;
}

constexpr double PURE = 0.002;  // beryllium share below which a cell is pure light (and 1 - this)

// Full EOS: pressure (Pa) and sound speed squared (m2/s2) from (rho g/cm3, beryllium mass
// share, specific internal energy J/kg); also the temperature (eV) and the cold energy.
struct State { double p, csq, t_ev, e_cold, p_cold; };
inline State FromRhoE(double rho, double fbe, double e) {
  State s;
  fbe = std::min(std::max(fbe, 0.0), 1.0);
  double rho_si = rho * U_RHO;
  if (fbe < PURE || fbe > 1 - PURE) {
    int m = fbe > 0.5 ? BE : LIGHT;
    double pc, ec, dpc; Cold(Mat(m), rho, &pc, &ec, &dpc);
    double eth = std::max(e - ec, 1e-3);
    double pth, gam, t;
    PureFromE(m, rho, eth, &pth, &gam, &t);
    s.p = pc + pth; s.csq = dpc + gam * pth / rho_si; s.t_ev = t; s.e_cold = ec; s.p_cold = pc;
    return s;
  }
  ColdMix c = ColdOfMix(rho, fbe);
  double eth = std::max(e - c.e, 1e-3), e_out, pth;
  double lt = MixTemperature(rho, fbe, eth, true, &e_out, &pth);
  s.p = c.p + pth; s.csq = c.k / rho_si + (5.0 / 3.0) * pth / rho_si;
  s.t_ev = std::pow(10.0, lt); s.e_cold = c.e; s.p_cold = c.p;
  return s;
}
// Specific internal energy (J/kg) and sound speed squared from (rho, fbe, p Pa).
inline void FromRhoP(double rho, double fbe, double p, double *e, double *csq) {
  fbe = std::min(std::max(fbe, 0.0), 1.0);
  double rho_si = rho * U_RHO;
  if (fbe < PURE || fbe > 1 - PURE) {
    int m = fbe > 0.5 ? BE : LIGHT;
    double pc, ec, dpc; Cold(Mat(m), rho, &pc, &ec, &dpc);
    double pth = std::max(p - pc, 1e-6 * std::max(p, 1e-6)), eth, gam;
    PureFromP(m, rho, pth, &eth, &gam);
    *e = ec + eth; if (csq) *csq = dpc + gam * pth / rho_si;
    return;
  }
  ColdMix c = ColdOfMix(rho, fbe);
  double pth = std::max(p - c.p, 1e-6 * std::max(p, 1e-6)), eth, p_out;
  MixTemperature(rho, fbe, pth, false, &eth, &p_out);
  *e = c.e + eth;
  if (csq) *csq = c.k / rho_si + (5.0 / 3.0) * pth / rho_si;
}

// Electron temperature (eV) from the electrons' energy (J/m3) in a cell of light and
// beryllium (partial densities in g/cm3): solve rho_b e_el,b(T) + rho_l e_el,l(T) = eps.
inline double ElectronTemperature(double rho_be, double rho_l, double eps_e, double lt_guess,
                                  double *z_be, double *z_l) {
  // a material with under 1e-4 of the mass is left out of the solve (its zbar still comes out)
  double tot = rho_be + rho_l;
  bool use_b = rho_be > 1e-4 * tot, use_l = rho_l > 1e-4 * tot;
  double lt = lt_guess;
  double ltarget = std::log10(std::max(eps_e, 1e-300));
  Fwd b{}, l{};
  for (int it = 0; it < 40; ++it) {
    double vb = 0, vl = 0, sb = 0, sl = 0;
    if (use_b) { b = Forward(BE, rho_be, lt); vb = rho_be * U_RHO * std::pow(10.0, b.lee); sb = vb * b.dlee; }
    if (use_l) { l = Forward(LIGHT, rho_l, lt); vl = rho_l * U_RHO * std::pow(10.0, l.lee); sl = vl * l.dlee; }
    double v = vb + vl;
    double slope = (sb + sl) / std::max(v, 1e-300);
    double step = (ltarget - std::log10(std::max(v, 1e-300))) / std::max(slope, 0.2);
    step = std::max(-2.0, std::min(2.0, step));
    lt += step;
    if (std::abs(step) < 1e-5) break;   // 2e-5 relative in T
  }
  // zbar at the final temperature (the last step is tiny: the table's value is close enough)
  *z_be = use_b ? b.z : Forward(BE, std::max(rho_be, 1e-30), lt).z;
  *z_l = use_l ? l.z : Forward(LIGHT, std::max(rho_l, 1e-30), lt).z;
  return std::pow(10.0, lt);
}
// Cold energy (J/kg) of a cell (single material or mixed).
inline double ColdEnergy(double rho, double fbe) {
  fbe = std::min(std::max(fbe, 0.0), 1.0);
  if (fbe < PURE || fbe > 1 - PURE) {
    double p, e, d; Cold(Mat(fbe > 0.5 ? BE : LIGHT), rho, &p, &e, &d); return e;
  }
  return ColdOfMix(rho, fbe).e;
}
// Electron energy (J/m3) at temperature T (eV).
inline double ElectronEnergy(double rho_be, double rho_l, double t_ev) {
  double lt = std::log10(std::max(t_ev, 1e-6));
  Fwd b = Forward(BE, std::max(rho_be, 1e-30), lt), l = Forward(LIGHT, std::max(rho_l, 1e-30), lt);
  return rho_be * U_RHO * std::pow(10.0, b.lee) + rho_l * U_RHO * std::pow(10.0, l.lee);
}

// --- collisions and transport ---------------------------------------------------------
inline double CoulombLog(double ne, double te_kev) {  // as burn1d.coulomb_log
  return std::max(2.0, 24.0 - std::log(std::sqrt(ne * 1e-6) / (te_kev * 1e3)));
}
inline double FermiEnergy(double ne) {
  return HBAR * HBAR * std::pow(3 * kPi * kPi * ne, 2.0 / 3.0) / (2 * ME);
}
// Electron collision time (s): Spitzer-Braginskii for ions of effective charge zeff, but no
// shorter than the time to cross one ion spacing (Lee-More minimum mean free path).
inline double TauE(double ne, double ni, double te_ev, double zeff) {
  double t = std::max(te_ev, 1e-3);
  double tau = 3.44e5 * std::pow(t, 1.5) / (ne * 1e-6 * std::max(zeff, 1e-3) * CoulombLog(ne, t * 1e-3));
  double a = std::cbrt(3.0 / (4 * kPi * std::max(ni, 1.0)));
  double v = std::sqrt(3 * t * QE / ME + 2 * FermiEnergy(ne) / ME);
  return std::max(tau, a / v);
}
// Braginskii coefficients for Z = 1, x = omega_ce tau_e.
inline double BragDelta(double x) { return x * x * x * x + 14.79 * x * x + 3.7703; }
inline double AlphaPerp(double x) { return 1.0 - (6.416 * x * x + 1.837) / BragDelta(x); }
inline double KappaPerp(double x) { return (4.664 * x * x + 11.92) / BragDelta(x); }
constexpr double KAPPA_PAR = 3.1616;
// beta_wedge / x: the Nernst coefficient without its 1/x singularity
inline double BetaWedgeOverX(double x) { return (1.5 * x * x + 3.053) / BragDelta(x); }

// D-T reactivity (m3/s), Bosch & Hale (1992), as burn1d.bosch_hale_dt.
inline double SigmaVDT(double t_kev) {
  if (t_kev < 0.2) return 0.0;
  double t = std::min(t_kev, 100.0);
  const double c1 = 1.17302e-9, c2 = 1.51361e-2, c3 = 7.51886e-2, c4 = 4.60643e-3, c5 = 1.35e-2,
      c6 = -1.0675e-4, c7 = 1.366e-5;
  double theta = t / (1 - t * (c2 + t * (c4 + t * c6)) / (1 + t * (c3 + t * (c5 + t * c7))));
  double xi = std::cbrt(34.3827 * 34.3827 / (4 * theta));
  return c1 * theta * std::sqrt(xi / (1124656.0 * t * t * t)) * std::exp(-3 * xi) * 1e-6;
}

// Alpha slowing down (confinement.alpha_range_factor and burn1d): path length in units of
// v0 tau_se, u_c = sqrt(E_c / E_0).
inline double AlphaRangeFactor(double uc) {
  if (uc < 1e-9) return 1.0;
  double a = uc;
  double lg = std::log((1 + a) * (1 + a) / (1 - a + a * a)) / 6;
  double at = (std::atan((2 - a) / (a * std::sqrt(3.0))) + kPi / 6) / std::sqrt(3.0);
  return 1.0 - a * (lg + at);
}
const double ALPHA_V0 = std::sqrt(2 * E_ALPHA / (A_HE * AMU));
const double ALPHA_LARMOR_TM = A_HE * AMU * ALPHA_V0 / (2 * QE);   // ~0.27 T m
struct Alpha { double range, tau_e, f_ions; };   // m, s (energy loss time), share to ions
inline Alpha AlphaSlowing(double ne, double te_ev, double sza) {
  double teff = std::sqrt(te_ev * te_ev + std::pow(0.6 * FermiEnergy(ne) / QE, 2));
  double tau_se = 6.27e8 * A_HE * std::pow(teff, 1.5) / (4 * ne * 1e-6 * CoulombLog(ne, teff * 1e-3));
  double ec = 14.8 * A_HE * teff * 1e-3 * std::pow(std::max(sza, 1e-6), 2.0 / 3.0);  // keV
  Alpha a;
  a.range = ALPHA_V0 * tau_se * AlphaRangeFactor(std::sqrt(std::min(ec / E_ALPHA_KEV, 0.999)));
  a.tau_e = tau_se / (2 * (1 + std::pow(2 * ec / E_ALPHA_KEV, 1.5)));
  a.f_ions = 1.0 / (1.0 + 32.0 / std::max(te_ev * 1e-3, 1e-6));
  return a;
}

// --- radiation (W/m3) ------------------------------------------------------------------
constexpr double BREMS = 5.35e-37;   // W m3 keV^-1/2, per n_e sum(n Z^2)
// Optically thin bremsstrahlung + recombination of every ion, and K-shell line emission
// of beryllium (van Regemorter excitation of the H- and He-like ions, coronal).
// n_be, n_l: ion densities; z_be, z_l: their mean charges; n_he: helium ash (bare).
inline double ThinRadiation(double ne, double te_ev, double n_be, double z_be, double n_l,
                            double z_l, double n_he) {
  double t = std::max(te_ev, 1e-2);
  double sq = std::sqrt(t * 1e-3);
  double brems = BREMS * ne * sq * (n_be * z_be * z_be + n_l * z_l * z_l + 4 * n_he);
  // recombination: n_e sqrt(T) sum n Z^2 (E_inf / T), E_inf the binding of the captured electron
  double rec = BREMS * ne * sq * (n_l * z_l * z_l * 13.6 / t + n_he * 4 * 54.4 / t);
  // beryllium: bare and H-like fractions from zbar
  double fz = z_be - std::floor(z_be);
  int zi = static_cast<int>(std::floor(z_be));
  const double ip[5] = {0, 9.32, 18.21, 153.9, 217.7};  // binding of the electron captured by charge q
  double rec_be = 0;
  if (zi >= 1) rec_be += (1 - fz) * zi * zi * ip[std::min(zi, 4)];
  if (zi + 1 <= 4) rec_be += fz * (zi + 1) * (zi + 1) * ip[zi + 1];
  rec += BREMS * ne * sq * n_be * rec_be / t;
  // lines: ions with K-shell electrons (zbar between 2 and 4)
  double line = 0;
  if (z_be > 2 && z_be < 4) {
    double h_like = (z_be >= 3) ? (4 - z_be) : (z_be - 2);   // Be3+ share
    double he_like = (z_be < 3) ? (3 - z_be) : 0.0;           // Be2+ share
    auto vr = [&](double de, double f) {   // van Regemorter, cm3/s -> m3/s; g = 0.2
      return 1.6e-5 * f * 0.2 * std::exp(-de / t) / (de * std::sqrt(t)) * 1e-6;
    };
    line = ne * n_be * (h_like * vr(163.3, 0.416) * 163.3 + he_like * vr(123.7, 0.55) * 123.7) * QE;
  }
  return brems + rec + line;
}
// Optically thick limit: a cell cannot lose more than a black body of depth `depth` (m).
inline double LimitedRadiation(double thin, double te_ev, double depth) {
  double tk = std::max(te_ev, 1e-3) * 11604.518;
  double bb = 4 * SIGMA_SB * tk * tk * tk * tk / depth;
  return thin / (1 + thin / bb);
}

// --- radiation transport along r: multigroup, column by column -----------------------
// Each cell emits its optically thin power eta (bremsstrahlung, recombination, lines) with a
// bremsstrahlung spectrum, exp(-E/T), split into photon groups. In each group its opacity is
// Kirchhoff's, alpha_g = eta_g / (4 sigma T^4 P_g) with P_g the group's share of the Planck
// spectrum, so opaque matter radiates as a black body and the opacity falls as E^-3 above T.
// Photons go along r. A cell of optical depth d lets out (1/2 - E3(d)) / (2d) of its own
// emission through each face and transmits 2 E3(d) of a flux crossing it; the rest is
// absorbed there. Two sweeps per group carry the outward and inward beams; the inward beam
// crosses the axis and goes out through the column again; what passes the outer edge leaves.
// In watts, so energy is conserved exactly (the spreading of rays in a cylinder is left out).
constexpr int NGROUP = 17;
inline const double *GroupEdgesKeV() {   // NGROUP + 1 edges, the last infinite
  static const double e[NGROUP + 1] = {0.0, 0.001, 0.003, 0.01, 0.03, 0.1, 0.2, 0.4, 0.7, 1.0,
                                       2.0, 4.0, 7.0, 10.0, 20.0, 50.0, 100.0, 1e30};
  return e;
}
// share of the Planck spectrum above x = E / T
inline double PlanckAbove(double x) {
  if (x <= 0) return 1.0;
  if (x > 700) return 0.0;
  const double norm = 15.0 / (kPi * kPi * kPi * kPi);
  if (x < 2.0) {   // 1 - integral_0^x t^3/(e^t - 1) dt
    double x2 = x * x, x3 = x2 * x;
    double below = x3 / 3 - x2 * x2 / 8 + x3 * x2 / 60 - x3 * x2 * x2 / 5040 + x3 * x3 * x3 / 272160
                   - x3 * x3 * x3 * x2 / 13305600;
    return std::max(0.0, 1.0 - norm * below);
  }
  double sum = 0;
  for (int n = 1; n <= 12; ++n) {
    double en = std::exp(-n * x);
    if (en < 1e-300) break;
    sum += en * (x * x * x / n + 3 * x * x / (n * n) + 6 * x / (n * n * n) + 6.0 / (n * n * n * n));
  }
  return norm * sum;
}
inline double ExpIntE1(double x) {   // Abramowitz & Stegun 5.1.53 and 5.1.56
  if (x <= 0) return 1e300;
  if (x <= 1) return -std::log(x) - 0.57721566 + x * (0.99999193 + x * (-0.24991055 + x * (0.05519968
                     + x * (-0.00976004 + x * 0.00107857))));
  return std::exp(-x) / x * (x * x + 2.334733 * x + 0.250621) / (x * x + 3.330657 * x + 1.681534);
}
inline double ExpIntE3(double x) {
  if (x <= 0) return 0.5;
  if (x > 700) return 0.0;
  return 0.5 * (std::exp(-x) * (1 - x) + x * x * ExpIntE1(x));
}
// One radial column, cell 0 at the axis: volumes (m3), widths (m), electron temperatures (eV)
// and thin emission eta (W/m3) -> net radiative loss per unit volume (W/m3; negative where
// the cell absorbs more than it emits). Returns the power that leaves the column (W).
inline double ColumnTransport(int n, const double *vol, const double *dr, const double *t_ev,
                              const double *eta, double *loss, std::vector<double> &work) {
  const double *edge = GroupEdgesKeV();
  work.assign(static_cast<size_t>(n) * 3, 0.0);
  double *emit = work.data(), *trans = emit + n, *absorbed = trans + n;   // per group, reused
  std::vector<double> leave(n, 0.0), gain(n, 0.0);
  double escaped = 0;
  for (int g = 0; g < NGROUP; ++g) {
    bool any = false;
    for (int i = 0; i < n; ++i) {
      emit[i] = 0; trans[i] = 1; absorbed[i] = 0;
      if (eta[i] <= 0) continue;
      double t = std::max(t_ev[i], 1e-3) * 1e-3;   // keV
      double lo = edge[g] / t, hi = edge[g + 1] / t;
      double e_g = eta[i] * (std::exp(-lo) - (hi > 700 ? 0.0 : std::exp(-hi)));   // W/m3 in group
      if (e_g <= 0) continue;
      double tk = std::max(t_ev[i], 1e-3) * 11604.518;
      double bb = 4 * SIGMA_SB * tk * tk * tk * tk * (PlanckAbove(lo) - PlanckAbove(hi));
      double d = (bb > 0) ? e_g / bb * dr[i] : 1e30;                                // optical depth
      double out_share = (d < 1e-6) ? 0.5 : (0.5 - ExpIntE3(d)) / (2 * d);         // per face
      emit[i] = e_g * vol[i] * out_share;   // W through each face
      trans[i] = 2 * ExpIntE3(d);
      leave[i] += 2 * emit[i];
      any = true;
    }
    if (!any) continue;
    // outward beam
    double beam = 0;
    for (int i = 0; i < n; ++i) {
      double through = beam * trans[i];
      absorbed[i] += beam - through;
      beam = through + emit[i];
    }
    escaped += beam;
    // inward beam, then across the axis and out again
    beam = 0;
    for (int i = n - 1; i >= 0; --i) {
      double through = beam * trans[i];
      absorbed[i] += beam - through;
      beam = through + emit[i];
    }
    for (int i = 0; i < n; ++i) {
      double through = beam * trans[i];
      absorbed[i] += beam - through;
      beam = through;
    }
    escaped += beam;
    for (int i = 0; i < n; ++i) gain[i] += absorbed[i];
  }
  for (int i = 0; i < n; ++i) loss[i] = (leave[i] - gain[i]) / vol[i];
  return escaped;
}

}  // namespace lp
#endif  // NANOFUSION_LINER_PHYSICS_HPP_
