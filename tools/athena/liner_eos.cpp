//========================================================================================
// Athena++ general EOS for the nanofusion liner (copy to src/eos/general/liner_eos.cpp,
// configure with --eos=general/liner_eos; needs the core patch athena-physics.patch, which
// lets the EOS know which cell it is evaluating, and liner_physics.hpp next to it).
//
// <hydro> eos_model = liner (default): beryllium and D-T with cold curves, Thomas-Fermi
//                     ionisation and degenerate electrons (liner_physics.hpp).
//                   = ideal: gamma-law gas (gamma), to compare with the adiabatic build.
// The beryllium share of a cell's mass is passive scalar 0.
//========================================================================================
#include <algorithm>
#include <cmath>
#include <string>

#include "../../athena.hpp"
#include "../../hydro/hydro.hpp"
#include "../../mesh/mesh.hpp"
#include "../../parameter_input.hpp"
#include "../../scalars/scalars.hpp"
#include "../eos.hpp"
#include "liner_physics.hpp"

namespace {
bool ideal_eos = false;
Real gamma_ideal = 5.0 / 3.0;

// Read the beryllium share of cell (k,j,i): from the conserved scalars (during
// conserved-to-primitive) or from the primitive ones (everywhere else).
void LinerSetCell(EquationOfState *peos, int k, int j, int i, bool cons) {
  MeshBlock *pmb = peos->GetBlock();
  Real f;
  if (cons) {
    Real d = pmb->phydro->u(IDN, k, j, i);
    f = (d > 0) ? pmb->pscalars->s(0, k, j, i) / d : 0.0;
  } else {
    f = pmb->pscalars->r(0, k, j, i);
  }
  peos->eos_comp_[0] = std::min(static_cast<Real>(1.0), std::max(static_cast<Real>(0.0), f));
}

// The Riemann solver asks for the energy and then the sound speed of the same two states:
// remember the last two (per thread; each MeshBlock is worked on by one thread at a time).
struct CacheEntry { double rho, f, p, e, csq; };
thread_local CacheEntry cache[2] = {{-1, -1, -1, 0, 0}, {-1, -1, -1, 0, 0}};
thread_local int cache_next = 0;
void FromRhoPCached(double rho, double f, double p, double *e, double *csq) {
  for (const CacheEntry &c : cache)
    if (c.rho == rho && c.p == p && c.f == f) { *e = c.e; *csq = c.csq; return; }
  lp::FromRhoP(rho, f, p, e, csq);
  cache[cache_next] = {rho, f, p, *e, *csq};
  cache_next ^= 1;
}
}  // namespace

Real EquationOfState::PresFromRhoEg(Real rho, Real egas) {
  if (ideal_eos) return (gamma_ideal - 1.0) * egas;
  lp::State s = lp::FromRhoE(rho, eos_comp_[0], egas / rho * lp::U_E);
  return s.p / lp::U_P;
}

Real EquationOfState::EgasFromRhoP(Real rho, Real pres) {
  if (ideal_eos) return pres / (gamma_ideal - 1.0);
  double e, csq;
  FromRhoPCached(rho, eos_comp_[0], pres * lp::U_P, &e, &csq);
  return rho * e / lp::U_E;
}

Real EquationOfState::AsqFromRhoP(Real rho, Real pres) {
  if (ideal_eos) return gamma_ideal * pres / rho;
  double e, csq;
  FromRhoPCached(rho, eos_comp_[0], pres * lp::U_P, &e, &csq);
  return csq / (lp::U_VEL * lp::U_VEL);
}

void EquationOfState::InitEosConstants(ParameterInput *pin) {
  std::string model = pin->GetOrAddString("hydro", "eos_model", "liner");
  ideal_eos = (model == "ideal");
  gamma_ideal = pin->GetOrAddReal("hydro", "gamma", 5.0 / 3.0);
  set_cell_ = LinerSetCell;
  if (!ideal_eos) { lp::Tables(lp::BE); lp::Tables(lp::LIGHT); }  // build the tables now
}
