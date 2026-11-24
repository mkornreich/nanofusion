"""Simulation engine, a port of org.nanofusion.simulation.Simulator.

The Java engine does all arithmetic in ``BigDecimal``: most operations are
exact, and selected ones are rounded to ``precision`` significant digits
(HALF_UP). Those roundings shape the results, so this port reproduces them
operation by operation with :mod:`decimal` and matches the Java output digit
for digit (see tests/test_engine.py).

Differences from the Java code, all in cases where Java crashed:

* Each particle interaction is computed once per pair instead of twice. The
  force sums are exact, so the result is identical.
* The Lorentz force works without an electric field (Java threw a
  NullPointerException).
* An overflowing exponent in log-time mode gives a zero step, matching Java's
  ``Math.exp`` returning Infinity.
"""

from __future__ import annotations

import math
import time as _time
from decimal import Decimal
from typing import Callable, Optional

import numpy as np

from . import constants as K
from .constants import EXACT, ZERO
from .model import Frame, FusionEvent, SimulationResult, SimulatorProperties, lattice_indices

# Offsets to neighbours at lattice (Manhattan) distance 1 or 2, i.e.
# Java LVector.isNearNeighbor.
_NEIGHBOUR_OFFSETS = [
    (dx, dy, dz)
    for dx in range(-2, 3)
    for dy in range(-2, 3)
    for dz in range(-2, 3)
    if 0 < abs(dx) + abs(dy) + abs(dz) < 3
]


class SimulationCancelled(Exception):
    pass


def neighbour_pairs(x_size: int, y_size: int, z_size: int) -> list[tuple[int, int]]:
    """All particle-id pairs (i < j) that interact, in a lattice of the given size."""
    pairs = []
    for k in range(x_size * y_size * z_size):
        x, y, z = k % x_size, (k // x_size) % y_size, k // (x_size * y_size)
        for dx, dy, dz in _NEIGHBOUR_OFFSETS:
            nx, ny, nz = x + dx, y + dy, z + dz
            if 0 <= nx < x_size and 0 <= ny < y_size and 0 <= nz < z_size:
                j = nx + ny * x_size + nz * x_size * y_size
                if j > k:
                    pairs.append((k, j))
    return pairs


def detect_fusion(result: SimulationResult) -> Optional[list[FusionEvent]]:
    """Find threshold crossings from stored float positions (for loaded files).

    Like the engine, it checks the positions each step starts from, so the
    last frame is not checked. Returns ``None`` if the lattice layout is
    unknown.
    """
    p = result.properties
    n = len(result.types)
    if p.particle_count != n:
        return None
    pairs = np.array(neighbour_pairs(p.lattice_x_size, p.lattice_y_size, p.lattice_z_size), dtype=int)
    events: list[FusionEvent] = []
    if len(pairs) == 0:
        return events
    threshold = float(p.threshold)
    for index, frame in enumerate(result.frames[:-1]):
        d = np.linalg.norm(frame.positions[pairs[:, 0]] - frame.positions[pairs[:, 1]], axis=1)
        for k in np.nonzero(d <= threshold)[0]:
            events.append(FusionEvent(index, int(pairs[k, 0]), int(pairs[k, 1]), Decimal(float(d[k]))))
    return events


class Simulator:
    def __init__(self, props: SimulatorProperties):
        props.validate()
        self.props = props
        self._ctx = K.rounding_context(props.precision)
        self.time: Decimal = ZERO
        self._started = _time.monotonic()
        self.fusion_events: list[FusionEvent] = []
        self._frame_count = 0
        self._build_lattice()
        self._build_pairs()

    # --- set-up -------------------------------------------------------------

    def _build_lattice(self) -> None:
        """Java ``Lattice(LVector size, int precision, PVector electricField)``."""
        p = self.props
        add = EXACT.add
        px: list[Decimal] = []
        py: list[Decimal] = []
        pz: list[Decimal] = []
        types: list[str] = []
        deut = True
        z = ZERO
        for _layer in range(p.lattice_z_size):
            y = ZERO
            deut = not deut
            for _column in range(p.lattice_y_size):
                x = ZERO
                deut = not deut
                for _row in range(p.lattice_x_size):
                    types.append(K.DEUTERIUM if deut else K.LITHIUM)
                    deut = not deut
                    px.append(x)
                    py.append(y)
                    pz.append(z)
                    x = add(x, K.INITIAL_LATTICE_SPACING)
                y = add(y, K.INITIAL_LATTICE_SPACING)
            z = add(z, K.INITIAL_LATTICE_SPACING)
        n = len(types)
        self.px, self.py, self.pz = px, py, pz
        self.vx = [ZERO] * n
        self.vy = [ZERO] * n
        self.vz = [ZERO] * n
        self.types = types
        self.is_li = [t == K.LITHIUM for t in types]

        mul = EXACT.multiply
        self._charge = {False: K.DEUTERIUM_CHARGE, True: K.LITHIUM_CHARGE}
        self._inverse_mass = {False: K.INVERSE_DEUTERIUM_MASS, True: K.INVERSE_LITHIUM_MASS}
        # Particle constructor: fieldForce = E * charge * e, lorentzMultiplier = charge * e / c
        self._field_force = None
        if p.electric_field is not None:
            self._field_force = {
                li: tuple(mul(mul(e, q), K.CHARGE) for e in p.electric_field)
                for li, q in self._charge.items()
            }
        self._lorentz = {li: mul(mul(q, K.CHARGE), K.INVERSE_C) for li, q in self._charge.items()}

    def _build_pairs(self) -> None:
        p = self.props
        xs, ys, zs = p.lattice_x_size, p.lattice_y_size, p.lattice_z_size
        ctx = self._ctx
        # CommonCalculatorUtility.getInteraction:
        #   potential = POTENTIAL.multiply(c_actor, mc).multiply(c_particle, mc)
        # Rounding makes this order-dependent, so keep both orders.
        pot = {
            (a, b): ctx.multiply(ctx.multiply(K.POTENTIAL, self._charge[a]), self._charge[b])
            for a in (False, True)
            for b in (False, True)
        }
        li = self.is_li
        # (i, j, potential of i acting on j, potential of j acting on i)
        self._pairs = [
            (i, j, pot[(li[i], li[j])], pot[(li[j], li[i])]) for i, j in neighbour_pairs(xs, ys, zs)
        ]

    # --- stepping -----------------------------------------------------------

    def _time_step(self, step: int) -> Decimal:
        ts = self.props.time_step
        if self.props.use_log_time:
            expnt = 1.0 * (step - self.props.transition_step) / float(self.props.damping_factor)
            try:
                e = math.exp(expnt)
            except OverflowError:
                e = math.inf
            ts = EXACT.multiply(ts, Decimal(1.0 / (1.0 + e)))
        return ts

    def step(self, iteration: int) -> Frame:
        """Advance one time step (Java ``Simulator.step``) and return the new frame."""
        ex_add, ex_sub, ex_mul = EXACT.add, EXACT.subtract, EXACT.multiply
        r_mul, r_div, r_add = self._ctx.multiply, self._ctx.divide, self._ctx.add
        sqrt = math.sqrt

        dt = self._time_step(iteration)
        self.time = ex_add(self.time, dt)
        t = self.time
        px, py, pz = self.px, self.py, self.pz
        vx, vy, vz = self.vx, self.vy, self.vz
        n = len(px)
        fx = [ZERO] * n
        fy = [ZERO] * n
        fz = [ZERO] * n
        threshold = self.props.threshold
        frame_index = self._frame_count  # distances are measured on the frame this step starts from

        # Coulomb forces between near neighbours.
        for i, j, pot_ij, pot_ji in self._pairs:
            xi, yi, zi = px[i], py[i], pz[i]
            xj, yj, zj = px[j], py[j], pz[j]
            dxi = ex_sub(xi, xj)
            dyi = ex_sub(yi, yj)
            dzi = ex_sub(zi, zj)
            # getSquareDistance: each square rounded, sum exact
            sq = ex_add(ex_add(r_mul(dxi, dxi), r_mul(dyi, dyi)), r_mul(dzi, dzi))
            # CommonCalculatorUtility.sqrt(BigDecimal): one Newton step from the double sqrt
            sf = float(sq)
            if sf == 0.0:
                raise ArithmeticError(f"Particles {i} and {j} are at the same position.")
            xf = sqrt(sf)
            xd = Decimal(xf)
            distance = ex_add(xd, Decimal(float(ex_sub(sq, ex_mul(xd, xd))) / (xf * 2.0)))
            denom = r_mul(distance, sq)
            um_ij = r_div(pot_ij, denom)  # i acting on j
            um_ji = um_ij if pot_ji is pot_ij else r_div(pot_ji, denom)
            # Force on j points from i to j; force on i from j to i.
            fx[j] = ex_add(fx[j], ex_mul(ex_sub(xj, xi), um_ij))
            fy[j] = ex_add(fy[j], ex_mul(ex_sub(yj, yi), um_ij))
            fz[j] = ex_add(fz[j], ex_mul(ex_sub(zj, zi), um_ij))
            fx[i] = ex_add(fx[i], ex_mul(dxi, um_ji))
            fy[i] = ex_add(fy[i], ex_mul(dyi, um_ji))
            fz[i] = ex_add(fz[i], ex_mul(dzi, um_ji))
            if distance <= threshold:
                self.fusion_events.append(FusionEvent(frame_index, i, j, distance))

        # Oscillating electric field: F = q E sin(omega t), same for every particle of a type.
        field_now = None
        if self._field_force is not None:
            wave = Decimal(math.sin(float(r_mul(K.OMEGA, t))))
            field_now = {li: tuple(ex_mul(c, wave) for c in ff) for li, ff in self._field_force.items()}
        lorentz = self.props.include_lorentz

        for k in range(n):
            li = self.is_li[k]
            fxk, fyk, fzk = fx[k], fy[k], fz[k]
            if field_now is not None:
                ffx, ffy, ffz = field_now[li]
                fxk, fyk, fzk = ex_add(fxk, ffx), ex_add(fyk, ffy), ex_add(fzk, ffz)
            vxk, vyk, vzk = vx[k], vy[k], vz[k]
            if lorentz:
                lm = self._lorentz[li]
                fxk = ex_add(fxk, ex_mul(lm, vzk))
                fyk = ex_add(fyk, ZERO)
                fzk = ex_add(fzk, ex_sub(ZERO, ex_mul(lm, vxk)))
            # Particle.update: explicit Euler, position first (old velocity), each result rounded.
            px[k] = r_add(px[k], ex_mul(vxk, dt))
            py[k] = r_add(py[k], ex_mul(vyk, dt))
            pz[k] = r_add(pz[k], ex_mul(vzk, dt))
            im = self._inverse_mass[li]
            vx[k] = r_add(vxk, ex_mul(ex_mul(fxk, im), dt))
            vy[k] = r_add(vyk, ex_mul(ex_mul(fyk, im), dt))
            vz[k] = r_add(vzk, ex_mul(ex_mul(fzk, im), dt))

        elapsed_ms = int((_time.monotonic() - self._started) * 1000)
        self._frame_count += 1
        return self.snapshot(elapsed_ms)

    # --- output -------------------------------------------------------------

    def snapshot(self, execution_time_ms: int = 0) -> Frame:
        px, py, pz = self.px, self.py, self.pz
        vx, vy, vz = self.vx, self.vy, self.vz
        types = self.types
        positions = np.array([px, py, pz], dtype=float).T
        particles = ",".join(
            '{"id":%d,"position":{"x":%s,"y":%s,"z":%s},'
            '"velocity":{"x":%s,"y":%s,"z":%s},"type":"%s"}'
            % (k, px[k], py[k], pz[k], vx[k], vy[k], vz[k], types[k])
            for k in range(len(px))
        )
        text = '{"particles":[%s],"executionTime":%d,"time":%s}' % (particles, execution_time_ms, self.time)
        return Frame(time=self.time, execution_time_ms=execution_time_ms, positions=positions, json=text)


ProgressCallback = Callable[[int, int], None]


def run_simulation(
    props: SimulatorProperties,
    progress: Optional[ProgressCallback] = None,
    cancelled: Optional[Callable[[], bool]] = None,
) -> SimulationResult:
    """Run a full simulation: the initial lattice plus one frame per iteration."""
    sim = Simulator(props)
    frames = [sim.snapshot(0)]
    total = props.iterations
    for ii in range(total):
        if cancelled is not None and cancelled():
            raise SimulationCancelled()
        frames.append(sim.step(ii))
        if progress is not None:
            progress(ii + 1, total)
    return SimulationResult(
        properties=props,
        types=list(sim.types),
        lattice_index=lattice_indices(props.lattice_x_size, props.lattice_y_size, props.lattice_z_size),
        frames=frames,
        fusion_events=sim.fusion_events,
    )
