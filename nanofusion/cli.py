"""Command-line simulation runner: ``python -m nanofusion run ...``."""

from __future__ import annotations

import argparse
import sys
import time
from decimal import Decimal, InvalidOperation

from .engine import run_simulation
from .fileio import save_simulation
from .model import SimulatorProperties


def _decimal(text: str) -> Decimal:
    try:
        value = Decimal(text)
    except InvalidOperation:
        raise argparse.ArgumentTypeError(f"not a number: {text!r}") from None
    if not value.is_finite():
        raise argparse.ArgumentTypeError(f"not a finite number: {text!r}")
    return value


def build_parser() -> argparse.ArgumentParser:
    defaults = SimulatorProperties()
    p = argparse.ArgumentParser(
        prog="nanofusion run",
        description="Run a lattice fusion simulation without the GUI and save it as JSON.",
    )
    p.add_argument("-o", "--output", help="file to write (JSON); omit to only print a summary")
    p.add_argument("--name", default=defaults.name)
    p.add_argument("--size", nargs=3, type=int, metavar=("X", "Y", "Z"),
                   default=[defaults.lattice_x_size, defaults.lattice_y_size, defaults.lattice_z_size],
                   help="lattice size in particles (default: 10 10 10)")
    p.add_argument("-n", "--iterations", type=int, default=defaults.iterations)
    p.add_argument("--precision", type=int, default=defaults.precision,
                   help="significant digits kept by rounded operations (default: 6)")
    p.add_argument("--time-step", type=_decimal, default=defaults.time_step, help="seconds (default: 2.0E-17)")
    p.add_argument("--efield", nargs=3, type=_decimal, metavar=("EX", "EY", "EZ"),
                   help="peak electric field in V/m, e.g. 1e12 0 0")
    p.add_argument("--threshold", type=_decimal, default=defaults.threshold,
                   help="fusion distance threshold in metres (default: 2.0E-15)")
    p.add_argument("--lorentz", action="store_true", help="include the Lorentz force")
    p.add_argument("--log-time", action="store_true", help="shrink the time step logistically")
    p.add_argument("--transition", type=int, default=0, help="log time: step where the time step halves")
    p.add_argument("--damping", type=_decimal, help="log time: damping factor (steps)")
    p.add_argument("--java-format", action="store_true",
                   help="write the original Java app's non-standard file layout")
    return p


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)
    props = SimulatorProperties(
        name=args.name,
        lattice_x_size=args.size[0],
        lattice_y_size=args.size[1],
        lattice_z_size=args.size[2],
        iterations=args.iterations,
        time_step=args.time_step,
        precision=args.precision,
        electric_field=tuple(args.efield) if args.efield else None,
        threshold=args.threshold,
        include_lorentz=args.lorentz,
        use_log_time=args.log_time,
        transition_step=args.transition,
        damping_factor=args.damping,
    )
    try:
        props.validate()
    except ValueError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2

    interactive = sys.stderr.isatty()

    def progress(done: int, total: int) -> None:
        if interactive:
            print(f"\rstep {done}/{total}", end="", file=sys.stderr, flush=True)

    start = time.monotonic()
    result = run_simulation(props, progress=progress)
    if interactive:
        print(file=sys.stderr)
    elapsed = time.monotonic() - start
    events = result.fusion_events or []
    print(f"{props.particle_count} particles, {props.iterations} steps in {elapsed:.2f} s")
    print(f"simulated time: {result.frames[-1].time} s")
    if events:
        first = events[0]
        print(f"fusion: YES ({len(events)} close approaches; first in frame {first.frame + 1}, "
              f"particles {first.particle_1} and {first.particle_2}, {first.distance:.3e} m)")
    else:
        print("fusion: NO")
    if args.output:
        save_simulation(result, args.output, java_format=args.java_format)
        print(f"saved {args.output}")
    return 0
