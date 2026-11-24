"""The Python engine must reproduce the original Java engine digit for digit.

Reference outputs in tests/fixtures/ come from the original Java code
(tests/java_reference/generate_fixtures.py).
"""

from __future__ import annotations

import gzip
import json
import re
import shlex
from decimal import Decimal
from pathlib import Path

import pytest

from nanofusion import constants as K
from nanofusion.engine import SimulationCancelled, Simulator, detect_fusion, run_simulation
from nanofusion.fileio import load_simulation, save_simulation
from nanofusion.model import SimulatorProperties

FIXTURES = Path(__file__).parent / "fixtures"
EXPECTED = json.loads((FIXTURES / "expected.json").read_text())
_EXECUTION_TIME = re.compile(r'"executionTime":\d+')


def props_from_args(spec: str) -> SimulatorProperties:
    """Same key=value parameters as tests/java_reference/ReferenceRun.java."""
    args: dict[str, str] = {}
    last = None
    for token in shlex.split(spec):
        if "=" in token:
            last, value = token.split("=", 1)
            args[last] = value
        else:
            args[last] += " " + token
    p = SimulatorProperties()
    p.name = args.get("name", p.name)
    p.lattice_x_size = int(args.get("x", p.lattice_x_size))
    p.lattice_y_size = int(args.get("y", p.lattice_y_size))
    p.lattice_z_size = int(args.get("z", p.lattice_z_size))
    p.iterations = int(args.get("iterations", p.iterations))
    p.precision = int(args.get("precision", p.precision))
    if "timeStep" in args:
        p.time_step = Decimal(args["timeStep"])
    if "threshold" in args:
        p.threshold = Decimal(args["threshold"])
    if any(k in args for k in ("ex", "ey", "ez")):
        p.electric_field = tuple(Decimal(args.get(k, "0")) for k in ("ex", "ey", "ez"))
    p.include_lorentz = args.get("lorentz") == "true"
    p.use_log_time = args.get("logTime") == "true"
    p.transition_step = int(args.get("transition", 0))
    if "damping" in args:
        p.damping_factor = Decimal(args["damping"])
    return p


def first_difference(expected: str, actual: str) -> str:
    for line_no, (a, b) in enumerate(zip(expected.splitlines(), actual.splitlines()), 1):
        if a != b:
            ta, tb = re.split(r"[,{}\[\]]", a), re.split(r"[,{}\[\]]", b)
            for x, y in zip(ta, tb):
                if x != y:
                    return f"line {line_no}: Java {x!r} vs Python {y!r}"
            return f"line {line_no} differs in length"
    return "different number of lines"


@pytest.mark.parametrize("case", sorted(EXPECTED))
def test_matches_java_engine(case: str, tmp_path: Path) -> None:
    props = props_from_args(EXPECTED[case]["args"])
    result = run_simulation(props)
    out = tmp_path / "out.json"
    save_simulation(result, out, java_format=True)
    with gzip.open(FIXTURES / f"{case}.json.gz", "rt", encoding="utf-8") as f:
        java = _EXECUTION_TIME.sub("", f.read())
    ours = _EXECUTION_TIME.sub("", out.read_text(encoding="utf-8"))
    assert ours == java, first_difference(java, ours)
    assert result.fused == EXPECTED[case]["fusion"]


def test_constants() -> None:
    assert K.INVERSE_C == Decimal("3.3E-9")  # 1/c rounded to 10 decimal places, as in Java
    assert K.DEUTERIUM_MASS == Decimal("2.01410178") * Decimal("1.660538921E-27")
    assert abs(float(K.INVERSE_DEUTERIUM_MASS) * float(K.DEUTERIUM_MASS) - 1) < 1e-15
    assert K.INVERSE_DEUTERIUM_MASS.as_tuple().exponent == -10
    assert float(K.OMEGA) == pytest.approx(2 * 3.141592653589793 * 3e8 / 3.6e-6)
    assert K.OMEGA.as_tuple().exponent == -18


def test_lattice_layout() -> None:
    sim = Simulator(SimulatorProperties(lattice_x_size=2, lattice_y_size=2, lattice_z_size=2, iterations=0))
    assert sim.types == ["D", "Li", "Li", "D", "Li", "D", "D", "Li"]
    assert (sim.px[1], sim.py[2], sim.pz[4]) == (K.INITIAL_LATTICE_SPACING,) * 3


def test_cancel() -> None:
    with pytest.raises(SimulationCancelled):
        run_simulation(SimulatorProperties(lattice_x_size=2, lattice_y_size=2, lattice_z_size=2,
                                           iterations=10), cancelled=lambda: True)


@pytest.mark.parametrize(
    "change, message",
    [
        ({"precision": 0}, "Precision"),
        ({"lattice_x_size": 0}, "Lattice X"),
        ({"use_log_time": True}, "damping"),
        ({"use_log_time": True, "damping_factor": Decimal(0)}, "damping"),
        ({"time_step": Decimal(0)}, "Time step"),
    ],
)
def test_validation(change: dict, message: str) -> None:
    props = SimulatorProperties(**change)
    with pytest.raises(ValueError, match=message):
        props.validate()


def test_lorentz_without_field_does_not_crash() -> None:
    # The Java engine threw NullPointerException here.
    result = run_simulation(SimulatorProperties(lattice_x_size=3, lattice_y_size=2, lattice_z_size=1,
                                                iterations=5, include_lorentz=True))
    assert len(result.frames) == 6


def test_detect_fusion_on_saved_file_matches_engine(tmp_path: Path) -> None:
    props = props_from_args(EXPECTED["fusion_threshold"]["args"])
    result = run_simulation(props)
    save_simulation(result, tmp_path / "f.json")
    loaded = load_simulation(tmp_path / "f.json")
    found = detect_fusion(loaded)
    key = lambda e: (e.frame, e.particle_1, e.particle_2)  # noqa: E731
    assert sorted(map(key, found)) == sorted(map(key, result.fusion_events))
