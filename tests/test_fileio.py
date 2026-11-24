from __future__ import annotations

import gzip
import json
from decimal import Decimal
from pathlib import Path

import numpy as np
import pytest

from nanofusion.__main__ import main as entry_point
from nanofusion.engine import run_simulation
from nanofusion.fileio import load_simulation, save_simulation
from nanofusion.model import SimulatorProperties

FIXTURES = Path(__file__).parent / "fixtures"
JAVA_SAMPLES = Path(__file__).parents[1] / "java" / "FusionSimulation"


def small_result():
    props = SimulatorProperties(name='Test "quoted" run', lattice_x_size=3, lattice_y_size=2,
                                lattice_z_size=2, iterations=6,
                                electric_field=(Decimal("1e12"), Decimal(0), Decimal("-2e11")),
                                include_lorentz=True)
    return run_simulation(props)


def test_saved_file_is_standard_json_and_round_trips(tmp_path: Path) -> None:
    result = small_result()
    path = tmp_path / "sim.json"
    save_simulation(result, path)
    data = json.loads(path.read_text())  # strict JSON
    assert set(data) == {"properties", "results"}
    assert len(data["results"]) == 7
    assert data["properties"]["name"] == 'Test "quoted" run'

    loaded = load_simulation(path)
    assert loaded.properties == result.properties
    assert loaded.types == result.types
    assert [f.time for f in loaded.frames] == [f.time for f in result.frames]
    for a, b in zip(loaded.frames, result.frames):
        np.testing.assert_array_equal(a.positions, b.positions)
    assert loaded.warning is None


def test_loads_original_java_format(tmp_path: Path) -> None:
    path = tmp_path / "java.json"
    with gzip.open(FIXTURES / "basic_4x4x4.json.gz", "rb") as f:
        path.write_bytes(f.read())
    loaded = load_simulation(path)
    assert len(loaded.frames) == 5
    assert loaded.properties.electric_field == (Decimal("1.0E+12"), Decimal(0), Decimal(0))
    assert str(loaded.frames[-1].time) == "8.0E-17"
    assert loaded.lattice_index[17].tolist() == [1, 0, 1]


def test_truncated_file_keeps_complete_frames(tmp_path: Path) -> None:
    result = small_result()
    path = tmp_path / "sim.json"
    save_simulation(result, path)
    text = path.read_text()
    cut = text.rfind('"particles"', 0, len(text) - 200) + 300  # inside the last-but-one frame
    path.write_text(text[:cut])
    loaded = load_simulation(path)
    assert 0 < len(loaded.frames) < 7
    assert "incomplete" in loaded.warning


def test_rejects_non_simulation(tmp_path: Path) -> None:
    path = tmp_path / "x.json"
    path.write_text("[1, 2, 3]")
    with pytest.raises(ValueError):
        load_simulation(path)


@pytest.mark.skipif(not (JAVA_SAMPLES / "simsample.json").exists(), reason="original samples not extracted")
def test_loads_original_sample_file() -> None:
    loaded = load_simulation(JAVA_SAMPLES / "simsample.json")
    assert len(loaded.frames) == 5 and len(loaded.types) == 64


def test_cli_run(tmp_path: Path, capsys: pytest.CaptureFixture) -> None:
    out = tmp_path / "cli.json"
    code = entry_point(["run", "--size", "2", "2", "2", "-n", "3", "--efield", "1e12", "0", "0", "-o", str(out)])
    assert code == 0
    assert "fusion: NO" in capsys.readouterr().out
    assert len(load_simulation(out).frames) == 4
