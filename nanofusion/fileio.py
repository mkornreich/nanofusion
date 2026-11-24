"""Saving and loading simulations.

Files written here are standard JSON with one lattice per line:

    {
    "properties": {...},
    "results": [
    {"particles": [...], "executionTime": 12, "time": 2.0E-17},
    ...
    ]
    }

The loader also reads the files written by the original Java application,
which used unquoted keys and a trailing comma (``properties: {...}``).
"""

from __future__ import annotations

import json
import os
import re
from decimal import Decimal
from pathlib import Path
from typing import Callable, Optional, TextIO, Union

import numpy as np

from .model import Frame, SimulationResult, SimulatorProperties, lattice_indices

PathLike = Union[str, os.PathLike]

_WS = re.compile(r"[\s,]*")
_TIME = re.compile(r'"time"\s*:\s*([-+0-9.eE]+)')
_KEY = re.compile(r'\s*(?:"([^"]+)"|([A-Za-z_]\w*))\s*:\s*')


def _encode(value: object) -> str:
    if isinstance(value, bool):
        return "true" if value else "false"
    if isinstance(value, Decimal):
        return str(value)  # same text as Java BigDecimal.toString()
    if isinstance(value, int):
        return str(value)
    if isinstance(value, str):
        return json.dumps(value)
    if isinstance(value, dict):
        return "{" + ",".join(f"{json.dumps(k)}:{_encode(v)}" for k, v in value.items()) + "}"
    raise TypeError(f"cannot encode {type(value).__name__}")


def properties_json(props: SimulatorProperties) -> str:
    return "{" + ",".join(f'"{k}":{_encode(v)}' for k, v in props.to_json_fields()) + "}"


def write_simulation(
    result: SimulationResult,
    out: TextIO,
    java_format: bool = False,
    progress: Optional[Callable[[int, int], None]] = None,
) -> None:
    """Write ``result`` as text to an open stream.

    ``java_format=True`` reproduces the original application's non-standard
    layout byte for byte (unquoted keys, trailing comma).
    """
    if not result.can_save:
        raise ValueError("Only simulations run in this session can be saved.")
    total = len(result.frames)
    out.write("{\n")
    out.write("properties: " if java_format else '"properties": ')
    out.write(properties_json(result.properties))
    out.write(",\n")
    out.write("results:\n[\n" if java_format else '"results": [\n')
    for index, frame in enumerate(result.frames):
        out.write(frame.json)
        if java_format or index < total - 1:
            out.write(",")
        out.write("\n")
        if progress is not None:
            progress(index + 1, total)
    out.write("]\n}" if java_format else "]\n}\n")


def save_simulation(
    result: SimulationResult,
    path: PathLike,
    java_format: bool = False,
    progress: Optional[Callable[[int, int], None]] = None,
) -> None:
    """Write ``result`` to ``path`` (atomically, via a ``.part`` file)."""
    if not result.can_save:
        raise ValueError("Only simulations run in this session can be saved.")
    tmp = Path(str(path) + ".part")
    with open(tmp, "w", encoding="utf-8", newline="\n") as out:
        write_simulation(result, out, java_format=java_format, progress=progress)
    os.replace(tmp, path)
    result.source_path = str(path)
    result.saved = True


def _frame_from_lattice(lattice: dict, raw_tail: str) -> Frame:
    particles = lattice["particles"]
    positions = np.array(
        [(p["position"]["x"], p["position"]["y"], p["position"]["z"]) for p in particles],
        dtype=float,
    )
    # Keep the time exactly as written ("8.0E-17"); json would turn it into a float.
    times = _TIME.findall(raw_tail)
    t = Decimal(times[-1]) if times else Decimal(repr(lattice.get("time", 0)))
    return Frame(
        time=t,
        execution_time_ms=int(lattice.get("executionTime", 0)),
        positions=positions,
    )


def load_simulation(
    path: PathLike,
    progress: Optional[Callable[[int, int], None]] = None,
) -> SimulationResult:
    """Read a simulation saved by this program or by the original Java app."""
    return load_simulation_text(Path(path).read_text(encoding="utf-8"), str(path), progress)


def load_simulation_text(
    text: str,
    source: Optional[str] = None,
    progress: Optional[Callable[[int, int], None]] = None,
) -> SimulationResult:
    """Parse the text of a saved simulation (see :func:`load_simulation`)."""
    decoder = json.JSONDecoder()
    exact = json.JSONDecoder(parse_float=Decimal)
    length = len(text)
    pos = _WS.match(text, 0).end()
    if text[pos : pos + 1] != "{":
        raise ValueError("Not a simulation file: expected '{' at the start.")
    pos += 1

    props_data: Optional[dict] = None
    lattices: list[dict] = []
    frames: list[Frame] = []
    warning: Optional[str] = None
    while warning is None:
        pos = _WS.match(text, pos).end()
        if pos >= length or text[pos] == "}":
            break
        key_match = _KEY.match(text, pos)
        if key_match is None:
            raise ValueError(f"Unexpected content at character {pos}.")
        key = key_match.group(1) or key_match.group(2)
        pos = key_match.end()
        if key == "properties":
            props_data, pos = exact.raw_decode(text, pos)
        elif key == "results":
            if text[pos : pos + 1] != "[":
                raise ValueError("'results' must be a list.")
            pos += 1
            while True:
                pos = _WS.match(text, pos).end()
                if pos < length and text[pos] == "]":
                    pos += 1
                    break
                try:
                    if pos >= length:
                        raise ValueError("end of file")
                    lattice, end = decoder.raw_decode(text, pos)
                except ValueError:
                    if not frames:
                        raise ValueError("The file ends before the first complete frame.") from None
                    warning = (f"The file is incomplete (it ends part-way through frame "
                               f"{len(frames) + 1}). Loaded the {len(frames):,} complete frames.")
                    break
                if not lattices:
                    lattices.append(lattice)  # keep the first one for particle types
                frames.append(_frame_from_lattice(lattice, text[max(pos, end - 120) : end]))
                pos = end
                if progress is not None:
                    progress(pos, length)
        else:
            _, pos = decoder.raw_decode(text, pos)

    if props_data is None:
        raise ValueError("The file has no 'properties' section.")
    props = SimulatorProperties.from_json(props_data)
    if not frames:
        raise ValueError("The file has no results.")
    particles = lattices[0]["particles"]
    types = [p.get("type", "D") for p in particles]
    n = len(types)
    if props.particle_count == n:
        index = lattice_indices(props.lattice_x_size, props.lattice_y_size, props.lattice_z_size)
    else:  # sizes missing or inconsistent: treat as a single row
        index = np.stack([np.arange(n), np.zeros(n, int), np.zeros(n, int)], axis=1)
    for frame in frames:
        if frame.positions.shape != (n, 3):
            raise ValueError("Frames have different particle counts.")
    if progress is not None:
        progress(length, length)
    return SimulationResult(
        properties=props,
        types=types,
        lattice_index=index,
        frames=frames,
        fusion_events=None,
        source_path=source,
        saved=True,
        warning=warning,
    )
