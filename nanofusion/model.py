"""Data model: simulation parameters and results."""

from __future__ import annotations

from dataclasses import dataclass, field
from decimal import Decimal
from typing import Optional

import numpy as np

from .constants import DEFAULT_THRESHOLD, INITIAL_LATTICE_SPACING, LITHIUM

Vector = tuple[Decimal, Decimal, Decimal]


@dataclass
class SimulatorProperties:
    """Parameters of one simulation run (Java ``SimulatorProperties``)."""

    name: str = "Simulation 1"
    lattice_x_size: int = 10
    lattice_y_size: int = 10
    lattice_z_size: int = 10
    iterations: int = 5
    time_step: Decimal = Decimal("2.0E-17")
    precision: int = 6
    electric_field: Optional[Vector] = None
    threshold: Decimal = DEFAULT_THRESHOLD
    include_lorentz: bool = False
    use_log_time: bool = False
    transition_step: int = 0
    damping_factor: Optional[Decimal] = None

    @property
    def particle_count(self) -> int:
        return self.lattice_x_size * self.lattice_y_size * self.lattice_z_size

    def validate(self) -> None:
        """Raise ``ValueError`` describing every invalid setting."""
        problems = []
        for label, value in (
            ("Lattice X size", self.lattice_x_size),
            ("Lattice Y size", self.lattice_y_size),
            ("Lattice Z size", self.lattice_z_size),
        ):
            if value < 1:
                problems.append(f"{label} must be at least 1.")
        if self.iterations < 0:
            problems.append("Iterations cannot be negative.")
        if self.precision < 1:
            problems.append("Precision must be at least 1 significant digit.")
        if not self.time_step.is_finite() or self.time_step <= 0:
            problems.append("Time step must be a positive number.")
        if not self.threshold.is_finite() or self.threshold < 0:
            problems.append("Threshold must be a non-negative number.")
        if self.electric_field is not None and not all(c.is_finite() for c in self.electric_field):
            problems.append("Electric field components must be numbers.")
        if self.use_log_time:
            if self.damping_factor is None or not self.damping_factor.is_finite() or self.damping_factor == 0:
                problems.append("Log time needs a non-zero damping factor.")
        if problems:
            raise ValueError("\n".join(problems))

    # --- JSON (same field names and order as the Java/Gson output) ---------

    def to_json_fields(self) -> list[tuple[str, object]]:
        fields: list[tuple[str, object]] = [
            ("name", self.name),
            ("latticeXSize", self.lattice_x_size),
            ("latticeYSize", self.lattice_y_size),
            ("latticeZSize", self.lattice_z_size),
            ("iterations", self.iterations),
            ("timeStep", self.time_step),
            ("precision", self.precision),
        ]
        if self.electric_field is not None:
            ex, ey, ez = self.electric_field
            fields.append(("electricField", {"x": ex, "y": ey, "z": ez}))
        fields += [
            ("threshold", self.threshold),
            ("includeLorentz", self.include_lorentz),
            ("useLogTime", self.use_log_time),
            ("transitionStep", self.transition_step),
        ]
        if self.damping_factor is not None:
            fields.append(("dampingFactor", self.damping_factor))
        return fields

    @classmethod
    def from_json(cls, data: dict) -> "SimulatorProperties":
        """Build from a parsed ``properties`` object (numbers as ``Decimal``).

        Older saved files omit fields that were added later; those keep their
        defaults.
        """
        defaults = cls()

        def dec(value, default):
            return default if value is None else Decimal(str(value))

        field_data = data.get("electricField")
        electric_field = None
        if field_data is not None:
            electric_field = tuple(Decimal(str(field_data.get(k, 0))) for k in ("x", "y", "z"))
        damping = data.get("dampingFactor")
        return cls(
            name=str(data.get("name", defaults.name)),
            lattice_x_size=int(data.get("latticeXSize", defaults.lattice_x_size)),
            lattice_y_size=int(data.get("latticeYSize", defaults.lattice_y_size)),
            lattice_z_size=int(data.get("latticeZSize", defaults.lattice_z_size)),
            iterations=int(data.get("iterations", defaults.iterations)),
            time_step=dec(data.get("timeStep"), defaults.time_step),
            precision=int(data.get("precision", defaults.precision)),
            electric_field=electric_field,
            threshold=dec(data.get("threshold"), defaults.threshold),
            include_lorentz=bool(data.get("includeLorentz", False)),
            use_log_time=bool(data.get("useLogTime", False)),
            transition_step=int(data.get("transitionStep", 0)),
            damping_factor=None if damping is None else Decimal(str(damping)),
        )


@dataclass
class Frame:
    """One lattice snapshot (Java ``Lattice``).

    ``positions`` holds float copies for plotting. ``json`` is the exact
    serialised lattice (Decimal digits as computed) when the frame came from
    the engine; it is ``None`` for frames read from a file.
    """

    time: Decimal
    execution_time_ms: int
    positions: np.ndarray  # shape (N, 3), metres
    json: Optional[str] = None


@dataclass(frozen=True)
class FusionEvent:
    """Two particles came within the fusion threshold during a step."""

    frame: int  # index into SimulationResult.frames of the positions that were measured
    particle_1: int
    particle_2: int
    distance: Decimal


@dataclass
class SimulationResult:
    properties: SimulatorProperties
    types: list[str]  # "D" or "Li" per particle id
    lattice_index: np.ndarray  # shape (N, 3) ints: row (x), column (y), layer (z)
    frames: list[Frame] = field(default_factory=list)
    fusion_events: Optional[list[FusionEvent]] = None  # None when unknown (loaded from file)
    source_path: Optional[str] = None
    saved: bool = False
    warning: Optional[str] = None  # e.g. the file was truncated
    report: Optional[dict] = None  # improved-model report (see nanofusion.improved)

    @property
    def is_lithium(self) -> np.ndarray:
        return np.array([t == LITHIUM for t in self.types], dtype=bool)

    @property
    def fused(self) -> Optional[bool]:
        if self.fusion_events is None:
            return None
        return len(self.fusion_events) > 0

    @property
    def can_save(self) -> bool:
        return bool(self.frames) and all(f.json is not None for f in self.frames)

    def bounds(self, pad: float | None = None) -> tuple[np.ndarray, np.ndarray]:
        """Min/max position over all frames, padded (default: one lattice spacing)."""
        if pad is None:
            pad = float(INITIAL_LATTICE_SPACING)
        lo = np.full(3, np.inf)
        hi = np.full(3, -np.inf)
        for frame in self.frames:
            lo = np.minimum(lo, frame.positions.min(axis=0))
            hi = np.maximum(hi, frame.positions.max(axis=0))
        return lo - pad, hi + pad


def lattice_indices(x_size: int, y_size: int, z_size: int) -> np.ndarray:
    """Lattice coordinates for particle ids in Java construction order."""
    ids = np.arange(x_size * y_size * z_size)
    return np.stack([ids % x_size, (ids // x_size) % y_size, ids // (x_size * y_size)], axis=1)
