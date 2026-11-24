"""Convert between the parameter form (strings, as typed) and SimulatorProperties."""

from __future__ import annotations

from decimal import Decimal, InvalidOperation
from typing import Any, Optional

from .constants import DEFAULT_THRESHOLD
from .model import SimulatorProperties

MAX_PARTICLES = 100_000


class FormError(ValueError):
    def __init__(self, errors: list[str]):
        super().__init__("\n".join(errors))
        self.errors = errors


def form_from_properties(props: SimulatorProperties) -> dict[str, Any]:
    field = props.electric_field
    return {
        "name": props.name,
        "iterations": str(props.iterations),
        "timeStep": str(props.time_step),
        "precision": str(props.precision),
        "sizeX": str(props.lattice_x_size),
        "sizeY": str(props.lattice_y_size),
        "sizeZ": str(props.lattice_z_size),
        "fieldX": "" if field is None else str(field[0]),
        "fieldY": "" if field is None else str(field[1]),
        "fieldZ": "" if field is None else str(field[2]),
        # Empty means "the original default" (the exact double 2.0E-15).
        "threshold": "" if props.threshold == DEFAULT_THRESHOLD else str(props.threshold),
        "lorentz": props.include_lorentz,
        "logTime": props.use_log_time,
        "transition": str(props.transition_step),
        "damping": "" if props.damping_factor is None else str(props.damping_factor),
    }


def default_form() -> dict[str, Any]:
    return form_from_properties(SimulatorProperties())


def properties_from_form(form: dict[str, Any]) -> SimulatorProperties:
    """Parse and validate. Raises :class:`FormError` listing every problem."""
    errors: list[str] = []

    def text(key: str) -> str:
        value = form.get(key, "")
        return "" if value is None else str(value).strip()

    def integer(key: str, label: str, lo: int, hi: int, default: int) -> int:
        raw = text(key)
        if not raw:
            errors.append(f"{label}: enter a whole number.")
            return default
        try:
            value = int(raw)
        except ValueError:
            errors.append(f"{label}: '{raw}' is not a whole number.")
            return default
        if not lo <= value <= hi:
            errors.append(f"{label} must be between {lo:,} and {hi:,}.")
        return value

    def number(key: str, label: str, default: Optional[Decimal]) -> Optional[Decimal]:
        raw = text(key)
        if not raw:
            return default
        try:
            value = Decimal(raw)
        except InvalidOperation:
            errors.append(f"{label}: '{raw}' is not a number.")
            return default
        if not value.is_finite():
            errors.append(f"{label}: '{raw}' is not a finite number.")
            return default
        return value

    def flag(key: str) -> bool:
        value = form.get(key, False)
        return value is True or str(value).lower() in ("true", "on", "1", "yes")

    iterations = integer("iterations", "Iterations", 0, 10_000_000, 5)
    precision = integer("precision", "Precision", 1, 200, 6)
    sx = integer("sizeX", "Lattice X", 1, 1000, 1)
    sy = integer("sizeY", "Lattice Y", 1, 1000, 1)
    sz = integer("sizeZ", "Lattice Z", 1, 1000, 1)
    if sx * sy * sz > MAX_PARTICLES:
        errors.append(f"The lattice has {sx * sy * sz:,} particles; the limit is {MAX_PARTICLES:,}.")
    time_step = number("timeStep", "Time step", None)
    if time_step is None and not text("timeStep"):
        errors.append("Time step: enter a value, e.g. 2.0E-17.")
    components = [number(f"field{a}", f"Electric field {a}", None) for a in "XYZ"]
    electric_field = None
    if any(c is not None for c in components):
        electric_field = tuple(c if c is not None else Decimal(0) for c in components)
    threshold = number("threshold", "Fusion threshold", DEFAULT_THRESHOLD)
    use_log_time = flag("logTime")
    transition = integer("transition", "Transition step", -10_000_000, 10_000_000, 0) if use_log_time else 0
    damping = number("damping", "Damping factor", None) if use_log_time else None
    if errors:
        raise FormError(errors)

    props = SimulatorProperties(
        name=text("name") or "Simulation",
        lattice_x_size=sx,
        lattice_y_size=sy,
        lattice_z_size=sz,
        iterations=iterations,
        time_step=time_step,
        precision=precision,
        electric_field=electric_field,
        threshold=threshold,
        include_lorentz=flag("lorentz"),
        use_log_time=use_log_time,
        transition_step=transition,
        damping_factor=damping,
    )
    try:
        props.validate()
    except ValueError as exc:
        raise FormError(str(exc).splitlines()) from None
    return props


# --- improved model -----------------------------------------------------------

IMPROVED_FIELDS = {
    # form name: (attribute, label, kind)
    "impShape": ("shape", "Crystal shape", "choice"),
    "impSize": ("size", "Crystal size", "float"),
    "impLithium": ("lithium", "Lithium isotope", "choice"),
    "impIntensity": ("intensity_w_cm2", "Peak intensity", "float"),
    "impWavelength": ("wavelength_um", "Wavelength", "float"),
    "impFwhm": ("fwhm_fs", "Pulse FWHM", "float"),
    "impPhase": ("phase_deg", "Phase at stripping", "float"),
    "impCentre": ("centre_fs", "Pulse centre", "float"),
    "impDuration": ("duration_fs", "Simulated time", "float"),
    "impFrames": ("frames", "Frames", "int"),
    "impTarget": ("target_radius_nm", "Extrapolate to radius", "float"),
    "impProfile": ("profile", "Density profile", "choice"),
    "impCore": ("core_density", "Core density", "float"),
    "impExpansion": ("expansion", "Edge expansion", "float"),
    "impGeometry": ("geometry", "Geometry", "choice"),
    "impGap": ("gap_nm", "Gap between crystals", "float"),
}


def form_from_improved(props) -> dict[str, Any]:
    out: dict[str, Any] = {"model": "improved", "name": props.name}
    for key, (attr, _, kind) in IMPROVED_FIELDS.items():
        value = getattr(props, attr)
        out[key] = value if kind == "choice" else f"{value:g}"
    return out


def improved_from_form(form: dict[str, Any]):
    from .improved import ImprovedProperties

    props = ImprovedProperties()
    errors: list[str] = []
    name = str(form.get("name") or "").strip()
    if name:
        props.name = name
    for key, (attr, label, kind) in IMPROVED_FIELDS.items():
        raw = str(form.get(key, "")).strip()
        if not raw:
            continue  # keep the default
        try:
            if kind == "choice":
                value: Any = raw
            elif kind == "int":
                value = int(raw)
            else:
                value = float(raw)
                if value != value or value in (float("inf"), float("-inf")):
                    raise ValueError
        except ValueError:
            errors.append(f"{label}: '{raw}' is not a valid number.")
            continue
        setattr(props, attr, value)
    if errors:
        raise FormError(errors)
    try:
        props.validate()
    except ValueError as exc:
        raise FormError(str(exc).splitlines()) from None
    return props


def all_defaults() -> dict[str, Any]:
    from .improved import ImprovedProperties

    improved = form_from_improved(ImprovedProperties(intensity_w_cm2=1e20))
    improved.pop("name")
    return {**default_form(), **improved, "model": "original"}
