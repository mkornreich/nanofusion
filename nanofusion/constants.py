"""Physical constants, mirroring org.nanofusion.simulation.Constants.

Every value is a :class:`~decimal.Decimal` built exactly as the Java
``BigDecimal`` original was, including the deliberately coarse ``INVERSE_C``
(``1/c`` rounded to 10 decimal places, i.e. ``3.3E-9``).
"""

from __future__ import annotations

import math
from decimal import (
    MAX_EMAX,
    MAX_PREC,
    MIN_EMIN,
    ROUND_DOWN,
    ROUND_HALF_UP,
    Context,
    Decimal,
    DivisionByZero,
    InvalidOperation,
    Overflow,
)

# Context for exact add/subtract/multiply (Java BigDecimal ops without a
# MathContext). Never use it for division: a non-terminating quotient would
# try to produce MAX_PREC digits.
EXACT = Context(
    prec=MAX_PREC,
    rounding=ROUND_HALF_UP,
    Emax=MAX_EMAX,
    Emin=MIN_EMIN,
    traps=[InvalidOperation, DivisionByZero, Overflow],
)


def rounding_context(precision: int) -> Context:
    """Equivalent of Java ``new MathContext(precision, RoundingMode.HALF_UP)``."""
    if precision < 1:
        raise ValueError("precision must be at least 1")
    return Context(
        prec=precision,
        rounding=ROUND_HALF_UP,
        Emax=MAX_EMAX,
        Emin=MIN_EMIN,
        traps=[InvalidOperation, DivisionByZero, Overflow],
    )


def divide_to_scale(dividend: Decimal, divisor: Decimal, scale: int) -> Decimal:
    """Java ``dividend.divide(divisor, scale, RoundingMode.HALF_UP)``.

    Truncating at high precision first and then rounding HALF_UP to the target
    scale gives the correctly rounded result (no double-rounding error).
    """
    wide = Context(prec=300, rounding=ROUND_DOWN, Emax=MAX_EMAX, Emin=MIN_EMIN)
    quotient = wide.divide(dividend, divisor)
    return quotient.quantize(Decimal(1).scaleb(-scale), rounding=ROUND_HALF_UP, context=wide)


_mul = EXACT.multiply

ZERO = Decimal(0)
ONE = Decimal(1)

SCALE = 10
CHARGE = Decimal("1.60217657E-19")
TIME_INTERVAL = Decimal("2.0E-17")
INITIAL_LATTICE_SPACING = Decimal("2.04E-10")
ATOMIC_MASS_UNIT = Decimal("1.660538921E-27")
DEUTERIUM_MASS = _mul(Decimal("2.01410178"), ATOMIC_MASS_UNIT)
INVERSE_DEUTERIUM_MASS = divide_to_scale(ONE, DEUTERIUM_MASS, 10)
DEUTERIUM_CHARGE = ONE
LITHIUM_MASS = _mul(Decimal("7.01600455"), ATOMIC_MASS_UNIT)
INVERSE_LITHIUM_MASS = divide_to_scale(ONE, LITHIUM_MASS, 10)
LITHIUM_CHARGE = Decimal("3.0")
K_E = Decimal("8.9875517873681765E+9")
POTENTIAL = _mul(_mul(K_E, CHARGE), CHARGE)
C = Decimal("3.0E8")
INVERSE_C = divide_to_scale(ONE, C, 10)
LAMBDA = Decimal("3.6E-6")
# new BigDecimal(Math.PI) is the exact binary value of the double, as is Decimal(math.pi).
OMEGA = divide_to_scale(_mul(_mul(Decimal("2.0"), Decimal(math.pi)), C), LAMBDA, 18)

# Default fusion threshold: Java `new BigDecimal(2.0E-15)` (exact double value).
DEFAULT_THRESHOLD = Decimal(2.0e-15)

DEUTERIUM = "D"
LITHIUM = "Li"
