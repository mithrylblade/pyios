"""Unit handling for inputs (age, height, weight) and oscillometry results."""
from __future__ import annotations

import re

# canonical unit -> multiplier to the package's internal standard unit
_AGE = {"years": 1.0, "months": 1.0 / 12.0, "weeks": 7.0 / 365.25, "days": 1.0 / 365.25}
_AGE_ALIASES = {
    "y": "years", "yr": "years", "yrs": "years", "year": "years", "years": "years",
    "mo": "months", "mos": "months", "month": "months", "months": "months",
    "w": "weeks", "wk": "weeks", "wks": "weeks", "week": "weeks", "weeks": "weeks",
    "d": "days", "day": "days", "days": "days",
}

_HEIGHT = {"cm": 1.0, "m": 100.0, "mm": 0.1, "in": 2.54}
_HEIGHT_ALIASES = {
    "cm": "cm", "centimetre": "cm", "centimetres": "cm", "centimeter": "cm", "centimeters": "cm",
    "m": "m", "metre": "m", "metres": "m", "meter": "m", "meters": "m",
    "mm": "mm", "millimetre": "mm", "millimetres": "mm", "millimeter": "mm", "millimeters": "mm",
    "in": "in", "inch": "in", "inches": "in",
}

_WEIGHT = {"kg": 1.0, "g": 0.001, "lb": 0.45359237}
_WEIGHT_ALIASES = {
    "kg": "kg", "kgs": "kg", "kilogram": "kg", "kilograms": "kg", "kilo": "kg", "kilos": "kg",
    "g": "g", "gram": "g", "grams": "g",
    "lb": "lb", "lbs": "lb", "pound": "lb", "pounds": "lb",
}

# pressure units -> kPa
_PRESSURE = {"kPa": 1.0, "hPa": 0.1, "cmH2O": 0.0980665}
_PRESSURE_ALIASES = {
    "kpa": "kPa", "hpa": "hPa", "mbar": "hPa", "millibar": "hPa",
    "cmh2o": "cmH2O", "cmh\u2082o": "cmH2O",
}


def _key(value) -> str:
    return re.sub(r"[^a-z0-9\u2082]", "", str(value).strip().lower())


def _lookup(value, aliases, table, kind):
    canon = aliases.get(_key(value))
    if canon is None:
        raise ValueError(
            f"Unsupported {kind} unit {value!r}. Supported: {', '.join(sorted(table))}."
        )
    return canon, table[canon]


def parse_age_unit(value):
    """Return (canonical name, factor to years)."""
    return _lookup(value, _AGE_ALIASES, _AGE, "age")


def parse_height_unit(value):
    """Return (canonical name, factor to cm)."""
    return _lookup(value, _HEIGHT_ALIASES, _HEIGHT, "height")


def parse_weight_unit(value):
    """Return (canonical name, factor to kg)."""
    return _lookup(value, _WEIGHT_ALIASES, _WEIGHT, "weight")


def parse_result_unit(value):
    """Return (canonical name, factor to kPa) for the pressure part of a result unit."""
    return _lookup(value, _PRESSURE_ALIASES, _PRESSURE, "result")


_REF_PRESSURE = re.compile(r"^(kPa|hPa|cmH2O)")


def reference_pressure_factor(unit: str):
    """kPa-per-unit factor for a reference equation unit string such as ``hPa*s/L``.

    Returns ``None`` when the quantity is not pressure-based (``Hz``, ``%``).
    """
    m = _REF_PRESSURE.match(str(unit))
    return _PRESSURE[m.group(1)] if m else None


def reference_pressure_name(unit: str):
    m = _REF_PRESSURE.match(str(unit))
    return m.group(1) if m else None
