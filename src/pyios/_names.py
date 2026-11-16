"""Flexible parsing of equation names, parameter names and sex codes."""
from __future__ import annotations

import re

import numpy as np
import pandas as pd

_DASHES = "\u2010\u2011\u2012\u2013\u2014\u2212"

_EQUATION_ALIASES = {"shulz": "schulz", "gochicoarangel": "gochicoa"}


def normalise_equation_key(name) -> str:
    """'Gochicoa_2023', 'gochicoa 2023', 'GOCHICOA-2023' -> 'gochicoa2023'."""
    key = re.sub(r"[^a-z0-9]", "", str(name).lower())
    for wrong, right in _EQUATION_ALIASES.items():
        if key.startswith(wrong):
            key = right + key[len(wrong):]
    return key


def resolve_equation(name, available):
    """Match a user-supplied equation name to one of the available reference ids."""
    key = normalise_equation_key(name)
    if key in available:
        return key
    matches = [a for a in available if a.startswith(key)] if key else []
    if len(matches) == 1:
        return matches[0]
    opts = ", ".join(sorted(available))
    if len(matches) > 1:
        raise ValueError(f"Equation {name!r} is ambiguous (matches {', '.join(matches)}).")
    raise ValueError(f"Unknown equation {name!r}. Available equations: {opts}.")


def normalise_param(param) -> str:
    """Return an upper-case canonical parameter name for matching.

    'r5', 'R5', 'R5-20', 'r5_r20', 'R5\u2013R20' -> 'R5-R20';
    '(R5-R20)/R20', 'R5-R20/R20' -> '(R5-R20)/R20'.
    """
    s = str(param).strip()
    for d in _DASHES:
        s = s.replace(d, "-")
    s = re.sub(r"\s+", "", s).replace("_", "-").upper()
    m = re.fullmatch(r"R(\d+)-(\d+)", s)
    if m:
        s = f"R{m.group(1)}-R{m.group(2)}"
    m = re.fullmatch(r"\(?R5-R(15|20)\)?/R(5|20)(\*100)?", s)
    if m:
        s = f"(R5-R{m.group(1)})/R{m.group(2)}"
    return s


_FEMALE = {"0", "f", "female", "w", "woman", "girl"}
_MALE = {"1", "m", "male", "man", "boy"}


def normalise_sex(value):
    """Return 'Male', 'Female' or None (missing / unrecognised)."""
    if value is None:
        return None
    if isinstance(value, (bool, np.bool_)):
        return None
    try:
        if pd.isna(value):
            return None
    except (TypeError, ValueError):
        return None
    if isinstance(value, (int, float, np.integer, np.floating)):
        if value == 0:
            return "Female"
        if value == 1:
            return "Male"
        return None
    s = str(value).strip().lower()
    if re.fullmatch(r"-?\d+\.0+", s):
        s = s.split(".")[0]
    if s in _FEMALE:
        return "Female"
    if s in _MALE:
        return "Male"
    return None


def normalise_sex_array(values) -> np.ndarray:
    """Vectorised :func:`normalise_sex` returning an object array."""
    cache = {}
    out = np.empty(len(values), dtype=object)
    for i, v in enumerate(values):
        try:
            hit = cache[v]
        except KeyError:
            hit = cache[v] = normalise_sex(v)
        except TypeError:
            hit = normalise_sex(v)
        out[i] = hit
    return out
