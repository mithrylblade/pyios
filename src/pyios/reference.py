"""Predicted values, % predicted, z-scores and limits of normal for oscillometry."""
from __future__ import annotations

import math
from collections import defaultdict
from dataclasses import asdict, dataclass
from typing import Dict, List, Optional

import numpy as np
import pandas as pd

from ._data import Model, Row, load_references
from ._names import normalise_sex_array, resolve_equation
from ._units import (
    parse_age_unit, parse_height_unit, parse_result_unit, parse_weight_unit,
    reference_pressure_name,
)

#: z-score corresponding to the 5th / 95th percentile (as used by Gochicoa et al.)
Z95 = 1.6449

_METRICS = ("predicted", "percent_predicted", "z_score", "lln", "uln")


@dataclass
class OscillometryResult:
    """Result of :meth:`OscillometryReference.calculate`. Unavailable values are ``None``."""
    predicted: Optional[float]
    percent_predicted: Optional[float]
    z_score: Optional[float]
    lln: Optional[float]
    uln: Optional[float]
    flag: Optional[str]
    param: str = ""
    reference: str = ""
    unit: str = ""

    def to_dict(self) -> dict:
        return asdict(self)

    def __getitem__(self, key):
        return getattr(self, key)


class _Flags:
    """Collects per-row flag messages; rows without a message stay ``None``."""

    def __init__(self, n: int):
        self.n = n
        self._msgs = defaultdict(list)

    def add(self, mask, text: str):
        for i in np.flatnonzero(mask):
            self._msgs[int(i)].append(text)

    def add_each(self, mask, fn):
        for i in np.flatnonzero(mask):
            self._msgs[int(i)].append(fn(int(i)))

    def add_all(self, text: str):
        for i in range(self.n):
            self._msgs[i].append(text)

    def result(self) -> np.ndarray:
        out = np.full(self.n, None, dtype=object)
        for i, msgs in self._msgs.items():
            out[i] = "; ".join(msgs)
        return out


def _g(x) -> str:
    return f"{x:.6g}"


def _range_flags(fb: _Flags, label: str, unit: str, values, lo, hi):
    lo = np.broadcast_to(np.asarray(lo, float), values.shape)
    hi = np.broadcast_to(np.asarray(hi, float), values.shape)
    with np.errstate(invalid="ignore"):
        low = values < lo
        high = values > hi
    fb.add_each(low, lambda i: f"{label} {_g(values[i])} {unit} below reference minimum {_g(lo[i])} {unit}")
    fb.add_each(high, lambda i: f"{label} {_g(values[i])} {unit} above reference maximum {_g(hi[i])} {unit}")


# ---------------------------------------------------------------- maths
def _fwd(y, tr: str, off: float):
    with np.errstate(invalid="ignore", divide="ignore"):
        if tr == "ln":
            return np.log(y + off)
        if tr == "ln_neg":
            return np.log(off - y)
    return y


def _inv(t, tr: str, off: float):
    with np.errstate(over="ignore", invalid="ignore"):
        if tr == "ln":
            return np.exp(t) - off
        if tr == "ln_neg":
            return off - np.exp(t)
    return t


def _linear_predictor(row: Row, age, h_cm, w, bmi, male):
    h = h_cm / 100.0 if row.height_unit == "m" else h_cm
    lin = np.full(age.shape, row.intercept, dtype=float)
    if row.sex_coef:
        lin += row.sex_coef * male
    if row.age_coef:
        lin += row.age_coef * (np.log(age) if row.age_transform == "ln" else age)
    if row.height_coef:
        lin += row.height_coef * h
    if row.inv_height_coef:
        lin += row.inv_height_coef / h
    if row.weight_coef:
        lin += row.weight_coef * w
    if row.bmi_coef:
        lin += row.bmi_coef * bmi
    if row.inv_bmi_coef:
        lin += row.inv_bmi_coef / bmi
    return lin


class OscillometryReference:
    """Reference equations for impulse oscillometry.

    Parameters
    ----------
    equation : str
        One of ``moitra2020``, ``oostveen2013``, ``schulz2013``, ``berger2021``,
        ``gochicoa2023``. Case, spaces, underscores and the year are optional
        (``"Gochicoa"``, ``"gochicoa_2023"``, ``"GOCHICOA 2023"`` all work).
    age_unit : {"years", "months", "weeks", "days"}
    height_unit : {"cm", "m", "mm", "in"}
    weight_unit : {"kg", "g", "lb"}
    result_unit : {"kPa", "hPa", "cmH2O"}
        Pressure unit of measured values given to the methods and of all returned
        results (e.g. ``kPa`` gives resistance in kPa*s/L).
    """

    def __init__(self, equation: str, age_unit: str = "years", height_unit: str = "cm",
                 weight_unit: str = "kg", result_unit: str = "kPa"):
        refs = load_references()
        self.equation = resolve_equation(equation, list(refs))
        self._ref = refs[self.equation]
        self.age_unit, self._age_factor = parse_age_unit(age_unit)
        self.height_unit, self._height_factor = parse_height_unit(height_unit)
        self.weight_unit, self._weight_factor = parse_weight_unit(weight_unit)
        self.result_unit, self._result_factor = parse_result_unit(result_unit)

    # ---------------------------------------------------------- introspection
    @staticmethod
    def available_equations() -> List[str]:
        return sorted(load_references())

    @property
    def parameters(self) -> List[str]:
        """Parameters available in the selected reference set."""
        return self._ref.parameters

    @property
    def metadata(self) -> dict:
        """Study metadata (population, country, ranges, device, DOI ...)."""
        return dict(self._ref.metadata)

    def __repr__(self) -> str:
        return (f"OscillometryReference(equation={self.equation!r}, age_unit={self.age_unit!r}, "
                f"height_unit={self.height_unit!r}, weight_unit={self.weight_unit!r}, "
                f"result_unit={self.result_unit!r})")

    # ---------------------------------------------------------------- public
    def calculate(self, age, height, weight, sex, param, value=None) -> OscillometryResult:
        """Evaluate a single subject.

        ``value`` is the measured value (in ``result_unit``). Without it only the
        predicted value and limits of normal are returned; ``percent_predicted``
        and ``z_score`` are ``None`` (and this is not flagged).
        """
        a = np.asarray([np.nan if age is None else age], float)
        h = np.asarray([np.nan if height is None else height], float)
        w = np.asarray([np.nan if weight is None else weight], float)
        s = normalise_sex_array([sex])
        v = None if value is None else np.asarray([value], float)
        for x in (age, height, weight, value):
            if x is not None and np.size(x) != 1:
                raise ValueError("calculate() takes scalar inputs; use compute() for many rows.")
        out, flags = self._evaluate(param, s, a, h, w, v, observed_given=value is not None)
        model = self._ref.find_model(param)
        unit = self._display_unit(model) if model else ""
        res = {k: (None if math.isnan(out[k][0]) else float(out[k][0])) for k in _METRICS}
        return OscillometryResult(flag=flags[0], param=model.param_key if model else str(param),
                                  reference=self.equation, unit=unit, **res)

    def compute(self, df: pd.DataFrame, param: str, sex_col: str = "sex", age_col: str = "age",
                height_col: str = "height", weight_col: str = "weight",
                param_col: Optional[str] = None, prefix: str = "") -> pd.DataFrame:
        """Evaluate every row of ``df``.

        Returns a DataFrame (same index as ``df``) with columns ``predicted``,
        ``percent_predicted``, ``z_score``, ``lln``, ``uln`` and ``flag`` (optionally
        prefixed). ``param_col`` names the column of measured values; if omitted and
        ``df`` has a column called ``param`` that column is used; if there is none,
        only predicted values and limits are returned.
        """
        for col in (sex_col, age_col, height_col, weight_col):
            if col not in df.columns:
                raise KeyError(f"Column {col!r} not found in DataFrame (columns: {list(df.columns)}).")
        if param_col is None and param in df.columns:
            param_col = param
        if param_col is not None and param_col not in df.columns:
            raise KeyError(f"param_col {param_col!r} not found in DataFrame.")
        num = lambda c: pd.to_numeric(df[c], errors="coerce").to_numpy(dtype=float)
        sex = normalise_sex_array(df[sex_col].tolist())
        value = None if param_col is None else num(param_col)
        out, flags = self._evaluate(param, sex, num(age_col), num(height_col), num(weight_col),
                                    value, observed_given=param_col is not None)
        res = pd.DataFrame({f"{prefix}{k}": out[k] for k in _METRICS}, index=df.index)
        res[f"{prefix}flag"] = pd.Series(flags, index=df.index, dtype=object)
        return res

    # --------------------------------------------------------------- internal
    def _scale(self, model: Model) -> float:
        """Multiplier converting reference-equation units to the requested result unit."""
        f = model.pressure_factor
        return 1.0 if f is None else f / self._result_factor

    def _display_unit(self, model: Model) -> str:
        name = reference_pressure_name(model.unit)
        return model.unit if name is None else self.result_unit + model.unit[len(name):]

    def _evaluate(self, param, sex, age, height, weight, value, observed_given):
        n = age.shape[0]
        out: Dict[str, np.ndarray] = {k: np.full(n, np.nan) for k in _METRICS}
        fb = _Flags(n)
        model = self._ref.find_model(param)
        if model is None:
            fb.add_all(f"parameter {str(param)!r} is not available in reference set "
                       f"{self.equation!r} (available: {', '.join(self._ref.parameters)})")
            return out, fb.result()

        age_y = age * self._age_factor
        h_cm = height * self._height_factor
        w_kg = weight * self._weight_factor
        with np.errstate(invalid="ignore", divide="ignore"):
            bmi = w_kg / (h_cm / 100.0) ** 2

        sex_ok = (sex == "Male") | (sex == "Female")
        male = (sex == "Male").astype(float)
        age_ok, h_ok, w_ok = age_y > 0, h_cm > 0, w_kg > 0
        bmi_ok = np.isfinite(bmi)

        fb.add(~sex_ok, "sex missing or not recognised")
        fb.add(~h_ok, "height missing or not positive")
        fb.add(~w_ok, "weight missing or not positive")
        if model.needs_age:
            fb.add(~age_ok, "age missing or not positive")
        else:
            fb.add(~age_ok, "age missing or not positive (not used by this equation, "
                            "but the age range could not be checked)")

        # ---- validity ranges
        amin, amax = np.full(n, np.nan), np.full(n, np.nan)
        for s in ("Male", "Female"):
            rng = model.age_range(s)
            if rng:
                m = sex == s
                amin[m], amax[m] = rng
        _range_flags(fb, "age", "years", np.where(age_ok, age_y, np.nan), amin, amax)
        lim = self._ref.limits
        _range_flags(fb, "height", "cm", np.where(h_ok, h_cm, np.nan), *lim["height"])
        _range_flags(fb, "weight", "kg", np.where(w_ok, w_kg, np.nan), *lim["weight"])
        _range_flags(fb, "BMI", "kg/m2", np.where(bmi_ok & w_ok & h_ok, bmi, np.nan), *lim["bmi"])

        ok = sex_ok & h_ok & (w_ok | (not model.needs_weight)) & (age_ok | (not model.needs_age))
        scale = self._scale(model)
        obs = None if value is None else value / scale   # measured value in reference units

        if not model.has_sd:
            fb.add_all(f"z-score not available: {self.equation} equations are quantile-based "
                       f"and provide no residual standard deviation")

        sexes = ["Both"] if "Both" in model.rows else ["Male", "Female"]
        for s in sexes:
            m = ok if s == "Both" else ok & (sex == s)
            if not m.any():
                continue
            a, h, w, b, ml = age_y[m], h_cm[m], w_kg[m], bmi[m], male[m]
            o = None if obs is None else obs[m]
            rows = model.rows[s]
            if model.method == "quantile":
                pred, lln, uln, z = self._quantile(model, rows, a, h, w, b, ml)
            else:
                pred, lln, uln, z = self._parametric(model, rows, a, h, w, b, ml, o)
            out["predicted"][m] = pred * scale
            out["lln"][m] = lln * scale
            out["uln"][m] = uln * scale
            if z is not None:
                out["z_score"][m] = z
            if value is not None:
                with np.errstate(invalid="ignore", divide="ignore"):
                    pp = 100.0 * value[m] / (pred * scale)
                out["percent_predicted"][m] = np.where(np.isfinite(pp), pp, np.nan)

        if observed_given:
            miss = ~np.isfinite(value)
            fb.add(miss & ok, "measured value missing: percent predicted and z-score not calculated")
            if model.dep_transform != "none" and model.has_sd:
                with np.errstate(invalid="ignore", divide="ignore"):
                    t = _fwd(obs, model.dep_transform, model.dep_offset)
                bad = np.isfinite(value) & ~np.isfinite(t) & ok
                fb.add(bad, "measured value outside the domain of the log-transformed "
                            "equation: z-score not calculated")
        return out, fb.result()

    @staticmethod
    def _parametric(model: Model, rows: List[Row], a, h, w, b, ml, obs):
        if len(rows) == 1:
            lin = _linear_predictor(rows[0], a, h, w, b, ml)
            sd = np.full(a.shape, rows[0].sd)
        else:                                  # segmented: segment k covers ages up to its age_max
            bounds = np.array([r.age_max for r in rows[:-1]])
            seg = np.searchsorted(bounds, a, side="left")
            lin, sd = np.empty(a.shape), np.empty(a.shape)
            for k, r in enumerate(rows):
                m = seg == k
                if m.any():
                    lin[m] = _linear_predictor(r, a[m], h[m], w[m], b[m], ml[m])
                    sd[m] = r.sd
        s = -1.0 if model.dep_transform == "ln_neg" else 1.0
        tr, off = model.dep_transform, model.dep_offset
        pred = _inv(lin, tr, off)
        lln = _inv(lin - s * Z95 * sd, tr, off)
        uln = _inv(lin + s * Z95 * sd, tr, off)
        z = None
        if obs is not None:
            with np.errstate(invalid="ignore", divide="ignore"):
                z = s * (_fwd(obs, tr, off) - lin) / sd
            z = np.where(np.isfinite(z), z, np.nan)
        return pred, lln, uln, z

    @staticmethod
    def _quantile(model: Model, rows: Dict[int, Row], a, h, w, b, ml):
        tr, off = model.dep_transform, model.dep_offset
        pred = _inv(_linear_predictor(rows[50], a, h, w, b, ml), tr, off)
        lln = _inv(_linear_predictor(rows[5], a, h, w, b, ml), tr, off)
        uln = _inv(_linear_predictor(rows[95], a, h, w, b, ml), tr, off)
        return pred, lln, uln, None
