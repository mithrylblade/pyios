"""Load the bundled reference-equation tables into ready-to-evaluate objects."""
from __future__ import annotations

import math
from dataclasses import dataclass, field
from functools import lru_cache
from importlib import resources
from typing import Dict, List, Optional, Tuple

import pandas as pd

from ._names import normalise_param
from ._units import reference_pressure_factor

_COEF_COLS = (
    "sex_coef", "age_coef", "height_coef", "weight_coef", "bmi_coef",
    "inv_height_coef", "inv_bmi_coef",
)


def _opt(value) -> Optional[float]:
    """NaN/NA or an explicit 0 both mean 'term not used'."""
    if value is None or pd.isna(value) or float(value) == 0.0:
        return None
    return float(value)


def _num(value) -> float:
    return math.nan if value is None or pd.isna(value) else float(value)


@dataclass
class Row:
    """One published equation (one sex / percentile / age segment)."""
    intercept: float
    sex_coef: Optional[float]
    age_coef: Optional[float]
    height_coef: Optional[float]
    weight_coef: Optional[float]
    bmi_coef: Optional[float]
    inv_height_coef: Optional[float]
    inv_bmi_coef: Optional[float]
    sd: float                      # rmse or rsd; NaN when not published
    age_transform: str
    height_unit: str               # unit the height coefficients were fitted in ('m' or 'cm')
    age_min: float
    age_max: float
    segment: Optional[int]
    percentile: Optional[int]


@dataclass
class Model:
    """All equations for one parameter of one reference set."""
    reference_id: str
    param_key: str
    parameter: str
    frequency_hz: Optional[float]
    unit: str
    method: str                    # linear | log-linear | quantile | segmented
    dep_transform: str             # none | ln | ln_neg
    dep_offset: float
    rows: Dict[str, object] = field(default_factory=dict)   # sex -> [Row] or {percentile: Row}
    ranges: Dict[str, Tuple[float, float]] = field(default_factory=dict)  # sex -> (age_min, age_max)
    needs_age: bool = False
    needs_weight: bool = False

    @property
    def pressure_factor(self):
        return reference_pressure_factor(self.unit)

    @property
    def has_sd(self) -> bool:
        return self.method != "quantile"

    def age_range(self, sex):
        return self.ranges.get(sex) or self.ranges.get("Both")


@dataclass
class ReferenceSet:
    reference_id: str
    metadata: dict
    models: Dict[str, Model]
    limits: Dict[str, Tuple[float, float]]    # height (cm), weight (kg), bmi (kg/m2)
    _upper: Dict[str, str] = field(default_factory=dict, repr=False)

    @property
    def parameters(self) -> List[str]:
        return sorted(self.models)

    def find_model(self, param) -> Optional[Model]:
        key = normalise_param(param)
        hit = self._upper.get(key)
        if hit is None and key == "AX":           # default AX = integrated from 5 Hz
            hit = self._upper.get("AX5")
        return self.models[hit] if hit else None


def _read_csv(name: str) -> pd.DataFrame:
    resource = resources.files("pyios") / "data" / name
    with resources.as_file(resource) as path:
        return pd.read_csv(path)


def _param_key(parameter: str, freq) -> str:
    if parameter == "AX" and not pd.isna(freq):
        return f"AX{int(freq)}"
    return parameter


def _make_row(rec: dict) -> Row:
    sd = _num(rec.get("rmse"))
    if math.isnan(sd):
        sd = _num(rec.get("rsd"))
    age_tf = rec.get("age_transform")
    return Row(
        intercept=float(rec["intercept"]),
        sd=sd,
        age_transform="none" if pd.isna(age_tf) else str(age_tf),
        height_unit=str(rec["height_unit"]),
        age_min=float(rec["age_min"]),
        age_max=float(rec["age_max"]),
        segment=None if pd.isna(rec["segment"]) else int(rec["segment"]),
        percentile=None if pd.isna(rec["percentile"]) else int(rec["percentile"]),
        **{c: _opt(rec.get(c)) for c in _COEF_COLS},
    )


def _limits(meta: dict) -> Dict[str, Tuple[float, float]]:
    hf = 100.0 if str(meta["input_height_unit"]).lower() == "m" else 1.0
    def pair(lo, hi, f=1.0):
        return (_num(meta[lo]) * f, _num(meta[hi]) * f)
    return {
        "height": pair("height_min", "height_max", hf),
        "weight": pair("weight_min", "weight_max"),
        "bmi": pair("bmi_min", "bmi_max"),
    }


@lru_cache(maxsize=1)
def load_references() -> Dict[str, ReferenceSet]:
    eq = _read_csv("reference_equations.csv")
    meta = _read_csv("reference_metadata.csv")
    if "output_unit" in eq.columns and "parameter_unit" not in eq.columns:
        eq = eq.rename(columns={"output_unit": "parameter_unit"})
    eq["paper_id"] = eq["reference_id"].str.split("_").str[0]
    eq["param_key"] = [_param_key(p, f) for p, f in zip(eq["parameter"], eq["frequency_hz"])]

    refs: Dict[str, ReferenceSet] = {}
    for rec in meta.to_dict("records"):
        pid = rec["reference_id"]
        models: Dict[str, Model] = {}
        sub = eq[eq["paper_id"] == pid]
        for key, grp in sub.groupby("param_key", sort=False):
            first = grp.iloc[0]
            freq = first["frequency_hz"]
            model = Model(
                reference_id=pid, param_key=key, parameter=str(first["parameter"]),
                frequency_hz=None if pd.isna(freq) else float(freq),
                unit=str(first["parameter_unit"]), method=str(first["equation_method"]),
                dep_transform=str(first["dependent_transform"]),
                dep_offset=float(first["dependent_offset"]),
            )
            for sex, sgrp in grp.groupby("sex", sort=False):
                rows = [(_make_row(r), r) for r in sgrp.to_dict("records")]
                if model.method == "quantile":
                    model.rows[sex] = {row.percentile: row for row, _ in rows}
                else:
                    ordered = sorted(rows, key=lambda t: t[0].segment or 0)
                    model.rows[sex] = [row for row, _ in ordered]
                allrows = [row for row, _ in rows]
                model.ranges[sex] = (min(r.age_min for r in allrows), max(r.age_max for r in allrows))
            flat = [r for v in model.rows.values()
                    for r in (v.values() if isinstance(v, dict) else v)]
            model.needs_age = model.method == "segmented" or any(r.age_coef for r in flat)
            model.needs_weight = any(r.weight_coef or r.bmi_coef or r.inv_bmi_coef for r in flat)
            models[key] = model
        refs[pid] = ReferenceSet(
            reference_id=pid, metadata=rec, models=models, limits=_limits(rec),
            _upper={k.upper(): k for k in models},
        )
    return refs
