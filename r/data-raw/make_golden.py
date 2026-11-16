#!/usr/bin/env python3
"""Generate tests/testthat/golden.csv: expected outputs of the *Python* pyios package.

The R test-suite (tests/testthat/test-golden.R) checks that the R port reproduces every value
and flag message. Re-run after changing equations or flag wording in either implementation:

    PYTHONPATH=src python3 r/data-raw/make_golden.py
"""
import csv
import math
import pathlib

from pyios import OscillometryReference as Ref

OUT = pathlib.Path(__file__).resolve().parents[1] / "tests" / "testthat" / "golden.csv"
COLS = ["equation", "age_unit", "height_unit", "weight_unit", "result_unit", "param", "sex",
        "age", "height", "weight", "value", "exp_predicted", "exp_percent_predicted",
        "exp_z_score", "exp_lln", "exp_uln", "exp_flag", "exp_param", "exp_unit"]
DEFAULT_UNITS = ("years", "cm", "kg", "kPa")
rows = []


def fmt(x):
    if x is None or (isinstance(x, float) and math.isnan(x)):
        return ""
    return format(float(x), ".17g")


def add(eq, param, sex, age, height, weight, value=None, units=DEFAULT_UNITS):
    au, hu, wu, ru = units
    ref = Ref(eq, age_unit=au, height_unit=hu, weight_unit=wu, result_unit=ru)
    r = ref.calculate(age, height, weight, sex, param, value=value)
    rows.append({
        "equation": eq, "age_unit": au, "height_unit": hu, "weight_unit": wu, "result_unit": ru,
        "param": param, "sex": "" if sex is None else sex, "age": fmt(age), "height": fmt(height),
        "weight": fmt(weight), "value": fmt(value),
        "exp_predicted": fmt(r.predicted), "exp_percent_predicted": fmt(r.percent_predicted),
        "exp_z_score": fmt(r.z_score), "exp_lln": fmt(r.lln), "exp_uln": fmt(r.uln),
        "exp_flag": r.flag or "", "exp_param": r.param, "exp_unit": r.unit,
    })


SUBJECTS = [(60, 170, 70), (10, 140, 35), (88, 160, 60), (30, 220, 150), (50, 175, 80)]

# 1. every equation x parameter x sex x subject (with a measured value on alternate subjects)
for eq in Ref.available_equations():
    ref = Ref(eq)
    for param in ref.parameters:
        for sex in ("male", "female"):
            for k, (a, h, w) in enumerate(SUBJECTS):
                value = None
                if k % 2 == 0:
                    p = ref.calculate(a, h, w, sex, param).predicted
                    value = 1.15 * p if p is not None else 0.3
                add(eq, param, sex, a, h, w, value)

# 2. unit variants
for eq in Ref.available_equations():
    add(eq, "R5", "female", 60, 170, 70, 0.35, units=("years", "cm", "kg", "hPa"))
    add(eq, "R5", "male", 60 * 12, 1.70, 70 / 0.45359237, 3.5, units=("months", "m", "lb", "hPa"))
    add(eq, "X5", "male", 60 * 365.25, 1700, 70000, -0.1, units=("days", "mm", "g", "kPa"))
    add(eq, "Fres", "female", 60, 170 / 2.54, 70, 12, units=("years", "in", "kg", "cmH2O"))
    add(eq, "R5", "male", 60, 170, 70, 3.5, units=("years", "cm", "kg", "cmH2O"))

# 3. parameter / equation / sex spelling variants and defaults
for p in ("r5-20", "R5_R20", "R5\u2013R20", " r5-r20 ", "R5-R20", "fres", "ax", "AX4", "AX5", "x5"):
    for eq in ("berger_2021", "Gochicoa", "OOSTVEEN 2013", "moitra", "Shulz"):
        add(eq, p, "male", 55, 172, 74, None)
for p in ("(R5-R20)/R20", "R5-R20/R20", "(R5-R20)/R5", "r5-20/r5"):
    for eq in ("schulz2013", "gochicoa2023"):
        add(eq, p, "female", 55, 165, 60, 5.0)
for s in ("0", "1", "f", "F", "FEMALE", "Male", "m", "M", "1.0", "0.0", "w", "x", ""):
    add("gochicoa", "R5", s, 50, 175, 80, 0.3)

# 4. missing / invalid inputs and unknown parameters
add("gochicoa", "R5", "male", None, 175, 80, 0.3)
add("gochicoa", "R5", "male", 50, 175, None, 0.3)
add("gochicoa", "R5", "male", 50, 0, 80, 0.3)
add("gochicoa", "R5", None, 50, 175, 80, 0.3)
add("berger", "R5", "male", None, 175, 80, 0.3)
add("berger", "Z5", "male", 50, 175, 80, 0.3)
add("schulz", "R10", "male", 60, 170, 70, 0.3)
add("schulz", "R5", "male", 30, 220, 70, 0.3)
add("moitra", "R5", "male", 50, 175, 80, -0.1)
add("moitra", "R5", "male", 50, 170, 100, 0.5)
add("oostveen", "X5", "female", 83, 165, 70, -0.1)
add("oostveen", "X5", "male", 50, 175, 80, 0.5)

# 5. Gochicoa age-segment boundaries and the gaps between segments
for a in (2.7, 2.69, 7.36, 7.365, 7.37, 19.1, 19.105, 19.11, 23.55, 23.555, 23.56, 90, 91):
    for p in ("R5", "R5-R20", "AX", "X10"):
        add("gochicoa2023", p, "male", a, 130 if a < 19 else 175, 30 if a < 19 else 80, None)

OUT.parent.mkdir(parents=True, exist_ok=True)
with open(OUT, "w", newline="", encoding="utf-8") as f:
    w = csv.DictWriter(f, fieldnames=COLS)
    w.writeheader()
    w.writerows(rows)
print(len(rows), "golden cases written to", OUT)
