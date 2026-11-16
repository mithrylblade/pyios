import math
import unittest

import numpy as np
import pandas as pd

from pyios import OscillometryReference as Ref

Z = 1.6449


def close(a, b, rel=1e-9, abs_=1e-12):
    return math.isclose(a, b, rel_tol=rel, abs_tol=abs_)


class TestHandComputedValues(unittest.TestCase):
    """Expected numbers are computed from the published coefficients by hand."""

    def test_moitra_male_r5(self):
        # lnR5 = -0.30 + 0.003*age + 0.01*weight - 0.83*height(m); RMSE 0.40
        lin = -0.30 + 0.003 * 50 + 0.01 * 80 - 0.83 * 1.75
        r = Ref("moitra").calculate(50, 175, 80, "male", "R5")
        self.assertTrue(close(r.predicted, math.exp(lin)))
        self.assertTrue(close(r.lln, math.exp(lin - Z * 0.40)))
        self.assertTrue(close(r.uln, math.exp(lin + Z * 0.40)))

    def test_moitra_x5_is_untransformed(self):
        pred = -0.23 - 0.002 * 50 - 0.002 * 80 + 0.15 * 1.75
        r = Ref("moitra").calculate(50, 175, 80, "male", "X5")
        self.assertTrue(close(r.predicted, pred))
        self.assertTrue(close(r.lln, pred - Z * 0.23))

    def test_moitra_ln_age_term(self):
        lin = 2.44 + 0.39 * math.log(50) + 0.001 * 80 - 0.55 * 1.75
        r = Ref("moitra").calculate(50, 175, 80, "male", "Fres")
        self.assertTrue(close(r.predicted, math.exp(lin)))

    def test_oostveen_r5_hpa_to_kpa(self):
        # ln(R5) = a + b*height(m) + c*age + d*weight in hPa*s/L ; RSD 0.2803
        lin = 5.327 - 3.032 * 1.75 - 0.00381 * 50 + 0.01390 * 80
        r = Ref("oostveen").calculate(50, 175, 80, "male", "R5")
        self.assertTrue(close(r.predicted, math.exp(lin) * 0.1))
        self.assertTrue(close(r.uln, math.exp(lin + Z * 0.2803) * 0.1))
        r2 = Ref("oostveen", result_unit="hPa").calculate(50, 175, 80, "male", "R5")
        self.assertTrue(close(r2.predicted, math.exp(lin)))

    def test_oostveen_x5_reflected_log(self):
        # ln(4 - X5) = a + b*h + c*age + d*w ; LLN of X is MORE negative than predicted
        lin = 2.683 - 0.703 * 1.75 + 0 * 50 + 0.00190 * 80
        sd = 0.0728
        r = Ref("oostveen").calculate(50, 175, 80, "male", "X5")
        self.assertTrue(close(r.predicted, (4 - math.exp(lin)) * 0.1))
        self.assertTrue(close(r.lln, (4 - math.exp(lin + Z * sd)) * 0.1))
        self.assertTrue(close(r.uln, (4 - math.exp(lin - Z * sd)) * 0.1))
        self.assertLess(r.lln, r.predicted)
        self.assertLess(r.predicted, r.uln)

    def test_berger_offsets_and_bmi(self):
        bmi = 80 / 1.75 ** 2
        # ln(R5-R20 + 1) = a + b*h + c*BMI  (cmH2O*s/L)
        lin = 0.37755 - 0.70125 * 1.75 + 0.03588 * bmi
        r = Ref("berger").calculate(50, 175, 80, "male", "R5-R20")
        self.assertTrue(close(r.predicted, (math.exp(lin) - 1) * 0.0980665))
        # ln(X5 + 4) = ...
        lin = 0.10509 + 0.82626 * 1.75 - 0.01775 * bmi
        r = Ref("berger").calculate(50, 175, 80, "male", "X5")
        self.assertTrue(close(r.predicted, (math.exp(lin) - 4) * 0.0980665))
        self.assertTrue(close(r.lln, (math.exp(lin - Z * 0.241371) - 4) * 0.0980665))

    def test_berger_ignores_age(self):
        a = Ref("berger").calculate(30, 175, 80, "male", "R5").predicted
        b = Ref("berger").calculate(70, 175, 80, "male", "R5").predicted
        self.assertEqual(a, b)

    def test_schulz_quantile_equations(self):
        p50 = 1.0719341 - 0.0001490 * 60 - 0.0059151 * 172 + 0.0029873 * 80
        p5 = 1.2018022 - 0.0025997 * 60 - 0.0048804 * 172        # no weight term
        p95 = 0.6514286 + 0.0035450 * 60 - 0.0026455 * 172
        r = Ref("schulz").calculate(60, 172, 80, "male", "Z5")
        self.assertTrue(close(r.predicted, p50))
        self.assertTrue(close(r.lln, p5))
        self.assertTrue(close(r.uln, p95))
        self.assertIsNone(r.z_score)
        self.assertIn("z-score not available", r.flag)

    def test_schulz_negative_intercepts(self):
        # men R5-R20 ratio 95th percentile intercept is -2.6268375 (sign checked on page image)
        p95 = -2.6268375 + 0.4546831 * 60 + 0.1793720 * 172
        r = Ref("schulz").calculate(60, 172, 80, "male", "(R5-R20)/R20")
        self.assertTrue(close(r.uln, p95))

    def test_gochicoa_segments_and_sex(self):
        def lin(inter, sexc, agec, hc, bc, rmse, age, sex, h=175.0, w=80.0):
            bmi = w / (h / 100) ** 2
            return inter + sexc * sex + agec * age + hc / h + bc / bmi, rmse
        # R5 third segment (19.11-90): age 50
        l, sd = lin(-0.28559, -0.015418, 0.0019316, 100.379581, -2.606237, 0.0724, 50, 1)
        r = Ref("gochicoa2023").calculate(50, 175, 80, "male", "R5")
        self.assertTrue(close(r.predicted, l))
        self.assertTrue(close(r.lln, l - Z * sd))
        # female: sex coefficient switched off
        l, sd = lin(-0.28559, -0.015418, 0.0019316, 100.379581, -2.606237, 0.0724, 50, 0)
        r = Ref("gochicoa2023").calculate(50, 175, 80, "female", "R5")
        self.assertTrue(close(r.predicted, l))
        # first segment (2.70-7.36), a child of 5
        l, sd = lin(0.446389, -0.015418, -0.054965, 100.379581, -2.606237, 0.143, 5, 1, 110, 19)
        r = Ref("gochicoa2023").calculate(5, 110, 19, "male", "R5")
        self.assertTrue(close(r.predicted, l))

    def test_gochicoa_segment_boundaries_and_gap(self):
        ref = Ref("gochicoa2023")
        seg1 = lambda a: -0.054965 * a + 0.446389
        seg2 = lambda a: -0.024709 * a + 0.2235
        def pred(a):
            return ref.calculate(a, 110, 19, "male", "R5").predicted
        # same covariates in both segments except intercept/age term
        base = -0.015418 + 100.379581 / 110 - 2.606237 / (19 / 1.1 ** 2)
        self.assertTrue(close(pred(7.36), seg1(7.36) + base))      # inclusive upper bound, seg 1
        self.assertTrue(close(pred(7.365), seg2(7.365) + base))    # in the 7.36-7.37 gap -> seg 2
        self.assertTrue(close(pred(7.37), seg2(7.37) + base))

    def test_gochicoa_ax_is_log(self):
        r = Ref("gochicoa2023").calculate(50, 175, 80, "male", "AX")
        self.assertGreater(r.predicted, 0)
        self.assertLess(r.lln, r.predicted)


class TestUnitsAndNames(unittest.TestCase):
    def test_equation_name_variants(self):
        for name in ["gochicoa", "gochicoa2023", "gochicoa_2023", "gochicoa 2023",
                     "GOCHICOA", "Gochicoa-2023", "GOCHICOA_2023"]:
            self.assertEqual(Ref(name).equation, "gochicoa2023")
        self.assertEqual(Ref("Shulz 2013").equation, "schulz2013")
        with self.assertRaises(ValueError):
            Ref("smith2020")
        with self.assertRaises(ValueError):
            Ref("gochicoa2019")

    def test_sex_variants(self):
        ref = Ref("gochicoa")
        female = [ref.calculate(50, 175, 80, s, "R5").predicted for s in
                  (0, "0", "f", "F", "female", "Female", "FEMALE", 0.0)]
        male = [ref.calculate(50, 175, 80, s, "R5").predicted for s in
                (1, "1", "m", "M", "male", "Male", "MALE", 1.0)]
        self.assertEqual(len(set(female)), 1)
        self.assertEqual(len(set(male)), 1)
        self.assertNotEqual(female[0], male[0])

    def test_invalid_sex_is_flagged(self):
        for bad in (None, "x", 2, float("nan"), ""):
            r = Ref("gochicoa").calculate(50, 175, 80, bad, "R5")
            self.assertIsNone(r.predicted)
            self.assertIn("sex", r.flag)

    def test_input_units(self):
        base = Ref("oostveen").calculate(50, 175, 80, "male", "R5")
        for kw in ({"height_unit": "m"}, {"height_unit": "in"}, {"weight_unit": "lb"},
                   {"age_unit": "months"}, {"age_unit": "days"}):
            h, w, a = 175, 80, 50
            if kw.get("height_unit") == "m": h = 1.75
            if kw.get("height_unit") == "in": h = 175 / 2.54
            if kw.get("weight_unit") == "lb": w = 80 / 0.45359237
            if kw.get("age_unit") == "months": a = 600
            if kw.get("age_unit") == "days": a = 50 * 365.25
            r = Ref("oostveen", **kw).calculate(a, h, w, "male", "R5")
            self.assertTrue(close(r.predicted, base.predicted, rel=1e-9), kw)

    def test_result_units(self):
        for eq in ("gochicoa", "oostveen", "berger", "moitra", "schulz"):
            k = Ref(eq, result_unit="kPa").calculate(60, 170, 70, "female", "R5")
            h = Ref(eq, result_unit="hPa").calculate(60, 170, 70, "female", "R5")
            self.assertTrue(close(h.predicted, k.predicted * 10), eq)
            self.assertTrue(close(h.lln, k.lln * 10), eq)
            self.assertTrue(close(h.uln, k.uln * 10), eq)
            self.assertEqual(h.unit, "hPa*s/L")
            self.assertEqual(k.unit, "kPa*s/L")

    def test_measured_value_unit_follows_result_unit(self):
        k = Ref("gochicoa").calculate(50, 175, 80, "male", "R5", value=0.30)
        h = Ref("gochicoa", result_unit="hPa").calculate(50, 175, 80, "male", "R5", value=3.0)
        self.assertTrue(close(k.percent_predicted, h.percent_predicted))
        self.assertTrue(close(k.z_score, h.z_score))

    def test_non_pressure_parameters_not_scaled(self):
        a = Ref("gochicoa", result_unit="kPa").calculate(50, 175, 80, "male", "Fres")
        b = Ref("gochicoa", result_unit="hPa").calculate(50, 175, 80, "male", "Fres")
        self.assertEqual(a.predicted, b.predicted)
        self.assertEqual(a.unit, "Hz")

    def test_bad_units_rejected(self):
        with self.assertRaises(ValueError):
            Ref("gochicoa", height_unit="furlong")
        with self.assertRaises(ValueError):
            Ref("gochicoa", result_unit="mmHg")

    def test_parameter_name_variants(self):
        ref = Ref("berger")
        a = ref.calculate(50, 175, 80, "male", "R5-R20").predicted
        for p in ("r5-r20", "R5_R20", "R5\u2013R20", "R5-20", " r5-20 "):
            self.assertEqual(ref.calculate(50, 175, 80, "male", p).predicted, a)
        self.assertIsNotNone(Ref("gochicoa").calculate(50, 175, 80, "male", "fres").predicted)

    def test_ax_default_and_explicit(self):
        r = Ref("oostveen")
        self.assertEqual(r.calculate(50, 175, 80, "male", "AX").predicted,
                         r.calculate(50, 175, 80, "male", "AX5").predicted)
        self.assertNotEqual(r.calculate(50, 175, 80, "male", "AX4").predicted,
                            r.calculate(50, 175, 80, "male", "AX5").predicted)
        self.assertIsNotNone(Ref("berger").calculate(50, 175, 80, "male", "AX").predicted)


class TestFlags(unittest.TestCase):
    def test_valid_has_null_flag(self):
        for eq in ("moitra", "oostveen", "berger", "gochicoa"):
            r = Ref(eq).calculate(50, 170, 70, "male", "R5", value=0.3)
            self.assertIsNone(r.flag, (eq, r.flag))

    def test_parameter_not_in_reference(self):
        r = Ref("schulz").calculate(60, 170, 70, "male", "R10")
        self.assertIsNone(r.predicted)
        self.assertIsNone(r.lln)
        self.assertIn("not available", r.flag)
        r = Ref("berger").calculate(50, 170, 70, "male", "Z5")
        self.assertIsNone(r.predicted)
        self.assertIn("Z5", r.flag)

    def test_out_of_range_values_flagged_but_computed(self):
        r = Ref("schulz").calculate(30, 170, 70, "male", "R5")
        self.assertIsNotNone(r.predicted)
        self.assertIn("age 30 years below reference minimum 45", r.flag)
        r = Ref("moitra").calculate(50, 170, 100, "male", "R5")        # BMI 34.6
        self.assertIn("BMI", r.flag)
        self.assertIn("above reference maximum 30", r.flag)
        r = Ref("berger").calculate(50, 220, 70, "male", "R5")
        self.assertIn("height 220 cm above reference maximum", r.flag)
        r = Ref("gochicoa").calculate(1.5, 80, 10, "male", "R5")
        self.assertIn("age 1.5 years below", r.flag)
        self.assertIn("height 80 cm below", r.flag)

    def test_sex_specific_age_limits(self):
        # Oostveen: men 18-84, women 19-81
        self.assertIsNone(Ref("oostveen").calculate(83, 175, 80, "male", "R5").flag)
        self.assertIn("age 83", Ref("oostveen").calculate(83, 165, 70, "female", "R5").flag)

    def test_missing_inputs(self):
        r = Ref("gochicoa").calculate(None, 175, 80, "male", "R5")
        self.assertIsNone(r.predicted)
        self.assertIn("age", r.flag)
        r = Ref("gochicoa").calculate(50, 175, None, "male", "R5")
        self.assertIn("weight", r.flag)
        r = Ref("gochicoa").calculate(50, 0, 80, "male", "R5")
        self.assertIn("height", r.flag)
        # Berger does not use age: still computed, but the flag says range could not be checked
        r = Ref("berger").calculate(None, 175, 80, "male", "R5")
        self.assertIsNotNone(r.predicted)
        self.assertIn("not used by this equation", r.flag)

    def test_multiple_messages_joined(self):
        r = Ref("schulz").calculate(30, 220, 70, "male", "R10")
        self.assertIn("not available", r.flag)
        r = Ref("schulz").calculate(30, 220, 70, "male", "R5")
        self.assertGreaterEqual(r.flag.count(";"), 2)

    def test_value_not_given_not_flagged(self):
        r = Ref("gochicoa").calculate(50, 175, 80, "male", "R5")
        self.assertIsNone(r.flag)
        self.assertIsNone(r.percent_predicted)
        self.assertIsNone(r.z_score)

    def test_measured_value_outside_log_domain(self):
        r = Ref("moitra").calculate(50, 175, 80, "male", "R5", value=-0.1)
        self.assertIsNone(r.z_score)
        self.assertIn("log-transformed", r.flag)


class TestStatistics(unittest.TestCase):
    SUBJECT = dict(age=60, height=170, weight=70)

    def test_limits_bracket_prediction_everywhere(self):
        for eq in Ref.available_equations():
            ref = Ref(eq)
            for param in ref.parameters:
                for sex in ("male", "female"):
                    r = ref.calculate(sex=sex, param=param, **self.SUBJECT)
                    self.assertIsNotNone(r.predicted, (eq, param, sex))
                    self.assertLess(r.lln, r.predicted, (eq, param, sex))
                    self.assertLess(r.predicted, r.uln, (eq, param, sex))

    def test_z_and_percent_roundtrip(self):
        for eq in ("moitra", "oostveen", "berger", "gochicoa"):
            ref = Ref(eq)
            for param in ref.parameters:
                p = ref.calculate(sex="male", param=param, **self.SUBJECT)
                at_pred = ref.calculate(sex="male", param=param, value=p.predicted, **self.SUBJECT)
                at_uln = ref.calculate(sex="male", param=param, value=p.uln, **self.SUBJECT)
                at_lln = ref.calculate(sex="male", param=param, value=p.lln, **self.SUBJECT)
                self.assertAlmostEqual(at_pred.z_score, 0.0, places=7, msg=(eq, param))
                self.assertAlmostEqual(at_uln.z_score, Z, places=6, msg=(eq, param))
                self.assertAlmostEqual(at_lln.z_score, -Z, places=6, msg=(eq, param))
                self.assertAlmostEqual(at_pred.percent_predicted, 100.0, places=7, msg=(eq, param))

    def test_schulz_has_no_z_but_percent(self):
        ref = Ref("schulz")
        p = ref.calculate(sex="female", param="R5", **self.SUBJECT)
        r = ref.calculate(sex="female", param="R5", value=p.predicted * 1.2, **self.SUBJECT)
        self.assertIsNone(r.z_score)
        self.assertAlmostEqual(r.percent_predicted, 120.0, places=7)

    def test_reactance_direction(self):
        # a more negative X5 than predicted must give a negative z-score
        ref = Ref("oostveen")
        p = ref.calculate(sex="male", param="X5", **self.SUBJECT)
        r = ref.calculate(sex="male", param="X5", value=p.predicted * 1.5, **self.SUBJECT)
        self.assertLess(r.z_score, 0)
        self.assertGreater(r.percent_predicted, 100)


class TestCompute(unittest.TestCase):
    def make_df(self):
        return pd.DataFrame({
            "sex": ["male", "F", 1, 0, "unknown", "male"],
            "age": [50, 62, 40, 33, 50, 50],
            "height": [175, 160, 180, 165, 170, 175],
            "weight": [80, 60, 85, 55, 70, np.nan],
            "R5": [0.30, 0.35, 0.28, 0.40, 0.3, 0.3],
        })

    def test_matches_calculate_row_by_row(self):
        df = self.make_df()
        for eq in Ref.available_equations():
            ref = Ref(eq)
            out = ref.compute(df, param="R5", param_col="R5")
            self.assertEqual(list(out.columns),
                             ["predicted", "percent_predicted", "z_score", "lln", "uln", "flag"])
            self.assertTrue(out.index.equals(df.index))
            for i, row in df.iterrows():
                r = ref.calculate(row.age, row.height, row.weight, row.sex, "R5", value=row.R5)
                for k in ("predicted", "percent_predicted", "z_score", "lln", "uln"):
                    got = out.loc[i, k]
                    want = getattr(r, k)
                    if want is None:
                        self.assertTrue(pd.isna(got), (eq, i, k))
                    else:
                        self.assertTrue(close(got, want), (eq, i, k, got, want))
                want_flag = r.flag
                got_flag = out.loc[i, "flag"]
                self.assertEqual(None if pd.isna(got_flag) else got_flag, want_flag, (eq, i))

    def test_defaults_and_param_col_default(self):
        df = self.make_df()
        ref = Ref("gochicoa")
        a = ref.compute(df, param="R5")                       # column named R5 used automatically
        b = ref.compute(df, param="R5", sex_col="sex", age_col="age", height_col="height",
                        weight_col="weight", param_col="R5")
        pd.testing.assert_frame_equal(a, b)

    def test_without_measured_values(self):
        df = self.make_df().drop(columns="R5")
        out = Ref("gochicoa").compute(df, param="R5")
        self.assertTrue(out["percent_predicted"].isna().all())
        self.assertTrue(out["z_score"].isna().all())
        self.assertTrue(out["predicted"].iloc[[0, 1, 2, 3]].notna().all())

    def test_custom_columns_units_and_prefix(self):
        df = pd.DataFrame({"Sex": ["M"], "Age_m": [600], "ht_m": [1.75], "wt_lb": [176.37], "r5": [3.0]})
        ref = Ref("gochicoa", age_unit="months", height_unit="m", weight_unit="lb", result_unit="hPa")
        out = ref.compute(df, param="R5", sex_col="Sex", age_col="Age_m", height_col="ht_m",
                          weight_col="wt_lb", param_col="r5", prefix="R5_")
        self.assertIn("R5_predicted", out.columns)
        want = Ref("gochicoa", result_unit="hPa").calculate(50, 175, 80.0, "male", "R5", value=3.0)
        self.assertAlmostEqual(out["R5_predicted"].iloc[0], want.predicted, places=3)

    def test_unknown_param_and_missing_columns(self):
        df = self.make_df()
        out = Ref("schulz").compute(df, param="R10")
        self.assertTrue(out["predicted"].isna().all())
        self.assertTrue(out["flag"].str.contains("not available").all())
        with self.assertRaises(KeyError):
            Ref("gochicoa").compute(df.drop(columns="age"), param="R5")
        with self.assertRaises(KeyError):
            Ref("gochicoa").compute(df, param="R5", param_col="nope")

    def test_nonstandard_index_and_text_numbers(self):
        df = pd.DataFrame({"sex": ["male", "female"], "age": ["50", "abc"], "height": [175, 165],
                           "weight": [80, 60]}, index=[10, 20])
        out = Ref("gochicoa").compute(df, param="R5")
        self.assertEqual(list(out.index), [10, 20])
        self.assertTrue(pd.notna(out.loc[10, "predicted"]))
        self.assertTrue(pd.isna(out.loc[20, "predicted"]))
        self.assertIn("age", out.loc[20, "flag"])

    def test_large_frame(self):
        rng = np.random.default_rng(0)
        n = 200_000
        df = pd.DataFrame({"sex": rng.integers(0, 2, n), "age": rng.uniform(20, 80, n),
                           "height": rng.uniform(150, 190, n), "weight": rng.uniform(50, 100, n)})
        df["R5"] = 0.3
        out = Ref("gochicoa").compute(df, param="R5")
        self.assertEqual(len(out), n)
        self.assertTrue(out["predicted"].notna().all())

    def test_does_not_modify_input(self):
        df = self.make_df()
        before = df.copy(deep=True)
        Ref("berger").compute(df, param="R5")
        pd.testing.assert_frame_equal(df, before)


class TestMetadata(unittest.TestCase):
    def test_metadata_and_parameters(self):
        ref = Ref("berger")
        self.assertEqual(ref.metadata["country"], "USA")
        self.assertIn("R5-R20", ref.parameters)
        self.assertEqual(Ref.available_equations(),
                         ["berger2021", "gochicoa2023", "moitra2020", "oostveen2013", "schulz2013"])


class TestRPackageData(unittest.TestCase):
    def test_r_package_data_in_sync(self):
        import pathlib
        root = pathlib.Path(__file__).resolve().parents[1]
        rdata = root / "r" / "inst" / "extdata"
        if not rdata.exists():
            self.skipTest("R package not present")
        for name in ("reference_equations.csv", "reference_metadata.csv"):
            self.assertEqual((rdata / name).read_bytes(),
                             (root / "src" / "pyios" / "data" / name).read_bytes(), name)


if __name__ == "__main__":
    unittest.main()
