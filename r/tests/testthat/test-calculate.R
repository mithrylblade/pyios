Z <- 1.6449
ref <- function(...) OscillometryReference(...)

test_that("Moitra equations (log-linear, X5 untransformed, ln age)", {
  lin <- -0.30 + 0.003 * 50 + 0.01 * 80 - 0.83 * 1.75
  r <- ref("moitra")$calculate(50, 175, 80, "male", "R5")
  expect_equal(r$predicted, exp(lin), tolerance = 1e-9)
  expect_equal(r$lln, exp(lin - Z * 0.40), tolerance = 1e-9)
  expect_equal(r$uln, exp(lin + Z * 0.40), tolerance = 1e-9)

  pred <- -0.23 - 0.002 * 50 - 0.002 * 80 + 0.15 * 1.75
  r <- ref("moitra")$calculate(50, 175, 80, "male", "X5")
  expect_equal(r$predicted, pred, tolerance = 1e-9)

  lin <- 2.44 + 0.39 * log(50) + 0.001 * 80 - 0.55 * 1.75
  expect_equal(ref("moitra")$calculate(50, 175, 80, "male", "Fres")$predicted, exp(lin),
               tolerance = 1e-9)
})

test_that("Oostveen equations: hPa->kPa and reflected log for reactance", {
  lin <- 5.327 - 3.032 * 1.75 - 0.00381 * 50 + 0.01390 * 80
  r <- ref("oostveen")$calculate(50, 175, 80, "male", "R5")
  expect_equal(r$predicted, exp(lin) * 0.1, tolerance = 1e-9)
  expect_equal(r$uln, exp(lin + Z * 0.2803) * 0.1, tolerance = 1e-9)
  r2 <- ref("oostveen", result_unit = "hPa")$calculate(50, 175, 80, "male", "R5")
  expect_equal(r2$predicted, exp(lin), tolerance = 1e-9)

  lin <- 2.683 - 0.703 * 1.75 + 0 * 50 + 0.00190 * 80
  sd <- 0.0728
  r <- ref("oostveen")$calculate(50, 175, 80, "male", "X5")
  expect_equal(r$predicted, (4 - exp(lin)) * 0.1, tolerance = 1e-9)
  expect_equal(r$lln, (4 - exp(lin + Z * sd)) * 0.1, tolerance = 1e-9)
  expect_equal(r$uln, (4 - exp(lin - Z * sd)) * 0.1, tolerance = 1e-9)
  expect_lt(r$lln, r$predicted)
  expect_lt(r$predicted, r$uln)
})

test_that("Berger equations: offsets, BMI, cmH2O, no age term", {
  bmi <- 80 / 1.75^2
  lin <- 0.37755 - 0.70125 * 1.75 + 0.03588 * bmi
  r <- ref("berger")$calculate(50, 175, 80, "male", "R5-R20")
  expect_equal(r$predicted, (exp(lin) - 1) * 0.0980665, tolerance = 1e-9)
  lin <- 0.10509 + 0.82626 * 1.75 - 0.01775 * bmi
  r <- ref("berger")$calculate(50, 175, 80, "male", "X5")
  expect_equal(r$predicted, (exp(lin) - 4) * 0.0980665, tolerance = 1e-9)
  expect_equal(r$lln, (exp(lin - Z * 0.241371) - 4) * 0.0980665, tolerance = 1e-9)
  expect_equal(ref("berger")$calculate(30, 175, 80, "male", "R5")$predicted,
               ref("berger")$calculate(70, 175, 80, "male", "R5")$predicted)
})

test_that("Schulz quantile equations, no z-score, flagged", {
  p50 <- 1.0719341 - 0.0001490 * 60 - 0.0059151 * 172 + 0.0029873 * 80
  p5  <- 1.2018022 - 0.0025997 * 60 - 0.0048804 * 172
  p95 <- 0.6514286 + 0.0035450 * 60 - 0.0026455 * 172
  r <- ref("schulz")$calculate(60, 172, 80, "male", "Z5")
  expect_equal(r$predicted, p50, tolerance = 1e-9)
  expect_equal(r$lln, p5, tolerance = 1e-9)
  expect_equal(r$uln, p95, tolerance = 1e-9)
  expect_true(is.na(r$z_score))
  expect_match(r$flag, "z-score not available")

  p95 <- -2.6268375 + 0.4546831 * 60 + 0.1793720 * 172
  expect_equal(ref("schulz")$calculate(60, 172, 80, "male", "(R5-R20)/R20")$uln, p95,
               tolerance = 1e-9)
})

test_that("Gochicoa: segments, sex coefficient, gaps between bands", {
  lin <- function(inter, sexc, agec, hc, bc, age, sex, h = 175, w = 80) {
    bmi <- w / (h / 100)^2
    inter + sexc * sex + agec * age + hc / h + bc / bmi
  }
  l <- lin(-0.28559, -0.015418, 0.0019316, 100.379581, -2.606237, 50, 1)
  r <- ref("gochicoa2023")$calculate(50, 175, 80, "male", "R5")
  expect_equal(r$predicted, l, tolerance = 1e-9)
  expect_equal(r$lln, l - Z * 0.0724, tolerance = 1e-9)
  l <- lin(-0.28559, -0.015418, 0.0019316, 100.379581, -2.606237, 50, 0)
  expect_equal(ref("gochicoa2023")$calculate(50, 175, 80, "female", "R5")$predicted, l,
               tolerance = 1e-9)

  g <- ref("gochicoa2023")
  base <- -0.015418 + 100.379581 / 110 - 2.606237 / (19 / 1.1^2)
  pred <- function(a) g$calculate(a, 110, 19, "male", "R5")$predicted
  expect_equal(pred(7.36),  -0.054965 * 7.36  + 0.446389 + base, tolerance = 1e-9)  # segment 1
  expect_equal(pred(7.365), -0.024709 * 7.365 + 0.2235   + base, tolerance = 1e-9)  # gap -> seg 2
  expect_equal(pred(7.37),  -0.024709 * 7.37  + 0.2235   + base, tolerance = 1e-9)
})

test_that("equation names are flexible", {
  for (nm in c("gochicoa", "gochicoa2023", "gochicoa_2023", "gochicoa 2023", "GOCHICOA",
               "Gochicoa-2023", "GOCHICOA_2023")) {
    expect_identical(ref(nm)$equation, "gochicoa2023")
  }
  expect_identical(ref("Shulz 2013")$equation, "schulz2013")
  expect_error(ref("smith2020"), "Unknown equation")
  expect_error(ref("gochicoa2019"), "Unknown equation")
  expect_identical(available_equations(),
                   c("berger2021", "gochicoa2023", "moitra2020", "oostveen2013", "schulz2013"))
})

test_that("sex codes are flexible and invalid codes are flagged", {
  g <- ref("gochicoa")
  p <- function(s) g$calculate(50, 175, 80, s, "R5")$predicted
  female <- vapply(list(0, "0", "f", "F", "female", "Female", "FEMALE", 0.0), p, numeric(1))
  male   <- vapply(list(1, "1", "m", "M", "male", "Male", "MALE", 1.0), p, numeric(1))
  expect_length(unique(female), 1L)
  expect_length(unique(male), 1L)
  expect_false(female[1] == male[1])
  for (bad in list(NULL, "x", 2, NA, "")) {
    r <- g$calculate(50, 175, 80, bad, "R5")
    expect_true(is.na(r$predicted))
    expect_match(r$flag, "sex")
  }
})

test_that("input and result units are converted", {
  base <- ref("oostveen")$calculate(50, 175, 80, "male", "R5")$predicted
  expect_equal(ref("oostveen", height_unit = "m")$calculate(50, 1.75, 80, "male", "R5")$predicted, base)
  expect_equal(ref("oostveen", height_unit = "in")$calculate(50, 175 / 2.54, 80, "male", "R5")$predicted, base)
  expect_equal(ref("oostveen", weight_unit = "lb")$calculate(50, 175, 80 / 0.45359237, "male", "R5")$predicted, base)
  expect_equal(ref("oostveen", age_unit = "months")$calculate(600, 175, 80, "male", "R5")$predicted, base)
  for (eq in available_equations()) {
    k <- ref(eq, result_unit = "kPa")$calculate(60, 170, 70, "female", "R5")
    h <- ref(eq, result_unit = "hPa")$calculate(60, 170, 70, "female", "R5")
    expect_equal(h$predicted, k$predicted * 10)
    expect_equal(h$lln, k$lln * 10)
    expect_equal(h$uln, k$uln * 10)
    expect_identical(h$unit, "hPa*s/L")
    expect_identical(k$unit, "kPa*s/L")
  }
  a <- ref("gochicoa", result_unit = "kPa")$calculate(50, 175, 80, "male", "Fres")
  b <- ref("gochicoa", result_unit = "hPa")$calculate(50, 175, 80, "male", "Fres")
  expect_identical(a$predicted, b$predicted)
  expect_identical(a$unit, "Hz")
  expect_error(ref("gochicoa", height_unit = "furlong"), "Unsupported height unit")
  expect_error(ref("gochicoa", result_unit = "mmHg"), "Unsupported result unit")
})

test_that("measured value is interpreted in the result unit", {
  k <- ref("gochicoa")$calculate(50, 175, 80, "male", "R5", value = 0.30)
  h <- ref("gochicoa", result_unit = "hPa")$calculate(50, 175, 80, "male", "R5", value = 3.0)
  expect_equal(k$percent_predicted, h$percent_predicted)
  expect_equal(k$z_score, h$z_score)
})

test_that("parameter names are flexible; AX defaults to the 5 Hz equation", {
  r <- ref("berger")
  a <- r$calculate(50, 175, 80, "male", "R5-R20")$predicted
  for (p in c("r5-r20", "R5_R20", "R5\u2013R20", "R5-20", " r5-20 ")) {
    expect_identical(r$calculate(50, 175, 80, "male", p)$predicted, a)
  }
  o <- ref("oostveen")
  expect_identical(o$calculate(50, 175, 80, "male", "AX")$predicted,
                   o$calculate(50, 175, 80, "male", "AX5")$predicted)
  expect_false(o$calculate(50, 175, 80, "male", "AX4")$predicted ==
                 o$calculate(50, 175, 80, "male", "AX5")$predicted)
  expect_false(is.na(ref("berger")$calculate(50, 175, 80, "male", "AX")$predicted))
})

test_that("flag is NA when valid and informative otherwise", {
  for (eq in c("moitra", "oostveen", "berger", "gochicoa")) {
    r <- ref(eq)$calculate(50, 170, 70, "male", "R5", value = 0.3)
    expect_true(is.na(r$flag), info = eq)
  }
  r <- ref("schulz")$calculate(60, 170, 70, "male", "R10")
  expect_true(is.na(r$predicted)); expect_true(is.na(r$lln))
  expect_match(r$flag, "not available")
  r <- ref("schulz")$calculate(30, 170, 70, "male", "R5")
  expect_false(is.na(r$predicted))
  expect_match(r$flag, "age 30 years below reference minimum 45", fixed = TRUE)
  expect_match(ref("moitra")$calculate(50, 170, 100, "male", "R5")$flag, "above reference maximum 30")
  expect_match(ref("berger")$calculate(50, 220, 70, "male", "R5")$flag, "height 220 cm above reference maximum")
  r <- ref("gochicoa")$calculate(1.5, 80, 10, "male", "R5")
  expect_match(r$flag, "age 1.5 years below"); expect_match(r$flag, "height 80 cm below")
  expect_true(is.na(ref("oostveen")$calculate(83, 175, 80, "male", "R5")$flag))
  expect_match(ref("oostveen")$calculate(83, 165, 70, "female", "R5")$flag, "age 83")
  r <- ref("gochicoa")$calculate(NULL, 175, 80, "male", "R5")
  expect_true(is.na(r$predicted)); expect_match(r$flag, "age")
  r <- ref("berger")$calculate(NULL, 175, 80, "male", "R5")
  expect_false(is.na(r$predicted)); expect_match(r$flag, "not used by this equation")
  r <- ref("gochicoa")$calculate(50, 175, 80, "male", "R5")
  expect_true(is.na(r$flag)); expect_true(is.na(r$percent_predicted)); expect_true(is.na(r$z_score))
  r <- ref("moitra")$calculate(50, 175, 80, "male", "R5", value = -0.1)
  expect_true(is.na(r$z_score)); expect_match(r$flag, "log-transformed")
  expect_error(ref("gochicoa")$calculate(c(1, 2), 175, 80, "male", "R5"), "scalar")
})

test_that("limits bracket the prediction; z and % round-trip", {
  subj <- list(age = 60, height = 170, weight = 70)
  for (eq in available_equations()) {
    r <- ref(eq)
    for (param in r$parameters) for (sex in c("male", "female")) {
      p <- r$calculate(subj$age, subj$height, subj$weight, sex, param)
      lbl <- paste(eq, param, sex)
      expect_false(is.na(p$predicted), label = lbl)
      expect_lt(p$lln, p$predicted); expect_lt(p$predicted, p$uln)
      if (eq != "schulz2013") {
        at <- function(v) r$calculate(subj$age, subj$height, subj$weight, sex, param, value = v)
        expect_equal(at(p$predicted)$z_score, 0, tolerance = 1e-7, label = lbl)
        expect_equal(at(p$uln)$z_score, Z, tolerance = 1e-6, label = lbl)
        expect_equal(at(p$lln)$z_score, -Z, tolerance = 1e-6, label = lbl)
        expect_equal(at(p$predicted)$percent_predicted, 100, tolerance = 1e-7, label = lbl)
      }
    }
  }
  s <- ref("schulz")
  p <- s$calculate(60, 170, 70, "female", "R5")
  r <- s$calculate(60, 170, 70, "female", "R5", value = p$predicted * 1.2)
  expect_true(is.na(r$z_score)); expect_equal(r$percent_predicted, 120, tolerance = 1e-7)
  o <- ref("oostveen")
  p <- o$calculate(60, 170, 70, "male", "X5")
  r <- o$calculate(60, 170, 70, "male", "X5", value = p$predicted * 1.5)
  expect_lt(r$z_score, 0); expect_gt(r$percent_predicted, 100)
})

test_that("metadata and parameters are exposed", {
  r <- ref("berger")
  expect_identical(r$metadata$country, "USA")
  expect_true("R5-R20" %in% r$parameters)
  expect_output(print(r), "OscillometryReference")
})
