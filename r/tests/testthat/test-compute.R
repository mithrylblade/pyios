make_df <- function() {
  data.frame(
    sex = c("male", "F", 1, 0, "unknown", "male"),
    age = c(50, 62, 40, 33, 50, 50),
    height = c(175, 160, 180, 165, 170, 175),
    weight = c(80, 60, 85, 55, 70, NA),
    R5 = c(0.30, 0.35, 0.28, 0.40, 0.3, 0.3),
    stringsAsFactors = FALSE
  )
}

test_that("compute() matches calculate() row by row for every equation", {
  df <- make_df()
  for (eq in available_equations()) {
    ref <- OscillometryReference(eq)
    out <- ref$compute(df, param = "R5", param_col = "R5")
    expect_identical(names(out), c("predicted", "percent_predicted", "z_score", "lln", "uln", "flag"))
    expect_identical(rownames(out), rownames(df))
    for (i in seq_len(nrow(df))) {
      r <- ref$calculate(df$age[i], df$height[i], df$weight[i], df$sex[i], "R5", value = df$R5[i])
      for (k in c("predicted", "percent_predicted", "z_score", "lln", "uln")) {
        expect_equal(out[[k]][i], r[[k]], tolerance = 1e-12, label = paste(eq, i, k))
      }
      expect_identical(out$flag[i], r$flag, label = paste(eq, i, "flag"))
    }
  }
})

test_that("default columns, param_col default and prefix", {
  df <- make_df()
  ref <- OscillometryReference("gochicoa")
  a <- ref$compute(df, param = "R5")
  b <- ref$compute(df, param = "R5", sex_col = "sex", age_col = "age", height_col = "height",
                   weight_col = "weight", param_col = "R5")
  expect_equal(a, b)
  out <- ref$compute(df, param = "R5", prefix = "R5_")
  expect_true("R5_predicted" %in% names(out))
  expect_true("R5_flag" %in% names(out))
})

test_that("compute() without measured values returns predictions only", {
  df <- make_df()[, c("sex", "age", "height", "weight")]
  out <- OscillometryReference("gochicoa")$compute(df, param = "R5")
  expect_true(all(is.na(out$percent_predicted)))
  expect_true(all(is.na(out$z_score)))
  expect_true(all(!is.na(out$predicted[1:4])))
})

test_that("custom columns and units", {
  df <- data.frame(Sex = "M", Age_m = 600, ht_m = 1.75, wt_lb = 176.37, r5 = 3.0)
  ref <- OscillometryReference("gochicoa", age_unit = "months", height_unit = "m",
                               weight_unit = "lb", result_unit = "hPa")
  out <- ref$compute(df, param = "R5", sex_col = "Sex", age_col = "Age_m", height_col = "ht_m",
                     weight_col = "wt_lb", param_col = "r5", prefix = "R5_")
  want <- OscillometryReference("gochicoa", result_unit = "hPa")$calculate(50, 175, 80.0, "male", "R5", value = 3.0)
  expect_equal(out$R5_predicted, want$predicted, tolerance = 1e-3)
})

test_that("unknown parameter, missing columns, text numbers, row names", {
  df <- make_df()
  out <- OscillometryReference("schulz")$compute(df, param = "R10")
  expect_true(all(is.na(out$predicted)))
  expect_true(all(grepl("not available", out$flag)))
  g <- OscillometryReference("gochicoa")
  expect_error(g$compute(df[, names(df) != "age"], param = "R5"), "not found")
  expect_error(g$compute(df, param = "R5", param_col = "nope"), "not found")

  d2 <- data.frame(sex = c("male", "female"), age = c("50", "abc"), height = c(175, 165),
                   weight = c(80, 60), stringsAsFactors = FALSE)
  rownames(d2) <- c("10", "20")
  out <- g$compute(d2, param = "R5")
  expect_identical(rownames(out), c("10", "20"))
  expect_false(is.na(out["10", "predicted"]))
  expect_true(is.na(out["20", "predicted"]))
  expect_match(out["20", "flag"], "age")
})

test_that("factor and numeric sex columns, and large inputs", {
  df <- data.frame(sex = factor(c("male", "female")), age = c(50, 50), height = c(175, 165),
                   weight = c(80, 60))
  out <- OscillometryReference("gochicoa")$compute(df, param = "R5")
  expect_false(anyNA(out$predicted))
  set.seed(1)
  n <- 100000
  big <- data.frame(sex = sample(0:1, n, TRUE), age = runif(n, 20, 80),
                    height = runif(n, 150, 190), weight = runif(n, 50, 100), R5 = 0.3)
  out <- OscillometryReference("gochicoa")$compute(big, param = "R5")
  expect_equal(nrow(out), n)
  expect_false(anyNA(out$predicted))
})

test_that("compute() does not modify its input", {
  df <- make_df()
  before <- df
  OscillometryReference("berger")$compute(df, param = "R5")
  expect_identical(df, before)
})
