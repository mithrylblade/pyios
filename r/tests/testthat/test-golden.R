# Parity with the Python package: every case in golden.csv was produced by pyios (Python).
golden <- utils::read.csv(test_path("golden.csv"), stringsAsFactors = FALSE, na.strings = "",
                          fileEncoding = "UTF-8")

test_that("R results match the Python implementation (values and flag messages)", {
  expect_gt(nrow(golden), 700)
  for (i in seq_len(nrow(golden))) {
    g <- golden[i, ]
    ref <- OscillometryReference(g$equation, age_unit = g$age_unit, height_unit = g$height_unit,
                                 weight_unit = g$weight_unit, result_unit = g$result_unit)
    value <- if (is.na(g$value)) NULL else g$value
    r <- ref$calculate(age = g$age, height = g$height, weight = g$weight, sex = g$sex,
                       param = g$param, value = value)
    label <- sprintf("case %d: %s %s %s", i, g$equation, g$param, g$sex)
    for (m in c("predicted", "percent_predicted", "z_score", "lln", "uln")) {
      expect_equal(r[[m]], g[[paste0("exp_", m)]], tolerance = 1e-8, label = paste(label, m))
    }
    expect_identical(r$flag, g$exp_flag, label = paste(label, "flag"))
    expect_identical(r$param, g$exp_param, label = paste(label, "param"))
    exp_unit <- if (is.na(g$exp_unit)) "" else g$exp_unit      # blank cell = unknown parameter
    expect_identical(r$unit, exp_unit, label = paste(label, "unit"))
  }
})
