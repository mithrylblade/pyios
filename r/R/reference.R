# Predicted values, % predicted, z-scores and limits of normal for oscillometry.

# z-score of the 5th / 95th percentile (as used by Gochicoa et al.)
Z95 <- 1.6449

.metrics <- c("predicted", "percent_predicted", "z_score", "lln", "uln")

# ------------------------------------------------------------------ flag helpers
# Append `text` to the flag of every row where `mask` is TRUE (NA counts as FALSE).
.add_flag <- function(flag, mask, text) {
  idx <- which(mask)
  if (length(idx) == 0L) return(flag)
  cur <- flag[idx]
  flag[idx] <- ifelse(is.na(cur), text, paste(cur, text, sep = "; "))
  flag
}

# 6 significant digits like Python's format(x, ".6g"). formatC() pads to width digits + 1 unless
# width is given, so the padding is switched off explicitly.
.g <- function(x) trimws(formatC(x, digits = 6, format = "g", width = 1))

.range_flags <- function(flag, label, unit, values, lo, hi) {
  lo <- rep_len(lo, length(values))
  hi <- rep_len(hi, length(values))
  low  <- !is.na(values) & !is.na(lo) & values < lo
  high <- !is.na(values) & !is.na(hi) & values > hi
  il <- which(low)
  ih <- which(high)
  if (length(il)) {
    flag <- .add_flag(flag, low, paste0(label, " ", .g(values[il]), " ", unit,
                                        " below reference minimum ", .g(lo[il]), " ", unit))
  }
  if (length(ih)) {
    flag <- .add_flag(flag, high, paste0(label, " ", .g(values[ih]), " ", unit,
                                         " above reference maximum ", .g(hi[ih]), " ", unit))
  }
  flag
}

# ------------------------------------------------------------------------ maths
.fwd <- function(y, tr, off) {
  if (tr == "ln") return(suppressWarnings(log(y + off)))
  if (tr == "ln_neg") return(suppressWarnings(log(off - y)))
  y
}

.inv <- function(t, tr, off) {
  if (tr == "ln") return(exp(t) - off)
  if (tr == "ln_neg") return(off - exp(t))
  t
}

.linear_predictor <- function(row, age, h_cm, w, bmi, male) {
  h <- if (row$height_unit == "m") h_cm / 100 else h_cm
  lin <- rep(row$intercept, length(age))
  if (!is.na(row$sex_coef)) lin <- lin + row$sex_coef * male
  if (!is.na(row$age_coef)) {
    lin <- lin + row$age_coef * (if (row$age_transform == "ln") log(age) else age)
  }
  if (!is.na(row$height_coef)) lin <- lin + row$height_coef * h
  if (!is.na(row$inv_height_coef)) lin <- lin + row$inv_height_coef / h
  if (!is.na(row$weight_coef)) lin <- lin + row$weight_coef * w
  if (!is.na(row$bmi_coef)) lin <- lin + row$bmi_coef * bmi
  if (!is.na(row$inv_bmi_coef)) lin <- lin + row$inv_bmi_coef / bmi
  lin
}

# Linear / log-linear / segmented models: limits are Z95 residual SDs either side in the
# transformed space.
.parametric <- function(model, rows, a, h, w, b, ml, obs) {
  if (length(rows) == 1L) {
    lin <- .linear_predictor(rows[[1L]], a, h, w, b, ml)
    sdv <- rep(rows[[1L]]$sd, length(a))
  } else {
    # segment k covers ages up to its age_max (ages in the gap between bands go to the next)
    bounds <- vapply(rows[-length(rows)], function(r) r$age_max, numeric(1))
    seg <- findInterval(a, bounds, left.open = TRUE)
    lin <- rep(NA_real_, length(a))
    sdv <- rep(NA_real_, length(a))
    for (k in seq_along(rows)) {
      m <- seg == (k - 1L)
      if (any(m)) {
        lin[m] <- .linear_predictor(rows[[k]], a[m], h[m], w[m], b[m], ml[m])
        sdv[m] <- rows[[k]]$sd
      }
    }
  }
  s <- if (model$dep_transform == "ln_neg") -1 else 1
  tr <- model$dep_transform
  off <- model$dep_offset
  pred <- .inv(lin, tr, off)
  lln <- .inv(lin - s * Z95 * sdv, tr, off)
  uln <- .inv(lin + s * Z95 * sdv, tr, off)
  z <- NULL
  if (!is.null(obs)) {
    z <- s * (.fwd(obs, tr, off) - lin) / sdv
    z[!is.finite(z)] <- NA_real_
  }
  list(pred = pred, lln = lln, uln = uln, z = z)
}

# Quantile regression (Schulz): P50 = predicted, P5 = LLN, P95 = ULN; no z-score.
.quantile <- function(model, rows, a, h, w, b, ml) {
  tr <- model$dep_transform
  off <- model$dep_offset
  list(pred = .inv(.linear_predictor(rows[["50"]], a, h, w, b, ml), tr, off),
       lln  = .inv(.linear_predictor(rows[["5"]],  a, h, w, b, ml), tr, off),
       uln  = .inv(.linear_predictor(rows[["95"]], a, h, w, b, ml), tr, off),
       z = NULL)
}

# ------------------------------------------------------------------- evaluation
# Multiplier converting reference-equation units into the requested result unit.
.scale <- function(state, model) {
  f <- model$pressure_factor
  if (is.na(f)) 1 else f / state$result_factor
}

.display_unit <- function(state, model) {
  nm <- .reference_pressure_name(model$unit)
  if (is.na(nm)) model$unit else paste0(state$result_unit, substring(model$unit, nchar(nm) + 1L))
}

.evaluate <- function(state, param, sex, age, height, weight, value, observed_given) {
  n <- length(age)
  out <- list(predicted = rep(NA_real_, n), percent_predicted = rep(NA_real_, n),
              z_score = rep(NA_real_, n), lln = rep(NA_real_, n), uln = rep(NA_real_, n))
  flag <- rep(NA_character_, n)
  ref <- state$ref

  model <- .find_model(ref, param)
  if (is.null(model)) {
    flag[] <- sprintf("parameter '%s' is not available in reference set '%s' (available: %s)",
                      as.character(param)[1], state$equation,
                      paste(ref$parameters, collapse = ", "))
    return(list(out = out, flag = flag))
  }

  age_y <- age * state$age_factor
  h_cm  <- height * state$height_factor
  w_kg  <- weight * state$weight_factor
  bmi   <- w_kg / (h_cm / 100)^2

  pos <- function(x) !is.na(x) & x > 0
  sex_ok <- sex %in% c("Male", "Female")
  male <- as.numeric(sex %in% "Male")
  age_ok <- pos(age_y)
  h_ok <- pos(h_cm)
  w_ok <- pos(w_kg)
  bmi_ok <- is.finite(bmi)

  flag <- .add_flag(flag, !sex_ok, "sex missing or not recognised")
  flag <- .add_flag(flag, !h_ok, "height missing or not positive")
  flag <- .add_flag(flag, !w_ok, "weight missing or not positive")
  if (model$needs_age) {
    flag <- .add_flag(flag, !age_ok, "age missing or not positive")
  } else {
    flag <- .add_flag(flag, !age_ok, paste0("age missing or not positive (not used by this ",
                                            "equation, but the age range could not be checked)"))
  }

  # ---- validity ranges
  amin <- rep(NA_real_, n)
  amax <- rep(NA_real_, n)
  for (s in c("Male", "Female")) {
    rng <- .age_range(model, s)
    if (!is.null(rng)) {
      m <- sex %in% s
      amin[m] <- rng[1L]
      amax[m] <- rng[2L]
    }
  }
  flag <- .range_flags(flag, "age", "years", ifelse(age_ok, age_y, NA_real_), amin, amax)
  lim <- ref$limits
  flag <- .range_flags(flag, "height", "cm", ifelse(h_ok, h_cm, NA_real_), lim$height[1L], lim$height[2L])
  flag <- .range_flags(flag, "weight", "kg", ifelse(w_ok, w_kg, NA_real_), lim$weight[1L], lim$weight[2L])
  flag <- .range_flags(flag, "BMI", "kg/m2", ifelse(bmi_ok & w_ok & h_ok, bmi, NA_real_),
                       lim$bmi[1L], lim$bmi[2L])

  ok <- sex_ok & h_ok & (w_ok | !model$needs_weight) & (age_ok | !model$needs_age)
  scale <- .scale(state, model)
  obs <- if (is.null(value)) NULL else value / scale     # measured value in reference units

  if (!model$has_sd) {
    flag <- .add_flag(flag, rep(TRUE, n), paste0(
      "z-score not available: ", state$equation, " equations are quantile-based and ",
      "provide no residual standard deviation"))
  }

  sexes <- if ("Both" %in% names(model$rows)) "Both" else c("Male", "Female")
  for (s in sexes) {
    m <- if (s == "Both") ok else ok & (sex %in% s)
    if (!any(m)) next
    a <- age_y[m]; h <- h_cm[m]; w <- w_kg[m]; b <- bmi[m]; ml <- male[m]
    o <- if (is.null(obs)) NULL else obs[m]
    rows <- model$rows[[s]]
    res <- if (model$method == "quantile") {
      .quantile(model, rows, a, h, w, b, ml)
    } else {
      .parametric(model, rows, a, h, w, b, ml, o)
    }
    out$predicted[m] <- res$pred * scale
    out$lln[m] <- res$lln * scale
    out$uln[m] <- res$uln * scale
    if (!is.null(res$z)) out$z_score[m] <- res$z
    if (!is.null(value)) {
      pp <- 100 * value[m] / (res$pred * scale)
      pp[!is.finite(pp)] <- NA_real_
      out$percent_predicted[m] <- pp
    }
  }

  if (observed_given) {
    miss <- !is.finite(value)
    flag <- .add_flag(flag, miss & ok,
                      "measured value missing: percent predicted and z-score not calculated")
    if (model$dep_transform != "none" && model$has_sd) {
      t <- .fwd(obs, model$dep_transform, model$dep_offset)
      bad <- is.finite(value) & !is.finite(t) & ok
      flag <- .add_flag(flag, bad, paste0("measured value outside the domain of the ",
                                          "log-transformed equation: z-score not calculated"))
    }
  }
  list(out = out, flag = flag)
}

# ------------------------------------------------------------------ public methods
.as_num_scalar <- function(x) if (is.null(x)) NA_real_ else suppressWarnings(as.numeric(x))

.calculate <- function(state, age, height, weight, sex, param, value = NULL) {
  for (x in list(age, height, weight, value)) {
    if (!is.null(x) && length(x) != 1L) {
      stop("calculate() takes scalar inputs; use compute() for many rows.", call. = FALSE)
    }
  }
  if (is.null(sex)) sex <- NA
  ev <- .evaluate(state, param, .normalise_sex(sex), .as_num_scalar(age), .as_num_scalar(height),
                  .as_num_scalar(weight), if (is.null(value)) NULL else .as_num_scalar(value),
                  observed_given = !is.null(value))
  model <- .find_model(state$ref, param)
  as_na <- function(v) if (is.nan(v)) NA_real_ else v
  data.frame(
    predicted = as_na(ev$out$predicted), percent_predicted = as_na(ev$out$percent_predicted),
    z_score = as_na(ev$out$z_score), lln = as_na(ev$out$lln), uln = as_na(ev$out$uln),
    flag = ev$flag, param = if (is.null(model)) as.character(param)[1] else model$param_key,
    reference = state$equation, unit = if (is.null(model)) "" else .display_unit(state, model),
    stringsAsFactors = FALSE
  )
}

.compute <- function(state, df, param, sex_col, age_col, height_col, weight_col,
                     param_col, prefix) {
  for (col in c(sex_col, age_col, height_col, weight_col)) {
    if (!(col %in% names(df))) {
      stop(sprintf("Column '%s' not found in data frame (columns: %s).", col,
                   paste(names(df), collapse = ", ")), call. = FALSE)
    }
  }
  if (is.null(param_col) && param %in% names(df)) param_col <- param
  if (!is.null(param_col) && !(param_col %in% names(df))) {
    stop(sprintf("param_col '%s' not found in data frame.", param_col), call. = FALSE)
  }
  num <- function(col) {
    x <- df[[col]]
    if (is.numeric(x)) as.numeric(x) else suppressWarnings(as.numeric(as.character(x)))
  }
  ev <- .evaluate(state, param, .normalise_sex(df[[sex_col]]), num(age_col), num(height_col),
                  num(weight_col), if (is.null(param_col)) NULL else num(param_col),
                  observed_given = !is.null(param_col))
  res <- data.frame(ev$out, stringsAsFactors = FALSE, check.names = FALSE)
  res$flag <- ev$flag
  names(res) <- paste0(prefix, names(res))
  rownames(res) <- row.names(df)
  res
}

#' Create an oscillometry reference
#'
#' @param equation One of "moitra2020", "oostveen2013", "schulz2013", "berger2021",
#'   "gochicoa2023". Case, spaces, underscores and the year are optional.
#' @param age_unit "years", "months", "weeks" or "days".
#' @param height_unit "cm", "m", "mm" or "in".
#' @param weight_unit "kg", "g" or "lb".
#' @param result_unit Pressure unit of measured values and all returned results:
#'   "kPa", "hPa" or "cmH2O".
#' @return An object with members `calculate()`, `compute()`, `parameters`, `metadata`.
#' @export
OscillometryReference <- function(equation, age_unit = "years", height_unit = "cm",
                                  weight_unit = "kg", result_unit = "kPa") {
  refs <- .load_references()
  eq <- .resolve_equation(equation, names(refs))
  au <- .parse_age_unit(age_unit)
  hu <- .parse_height_unit(height_unit)
  wu <- .parse_weight_unit(weight_unit)
  ru <- .parse_result_unit(result_unit)
  state <- list(ref = refs[[eq]], equation = eq,
                age_unit = au$name, age_factor = au$factor,
                height_unit = hu$name, height_factor = hu$factor,
                weight_unit = wu$name, weight_factor = wu$factor,
                result_unit = ru$name, result_factor = ru$factor)
  obj <- list(
    equation = eq, age_unit = au$name, height_unit = hu$name,
    weight_unit = wu$name, result_unit = ru$name,
    parameters = state$ref$parameters,
    metadata = state$ref$metadata,
    calculate = function(age, height, weight, sex, param, value = NULL) {
      .calculate(state, age, height, weight, sex, param, value)
    },
    compute = function(df, param, sex_col = "sex", age_col = "age", height_col = "height",
                       weight_col = "weight", param_col = NULL, prefix = "") {
      .compute(state, df, param, sex_col, age_col, height_col, weight_col, param_col, prefix)
    }
  )
  class(obj) <- "OscillometryReference"
  obj
}

#' List the available reference equation sets
#' @return Character vector of equation ids.
#' @export
available_equations <- function() sort(names(.load_references()), method = "radix")

#' @export
print.OscillometryReference <- function(x, ...) {
  cat(sprintf(paste0("<OscillometryReference equation='%s' age_unit='%s' height_unit='%s' ",
                     "weight_unit='%s' result_unit='%s'>\n"),
              x$equation, x$age_unit, x$height_unit, x$weight_unit, x$result_unit))
  invisible(x)
}
