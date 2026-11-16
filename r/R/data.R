# Load the bundled reference-equation tables into ready-to-evaluate structures.

.pkg_env <- new.env(parent = emptyenv())

.coef_cols <- c("sex_coef", "age_coef", "height_coef", "weight_coef", "bmi_coef",
                "inv_height_coef", "inv_bmi_coef")

# NA or an explicit 0 both mean "term not used"
.opt <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  if (length(v) != 1L || is.na(v) || v == 0) NA_real_ else v
}

.num <- function(v) {
  v <- suppressWarnings(as.numeric(v))
  if (length(v) == 0L) NA_real_ else v
}

.read_extdata <- function(name, ...) {
  path <- system.file("extdata", name, package = "rios", mustWork = TRUE)
  utils::read.csv(path, stringsAsFactors = FALSE, fileEncoding = "UTF-8", ...)
}

# One published equation (one sex / percentile / age segment), r = one-row data.frame
.make_row <- function(r) {
  sd <- .num(r$rmse)
  if (is.na(sd)) sd <- .num(r$rsd)
  age_tf <- r$age_transform
  row <- list(
    intercept = .num(r$intercept),
    sd = sd,
    age_transform = if (is.na(age_tf)) "none" else as.character(age_tf),
    height_unit = as.character(r$height_unit),
    age_min = .num(r$age_min),
    age_max = .num(r$age_max),
    segment = .num(r$segment),
    percentile = .num(r$percentile)
  )
  for (cc in .coef_cols) row[[cc]] <- .opt(r[[cc]])
  row
}

.flat_rows <- function(model) {
  out <- list()
  for (sx in names(model$rows)) out <- c(out, unname(model$rows[[sx]]))
  out
}

.limits <- function(meta) {
  hf <- if (tolower(meta$input_height_unit) == "m") 100 else 1
  list(
    height = c(.num(meta$height_min), .num(meta$height_max)) * hf,
    weight = c(.num(meta$weight_min), .num(meta$weight_max)),
    bmi    = c(.num(meta$bmi_min), .num(meta$bmi_max))
  )
}

.param_key <- function(parameter, freq) {
  ifelse(parameter == "AX" & !is.na(freq), paste0("AX", as.integer(freq)), parameter)
}

.build_model <- function(pid, key, grp) {
  first <- grp[1L, ]
  freq <- .num(first$frequency_hz)
  model <- list(
    reference_id = pid, param_key = key, parameter = as.character(first$parameter),
    frequency_hz = freq, unit = as.character(first$parameter_unit),
    method = as.character(first$equation_method),
    dep_transform = as.character(first$dependent_transform),
    dep_offset = .num(first$dependent_offset),
    rows = list(), ranges = list(), needs_age = FALSE, needs_weight = FALSE
  )
  model$pressure_factor <- .reference_pressure_factor(model$unit)
  model$has_sd <- model$method != "quantile"

  for (sx in unique(grp$sex)) {
    sgrp <- grp[grp$sex == sx, ]
    rows <- lapply(seq_len(nrow(sgrp)), function(i) .make_row(sgrp[i, ]))
    if (model$method == "quantile") {
      names(rows) <- as.character(vapply(rows, function(r) r$percentile, numeric(1)))
    } else {
      seg <- vapply(rows, function(r) if (is.na(r$segment)) 0 else r$segment, numeric(1))
      rows <- rows[order(seg)]
    }
    model$rows[[sx]] <- rows
    model$ranges[[sx]] <- c(min(vapply(rows, function(r) r$age_min, numeric(1))),
                            max(vapply(rows, function(r) r$age_max, numeric(1))))
  }
  flat <- .flat_rows(model)
  has_coef <- function(nm) any(vapply(flat, function(r) !is.na(r[[nm]]), logical(1)))
  model$needs_age <- model$method == "segmented" || has_coef("age_coef")
  model$needs_weight <- has_coef("weight_coef") || has_coef("bmi_coef") || has_coef("inv_bmi_coef")
  model
}

.load_references <- function() {
  if (!is.null(.pkg_env$refs)) return(.pkg_env$refs)
  eq <- .read_extdata("reference_equations.csv", na.strings = c("NA", ""))
  meta <- .read_extdata("reference_metadata.csv", colClasses = "character",
                        na.strings = c("NA", ""))
  if ("output_unit" %in% names(eq) && !("parameter_unit" %in% names(eq))) {
    names(eq)[names(eq) == "output_unit"] <- "parameter_unit"
  }
  eq$paper_id <- sub("_.*$", "", eq$reference_id)
  eq$param_key <- .param_key(eq$parameter, .num(eq$frequency_hz))

  refs <- list()
  for (i in seq_len(nrow(meta))) {
    rec <- as.list(meta[i, , drop = FALSE])
    pid <- rec$reference_id
    for (nm in c("year", "age_min", "age_max", "height_min", "height_max", "weight_min",
                 "weight_max", "bmi_min", "bmi_max")) rec[[nm]] <- .num(rec[[nm]])
    sub <- eq[eq$paper_id == pid, ]
    models <- list()
    for (key in unique(sub$param_key)) {
      models[[key]] <- .build_model(pid, key, sub[sub$param_key == key, ])
    }
    refs[[pid]] <- list(
      reference_id = pid, metadata = rec, models = models, limits = .limits(rec),
      parameters = sort(names(models), method = "radix"),
      upper = setNames(as.list(names(models)), toupper(names(models)))
    )
  }
  .pkg_env$refs <- refs
  refs
}

# Look up a model by (flexible) parameter name; NULL if the set does not have it.
.find_model <- function(ref, param) {
  key <- .normalise_param(param)
  hit <- ref$upper[[key]]
  if (is.null(hit) && key == "AX") hit <- ref$upper[["AX5"]]   # default AX = integrated from 5 Hz
  if (is.null(hit)) NULL else ref$models[[hit]]
}

.age_range <- function(model, sex) {
  r <- model$ranges[[sex]]
  if (is.null(r)) r <- model$ranges[["Both"]]
  r
}
