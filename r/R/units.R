# Unit handling for inputs (age, height, weight) and oscillometry results.

.age_units   <- c(years = 1, months = 1 / 12, weeks = 7 / 365.25, days = 1 / 365.25)
.age_aliases <- c(y = "years", yr = "years", yrs = "years", year = "years", years = "years",
                  mo = "months", mos = "months", month = "months", months = "months",
                  w = "weeks", wk = "weeks", wks = "weeks", week = "weeks", weeks = "weeks",
                  d = "days", day = "days", days = "days")

.height_units   <- c(cm = 1, m = 100, mm = 0.1, "in" = 2.54)
.height_aliases <- c(cm = "cm", centimetre = "cm", centimetres = "cm",
                     centimeter = "cm", centimeters = "cm",
                     m = "m", metre = "m", metres = "m", meter = "m", meters = "m",
                     mm = "mm", millimetre = "mm", millimetres = "mm",
                     millimeter = "mm", millimeters = "mm",
                     "in" = "in", inch = "in", inches = "in")

.weight_units   <- c(kg = 1, g = 0.001, lb = 0.45359237)
.weight_aliases <- c(kg = "kg", kgs = "kg", kilogram = "kg", kilograms = "kg",
                     kilo = "kg", kilos = "kg",
                     g = "g", gram = "g", grams = "g",
                     lb = "lb", lbs = "lb", pound = "lb", pounds = "lb")

# pressure units -> kPa
.pressure_units   <- c(kPa = 1, hPa = 0.1, cmH2O = 0.0980665)
.pressure_aliases <- c(kpa = "kPa", hpa = "hPa", mbar = "hPa", millibar = "hPa",
                       cmh2o = "cmH2O")

.unit_key <- function(value) {
  s <- tolower(trimws(as.character(value)[1]))
  s <- gsub("\u2082", "2", s, fixed = TRUE)          # subscript two in cmH2O
  gsub("[^a-z0-9]", "", s)
}

# Returns list(name = canonical unit, factor = multiplier to the internal standard unit)
.parse_unit <- function(value, aliases, table, kind) {
  canon <- unname(aliases[.unit_key(value)])
  if (length(canon) != 1L || is.na(canon)) {
    stop(sprintf("Unsupported %s unit '%s'. Supported: %s.", kind,
                 as.character(value)[1], paste(sort(names(table)), collapse = ", ")),
         call. = FALSE)
  }
  list(name = canon, factor = unname(table[canon]))
}

.parse_age_unit    <- function(value) .parse_unit(value, .age_aliases, .age_units, "age")
.parse_height_unit <- function(value) .parse_unit(value, .height_aliases, .height_units, "height")
.parse_weight_unit <- function(value) .parse_unit(value, .weight_aliases, .weight_units, "weight")
.parse_result_unit <- function(value) .parse_unit(value, .pressure_aliases, .pressure_units, "result")

# Pressure part of a reference unit string such as "hPa*s/L"; NA for Hz or %.
.reference_pressure_name <- function(unit) {
  m <- regmatches(unit, regexpr("^(kPa|hPa|cmH2O)", unit))
  if (length(m) == 1L) m else NA_character_
}

.reference_pressure_factor <- function(unit) {
  nm <- .reference_pressure_name(unit)
  if (is.na(nm)) NA_real_ else unname(.pressure_units[nm])
}
