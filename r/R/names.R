# Flexible parsing of equation names, parameter names and sex codes.

.equation_aliases <- c(shulz = "schulz", gochicoarangel = "gochicoa")

# 'Gochicoa_2023', 'gochicoa 2023', 'GOCHICOA-2023' -> 'gochicoa2023'
.normalise_equation_key <- function(name) {
  key <- gsub("[^a-z0-9]", "", tolower(as.character(name)[1]))
  for (wrong in names(.equation_aliases)) {
    if (startsWith(key, wrong)) {
      key <- paste0(.equation_aliases[[wrong]], substring(key, nchar(wrong) + 1L))
    }
  }
  key
}

.resolve_equation <- function(name, available) {
  key <- .normalise_equation_key(name)
  if (key %in% available) return(key)
  matches <- if (nzchar(key)) available[startsWith(available, key)] else character(0)
  if (length(matches) == 1L) return(matches)
  if (length(matches) > 1L) {
    stop(sprintf("Equation '%s' is ambiguous (matches %s).", name,
                 paste(matches, collapse = ", ")), call. = FALSE)
  }
  stop(sprintf("Unknown equation '%s'. Available equations: %s.", name,
               paste(sort(available, method = "radix"), collapse = ", ")), call. = FALSE)
}

# Upper-case canonical parameter name used for matching:
# 'r5', 'R5-20', 'r5_r20', 'R5<en dash>R20' -> 'R5-R20'; 'R5-R20/R20' -> '(R5-R20)/R20'
.normalise_param <- function(param) {
  s <- trimws(as.character(param)[1])
  for (d in c("\u2010", "\u2011", "\u2012", "\u2013", "\u2014", "\u2212")) {
    s <- gsub(d, "-", s, fixed = TRUE)
  }
  s <- toupper(gsub("_", "-", gsub("[[:space:]]+", "", s), fixed = TRUE))
  s <- sub("^R([0-9]+)-([0-9]+)$", "R\\1-R\\2", s)
  s <- sub("^\\(?R5-R(15|20)\\)?/R(5|20)(\\*100)?$", "(R5-R\\1)/R\\2", s)
  s
}

.female_codes <- c("0", "f", "female", "w", "woman", "girl")
.male_codes   <- c("1", "m", "male", "man", "boy")

# One value -> "Male", "Female" or NA
.normalise_sex_one <- function(v) {
  if (is.null(v) || length(v) != 1L || is.na(v)) return(NA_character_)
  if (is.logical(v)) return(NA_character_)
  if (is.numeric(v)) {
    if (v == 0) return("Female")
    if (v == 1) return("Male")
    return(NA_character_)
  }
  s <- tolower(trimws(as.character(v)))
  if (grepl("^-?[0-9]+\\.0+$", s)) s <- sub("\\..*$", "", s)
  if (s %in% .female_codes) return("Female")
  if (s %in% .male_codes) return("Male")
  NA_character_
}

# Vectorised version returning a character vector of "Male"/"Female"/NA
.normalise_sex <- function(x) {
  if (is.factor(x)) x <- as.character(x)
  if (is.null(x)) return(NA_character_)
  u <- unique(x)
  mapped <- vapply(seq_along(u), function(i) .normalise_sex_one(u[[i]]), character(1))
  mapped[match(x, u)]
}
