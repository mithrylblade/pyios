\name{OscillometryReference}
\alias{OscillometryReference}
\alias{print.OscillometryReference}
\title{Create an oscillometry reference}
\description{
Creates an object that calculates predicted values, percent predicted, z-scores and lower and
upper limits of normal (LLN, ULN) for impulse oscillometry parameters using one of five
published reference equation sets.
}
\usage{
OscillometryReference(equation, age_unit = "years", height_unit = "cm",
                      weight_unit = "kg", result_unit = "kPa")
}
\arguments{
  \item{equation}{One of \code{"moitra2020"}, \code{"oostveen2013"}, \code{"schulz2013"},
    \code{"berger2021"}, \code{"gochicoa2023"}. Case, spaces, underscores and the year are
    optional, so \code{"Gochicoa"}, \code{"gochicoa_2023"} and \code{"GOCHICOA 2023"} are
    equivalent.}
  \item{age_unit}{Unit of supplied ages: \code{"years"}, \code{"months"}, \code{"weeks"} or
    \code{"days"}.}
  \item{height_unit}{Unit of supplied heights: \code{"cm"}, \code{"m"}, \code{"mm"} or
    \code{"in"}.}
  \item{weight_unit}{Unit of supplied weights: \code{"kg"}, \code{"g"} or \code{"lb"}.}
  \item{result_unit}{Pressure unit of measured values passed in and of all returned results:
    \code{"kPa"}, \code{"hPa"} or \code{"cmH2O"}. Non-pressure parameters (Fres in Hz, ratios
    in percent) are never scaled.}
}
\value{
An object of class \code{OscillometryReference}, a list with members:
\describe{
  \item{\code{calculate(age, height, weight, sex, param, value = NULL)}}{Evaluates one subject
    and returns a one-row data frame with \code{predicted}, \code{percent_predicted},
    \code{z_score}, \code{lln}, \code{uln}, \code{flag}, \code{param}, \code{reference} and
    \code{unit}. \code{value} is the optional measured result (in \code{result_unit});
    without it \code{percent_predicted} and \code{z_score} are \code{NA}.}
  \item{\code{compute(df, param, sex_col = "sex", age_col = "age", height_col = "height",
    weight_col = "weight", param_col = NULL, prefix = "")}}{Evaluates every row of a data
    frame and returns a data frame (same row names) with columns \code{predicted},
    \code{percent_predicted}, \code{z_score}, \code{lln}, \code{uln} and \code{flag}.
    \code{param_col} names the column of measured values; if \code{NULL} and \code{df} has a
    column named \code{param}, that column is used.}
  \item{\code{parameters}}{Parameters available in the selected reference set.}
  \item{\code{metadata}}{Study metadata: population, country, validity ranges, device, DOI.}
}
Unavailable values are \code{NA}. \code{flag} is \code{NA} when everything is valid, and
otherwise a message such as "age 30 years below reference minimum 45 years". Sex may be given
as 0/"0"/"f"/"female" (female) or 1/"1"/"m"/"male" (male), in any case.
}
\examples{
reference <- OscillometryReference(equation = "gochicoa2023")
reference$calculate(age = 50, height = 175, weight = 80, sex = "male", param = "R5",
                    value = 0.30)

df <- data.frame(sex = c("male", "female"), age = c(50, 62), height = c(175, 160),
                 weight = c(80, 60), R5 = c(0.30, 0.35))
reference$compute(df, param = "R5", param_col = "R5")
}
