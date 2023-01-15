#' Data Input: load competing risk data
#'
#' The `load_data` function loads data from `data_path` and extracts
#' the competing risk data consisting of observations (T\eqn{_i}, \eqn{\delta_i}).
#'
#' The competing risk data should consist of the observed sample data
#' (T\eqn{_i}, \eqn{\delta_i}).
#' \itemize{
#' \item `competing risks`:
#' \itemize{
#' \item \eqn{X} (time to death)
#' \item \eqn{Y} (censoring time)
#' }
#' \item `sample data`:
#' \itemize{
#' \item `time to event` \eqn{T =} min(X,Y)
#' \item `event-indicator` \eqn{\delta} = 1 if \eqn{X \le Y}; \eqn{\delta} = 0 if \eqn{X > Y}
#' }
#' }
#'
#' @param data_path Path to file. Files with extension '.csv' are supported. If
#' `NA` is passed, the package data `Melanoma` will be loaded.
#'
#' @importFrom magrittr "%>%"
#' @export
#'
#' @examples \donttest{## load data set:
#' # df <- load_data()
#' }
load_data <- function(data_path = NA) {
  if (is.na(data_path)) {
    sample_data <- copulagraphicr::Melanoma %>%
      dplyr::select("time", "delta") %>%
      dplyr::rename("T" = "time")
    ## alternatively call
    # sample_data <- riskRegression::Melanoma %>%
    #   dplyr::filter(status != 2) %>%
    #   dplyr::select("time", "delta") %>%
    #   dplyr::rename("T" = "time")
  } else {
    sample_data <- data.table::fread(data_path) %>%
    dplyr::select("time", "delta") %>%
    dplyr::rename("T" = "time")
  }
  return(sample_data)
}

#' Melanoma data
#'
#' The study observed a total of 205 patients during the time from 1962 until 1977.
#' At the end 134 were still alive, while 57 had died from cancer and 14 died from other causes.
#' The original dataset contains 12 variables of which we will only consider the following
#' \itemize{
#' \item `time`: time-to-event in days from operation (event: "status change"),
#' \item `status`: numeric with values (0 = censored, 1 = death.melanoma, 2 = death other).
#' }
#' Note that the Copula-Graphic estimator is designed for the case of two competing risks.
#' Hence, only patients with status 0 or 1 are considered in the Melanoma package data.
#' The two competing events are A = "died from Melanoma" and B = "censored". The sample data
#' contains the respective time-to-event \eqn{T} and the information if A or B was observed, captured by \eqn{\delta}.
#' @name Melanoma
#' @docType data
#' @references The data used is taken from riskRegression::Melanoma
#' @keywords data
#' @examples
#' # library(riskRegression)
#' # df <- load(Melanoma)
NULL

