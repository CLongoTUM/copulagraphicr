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
#' @param data_path Path to file. Files with extension '.csv' are supported.
#'
#' @importFrom magrittr "%>%"
#' @export
#'
#' @examples \donttest{# load data set:
#' df <- load_data()
#' }
load_data <- function(data_path = system.file("inst/Melanoma.csv", package = "copulagraphicr")) {
  sample_data <- data.table::fread(data_path) %>%
    dplyr::select("time", "delta") %>%
    dplyr::rename("T" = "time")
  return(sample_data)
}



