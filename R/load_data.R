#' The "load_data" function loads data from data_path and extracts
#' the competing risk data consisting of observations (T_i, delta_i).
#'
#' @param data_path Path to file. Files with extension '.csv' are supported.
#'
#' @export
#'
#' @examples # load data set:
#' df <- load_data()
load_data <- function(data_path = system.file("inst/Melanoma.csv", package = "copulagraphicr")) {

  ## competing risk data
  # variables: time; delta

  ## variable description
  # time: time in days from operation, stopped after "risk-event" occurred
  # risk-events for example data: X (time to death), Y (time to censoring)
  # delta: indicator, in risk-event was observed or censored, 1 if X =< Y; 0 if X > Y

  # observations (T_i, delta_i)
  sample_data <- data.table::fread(data_path) %>%
    select("time", "delta")
}



