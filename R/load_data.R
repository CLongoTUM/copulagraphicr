library("riskRegression")
library("dplyr")
library("lubridate")
library("magrittr")
load_data <- function() {

  load(file = paste0("C:/Users/", Sys.getenv("USERNAME"), "/OneDrive/Documents/R (git)/Projects/projectr/inst/Melanoma.rda"))
  # system.file("Melanoma.rda")

  ## competing risk data
  df <- Melanoma

  # time: time in days from operation, stopped after "event" occurred
  # status: 1 = died from melanoma; 0 = censored

  # risk: X (time to death), Y (time to censoring)
  # observations (T_i, delta_i), delta = 1 if X =< Y
  sample_data <- data.table::data.table(T = df$time, delta = df$status)
}


