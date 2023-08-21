#' A S4 class containing the copula graphic object
#'
#' @slot sample_data
#'
setClass(
  "copula_graphic",
  slots = c(
    sample_data = "data.frame",
    empirical_quantities = "list",
    param_frank = "list",
    est_frank = "list",
    surv_frank = "data.frame",
    param_gumbel = "list",
    est_gumbel = "list",
    surv_gumbel = "data.frame",
    param_joe = "list",
    est_joe = "list",
    surv_joe = "data.frame",
    param_ind = "list",
    est_ind = "list",
    surv_ind = "data.frame",
    param_kaplan_meier = "list",
    kaplan_meier = "data.frame"
  ),
  prototype = prototype(
    sample_data = copulagraphicr::load_data(data_path = NA),
    empirical_quantities = list(
      t_grid = NA_real_,
      k = NA_real_,
      p_1 = NA_real_,
      p_2 = NA_real_
    ),
    param_frank = list(
      copula = copulagraphicr::C_Frank,
      tau = NA
    ),
    est_frank = list(),
    surv_frank = data.frame(),
    param_gumbel = list(
      copula = copulagraphicr::C_Gumbel,
      tau = NA
    ),
    est_gumbel = list(),
    surv_gumbel = data.frame(),
    param_joe = list(
      copula = copulagraphicr::C_Joe,
      tau = NA
    ),
    est_joe = list(),
    surv_joe = data.frame(),
    param_ind = list(),
    est_ind = list(),
    surv_ind = data.frame(),
    param_kaplan_meier = list(),
    kaplan_meier = data.frame()
  )
)
