#' Case Study: Copula-Graphic Estimator applied to Melanoma data
#'
#' @param sample_data Competing risk data. If `NA` is passed, the function will
#' load data from the `data_path` or use the package data `Melanoma`.
#' @param data_path Path to file with competing risk data. If `NA` is passed as input
#' argument for `sample_data`, the function will use the `load_data` function to load
#' the competing risk data from `data_path`. If `NA` is passed as input argument for
#' `data_path`, the package data `Melanoma` is loaded.
#' @param tau Kendall's \eqn{\tau}. This parameter is used to specify Kendall's \eqn{\tau}
#' for the copulas that are used to calculate the Copula-Graphic estimator.
#' @param error_A Error tolerance for the first iteration step. The estimate for
#' \eqn{G(t[i])} is accepted if \eqn{|mu_C(A(t[i])) -k(t[i])| <} `error_A`.
#' @param error_B Error tolerance for the second iteration step. The estimated
#' pair \eqn{(F(t[i]), G(t[i]))} is accepted if
#' \eqn{|mu_C(B_t[i])) - p_1(t[i])| <} `error_B`.
#' @param p_1 Empirical estimate est \eqn{p_1(t[i])} for \eqn{P(X \le t[i], \delta = 1)}.
#' See documentation of function \link[copulagraphicr]{estimate_quantities}. If `NA` is passed,
#' the function `estimate_quantities`is called to calculate `p_1`.
#' @param k Empirical estimate est \eqn{k(t[i])} for \eqn{P(X > t[i], Y > t[i])}.
#' See documentaion of function \link[copulagraphicr]{estimate_quantities}. If `NA` is passed,
#' the function `estimate_quantities`is called to calculate `k`.
#' @param t_grid Time grid of consisting of the unique event times \eqn{T}.
#'
#' @export
#'
#' @examples ## use existing data from the einvironment as input data
#' # CG_results(sample_data = sample_data, load_data = TRUE, data_path = NA)
#'
#' ## use package data as input data
#' # CG_results(sample_data = NA, load_data = TRUE, data_path = NA)
#'
#' ## load data from data_path -> pass input for data_path
#' # CG_results(sample_data = NA, load_data = TRUE, data_path = data_path)
CG_results <- function(sample_data = NA, data_path = NA, tau = 0.5, error_A = 1e-6, error_B = 1e-5,
                       k = NA, p_1 = NA, t_grid = NA) {

## if no sample_data is loaded, load sample data
if (is.na(sample_data[[1]][1])) {
  if (is.na(data_path)) {
    sample_data <- copulagraphicr::load_data()
  } else if (!is.na(data_path)) {
    sample_data <- copulagraphicr::load_data(data_path = data_path)
  }
}

## if no estimates are loaded, calculate empirical quantities
if (is.na(k[1]) | is.na(p_1[1]) | is.na(t_grid[1])) {
  est <- copulagraphicr::estimate_quantities(sample_data = sample_data)

  t_grid <- est[[1]]
  k <- est[[2]]
  p_1 <- est[[3]]
  p_2 <- est[[4]]
}

## calculate Copula-Graphic estimator
est_Frank <- copulagraphicr::CG_fit(k = k, p_1 = p_1, t_grid = t_grid, copula = copulagraphicr::C_Frank, error_A = error_A, error_B = error_B, tau = tau)
est_Gumbel <- copulagraphicr::CG_fit(k = k, p_1 = p_1, t_grid = t_grid, copula = copulagraphicr::C_Gumbel, error_A = error_A, error_B = error_B, tau = tau)
est_Clayton <- copulagraphicr::CG_fit(k = k, p_1 = p_1, t_grid = t_grid, copula = copulagraphicr::C_Clayton, error_A = error_A, error_B = error_B, tau = tau)
est_Joe <- copulagraphicr::CG_fit(k = k, p_1 = p_1, t_grid = t_grid, copula = copulagraphicr::C_Joe, error_A = error_A, error_B = error_B, tau = tau)
est_Ind <- copulagraphicr::CG_fit(k = k, p_1 = p_1, t_grid = t_grid, copula = copulagraphicr::C_Independence, error_A = error_A, error_B = error_B, tau = tau)
if (tau == 0.5) {
  est_Gamma <- copulagraphicr::CG_fit(k = k, p_1 = p_1, t_grid = t_grid, copula = copulagraphicr::C_Gamma_Frailty, error_A = error_A, error_B = error_B, tau = tau)
  S_Gamma <- 1-est_Gamma[[1]]
}


S_Frank <- 1 - est_Frank[[1]]
S_Gumbel <- 1 - est_Gumbel[[1]]
S_Clayton <- 1 - est_Clayton[[1]]
S_Joe <- 1 -est_Joe[[1]]
S_Ind <- 1 - est_Ind[[1]]

## Benchmark: survival fit, Kaplan-Meier for censored data
d <- rep(0, nrow(sample_data))
x <- sample_data$T
e <- sample_data$delta
fit <- survival::survfit(survival::Surv(d, x, e) ~ 1)

## Plot
y_lim_lower <- round(min(S_Frank, S_Gumbel, S_Clayton, S_Joe, S_Ind) - 0.15, 1)

grDevices::dev.new()
coord <- graphics::par("usr")
graphics::plot(t_grid[1:length(S_Joe)],S_Joe,type="s",ylim = c(y_lim_lower,1), xlab="time", ylab="survival function", col = "white")
graphics::lines(fit, col = "red", conf.int = FALSE, type="s")
graphics::lines(t_grid[1:length(S_Frank)], S_Frank, col = "orange", type="s")
graphics::lines(t_grid[1:length(S_Gamma)], S_Gamma, col = "green", type="s")
graphics::lines(t_grid[1:length(S_Gumbel)], S_Gumbel, col = "purple", type="s")
graphics::lines(t_grid[1:length(S_Clayton)], S_Clayton, col = "blue", type="s")
graphics::lines(t_grid[1:length(S_Joe)], S_Joe, col = "brown", type="s")
graphics::legend(0, y_lim_lower+0.25, legend = c("Independence", "Frank's", "Gamma Frailty", "Gumbel", "Clayton", "Joe"),
       c("red", "orange", "green", "purple", "blue", "brown"), cex = 1, box.lty = 0, lty = 1,
       fill = 0, border = 0, x.intersp = 0.2, y.intersp = 0.8)

return(list(S_Frank, S_Gumbel, S_Clayton, S_Joe, S_Ind))
}
