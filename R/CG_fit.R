#' Copula-Graphic Estimator
#'
#' @param p_1 Empirical estimate est \eqn{p_1(t[i])} for \eqn{P(X \le t[i], \delta = 1)}.
#' See documentation of function \link[copulagraphicr]{estimate_quantities}.
#' @param k Empirical estimate est \eqn{k(t[i])} for \eqn{P(X > t[i], Y > t[i])}.
#' See documentaion of function \link[copulagraphicr]{estimate_quantities}.
#' @param error_A Error tolerance for the first iteration step. The estimate for
#' \eqn{G(t[i])} is accepted if \eqn{|mu_C(A(t[i])) -k(t[i])| <} `error_A`.
#' @param error_B Error tolerance for the second iteration step. The estimated
#' pair \eqn{(F(t[i]), G(t[i]))} is accepted if
#' \eqn{|mu_C(B_t[i])) - p_1(t[i])| <} `error_B`.
#' @param t_grid Time grid of consisting of the unique event times \eqn{T}.
#' @param copula The Copula to be assumed for calculating the Copula-Graphic Estimator.
#' @param theta Parameter of the following copulas:
#' \itemize{
#' \item Gamma Frailty Copula, `theta` \eqn{\ge 1}
#' \item Gumbel Copula, `theta` \eqn{\ge 1}
#' \item Frank Copula, `theta` \eqn{=/= 0}
#' \item Clayton Copula, `theta` \eqn{\ge -1}
#' \item Joe's Copula, `theta` \eqn{\ge 1}
#' }
#' If set to `NA`, `theta` is chosen such that
#' Kendall's \eqn{\tau = 0.5}. Note that the parameter `theta` is only used, if `NA` value for
#' `tau` is passed.
#' @param tau Kendall's \eqn{\tau}.
#'
#' @export
#'
#' @examples
#' ## load data set:
#' # sample_data <- copulagraphicr::load_data()
#'
#' ## estimate empirical quantities:
#' # est <- copulagraphicr::estimate_quantities(sample_data = sample_data)
#' # t_grid <- est[[1]]
#' # k <- est[[2]]
#' # p_1 <- est[[3]]
#' # p_2 <- est[[4]]
#'
#' ## Copula-Graphic Estimates for different input copulae:
#' # est_Ind <- copulagraphicr::CG_fit(k = k, p_1 = p_1, t_grid = t_grid, copula = copulagraphicr::C_Independence)
#' # est_Clayton <- copulagraphicr::CG_fit(k = k, p_1 = p_1, t_grid = t_grid, copula = copulagraphicr::C_Gumbel, tau = 0.5)
#' # est_Frank <- copulagraphicr::CG_fit(k = k, p_1 = p_1, t_grid = t_grid, copula = copulagraphicr::C_Frank, tau = 0.5)
CG_fit <- function(p_1,
                   k,
                   t_grid,
                   copula,
                   error_A = 1e-7,
                   error_B = 1e-5,
                   theta = NA,
                   tau = NA) {
  F_hat <- c(0)
  G_hat <- c(0)
  mu_C_B <- c(0)
  mu_C_A <- c(0)
  copula <- copula
    for (i in 2:(length(t_grid))) {
    solve <- copulagraphicr::solve_G_F(F_t_i_minus_1 = F_hat[i-1],
                                       G_t_i_minus_1 = G_hat[i-1],
                                       mu_C_B_t_i_minus_1 = mu_C_B[i-1],
                                       copula = copula,
                                       error_A,
                                       error_B,
                                       i = i,
                                       p_1 = p_1,
                                       k = k,
                                       theta = theta,
                                       tau = tau)
    F_hat[i] <- solve[[1]]
    G_hat[i] <- solve[[2]]
    mu_C_B[i] <- solve[[3]]
    mu_C_A[i] <- solve[[4]]
    print(paste0("Currently evaluating ", i, "-th time point."))
    }
  return(list(F_hat, G_hat))
}


# est_ind <- copulagraphicr::CG_fit(k = k, p_1 = p_1, t_grid = t_grid, copula = copulagraphicr::C_Independence, error_A = 1e-6, error_B = 1e-3)
# est_Clayton <- copulagraphicr::CG_fit(k = k, p_1 = p_1, t_grid = t_grid, copula = copulagraphicr::C_Gumbel, tau = 0.5, error_A = 1e-6, error_B = 1e-3)
# est_Frank <- copulagraphicr::CG_fit(k = k, p_1 = p_1, t_grid = t_grid, copula = copulagraphicr::C_Frank, error_A = 1e-6, error_B = 1e-3)
#
# plot(t_grid, 1-est_ind[[1]], type = "s", ylim = c(0,1))
# lines(t_grid, 1-est_Clayton[[1]], type = "s", ylim = c(0,1), col = "red")
# lines(t_grid, 1-est_Frank[[1]], type = "s", ylim = c(0,1), col = "orange")
