#' Iterative scheme for calculation of Copula-Graphic Estimator
#'
#' @description Two step iterative bisection scheme to solve the following system of non-linear equations:
#' \itemize{
#' \item (1) \eqn{mu_C(A_t[i]) = k(t[i]),}
#' \item (2) \eqn{mu_C(B_t[i]) =  p_1(t[i]).}
#' }
#' Notation:
#' \itemize{
#' \item \eqn{A_t[i] = { (x,y) : F(t[i]) \le x \le 1, G(t[i]) \le y \le 1 } }
#' \item \eqn{B_t[i] = { (x,y) : 0 \le F(t[i]) \le 1, GF^{-1}(t[i]) \le y \le 1 } }
#' \item \eqn{k(t[i]) = 1/n \sum_k 1_(X_k > t[i], Y_k > t[i])}
#' \item \eqn{p_1(t[i]) = 1/n \sum_k 1_(X_k \le t[i], \delta_k = 1)}
#' }
#' @description It can be shown that
#' \itemize{
#' \item \eqn{ mu_C(A_t[i]) = 1 - G(t[i]) - F(t[i]) + C(F(t[i]), G(t[i]))} and
#' \item \eqn{mu_C(B_t[i]) = mu_C(B_t[i-1]) + F(t[i]) - F(t[i-1]) + C(F(t[i-1]), G(t[i])) - C(F(t[i]), G(t[i]))}.
#' }
#' Note: in the i-th iteration step, the variables \eqn{F(t[i-1])} and \eqn{mu_C(B_t[i-1])} are known input parameters.
#' Hence, there exists a unique solution, since there are two equations and two unknown variables.
#' @description In the first step a candidate \eqn{F_t_i} for \eqn{F(t[i])} is proposed. Then the
#' first equation is solved for the corresponding candidate \eqn{G_t_i} for \eqn{G(t[i])}.
#' The equation is solved using a bisection algorithm and the candidate \eqn{G_t_i} is accepted, if
#' the error \eqn{|mu_C(A_t[i]) - k(t[i])|} is smaller than `error_A`. Note that for a bad initial choice of
#' \eqn{F_t_i}, it may happen that there is no solution. In this case, the algorithm stops after \eqn{iter_max_A}
#' steps and a new candidate \eqn{F_t_i} is selected.
#' @description In the second step the pair \eqn{(F_t_i, G_t_i)} is proposed as a solution for
#' the second equation. If \eqn{|mu_C(B_t[i]) - p_1(t[i])|} is smaller than `error_B`, the pair is accepted
#' as a solution. Else a new candidate for \eqn{F_t_i} is selected.
#'
#' @param F_t_i_minus_1 The estimate for \eqn{F(t[i-1])}. This input
#' parameter is given by the previous iteration step.
#' @param G_t_i_minus_1 The estimate for \eqn{G(t[i-1])}. This input
#' parameter is given by the previous iteration step.
#' @param mu_C_B_t_i_minus_1 The estimate for \eqn{mu_C(B_t[i-1]))}.
#' This input parameter is given by the previous iteration step.
#' @param copula Copula function to be applied for the Copula-Graphic Estimator.
#' @param error_A Error tolerance for the first iteration step. The estimate for
#' \eqn{G(t[i])} is accepted if \eqn{|mu_C(A(t[i])) -k(t[i])| <} `error_A`.
#' @param error_B Error tolerance for the second iteration step. The estimated
#' pair \eqn{(F(t[i]), G(t[i]))} is accepted if
#' \eqn{|mu_C(B_t[i])) - p_1(t[i])| <} `error_B`.
#' @param i Iteration parameter.
#' @param p_1 Empirical estimate est \eqn{p_1(t[i])} for \eqn{P(X \le t[i], \delta = 1)}.
#' See documentation of function \link[copulagraphicr]{estimate_quantities}.
#' @param k Empirical estimate est \eqn{k(t[i])} for \eqn{P(X > t[i], Y > t[i])}.
#' See documentaion of function \link[copulagraphicr]{estimate_quantities}.
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
#' @examples \donttest{## load data set:
#' # sample_data <- copulagraphicr::load_data()
#'
#' ## estimate empirical quantities:
#' # est <- copulagraphicr::estimate_quantities(sample_data = sample_data)
#' # t_grid <- est[[1]]
#' # k <- est[[2]]
#' # p_1 <- est[[3]]
#' # p_2 <- est[[4]]
#' }
#'
#' # F_hat <- c(0)
#' # G_hat <- c(0)
#' # mu_C_B <- c(0)
#' # mu_C_A <- c(0)
#' #
#' # solve <- copulagraphicr::solve_G_F(F_t_i_minus_1 = F_hat[i-1],
#' #                      G_t_i_minus_1 = G_hat[i-1],
#' #                      mu_C_B_t_i_minus_1 = mu_C_B[i-1],
#' #                      copula = copulagraphicr::C_Independence,
#' #                      i = 2,
#' #                      p_1 = p_1,
#' #                      k = k)
#' # F_hat[i] <- solve[[1]]
#' # G_hat[i] <- solve[[2]]
#' # mu_C_B[i] <- solve[[3]]
#' # mu_C_A[i] <- solve[[4]]
#' @export
#'
solve_G_F <- function(F_t_i_minus_1,
                      G_t_i_minus_1,
                      mu_C_B_t_i_minus_1,
                      copula,
                      error_A = 1e-7,
                      error_B = 1e-5,
                      i,
                      p_1,
                      k,
                      theta = NA,
                      tau = NA) {
  # bisection algorithm to find a solution (F(t_i), G(t_i))

  ## The bisection chooses a parameter in the interval [lower_F,upper_F].
  ## In every bisection step the parameter interval length is divided by 2.
  ## Hence, after n steps, the interval length is |upper_F - lower_F| = o(2^(-n)).
  ## Assuming that if, |upper_F - lower_F| < error_A / 10^3, making the interval
  ## any smaller, won't improve the precision of the estimate, we choose the
  ## max iteration as iter_max_A <- log(error_A/10^3)/log(1/2), which is a
  ## solution for: |upper_F - lower_F| < 2^(-iter_max) < error_A / 10^3.
  ## Similary, we choose iter_max_B <- log(error_B/10^3)/log(1/2).
  iter_max_A <- log(error_A/10^3)/log(1/2)
  iter_max_B <- log(error_B/10^3)/log(1/2)

  # initialization
  upper_F <- 1 # F_t =< 1
  lower_F <- F_t_i_minus_1 # F_t is non-decreasing
  F_t_i = mean(c(lower_F, upper_F))
  error_B_0 <- error_B

  eps_B <- 1 # error term for |mu_C(B_t) - p_1(t)|
  iter_B <- 0
  ## error tolerance: |mu_C(B_t) - p_1(t)| < error_B
  while (eps_B > error_B) {
    upper_G <- 1
    lower_G <- G_t_i_minus_1 # G is non-decreasing, so G_hat[i] >= G_hat[i-1]
    iter_B <- iter_B + 1
    if (iter_B > iter_max_B) {
      ## after iter = n iterations, the interval [G_lower, G_upper] has length 2^(-n)
      ## further iterations will only slightly improve the "solution"
      ## it may happen that the found "solution" does not satisfy the given precision
      error_B <- error_B + error_B / 10
      if (error_B > error_B_0 * 100) {
        print("No solution found")
        break
      }
      print(paste0("In iteration ", i, " Error B was adjusted to: ", error_B,
                   "Iteration_B count: ", iter_B, " Eps B value is: ", eps_B))
      iter_B <- 0
    }

    ## Step 1: given F(t_i) find corresponding G(t_i), solving mu_C(A_t_i) = est k(t_i)
    eps_A <- 1 # error term for |mu_C(A_t_i) - k(t_i)|
    iter_A <- 0
    ## error tolerance: |mu_C(A_t_i) - k(t_i)| < error_A
    while (eps_A > error_A) {
      mu_C_A_i <- 1 - mean(c(lower_G, upper_G)) - F_t_i + copula(x = F_t_i, y = mean(c(lower_G, upper_G)), theta = theta, tau = tau)
      if (mu_C_A_i > k[i]) {
        # P(X > t_i, Y > t_i) < mu_C_A_i, so prob. to survive was too large / die was too low
        lower_G <- mean(c(lower_G, upper_G)) # increase prob. to die by choosing higher F
      } else if (mu_C_A_i <= k[i]) {
        # P(X > t_i, Y > t_i) > mu_C_A_i, so prob. to survive was too low
        upper_G <- mean(c(lower_G, upper_G)) # increase prob. to survive
      }
      eps_A <- abs(mu_C_A_i - k[i])

      ## it may happen that the estimate for F_t_i is too far off and eps_A is bounded from below
      iter_A <- iter_A + 1
      if (iter_A > iter_max_A) {
        break
      }
    }
    G_t_i <- mean(c(lower_G, upper_G))

    ## Step 2: check if the pair (F(t_i), G(t_i)) satisfies mu_C(B_t) = est p_1(t)
    mu_C_B_i <- mu_C_B_t_i_minus_1 + F_t_i - F_t_i_minus_1 + copula(x = F_t_i_minus_1, y = G_t_i,  theta = theta, tau = tau) - copula(x = F_t_i, y = G_t_i, theta = theta, tau = tau)

    if (mu_C_B_i > p_1[i]) {
      # F(t_i) was too large
      upper_F <- mean(c(lower_F, upper_F)) # increase prob. to survive
    } else if (mu_C_B_i <= p_1[i]) {
      # F(t_i) was too small
      lower_F <- mean(c(lower_F, upper_F)) # increase prob. to die by choosing higher F
    }
    eps_B <- abs(mu_C_B_i - p_1[i])
    F_t_i <- mean(c(lower_F, upper_F))
  }
  return(list(F_t_i, G_t_i, mu_C_B_i, mu_C_A_i))
}
