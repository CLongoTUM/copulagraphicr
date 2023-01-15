#' Estimate empirical quantities `k`, `p_1` and `p_2`
#'
#' In order to calculate the Copula-Graphic estimator the following quantities
#' are needed
#' \deqn{k(t[i]) = pr(X > t[i], Y > t[i]),
#' p_1(t[i]) = pr(X \le t[i], \delta = 1) and p_2(t[i]) = pr(Y \le t[i], \delta = 0).}
#' These quantities can be estimated by their empirical estimates. Given the sample
#' data (T\eqn{_i}, \eqn{\delta_i}), we have
#' \itemize{
#' \item \eqn{est k(t[i]) = 1/n \sum_k 1_(T_k > t[i]) = 1/n \sum_k 1_(X_k > t[i], Y_k > t[i])},
#' \item \eqn{est p_1(t[i]) = 1/n \sum_k 1_(X_k \le t[i], \delta_k = 1) = 1/n \sum_k 1_(X_k \le t[i], X_k \le Y_k)},
#' \item \eqn{est p_2(t[i]) = 1/n \sum_k 1_(X_k \le t[i], \delta_k = 0) = 1/n \sum_k 1_(X_k \le t[i], X_k > Y_k)}.
#'}
#'
#' @param sample_data Sample data consisting of observed competing risk data (T\eqn{_1}, \eqn{\delta_1}), ..., (T\eqn{_n}, \eqn{\delta_n}).
#' The dataset should be loaded using the `load_data` function.
#'
#' @export
#'
#' @examples \donttest{## load data set:
#' # data_path <- system.file("inst/Melanoma.csv", package = "copulagraphicr")
#' # sample_data <- copulagraphicr::load_data(data_path = data_path)
#'
#' ## estimate empirical quantities:
#' # est <- copulagraphicr::estimate_quantities(sample_data = sample_data)
#' # t_grid <- est[[1]]
#' # k <- est[[2]]
#' # p_1 <- est[[3]]
#' # p_2 <- est[[4]]
#' }
estimate_quantities <- function(sample_data) {
  n <- nrow(sample_data)
  t_grid <- c(0, unique(sample_data$T))
  k <- c()
  p_1 <- c()
  p_2 <- c()

  for (i in seq_along(t_grid)) {
    k[i] <- 1 / n * sum(sample_data$T > t_grid[i])
    p_1[i] <-   1 / n * sum(sample_data$T <= t_grid[i] & sample_data$delta == 1)
    p_2[i] <- 1 / n * sum(sample_data$T <= t_grid[i] & sample_data$delta == 0)
  }
  return(list(t_grid, k, p_1, p_2))
}
