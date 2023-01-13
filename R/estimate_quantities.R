#' Estimate empirical quantities `k`, `p_1` and `p_2`
#'
#' In order to calculate the Copula-Graphic estimator the following quantities
#' are needed
#' \deqn{k(t) = pr(X>t, Y>t),
#' p_1(t) = pr(X \le t, \delta = 1) and p_2(t) = pr(Y \le t, \delta = 0).}
#' These quantities can be estimated by their empirical estimates
#' \itemize{
#' \item \eqn{est k(t_i) = 1/n \sum 1_(X_j > t_i, Y_j > t_i)},
#' \item \eqn{est p_1(t_i) = 1/n \sum 1_(X_j \le t_i, X_j \le Y_j)},
#' \item \eqn{est p_2(t_i) = 1/n \sum 1_(X_j \le t_i, X_j > Y_j)}.
#'}
#'
#' @param sample_data Sample data consisting of observed competing risk data (T\eqn{_i}, \eqn{\delta_i}).
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
