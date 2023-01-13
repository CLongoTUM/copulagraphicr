#' Visualization of competing risk data
#'
#' @param sample_data Sample data consisting of observed competing risk data (T\eqn{_i}, \eqn{\delta_i}).
#' The dataset should be loaded using the `load_data` function.
#'
#' @export
#'
#' @examples \donttest{## load data set:
#' # data_path <- system.file("inst/Melanoma.csv", package = "copulagraphicr")
#' # sample_data <- copulagraphicr::load_data(data_path = data_path)
#' # copulagraphicr::visualize_data(sample_data = sample_data)
#' }
visualize_data <- function(sample_data) {
  grDevices::dev.new()
  plot(sample_data)
  coord <- graphics::par("usr")
  graphics::points(sample_data[sample_data$delta == 1, ], col = "red")
  graphics::points(sample_data[sample_data$delta == 0, ], col = "blue")
  graphics::legend(x = 0.65 * coord[2], y = 0.8 * coord[4], legend = c("T = X, delta = 1", "T = Y, delta = 0"),
         col = c("red", "blue"), pch = 1, cex = 1, box.lty = 0)
  graphics::title("Competing risk data, observations (T = min(X,Y), delta)")
}
