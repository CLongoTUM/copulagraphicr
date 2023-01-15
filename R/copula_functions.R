#' Independence Copula
#'
#' @param x It is required that \eqn{0 \le x \le 1}
#' @param y It is required that \eqn{0 \le y \le 1}
#' @param theta Placeholder argument. Not applicable for independence copula.
#' @param tau Placeholder argument. Kendall's \eqn{\tau} is zero for the independence copula.
#'
#' @export
#'
#' @examples x <- runif(10^2)
#' y <- runif(10^2)
#' sample <- C_Independence(x, y)
C_Independence <- function(x, y, theta = NA, tau = 0) {
  if (sum(x[x < 0 | x > 1]) + sum(y[y < 0 | y > 1]) == 0) {
    C <- x*y
    return(C)
  } else {
    print("Input in [0,1] x [0,1] is required.")
  }
}


#' Gamma frailty Copula
#'
#' @param x It is required that \eqn{0 \le x \le 1}
#' @param y It is required that \eqn{0 \le y \le 1}
#' @param theta Parameter of the Gamma Frailty Copula.
#' It is required that `theta` \eqn{\ge 1}.
#' @param tau Placeholder argument. It is required that each copula has the input
#' argument `tau` (Kendall's \eqn{\tau}. Unfortunately the closed form solution for
#' Kendall's \eqn{\tau} is unknown to the author.
#'
#' @export
#'
#' @examples x <- runif(10^2)
#' y <- runif(10^2)
#' sample <- C_Gamma_Frailty(x, y)
C_Gamma_Frailty <- function(x, y, theta = 3, tau = NA) {
  if (sum(x[x < 0 | x > 1]) + sum(y[y < 0 | y > 1]) == 0) {
    if (is.na(theta)) {
      theta <- 3
    }
    if (theta >= 1) {
      C <- x+y-1+((1/(1-x))^(theta-1)+(1/(1-y))^(theta-1)-1)^(-1/(theta-1))
      return(C) } else {
        print("Error. Parameter theta >= 1 is required")
      }
  } else {
    print("Input in [0,1] x [0,1] is required.")
  }
}



#' Gumbel Copula
#'
#' This function evaluates the Gumbel Copula distribution at (`x`, `y`). One can
#' define the parameter `theta` of the Gumbel Copula. Alternatively one can
#' set the target value for Kendall's \eqn{\tau} and the corresponding value of `theta`
#' will be set automatically.
#'
#'
#' @param x It is required that \eqn{0 \le x \le 1}
#' @param y It is required that \eqn{0 \le y \le 1}
#' @param theta Parameter of the Gumbel Copula.
#' It is required that `theta` \eqn{\ge 1}. If set to `NA`, `theta` is chosen such that
#' Kendall's \eqn{\tau = 0.5}. Note that the parameter `theta` is only used, if `NA` value for
#' `tau` is passed.
#' @param tau Kendall's \eqn{\tau} for the Gumbel Copula. If set to `NA`, the
#' value of Kendall's \eqn{\tau} is defined by `theta`.
#' \itemize{
#' \item `NA`: If no value for `tau` is passed to the function,
#' the parameter `theta` will be used.
#' \item `tau`: If a value for `tau` is assigned, the parameter `theta`
#' will be chosen such that the passed value for the Kendall's \eqn{\tau} is achieved.
#' }
#'
#' @export
#'
#' @examples x <- runif(10^2)
#' y <- runif(10^2)
#' sample <- C_Gumbel(x, y)
C_Gumbel <- function(x, y, theta = VineCopula::BiCopTau2Par(family = 4, tau = 0.5, check.taus = TRUE), tau = NA) {
  if (is.na(theta)) {
    theta <- VineCopula::BiCopTau2Par(family = 4, tau = 0.5, check.taus = TRUE)
  }
  if (!is.na(tau)) {
    theta <- VineCopula::BiCopTau2Par(family = 4, tau = tau, check.taus = TRUE)
  }
  if (sum(x[x < 0 | x > 1]) + sum(y[y < 0 | y > 1]) == 0) {
    if (theta >= 1) {
    C <- exp(-( (-log(x))^theta + (-log(y))^(theta))^(1/theta))
      return(C) } else {
      print("Error. Parameter theta >= 1 is required")
    }
  } else {
    print("Input in [0,1] x [0,1] is required.")
  }
}



# BiCopPar2Tau(family = 5, par = 5.747565, par2 = 0, obj = NULL, check.pars = TRUE)
# setTheta(copFrank, iTau(copFrank, 0.5)) -> 5.736283
#' Frank's Copula
#'
#' This function evaluates the Frank Copula distribution at (`x`, `y`). One can
#' define the parameter `theta` of the Frank Copula. Alternatively one can
#' set the target value for Kendall's \eqn{\tau} and the corresponding value of `theta`
#' will be set automatically.
#'
#' @param x It is required that \eqn{0 \le x \le 1}
#' @param y It is required that \eqn{0 \le y \le 1}
#' @param theta Parameter of the Frank Copula.
#' It is required that `theta` \eqn{=/= 0}. If set to `NA`, `theta` is chosen such that
#' Kendall's \eqn{\tau = 0.5}. Note that the parameter `theta` is only used, if `NA` value for
#' `tau` is passed.
#' @param tau Kendall's \eqn{\tau} for the Frank Copula. If set to `NA`, the
#' value of Kendall's \eqn{\tau} is defined by `theta`.
#' \itemize{
#' \item `NA`: If no value for `tau` is passed to the function,
#' the parameter `theta` will be used.
#' \item `tau`: If a value for `tau` is assigned, the parameter `theta`
#' will be chosen such that the passed value for the Kendall's \eqn{\tau} is achieved.
#' }
#'
#' @export
#'
#' @examples x <- runif(10^2)
#' y <- runif(10^2)
#' sample <- C_Frank(x, y)
C_Frank <- function(x, y, theta = VineCopula::BiCopTau2Par(family = 5, tau = 0.5, check.taus = TRUE), tau = NA) {
  if (is.na(theta)) {
    theta <- VineCopula::BiCopTau2Par(family = 5, tau = 0.5, check.taus = TRUE)
  }
  if (!is.na(tau)) {
    theta <- VineCopula::BiCopTau2Par(family = 5, tau = tau, check.taus = TRUE)
  }
  if (sum(x[x < 0 | x > 1]) + sum(y[y < 0 | y > 1]) == 0) {
    if (theta != 0) {
      C <- -1/theta*log(1+(exp(-theta*x)-1)*(exp(-theta*y)-1)/(exp(-theta)-1))
      return(C) } else {
      print("Error. Parameter theta =/= 0 is required")
    }
  } else {
    print("Input in [0,1] x [0,1] is required.")
  }
}



# setTheta(copClayton, iTau(copClayton, 0.5)) -> theta = 2
# BiCopPar2Tau(family = 3, par = 2, par2 = 0, obj = NULL, check.pars = TRUE)
#' Clayton copula
#'
#' This function evaluates the Clayton Copula distribution at (`x`, `y`). One can
#' define the parameter `theta` of the Clayton Copula. Alternatively one can
#' set the target value for Kendall's \eqn{\tau} and the corresponding value of `theta`
#' will be set automatically.
#'
#' @param x It is required that \eqn{0 \le x \le 1}
#' @param y It is required that \eqn{0 \le y \le 1}
#' @param theta Parameter of the Clayton Copula.
#' It is required that `theta` \eqn{\ge -1}. If set to `NA`, `theta` is chosen such that
#' Kendall's \eqn{\tau = 0.5}. Note that the parameter `theta` is only used, if `NA` value for
#' `tau` is passed.
#' @param tau Kendall's \eqn{\tau} for the Clayton Copula. If set to `NA`, the
#' value of Kendall's \eqn{\tau} is defined by `theta`.
#' \itemize{
#' \item `NA`: If no value for `tau` is passed to the function,
#' the parameter `theta` will be used.
#' \item `tau`: If a value for `tau` is assigned, the parameter `theta`
#' will be chosen such that the passed value for the Kendall's \eqn{\tau} is achieved.
#' }
#'
#' @export
#'
#' @examples x <- runif(10^2)
#' y <- runif(10^2)
#' sample <- C_Clayton(x, y)
C_Clayton <- function(x, y, theta = VineCopula::BiCopTau2Par(family = 3, tau = 0.5, check.taus = TRUE), tau = NA) {
  if (is.na(theta)) {
    theta <- VineCopula::BiCopTau2Par(family = 3, tau = 0.5, check.taus = TRUE)
  }
  if (!is.na(tau)) {
    theta <- VineCopula::BiCopTau2Par(family = 3, tau = tau, check.taus = TRUE)
  }
  if (sum(x[x < 0 | x > 1]) + sum(y[y < 0 | y > 1]) == 0) {
    if (theta >= -1) {
      C <- (max(x^(-theta)+ y^(-theta),0))^(-1/theta)
      return(C) } else {
      print("Error. Parameter theta >= -1 is required")
    }
  } else {
    print("Input in [0,1] x [0,1] is required.")
  }
}


# BiCopPar2Tau(family = 6, par = 2.85625, par2 = 0, obj = NULL, check.pars = TRUE)
#' Joe copula
#'
#' This function evaluates the Joe Copula distribution at (`x`, `y`). One can
#' define the parameter `theta` of the Joe Copula. Alternatively one can
#' set the target value for Kendall's \eqn{\tau} and the corresponding value of `theta`
#' will be set automatically.
#'
#' @param x It is required that \eqn{0 \le x \le 1}
#' @param y It is required that \eqn{0 \le y \le 1}
#' @param theta Parameter of the Joe Copula.
#' It is required that `theta` \eqn{\ge 1}. If set to `NA`, `theta` is chosen such that
#' Kendall's \eqn{\tau = 0.5}. Note that the parameter `theta` is only used, if `NA` value for
#' `tau` is passed.
#' @param tau Kendall's \eqn{\tau} for the Joe Copula. If set to `NA`, the
#' value of Kendall's \eqn{\tau} is defined by `theta`.
#' \itemize{
#' \item `NA`: If no value for `tau` is passed to the function,
#' the parameter `theta` will be used.
#' \item `tau`: If a value for `tau` is assigned, the parameter `theta`
#' will be chosen such that the passed value for the Kendall's \eqn{\tau} is achieved.
#' }
#'
#' @export
#'
#' @examples x <- runif(10^2)
#' y <- runif(10^2)
#' sample <- C_Joe(x, y)
C_Joe <- function(x, y, theta = VineCopula::BiCopTau2Par(family = 6, tau = 0.5, check.taus = TRUE), tau = NA) {
  if (is.na(theta)) {
    theta <- VineCopula::BiCopTau2Par(family = 6, tau = 0.5, check.taus = TRUE)
  }
  if (!is.na(tau)) {
    theta <- VineCopula::BiCopTau2Par(family = 6, tau = tau, check.taus = TRUE)
  }
  if (sum(x[x < 0 | x > 1]) + sum(y[y < 0 | y > 1]) == 0) {
    if (theta >= 1) {
    C <- 1-((1-x)^theta+(1-y)^theta-(1-x)^theta*(1-y)^theta)^(1/theta)
    return(C) } else {
      print("Error. Parameter theta >= 1 is required")
    }
  } else {
    print("Input in [0,1] x [0,1] is required.")
  }
}
