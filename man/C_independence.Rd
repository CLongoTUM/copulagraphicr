% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/copula_functions.R
\name{C_Independence}
\alias{C_Independence}
\title{Independence Copula}
\usage{
C_Independence(x, y, theta = NA, tau = 0)
}
\arguments{
\item{x}{It is required that \eqn{0 \le x \le 1}}

\item{y}{It is required that \eqn{0 \le y \le 1}}

\item{theta}{Placeholder argument. Not applicable for independence copula.}

\item{tau}{Placeholder argument. Kendall's \eqn{\tau} is zero for the independence copula.}
}
\description{
Independence Copula
}
\examples{
x <- runif(10^2)
y <- runif(10^2)
sample <- C_Independence(x, y)
}
