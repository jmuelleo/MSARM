% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/eta_function.R
\name{eta_function}
\alias{eta_function}
\title{eta_function}
\usage{
eta_function(alpha, Y_T, N, K)
}
\arguments{
\item{alpha}{alpha vector (intercept, coefficients and error variance) from the last iteration}

\item{Y_T}{Data}

\item{N}{Number of underlying regimes}

\item{K}{Lag-order of the AR model}
}
\value{
returns a matrix of the dimensions (T x N) containing the values of the conditional density function
}
\description{
The eta-vector, defined as in Hamilton (1994, page 691) is computed
}
\examples{
eta_function(alpha,Y_T,N,K)
}
