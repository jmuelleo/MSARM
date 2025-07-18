% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/Log_Likelihood_function.R
\name{Log_Likelihood_function}
\alias{Log_Likelihood_function}
\title{Log_Likelihood_function}
\usage{
Log_Likelihood_function(alpha, P, Y_T, N, K)
}
\arguments{
\item{alpha}{alpha vector from the previous iteration}

\item{P}{P matrix from the previous iteration}

\item{Y_T}{Data}

\item{N}{Number of underlying regimes}

\item{K}{Lag-order of the AR model}
}
\value{
Returns the value of the log-likelihood function at a specific theta
}
\description{
Log_Likelihood_function allows its user to compute the value of the log-likelihood function based on Hamilton (1994, page 692)
}
\examples{
Log_Likelihood_function(alpha,P,Y_T,N,K)
}
