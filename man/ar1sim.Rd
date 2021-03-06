% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/demos.R
\name{ar1sim}
\alias{ar1sim}
\title{ar1sim - a demonstration that simulating data using arima.sim and from
the analogous error covariance matrix yields assymptotically equivalent
results}
\usage{
ar1sim(B = 20, n = 20, phi = 0.5, sigma2 = 0.75)
}
\description{
ar1sim - a demonstration that simulating data using arima.sim and from
the analogous error covariance matrix yields assymptotically equivalent
results
}
\examples{
Out <- ar1sim(20000)
}
\author{
Stephen Tueller \email{stueller@rti.org}
}
