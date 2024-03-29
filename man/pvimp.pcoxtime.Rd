% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/posthocfuns.R
\name{pvimp.pcoxtime}
\alias{pvimp.pcoxtime}
\title{Permutation variable importance}
\usage{
pvimp.pcoxtime(
  object,
  newdata,
  nrep = 50,
  parallelize = FALSE,
  nclusters = 1,
  estimate = c("mean", "quantile"),
  probs = c(0.025, 0.5, 0.975),
  seed = NULL,
  ...
)
}
\arguments{
\item{object}{fitted \code{\link[pcoxtime]{pcoxtime}}.}

\item{newdata}{data frame containing the variables appearing on the right hand side of \code{\link[pcoxtime]{pcoxtime}} formula.}

\item{nrep}{number of replicates for permutations. Default is \code{nrep = 50}.}

\item{parallelize}{whether to run in parallel. Default is \code{FALSE}.}

\item{nclusters}{number of cores to use if \code{parallelize = TRUE}.}

\item{estimate}{character string specify which summary statistic to use for the estimates. Default is \code{"mean"}.}

\item{probs}{numeric vector of probabilities with values in \code{[0,1]}.}

\item{seed}{a single value for for random number generation.}

\item{...}{for future implementation.}
}
\value{
a named vector of variable scores (\code{estimate = "mean"}) or a data frame (\code{estimate = "quantile"}).
}
\description{
Computes the relative importance based on random permutation of focal variable for pcoxtime model.
}
\details{
Given predictors \code{x_1, x_2, ..., x_n} used to predict the survival outcome, \code{y}. Suppose, for example, \code{x_1} has low predictive power for the response. Then, if we randomly permute the observed values for \code{x_1}, then the prediction for \code{y} will not change much. Conversely, if any of the predictors highly predicts the response, the permutation of that specific predictor will lead to a considerable change in the predictive measure of the model. In this case, we conclude that this predictor is important. In our implementation, Harrel's concordance index is used to measure the prediction accuracy.
}
\keyword{internal}
