% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/helperfuns.R
\name{unscale}
\alias{unscale}
\title{Unscale scaled predictors.}
\usage{
unscale(x, center, scale)
}
\arguments{
\item{x}{a model matrix or data frame with numeric variables.}

\item{center}{a vector of means of predictors.}

\item{scale}{a vector of standard deviation of predictors.}
}
\value{
A matrix (model.matrix) with the columns on the original scale.
}
\description{
Unscale scaled predictors.
}
\keyword{internal}
