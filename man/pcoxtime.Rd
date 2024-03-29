% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/pcoxtime.R
\docType{package}
\name{pcoxtime}
\alias{pcoxtime}
\title{Fit penalized Cox model}
\usage{
pcoxtime(
  formula,
  data,
  alpha = 1,
  lambda = 1,
  maxiter = 1e+05,
  tol = 1e-08,
  quietly = FALSE,
  lambmax = FALSE,
  origin_scale = TRUE,
  contrasts.arg = NULL,
  xlevs = NULL,
  na.action = na.omit,
  ...
)
}
\arguments{
\item{formula}{object of class formula describing
the model. The response is specified similar to
\code{\link[survival]{Surv}} function from package 
\strong{survival}. The terms (predictors) are specified
on the right of "~" in the formula.}

\item{data}{optional data frame containing
variables specified in the formula.}

\item{alpha}{elasticnet mixing parameter, with 
\eqn{0 \le\alpha\le 1}. See details}

\item{lambda}{tuning parameter for the lasso
penalization, with \eqn{\lambda \ge 0}. 
\eqn{\lambda = 0} fits unpenalized Cox
model. See details}

\item{maxiter}{maximum number of iterations to
convergence. Default is \eqn{1e4}. Consider
increasing it if the model does not converge.}

\item{tol}{convergence threshold for proximal
gradient gradient descent. Each proximal update
continues until the relative change in all the 
coefficients 
(i.e. \eqn{\sqrt{\sum(\beta_{k+1} - \beta_k)^2}}/stepsize)
is less than tol. The default value is \eqn{1e-8}.}

\item{quietly}{logical. If TRUE, iteration progress 
printed.}

\item{lambmax}{logical. Sufficiently large, 
\eqn{\lambda_{\max}}, that sets \eqn{\beta = 0}
for regularization path. If TRUE, \eqn{\lambda_{\max}}
is returned.}

\item{origin_scale}{logical. If TRUE (default), 
the estimated coefficients are returned on the 
original covariate scale. Otherwise, FALSE, 
coefficients are standardized.}

\item{contrasts.arg}{an optional list. See
the contrasts.arg of
\code{\link[stats]{model.matrix.default}}.}

\item{xlevs}{a named list of character vectors
giving the full set of levels to be assumed
for each factor. See \code{\link[stats]{model.frame}}.}

\item{na.action}{a function which indicates
what should happen when the data contain NAs.
See \code{\link[stats]{model.frame}}.}

\item{...}{additional arguments not implemented.}
}
\value{
An S3 object of class \code{\link[pcoxtime]{pcoxtime}}:
\item{coef}{a named vector of coefficients. If any of the coefficients violates KKT conditions, the model will print a warning but still return coefficient estimates.}
\item{min.nloglik}{estimated log-likelihood at convergence.}
\item{min.dev}{the deviation satisfying the \code{tol} stopping creteria.}
\item{iter.dev}{deviations between previous and current coefficient estimate at each iteration.}
\item{convergence}{convergence message containing the number of iterations}
\item{n}{the number of observations used in the fit.}
\item{n.risk}{the number of individuals at risk at time \code{t}.}
\item{n.event}{the number of events that occur at time \code{t}.}
\item{n.censor}{the number of subjects who exit the risk set, without an event, at time \code{t}.}
\item{time}{time points defined by the risk set.}
\item{Y}{Surv object defining the event times and event status.}
\item{data}{data frame used.}
\item{timevarlabel, eventvarlabel}{time and event variables, respectively.}
\item{predictors}{a vector of predictors/covariates in the model.}
\item{lambda, alpha}{lambda and alpha used, respectively.}
\item{formula}{model formula used in the fit.}
\item{means}{vector of column means of the X matrix. Subsequent survival curves are adjusted to this value.}
\item{assign, xlevels, terms}{See \code{\link[stats]{model.frame}} for \code{assign}, \code{xlevels}, \code{contrasts} and \code{terms}.}
}
\description{
Fits a Cox model with either lasso, ridge or elasticnet penalty 
for both time-independent and time-dependent (varying) covariates survival data.
}
\details{
The algorithm estimates the coefficients based on observed
survival data, with either time-independent or time-dependent
covariates, through penalized partial log-likelihood

\deqn{\mathrm{pen} ~ \mathit{\ell(\beta)_{\alpha, \lambda}} = -\mathit{\ell(\beta)} + P_{\alpha, \lambda}(\beta)}

using elasticnet (which combines both lasso and ridge) penalty

\deqn{\lambda\left(\alpha\sum_{i=1}^p|\beta_i| + 0.5(1 - \alpha)\sum_{i=1}^p\beta_i^2 \right)}.

\code{alpha = 1} (\eqn{\alpha}) is the lasso penalty, 
and \code{alpha = 0} is the ridge penalty. \code{lambda = 0} fits 
the standard Cox proportional hazard model.

User can provide a particular lambda. Typical usage is to 
use the \code{\link[pcoxtime]{pcoxtimecv}} to select the optimal 
lambda first.

The routine to handle time-dependent covariates is similar
to that implemented in \code{\link[survival]{coxph}}: if there
are tied event times, Breslow approximation is used.
}
\examples{

# Time-independent covariates
if (packageVersion("survival")>="3.2.9") {
	data(cancer, package="survival")
} else {
	data(veteran, package="survival")
}
## Fit unpenalized Cox using pcoxtime
lam <- 0 # Should fit unpenalized Cox model
pfit1 <- pcoxtime(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
	, data = veteran
	, lambda = lam 
	, alpha = 1
)
print(pfit1)

## fit survival::coxph
cfit1 <- coxph(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
	, data = veteran
	, method = 'breslow'
	, ties = "breslow"
)
print(cfit1)

## Penalized Cox model (pcoxtime) 
lam <- 0.1
alp <- 0.5
pfit2 <- pcoxtime(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
	, data = veteran
	, lambda = lam 
	, alpha = alp
)
print(pfit2)

# Time-varying covariates
data(heart, package="survival")
lam <- 0.1
alp <- 0.8
pfit2 <- pcoxtime(Surv(start, stop, event) ~ age + year + surgery + transplant
	, data = heart
	, lambda = lam
	, alpha = alp
)
print(pfit2)

}
\seealso{
\code{\link[survival]{coxph}}, \code{\link[pcoxtime]{pcoxtimecv}}
}
