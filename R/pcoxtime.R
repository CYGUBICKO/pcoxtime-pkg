#' Fit penalized Cox model
#'
#' Fits a Cox model with either lasso, ridge or elasticnet penalty 
#' for both time-independent and time-dependent (varying) covariates survival data.
#'
#' @details
#' The algorithm estimates the coefficients based on observed
#' survival data, with either time-independent or time-dependent
#' covariates, through penalized partial log-likelihood
#'
#' \deqn{\mathrm{pen} ~ \mathit{\ell(\beta)_{\alpha, \lambda}} = -\mathit{\ell(\beta)} + P_{\alpha, \lambda}(\beta)}
#' 
#' using elasticnet (which combines both lasso and ridge) penalty
#'
#' \deqn{\lambda\left(\alpha\sum_{i=1}^p|\beta_i| + 0.5(1 - \alpha)\sum_{i=1}^p\beta_i^2 \right)}.
#'
#' \code{alpha = 1} (\eqn{\alpha}) is the lasso penalty, 
#' and \code{alpha = 0} is the ridge penalty. \code{lambda = 0} fits 
#' the standard Cox proportional hazard model.
#'
#' User can provide a particular lambda. Typical usage is to 
#' use the \code{\link[pcoxtime]{pcoxtimecv}} to select the optimal 
#' lambda first.
#'
#' The routine to handle time-dependent covariates is similar
#' to that implemented in \code{\link[survival]{coxph}}: if there
#' are tied event times, Breslow approximation is used.
#'
#' @param formula object of class formula describing
#' the model. The response is specified similar to
#' \code{\link[survival]{Surv}} function from package 
#' \strong{survival}. The terms (predictors) are specified
#' on the right of "~" in the formula.
#' @param data optional data frame containing
#' variables specified in the formula.
#' @param alpha elasticnet mixing parameter, with 
#' \eqn{0 \le\alpha\le 1}. See details
#' @param lambda tuning parameter for the lasso
#' penalization, with \eqn{\lambda \ge 0}. 
#' \eqn{\lambda = 0} fits unpenalized Cox
#' model. See details
#' @param maxiter maximum number of iterations to
#' convergence. Default is \eqn{1e4}. Consider
#' increasing it if the model does not converge.
#' @param tol convergence threshold for proximal
#' gradient gradient descent. Each proximal update
#' continues until the relative change in all the 
#' coefficients 
#' (i.e. \eqn{\sqrt{\sum(\beta_{k+1} - \beta_k)^2}}/stepsize)
#' is less than tol. The default value is \eqn{1e-8}.
#' @param quietly logical. If TRUE, iteration progress 
#' printed.
#' @param lambmax logical. Sufficiently large, 
#' \eqn{\lambda_{\max}}, that sets \eqn{\beta = 0}
#' for regularization path. If TRUE, \eqn{\lambda_{\max}}
#' is returned. 
#' @param origin_scale logical. If TRUE (default), 
#' the estimated coefficients are returned on the 
#' original covariate scale. Otherwise, FALSE, 
#' coefficients are standardized.
#' @param contrasts.arg an optional list. See
#' the contrasts.arg of
#' \code{\link[stats]{model.matrix.default}}.
#' @param xlevs a named list of character vectors
#' giving the full set of levels to be assumed
#' for each factor. See \code{\link[stats]{model.frame}}.
#' @param na.action a function which indicates
#' what should happen when the data contain NAs.
#' See \code{\link[stats]{model.frame}}.
#' @param ... additional arguments not implemented.
#'
#' @return An S3 object of class \code{\link[pcoxtime]{pcoxtime}}:
#' \item{coef}{a named vector of coefficients. If any of the coefficients violates KKT conditions, the model will print a warning but still return coefficient estimates.}
#' \item{min.nloglik}{estimated log-likelihood at convergence.}
#' \item{min.dev}{the deviation satisfying the \code{tol} stopping creteria.}
#' \item{iter.dev}{deviations between previous and current coefficient estimate at each iteration.}
#' \item{convergence}{convergence message containing the number of iterations}
#' \item{n}{the number of observations used in the fit.}
#' \item{n.risk}{the number of individuals at risk at time \code{t}.}
#' \item{n.event}{the number of events that occur at time \code{t}.}
#' \item{n.censor}{the number of subjects who exit the risk set, without an event, at time \code{t}.}
#' \item{time}{time points defined by the risk set.}
#' \item{Y}{Surv object defining the event times and event status.}
#' \item{data}{data frame used.}
#'	\item{timevarlabel, eventvarlabel}{time and event variables, respectively.}
#' \item{predictors}{a vector of predictors/covariates in the model.}
#' \item{lambda, alpha}{lambda and alpha used, respectively.}
#' \item{formula}{model formula used in the fit.}
#' \item{means}{vector of column means of the X matrix. Subsequent survival curves are adjusted to this value.}
#' \item{}{See \code{\link[stats]{model.frame}} for \code{assign}, \code{xlevels}, \code{contrasts} and \code{terms}.}
#'
#' @seealso
#' \code{\link[survival]{coxph}}, \code{\link[pcoxtime]{pcoxtimecv}}
#'
#' @export
#' @import Rcpp
#' @importFrom stats .getXlevels aggregate approx as.formula coef coefficients delete.response model.extract model.frame model.matrix na.omit na.pass predict setNames terms weighted.mean quantile formula reorder
#' @importFrom grDevices rainbow
#'
#' @examples
#'
#' # Time-independent covariates
#' if (packageVersion("survival")>="3.2.9") {
#' 	data(cancer, package="survival")
#' } else {
#' 	data(veteran, package="survival")
#' }
#' ## Fit unpenalized Cox using pcoxtime
#' lam <- 0 # Should fit unpenalized Cox model
#' pfit1 <- pcoxtime(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, lambda = lam 
#'		, alpha = 1
#'	)
#' print(pfit1)
#'
#' ## fit survival::coxph
#' cfit1 <- coxph(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, method = 'breslow'
#'		, ties = "breslow"
#'	)
#' print(cfit1)
#'
#' ## Penalized Cox model (pcoxtime) 
#' lam <- 0.1
#' alp <- 0.5
#' pfit2 <- pcoxtime(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, lambda = lam 
#'		, alpha = alp
#'	)
#' print(pfit2)
#'
#' # Time-varying covariates
#' data(heart, package="survival")
#' lam <- 0.1
#' alp <- 0.8
#' pfit2 <- pcoxtime(Surv(start, stop, event) ~ age + year + surgery + transplant
#' 	, data = heart
#' 	, lambda = lam
#'		, alpha = alp
#'	)
#' print(pfit2)

#'
#' @docType package
#' @name pcoxtime
#' @useDynLib pcoxtime, .registration=TRUE

pcoxtime <- function(formula, data, alpha = 1, lambda = 1
	, maxiter = 1e5, tol = 1e-8, quietly = FALSE, lambmax = FALSE
	, origin_scale = TRUE, contrasts.arg = NULL, xlevs = NULL
	, na.action = na.omit, ...) {
	
	# Reset alpha if > 1
	if (alpha > 1){alpha_old <- alpha; alpha <- 1} else{alpha_old <- alpha}

	# survival package data format
	if (missing(formula)) stop("a formula argument is required")
   sobj <- if (missing(data)) riskset(formula = formula, na.action = na.action) else riskset(formula = formula, data = data, na.action = na.action)

	sobj <- riskset(formula = formula, data = data, contrasts.arg = contrasts.arg, xlevs = xlevs, na.action = na.action)
	Y <- sobj[["Y"]]
	storage.mode(Y) <- "double"
	X <- sobj[["X"]]
	storage.mode(X) <- "double"
	p <- NCOL(X)
	events <- sobj[["events"]]
	eventvarlabel <- sobj[["eventvarlabel"]]
	times <- sobj[["times"]]
	timevarlabel <- sobj[["timevarlabel"]]
	scale_sd <- sobj[["scale_sd"]]
	scale_mu <- sobj[["scale_mu"]]
	predictors <-  sobj[["term.labels"]]
	assign <- sobj[["assign"]]
	xlevels <- sobj[["xlevels"]]
	contrasts <- sobj[["contrasts"]]
	terms <- sobj[["terms"]]

	if(!is.null(colnames(X))){xnames = colnames(X)}else{xnames = paste0("X", 1:p)}

	# Initialise beta values
	beta0 <- as.double(rep(0, p))
	
	# Proximal iterations and updates
	model_est <- proxiterate(Y = Y
		, X = X
		, beta0 = beta0
		, lambda = as.double(lambda)
		, alpha = as.double(alpha)
		, p = p
		, maxiter = maxiter
		, tol = tol
		, xnames = xnames
		, lambmax = lambmax
	)
	
	if (lambmax) return(model_est)

	if(!quietly) message(model_est[["message"]])

	if (alpha_old > 1) warning("alpha >1; set to default: alpha = 1")

	if (origin_scale){
		beta_hat <- model_est[["beta_hat"]]/scale_sd
	} else {
		beta_hat <- model_est[["beta_hat"]]
	}
	
	## Check KKT conditions
	kkt_violations <- pcoxKKTcheck(model_est[["grad.opt"]]
		, beta_hat
		, lambda
		, alpha
	)

	if (any(kkt_violations)){
		warning(paste0("There are ", sum(kkt_violations), " coefficients which violated KKT conditons. Consider increasing maxiter."))
	}

	beta_hat <- as.matrix(beta_hat, ncol = 1)
	colnames(beta_hat) <- "s0"
	
	# Back transform scaled predictors
	Xmatrix <- unscale(X, scale_mu, scale_sd)
	means <- colMeans(Xmatrix)
	temp_df <- cbind(times, events)
	colnames(temp_df) <- c(timevarlabel, eventvarlabel)
	data <- cbind(temp_df, Xmatrix)

	# Some posthocs
	n.risk <- sort(unique(drop(model_est[["n.risk"]])), decreasing = TRUE)
	temp_df <- data.frame(temp_df[order(temp_df[, timevarlabel]), ])
	temp_df$censorvar <- 1 - temp_df[, eventvarlabel]
	n.event <- as.vector(t(aggregate(as.formula(paste0(eventvarlabel, "~", timevarlabel)), temp_df, FUN = sum)[2]))
	n.censor <- as.vector(t(aggregate(as.formula(paste0("censorvar", "~", timevarlabel)), temp_df, FUN = sum)[2]))
	time <- unique(temp_df[, timevarlabel])
	n <- length(times)
	
	result <- list(coef = beta_hat
		, min.nloglik = model_est[["min.nloglik"]]
		, min.dev = model_est[["min.dev"]]
		, iter.dev = drop(model_est[["deviances"]])
		, convergence = model_est[["message"]]
		, n = n, n.risk = n.risk, n.event = n.event
		, n.censor = n.censor, time = time
		, Y = Y, data = data, timevarlabel = timevarlabel
		, eventvarlabel = eventvarlabel, predictors = predictors
		, lambda = lambda, alpha = alpha
		, formula = formula, means = means, assign = assign
		, xlevels = xlevels, contrasts = contrasts, terms = terms
	)
	result$call <- match.call()
	class(result) <- "pcoxtime"
	result
}



