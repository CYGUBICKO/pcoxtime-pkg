#' Print coefficients from a pcoxtime object
#'
#' This function prints a summary of the pcoxtime object.
#'
#' @details The call that produced \code{\link[pcoxtime]{pcoxtime}} is printed, followed by coefficient estimates with their corresponding exponentiated values.
#' Depending on the number of coefficients, \code{nprint} can be used to specify the number of coefficients to print out.
#'
#' @param x fitted \code{\link[pcoxtime]{pcoxtime}} model object
#' @param ... for future implementations
#' @param nprint number of coefficients to print out
#'
#' @return A two column output, the first column is the coefficient estimate and the second column is the exponent of the coefficient estimate. Additonal summary about the number of nonzero coefficents, the number of observations and the number of event of interest are also printed.
#'
#' @method print pcoxtime
#' @export
print.pcoxtime <- function(x, ..., nprint = 10){
	cat("Call:\n")
	print(x$call)
	cat("\n")
	if (length(coef(x))<=nprint){
		print(cbind(coef = coef(x), "exp(coef)" = exp(coef(x))))
	}
	cat(sum(coef(x)!=0), "out of", length(coef(x)), "coefficients are nonzero")
	cat("\n")
	cat("n =", x$n,", number of events =", sum(x$n.event), "\n")
}

#' Print a short summary of survival function
#'
#' Print the number of observations and number of events.
#'
#' @details Provide a summary of \code{\link[pcoxtime]{pcoxsurvfit.pcoxtime}} object.
#'
#' @param x the result of a call to the \code{\link[pcoxtime]{pcoxsurvfit.pcoxtime}} function.
#' @param ... for future implementations
#'
#' @return The call to the \code{\link[pcoxtime]{pcoxsurvfit.pcoxtime}} and the summary of the survival function.
#'
#' @method print pcoxsurvfit
#' @export
print.pcoxsurvfit <- function(x, ...){
	cat("Call:\n")
	print(x$call)
	out <- data.frame(cbind(n = x$n, events = sum(x$events)))
	print(out, row.names = FALSE, ...)
	cat("\n")
}

#' Print baseline hazard function data frame
#'
#' Print the head of baseline hazard function data frame.
#'
#' @details Provide a summary of \code{\link[pcoxtime]{pcoxbasehaz.pcoxtime}} object.
#'
#' @param x the result of a call to the \code{\link[pcoxtime]{pcoxbasehaz.pcoxtime}} function.
#' @param n number of rows to print. Default is 5.
#' @param ... for future implementations
#'
#' @return The call to the \code{\link[pcoxtime]{pcoxbasehaz.pcoxtime}} and the head of baseline hazard function data frame.
#'
#' @method print pcoxbasehaz
#' @export
print.pcoxbasehaz <- function(x, n=5, ...){
	cat("Call:\n")
	print(x$call)
	out <- data.frame(time=x$time, hazard=x$hazard, surv=x$surv)
	print(head(out, n=n), row.names = FALSE, ...)
	cat("\n")
}

#' Print cross-validated pcoxtime object
#'
#' Print the summary of the result of cross-validation for a pcoxtime object.
#'
#' @details
#' A summary of optimal lambda and alpha for training pcoxtime model.
#'
#' @param x \code{\link[pcoxtime]{pcoxtimecv}} object
#' @param ... for future implementations
#'
#' @return The call to the \code{\link[pcoxtime]{pcoxtimecv}} and the summary of the optimal alpha and lambdas.
#'
#' @method print pcoxtimecv
#' @export
print.pcoxtimecv <- function(x, ...){
	cat("Call:\n")
	print(x$call)
	cat("\nOptimal parameter values\n")
	out <- data.frame(cbind(lambda.min = x$lambda.min, lambda.1se = x$lambda.1se, alpha.optimal = x$alpha.optimal))
	print(out, row.names = FALSE, ...)
	cat("\n")
}

#' Extract coefficient estimates of pcoxtimecv object
#' 
#' This function extracts cross-validation estimates for a particular lambda.
#'
#' @details Extract the coefficient estimates for optimal lambda-alpha pair or based on specified the value of lambda for an optimal alpha. If the value of lambda specified is not exact (not in lambdas), the nearest value is used, based on \code{nearest <- function(values, value){values[which(abs(values-value)==min(abs(values-value)))]}}. It requires that \code{\link[pcoxtime]{pcoxtimecv}} is run with \code{refit = TRUE}.
#'
#' @param object \code{\link[pcoxtime]{pcoxtimecv}} object
#' @param lambda the value of lambda for which to return the coefficient estimates. It can be any of the character string, "min", "optimal" or "best" for optimal lambda; "1se" for 1 standard error lambda; or any numeric value for lambda. See details.
#' @param ... for future implementations
#'
#' @return A data frame of coefficient estimates.
#'
#' @method coef pcoxtimecv
#' @export
#' @importFrom utils head
coef.pcoxtimecv <- function(object, lambda, ...){
	betas <- object$fit$beta
	if (is.null(betas))stop("Run pcoxtimecv with refit = TRUE to extract coefficients")
	nearest <- function(values, value){
		values[which(abs(values-value)==min(abs(values-value)))]
	}
	if (missing(lambda)) lambda <- "min"
	if (lambda=="min" | lambda=="optimal" | lambda=="best") {
		lambda <- object$lambda.min
	} else if (lambda=="1se") {
		lambda <- object$lambda.1se
	}
	values <- unique(betas$lambda)
	lambda <- nearest(values, lambda)
	betas <- betas[betas$lambda==lambda, ]
	return(betas)
}

#' Extract coefficient estimates of pcoxtimecv object
#' 
#' @return A vector of coefficient estimates.
#'
#' @method coefficients pcoxtimecv
#' @rdname coef.pcoxtimecv
#' @export 
coefficients.pcoxtimecv <- function(object, lambda, ...){
	return(coef.pcoxtimecv(object, lambda, ...))
}


#' Extract coefficient estimates of pcoxtime object
#'
#' This function extracts the estimates for all the coefficients.
#'
#' @details The call that produced \code{\link[pcoxtime]{pcoxtime}} is printed, followed by coefficient estimates.
#'
#' @param object fitted \code{\link[pcoxtime]{pcoxtime}} model object
#' @param ... for future implementations
#'
#' @return A vector of coefficient estimates.
#'
#' @method coef pcoxtime
#' @export
coef.pcoxtime <- function(object, ...){
	return(drop(object$coef))
}

#' Extract coefficient estimates of pcoxtime object
#'
#' @return A vector of coefficient estimates.
#'
#' @method coefficients pcoxtime
#' @rdname coef.pcoxtime
#' @export
coefficients.pcoxtime <- function(object, ...){
	return(drop(object$coef))
}

#' @export
pcoxsurvfit <- function(fit, newdata, ...) UseMethod("pcoxsurvfit")

#' @export
pcoxbasehaz <- function(fit, centered = TRUE) UseMethod("pcoxbasehaz")

#' @export
concordScore <- function(fit, newdata = NULL, stats = FALSE, reverse = TRUE, ...) UseMethod("concordScore")

#' @export
extractoptimal <- function(object, what=c("optimal", "cvm", "coefs"), ...) UseMethod("extractoptimal") 

#' @export
varimp <- function(object, newdata, type=c("coef", "perm", "model")
	, relative=TRUE, nrep=50, parallelize=FALSE, nclusters=1
	, estimate=c("mean", "quantile"), probs=c(0.025, 0.5, 0.975)
	, seed=NULL, ...) UseMethod("varimp")

