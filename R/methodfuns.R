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
#' @export print.pcoxtime
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
#' @export print.pcoxsurvfit
print.pcoxsurvfit <- function(x, ...){
	if (!inherits(x, "pcoxbasehaz")){
		cat("Call:\n")
		print(x$call)
		out <- data.frame(cbind(n = x$n, events = sum(x$events)))
		print(out, row.names = FALSE, ...)
		cat("\n")
	}
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
#' @export print.pcoxtimecv
print.pcoxtimecv <- function(x, ...){
	cat("Call:\n")
	print(x$call)
	cat("\nOptimal parameter values\n")
	out <- data.frame(cbind(lambda.min = x$lambda.min, lambda.1se = x$lambda.1se, alpha.optimal = x$alpha.optimal))
	print(out, row.names = FALSE, ...)
	cat("\n")
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
#' @export coef.pcoxtime
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
#' @export coefficients.pcoxtime
coefficients.pcoxtime <- function(object, ...){
	return(drop(object$coef))
}

#' @export
pcoxsurvfit <- function(fit, newdata, ...) UseMethod("pcoxsurvfit")

#' @export
pcoxbasehaz <- function(fit, centered = TRUE) UseMethod("pcoxbasehaz")

#' @export
concordScore <- function(fit, newdata = NULL, stats = FALSE) UseMethod("concordScore")
