#' Compute the risk set indicator
#'
#' Identify set of individuals at risk of 
#' experiencing the event before event time 
#' (failure time), \eqn{t_i}. The risk set, 
#' \eqn{R_i}, is the set of individuals, 
#' \eqn{j}, who had not experienced the 
#' event or had been censored by time \eqn{t_i}. 
#' This function identifies this set and 
#' computes event times satisfying this.
#'
#' @details 
#' Let \eqn{t_1 < t_2 <, ..., t_m}, such
#' that \eqn{m < n} if there are not ties, 
#' otherwise \eqn{m = n}. If covariates are 
#' time-independent the risk set, \eqn{R_i},
#' is the set of individuals who are still at risk
#' at time \eqn{t_i}, i.e., individuals with
#' event/censoring time \eqn{y_j\ge t_i}.
#' For time-dependent covariates risk set at time 
#' \eqn{t_i} is now defined as 
#' \eqn{R(t_i) = \{j : (y^{stop}_{j} \ge t_i) \wedge (y^{start}_{j} < t_i)\}}.
#' The first condition, \eqn{(y^{stop}_{j} \ge t_i)}, 
#' ensures that individual \eqn{j} either experienced 
#' the event or was censored at a later time point than 
#' \eqn{t_i}, while the second condition, 
#' \eqn{(y^{start}_{j} < t_i)}, ensures the start 
#' time was observed before the event.
#'
#' @param formula Object of class formula describing 
#' the model. The response and terms are specified 
#' similar to \code{\link[survival]{Surv}} function.
#' @param data optional data frame containing 
#' variables specified in the formula.
#' @param contrasts.arg an optional list. See 
#' the contrasts.arg of
#' \code{[stats]{model.matrix.default}}.
#' @param xlevs a named list of character vectors 
#' giving the full set of levels to be assumed 
#' for each factor. See \code{[stats]{model.frame}}.
#' @param scaleX logical. If TRUE (default), 
#' predictors are scaled/standardized. This is 
#' used internally.
#' @param na.action a function which indicates 
#' what should happen when the data contain NAs. 
#' See \code{[stats]{model.frame}}.
#'
#' @return A list of survival objects:
#' \item{Y}{Surv object defining the event times and event status.}
#' \item{X}{model matrix of model terms.}
#' \item{events}{observed events.}
#' \item{times}{event times defined by risk set condition.}
#'	\item{timevarlabel, eventvarlabel}{time and event variables, respectively.}
#' \item{scale_sd, scale_mu}{standard deviation and mean of each of the variable used in standardization.}
#' 
#' @keywords internal

riskset <- function(formula, data, contrasts.arg = NULL, xlevs = NULL, scaleX = TRUE, na.action = na.omit){
	
	call <- match.call()
	m <- match.call(expand.dots = FALSE)
	temp <- c("", "formula", "data")
	m <- m[match(temp, names(m), nomatch = 0)]
	if (missing(data)) {
		data <- sys.parent()
	}
	if (missing(formula)) stop("a formula argument is required")
	Terms <- terms(formula, data = data)
	m$formula <- Terms
	m$na.action <- na.action
	m[[1]] <- as.name("model.frame")
	m <- eval(m, sys.parent())
	Terms2 <- terms(m)
	y <- model.extract(m, "response")
	term.labels <- attr(Terms, "term.labels")
	xlevels <- .getXlevels(Terms, m)
	
	if(!inherits(y, "Surv")) stop("formula: must be a survival formula. ?survival")
	ynames <- deparse(formula[[2]])
	N <- NROW(y)
	p <- NCOL(y)
  
	X <- model.matrix(Terms, m, contr = contrasts.arg, xlev = xlevs)
	contrasts <- attr(X, "contrasts")
	xnames <- colnames(X)
	assign <- setNames(attr(X, "assign"), xnames)[-1]
	if (scaleX){
		X <- scale(X[, -1, drop = FALSE])
		scale_mu <- attr(X, "scaled:center")
		scale_sd <- attr(X, "scaled:scale")
	} else {
		X <- X[, -1, drop = FALSE]
		scale_mu <- 0
		scale_sd <- 1
	}
	
	eventvarlabel <- trimws(gsub(".*\\,|\\)", "", ynames))
	# eventvarlabel <- ynames[length(ynames)]
	
	if (p == 2) {
		timevarlabel <- gsub(".*\\(|\\,.*", "", ynames)
		endtime <- y[,1]
		events <- y[,2]
	} else {
		starttime <- y[,1]
		endtime <- y[,2]
		events <- y[,3]
		timevarlabel <- trimws(strsplit(ynames, "\\,", perl = TRUE)[[1]][2])
	}
	
	out <- list(Y = y, X = X, events = events, times = endtime
		, timevarlabel = timevarlabel, eventvarlabel = eventvarlabel
		, scale_sd = scale_sd, scale_mu = scale_mu, assign = assign
		, term.labels = term.labels, xnames = xnames[-1]
		, xlevels = xlevels, contrasts = contrasts, terms = Terms2
	)
	return(invisible(out))
}

#' Extract response from the formula
#'
#' @param model formula.
#'
#' @return A character string. The name of the response variable
#'
#' @keywords internal

getResponse <- function(formula) {
	tt <- terms(formula)
	vars <- as.character(attr(tt, "variables"))[-1]
	response <- attr(tt, "response")
	vars[response]
}

#' Parse formula and return response variable.
#'
#' @param formula Object of class formula describing 
#' the model. 
#' @param data optional data frame containing.
#' variables specified in the formula.
#'
#' @return A character string. The name of the response variable
#'
#' @keywords internal

parseFormula <- function(formula, data, env = parent.frame()) {
	f <- as.formula(formula)
	response <- data.frame(eval(f[[2]], envir = data, enclos = env))
	colnames(response) <- deparse(f[[2]])
	return(response)
}

#' Unscale scaled predictors.
#'
#' @param x a model matrix or data frame with numeric variables.
#' @param center a vector of means of predictors.
#' @param scale a vector of standard deviation of predictors.
#'
#' @return A matrix (model.matrix) with the columns on the original scale.
#'
#' @keywords internal

unscale <- function(x, center, scale) {
	X <- sweep(x, 2, scale, "*")
	X <- sweep(X, 2, center, "+")
	return(X)
}
