#' Simulate survival data for time-dependent covariates
#' 
#' This function uses the permutation algorithm to generate a dataset based on user specified list of covariates, of which some can be time-dependent. User can also specify distribution of event and censoring times.
#'
#' @details 
#' This function is a wrapper to the permutation algorithm implemented in \code{\link[PermAlgo]{permalgorithm}}. The user can fix the pairwise correlation between any 2 predictors (only time-dependent covariates) by specify \code{-1 <= rho < 1}.
#'
#' @param nSubjects number of subjects to simulate, default is \code{100}.
#' @param maxTime a non-zero integer specifying the maximum length of follow-up time, default is \code{365}.
#' @param pfixed number of time-independent (fixed) covariates. They are randomly drawn from a normal distribution with \code{mean = 0} and \code{sd = 1}. The values are replicated for each subject upto the respective follow-up time.
#' @param ptdc number of time-dependent covariates. By default, these are drawn from a normal distribution with \code{mean = 0} and \code{sd = 1} but the user can by-pass this option and specify a matrix of time-dependent covariates via \code{tdcmat}.
#' @param pbin optional. Number of binary covariates. This are treated as fixed covariates and are drawn form a binomial distribution with \code{p = 0.5}.
#' @param betas a vector of 'true' effect sizes (regression coefficients) representing the magnitude of the relationship between the respective covariate and the risk of event. If \code{NULL}, the algorithm generates \code{betas} from a uniform distribution and then converts them to log hazard, i.e., \code{log(runif((pfixed+ptdc+pbin), 0, 2))}. The length of \code{betas} must be the same the total of covariates to generate.
#' @param tdcmat specify own time-dependent covariates. If specified (a matrix with nSubjects*maxTime rows), \code{ptdc} is ignored. This is important in mechanistic simulation of the time-dependent covariates.
#' @param xmat specify an entire matrix for all the covariates. If specified (a matrix with nSubjects*maxTime rows), all the previous specifications for number of covariates and \code{tdcmat} are ignored. This is important in mechanistic simulation of all the covariates or some specific distributional assumptions are required.
#' @param rho specify the pairwise correlation between the time-independent covariates. The default \code{rho = 1} means no pairwise correlation between the covariates.
#' @param eventRandom a non-negative integers of length \code{nSubjects} which represent the subject's event times or a random generating function with \code{n} option specified. If \code{NULL}, the algorithm generates \code{nSubjects} random deviates from exponential distribution with \code{rate = rate}. See \code{rate} option.
#' @param rate the rate for the exponential random deviates for \code{eventRandom}.
#' @param censorRandom a non-negative integers of length \code{nSubjects} which represent the subject's censoring times or a random generating function with \code{n} option specified. If \code{NULL}, the algorithm generates \code{nSubjects} random numbers based on uniform distribution, i.e., \code{runif(nSubjects, 1, maxTime)}.
#' @param groupByD see \code{\link[PermAlgo]{permalgorithm}}.
#' @param x logical. Whether to return matrix of generated covariates in addition to the entire dataset.
#'
#' @seealso \code{\link[PermAlgo]{permalgorithm}}.
#'
#' @return a list of dataset, betas (and matrix of covariates). The covariates have a suffix depending on their type, xbin* for binary, xtf* for time-independent (fixed) and xtd* for time-dependent covariates.
#' \itemize{
#' \item{data}{simulated \code{data.frame} with the following columns}
#'  \itemize{
#'		\item{Id}{subject id. Identifies each of the nSubjects individuals}
#'		\item{Event}{event indicator. \code{Event = 1} if the event occurs otherwise \code{0}.}
#'		\item{Fup}{individual max follow-up time}
#' 	\item{Start}{start of each time interval}
#' 	\item{Stop}{end of each time interval}
#'		\item{x*}{all generated covariates}
#' }
#' \item{betas}{named vector of coefficients specified in the function call. Otherwise, internally generated.}
#' \item{xmat}{if \code{x = TRUE}, matrix of covariates}
#' }
#'
#' @references
#' Sylvestre M.-P., Abrahamowicz M. (2008) \emph{Comparison of algorithms to generate event times conditional on time-dependent covariates}. Statistics in Medicine 27(14):2618--34
#'
#' @examples
#'
#' \dontrun{
#'		library(PermAlgo)
#' 	library(survival)
#' 	library(ggplot2)
#' 	pcoxtheme()
#' 
#' 	set.seed(123407)
#'		# Simulate with default values
#'		df <- simtdc()
#'		head(df$data)
#' 	# Simulate for a number of times to check stability of the estimates
#' 	nrep <- 500
#' 	betas <- log(runif(6, 0, 2))
#' 	beta_list <- list()
#' 	true_list <- list()
#' 	for (i in 1:nrep){
#' 		sim <- simtdc(pfixed = 3, ptdc = 2, pbin = 1, betas = betas)
#' 		df <- sim$data
#' 		vnames <- colnames(df)[!colnames(df) %in% c("Id", "Fup")]
#' 		df <- df[ ,vnames]
#' 		# Estimate coefficients using coxph
#' 		mod <- coxph(Surv(Start, Stop, Event) ~ ., df)
#' 		beta_list[[i]] <- coef(mod)
#' 		true_list[[i]] <- sim$betas
#' 	}
#' 	beta_df <- data.frame(do.call("rbind", beta_list))
#' 	beta_df <- stack(beta_df)
#' 
#' 	true_df <- data.frame(ind = names(true_list[[1]]), values = true_list[[1]])
#' 	p1 <- (ggplot(beta_df, aes(x = values))
#' 		+ geom_histogram(alpha = 0.3)
#' 		+ geom_vline(data = true_df, aes(xintercept = values), col = "blue")
#' 		+ facet_wrap(~ind, scales = "free")
#' 		+ labs(x = "Beta estimate", y = "")
#' 	)
#' 	print(p1)
#' }
#'
#' @export
#' @importFrom stats rbinom rexp rnorm runif

simtdc <- function(nSubjects = 100, maxTime = 365, pfixed = 2
	, ptdc = 2, pbin = NULL, betas = NULL, tdcmat = NULL
	, xmat = NULL, rho = 1, eventRandom = NULL, rate = 0.012 
	, censorRandom = NULL, groupByD = FALSE, x = FALSE) {
	if (is.null(xmat)){
		xmat <- genX(nSubjects = nSubjects, maxTime = maxTime
			, pfixed = pfixed, ptdc = ptdc, pbin = pbin
			, rho = rho, tdcmat = tdcmat
		)
	} 
	if (!is.null(betas) & (length(betas)!=NCOL(xmat))) stop("Number of coefficients must be the same as number of covariates")
	if (is.null(betas)) betas <- log(runif(NCOL(xmat), 0, 2))
	if (is.null(eventRandom)) eventRandom <- round(rexp(nSubjects, rate)+1,0)
	if (is.null(censorRandom)) censorRandom <- round(runif(nSubjects, 1, maxTime),0)
	
	df <- PermAlgo::permalgorithm(numSubjects = nSubjects
		, maxTime = maxTime, Xmat = xmat, XmatNames = colnames(xmat)
		, eventRandom = eventRandom, censorRandom = censorRandom
		, betas = betas, groupByD = groupByD
	)
	names(betas) <- colnames(xmat)
	if (x){
		out <- list(data = df, betas = betas, xmat = xmat)
	} else {
		out <- list(data = df, betas = betas)
	}
	return(out)
}

#' Generate covariates for which some are time-dependent and pairwise correlated.
#' @rdname simtdc
#' @keywords internal

genX <- function(nSubjects = 100, maxTime = 365, pfixed = 2, ptdc = 2
	, pbin = NULL, tdcmat = NULL, rho = 1){
	if(rho < -1 | rho > 1) stop("choose rho between -1 and 1")
	if (any(c(pfixed, ptdc, pbin) < 1)) stop("Number of covariates should be at least 1.")
	xfixed <- replicate(pfixed, rep(rnorm(nSubjects), each = maxTime))
	if(abs(rho)<1){
		b <- sqrt(rho/(1-rho))
		z <- rep(rnorm(nSubjects), each = maxTime)
		x0 <- b*matrix(z,nrow=nSubjects*maxTime,ncol=pfixed,byrow=FALSE) + xfixed
	}
	if (abs(rho)==1) {x0 = xfixed}
	colnames(x0) <- paste0("xtf", 1:pfixed)
	xfixed <- x0
	if(!is.null(pbin)){
		xbin <- replicate(pbin, rep(rbinom(nSubjects, 1, 0.5), each = maxTime))
		colnames(xbin) <- paste0("xbin", 1:pbin)
		xfixed <- cbind(xbin, xfixed)
	}
	if (is.null(tdcmat)){
		xtdc <- replicate(ptdc, do.call("c", lapply(1:nSubjects, function(i)rnorm(maxTime))))
	} else {
		xtdc <- tdcmat
	}
	colnames(xtdc) <- paste0("xtd", 1:NCOL(xtdc))
	x <- cbind(xfixed, xtdc)
	return(x)
}

