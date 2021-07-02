#' Cross-validation for pcoxtime
#'
#' Performs \code{k}-fold cross-validation for pcoxtime, plots 
#' solution path plots, and returns optimal value of lambda
#' (and optimal alpha if more than one is given).
#'
#' @details
#' The function fits \code{\link[pcoxtime]{pcoxtime}} \code{folds + 1} (if \code{refit = FALSE}) or \code{folds + 2} times  (if \code{refit = FALSE}). In the former case, the solution path to display using \code{\link[pcoxtime]{plot.pcoxtimecv}} is randomly picked from all the cross-validation runs. However, in the later case, the solution path plot is based on the model refitted using the optimal parameters. In both cases, the function first runs  \code{\link[pcoxtime]{plot.pcoxtimecv}} to compute the lambda sequence and then perform cross-validation on \code{nfolds}.
#'
#' If more than one \code{alphas} is specified, say code{(0.2, 0.5, 1)}, the \code{pcoxtimecv} will search (experimental) for optimal values for alpha with respect to the corresponding lambda values. In this case, optimal alpha and lambda sequence will be returned, i.e., the \code{(alphas, lambdas)} pair that corresponds to the lowest predicted cross-validated error (likelihood deviance).
#'
#' For data sets with a very large number of predictors, it is recommended to only calculate partial paths by lowering the value of \code{lamfract}. In other words, for \code{p > n} problems, the near \code{lambda = 0} solution is poorly behaved and this may account for over \code{99\%} of the function's runtime. We therefore recommend always specifying \code{lamfract < 1} and increase if the optimal lambda suggests lower values. 
#'
#'
#' @param formula object of class formula describing
#' the model. The response is specified similar to
#' \code{\link[survival]{Surv}} function from package
#' \strong{survival}. The terms (predictors) are specified
#' on the right of "~" in the formula.
#' @param data optional data frame containing
#' variables specified in the formula.
#' @param alphas elasticnet mixing parameter, with
#' \code{0 <= alphas <= 1}. If a vector of  \code{alphas} is supplied, cross-validation will be performed for each of the \code{alphas} and optimal value returned. The default is \code{1}.
#' @param lambdas optional user-supplied sequence. If \code{lambdas = NULL} (default -- highly recommended), the algorithm chooses its own sequence.
#' @param nlambdas the default number of lambdas values. Default is \code{100}.
#'	@param lammin_fract smallest value of \code{lambda}, as fraction of maximum \code{lambda}. If \code{NULL}, default, it depends on the number of observations (n) relative to the number of variables (p). If \code{n > p}, the default is \code{0.0001}, otherwise \code{0.01}. Increasing this value may lead to faster convergence. 
#' @param lamfract proportion of regularization path to consider. If \code{lamfract = 1}, complete regularization path is considered. However, if \code{0.5 <= lamfract <1}, only a proportion of the \code{nlambdas} considered. Choosing a smaller \code{lamfract} reduces computational time and potentially stable estimates for model with large number of predictors. See details.
#' @param nfolds number of folds. Default is \code{10}. The smallest allowable is \code{nfolds = 3}.
#' @param foldids an optional sequence of values between \code{1} and {nfolds} specifying what fold each observation is in. This is important when comparing performance across models. If specified, \code{nfolds} can be missing.
#' @param devtype loss to use for cross-validation. Currently, two options are available but versions will implement \code{\link[pcoxtime]{concordScore.pcoxtime}} loss too. The two are, default \code{(devtype = "vv")} Verweij Van Houwelingen partial-likelihood deviance and basic cross-validated parial likelihood \code{devtype = "basic"}. See Dai, B., and Breheny, P. (2019) for details.
#' @param refit logical. Whether to return solution path based on optimal lambda and alpha picked by the model. Default is \code{refit = FALSE}.
#' @param maxiter maximum number of iterations to convergence. Default is \eqn{1e5}. Consider increasing it if the model does not converge.
#' @param tol convergence threshold for proximal gradient gradient descent. Each proximal update continues until the relative change in all the coefficients (i.e. \eqn{\sqrt{\sum(\beta_{k+1} - \beta_k)^2}}/stepsize) is less than tol. The default value is \eqn{1e-8}.
#' @param quietly logical. If TRUE, refit progress is printed.
#' @param seed random seed. Default is \code{NULL}, which generated the seed internally.
#' @param nclusters number of cores to use to run the cross-validation in parallel. Default is \code{nclusters = 1} which runs serial.
#' @param na.action a function which indicates what should happen when the data contain NAs.
#' @param ... additional arguments not implemented.
#'
#' @return An S3 object of class \code{\link[pcoxtime]{pcoxtimecv}}:
#' \item{lambda.min}{the value of lambda that gives minimum cross-validated error.}
#' \item{lambda.1se}{largest value of lambda such that error is within \code{1} standard error of the minimum.}
#' \item{alpha.optimal}{optimal alpha corresponding to \code{lambda.min}.}
#' \item{lambdas.optimal}{the sequence of lambdas containing \code{lambda.min}.}
#' \item{foldids}{the fold assignment used.}
#' \item{dfs}{list of data frames containing mean cross-validated error summaries and estimated coefficients in each fold.}
#' \item{fit}{if \code{refit = TRUE}, summaries corresponding to the optimal \code{alpha} and \code{lambdas}. This is used to plot solution path}.
#'
#' @seealso
#'\code{\link[pcoxtime]{plot.pcoxtimecv}}, \code{\link[pcoxtime]{pcoxtime}}
#'
#' @references 
#' Dai, B., and Breheny, P. (2019). \emph{Cross validation approaches for penalized Cox regression}. \emph{arXiv preprint arXiv:1905.10432}.
#'
#'Simon, N., Friedman, J., Hastie, T., Tibshirani, R. (2011) \emph{Regularization Paths for Cox's Proportional Hazards Model via Coordinate Descent, Journal of Statistical Software, Vol. 39(5) 1-13} \url{https://www.jstatsoft.org/v39/i05/}
#'
#' @examples
#'
#' # Time-independent covariates
#' if (packageVersion("survival")>="3.2.9") {
#'    data(cancer, package="survival")
#' } else {
#'    data(veteran, package="survival")
#' }
#' \donttest{
#' cv1 <- pcoxtimecv(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, alphas = 1
#'		, refit = FALSE
#'		, lamfract = 0.6
#'	)
#' print(cv1)
#'
#' # Train model using optimal alpha and lambda
#' fit1 <- pcoxtime(Surv(time, status) ~ factor(trt) + karno + diagtime + age + prior
#'		, data = veteran
#'		, alpha = cv1$alpha.optimal
#'		, lambda = cv1$lambda.min
#'	)
#' print(fit1)
#' # Time-varying covariates
#' data(heart, package="survival")
#' cv2 <- pcoxtimecv(Surv(start, stop, event) ~ age + year + surgery + transplant
#' 	, data = heart
#'		, alphas = 1
#'		, refit = FALSE
#'		, lamfract = 0.6
#'	)
#' print(cv2)
#'
#' # Train model
#' fit2 <- pcoxtime(Surv(start, stop, event) ~ age + year + surgery + transplant
#' 	, data = heart
#'		, alpha = cv2$alpha.optimal
#' 	, lambda = cv2$lambda.min
#'	)
#' print(fit2)
#' }
#' @export
#' @import parallel
#' @import foreach
#' @import doParallel 

pcoxtimecv <- function(formula, data, alphas = 1, lambdas = NULL
	, nlambdas = 100, lammin_fract = NULL, lamfract = 0.6, nfolds = 10
	, foldids = NULL, devtype = "vv", refit = FALSE, maxiter = 1e5
	, tol = 1e-8, quietly = FALSE, seed = NULL, nclusters = 1
	, na.action = na.omit, ...) {
	
	if(is.null(seed)){seed = 1254}
	set.seed(seed)
	
	# survival package data format
	if (missing(formula)) stop("a formula argument is required")
	sobj <- if (missing(data)) riskset(formula = formula, na.action = na.action) else riskset(formula = formula, data = data, na.action = na.action) 
	Y <- sobj[["Y"]]
	storage.mode(Y) <- "double"
	X <- sobj[["X"]]
	N <- NROW(X)
	storage.mode(X) <- "double"
	p <- NCOL(X)
	events <- sobj[["events"]]
	eventvarlabel <- sobj[["eventvarlabel"]]
	times <- sobj[["times"]]
	timevarlabel <- sobj[["timevarlabel"]]
	scale_sd <- sobj[["scale_sd"]]
	scale_mu <- sobj[["scale_mu"]]
	varnames <- colnames(data)
	
	if(lamfract<0.5 | lamfract > 1)stop("Choose lamfract between 0.5 and 1")
	if(!is.null(colnames(X))){xnames = colnames(X)}else{xnames = paste0("X", 1:p)}
	if (any(alphas > 1) | any(alphas < 0))stop("Choose alphas between 0 to 1.")
	# Reset alpha if > 1
	if (length(alphas)==1 && any(alphas > 1)){alphas_old <- alphas; alphas <- 1} else{alphas_old <- alphas}
	if(!is.null(lambdas) && length(lambdas)<2)stop("Need more than one value of lambda for cross-validation.")
	if(nclusters < 1)stop("Number of clusters must be at least 1.")
	if(!is.null(lambdas)) nlambdas <- length(lambdas)
	if (is.null(foldids)){
		foldids <- sample(rep(seq(nfolds), length.out = N))
	} else {
		nfolds = max(foldids)
	}
	if (nfolds < 3)stop("Number of folds should be at least 3: nfolds = 10 recommended")
	if (!devtype %in% c("basic", "vv", "C"))stop("Likelihood deviance can only be basic or vv. See details.")
	
	# Initialise beta values
	beta0 <- as.double(rep(0, p))
	
	old_lambdas <- lambdas	
	
	# Perform CV
	## Setup parallel because serial takes a lot of time. Otherwise you can turn it off
	nn <- min(parallel::detectCores(), nclusters)
	if (nn < 2){
		foreach::registerDoSEQ()
	} else{
		cl <-  parallel::makeCluster(nn)
		doParallel::registerDoParallel(cl)
		on.exit(parallel::stopCluster(cl))
	}

	if (is.null(old_lambdas)) {
		names(alphas) <- alphas
		maxlambda_alpha <- lapply(alphas, function(al){
			lmax <- proxiterate(Y = Y, X = X, beta0 = beta0, lambda = 1
				, alpha = al, p = p, maxiter = maxiter, tol = tol
				, xnames = xnames, lambmax = TRUE
			)[["max.grad"]]
			if (!is.null(lammin_fract)){
				if (lammin_fract < 0.0001) stop("lammin_fract too small. Consider increasing it.")
				if (lammin_fract > lmax) stop("lammin_fract too large. Consider reducing it.")
				eps <- lammin_fract
			} else {
				if (N > p){eps <- 0.0001} else {eps <- 0.01}
			}
			lmin <- eps*lmax
			lminus <- nlambdas - 1
			lambdas <- lmax*((lmin/lmax)^((0:lminus)/lminus))
			lambdas <- lambdas[1:floor(lamfract*nlambdas)]
			return(lambdas)
		})
	}
	
	tunegrid_df <- expand.grid(folds = 1:nfolds,alpha = alphas)
	fold <- NULL
	cvraw <- foreach (fold = seq(NROW(tunegrid_df)), .combine="rbind", .packages = "pcoxtime"
		, .export = "onefold", .multicombine = TRUE) %dopar% {

		alp <- tunegrid_df[["alpha"]][fold]
		if (is.null(old_lambdas)){
			lam <- maxlambda_alpha[[as.character(alp)]]
		} else {
			lam <- old_lambdas 
		}
		ff <- tunegrid_df[["folds"]][fold]
		index <- which(foldids == ff)
		Y_train <- Y[-index, , drop = FALSE]
		X_train <- X[-index, , drop = FALSE]

		if (devtype == "basic") {
			Y_test <- Y[index, , drop = FALSE]
			X_test <- X[index, , drop = FALSE]
		} else {
			Y_test <- Y
			X_test <- X
		}
		cvraw <- onefold(Y_train = Y_train, X_train = X_train, Y_test = Y_test
			, X_test = X_test, beta0 = beta0, alpha = alp, lambdas = lam
			, devtype = devtype, p = p, tol = tol, xnames = xnames, lambmax = FALSE
			, maxiter = maxiter
		)
		list(alpha = unname(alp), cvraw = cvraw)
	}
	rownames(cvraw) <- NULL
	cvraw <- split(cvraw[,"cvraw"], unlist(cvraw[, "alpha"]))
	temp_alphas <- names(cvraw)
	names(temp_alphas) <- temp_alphas
	out <- lapply(temp_alphas, function(alp){
		cvraw <- do.call(rbind, cvraw[[alp]])
		perfold <- nfolds - apply(is.na(cvraw),2,sum)
		weights <- as.vector(tapply(events,foldids,sum))
		cvraw <- cvraw/weights
		cvm <- apply(cvraw,2,weighted.mean,w=weights,na.rm=TRUE)
		cvsd <- sqrt(apply(scale(cvraw,cvm,FALSE)^2,2,weighted.mean,w=weights,na.rm=TRUE)/(perfold-1))
		if (is.null(old_lambdas)){
			lam <- maxlambda_alpha[[alp]]
		} else {
			lam <- old_lambdas 
		}
		lamin <- getmin(lam,cvm,cvsd)
		cvm_df <- data.frame(lambda = lam, alpha = as.numeric(unname(alp))
			, cvm = cvm, cvsd = cvsd, cvlo = cvm - cvsd, cvup = cvm + cvsd
		)
		min_metrics_df <- data.frame(lambda.min = lamin$lambda.min
			, lambda.1se = lamin$lambda.1se, cv.min = lamin$cv.min
			, alpha = as.numeric(unname(alp))
		)
		res <- list(cvm_df = cvm_df, min_metrics_df = min_metrics_df)
		return(res)	
	})
	out <- do.call("rbind", out)
	cvm_df <- do.call("rbind", out[, "cvm_df"])
	rownames(cvm_df) <- NULL
	min_metrics_df <- do.call("rbind", out[, "min_metrics_df"])
	rownames(min_metrics_df) <- NULL
	
	### Min metrics
	min_lambdas <- min_metrics_df[which.min(min_metrics_df$cv.min), ]

	### Min lambda: optimal lambda
	lambda.min <- min_lambdas$lambda.min
	### 1 std error
	lambda.1se <- min_lambdas$lambda.1se
	### Optimal alpha
	alpha.optimal <- min_lambdas$alpha
	lambdas.optimal <- cvm_df$lambda[cvm_df$alpha==alpha.optimal]

	#### Fit the final model with all the dataset and all opt lambdas and alpha
	if (refit){
		if (!quietly) {cat("Progress: Refitting with optimal lambdas...", "\n")}

		i <- NULL
		beta_refit_df <- foreach (i = 1:length(lambdas.optimal), .combine = "rbind", .packages = "pcoxtime", .export = "proxiterate") %dopar% {
			lam <- lambdas.optimal[[i]]
			fitobj <- proxiterate(Y = Y
				, X = X
				, beta0 = beta0
				, lambda = as.double(lam)
				, alpha =	as.double(alpha.optimal)
				, p = p
				, maxiter = maxiter
				, tol = tol
				, xnames = xnames
				, lambmax = FALSE
			)
			beta_hat <- fitobj$beta_hat/scale_sd
			l1_norm <- sum(abs(beta_hat))
			nzero <- length(beta_hat[beta_hat!=0])
			beta_est <- data.frame(term = names(beta_hat)
				, estimate = unname(beta_hat)
				, alpha = rep(alpha.optimal, p), lambda = rep(lam, p)
				, l1_norm = rep(l1_norm, p), nzero = rep(nzero, p)
			)
			beta_est
		}
	} else {
		beta_refit_df <- NULL
	}
	
	result <- list(lambda.min = lambda.min
		, lambda.1se = lambda.1se
		, alpha.optimal = alpha.optimal
		, lambdas.optimal = lambdas.optimal, foldids = foldids
		, dfs = list(min_metrics_df = min_metrics_df, cvm_df = cvm_df)
		, fit = list(beta = beta_refit_df)
	)
	result$call <- match.call()
	class(result) <- "pcoxtimecv"
	return(result)
}

## Iterate over folds 
onefold <- function(Y_train, X_train, Y_test, X_test, beta0
	, alpha, lambdas, devtype, p, tol, xnames, lambmax, maxiter){
 	
	tryCatch({	
		fitobj <- lambdaiterate(Y = Y_train
			, X = X_train
			, beta0 = as.double(beta0)
			, lambdas = lambdas
			, alpha =	as.double(alpha)
			, p = p
			, maxiter = maxiter
			, tol = tol
			, xnames = xnames
			, lambmax = lambmax
		)
		
		beta_hat <- fitobj$betaL
		plminusk <-  2*fitobj[["min.nloglikL"]]
		
		plfull <- lapply(1:length(lambdas), function(l){
			plfull <- 2*nloglik(Y = Y_test, X = X_test
				, beta0 = as.vector(beta_hat[,l])
				, alpha = alpha
				, lambda = lambdas[l]
			)[["nll.est"]]
		})
		plfull <- as.vector(unlist(plfull))

		if (devtype == "basic"){
			cvraw <- plfull	
		} else {
			cvraw <- plfull - plminusk
		}
		return(cvraw)
	}, error = function(e){
		cat("Possible non-convergence for some lambdas in one of the folds. \nConsider increasing maxiter or reducing tol")
	})
}


## These function is directly copied from the
## glmnet package:
##        Jerome Friedman, Trevor Hastie, Robert Tibshirani (2010).
##        Regularization Paths for Generalized Linear Models via
##        Coordinate Descent.
##        Journal of Statistical Software, 33(1), 1-22.
##        URL http://www.jstatsoft.org/v33/i01/.
## The reason it is copied here is because it is an internal function
## and hence not exported into the global environment.

getmin <- function(lambda, cvm, cvsd) {
	cvmin <- min(cvm)
	idmin <- cvm <= cvmin
	lambda.min <- max(lambda[idmin])
	idmin <- match(lambda.min, lambda)
	semin <- (cvm + cvsd)[idmin]
	idmin <- cvm <= semin
	lambda.1se <- max(lambda[idmin])
	list(lambda.min = lambda.min, lambda.1se = lambda.1se, cv.min = cvmin)
}

