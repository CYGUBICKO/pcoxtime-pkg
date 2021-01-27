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
#' @param refit logical. Whether to return solution path based on optimal lambda and alpha picked by the model. Default is \code{refit = TRUE}.
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
#' data(veteran, package="survival")
#' \dontrun{
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

pcoxtimecv <- function(formula = formula(data), data = sys.parent()
	, alphas = 1, lambdas = NULL, nlambdas = 100, lammin_fract = NULL
	, lamfract = 0.6, nfolds = 10, foldids = NULL, devtype = "vv"
	, refit = TRUE, maxiter = 1e5, tol = 1e-8, quietly = FALSE
	, seed = NULL, nclusters = 1, na.action = na.omit, ...) {
	
	if(is.null(seed)){seed = 1254}
	set.seed(seed)
	
	# survival package data format
	sobj <- riskset(formula = formula, data = data, na.action = na.action)
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
	if (any(alphas > 1))stop("Choose alphas between 0 to 1.")
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
	if (nn < 3){
		foreach::registerDoSEQ()
	} else{
		cl <-  parallel::makeCluster(nn)
		doParallel::registerDoParallel(cl)
		on.exit(parallel::stopCluster(cl))
	}
	
	if (is.null(old_lambdas)){
		## Generate a sequence of lambdas based on a given alpha
		maxlambda_df <- lapply(alphas, function(al){
			lmax <- proxiterate(Y = Y, X = X, beta0 = beta0, lambda = 1
				, alpha = al, p = p, maxiter = maxiter, tol = tol
				, xnames = xnames, lambmax = TRUE
			)[["max.grad"]]
			if (!is.null(lammin_fract)){
				if (lammin_fract < 0.0001) stop("lammin_fract too small. Consider increasing it.")
				if (lammin_fract > lmax) stop("lammin_fract too large. Consider reducing it.")
				eps <- lammin_fract
			} else {
				if (N >= p){eps <- 0.0001} else {eps <- 0.01}
			}
			lmin <- eps*lmax
			lminus <- nlambdas - 1
			lambdas <- lmax*((lmin/lmax)^((0:lminus)/lminus))
			#lambdas <- exp(seq(log(lmax), log(lmin),  length.out = nlambdas))
			ll_df <- expand.grid(folds = 1:nfolds, lambda = lambdas[1:floor(lamfract*nlambdas)]
				, alpha = al, cvraw = 0, weights = 0
			)
		})
		tunegrid_df <- do.call("rbind", maxlambda_df)
	} else {
		tunegrid_df <- expand.grid(folds = 1:nfolds, lambda = lambdas[1:floor(lamfract*nlambdas)]
			, alpha = alphas, cvraw = 0, weights = 0
		)
	}
	fold = NULL
	out <- foreach (fold = seq(NROW(tunegrid_df)), .combine = "rbind"
		, .packages = "pcoxtime", .export = "cvfitfun") %dopar% {
		fitcv <- cvfitfun(fold = fold, Y = Y, X = X,  events = events, beta0 = beta0
			, foldids = foldids, tunegrid_df = tunegrid_df, devtype = devtype, p = p
			, maxiter = maxiter, tol = tol, xnames = xnames, lambmax = FALSE, scale_sd = scale_sd
		)
		beta_df <- data.frame(term = names(fitcv$beta_hat) 
			, fold = rep(fitcv$fold, p), estimate = unname(fitcv$beta_hat) 
			, alpha = rep(fitcv$alpha, p), lambda = rep(fitcv$lambda, p)
			, l1_norm = rep(fitcv$l1_norm, p), nzero = rep(fitcv$ncoefs, p), kkt_pass = fitcv$kkt_pass
		)
		tune_df <- data.frame(folds = fitcv$fold, lambda = fitcv$lambda
			, alpha = fitcv$alpha, cvraw = fitcv$cvraw, weights = fitcv$weights, kkt_pass = fitcv$kkt_pass
		)
		list(tune_df = tune_df, beta_df = beta_df)
	}

	tune_df <- do.call("rbind", out[,1])
	rownames(tune_df) <- NULL
	tune_df <- tune_df[!is.na(tune_df$cvraw)&!is.nan(tune_df$cvraw)&!is.infinite(tune_df$cvraw), ]
	cvm_df <- extractcv(tune_df)
	beta_df <- do.call("rbind", out[, 2])
	rownames(beta_df) <- NULL

	## Optimal parameter values
	min_lambdas_df <- sapply(split(cvm_df, cvm_df$alpha), function(d){
		min_df <- getmin(d$lambda, d$cvm, d$cvsd)
		min_df$alpha <- unique(d$alpha)
		return(unlist(min_df))
	}, simplify = FALSE)
	min_metrics_df <- data.frame(do.call("rbind", min_lambdas_df))
	rownames(min_metrics_df) <- NULL
	
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
	
	out <- list(lambda.min = lambda.min
		, lambda.1se = lambda.1se
		, alpha.optimal = alpha.optimal
		, lambdas.optimal = lambdas.optimal, foldids = foldids
		, dfs = list(min_metrics_df = min_metrics_df, cvm_df = cvm_df, beta = beta_df)
		, fit = list(beta = beta_refit_df)
	)
	out$call <- match.call()
	class(out) <- "pcoxtimecv"
	return(out)
}

cvfitfun <- function(fold, Y, X, events, beta0
	, foldids, tunegrid_df, devtype, p, tol
	, xnames, lambmax, scale_sd, maxiter){
 	
	alp <- tunegrid_df[["alpha"]][fold]
	lam <- tunegrid_df[["lambda"]][fold]
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
	tryCatch({	
		fitobj <- proxiterate(Y = Y_train
			, X = X_train
			, beta0 = as.double(beta0)
			, lambda = as.double(lam)
			, alpha =	as.double(alp)
			, p = p
			, maxiter = maxiter
			, tol = tol
			, xnames = xnames
			, lambmax = lambmax
		)
		
		beta_hat <- fitobj$beta_hat
		plminusk <-  2*fitobj[["flog"]]
		
		plfull <- 2*nloglik(Y = Y_test, X = X_test
			, beta0 = as.double(beta_hat)
			, alpha = as.double(alp)
			, lambda = as.double(lam)
		)[["flog"]]
		
		if (devtype == "basic"){
			cvraw <- plfull	
		} else {
			cvraw <- plfull - plminusk
		}

		## Check KKT conditions
		kkt_violations <- pcoxKKTcheck(fitobj[["grad.opt"]]
			, beta_hat
			, lam
			, alp
		)
		kkt_pass <- TRUE
		if (any(kkt_violations)){
			kkt_pass <- FALSE
		}

		## Nonzero coefficients
		ncoefs <- length(beta_hat[beta_hat != 0])
		
		## Compute weights used in ll deviance. Borrowed from glmnet code
		weights <- sum(events[index])

		beta_hat <- beta_hat #/scale_sd
		l1_norm <- sum(abs(beta_hat), na.rm = TRUE)
		
		out <- list(fold = ff, alpha = alp, lambda = lam, beta_hat = beta_hat
			, cvraw = cvraw, ncoefs = ncoefs, weights = weights, l1_norm = l1_norm, kkt_pass = kkt_pass
		)
		return(out)
	}, error = function(e){
		cat("Possible non-convergence for lambda = ", lam, "in one of the folds. \nConsider increasing maxiter or reducing tol")
	})
}


## Compute Partial likelihood-deviance for each alpha
## Adopted from glmnet
extractcv <- function(cvraw_df){
	cv_df <- sapply(split(cvraw_df, cvraw_df$alpha), function(d1){
		cv_df <- sapply(split(d1, d1$lambda), function(d){
			lambda <- unique(d$lambda)
			alpha <- unique(d$alpha)
			nfolds <- length(unique(d$folds)) - (sum(is.na(d$cvraw)) + sum(is.nan(d$cvraw)))
			cvm <- weighted.mean(d$cvraw/d$weights, w = d$weights)
			cvsd <- sqrt(weighted.mean(scale(d$cvraw/d$weights, cvm, FALSE)^2
				, w = d$weights, na.rm = TRUE)/(nfolds-1)
			)
			cve_df <- data.frame(lambda = lambda, alpha = alpha, cvm = cvm
				, cvsd = cvsd ,cvlo = cvm - cvsd, cvup = cvm + cvsd
			)
			return(cve_df)
		}, simplify = FALSE)
		cv_df <- do.call("rbind", cv_df)
	}, simplify = FALSE)
	cv_df <- do.call("rbind", cv_df)
	rownames(cv_df) <- NULL
	return(cv_df)
}


## Minimum cve: copied from glmnet because it's internal function
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
