# include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]

using namespace Rcpp;
using namespace arma;

// Compute negative log-likelihood of time-dependant covariates
// [[Rcpp::export]]
Rcpp::List nloglik(const arma::mat& Y, const arma::mat& X, const arma::vec& beta0
  , const double lambda, const double alpha) {
  
  // Number of observations
  int N = X.n_rows;
  
  // Number of time variables (defined by Surv)
  int p = Y.n_cols;
  
  // Linear predictor
  arma::vec eta = X * beta0;
  // Relative hazard
  arma::vec relhaz = arma::exp(eta);
  
  // Risk score for each patient
  arma::vec risk = zeros<vec>(N);
  // Individuals at risk
  arma::vec n_risk = zeros<vec>(N);
  // Residuals (events - P_mat %*% events)
  arma::vec res_est(N);

  // Sum of log-likelihood
  double ll = 0.0;
  double nll;
  double flog;
  
  arma::vec starttime;
  arma::vec endtime;
  arma::vec events;
  Rcpp::NumericVector cond;
  arma::rowvec P;

  if (p == 2) {
    endtime = Y.col(0);
    events = Y.col(1);
    for (int i = 0; i < N; i++){
      cond = (endtime(i)) >= (endtime);
      risk += cond * relhaz(i);
      n_risk += cond;
    }

    // Compute log-likelihood
    for (int i = 0; i < N; i++){
      cond = (endtime(i)) >= (endtime);
      P = Rcpp::as<Rcpp::NumericVector>(wrap(relhaz(i)/risk)) * cond;
      res_est[i] = events(i) - as<double>(wrap(P * events));
      if (events[i] == 1){
        ll += log(as<double>(wrap(P(i))));
      }
    }
  } else {
    starttime = Y.col(0);
    endtime = Y.col(1);
    events = Y.col(2);
    for (int i = 0; i < N; i++){
      cond = ((endtime(i)) >= (endtime)) && (starttime(i) < endtime);
      risk += cond * relhaz(i);
      n_risk += cond;
    }
    
    // Compute log-likelihood
    for (int i = 0; i < N; i++){
      cond = ((endtime(i)) >= (endtime)) && (starttime(i) < endtime);
      P = Rcpp::as<Rcpp::NumericVector>(wrap(relhaz(i)/risk)) * cond;
      res_est[i] = events(i) - as<double>(wrap(P * events));
      if (events[i] == 1){
        ll += log(as<double>(wrap(P(i))));
      }
    }
  }

  // Penalized negative log-likelihood
  nll = -ll + lambda*(alpha*sum(abs(beta0)) + 0.5*(1-alpha)*sum(square(beta0)));
  flog = -ll; 
  return Rcpp::List::create(Rcpp::Named("res.est", wrap(res_est))
    , Rcpp::Named("nll.est", wrap(nll))
    , Rcpp::Named("n.risk", wrap(n_risk))
    , Rcpp::Named("flog", wrap(flog))
  );
}


// Compute the gradient 
// [[Rcpp::export]]
Rcpp::NumericVector gradient(const arma::mat& X, const arma::vec& beta0,  const arma::vec& res_est, const double lambda, const double alpha){
  // Number of observations
  double N = X.n_rows;
  Rcpp::NumericVector grad;
  grad = -1/N * (X.t() * res_est) + lambda*(1-alpha)*beta0;
  return(grad);
}

//  Defines operator for proximal gradient descent update
// [[Rcpp::export]]
Rcpp::NumericVector proxupdate(const arma::vec& beta0, const  arma::vec& grad, const double gamma, const double lambda, const double alpha){
  arma::vec z = beta0 - gamma * grad;
  Rcpp::NumericVector z2 = Rcpp::as<Rcpp::NumericVector>(wrap(z));
  double theta = gamma * lambda * alpha;
  Rcpp::NumericVector betak = Rcpp::as<Rcpp::NumericVector>(sign(z2)) * (pmax(abs(z2) - theta, 0.0));
  return(betak);
}


// Barzilia-Borwien step size adjustment
// [[Rcpp::export]]
double bbstep(const Rcpp::NumericVector& beta, const Rcpp::NumericVector& beta_prev, const Rcpp::NumericVector& grad, const Rcpp::NumericVector& grad_prev){
  Rcpp::NumericVector sk = beta - beta_prev;
  Rcpp::NumericVector yk = grad - grad_prev;
  double gamma = std::max(sum(sk*yk)/sum(pow(yk, 2)), sum(pow(sk,2))/sum(sk*yk));
  return(gamma);
}


// Proximal iteration updates
// [[Rcpp::export]]
Rcpp::List proxiterate(const arma::mat& Y, const arma::mat& X
  , const arma::vec& beta0, const double lambda, const double alpha
  , const int p, const int maxiter, const double tol
  , const CharacterVector& xnames, bool lambmax = false){
  
  // Intialize Beta
  Rcpp::NumericMatrix beta(p, maxiter+1);
  rownames(beta) = xnames;
  
  // Initialize gradients
  Rcpp::NumericMatrix grads(p, maxiter+1);
  rownames(grads) = xnames;
  
  // Negative log-likelihoods
  arma::vec nlls(maxiter+1);

  // Beta deviance
  arma::vec deviances(maxiter+1);

  // Step sizes
  arma::vec gamma(maxiter+1);
  
  // Barzilia-Borwien step initialization
  List nll_init = nloglik(Y, X, beta0, lambda, alpha);
  arma::vec res_est_init = nll_init["res.est"];
  Rcpp::NumericVector grad_init = gradient(X, beta0, res_est_init, lambda, alpha);
  
  // Solution path using lambda which give all B_hat = 0
  double gtol = max(abs(grad_init));
  if (alpha == 0.0){ gtol = gtol/1e-3;} else {gtol = gtol/alpha; }
  if (lambmax){
    return Rcpp::List::create(Named("max.grad", gtol));
  }

  // Initialization (0.01)
  beta(_, 0) = proxupdate(beta0, grad_init, 0.001*gtol, lambda, alpha);
  List nll1 = nloglik(Y, X, beta(_, 0), lambda, alpha);
  grads(_, 0) = gradient(X, beta(_, 0), nll1["res.est"], lambda, alpha);
  gamma(0) = bbstep(beta(_, 0), Rcpp::as<Rcpp::NumericVector>(wrap(beta0)), grads(_,0), grad_init);
  nlls(0) = nll1["nll.est"];
  
  // Convergence message
  std::string message = "Model did not converge after " + std::to_string(maxiter) + " iterations, consider increasing maxiter...";
  
  // Proximal updates and iterations
  Rcpp::NumericVector beta_hat(p);
  double devcheck;
  List nll;
  int iter = 0;
  for ( int k = 0; k < maxiter; k++){
    if (lambda >= gtol){
      double gamma2 = 1;
      nll = nloglik(Y, X, beta(_, k), lambda, alpha);
      grads(_, k+1) = gradient(X, beta(_, k), nll["res.est"], lambda, alpha);
      beta(_, k+1) = proxupdate(beta(_, k), grads(_, k+1), gamma2, lambda, alpha);
      
      // Check for convergence
      // devcheck = sqrt(sum(pow(beta(_,k+1) - beta(_, k), 2)))/gamma2;
		devcheck = max(abs(beta(_,k+1) - beta(_, k)))/gamma2;
    } else {
      beta(_, k+1) = proxupdate(beta(_, k), grads(_, k), gamma(k), lambda, alpha);
      nll = nloglik(Y, X, beta(_, k+1), lambda, alpha);
      grads(_, k+1) = gradient(X, beta(_, k+1), nll["res.est"], lambda, alpha);
      gamma(k+1) = bbstep(beta(_, k+1), Rcpp::as<Rcpp::NumericVector>(wrap(beta(_, k))), grads(_,k+1), grads(_,k));
      // Check for convergence
		// devcheck = sqrt(sum(pow((beta(_,k+1) - beta(_, k)), 2)))/gamma(k);
		devcheck = max(abs(beta(_,k+1) - beta(_, k)))/gamma(k);
    }
    deviances(k) = devcheck;
    nlls(k+1) = nll["nll.est"];
	 // devcheck = abs(nlls(k+1) - nlls(k))/(abs(nlls(k)) + tol);
    iter += 1;
    if (devcheck <= tol || lambda >= gtol){
      beta_hat = beta(_, k+1);
      message = "Model converged after " + std::to_string(k+2) + " iterations";
      break;
    }
  }
  
  beta_hat.attr("names") = xnames;
  Rcpp::NumericVector grad_opt = grads(_, iter);
  grad_opt.attr("names") = xnames;
  Rcpp::NumericVector n_risk = nll["n.risk"];
  n_risk.attr("dim") = R_NilValue;
  
  return Rcpp::List::create(Rcpp::Named("beta_hat", wrap(beta_hat))
    , Rcpp::Named("grad.opt", wrap(grad_opt))
    , Rcpp::Named("n.risk", wrap(n_risk))
    , Rcpp::Named("min.dev", wrap(devcheck))
    , Rcpp::Named("deviances", wrap(deviances.rows(0, iter)))
    , Rcpp::Named("min.nloglik", wrap(nll["nll.est"]))
    , Rcpp::Named("nlogliks", wrap(nlls.rows(0, iter)))
    , Rcpp::Named("flog", wrap(nll["flog"]))
    , Rcpp::Named("message", wrap(message))
  );
}

// Check KKT violations
// [[Rcpp::export]]
LogicalVector pcoxKKTcheck(const Rcpp::NumericVector& grad, const Rcpp::NumericVector& beta0, const double lambda, const double alpha) {
  LogicalVector index1 = abs(grad) >= (alpha * lambda);
  LogicalVector index2 = beta0 == 0.0;
  LogicalVector index = index1 & index2;
  return(index);
}

// Iterate over a vector of lambda
// [[Rcpp::export]]
Rcpp::List lambdaiterate(const arma::mat& Y, const arma::mat& X
	, const arma::vec& beta0, const arma::vec& lambdas, const double alpha
	, const int p, const int maxiter, const double tol
	, const CharacterVector& xnames, bool lambmax = false){
  
	// Number of lambdas
	int nlambda = lambdas.n_elem;
  
	// Intialize Beta-lambda matrix
	Rcpp::NumericMatrix betaL(p, nlambda);
	rownames(betaL) = xnames;

	// Min Negative log-likelihoods
	Rcpp::NumericVector min_nloglikL(nlambda);

	// Cross-validation Min Negative log-likelihoods
 	Rcpp::NumericVector flogL(nlambda);
	
	for (int l = 0; l < nlambda; l++) {
		List out = proxiterate(Y, X, beta0, lambdas(l), alpha, p, maxiter, tol, xnames, lambmax = false);
		betaL(_, l) = Rcpp::as<Rcpp::NumericVector>(wrap(out["beta_hat"]));
		min_nloglikL[l] = as<double>(wrap(out["min.nloglik"]));
		flogL[l] = as<double>(wrap(out["flog"]));
	}

	return Rcpp::List::create(Rcpp::Named("betaL", wrap(betaL))
		, Rcpp::Named("min.nloglikL", wrap(min_nloglikL))
		, Rcpp::Named("flogL", wrap(flogL))
	);

}

//// Compute predicted relative hazard for every individual
//// [[Rcpp::export]]
//NumericVector predictedHazard(const arma::mat& Y, const arma::mat& X, const arma::vec& beta_hat) {
//  
//  // Number of observations
//  int N = X.n_rows;
//  
//  // Number of time variables (defined by Surv)
//  int p = Y.n_cols;
//  
//  // Linear predictor
//  arma::vec eta = vectorise(X * beta_hat);
//  // Relative hazard
//  arma::vec relhaz = arma::exp(eta);
//  
//  // Risk score for each patient
//  arma::vec risk = zeros<vec>(N);
//  
//  arma::vec starttime;
//  arma::vec endtime;
//  arma::vec events;
//  NumericVector cond;
//  
//  if (p == 2) {
//    endtime = Y.col(0);
//    events = Y.col(1);
//    for (int i = 0; i < N; i++){
//      cond = (endtime(i)) >= (endtime);
//      risk += cond * relhaz(i);
//    }
//  } else {
//    starttime = Y.col(0);
//    endtime = Y.col(1);
//    events = Y.col(2);
//    for (int i = 0; i < N; i++){
//      cond = ((endtime(i)) >= (endtime)) && (starttime(i) < endtime);
//      risk += cond * relhaz(i);
//    }
//  }
//  
//  return(Rcpp::as<Rcpp::NumericVector>(wrap(risk)));
//}

