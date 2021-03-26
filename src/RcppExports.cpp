// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

// nloglik
Rcpp::List nloglik(const arma::mat& Y, const arma::mat& X, const arma::vec& beta0, const double lambda, const double alpha);
RcppExport SEXP _pcoxtime_nloglik(SEXP YSEXP, SEXP XSEXP, SEXP beta0SEXP, SEXP lambdaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(nloglik(Y, X, beta0, lambda, alpha));
    return rcpp_result_gen;
END_RCPP
}
// gradient
Rcpp::NumericVector gradient(const arma::mat& X, const arma::vec& beta0, const arma::vec& res_est, const double lambda, const double alpha);
RcppExport SEXP _pcoxtime_gradient(SEXP XSEXP, SEXP beta0SEXP, SEXP res_estSEXP, SEXP lambdaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type res_est(res_estSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(gradient(X, beta0, res_est, lambda, alpha));
    return rcpp_result_gen;
END_RCPP
}
// proxupdate
Rcpp::NumericVector proxupdate(const arma::vec& beta0, const arma::vec& grad, const double gamma, const double lambda, const double alpha);
RcppExport SEXP _pcoxtime_proxupdate(SEXP beta0SEXP, SEXP gradSEXP, SEXP gammaSEXP, SEXP lambdaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::vec& >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type grad(gradSEXP);
    Rcpp::traits::input_parameter< const double >::type gamma(gammaSEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(proxupdate(beta0, grad, gamma, lambda, alpha));
    return rcpp_result_gen;
END_RCPP
}
// bbstep
double bbstep(const NumericVector& beta, const NumericVector& beta_prev, const NumericVector& grad, const NumericVector& grad_prev);
RcppExport SEXP _pcoxtime_bbstep(SEXP betaSEXP, SEXP beta_prevSEXP, SEXP gradSEXP, SEXP grad_prevSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type beta(betaSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type beta_prev(beta_prevSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type grad(gradSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type grad_prev(grad_prevSEXP);
    rcpp_result_gen = Rcpp::wrap(bbstep(beta, beta_prev, grad, grad_prev));
    return rcpp_result_gen;
END_RCPP
}
// proxiterate
Rcpp::List proxiterate(const arma::mat& Y, const arma::mat& X, const arma::vec& beta0, const double lambda, const double alpha, const int p, const int maxiter, const double tol, const CharacterVector& xnames, bool lambmax);
RcppExport SEXP _pcoxtime_proxiterate(SEXP YSEXP, SEXP XSEXP, SEXP beta0SEXP, SEXP lambdaSEXP, SEXP alphaSEXP, SEXP pSEXP, SEXP maxiterSEXP, SEXP tolSEXP, SEXP xnamesSEXP, SEXP lambmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type xnames(xnamesSEXP);
    Rcpp::traits::input_parameter< bool >::type lambmax(lambmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(proxiterate(Y, X, beta0, lambda, alpha, p, maxiter, tol, xnames, lambmax));
    return rcpp_result_gen;
END_RCPP
}
// pcoxKKTcheck
LogicalVector pcoxKKTcheck(const NumericVector& grad, const NumericVector& beta0, const double lambda, const double alpha);
RcppExport SEXP _pcoxtime_pcoxKKTcheck(SEXP gradSEXP, SEXP beta0SEXP, SEXP lambdaSEXP, SEXP alphaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const NumericVector& >::type grad(gradSEXP);
    Rcpp::traits::input_parameter< const NumericVector& >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< const double >::type lambda(lambdaSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    rcpp_result_gen = Rcpp::wrap(pcoxKKTcheck(grad, beta0, lambda, alpha));
    return rcpp_result_gen;
END_RCPP
}
// lambdaiterate
Rcpp::List lambdaiterate(const arma::mat& Y, const arma::mat& X, const arma::vec& beta0, const arma::vec& lambdas, const double alpha, const int p, const int maxiter, const double tol, const CharacterVector& xnames, bool lambmax);
RcppExport SEXP _pcoxtime_lambdaiterate(SEXP YSEXP, SEXP XSEXP, SEXP beta0SEXP, SEXP lambdasSEXP, SEXP alphaSEXP, SEXP pSEXP, SEXP maxiterSEXP, SEXP tolSEXP, SEXP xnamesSEXP, SEXP lambmaxSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< const arma::mat& >::type Y(YSEXP);
    Rcpp::traits::input_parameter< const arma::mat& >::type X(XSEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type beta0(beta0SEXP);
    Rcpp::traits::input_parameter< const arma::vec& >::type lambdas(lambdasSEXP);
    Rcpp::traits::input_parameter< const double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< const int >::type p(pSEXP);
    Rcpp::traits::input_parameter< const int >::type maxiter(maxiterSEXP);
    Rcpp::traits::input_parameter< const double >::type tol(tolSEXP);
    Rcpp::traits::input_parameter< const CharacterVector& >::type xnames(xnamesSEXP);
    Rcpp::traits::input_parameter< bool >::type lambmax(lambmaxSEXP);
    rcpp_result_gen = Rcpp::wrap(lambdaiterate(Y, X, beta0, lambdas, alpha, p, maxiter, tol, xnames, lambmax));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_pcoxtime_nloglik", (DL_FUNC) &_pcoxtime_nloglik, 5},
    {"_pcoxtime_gradient", (DL_FUNC) &_pcoxtime_gradient, 5},
    {"_pcoxtime_proxupdate", (DL_FUNC) &_pcoxtime_proxupdate, 5},
    {"_pcoxtime_bbstep", (DL_FUNC) &_pcoxtime_bbstep, 4},
    {"_pcoxtime_proxiterate", (DL_FUNC) &_pcoxtime_proxiterate, 10},
    {"_pcoxtime_pcoxKKTcheck", (DL_FUNC) &_pcoxtime_pcoxKKTcheck, 4},
    {"_pcoxtime_lambdaiterate", (DL_FUNC) &_pcoxtime_lambdaiterate, 10},
    {NULL, NULL, 0}
};

RcppExport void R_init_pcoxtime(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
