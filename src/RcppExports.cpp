// Generated by using Rcpp::compileAttributes() -> do not edit by hand
// Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

#include <RcppArmadillo.h>
#include <Rcpp.h>

using namespace Rcpp;

#ifdef RCPP_USE_GLOBAL_ROSTREAM
Rcpp::Rostream<true>&  Rcpp::Rcout = Rcpp::Rcpp_cout_get();
Rcpp::Rostream<false>& Rcpp::Rcerr = Rcpp::Rcpp_cerr_get();
#endif

// threshold_beta
arma::vec threshold_beta(arma::vec betahat, arma::vec alambda_beta, bool softbeta);
RcppExport SEXP _regcell_threshold_beta(SEXP betahatSEXP, SEXP alambda_betaSEXP, SEXP softbetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type alambda_beta(alambda_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type softbeta(softbetaSEXP);
    rcpp_result_gen = Rcpp::wrap(threshold_beta(betahat, alambda_beta, softbeta));
    return rcpp_result_gen;
END_RCPP
}
// threshold_delta
arma::mat threshold_delta(arma::mat deltahat, arma::mat alambda_delta, bool softdelta);
RcppExport SEXP _regcell_threshold_delta(SEXP deltahatSEXP, SEXP alambda_deltaSEXP, SEXP softdeltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::mat >::type deltahat(deltahatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type alambda_delta(alambda_deltaSEXP);
    Rcpp::traits::input_parameter< bool >::type softdelta(softdeltaSEXP);
    rcpp_result_gen = Rcpp::wrap(threshold_delta(deltahat, alambda_delta, softdelta));
    return rcpp_result_gen;
END_RCPP
}
// reg_beta
List reg_beta(arma::vec y, arma::mat x, arma::vec betahat, double intercept, arma::vec lambdavec_beta, bool softbeta, double maxiterbeta);
RcppExport SEXP _regcell_reg_beta(SEXP ySEXP, SEXP xSEXP, SEXP betahatSEXP, SEXP interceptSEXP, SEXP lambdavec_betaSEXP, SEXP softbetaSEXP, SEXP maxiterbetaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< double >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type lambdavec_beta(lambdavec_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type softbeta(softbetaSEXP);
    Rcpp::traits::input_parameter< double >::type maxiterbeta(maxiterbetaSEXP);
    rcpp_result_gen = Rcpp::wrap(reg_beta(y, x, betahat, intercept, lambdavec_beta, softbeta, maxiterbeta));
    return rcpp_result_gen;
END_RCPP
}
// reg_delta
List reg_delta(arma::vec y, arma::mat x, arma::vec betahat, arma::mat deltahat, double lambda_delta, arma::mat cellweight, double alpha, bool softdelta, double maxiterdelta);
RcppExport SEXP _regcell_reg_delta(SEXP ySEXP, SEXP xSEXP, SEXP betahatSEXP, SEXP deltahatSEXP, SEXP lambda_deltaSEXP, SEXP cellweightSEXP, SEXP alphaSEXP, SEXP softdeltaSEXP, SEXP maxiterdeltaSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type deltahat(deltahatSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_delta(lambda_deltaSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type cellweight(cellweightSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< bool >::type softdelta(softdeltaSEXP);
    Rcpp::traits::input_parameter< double >::type maxiterdelta(maxiterdeltaSEXP);
    rcpp_result_gen = Rcpp::wrap(reg_delta(y, x, betahat, deltahat, lambda_delta, cellweight, alpha, softdelta, maxiterdelta));
    return rcpp_result_gen;
END_RCPP
}
// reg_beta_delta
List reg_beta_delta(arma::vec y, arma::mat x, arma::vec betahat, double intercept, arma::mat deltahat, arma::vec rowweight, arma::mat colweight, double lambda_beta, bool adabeta, bool softbeta, double lambda_delta, bool adadelta, bool softdelta, double alpha, double maxiter);
RcppExport SEXP _regcell_reg_beta_delta(SEXP ySEXP, SEXP xSEXP, SEXP betahatSEXP, SEXP interceptSEXP, SEXP deltahatSEXP, SEXP rowweightSEXP, SEXP colweightSEXP, SEXP lambda_betaSEXP, SEXP adabetaSEXP, SEXP softbetaSEXP, SEXP lambda_deltaSEXP, SEXP adadeltaSEXP, SEXP softdeltaSEXP, SEXP alphaSEXP, SEXP maxiterSEXP) {
BEGIN_RCPP
    Rcpp::RObject rcpp_result_gen;
    Rcpp::RNGScope rcpp_rngScope_gen;
    Rcpp::traits::input_parameter< arma::vec >::type y(ySEXP);
    Rcpp::traits::input_parameter< arma::mat >::type x(xSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type betahat(betahatSEXP);
    Rcpp::traits::input_parameter< double >::type intercept(interceptSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type deltahat(deltahatSEXP);
    Rcpp::traits::input_parameter< arma::vec >::type rowweight(rowweightSEXP);
    Rcpp::traits::input_parameter< arma::mat >::type colweight(colweightSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_beta(lambda_betaSEXP);
    Rcpp::traits::input_parameter< bool >::type adabeta(adabetaSEXP);
    Rcpp::traits::input_parameter< bool >::type softbeta(softbetaSEXP);
    Rcpp::traits::input_parameter< double >::type lambda_delta(lambda_deltaSEXP);
    Rcpp::traits::input_parameter< bool >::type adadelta(adadeltaSEXP);
    Rcpp::traits::input_parameter< bool >::type softdelta(softdeltaSEXP);
    Rcpp::traits::input_parameter< double >::type alpha(alphaSEXP);
    Rcpp::traits::input_parameter< double >::type maxiter(maxiterSEXP);
    rcpp_result_gen = Rcpp::wrap(reg_beta_delta(y, x, betahat, intercept, deltahat, rowweight, colweight, lambda_beta, adabeta, softbeta, lambda_delta, adadelta, softdelta, alpha, maxiter));
    return rcpp_result_gen;
END_RCPP
}

static const R_CallMethodDef CallEntries[] = {
    {"_regcell_threshold_beta", (DL_FUNC) &_regcell_threshold_beta, 3},
    {"_regcell_threshold_delta", (DL_FUNC) &_regcell_threshold_delta, 3},
    {"_regcell_reg_beta", (DL_FUNC) &_regcell_reg_beta, 7},
    {"_regcell_reg_delta", (DL_FUNC) &_regcell_reg_delta, 9},
    {"_regcell_reg_beta_delta", (DL_FUNC) &_regcell_reg_beta_delta, 15},
    {NULL, NULL, 0}
};

RcppExport void R_init_regcell(DllInfo *dll) {
    R_registerRoutines(dll, NULL, CallEntries, NULL, NULL);
    R_useDynamicSymbols(dll, FALSE);
}
