# Generated by using Rcpp::compileAttributes() -> do not edit by hand
# Generator token: 10BE3573-1514-4C36-9D1C-5A225CD40393

threshold_vec <- function(xvec, lambdavec, soft = TRUE) {
    .Call(`_regcell_threshold_vec`, xvec, lambdavec, soft)
}

threshold_mat <- function(xmat, lambdamat, soft = TRUE) {
    .Call(`_regcell_threshold_mat`, xmat, lambdamat, soft)
}

threshold_svd <- function(x, lambdavec, soft = TRUE) {
    .Call(`_regcell_threshold_svd`, x, lambdavec, soft)
}

rob_pca <- function(x, xc, delta, lambda, maxiter = 1000L) {
    .Call(`_regcell_rob_pca`, x, xc, delta, lambda, maxiter)
}

reg_beta <- function(yclean, xclean, betahat, intercept, alambdavec_beta, softbeta = TRUE, tol = 0.001, maxiterbeta = 5) {
    .Call(`_regcell_reg_beta`, yclean, xclean, betahat, intercept, alambdavec_beta, softbeta, tol, maxiterbeta)
}

reg_delta <- function(yclean, x, betahat, deltahat, alambdamat_delta, alpha, softdelta = TRUE, maxiterdelta = 5) {
    .Call(`_regcell_reg_delta`, yclean, x, betahat, deltahat, alambdamat_delta, alpha, softdelta, maxiterdelta)
}

reg_zeta <- function(y, xclean, betahat, alambdavec_zeta, softzeta = TRUE) {
    .Call(`_regcell_reg_zeta`, y, xclean, betahat, alambdavec_zeta, softzeta)
}

reg_beta_delta <- function(y, x, betahat, intercept, deltahat, zetahat, lambda_beta, softbeta, lambda_delta, softdelta, lambda_zeta, softzeta, alpha, tol = 0.001, maxiter = 30) {
    .Call(`_regcell_reg_beta_delta`, y, x, betahat, intercept, deltahat, zetahat, lambda_beta, softbeta, lambda_delta, softdelta, lambda_zeta, softzeta, alpha, tol, maxiter)
}

