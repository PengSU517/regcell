// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
// Rcpp::sourceCpp("src/srlm.cpp")
// Rcpp::compileAttributes()

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"

using namespace arma;
using namespace Rcpp;


// [[Rcpp::export]]
arma::vec threshold_vec(arma::vec betahat, arma::vec alambda_beta, bool soft = true){   //should be a vector
  double p = betahat.n_rows;
  arma::vec diff = abs(betahat) - alambda_beta;
  arma::vec value(p);
  if(soft){
    value = sign(betahat)%(diff>0)%diff;
  }else{
    value = betahat%(diff>0);
  }
  return value;
}


// [[Rcpp::export]]
arma::mat threshold_mat(arma::mat deltahat, arma::mat alambda_delta, bool soft = true){   //should be a vector
  double n = deltahat.n_rows;
  double p = deltahat.n_cols;
  arma::mat value(n,p);
  arma::mat diff = abs(deltahat) - alambda_delta;
  if(soft){
    value = sign(deltahat)%(diff>0)%diff;
  }else{
    value = deltahat%(diff>0);
  }
  return value;
}


// [[Rcpp::export]]
arma::mat threshold_svd(arma::mat x, arma::vec lambdavec, bool soft = true){   //should be a vector
  double n = x.n_rows;
  double p = x.n_cols;
  double dim = lambdavec.n_elem;

  arma::mat U;
  arma::vec s;
  arma::mat V;

  svd(U, s, V, x);
  arma::vec snew(dim);
  snew = threshold_vec(s, lambdavec);//什么毛病

  arma::mat S(n,p);
  S.diag() =  snew;
  arma::mat xc = U*S*trans(V);

  return xc;
}


// [[Rcpp::export]]
List rob_pca(arma::mat x, arma::mat xc, arma::mat delta, double lambda, int maxiter = 10){
  unsigned int n = x.n_rows;
  unsigned int p = x.n_cols;
  unsigned int dim = std::min(n,p);

  arma::mat zeta(n,p);
  zeta.fill(0);

  double L = 1;
  arma::mat lmat(n,p);
  lmat.fill(L);

  double stepsize = 1/L;
  arma::mat stepmat(n,p);
  stepmat.fill(stepsize);

  arma::vec stepvec(dim);
  stepvec.fill(stepsize);

  arma::mat lambdamat(n,p);
  lambdamat.fill(lambda);

  int k = 1;
  while(k < maxiter){
    delta = threshold_mat(x - xc + stepmat%zeta, stepmat%lambdamat);
    xc = threshold_svd(x - delta + stepmat%zeta, stepvec);
    zeta = zeta + lmat%(x - xc - delta);
    if(abs(x-xc-delta).max() < pow(10, -8)){break;}
    k = k+1;
  }

  List results = List::create(
    Named("xc") = xc,
    Named("delta") = delta
  );

  return results;
}

