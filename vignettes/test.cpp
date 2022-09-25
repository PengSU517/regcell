// -*- mode: C++; c-indent-level: 4; c-basic-offset: 4; indent-tabs-mode: nil; -*-
// we only include RcppArmadillo.h which pulls Rcpp.h in for us
// Rcpp::sourceCpp("vignettes/test.cpp")
// Rcpp::compileAttributes()

// [[Rcpp::plugins(cpp11)]]
// [[Rcpp::depends(RcppArmadillo)]]
#include "RcppArmadillo.h"

using namespace arma;
using namespace Rcpp;





// [[Rcpp::export]]
arma::vec diff_func(arma::vec res, double b = 4.6845){   //should be a vector
  // rho function Tukey's biweight

  unsigned int n = res.n_elem;

  arma::vec bvec(n);
  bvec.fill(b);

  arma::vec value(n);
  arma::vec resb = res/bvec;

  value = (abs(res) <= bvec)%(1 - 6*pow(resb,2) + 5*pow(resb,4));

  return value;

}





