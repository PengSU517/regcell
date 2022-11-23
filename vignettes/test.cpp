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
arma::vec rho_func(arma::vec res, double b = 4.6845){   //should be a vector
  // rho function Tukey's biweight

  unsigned int n = res.n_elem;

  arma::vec bvec(n);
  bvec.fill(b);

  arma::vec value(n);
  value = ((abs(res) <= bvec)%(1-pow(  (1 - pow(res/bvec,2.0)),3.0))%pow(bvec,2.0)/6) + (abs(res) > bvec)%pow(bvec,2.0)/6;
  return value;
}


// [[Rcpp::export]]
arma::vec psi_func(arma::vec res, double b = 4.6845){
  // rho function Tukey's biweight

  unsigned int n = res.n_elem;

  arma::vec bvec(n);
  bvec.fill(b);

  arma::vec value(n);
  value = (abs(res) <= bvec)%res%(pow((1 - pow(res/bvec,2)),2));

  return value;
}


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


// [[Rcpp::export]]

List TukeyM(arma::vec y, arma::mat x, arma::vec betahat, double sigmahat, double maxiter = 30){
  //sparse M-regressioner with Tukey's biweight loss and atan penalty
  //y: response
  //x: design matrix
  //betahat: initial estimate
  //sigmahat:estimate of sigma (keeping fixed)
  //maxiter: maximum iteration time.


  //initialize variables
  double n = x.n_rows; //sample size
  double p = x.n_cols; //dimension size

  arma::mat onevec(n,1); //constant vector (first column of the design matrix)
  onevec.fill(1);

  arma::vec sigmahatvec(n);//a vector filled with sigma, just for computation convenience
  sigmahatvec.fill(sigmahat);

  arma::mat sigmahatmat(p+1, p+1);//a matrix filled with sigma, just for computation convenience
  sigmahatmat.fill(sigmahat);

  arma::mat xc = join_rows(onevec, x);// complete design matrix with the first column constant


  arma::mat stepsize(p+1,p+1); //stepsize for proximal gradient descent algorithm
  stepsize.fill(1);//stepsize setting need to be improved

  arma::mat shrinksize(p+1,p+1);
  shrinksize.fill(1);

  arma::vec res = y - xc*betahat;
  arma::vec betanew(p+1);

  arma::mat weightmat = repmat(diff_func(res/sigmahatvec),1,p+1);

  double L = arma::eig_sym(trans(xc)*(weightmat%xc)).max()/pow(sigmahat,2);
  double t = 1/L;  //stepsize could be improved

  arma::vec steptvec(p+1);// a constant vector filled with 1/pi
  steptvec.fill(t);

  double k = 1;//iterator
  while(k <= maxiter)//start iteration
  {
    res = y - xc*betahat;//residuals
    betanew = betahat + steptvec%(trans(xc) * psi_func(res/sigmahatvec));
    if((abs(betanew-betahat)).max() < pow(10, -4)) {break;};
    stepsize = stepsize*shrinksize;
    betahat = betanew;
    k++;
  }
  double loss = 2*accu(rho_func(y-xc*betahat)/sigmahatvec);

  List results = List::create(
    Named("betahat") = betahat,
    Named("sigmahat") = sigmahat,
    Named("loss") = loss,
    Named("t") = t,
    Named("L") = L,
    Named("weightmat") = weightmat,
    Named("k") = k-1
  );
  return(results);

}


// [[Rcpp::export]]
List TukeyM_iter(arma::vec y, arma::mat x, arma::mat ximp, arma::vec betahat, double sigmahat, String tech = "row",
                 double b = 4.6845,  double maxiter = 30){

  //initialize variables
  double n = x.n_rows; //sample size
  double p = x.n_cols; //dimension size
  arma::mat xtilde = ximp;

  arma::mat onevec = ones(n,1); //constant vector (first column of the design matrix)
  arma::mat xc    = join_rows(onevec, x);
  arma::mat xtildec(n,p+1);
  xtildec = join_rows(onevec, xtilde);

  List outputs = TukeyM(y, xtilde, betahat, sigmahat);
  betahat = as<arma::vec>(outputs["betahat"]);

  int i,j,k;

  if(tech=="row"){
    arma::vec res_ori(n);
    arma::vec res_tilde(n);

    arma::vec indexn(n);
    arma::vec indexp(p);
    for(k = 1; k < maxiter; k++){
      xtildec = join_rows(onevec, xtilde);
      res_ori =   y -      xc*betahat;
      res_tilde = y - xtildec*betahat;

      indexn = ((abs(res_ori) - abs(res_tilde)) <  zeros(n))%ones(n);
      //need to improve this condition!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      if(accu(indexn)==0){break;}
      indexp = (abs(betahat.rows(1,p)) > zeros(p))%ones(p);

      for(i=0; (i<n); i++){
        for(j=0; (j<p); j++){
          if((indexn(i)==1)&&(indexp(j)==1)){
            xtilde(i,j) = x(i,j);
          }
        }

      }
      outputs = TukeyM(y, xtilde, betahat, sigmahat);
      betahat = as<arma::vec>(outputs["betahat"]);

    }
  }


  if(tech == "cell"){
    double resori, restilde, xdiff, diff;
    for(k = 1; k < maxiter; k++){
      xtildec = join_rows(onevec, xtilde);
      for(i=0; (i<n); i++){
        restilde =   as_scalar(y[i] - xtildec.row(i)*betahat) ;
        for(j=0; (j<p); j++){
          xdiff = as_scalar(xtildec(i,(j+1))- xc(i,(j+1)));
          if((xdiff!=0)&&(as_scalar(betahat(j+1))!=0)){
            diff = xdiff*as_scalar(betahat(j+1));
          }else{
            diff = 0;
          }
          resori = restilde + diff;
          //change it to restilde
          //如果没有imputation可以跳过计算
          if(abs(resori) < abs(restilde)){
            xtilde(i,j) = x(i,j);
          }
        }
      }
      outputs = TukeyM(y, xtilde, betahat, sigmahat);
      arma::vec betanew = as<arma::vec>(outputs["betahat"]);
      if((abs(betanew-betahat)).max() < pow(10, -4)) {break;};
      betahat = betanew;
    }
  }

  double loss = as<double>(outputs["loss"]);
  double L = as<double>(outputs["L"]);
  double t = as<double>(outputs["t"]);
  arma::mat weightmat = as<arma::mat>(outputs["weightmat"]);

  List results = List::create(
    Named("betahat") = betahat,
    Named("sigmahat") = sigmahat,
    Named("xtilde") = xtilde,
    Named("t") = t,
    Named("L") = L,
    Named("k") = k,
    Named("weightmat") = weightmat,
    Named("loss") = loss
  );
  return(results);

}


