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
arma::vec soft_threshold_beta(arma::vec betahat, arma::vec alambda_beta){   //should be a vector
  arma::vec diff = abs(betahat) - alambda_beta;
  arma::vec value = sign(betahat)%(diff>0)%diff;
  return value;
}


// [[Rcpp::export]]
arma::mat soft_threshold_delta(arma::mat deltahat, arma::mat alambda_delta){   //should be a vector
  arma::mat diff = abs(deltahat) - alambda_delta;
  arma::mat value = sign(deltahat)%(diff>0)%diff;
  return value;
}



// [[Rcpp::export]]
List reg_beta(arma::vec y, arma::mat x, arma::vec betahat, double intercept, double lambda_beta, bool adabeta, double maxiter = 30){

  //initialize variables
  unsigned int n = x.n_rows; //sample size
  unsigned int p = x.n_cols; //dimension size

  arma::vec lambdavec_beta(p);//a vector filled with lambda, just for computation convenience
  lambdavec_beta.fill(lambda_beta);

  arma::vec gammavec(p);//a vector filled with gamma, just for computation convenience
  gammavec.fill(0.05);
  arma::vec twoerpi(p);// a constant vector filled with 1/pi
  twoerpi.fill(2/M_PI);

  arma::vec alambdavec_beta = lambdavec_beta;//-----------------------whether using adaptive penalties
  // compute adaptive penaltive for each variables//-----------whether penalties should be related with delta

  arma::mat weightmat = eye(n,n);
  double L = arma::eig_sym(trans(x)*weightmat*x).max();
  arma::vec steptvec(p);// a constant vector filled with 1/pi
  steptvec.fill(1/L);

  arma::vec interceptvec(n);
  arma::vec res(n);


  arma::vec betapos(p);
  arma::vec betanew(p);
  double loss;

  double k = 1;//iterator
  while(k <= maxiter)//start iteration
  {
    interceptvec.fill(intercept);
    res = y - interceptvec - x*betahat;//residual
    betapos = betahat + steptvec%(trans(x) * res);

    if(adabeta){
      alambdavec_beta = gammavec%(gammavec + twoerpi)/(square(gammavec) +square(betahat))%lambdavec_beta;
    }

    betanew = soft_threshold_beta(betapos, steptvec%alambdavec_beta);
    if((abs(betanew-betahat)).max() < pow(10, -4)) {break;};
    betahat = betanew;
    res = y - x*betahat;//residuals
    intercept  = mean(res);
    k++;
  }
  loss = accu(pow(y-x*betahat,2));//eps and x have been std


  List results = List::create(
    Named("betahat") = betahat,
    Named("interceptvec") = interceptvec,
    Named("loss") = loss
    );
  return(results);

}





// [[Rcpp::export]]
List reg_delta(arma::vec y, arma::mat x, arma::vec betahat, arma::mat deltahat,
               double lambda_delta, arma::mat cellweight, double alpha, double maxiter = 30){

  //initialize variables
  unsigned int n = x.n_rows; //sample size
  unsigned int p = x.n_cols; //dimension size

  arma::mat lambdamat_delta(n,p);//a vector filled with lambda,
  lambdamat_delta.fill(lambda_delta);
  arma::mat alphamat(p,p);
  alphamat.fill(alpha);

  arma::mat one_minus_alphahat(p,p);
  one_minus_alphahat.fill((1-alpha));

  arma::mat alphamatnp(n,p);
  alphamatnp.fill(alpha);

  arma::mat alambdamat_delta = lambdamat_delta%cellweight;//------whether using adaptive penalties

  arma::mat Sigmahat = eye(p,p);
  arma::mat xtx = alphamat%(betahat*trans(betahat))+ one_minus_alphahat%(Sigmahat);
  arma::mat ytx = x*xtx - alphamatnp%(y*trans(betahat));

  double L = arma::eig_sym(xtx).max();
  arma::mat stepmat(n,p);
  stepmat.fill(1/L);
  arma::mat deltapos(n,p);
  arma::mat deltanew(n,p);
  double k = 1;//iterator
  while(k <= maxiter)//start iteration
  {
    deltapos = deltahat + stepmat%(ytx - deltahat*xtx);
    deltanew = soft_threshold_delta(deltapos, stepmat%alambdamat_delta);//------------------------------!!!!!!!!!!!
    if((abs(deltanew-deltahat)).max() < pow(10, -4)) {break;};
    deltahat = deltanew;
    k++;
  }

  List results = List::create(
    Named("deltahat") = deltahat
  );
  return(results);


}


// [[Rcpp::export]]
List reg_beta_delta(arma::vec y, arma::mat x, arma::vec betahat, double intercept, arma::mat deltahat, arma::mat cellweight,
                    double lambda_beta, bool adabeta, double lambda_delta,double alpha, double maxiter = 30){

  //initialize variables
  unsigned int n = x.n_rows; //sample size
  unsigned int p = x.n_cols; //dimension size

  List outputs_beta;
  List outputs_delta;
  arma::vec betaget(p);
  arma::vec interceptvec(n);
  interceptvec.fill(intercept);
  arma::mat deltaget(n,p);
  double loss;
  double k = 1;//iterator
  while(k <= maxiter)//start iteration
  {
    outputs_delta = reg_delta(y-interceptvec, x, betahat, deltahat, lambda_delta, cellweight, alpha);
    deltaget = as<arma::mat>(outputs_delta["deltahat"]);
    outputs_beta = reg_beta(y, x-deltahat, betahat, intercept, lambda_beta, adabeta);
    betaget = as<arma::vec>(outputs_beta["betahat"]);
    interceptvec = as<arma::vec>(outputs_beta["interceptvec"]);

    if((abs(betaget-betahat)).max() < pow(10, -4)) {break;};
    betahat = betaget;
    deltahat = deltaget;
    k++;
  }

  arma::mat lambdamat_delta(n,p);
  lambdamat_delta.fill(lambda_delta);
  loss = alpha * as<double>(outputs_beta["loss"]) + (1 - alpha)*accu(pow(x-deltahat,2)) + accu(lambdamat_delta%abs(deltahat));
  intercept= as_scalar(interceptvec[0]);


  List results = List::create(
    Named("intercept") = intercept,
    Named("betahat") = betahat,
    Named("deltahat") = deltahat,
    Named("loss") = loss
  );
  return(results);


}


// [[Rcpp::export]]

List lambdamax_beta(arma::vec y, arma::mat x, arma::vec betahat, double intercept, bool adabeta, double lambdamax = 1){
  //initialize variables
  double n = x.n_rows; //sample size
  double p = x.n_cols; //dimension size
  arma::vec betaget(p);
  arma::vec iovec(2, fill::zeros);
  int ind = 0;
  bool label = true;
  int k = 0;
  while((label)){
    List outputs;
    outputs = reg_beta(y, x, betahat,intercept,lambdamax,adabeta);
    betaget = as<arma::vec>(outputs["betahat"]);
    ind = accu(abs(betaget))==0;// 1 means zero model; 0 means others
    label = iovec.min()<=5;
    lambdamax = (pow(0.5, (ind-0.5)*2*1/pow(2,iovec.min()))*lambdamax);
    iovec[ind] = iovec[ind]+1;
    k++;

  }

  List results = List::create(
    Named("iovec") = iovec,
    Named("ind") = ind,
    Named("label") = label,
    Named("lambdamax") = lambdamax,
    Named("betaget") = betaget,
    Named("k") = k
  );
  return(results);

}


