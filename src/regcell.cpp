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
arma::vec threshold_beta(arma::vec betahat, arma::vec alambda_beta, bool softbeta){   //should be a vector
  double p = betahat.n_rows;
  arma::vec diff = abs(betahat) - alambda_beta;
  arma::vec value(p);
  if(softbeta){
    value = sign(betahat)%(diff>0)%diff;
  }else{
    value = betahat%(diff>0);
  }
  return value;
}


// [[Rcpp::export]]
arma::mat threshold_delta(arma::mat deltahat, arma::mat alambda_delta, bool softdelta){   //should be a vector
  double n = deltahat.n_rows;
  double p = deltahat.n_cols;
  arma::mat value(n,p);
  arma::mat diff = abs(deltahat) - alambda_delta;
  if(softdelta){
    value = sign(deltahat)%(diff>0)%diff;
  }else{
    value = deltahat%(diff>0);
  }
  return value;
}



// [[Rcpp::export]]
List reg_beta(arma::vec y, arma::mat x, arma::vec betahat, double intercept,
              arma::vec lambdavec_beta, bool softbeta = true, double maxiterbeta = 100){

  //initialize variables
  unsigned int n = x.n_rows; //sample size
  unsigned int p = x.n_cols; //dimension size

  arma::vec alambdavec_beta = lambdavec_beta;//------------whether using adaptive penalties

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
  while(k <= maxiterbeta)//start iteration
  {
    interceptvec.fill(intercept);
    res = y - interceptvec - x*betahat;//residual
    betapos = betahat + steptvec%(trans(x) * res);

    betanew = threshold_beta(betapos, steptvec%alambdavec_beta, softbeta = softbeta);
    if((abs(betanew-betahat)).max() < pow(10, -6)) {break;};
    betahat = betanew;
    //res = y - x*betahat;//residuals
    intercept  = intercept + mean(res);
    k++;
  };
  loss = accu(pow(y-x*betahat,2));//eps and x have been std


  List results = List::create(
    Named("betahat") = betahat,
    Named("intercept") = intercept,
    Named("loss") = loss
    );
  return(results);

}





// [[Rcpp::export]]
List reg_delta(arma::vec y, arma::mat x, arma::vec betahat, arma::mat deltahat,
               double lambda_delta, arma::mat cellweight, double alpha, bool softdelta = true, double maxiterdelta = 100){

  //initialize variables
  unsigned int n = x.n_rows; //sample size
  unsigned int p = x.n_cols; //dimension size

  arma::mat lambdamat_delta(n,p);//a vector filled with lambda,
  lambdamat_delta.fill(lambda_delta);
  arma::mat alphamatnp(n,p);
  alphamatnp.fill(alpha);

  arma::mat one_minus_alphamatnp(n,p);
  one_minus_alphamatnp.fill(1-alpha);

  arma::mat alambdamat_delta = lambdamat_delta%cellweight;//------whether using adaptive penalties

  double L = alpha*as_scalar(trans(betahat)*betahat) + (1-alpha);
  //stepwise should be related with the second order derivative
  // only unit matrix could compute like this
  arma::mat stepmat(n,p);
  stepmat.fill(0.9/L);

  arma::mat xclean(n,p);
  arma::mat gradient(n,p);
  arma::mat deltapos(n,p);
  arma::mat deltanew(n,p);
  arma::vec resclean(n);

  double k = 1;//iterator
  while(k <= maxiterdelta)//start iteration
  {

    xclean = x - deltahat;
    resclean = y - xclean*betahat;
    gradient = alphamatnp%(resclean*trans(betahat)) - (one_minus_alphamatnp%xclean);
    deltapos = deltahat - stepmat%gradient;
    deltanew = threshold_delta(deltapos, stepmat%alambdamat_delta, softdelta = softdelta);//------------------------------!!!!!!!!!!!c
    if((abs(deltanew-deltahat)).max() < pow(10, -6)) {break;};
    deltahat = deltanew;
    k++;
  };

  List results = List::create(
    Named("deltahat") = deltahat,
    Named("xclean") = xclean,
    Named("resclean") = resclean,
    Named("gradient") = gradient,
    Named("deltapos") = deltapos,
    Named("alambdamat_delta") = alambdamat_delta
  );
  return(results);


}


// [[Rcpp::export]]
List reg_beta_delta(arma::vec y, arma::mat x,
                    arma::vec betahat, double intercept, arma::mat deltahat,
                    arma::vec rowweight, arma::mat colweight,
                    double lambda_beta, bool adabeta, bool softbeta,
                    double lambda_delta,bool adadelta, bool softdelta,
                    double alpha, double maxiter = 20){

  //initialize variables
  unsigned int n = x.n_rows; //sample size
  unsigned int p = x.n_cols; //dimension size

  List outputs_beta;
  List outputs_delta;

  arma::vec betaget(p);
  arma::vec interceptvec(n);
  arma::mat deltaget(n,p);
  if(!adadelta){
    rowweight.fill(1);
    colweight.fill(1);
  };

  arma::mat alphamatnp(n,p);
  alphamatnp.fill(alpha);

  arma::mat one_minus_alphahatnp(n,p);
  one_minus_alphahatnp.fill((1-alpha));

  arma::mat cellweight;
  cellweight = rowweight*trans(ones(p))%alphamatnp + colweight%one_minus_alphahatnp;

  arma::vec lambdavec_beta(p);
  lambdavec_beta.fill(lambda_beta);

  arma::vec lambdavec_betanew = lambdavec_beta;

  arma::vec lambdavec_delta(p);
  lambdavec_delta.fill(lambda_delta);

  double regloss;
  double scaleloss;
  double penaltyloss;
  double k = 1;//iterator
  while(k <= maxiter)//start iteration
  {
    if(adabeta){//adabeta only useful when not all important variables are contaminated
      cellweight = alphamatnp%(rowweight*trans(betahat)) + one_minus_alphahatnp%colweight;
      lambdavec_betanew = lambdavec_delta%trans(sum(abs(deltahat))) + lambdavec_beta;
    };

    interceptvec.fill(intercept);
    // reg_delta(arma::vec y, arma::mat x, arma::vec betahat, arma::mat deltahat,
    // double lambda_delta, arma::mat cellweight, double alpha
    outputs_delta= reg_delta(y = y - interceptvec, x = x, betahat = betahat, deltahat = deltahat,
                             lambda_delta = lambda_delta, cellweight = cellweight, alpha = alpha, softdelta = softdelta);

    deltaget = as<arma::mat>(outputs_delta["deltahat"]);
    //arma::vec y, arma::mat x, arma::vec betahat, double intercept, arma::vec lambdavec_beta, double maxiter = 20
    outputs_beta = reg_beta(y = y, x = x-deltaget, betahat = betahat, intercept = intercept,
                            lambdavec_beta = lambdavec_betanew, softbeta = softbeta);//找bug找了两天，我吐了
    //但是即使没有及时更新也应该收敛才对吧
    betaget = as<arma::vec>(outputs_beta["betahat"]);
    intercept = double(outputs_beta["intercept"]);

    if((abs(betaget-betahat)).max() < pow(10, -6)) {break;};
    betahat = betaget;
    deltahat = deltaget;
    k++;
  };

  arma::mat lambdamat_delta(n,p);
  lambdamat_delta.fill(lambda_delta);
  regloss = as<double>(outputs_beta["loss"]);
  scaleloss = accu(pow(x-deltahat,2));
  penaltyloss = accu(lambdamat_delta%abs(deltahat));
  //-------------------------------------------------------------有大问题

  List results = List::create(
    Named("intercept") = intercept,
    Named("betahat") = betahat,
    Named("deltahat") = deltahat,
    Named("rowweight") = rowweight,
    Named("colweight") = colweight,
    Named("cellweight") = cellweight,
    Named("lambda_beta") = lambda_beta,
    Named("lambda_delta") = lambda_delta,
    Named("k") = k,
    Named("regloss") = regloss,
    Named("scaleloss") = scaleloss,
    Named("penaltyloss") = penaltyloss
  );
  return(results);


}


