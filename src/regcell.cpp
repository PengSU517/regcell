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
List rob_pca(arma::mat x, arma::mat xc, arma::mat delta, double lambda, int maxiter = 100){
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




// [[Rcpp::export]]
List reg_beta(arma::vec y, arma::mat x, arma::vec betahat, double intercept,
              arma::vec lambdavec_beta, bool softbeta = true, double maxiterbeta = 2){

  //initialize variables
  unsigned int n = x.n_rows; //sample size
  unsigned int p = x.n_cols; //dimension size

  arma::vec alambdavec_beta = lambdavec_beta;//------------whether using adaptive penalties

  //arma::mat weightmat = eye(n,n);
  double L = arma::eig_sym(trans(x)*x).max();//直接用奇异值分解不好吗
  arma::vec steptvec(p);// a constant vector filled with 1/pi
  steptvec.fill(0.9/L);

  arma::vec interceptvec(n);
  arma::vec res(n);


  arma::vec betapos(p);
  arma::vec betanew(p);


  double k = 1;//iterator
  while(k <= maxiterbeta)//start iteration
  {
    interceptvec.fill(intercept);
    res = y - interceptvec - x*betahat;//residual
    betapos = betahat + steptvec%(trans(x) * res);

    betanew = threshold_vec(betapos, steptvec%alambdavec_beta, softbeta);
    if((abs(betanew-betahat)).max() < pow(10, -6)) {break;};
    betahat = betanew;
    // 收敛速度太慢
    //res = y - x*betahat;//residuals
    intercept  = intercept + mean(res);
    k++;
  };

  List results = List::create(
    Named("betahat") = betahat,
    Named("intercept") = intercept
    );
  return(results);

}





// [[Rcpp::export]]
List reg_delta(arma::vec y, arma::mat x, arma::vec betahat, arma::mat deltahat,
               double lambda_delta, arma::mat cellweight, double alpha, bool softdelta = true, double maxiterdelta = 2){

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
  //在pca时步长如何选择是个问题
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
    gradient = alphamatnp%(resclean*trans(betahat)) - (one_minus_alphamatnp%xclean);//这里是不是有问题
    deltapos = deltahat - stepmat%gradient;
    deltanew = threshold_mat(deltapos, stepmat%alambdamat_delta, softdelta);
    //------------------------------!!!!!!!!!!!关于加罚的选择
    if((abs(deltanew-deltahat)).max() < pow(10, -6)) {break;};
    deltahat = deltanew;
    k++;
  };

  List results = List::create(
    Named("deltahat") = deltahat,
    Named("xclean") = xclean,
    Named("resclean") = resclean,
    Named("gradient") = gradient,
    Named("stepmat") = stepmat,
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
                    double lambda_delta,bool adadelta, bool adarow, bool softdelta,
                    double alpha, double maxiter = 20){

  //initialize variables
  unsigned int n = x.n_rows; //sample size
  unsigned int p = x.n_cols; //dimension size

  List outputs_beta;
  List outputs_delta;

  arma::vec interceptvec(n);
  interceptvec.fill(intercept);

  arma::vec ycenter = y - interceptvec;
  arma::mat xclean = x - deltahat;

  arma::vec betaget(p);
  arma::mat deltaget(n,p);

  if(!adadelta){
    colweight.fill(1);
  };

  if(!adarow){
    rowweight.fill(1);
  };

  arma::mat alphamatnp(n,p);
  alphamatnp.fill(alpha);

  arma::mat one_minus_alphahatnp(n,p);
  one_minus_alphahatnp.fill((1-alpha));

  arma::mat cellweight;
  //cellweight = rowweight*trans(ones(p))%alphamatnp + colweight%one_minus_alphahatnp;
  cellweight = colweight;

  arma::vec lambdavec_beta(p);
  lambdavec_beta.fill(lambda_beta);

  arma::vec lambdavec_betanew = lambdavec_beta;

  arma::vec lambdavec_delta(p);
  lambdavec_delta.fill(lambda_delta);

  double k = 1;//iterator
  while(k <= maxiter)//start iteration
  {
    if(adabeta){//adabeta only useful when not all important variables are contaminated
      cellweight = alphamatnp%(rowweight*trans(betahat)) + one_minus_alphahatnp%colweight;
      lambdavec_betanew = lambdavec_delta%trans(sum(abs(deltahat))) + lambdavec_beta;
    };

    // reg_delta(arma::vec y, arma::mat x, arma::vec betahat, arma::mat deltahat,
    // double lambda_delta, arma::mat cellweight, double alpha
    outputs_delta= reg_delta(ycenter, x, betahat, deltahat,lambda_delta, cellweight, alpha, softdelta);

    //调用函数的时候不能申明参数位置？必须按顺序吗？
    deltaget = as<arma::mat>(outputs_delta["deltahat"]);
    //arma::vec y, arma::mat x, arma::vec betahat, double intercept, arma::vec lambdavec_beta, double maxiter = 20
    outputs_beta = reg_beta(y, xclean, betahat, intercept, lambdavec_betanew, softbeta);//找bug找了两天，我吐了
    //但是即使没有及时更新也应该收敛才对吧
    betaget = as<arma::vec>(outputs_beta["betahat"]);

    if((abs(betaget-betahat)).max() < pow(10, -6)) {break;};

    betahat = betaget;
    deltahat = deltaget;
    intercept = double(outputs_beta["intercept"]);
    interceptvec.fill(intercept);
    ycenter = y - interceptvec;
    xclean = x - deltahat;
    k++;

  };

  arma::mat lambdamat_delta(n,p);
  lambdamat_delta.fill(lambda_delta);
//
//   arma::mat loss1 = x;
//   arma::mat loss2 = deltahat;
//   arma::mat loss3 = (x-deltahat);
//   arma::mat loss4 = pow((x-deltahat),2);

  double regloss = accu(pow(y - interceptvec - (x - deltahat)*betahat,2));
  double scaleloss = accu(pow(x-deltahat,2));
  double penaltyloss = accu(lambdamat_delta%abs(deltahat));
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


