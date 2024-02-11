require("cellWise")
require("robustbase")



lambdamax_beta = function(y, x, betahat, intercept,
                          deltahat, zetahat,
                          softbeta,
                          lambda_delta, softdelta,
                          lambda_zeta, softzeta,
                          alpha, tol, maxiter){

  n = dim(x)[1]
  p = dim(x)[2]

  lambdamax = 2*max(t(x - deltahat)%*%(y-intercept - zetahat))


  iovec = rep(5,2)
  ind = 0
  k = 0
  while((min(iovec)<10)|(ind==FALSE)){
    outputs = reg_beta_delta(y = y, x = x, betahat = betahat, intercept = intercept,
                             deltahat = deltahat, zetahat = zetahat,
                             lambda_beta = lambdamax, softbeta = softbeta,
                             lambda_delta = lambda_delta, softdelta = softdelta,
                             lambda_zeta = lambda_zeta, softzeta = softzeta,
                             alpha = alpha, tol = tol, maxiter = maxiter)

    betaget = outputs$betahat
    ind = (sum(abs(betaget))==0)
    lambdamax = ((1/2)^((ind-0.5)*2*1/(2^min(iovec))))*lambdamax
    iovec[ind+1] = iovec[ind+1]+1
    k = k+1
  }
  result = list(iovec = iovec, lambdamax = lambdamax, k = k)
  return(result)

}




#' robust sparse regression under cellwise contamination (with a grid of lambdas)
#'
#' @param y response
#' @param x design matrix
#' @param betahat initial estimate of beta
#' @param intercept initial estimate of intercept
#' @param softbeta whether to use soft/hard threshold for beta
#' @param softdelta whether to use soft/hard threshold for delta(outliers in x)
#' @param softzeta whether to use soft/hard threshold for zeta(outliers in y)
#' @param lambda_delta tuning parameter of delta
#' @param lambda_zeta tuning parameter of zeta
#' @param alpha the importance factor of the regression loss (between 0-1, by default is 0.5)
#' @param penal the penalty parameter for model selection (by default is 1, equivalent to BIC )
#' @param penaldelta the penalty of number of detected outliers (for debug, by default is 0)
#' @param tol the tolerance of convergence
#' @param maxiter number of iterations
#'
#' @return
#' betahat: the estimated beta
#'
#' intercept_hat: the estiamted intercept
#'
#' betahat_opt: the estimated beta with post-cellwise-robust regression
#'
#' intercept_opt: the estimated intercept with post-cellwise-robust regression

#' @export
#'
#' @examples
#'
#' data = genevar()
#' y = data$y
#' x = data$x
#' fit = sregcell(y,x)
#'
sregcell = function(y,x, betahat = NULL, intercept = NULL,
                    softbeta = TRUE, softdelta = TRUE, softzeta = TRUE,
                    lambda_delta = 2.56, lambda_zeta = 2.56, alpha = 0.5,
                    penal = 1, penaldelta = 1,
                    tol = 1e-3, maxiter = 100){

  # {
  #   betahat = NULL
  #   intercept = NULL
  #
  #   softbeta = TRUE
  #   softdelta = TRUE
  #   softzeta = TRUE
  #
  #   lambda_delta = 2.56
  #   lambda_zeta = 2.56
  #   alpha = 0.5
  #   penal = 1
  #   penaldelta = 1
  #
  #   maxiter = 100
  # }

  n = dim(x)[1]
  p = dim(x)[2]


  deltahat = threshold_mat(x,matrix(lambda_delta, n,p))

  xc = x - deltahat
  zetahat = rep(0,n)

  if(is.null(betahat)){
    betahat = rep(0,p)
    intercept = 0
  }

  { #lambda grid
    length = 50
    # grid = get_lambda_grid(y,xc, betahat,intercept, softbeta, length = length)
    lambdamax = lambdamax_beta(y = y, x = x, betahat = betahat, intercept = intercept,
                               deltahat = deltahat, zetahat = zetahat,
                               softbeta = softbeta,
                               lambda_delta = lambda_delta, softdelta = softdelta,
                               lambda_zeta = lambda_zeta, softzeta = softzeta,
                               alpha = alpha, tol = tol, maxiter = maxiter)$lambdamax

    lmin = lambdamax/10^3
    grid = c(exp(seq(log(lambdamax),log(lmin),length = length)))
    if(p < n){ grid = c(grid,0)}
  }

  allfits = lapply(grid, function(lambda){reg_beta_delta(y = y, x = x, betahat = betahat, intercept = intercept,
                                                         deltahat = deltahat, zetahat = zetahat,
                                                         lambda_beta = lambda, softbeta = softbeta,
                                                         lambda_delta = lambda_delta, softdelta = softdelta,
                                                         lambda_zeta = lambda_zeta, softzeta = softzeta,
                                                         alpha = alpha, tol = tol, maxiter = maxiter)})


  regloss = unlist(lapply(allfits, function(fit) fit$regloss))
  scaleloss = unlist(lapply(allfits, function(fit) fit$scaleloss))
  penaltylossx = unlist(lapply(allfits, function(fit) fit$penaltylossx))
  penaltylossy = unlist(lapply(allfits, function(fit) fit$penaltylossy))
  activebeta = unlist(lapply(allfits, function(fit) sum(as.logical(fit$betahat))))
  activedelta = unlist(lapply(allfits, function(fit)  sum(as.logical(fit$deltahat))))
  activedeltabeta = unlist(lapply(allfits, function(fit)  sum((fit$deltahat!=0)[,as.logical(fit$betahat)])))
  betamat = as.data.frame(t(as.data.frame(lapply(allfits, function(fit) fit$betahat))))
  interceptvec = unlist(lapply(allfits, function(fit) fit$intercept))
  colactivedelta = unlist(lapply(allfits, function(fit)  max(colSums(fit$deltahat!=0)*(fit$betahat!=0))))


  ic = regloss + 2*penaltylossy + penal*log(n)*activebeta + penaldelta*activedeltabeta
  icadj = ic
  icadj[colactivedelta>(length(y)*0.3)] = NA


  label = which.min(icadj)
  result_opt = allfits[[label]]
  betahat = result_opt$betahat
  intercept_hat = result_opt$intercept


  { # post regression
    betahat_post = rep(0,p)
    intercept_hat_post = intercept_hat

    labels = which(betahat!=0)
    if(length(labels)!=0){
      xpost = t(t(x[,labels]))
      deltahat_post = t(t(deltahat[,labels]))

      fit_post = reg_beta_delta(y = y, x = xpost, betahat = betahat[labels], intercept = intercept,
                                deltahat = deltahat_post, zetahat = zetahat,
                                lambda_beta = 0, softbeta = TRUE,
                                lambda_delta = lambda_delta, softdelta = softdelta,
                                lambda_zeta = lambda_zeta, softzeta = softzeta,
                                alpha = alpha, tol = tol, maxiter = maxiter)

      betahat_post[labels] = fit_post$betahat
      intercept_hat_post = fit_post$intercept
    }

  }


  return(list(allfits = allfits, result_opt = result_opt,
              grid = grid, activebeta = activebeta, ic = ic,
              colactivedelta = colactivedelta, betamat = betamat, label = label,
              k = result_opt$k, k_post = fit_post$k,colactivedelta = colactivedelta,
              betahat = betahat, intercept_hat = intercept_hat,
              betahat_post = betahat_post, intercept_hat_post = intercept_hat_post))
}


#' Cellwise regularized robust sparse regression (with a specific lambda)
#'
#' @param y response
#' @param x design matrix
#' @param softbeta whether to use soft/hard threshold for beta
#' @param softdelta whether to use soft/hard threshold for delta(outliers in x)
#' @param softzeta whether to use soft/hard threshold for zeta(outliers in y)
#' @param lambda_delta tuning parameter of delta
#' @param lambda_zeta tuning parameter of zeta
#' @param lambda tuning parameter of beta
#' @param alpha the importance factor of the regression loss (between 0-1, by default is 0.5)
#' @param tol the tolerence of convergence
#' @param maxiter the number of interations
#'
#' @return
#' fit
#' @export
#'
#' @examples
#'
#'
#' data = genevar()
#' y = data$y
#' x = data$x
#' fit = sregcell_lambda(y,x, lambda = 1)
#'
sregcell_lambda = function(y,x, softbeta = TRUE, softdelta = TRUE, softzeta = TRUE,
                           lambda_delta = 2.56, lambda_zeta = 2.56, lambda = 0, alpha = 0.5,
                           tol = 1e-3, maxiter = 100){

  n = dim(x)[1]
  p = dim(x)[2]

  delta = threshold_mat(x,matrix(lambda_delta, n,p))
  xc = x - delta
  intercept = 0
  zetahat = rep(0,n)
  betahat = rep(0,p)


  result = reg_beta_delta(y = y, x = x, betahat = betahat, intercept = intercept,
                          deltahat = delta, zetahat = zetahat,
                          lambda_beta = lambda, softbeta = softbeta,
                          lambda_delta = lambda_delta, softdelta = softdelta,
                          lambda_zeta = lambda_zeta, softzeta = softzeta,
                          alpha = alpha, tol = tol, maxiter = maxiter)
  return(result)
}






### STD
quantscale = function(xvec, df = Inf, tau = 0.95){
  (quantile(xvec, tau) - quantile(xvec, 1-tau))/(qt(tau, df = df) - qt(1-tau, df = df))
}

qnscale = function(xvec, df = Inf){
  d = 1/(sqrt(2)*qt(5/8, df = df))
  robustbase::Qn(xvec, constant = d)
}

robstd = function(x, centerf = median, scalef = qnscale, df = df){
  center  =  apply(x, 2, centerf)
  x = sweep(x, 2, center, check.margin = FALSE)
  scale = apply(x, 2, scalef, df = df)
  x = sweep(x, 2, scale, "/", check.margin = FALSE)
  return(list(center = center, scale = scale, x = x))
}


#' robust sparse regression under cellwise contamination with a grid of lambdas (standardize predictors first)
#'
#' @param y response
#' @param x design matrix
#' @param scale.method method we used to obtain robust scales, by default is qn
#' @param df degrees of freedom of the assumed distribution
#' @param softbeta whether to use soft/hard threshold for beta
#' @param softdelta whether to use soft/hard threshold for delta(outliers in x)
#' @param softzeta whether to use soft/hard threshold for zeta(outliers in y)
#' @param lambda_delta tuning parameter of delta
#' @param lambda_zeta tuning parameter of zeta
#' @param prob probability of quantiles, by default is 0.995
#' @param the importance factor of the regression loss (between 0-1, by default is 0.5)
#' @param penal the penalty parameter for model selection (by default is 1, equivalent to BIC )
#' @param penaldelta the penalty of number of detected outliers (for debug, by default is 0)
#' @param tol tolerence of convergence
#' @param maxiter the number of interations
#'
#' @return
#' betahat: the estimated beta
#'
#' intercept_hat: the estiamted intercept
#'
#' betahat_opt: the estimated beta with post-cellwise-robust regression
#'
#' intercept_opt: the estimated intercept with post-cellwise-robust regression
#' @export
#'
#' @examples
#' data = genevar()
#' y = data$y
#' x = data$x
#' fit = sregcell_std(y,x)
sregcell_std = function(y,x,
                        scale.method = qnscale,  df = Inf,
                        softbeta = TRUE, softdelta = TRUE, softzeta = TRUE,
                        lambda_delta = NULL, lambda_zeta = 1, prob = 0.995,
                        alpha = 0.5,
                        penal = 1, penaldelta = 0,
                        tol = 1e-3, maxiter = 30){


  x = as.matrix(x)
  n = dim(x)[1]
  p = dim(x)[2]

  if(is.null(lambda_delta)){
    lambda_delta = qt(p = prob, df = df)
  }

  { # standardization
    stdrst = robstd(x, centerf = median, scalef = scale.method, df = df)
    mux = stdrst$center
    scalex = stdrst$scale
    xstd = stdrst$x
  }

  { set.seed(2)
    fit0 = regcell::Rlars(y,xstd)
    sigmahat = fit0$sigmahat # tend to overestimate

    yscl = (y-median(y))/sigmahat
    intercept = (fit0$betahat[1] - median(y))/sigmahat
    betahat = fit0$betahat[-1]/sigmahat

    # fit00 = regcell::Rlars(yscl,xstd)
    # intercept = fit00$betahat[1]*sigmahat + median(y)
    # betahat = fit00$betahat[-1]*sigmahat
  }

  rst = sregcell(y = yscl, x = xstd, betahat = betahat, intercept = intercept,
                 softbeta = softbeta, softdelta = softdelta, softzeta = softzeta,
                 lambda_delta = lambda_delta, lambda_zeta = lambda_zeta, alpha = alpha,
                 penal = penal, penaldelta = penaldelta,
                 tol = tol, maxiter = maxiter)

  #de-standardization
  betahat = as.numeric((1/scalex)*rst$betahat*sigmahat)
  intercept_hat = (rst$intercept_hat - sum((1/scalex)*mux*rst$betahat))*sigmahat + median(y)
  label = rst$label

  betahat_post = (1/scalex)*rst$betahat_post*sigmahat
  intercept_hat_post = (rst$intercept_hat_post - sum((1/scalex)*mux*rst$betahat_post))*sigmahat + median(y)

  return(list(betahat = betahat, intercept_hat = intercept_hat, label = label,
              betahat_post = betahat_post, intercept_hat_post = intercept_hat_post,
              k = rst$k, k_post = rst$k_post, colactivedelta = rst$colactivedelta,
              sigmahat = sigmahat, mux = mux, scalex = scalex))
}



