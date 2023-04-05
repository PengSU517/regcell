require("cellWise")
require("purrr")
require("robustbase")
require("robcovsel")

lambdamax_beta = function(y, x, betahat, intercept,
                          deltahat, zetahat,
                          softbeta,
                          lambda_delta, softdelta,
                          lambda_zeta, softzeta,
                          alpha, maxiter){

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
                             alpha = alpha, maxiter = maxiter)

    betaget = outputs$betahat
    ind = (sum(abs(betaget))==0)
    lambdamax = ((1/2)^((ind-0.5)*2*1/(2^min(iovec))))*lambdamax
    iovec[ind+1] = iovec[ind+1]+1
    k = k+1
  }
  result = list(iovec = iovec, lambdamax = lambdamax, k = k)
  return(result)

}

# lambdamax2_beta = function(y,x,betahat,intercept, softbeta){
#
#   n = dim(x)[1]
#   p = dim(x)[2]
#
#   lambdamax = 20*max(t(x)%*%(y-mean(y)))
#   outputs = reg_beta(y = y, x = x, betahat = betahat, intercept = intercept,
#                           alambdavec_beta = rep(lambdamax,p), softbeta = softbeta, maxiterbeta = 100)
#   lambdamax = max(abs(outputs$mgradient))*1.01
#   return(list(lambdamax = lambdamax))
#
# }


lambdamax3_beta = function(y, x, betahat, intercept,
                           deltahat, zetahat,
                           softbeta,
                           lambda_delta, softdelta,
                           lambda_zeta, softzeta,
                           alpha, maxiter){

  n = dim(x)[1]
  p = dim(x)[2]

  lambdamax = 2*max(t(x - deltahat)%*%(y-intercept - zetahat))
  outputs = reg_beta_delta(y = y, x = x, betahat = betahat, intercept = intercept,
                           deltahat = deltahat, zetahat = zetahat,
                           lambda_beta = lambdamax, softbeta = softbeta,
                           lambda_delta = lambda_delta, softdelta = softdelta,
                           lambda_zeta = lambda_zeta, softzeta = softzeta,
                           alpha = alpha, maxiter = maxiter)
  #outputs$mgradient
  lambdamax = max(abs(outputs$mgradient))*1.01
  return(list(lambdamax = lambdamax))
}


## sparse robust regression with imputed design matrix
sregcell = function(y,x, betahat = NULL, intercept = NULL,
                    softbeta = TRUE, softdelta = TRUE, softzeta = TRUE,
                    lambda_delta = 2.56, lambda_zeta = 2.56, alpha = 0.5,
                    penal = 1, penaldelta = 1,
                    maxiter = 100){


  n = dim(x)[1]
  p = dim(x)[2]


  deltahat = threshold_mat(x,matrix(lambda_delta, n,p))##这个没问题吧
  # 根据deltahat设置adaptive lambda_deltavec
  colSums(deltahat!=0)
  colSums(data$outlier!=0)

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
                                alpha = alpha, maxiter = maxiter)$lambdamax

    lmin = lambdamax/10^2
    grid = c(exp(seq(log(lambdamax),log(lmin),length = length)))
    if(p < n){ grid = c(grid,0)}
  }

  allfits = lapply(grid, function(lambda){reg_beta_delta(y = y, x = x, betahat = betahat, intercept = intercept,
                                                      deltahat = deltahat, zetahat = zetahat,
                                                      lambda_beta = lambda, softbeta = softbeta,
                                                      lambda_delta = lambda_delta, softdelta = softdelta,
                                                      lambda_zeta = lambda_zeta, softzeta = softzeta,
                                                      alpha = alpha, maxiter = maxiter)})


  regloss = unlist(lapply(allfits, function(fit) fit$regloss))
  scaleloss = unlist(lapply(allfits, function(fit) fit$scaleloss))
  penaltylossx = unlist(lapply(allfits, function(fit) fit$penaltylossx))
  penaltylossy = unlist(lapply(allfits, function(fit) fit$penaltylossy))
  activebeta = unlist(lapply(allfits, function(fit) sum(as.logical(fit$betahat))))
  activedelta = unlist(lapply(allfits, function(fit)  sum(as.logical(fit$deltahat))))
  activedeltabeta = unlist(lapply(allfits, function(fit)  sum((fit$deltahat!=0)[,as.logical(fit$betahat)])))
  betamat = as.data.frame(t(as.data.frame(lapply(allfits, function(fit) fit$betahat))))
  colactivedelta = unlist(lapply(allfits, function(fit)  max(colSums(fit$deltahat!=0)*(fit$betahat!=0))))


  ic = regloss + 2*penaltylossy + penal*log(n)*activebeta + penaldelta*activedeltabeta
  icadj = ic
  icadj[colactivedelta>30] = NA
  # penalize active delta
  # cs-lasso is good but scad is bad

  # plot(ic)
  # plot(activebeta, ic)


  # ic = alpha*regloss + (1-alpha)*scaleloss + lambda_delta*penaltyloss + 2*log(n)*activeseq
  # ic = alpha*regloss + (1-alpha)*scaleloss + lambda_delta*(penaltylossx + penaltylossy) + log(n)*activebeta
  # the original bic selection还是很差，prediction也不是很好， 在zeta = 0.5时

  # ic = alpha*regloss + (1-alpha)*scaleloss + 2*(1-alpha)*penaltylossx + 2*(alpha)*penaltylossy + penal*log(n)*activebeta
  # the classic bic
  # zeta = 0.5 lasso非常差
  # zeta= 1 有很多outlier

  # ic = regloss + 2*penaltylossy + penal*log(n)*activebeta
  # 这样的话scad非常差
  # ic = alpha*regloss + 2*(alpha)*penaltylossy + penal*log(n)*activebeta
  # only consider regression loss
  # zeta =1, lasso and post-lasso is very good but with a few outliers, scad is worse with little contamination

  # ic = alpha*regloss + (1-alpha)*scaleloss + 2*(1-alpha)*activedelta + 2*(alpha)*penaltylossy + log(n)*activebeta
  # penalize active delta
  # cs-lasso is good but scad is bad

  # ic = 2*(regloss + 2*penaltylossy) + 10*log(n)*activebeta

  # penalize how many non-zero delta in active variables
  # penalize active delta
  # cs-lasso is good but scad is bad


  # ic = alpha*regloss + (1-alpha)*scaleloss + 2*(alpha)*penaltylossy + (activebeta + activedelta)
  # penalize on delta


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
                                alpha = alpha, maxiter = maxiter)

      betahat_post[labels] = fit_post$betahat
      intercept_hat_post = fit_post$intercept
    }

  }


  return(list(allfits = allfits, result_opt = result_opt,
              grid = grid, activebeta = activebeta, ic = ic,  colactivedelta = colactivedelta, betamat = betamat, label = label,
              betahat = betahat, intercept_hat = intercept_hat,
              betahat_post = betahat_post, intercept_hat_post = intercept_hat_post))
}


# sregcell_lambda = function(y,x, softbeta = TRUE, softdelta = TRUE, softzeta = TRUE,
#                            lambda_delta = 2.56, lambda_zeta = 2.56, lambda = 0, alpha = 0.5, maxiter = 100){
#
#   n = dim(x)[1]
#   p = dim(x)[2]
#
#   delta = threshold_mat(x,matrix(lambda_delta, n,p))##这个没问题吧
#   xc = x - delta
#   intercept = 0
#   zetahat = rep(0,n)
#   betahat = rep(0,p)
#
#
#   result = reg_beta_delta(y = y, x = x, betahat = betahat, intercept = intercept,
#                           deltahat = delta, zetahat = zetahat,
#                           lambda_beta = lambda, softbeta = softbeta,
#                           lambda_delta = lambda_delta, softdelta = softdelta,
#                           lambda_zeta = lambda_zeta, softzeta = softzeta,
#                           alpha = alpha, maxiter = maxiter)
#   return(result)
# }






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


sregcell_std = function(y,x,
                        scale.method = qnscale,  df = Inf,
                        softbeta = TRUE, softdelta = TRUE, softzeta = TRUE,
                        lambda_delta = NULL, lambda_zeta = 1, prob = 0.995,
                        alpha = 0.5,
                        penal = 1, penaldelta = 0,
                        maxiter = 100){


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

  { # initialization #########could be improved
    # gausscor = robcovsel::paircorxyf(x,y)
    # thresh = sort(abs(gausscor), decreasing = T)[min(10,p)]
    # indx = which(abs(gausscor)>=thresh)
    # xscreen = x[,indx]
    #
    # ximp = purrr::quietly(cellWise::DDC)(xscreen)$result$Ximp
    # fit000 = suppressWarnings(robustbase::lmrob(y~ximp))
    # sigmahat = fit000$scale

    # fit00 = robcovsel::covlasso(xstd,y)
    # sigmahat = sqrt(fit00$sigma2hat_opt)

    set.seed(2)
    fit0 = regcell::Rlars(y,xstd)
    sigmahat = fit0$sigmahat # tend to overestimate

    yscl = y/sigmahat

    # fit00 = regcell::Rlars(yscl,xstd)
    # intercept = fit00$betahat[1]
    # betahat = fit00$betahat[-1]
  }

  rst = sregcell(y = yscl, x = xstd, betahat = rep(0,p), intercept = 0,
                 softbeta = softbeta, softdelta = softdelta, softzeta = softzeta,
                 lambda_delta = lambda_delta, lambda_zeta = lambda_zeta, alpha = alpha,
                 penal = penal, penaldelta = penaldelta,
                 maxiter = 20)

  #de-standardization
  betahat = as.numeric((1/scalex)*rst$betahat*sigmahat)
  intercept_hat = (rst$intercept_hat - sum((1/scalex)*mux*rst$betahat))*sigmahat
  label = rst$label

  betahat_post = (1/scalex)*rst$betahat_post*sigmahat
  intercept_hat_post = (rst$intercept_hat_post - sum((1/scalex)*mux*rst$betahat_post))*sigmahat

  return(list(betahat = betahat, intercept_hat = intercept_hat, label = label,
              betahat_post = betahat_post, intercept_hat_post = intercept_hat_post,
              sigmahat = sigmahat, mux = mux, scalex = scalex))
}



