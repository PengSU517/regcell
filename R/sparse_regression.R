# lambdamax_beta = function(y,x,betahat,intercept, softbeta){
#
#   n = dim(x)[1]
#   p = dim(x)[2]
#
#   lambdamax = max(t(x)%*%(y-mean(y)))
#
#   iovec = rep(5,2)
#   ind = 0
#   k = 0
#   while((min(iovec)<10)|(ind=FALSE)){
#     outputs_beta = reg_beta(y = y, x = x, betahat = betahat, intercept = intercept,
#                             alambdavec_beta = rep(lambdamax,p), softbeta = softbeta)
#
#     betaget = outputs_beta$betahat
#     ind = (sum(abs(betaget))==0)
#     lambdamax = ((1/2)^((ind-0.5)*2*1/(2^min(iovec))))*lambdamax
#     iovec[ind+1] = iovec[ind+1]+1
#     k = k+1
#   }
#   result = list(iovec = iovec, lambdamax = lambdamax, k = k)
#   return(result)
#
# }

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

  lambdamax = 20*max(t(x - deltahat)%*%(y-intercept - zetahat))
  outputs = reg_beta_delta(y = y, x = x, betahat = betahat, intercept = intercept,
                           deltahat = deltahat, zetahat = zetahat,
                           lambda_beta = lambdamax, softbeta = softbeta,
                           lambda_delta = lambda_delta, softdelta = softdelta,
                           lambda_zeta = lambda_zeta, softzeta = softzeta,
                           alpha = alpha, maxiter = maxiter)
  lambdamax = max(abs(outputs$mgradient))*1.01
  return(list(lambdamax = lambdamax))
}


## sparse robust regression with imputed design matrix
sregcell = function(y,x, betahat = NULL, intercept = NULL,
                    softbeta = TRUE, softdelta = TRUE, softzeta = TRUE,
                    lambda_delta = 2.56, lambda_zeta = 2.56, alpha = 0.5,
                    penal_beta = 2, penal_cell = 2, maxiter = 100){


  n = dim(x)[1]
  p = dim(x)[2]


  deltahat = threshold_mat(x,matrix(lambda_delta, n,p))##这个没问题吧
  xc = x - deltahat
  zetahat = rep(0,n)

  if(is.null(betahat)){
    betahat = rep(0,p)
    intercept = 0
  }

  { #lambda grid
    length = 50
    # grid = get_lambda_grid(y,xc, betahat,intercept, softbeta, length = length)
    lambdamax = lambdamax3_beta(y = y, x = x, betahat = betahat, intercept = intercept,
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


  regloss = unlist(lapply(1:length, function(i) allfits[[i]]$regloss))
  scaleloss = unlist(lapply(1:length, function(i) allfits[[i]]$scaleloss))
  penaltyloss = unlist(lapply(1:length, function(i) allfits[[i]]$penaltyloss))
  penaltylossx = unlist(lapply(1:length, function(i) allfits[[i]]$penaltylossx))
  penaltylossy = unlist(lapply(1:length, function(i) allfits[[i]]$penaltylossy))
  activebeta = unlist(lapply(1:length, function(i)  sum(as.logical(allfits[[i]]$betahat))))#为什么会减一啊
  loss = unlist(lapply(1:length, function(i) allfits[[i]]$loss))

  activedelta = unlist(lapply(1:length, function(i)  sum(as.logical(allfits[[i]]$delta))))
  activezeta = unlist(lapply(1:length, function(i)  sum(as.logical(allfits[[i]]$zeta))))

  # penal_beta = log(n)
  # penal_cell = 1
  # plot(ic)
  # plot(activezeta)
  # plot(activedelta)
  # plot(penaltyloss)
  # plot(penaltylossx)
  # plot(penaltylossy)
  # plot(regloss + penaltylossy)
  # plot(activebeta + activedelta)
  # plot(regloss + penaltylossy + log(n)*(activebeta + activedelta))

  # ic = alpha*regloss + (1-alpha)*scaleloss + penal_beta*(activebeta) + penal_cell*(penaltyloss)
  ic = alpha*regloss + (1-alpha)*scaleloss + penaltylossy + penal_beta*(activebeta + activedelta)
  result_opt = allfits[[which.min(ic)]]
  # result_opt$betahat

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
quantscale = function(xvec, df = Inf, tau = 0.9){
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


sregcell_std = function(y,x, scale.method = qnscale,
                        softbeta = TRUE, softdelta = TRUE, softzeta = TRUE,
                        df = Inf,
                        lambda_delta = NULL, lambda_zeta = 2.58,
                        alpha = 0.5,
                        penal_beta = 2, penal_cell = 2,
                        maxiter = 100){


  n = dim(x)[1]
  p = dim(x)[2]

  if(is.null(lambda_delta)){
    lambda_delta = qt(0.995, df = df)
  }

  { # standardization
    stdrst = robstd(x, centerf = median, scalef = scale.method, df = df)
    mux = stdrst$center
    scalex = stdrst$scale
    xstd = stdrst$x

  }

  { # initialization #####################################could be improved
    fit0 = regcell::Rlars(y,x)
    sigmahat = fit0$sigmahat
    yscl = y/sigmahat
    fit00 = regcell::Rlars(yscl,xstd)
    intercept = fit00$betahat[1]
    betahat = fit00$betahat[-1]
  }

  rst = sregcell(y = yscl, x = xstd, betahat = betahat, intercept = intercept,
                 softbeta = softbeta, softdelta = softdelta, softzeta = softzeta,
                 lambda_delta = lambda_delta, lambda_zeta = lambda_zeta, alpha = alpha,
                 penal_beta = penal_beta, penal_cell = penal_cell, maxiter = 100)

  #de-standardization
  betahat = as.numeric((1/scalex)*rst$betahat*sigmahat)
  intercept_hat = (rst$intercept_hat - sum((1/scalex)*mux*rst$betahat))*sigmahat

  betahat_post = (1/scalex)*rst$betahat_post*sigmahat
  intercept_hat_post = (rst$intercept_hat_post - sum((1/scalex)*mux*rst$betahat_post))*sigmahat

  return(list(betahat = betahat, intercept_hat = intercept_hat,
              betahat_post = betahat_post, intercept_hat_post = intercept_hat_post,
              sigmahat = sigmahat, mux = mux, scalex = scalex))
}



