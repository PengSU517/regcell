require(cellWise)


lambdamax_beta = function(y,x,betahat,intercept,lambdamax = 1){

  n = dim(x)[1]
  p = dim(x)[2]
  iovec = rep(0,2)
  ind = 0
  k = 0
  while(min(iovec)<5){
    outputs_beta = reg_beta(y = y, x = x, betahat = betahat, intercept = intercept, lambdavec_beta = rep(lambdamax,p))
    betaget = outputs_beta$betahat
    ind = (sum(abs(betaget))==0)
    lambdamax = ((1/2)^((ind-0.5)*2*1/(2^min(iovec))))*lambdamax
    iovec[ind+1] = iovec[ind+1]+1
    k = k+1
  }
  result = list(iovec = iovec, lambdamax = lambdamax, k = k)
  return(result)

}





#### get lambda grid
get_lambda_grid = function(y,ximp,betahat,intercept,adabeta, length){
  n = dim(ximp)[1]
  p = dim(ximp)[2]

  lambdamax = lambdamax_beta(y,ximp,betahat, intercept)$lambdamax

  lmin = lambdamax/10^3
  grid = c(exp(seq(log(lambdamax),log(lmin),length = length)))
  if(p <n){ grid = c(grid,0)}
  return(grid)
}


## sparse robust regression with imputed design matrix
sregcell = function(y,x, initial = "rlars", crit = "bic",
                    adadelta = TRUE, adabeta = FALSE, softbeta = TRUE, softdelta = TRUE,
                    lambda_delta = 2.56, alpha = 0.5, maxiter = 100){

  {
    ximp <- suppressMessages(cellWise::DDC(x)$Ximp)
    #if(initial == "ddc"  ){fit0 = slm(y, ximp)}
    fit0 = Rlars(y, x)

    n = dim(x)[1]
    p = dim(x)[2]

    sigmahat = fit0$sigmahat
    intercept = fit0$betahat[1]#/sigmahat
    #betahat  = fit0$betahat[-1]#/sigmahat
    betahat = rep(0,p)
    ynew = y#/sigmahat

    rowweight = 1/abs(fit0$res/fit0$sigmahat/lambda_delta)
    colweight = 1/abs(x/lambda_delta)
  }

  length = 100
  grid = get_lambda_grid(ynew,ximp, betahat,intercept, length = length)
  allfits = lapply(grid, function(lambda){reg_beta_delta(y = ynew, x = x, betahat = betahat, intercept = intercept,
                                                      deltahat = x-ximp, rowweight = rowweight, colweight = colweight,
                                                      lambda_beta = lambda, adabeta = adabeta, softbeta = softbeta,
                                                      lambda_delta = lambda_delta, adadelta = adadelta, softdelta = softdelta,
                                                      alpha = alpha, maxiter = maxiter)})


  regloss = unlist(lapply(1:length, function(i) allfits[[i]]$regloss))
  scaleloss = unlist(lapply(1:length, function(i) allfits[[i]]$scaleloss))
  penaltyloss = unlist(lapply(1:length, function(i) allfits[[i]]$penaltyloss))
  activeseq = unlist(lapply(1:length, function(i)  sum(as.logical(allfits[[i]]$betahat[-1]))))
  ic = alpha*regloss + (1-alpha)*scaleloss + penaltyloss + log(n)*activeseq

  result_opt = allfits[[which.min(ic)]]


  # plot(activeseq)
  # plot(activeseq,regloss)
  # plot(activeseq,scaleloss)
  # plot(activeseq, penaltyloss)
  # plot(activeseq, ic)
  #
  # max(abs(x - allfits[[1]]$deltahat))
  #
  #
  # fit1 = reg_beta_delta(y = ynew, x = x, betahat = betahat, intercept = intercept,
  #                       deltahat = x-ximp, rowweight = rowweight, colweight = colweight,
  #                       lambda_beta = grid[1], adabeta = adabeta, softbeta = softbeta,
  #                       lambda_delta = lambda_delta, adadelta = adadelta, softdelta = softdelta,
  #                       alpha = alpha, maxiter = maxiter)
  #
  # fit1$betahat
  # fit1$deltahat
  # fit1$penaltyloss

  #result_opt$lambda = grid[which.min(ic)]
  #result_opt$intercept = result_opt$intercept*sigmahat
  #result_opt$betahat = result_opt$betahat*sigmahat

  finallist = list(fits = allfits, result_opt = result_opt)
  return(finallist)
}


sregcell_lambda = function(y,x, initial = "rlars", crit = "bic",
                           adadelta = TRUE, adabeta = FALSE, softbeta = TRUE, softdelta = TRUE,
                           lambda_delta = 2.56, lambda = 5*log(length(y)), alpha = 0.5, maxiter = 100 ){

  {
    ximp <- suppressMessages(cellWise::DDC(x)$Ximp)
    #if(initial == "ddc"  ){fit0 = slm(y, ximp)}
    fit0 = Rlars(y, x)

    n = dim(x)[1]
    p = dim(x)[2]

    sigmahat = fit0$sigmahat
    intercept = fit0$betahat[1]#/sigmahat
    #betahat  = fit0$betahat[-1]#/sigmahat
    betahat = rep(0,p)
    ynew = y#/sigmahat

    rowweight = 1/abs(fit0$res/fit0$sigmahat/lambda_delta)
    colweight = 1/abs(x/lambda_delta)
  }


  result = reg_beta_delta(y = ynew, x = x, betahat = betahat, intercept = intercept,
                          deltahat = x-ximp, rowweight = rowweight, colweight = colweight,
                          lambda_beta = lambda, adabeta = adabeta, softbeta = softbeta,
                          lambda_delta = lambda_delta, adadelta = adadelta, softdelta = softdelta,
                          alpha = alpha, maxiter = maxiter)

  #result$intercept = result$intercept*sigmahat
  #result$betahat = result$betahat*sigmahat
  #result$loss

  return(result)
}
