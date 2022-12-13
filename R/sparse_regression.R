require(cellWise)

#### get lambda grid
get_lambda_grid = function(y,x,betahat,intercept,adabeta){
  n = dim(x)[1]
  p = dim(x)[2]
  lambdamax = lambdamax_beta(y, x, betahat, intercept, adabeta)$lambdamax####可以直接确定lambda
  lmin = lambdamax/10^3
  grid = c(exp(seq(log(lambdamax),log(lmin),length = 100)))
  if(p <n){ grid = c(grid,0)}
  return(grid)
}


## sparse robust regression with imputed design matrix
sregcell = function(y,x, initial = "rlars", crit = "bic",
                    adadelta = TRUE, adabeta = FALSE,
                    lambda = NULL, lambda_delta = 2.56/2, alpha = 0.5){

  ximp <- suppressMessages(cellWise::DDC(x)$Ximp)
  if(initial == "ddc"  ){fit0 = slm(y, ximp)}
  if(initial == "rlars"){fit0 = Rlars(y, x)}

  invstdres2 = 1/abs(fit0$res/fit0$sigmahat/lambda_delta)^2

  n = dim(x)[1]
  p = dim(x)[2]

  if(adadelta){
    cellweight =  matrix(invstdres2,nrow = n, ncol = p)##可能不合适，lambda=0时会过拟合
  }else{
    cellweight =  matrix(1,nrow = n, ncol = p)
  }


  sigmahat = fit0$sigmahat
  intercept = fit0$betahat[1]/sigmahat
  betahat  = fit0$betahat[-1]/sigmahat
  ynew = y/sigmahat

  grid = get_lambda_grid(ynew,ximp, betahat,intercept,adabeta)
  fits = lapply(grid, function(lambda){reg_beta_delta(y = ynew, x = x, betahat = betahat, intercept = intercept,
                                                      deltahat = x-ximp, cellweight = cellweight, adabeta = adabeta,
                                                      lambda_beta = lambda,lambda_delta = lambda_delta,alpha = alpha)})

  if(crit == "aic"){ic = unlist(lapply(fits, function(fit) (fit$loss) + 2*sum(as.logical(fit$betahat[-1]))))}
  if(crit == "bic"){ic = unlist(lapply(fits, function(fit) (fit$loss) + (log(n))*sum(as.logical(fit$betahat[-1]))))}
  result = fits[[which.min(ic)]]
  result$lambda_opt = grid[which.min(ic)]
  result$intercept = result$intercept*sigmahat
  result$betahat = result$betahat*sigmahat
  ##result = reg_beta_delta(y = ynew, x = x, betahat = betahat, intercept = intercept, deltahat = x-ximp,lambda_beta = lambda)

  return(result)
}


sregcell_lambda = function(y,x, initial = "rlars", crit = "bic", lambda = 1){

  ximp <- suppressMessages(cellWise::DDC(x)$Ximp)
  if(initial == "ddc"  ){fit0 = slm(y, ximp)}
  if(initial == "rlars"){fit0 = Rlars(y, x)}

  sigmahat = fit0$sigmahat
  intercept = fit0$betahat[1]/sigmahat
  betahat  = fit0$betahat[-1]/sigmahat
  ynew = y/sigmahat

  result = reg_beta_delta(y = ynew, x = x, betahat = betahat, intercept = intercept, deltahat = x-ximp, lambda_beta = lambda)

  result$betahat
  result$loss
  result$intercept = result$intercept*sigmahat
  result$betahat = result$betahat*sigmahat
  return(result)
}
