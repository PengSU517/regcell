require(cellWise)

#### get lambda grid
get_lambda_gridm = function(y,x,betahat, sigmahat){
  n = dim(x)[1]
  p = dim(x)[2]
  lambdamax = lambdamax_matan(y, x, betahat, sigmahat)$lambdamax
  lmin = lambdamax/10^3
  grid = c(exp(seq(log(lambdamax),log(lmin),length = 100)))
  if(p <n){ grid = c(grid,0)}
  return(grid)
}


#### get lambda grid
get_lambda_grids = function(y,x,betahat, sigmahat){
  n = dim(x)[1]
  p = dim(x)[2]
  lambdamax = lambdamax_satan(y, x, betahat, sigmahat)$lambdamax
  lmin = lambdamax/10^3
  grid = c(exp(seq(log(lambdamax),log(lmin),length = 100)))
  if(p <n){ grid = c(grid,0)}
  return(grid)
}

## sparse robust regression with imputed design matrix
srlms = function(y,x, initial = "ddc", iter = TRUE, crit = "bic"){

  ximp <- suppressMessages(cellWise::DDC(x)$Ximp)
  if(initial == "ddc"  ){fit0 = slm(y, ximp)}
  if(initial == "rlars"){fit0 = Rlars(y, x)}

  betahat  = fit0$betahat
  sigmahat = fit0$sigmahat

  n = dim(x)[1]
  p = dim(x)[2]
  grid = get_lambda_grids(y,ximp, betahat, sigmahat)
  if(iter){
    fit = lapply(grid, function(lambda){TukeyS_atan_iter(y, x, ximp, betahat, sigmahat, lambda)})
    }else{
      fit = lapply(grid, function(lambda){TukeyS_atan(y, x, betahat, sigmahat, lambda)})
    }
  if(crit == "ebic"){
    bic = unlist(lapply(fit, function(x) (x$logloss) + (log(n)+log(p))*sum(as.logical(x$betahat[-1]))))
  }else{
    bic = unlist(lapply(fit, function(x) (x$logloss) + (log(n))*sum(as.logical(x$betahat[-1]))))
  }

  result = fit[[which.min(bic)]]
  return(result)
}


## sparse robust regression with imputed design matrix
srlmm = function(y,x, initial = "ddc", iter = TRUE, tech = "row", crit = "bic"){

  ximp <- suppressMessages(cellWise::DDC(x)$Ximp)
  if(initial == "ddc"  ){fit0 = slm(y, ximp)}
  if(initial == "rlars"){fit0 = Rlars(y, x)}

  betahat  = fit0$betahat
  sigmahat = fit0$sigmahat

  n = dim(x)[1]
  p = dim(x)[2]
  grid = get_lambda_gridm(y,ximp, betahat, sigmahat)
  if(iter){
    fit = lapply(grid, function(lambda){TukeyM_atan_iter(y, x, ximp, betahat, sigmahat, tech, lambda)})
  }else{
    fit = lapply(grid, function(lambda){TukeyM_atan(y, ximp, betahat, sigmahat, lambda)})
  }
  if(crit == "ebic"){
    bic = unlist(lapply(fit, function(x) 2*(x$loss) + (log(n)+log(p))*sum(as.logical(x$betahat[-1]))))
  }else{
    bic = unlist(lapply(fit, function(x) 2*(x$loss) + (log(n))*sum(as.logical(x$betahat[-1]))))
  }

  result = fit[[which.min(bic)]]
  result$lambda_opt = grid[which.min(bic)]

  return(result)
}


