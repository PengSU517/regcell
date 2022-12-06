require(cellWise)

#### get lambda grid
get_lambda_gridm = function(y,x,betahat, sigmahat,reg){
  n = dim(x)[1]
  p = dim(x)[2]
  lambdamax = lambdamax_matan(y, x, betahat, sigmahat,reg)$lambdamax
  lmin = lambdamax/10^3
  grid = c(exp(seq(log(lambdamax),log(lmin),length = 50)))
  if(p <n){ grid = c(grid,0)}
  return(grid)
}


## sparse robust regression with imputed design matrix
srlmm = function(y,x, initial = "rlars", iter = TRUE, tech = "row", reg = "robust", crit = "bic", lambda = NULL){

  ximp <- suppressMessages(cellWise::DDC(x)$Ximp)
  if(initial == "ddc"  ){fit0 = slm(y, ximp)}
  if(initial == "rlars"){fit0 = Rlars(y, x)}

  betahat  = fit0$betahat
  sigmahat = fit0$sigmahat

  n = dim(x)[1]
  p = dim(x)[2]

  if(is.null(lambda)){
    grid = get_lambda_gridm(y,ximp, betahat, sigmahat, reg)
    if(iter){
      fit = lapply(grid, function(lambda){TukeyM_atan_iter(y, x, ximp, betahat, sigmahat, tech, lambda,reg)})
    }else{
      fit = lapply(grid, function(lambda){TukeyM_atan(y, ximp, betahat, sigmahat, lambda,reg)})
    }

    if(crit == "aic"){
      ic = unlist(lapply(fit, function(x) 2*(x$loss) + 2*sum(as.logical(x$betahat[-1]))))
    }else{
      if(crit == "bic"){
        ic = unlist(lapply(fit, function(x) 2*(x$loss) + (log(n))*sum(as.logical(x$betahat[-1]))))
        #unlist(lapply(fit, function(x) 2*(x$loss)))
      }else{
        print("Unappropriate IC")
      }
    }

    result = fit[[which.min(ic)]]
    result$lambda_opt = grid[which.min(ic)]
  }else{
    if(iter){
      result = TukeyM_atan_iter(y, x, ximp, betahat, sigmahat, tech, lambda,reg)
    }else{
      result = TukeyM_atan(y, ximp, betahat, sigmahat, lambda,reg)
    }
  }

  xtilde = result$xtilde
  result$flagger = (x!=xtilde)

  return(result)
}



## robust regression with imputed design matrix
rlmm = function(y,x, initial = "rlars", iter = TRUE, tech = "row"){
  ximp <- suppressMessages(cellWise::DDC(x)$Ximp)
  if(initial == "ddc"  ){fit0 = slm(y, ximp)}
  if(initial == "rlars"){fit0 = Rlars(y, x)}

  betahat  = fit0$betahat
  sigmahat = fit0$sigmahat

  n = dim(x)[1]

  p = dim(x)[2]

  if(iter){
    result = TukeyM_iter(y, x, ximp, betahat, sigmahat, tech)
  }else{
    result = TukeyM(y, ximp, betahat, sigmahat)
  }

  xtilde = result$xtilde
  if(is.null(xtilde)){xtilde = ximp}
  result$flagger = (x!=xtilde)

  return(result)
}


