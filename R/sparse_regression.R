require(cellWise)


lambdamax_beta = function(y,x,betahat,intercept, softbeta){

  n = dim(x)[1]
  p = dim(x)[2]

  lambdamax = max(t(x)%*%(y-mean(y)))

  iovec = rep(5,2)
  ind = 0
  k = 0
  while(min(iovec)<10){
    outputs_beta = reg_beta(y = y, x = x, betahat = betahat, intercept = intercept, alambdavec_beta = rep(lambdamax,p), softbeta = softbeta)
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
get_lambda_grid = function(y,ximp,betahat,intercept,softbeta, length){
  n = dim(ximp)[1]
  p = dim(ximp)[2]

  lambdamax = lambdamax_beta(y,ximp,betahat, intercept, softbeta)$lambdamax

  lmin = lambdamax/10^3
  grid = c(exp(seq(log(lambdamax),log(lmin),length = length)))
  if(p <n){ grid = c(grid,0)}
  return(grid)
}


## sparse robust regression with imputed design matrix
sregcell = function(y,x, crit = "bic", method.weight = "rpca",
                    softbeta = TRUE, softdelta = TRUE,
                    lambda_delta = 2.56, alpha = 0.5, maxiter = 100){

  n = dim(x)[1]
  p = dim(x)[2]

  if(method.weight=="rpca"){
    fitpca = rob_pca(x, xc = matrix(0,n,p), delta = matrix(0,n,p),lambda = 1/sqrt(min(n,p)), maxiter = 100)
    cellweight = 0.1/(abs(fitpca$delta)+0.1)
    xc = fitpca$xc
    delta = fitpca$delta
    }
  if(method.weight=="ddc"){
    fitddc = cellWise::DDC(x)
    cellweight = 0.1/(abs(x - fitddc$Ximp)+0.1)
    xc = fitddc$Ximp
    delta = x - xc
  }
  if(method.weight=="equal"){
    cellweight = matrix(1,n,p)
    xc = threshold_mat(x,matrix(lambda_delta, n,p))
    delta = x - xc
  }

  intercept = 0
  betahat = rep(0,p)

  length = 100
  grid = get_lambda_grid(y,xc, betahat,intercept, softbeta, length = length)
  allfits = lapply(grid, function(lambda){reg_beta_delta(y = y, x = x, betahat = betahat, intercept = 0,
                                                      deltahat = delta, cellweight = cellweight,
                                                      lambda_beta = lambda, softbeta = softbeta,
                                                      lambda_delta = lambda_delta, softdelta = softdelta,
                                                      alpha = alpha, maxiter = maxiter)})


  regloss = unlist(lapply(1:length, function(i) allfits[[i]]$regloss))
  scaleloss = unlist(lapply(1:length, function(i) allfits[[i]]$scaleloss))
  penaltyloss = unlist(lapply(1:length, function(i) allfits[[i]]$penaltyloss))
  activeseq = unlist(lapply(1:length, function(i)  sum(as.logical(allfits[[i]]$betahat[-1]))))
  ic = alpha*regloss + (1-alpha)*scaleloss + lambda_delta*penaltyloss + 2*log(n)*activeseq

  result_opt = allfits[[which.min(ic)]]

  finallist = list(fits = allfits, result_opt = result_opt)
  return(finallist)
}


sregcell_lambda = function(y,x, adadelta = TRUE, softbeta = TRUE, softdelta = TRUE,
                           lambda_delta = 2.56, lambda = 5*log(length(y)), alpha = 0.5, maxiter = 100 ){

  n = dim(x)[1]
  p = dim(x)[2]

  if(method.weight=="rpca"){
    fitpca = rob_pca(x, xc = matrix(0,n,p), delta = matrix(0,n,p),lambda = 1/sqrt(min(n,p)), maxiter = 100)
    cellweight = 0.1/(abs(fitpca$delta)+0.1)
  }
  if(method.weight=="ddc"){
    fitddc = cellWise::DDC(x)
    cellweight = 0.1/(abs(x - fitddc$Ximp)+0.1)
  }
  if(method.weight=="equal"){
    cellweight = matrix(1,n,p)
  }

  intercept = 0
  betahat = rep(0,p)

  result = reg_beta_delta(y = y, x = x, betahat = betahat, intercept = intercept,
                          deltahat = fitpca$delta, cellweight = cellweight,
                          lambda_beta = lambda,  softbeta = softbeta,
                          lambda_delta = lambda_delta,  softdelta = softdelta,
                          alpha = alpha, maxiter = maxiter)
  return(result)
}






sregcell_step = function(y,x, adadelta = TRUE, softbeta = TRUE, softdelta = TRUE,
                           lambda_delta = 2.56, alpha = 0.5, maxiter = 100 ){

  {
    n = dim(x)[1]
    p = dim(x)[2]
    fitpca = rob_pca(x, xc = matrix(0,n,p), delta = matrix(0,n,p),lambda = 1/sqrt(min(n,p)), maxiter = 100)
    if(adadelta){cellweight = 0.3/(abs(fitpca$delta)+0.1)}else{cellweight = matrix(1,n,p)}
  }

  {
    lambda = max(t(x)%*%(y-mean(y)))
    intercept = 0
    betahat = rep(0,p)
    deltahat = fitpca$delta
  }

  k = 1
  while(TRUE){
    result = reg_beta_delta(y = y, x = x, betahat = betahat, intercept = intercept,
                            deltahat = deltahat, cellweight = cellweight,
                            lambda_beta = lambda,  softbeta = softbeta,
                            lambda_delta = lambda_delta,  softdelta = softdelta,
                            alpha = alpha, maxiter = maxiter)
    intercept = result$intercept
    betahat = result$betahat

    k = k+1
  }




  return(result)
}
















