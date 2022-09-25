require(robustbase)
require(robustHD)

## robust lars ###tried to use AIC
Rlars = function(y, x){
  fit = suppressWarnings(robustHD::rlars(x,y))
  #winsorize = F, prob = 0.99,
  betahat = fit$coefficients[,fit$crit$best]
  sigmahat = fit$scale[fit$crit$best]
  return(list(betahat = betahat, sigmahat = sigmahat))

}

### sparse linear model
slm = function(y,x, type = "stepwise", crit = "logbic"){
  n = dim(x)[1]
  p = dim(x)[2]
  fit = lars::lars(x, y, type = type, max.steps = min(n/4, p),intercept = T)

  if(crit == "logbic"){
    bic = n*log(fit$RSS/n) + fit$df*log(n)##
  }else{
    bic = fit$RSS + fit$df*log(n)##
  }

  ## doesn't work in high dimensional cases!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  betahat = fit$beta[which.min(bic),]
  beta0 = mean(y - x%*%betahat)
  betahat = c(0, betahat)
  sigmahat = sqrt(fit$RSS[which.min(bic)]/n)
  return(list(betahat = betahat, sigmahat = sigmahat))
}

