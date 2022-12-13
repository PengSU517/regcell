require(robustbase)
require(robustHD)


## robust lars ###tried to use AIC
Rlars = function(y, x){
  fit = suppressWarnings(robustHD::rlars(x,y))
  #winsorize = F, prob = 0.99,
  betahat = fit$coefficients[,fit$crit$best]
  sigmahat = fit$scale[fit$crit$best]
  res = fit$residuals[,fit$crit$best]
  return(list(betahat = betahat, sigmahat = sigmahat, res = res))

}

### sparse linear model
slm = function(y,x, type = "stepwise"){
  n = dim(x)[1]
  p = dim(x)[2]
  fit = lars::lars(x, y, type = type, max.steps = min(n/4, p),intercept = T)

  bic = fit$RSS + fit$df*log(n)##

  betahat = fit$beta[which.min(bic),]
  beta0 = mean(y - x%*%betahat)
  betahat = c(0, betahat)
  sigmahat = sqrt(fit$RSS[which.min(bic)]/n)
  return(list(betahat = betahat, sigmahat = sigmahat))
}




