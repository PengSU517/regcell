require(robustbase)
require(robustHD)



#' robust Lars
#'
#' @param y response
#' @param x design matrix
#'
#' @return
#' betahat: estimated beta
#'
#' sigmahat: estimated scale
#'
#' res: residuals
#' @export
#'
#' @examples
#'
#' data = genevar()
#' y = data$y
#' x = data$x
#' fit = Rlars(y,x)
Rlars = function(y, x){
  fit = suppressWarnings(robustHD::rlars(x,y))
  #winsorize = F, prob = 0.99,
  betahat = fit$coefficients[,fit$crit$best]
  sigmahat = fit$scale[fit$crit$best]
  res = fit$residuals[,fit$crit$best]
  return(list(betahat = betahat, sigmahat = sigmahat, res = res))

}

### sparse linear model
slm = function(y,x, type = "lasso"){
  n = dim(x)[1]
  p = dim(x)[2]
  fit = lars::lars(x, y, type = type, #max.steps = min(n/4, p),
                   intercept = T)

  bic = fit$RSS + fit$df*log(n)##

  betahat = fit$beta[which.min(bic),]
  beta0 = mean(y - x%*%betahat)
  betahat = c(beta0, betahat)
  sigmahat = sqrt(fit$RSS[which.min(bic)]/n)
  return(list(betahat = betahat, sigmahat = sigmahat))
}


#' Lasso with cross validation
#'
#' @param y response
#' @param x the design matrix
#'
#' @return
#' betahat: the estimated regression coefficient vector
#' @export
#'
#' @examples
#' data = genevar()
#' y = data$y
#' x = data$x
#' fit = lassocv(y,x)
#'
lassocv = function(y,x){
  fit = glmnet::cv.glmnet(x, y)
  betatilde = as.numeric(coef(fit, s = "lambda.min"))
  return(list(betahat = betatilde))
}








