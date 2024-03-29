#' generate simulated datasets
#'
#' @param n number of observations
#' @param p number of candidate predictors
#' @param pr number of active predictors
#' @param e contamination rate
#' @param r correlation among predictors
#' @param df degrees of freedom used to generate predictors
#' @param beta regression coefficients
#' @param intercept intercept of the regression model
#' @param sigma scale of the error item
#' @param mux averages of predictors
#' @param scalex scales of predictors
#' @param gamma magnitudes of outliers
#' @param outtype type of outliers: cellwise, rowwise
#'
#' @return
#' x: the generated design matrix
#'
#' xc: the generated design matrix without contamination
#'
#' y: the generated response
#'
#' ynew: another generated response (as an independent test set)
#'
#' outlier: the generated outlier in the design matrix
#'
#' erroroutlier: the generated outlier in response
#'
#' @export
#'
#' @examples
#' data = genevar()
#' y = data$y
#' x = data$x
#'
genevar = function(n = 100, p = 10, pr = 5, e = 0, r = 0.5, df = Inf,
                   beta = NULL, intercept = 0, sigma = 1,  mux = rep(0,p), scalex = 1,  gamma = 6,
                   outtype = "cellwise"){

  {
    if(is.null(beta)){
      beta = c(rep(1,pr), rep(0, p-pr))
    }

    # mux = rep(0,p)
    sigmax = diag(rep(scalex^2,p))
    for (i in 1:p) {for (j in 1:p) {
      if (i != j) sigmax[i,j] = sqrt(sigmax[i,i]*sigmax[j,j])*r^abs(i-j)}}
  }

  {
    meanx = matrix(rep(mux, each = n), ncol = p)
    xc = mvtnorm::rmvt(n = n, sigma = sigmax, df = df) + meanx
    errorc = rnorm(n,0,sigma)

    errornew = rnorm(n,0,sigma)
    ynew = intercept + xc%*%beta + errornew

    if(outtype=="cellwise_x"){

      outlierlabel = apply(matrix(0, nrow = n, ncol = p), 2,
                           function(xvec) {xvec[sample(x = 1:n, size = e*n)] = 1; return(xvec)})
      outliervalue = rnorm(n = n*p, mean = gamma, sd = 1)
      outliersign = sample(c(-1,1), size = n*p, replace = T)
      outlier = matrix(outliervalue, nrow = n, ncol=p)*outlierlabel*outliersign

      erroroutlier = rep(0,n)

    }

    if(outtype=="cellwise"){
      outlierlabel = apply(matrix(0, nrow = n, ncol = p), 2,
                           function(xvec) {xvec[sample(x = 1:n, size = e*n)] = 1; return(xvec)})
      outliervalue = rnorm(n = n*p, mean = gamma, sd = 1)
      outliersign = sample(c(-1,1), size = n*p, replace = T)
      outlier = matrix(outliervalue, nrow = n, ncol=p)*outlierlabel*outliersign

      erroroutlierlabel = sample(c(rep(0,(1-e)*n), rep(1,e*n)))
      erroroutliervalue = rnorm(n = n, mean = gamma, sd = 1)
      erroroutliersign = sample(c(-1,1), size = n, replace = T)
      erroroutlier = erroroutlierlabel*erroroutliervalue*erroroutliersign

    }

    if(outtype=="rowwise"){
      xvec = rep(0,n)
      xvec[sample(x = 1:n, size = e*n)] = 1
      outlierlabel = matrix( rep(xvec,p), nrow = n, ncol = p)
      outliervalue = rnorm(n = n*p, mean = gamma, sd = 1)
      outliersign = sample(c(-1,1), size = n*p, replace = T)
      outlier = matrix(outliervalue, nrow = n, ncol=p)*outlierlabel*outliersign

      erroroutlierlabel = xvec
      erroroutliervalue = rnorm(n = n, mean = gamma, sd = 1)
      erroroutliersign = sample(c(-1,1), size = n, replace = T)
      erroroutlier = erroroutlierlabel*erroroutliervalue*erroroutliersign
    }

    x = xc + outlier#*sign(xc - meanx)
    error = errorc + erroroutlier#*sign(errorc)
    y = intercept + xc%*%beta + error



  }
  return(list(x = x, xc = xc, y = y, ynew = ynew, beta = beta,
              intercept = intercept, sigma = sigma, mux = mux, sigmax = sigmax,
              outlier = outlier, erroroutlier = erroroutlier))


}



