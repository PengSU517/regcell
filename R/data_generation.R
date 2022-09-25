genevar = function(n = 100, p = 20, e = 0,r = 0.5, df = 0,
                   beta = c(1,2,1,2,1), gamma = 6, errorsigma = 1,
                   outtype = "cellwise"){

  {
    beta = c(beta, rep(0,(p-length(beta))))
    mu = rep(0,p)
    sigma = diag(rep(1^2,p))
    for (i in 1:p) {for (j in 1:p) {
      if (i !=j)sigma[i,j] = sqrt(sigma[i,i]*sigma[j,j])*r^abs(i-j)}}
  }

  {
    mean = matrix(rep(mu, each = n), ncol = p)
    xc = mvtnorm::rmvt(n = n, sigma = sigma, df = df)

    ##creating response
    error = rnorm(n,0,errorsigma)
    y = xc%*%beta + error

    errornew = rnorm(n,0,errorsigma)
    ynew = xc%*%beta + errornew


    #creating outliers
    if(outtype=="cellwise"){
      outlierlabel = apply(matrix(0, nrow = n, ncol = p), 2,
                           function(xvec) {xvec[sample(x = 1:n, size = e*n)] = 1; return(xvec)})
    }
    if(outtype=="rowwise"){
      xvec = rep(0,n)
      xvec[sample(x = 1:n, size = e*n)] = 1
      outlierlabel = matrix( rep(xvec,p), nrow = n, ncol = p)
    }

    outliervalue = rnorm(n = n*p, mean = gamma, sd = 1)
    outliersign = sample(c(-1,1), size = n*p, replace = T)
    outlier = matrix(outliervalue*outliersign, nrow = n, ncol=p)
    x = xc*(1-outlierlabel)+outlier*outlierlabel

  }
  return(list(x = x, xc = xc, y = y, ynew = ynew, beta = beta, outlierlabel = outlierlabel,
              mu = mu, Sigma = sigma, sigma = errorsigma))


}



