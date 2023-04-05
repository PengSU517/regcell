library(doParallel)
registerDoParallel(cores=10)
getDoParWorkers()

###### data generation settings

{
  ms = 1:200
  ns = c(200)
  ps = c(50)
  prs = c(10)
  es = c(0, 0.02, 0.05)
  rs = c(0.5)
  gammas = c(0, 2, 4, 6, 8)
  # sigmas = c(2)
  dfs = c(4, Inf)
  outtypes = c("mixed")
}

# {
#   m = 1
#   n = 100
#   p = c(20)
#   pr = c(10)
#   e = c(0.05)
#   r = c(0.5)
#   gamma = c(6)
#   df = c(Inf)
#   outtype = c("mixed")
# }




##### methods
mtds = list(
  sss = function(y, x){shootings::sparseshooting(x,y)$coef},
  rlars = function(y, x){regcell::Rlars(y, x)$betahat},
  mmlasso = function(y,x){mmlasso::mmlasso(x,y)$coef.MMLasso.ad},
  slts = function(y,x){robustHD::sparseLTS(x,y)$coefficients},

  # cell_scad = function(y,x, penal, penaldelta){
  #   fit = regcell::sregcell_std(y = y, x = x, softbeta = FALSE, lambda_zeta = 1, penal = 0.3, penaldelta = 0)
  #   return(c(fit$intercept_hat, fit$betahat))
  # },

  cell_lasso_post = function(y,x){
    fit = regcell::sregcell_std(y = y, x = x, softbeta = TRUE, lambda_zeta = 1, penal = 1, penaldelta =0)
    return(c(fit$intercept_hat_post, fit$betahat_post))
  },

  # cell_lasso_post2 = function(y,x){
  #   fit = regcell::sregcell_std(y = y, x = x, softbeta = TRUE, lambda_zeta = 1, penal = 2, penaldelta =0)
  #   return(c(fit$intercept_hat_post, fit$betahat_post))
  # },

  lasso  = function(y,x){
    return(regcell::slm(y,x, type = "lasso")$betahat)
  }

)

{
  result <- foreach(m = ms,
                    .packages = c("lars", "robustHD", "robustbase" , "mmlasso",
                                  "shootings", "cellWise", "regcell", "purrr", "robcovsel"))%:%
    foreach(n = ns)%:%
    foreach(p = ps)%:%
    foreach(pr = prs)%:%
    foreach(e = es)%:%
    foreach(r = rs)%:%
    foreach(gamma = gammas)%:%
    # foreach(sigma = sigmas)%:%
    foreach(df = dfs)%:%
    foreach(outtype = outtypes)%dopar% {

      {
        seed = m
        set.seed(seed = seed)
        beta = c(rep(1,pr),rep(0,p-pr))#*c(1,-1)#正负号影响了信噪比
        dataset = regcell::genevar(n = n, p = p, e = e, r = r, beta = beta,intercept = 1,
                                   gamma = gamma, df = df, outtype = outtype, sigma = 3,
                                   mux = rep(0,p), scalex = 1)##gamma和scalex应该有关系
        x  = dataset$x
        xc = dataset$xc
        y  = dataset$y
        ynew = dataset$ynew
      }

      rst = list()
      if(((e!=0)&(gamma!=0))|((e==0)&(gamma==0))){
        for (mtd in 1:length(mtds)) {
          rst[[mtd]]=rep(NA, 27)
          try({
            timing  = system.time({betahat = mtds[[mtd]](y,x)})["elapsed"]
            mspe = mean((ynew - cbind(1,xc)%*%betahat)^2)
            mape = mean(abs(ynew - cbind(1,xc)%*%betahat))
            rtmspe = sqrt(mean(sort((ynew - cbind(1,xc)%*%betahat)^2)[1:(0.9*n)]))##trimed mean 要先排序啊

            tpf<-function(betahat, beta){sum((as.logical(betahat)==as.logical(beta))[1:pr])}
            tnf<-function(betahat,beta){sum((as.logical(betahat)==as.logical(beta))[-(1:pr)])}
            tp = tpf(betahat[-1], beta)
            tn = tnf(betahat[-1], beta)

            rst[[mtd]] = c(m = m, n = n, p = p, pr = pr, e = e, r = r,
                           gamma = gamma, outtype = outtype, df = df,
                           method = names(mtds)[mtd], seed =seed,
                           MSPE= mspe, MAPE = mape, RTMSPE = rtmspe,
                           TP = tp, FP = (p-pr)-tn, Time = timing, betahat[2:11])
          }, TRUE
          )

        }
      }else{
        rst[[1]] = rep(NA,27)
      }
      rst
    }

  save(result, file = "result_0404_n200_sigma3.RData")
}



