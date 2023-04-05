library(doParallel)
registerDoParallel(cores=10)
getDoParWorkers()

###### data generation settings

{
  ms = 1:200
  ns = c(100, 200)
  ps = c(20, 50, 100, 200)
  prs = c(5, 10)
  es = c(0.05)
  sigmas = c(1,2)
  rs = c(0.5)
  gammas = c(6)
  dfs = Inf
  penals = seq(1,5, 1)
  penaldeltas = 0
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
#   penal = 1
#   penaldelta = 1
#   outtype = c("mixed")
# }




##### methods
mtds = list(

  cell_lasso_post = function(y,x, penal, penaldelta){
    fit = regcell::sregcell_std(y = y, x = x, softbeta = TRUE, lambda_zeta = 1, penal = penal, penaldelta = penaldelta)
    return(c(fit$intercept_hat_post, fit$betahat_post))
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
    foreach(sigma = sigmas)%:%
    foreach(r = rs)%:%
    foreach(penal = penals)%:%
    foreach(penaldelta = penaldeltas)%:%
    foreach(gamma = gammas)%:%
    foreach(df = dfs)%:%
    foreach(outtype = outtypes)%dopar% {

      {
        seed = m
        set.seed(seed = seed)
        beta = c(rep(1,pr),rep(0,p-pr))#*c(1,-1)#正负号影响了信噪比
        dataset = regcell::genevar(n = n, p = p, e = e, r = r, beta = beta,intercept = 1,
                                   gamma = gamma, df = df, outtype = outtype, sigma = sigma,
                                   mux = rep(0,p), scalex = 1)##gamma和scalex应该有关系
        x  = dataset$x
        xc = dataset$xc
        y  = dataset$y
        ynew = dataset$ynew
      }

      rst = list()
      if(((e!=0)&(gamma!=0))|((e==0)&(gamma==0))){
        for (mtd in 1:length(mtds)) {
          rst[[mtd]]=rep(NA, 30)
          try({
            timing  = system.time({betahat = mtds[[mtd]](y,x, penal, penaldelta)})["elapsed"]
            mspe = mean((ynew - cbind(1,xc)%*%betahat)^2)
            mape = mean(abs(ynew - cbind(1,xc)%*%betahat))
            rtmspe = sqrt(mean(sort((ynew - cbind(1,xc)%*%betahat)^2)[1:(0.9*n)]))##trimed mean 要先排序啊

            tpf<-function(betahat, beta){sum((as.logical(betahat)==as.logical(beta))[1:pr])}
            tnf<-function(betahat,beta){sum((as.logical(betahat)==as.logical(beta))[-(1:pr)])}
            tp = tpf(betahat[-1], beta)
            tn = tnf(betahat[-1], beta)

            rst[[mtd]] = c(m = m, n = n, p = p, pr = pr, e = e, sigma = sigma, r = r,
                           gamma = gamma, outtype = outtype, df = df, penal = penal, penaldelta = penaldelta,
                           method = names(mtds)[mtd], seed =seed,
                           MSPE= mspe, MAPE = mape, RTMSPE = rtmspe,
                           TP = tp, FP = (p-pr)-tn, Time = timing, betahat[2:11])
          }, TRUE
          )

        }
      }else{
        rst[[1]] = rep(NA,30)
      }
      rst
    }

  save(result, file = "result_debug_0330.RData")
}



