library(doParallel)
registerDoParallel(cores=10)
getDoParWorkers()

###### data generation settings
{
  ms = 1:20
  ns = 100
  ps = 20 # ps = 200 in high dimensional settings
  prs = c(5,10)
  es = c(0, 0.02, 0.05, 0.1)
  rs = c(0.5)
  gammas = c(0,2,4,6,8,10)
  dfs = c(3, 5, 0)
  outtypes = c("cellwise")
}


{
  m = 9
  n = 100
  p = 200 # ps = 200 in high dimensional settings
  pr = 5
  e = c(0.05)
  r = c(0.5)
  gamma = c(10)
  df = c(0)
  outtype = c("cellwise")
}



##### methods
{
  slmf = function(y,x){srlmcell::slm(y,x)$betahat}
  sssf = function(y, x){shootings::sparseshooting(x,y)$coef_ln}
  rlarsf = function(y, x){srlmcell::Rlars(y, x)$betahat}
  sreg_ada_f = function(y,x){fit = regcell::sregcell(y = y, x = x, adadelta = T,adabeta = T); c(fit$intercept, fit$betahat)}
  sreg_adadelta_f = function(y,x){fit = regcell::sregcell(y = y, x = x, adadelta = T,adabeta = F); c(fit$intercept, fit$betahat)}
  sreg_adabeta_f = function(y,x){fit = regcell::sregcell(y = y, x = x, adadelta = F,adabeta = T); c(fit$intercept, fit$betahat)}
  sreg_f = function(y,x){fit = regcell::sregcell(y = y, x = x, adadelta = F,adabeta = F); c(fit$intercept, fit$betahat)}
  sreg_ada_enet_f = function(y,x){fit = regcell::sregcell(y = y, x = x, adadelta = T,adabeta = T,alpha = 0.9); c(fit$intercept, fit$betahat)}
  sreg_adadelta_enet_f = function(y,x){fit = regcell::sregcell(y = y, x = x, adadelta = T,adabeta = F,alpha = 0.9); c(fit$intercept, fit$betahat)}
  sreg_adabeta_enet_f = function(y,x){fit = regcell::sregcell(y = y, x = x, adadelta = F,adabeta = T,alpha = 0.9); c(fit$intercept, fit$betahat)}
  sreg_enet_f = function(y,x){fit = regcell::sregcell(y = y, x = x, adadelta = F,adabeta = F,alpha = 0.9); c(fit$intercept, fit$betahat)}

  mtds = list(
    rlars = rlarsf,
    sreg_ada = sreg_ada_f,
    sreg_adadelta = sreg_adadelta_f,
    sreg_adabeta = sreg_adabeta_f,
    sreg = sreg_f,
    sreg_enet_ada = sreg_ada_enet_f,
    sreg_enet_adadelta = sreg_adadelta_enet_f,
    sreg_enet_adabeta = sreg_adabeta_enet_f,
    sreg_enet = sreg_enet_f
  )
}



{
  result <- foreach(m = ms,
                    .packages = c("lars", "robustHD","mmlasso","shootings", "cellWise", "regcell","srlmcell"))%:%
    foreach(n = ns)%:%
    foreach(p = ps)%:%
    foreach(pr = prs)%:%
    foreach(e = es)%:%
    foreach(r = rs)%:%
    foreach(gamma = gammas)%:%
    foreach(df = dfs)%:%
    foreach(outtype = outtypes)%dopar% {

      {
        seed = m
        set.seed(seed = seed)
        beta = c(rep(1,pr),rep(0,p-pr))
        dataset = srlmcell::genevar(n = n, p = p, e = e, r = r, beta = beta, gamma = gamma, df = df, outtype = outtype)
        x  = dataset$x
        xc = dataset$xc
        y  = dataset$y
        ynew = dataset$ynew
        #outlierlabel = dataset$outlierlabel
      }

      rst = list()
      if(((e!=0)&(gamma!=0))|((e==0)&(gamma==0))){
        for (mtd in 1:length(mtds)) {
          timing  = system.time({betahat = mtds[[mtd]](y,x)})["elapsed"]
          mape = mean(abs(ynew - cbind(1,xc)%*%betahat))
          rtmspe = sqrt(mean(((ynew - cbind(1,xc)%*%betahat)^2)[1:(0.8*n)]))

          tpf<-function(betahat, beta){sum((as.logical(betahat)==as.logical(beta))[1:pr])}
          tnf<-function(betahat,beta){sum((as.logical(betahat)==as.logical(beta))[-(1:pr)])}
          tp = tpf(betahat[-1], beta)
          tn = tnf(betahat[-1], beta)

          rst[[mtd]] = c(m = m, n = n, p = p, pr = pr, e = e, r = r, gamma = gamma, outtype = outtype, df = df,
                         method = names(mtds)[mtd], seed =seed,
                         MAPE = mape, RTMSPE = rtmspe, TP = tp, FP = (p-pr)-tn, betahat[2:6], Time = timing)
        }
      }else{
        rst[[1]] = rep(NA,21)
      }
      rst
    }

  save(result, file = "result_simu_p20_1213.RData")
}



