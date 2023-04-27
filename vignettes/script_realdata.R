library(doParallel)
registerDoParallel(cores=10)
getDoParWorkers()


mtds = list(
  sss = function(y, x){shootings::sparseshooting(x,y)$coef},
  rlars = function(y, x){regcell::Rlars(y, x)$betahat},
  mmlasso = function(y,x){mmlasso::mmlasso(x,y)$coef.MMLasso.ad},
  slts = function(y,x){robustHD::sparseLTS(x,y)$coefficients},

  cell_lasso_post = function(y,x){
    fit = regcell::sregcell_std(y = y, x = x, softbeta = TRUE, lambda_zeta = 1, penal = 1, penaldelta =0)
    return(c(fit$intercept_hat_post, fit$betahat_post))
  },

  lasso  = function(y,x){
    return(regcell::slm2(y,x)$betahat)
  }

)

load("datascreen_corh.RData")

# load("vignettes/datascreen_corh.RData")
y = datascreen[,1]
x = as.matrix(robustHD::robStandardize(datascreen[,-1]))


result <- foreach(mtd = 1:length(mtds),
                  .packages = c("lars", "robustHD", "robustbase" , "mmlasso",
                                "shootings", "cellWise", "regcell", "purrr", "robcovsel"))%:%
  foreach(p = c(100))%:%
  foreach(obs = 1:84)%dopar%{

    x = x[,1:p]
    ytrain = y[-obs]
    xtrain = as.matrix(x[-obs,])
    ytest = y[obs]
    xtest = x[obs,]

    betahat = mtds[[mtd]](ytrain,xtrain)

    res = ytest - betahat[1] - sum(xtest*betahat[-1])
    rst = c(mtd = names(mtds)[mtd], n = obs, p = p, res = res)

    rst
  }


save(result, file = "result_real_data_corh_0425.RData")






