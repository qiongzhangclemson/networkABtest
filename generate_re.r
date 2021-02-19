rm(list=ls())
require(Matrix)
require(tidyverse)

source("evaluate.R")
source("opt.R")
require(expm)


#cases <- list(
#  n = c(1000, 2000),
#  p = c(5, 20, 50),
#  net.p = c(0.0001, 0.0005) # network density
#)

cases <- list(
  n = c(50, 100),
  p = c(5, 10),
  net.p = c(0.01) # network density
)


cases.run <- expand.grid(cases)



for(i in 1:nrow(cases.run)) {
  n <- cases.run$n[i]
  p <- cases.run$p[i]
  net.p <- cases.run$net.p[i]
  rho_design <- 0.5
  
  # generate data
  set.seed(2)
  Z <- matrix(sample(c(-1, 1), n*p, replace=TRUE), n, p)
  W <- sample(c(0, 1), n*n, replace=TRUE, prob=c(1-net.p/2,net.p/2))
  W <- matrix(W, n, n)
  W <- 1*((W+t(W))>0)
  nedge.total <- sum(W)
  diag(W) <- 0
  
  # generate correlated data
  # Lmat <- diag(apply(W, 1, sum)) - W
  # Z <- eigen(Lmat)$vectors[,1:p]
  
  # optimal design without network
  x.opt0 <- opt.design(W, Z, rho=rho_design, network=FALSE)
  
  # optimal with linear network effect
  x.opt.linear <- opt.design.linear(W, Z, covariates = FALSE)$x

  
  #time1 <-system.time(x.opt1 <- opt.design(W, Z, rho=rho_design))
  
  
  # modified optimal design  (T1 as constraint)
  time2 <-system.time(x.opt <- opt.design.modif2(alpha=0.001, W, Z, rho=rho_design))
  
  for(rho_true in c(0.1, 0.5, 0.9)) {
    x.opt.true <- opt.design.modif2(alpha=0.001, W, Z, rho=rho_true)
    
    var_re_opt0 <- evaluate(x.opt0, W, Z, rho=rho_true)
    var_re_opt <- evaluate(x.opt, W, Z, rho=rho_true)
    var_re_opt_true <- evaluate(x.opt.true, W, Z, rho=rho_true)
    var_re_linear <- evaluate(x.opt.linear, W, Z, rho=rho_true)
    
    # attach results
    re.vector <- cbind(n, p, net.p, rho_true, rho_design, 
                     rbind(var_re_opt0, var_re_opt, var_re_linear))
    re.vector <- rbind(re.vector, c(n, p, net.p, rho_true, rho_true, var_re_opt_true))
    re.vector <- data.frame(re.vector)
    names(re.vector) <- c("n", "p", "network_density", "network_cor", 
                          "design_network_cor", "Var_Re", "Obj","Rand", "T1", "T2")
    re.vector$method <- c("Optim_no_network","Optim", "Optim_linear", "Optim_true")
    print(re.vector)
    save(re.vector, file=paste("n=",n, "p=",p,"net.p=",net.p, 
                               "rho_true", rho_true, "simu.RData"))
  }
}





