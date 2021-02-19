rm(list=ls())
require(Matrix)
require(tidyverse)
require(mclcar)
require(CARBayes)
require(spam)
require(spdep)
require(spatialreg)
require(expm)


source("opt.R")
source("gen_outcome.R")

set.seed(1)
rho_true_list <- c(0.1, 0.5, 0.9)
rho_true_list <- 0.5
nrep <- 100
cases <- list(
  n = c(50),
  net.p = 0.01, # network density
  p = 10,
  beta = 1,
  mymodel=6,
  rho_design=0.5,
  alpha=0.001
)

# mymodel: the true model
# 1: car model
# 2: linear model
# 3: sar model
# 4: car model with binary response
# 5: car model with gamma response
# 6: car model with rho varying from different links

cases.run <- expand.grid(cases)

for(i in 1:nrow(cases.run)) {
  n <- cases.run$n[i]
  net.p <- cases.run$net.p[i]
  beta <- cases.run$beta[i]
  p <- cases.run$p[i]
  mymodel <- cases.run$mymodel[i]
  rho_design <- cases.run$rho_design[i]
  alpha <- cases.run$alpha[i]
  
  
    
  # generate network data
  W <- gen_network(n, net.p)
  Z <- matrix(sample(c(-1, 1), n*p, replace=TRUE), n, p)
  
  
  # optimal design
  #x.opt <- opt.design0(W, alpha=0.6)
  # optimal design without network
  x.opt0 <- opt.design(W, Z, rho=rho_design, network=FALSE)
  
  # modified optimal design  (T1 as constraint)
  x.opt <- opt.design.modif2(alpha=alpha, W, Z, rho=rho_design)
  

  
  # optimal with linear network effect
  x.opt.linear <- opt.design.linear(W, Z, covariates = TRUE)$x
  
  
  # random design 
  rand.designs <- replicate(nrep, sample(c(rep(1, n/2), rep(-1, n/2))))
  
  collect.results <- NULL
  for(rho_true in rho_true_list) {
    pars.true <- c(rho_true, beta)
    betas <- NULL
    
    

    # opt design results
    betas  <- gen_outcome(x.opt, W, pars.true, Z=Z, model=mymodel)
    betas  <- cbind(betas, gen_outcome(x.opt0, W, pars.true, Z=Z, model=mymodel))
    betas  <- cbind(betas, gen_outcome(x.opt.linear, W, pars.true, Z=Z, model=mymodel))
    
    #betas <- matrix(betas, ncol=1)
    #results <-cbind(apply(betas, 2, mean)-beta, 
    #                apply(betas, 2, var))
    #results <- cbind(results, results[,1]^2+results[,2])
    #print(results)
    
    # random design results
    for(j in 1:ncol(rand.designs)) {
      try(betas <- cbind(betas, gen_outcome(rand.designs[,j], W, pars.true, Z=Z, 
                                                  model=mymodel)), TRUE)
      print(j)
    }
    est.beta <- beta
    if(mymodel==4) {
      est.beta <- mean(apply(betas, 2, mean))
    }
    
    results <-cbind(apply(betas, 2, mean)-est.beta, 
                                apply(betas, 2, var))
    results <- cbind(results, results[,1]^2+results[,2])
    
        
    results <- cbind(n, net.p, rho_true, results)
    results  <- data.frame(results)
    names(results) <- c("n", "network_density", "network_cor",
                          "bias", "var", "mse")
    results$method <- c("Optim", "Optim_No_Network", "Optim_linear", paste("Rand", 1:(nrow(results)-3)))
    results$quantile <- rank(results$mse)/((3+nrep)+1)
    results$rho_true <- rho_true
    results$model <- mymodel
    results$rho_design <- rho_design
    results$alpha <- alpha
    collect.results <- rbind(collect.results, results)
  }
  save(collect.results, file=paste("n=",n, "net.p=",net.p, "model=",mymodel,
                           "alpha", alpha,"rho_design", rho_design, 
                           "evaluate_simu_robust.RData"))
  re <- collect.results %>%
    filter(method %in% c("Optim","Optim_No_Network","Optim_linear"))
  print(re)
}





