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
load("HR.Rdata")
W0 <- W
Z0 <- X

set.seed(1)
nrep <- 10
cases <- list(
  p = rep(c(5, 20), each=25),
beta = 1, #treatment effect
  rho_design = 0.5,
  alpha = 0.001
)

cases.run <- expand.grid(cases)

for(i in 1:nrow(cases.run)) {
  beta <- cases.run$beta[i]
  p <- cases.run$p[i]
  rho_design <- cases.run$rho_design[i]
  alpha <- cases.run$alpha[i]
  
  
    
  # generate network data
  id <- sample(1:nrow(W0), 2000)
  W <- as.matrix(W0[id, id])
  Z <- as.matrix(Z0[id, ])
  id1 <- apply(W, 1, sum)!=0 # remove isolated users
  W <- W[id1, id1]
  Z <- Z[id1,]
  Z <- Z[ , apply(Z, 2, sd)>0]
  Z <- Z[,1:p]
  n <- nrow(W)
  net.p <- mean(W)
  
  # optimal design without network
  x.opt0 <- opt.design(Z)
  
  # optimal design
  x.opt <- opt.design(Z, W, rho=rho_design)
  
  
  
  # random design 
  rand.designs <- replicate(nrep, sample(c(rep(1, n/2), rep(-1, n/2))))
  
  collect.results <- NULL
  pars.true <- c(0.5, beta)
  betas <- NULL
    
    

  # opt design results
  betas  <- gen_outcome(x.opt, W, pars.true, Z=Z)
  betas  <- cbind(betas, gen_outcome(x.opt0, W, pars.true, Z=Z))

  for(j in 1:ncol(rand.designs)) {
      try(betas <- cbind(betas, gen_outcome(rand.designs[,j], W, pars.true, Z=Z)), TRUE)
      print(j)
  }
    
  est.beta <- beta
  results <-cbind(apply(betas, 2, mean)-est.beta, apply(betas, 2, var))
  results <- cbind(results, results[,1]^2+results[,2])
    
        
  results <- cbind(n, net.p, results)
  results  <- data.frame(results)
  names(results) <- c("n", "network_density", "bias", "var", "mse")
  results$method <- c("Optim", "Optim_No_Network", paste("Rand", 1:(nrow(results)-2)))
  results$quantile <- rank(results$mse)/((2+nrep)+1)
  results$rho_design <- rho_design
  results$alpha <- alpha
  collect.results <- rbind(collect.results, results)
  
  save(collect.results, file=paste("case", i, "p=", p, "robust.RData"))
  
  re <- collect.results %>%
      filter(method %in% c("Optim","Optim_No_Network"))
  print(re)
}





