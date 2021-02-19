rm(list=ls())
require(Matrix)
require(tidyverse)

source("evaluate.R")
source("opt.R")
require(expm)


cases <- list(
  n = c(1000),
  p = c(5, 20),
  net.p = c(0.0003), # network density
  alpha = c(0.00001, 0.0001, 0.001, 0.01)
  #rho = c(0.1, 0.5, 0.9)
)

cases.run <- expand.grid(cases)
cases.run <- cases.run %>%
  filter(
    p<n
  )

results <- NULL
for(i in 1:nrow(cases.run)) {
  n <- cases.run$n[i]
  p <- cases.run$p[i]
  net.p <- cases.run$net.p[i]
  alpha <- cases.run$alpha[i]
  rho_design <- 0.5
  
  # generate data
  set.seed(2)
  Z <- matrix(sample(c(-1, 1), n*p, replace=TRUE), n, p)
  W <- sample(c(0, 1), n*n, replace=TRUE, prob=c(1-net.p/2,net.p/2))
  W <- matrix(W, n, n)
  W <- 1*((W+t(W))>0)
  nedge.total <- sum(W)
  diag(W) <- 0
  
  

  # modified optimal design  (T1 as constraint)
  time1 <-system.time(x.opt <- opt.design.modif2(alpha=alpha, W, Z, rho=rho_design))[3]

  
  for(rho_true in c(0.1, 0.5, 0.9)) {
    var_re_opt <- evaluate(x.opt, W, Z, rho=rho_true)
    
    # attach results
    re.vector <- cbind(n, p, alpha, rho_true, var_re_opt[1], var_re_opt[4])
    results <- rbind(results, re.vector)
  }
}

results <- data.frame(results)
names(results) <- c("n", "p", "alpha", "network_cor", 
                      "Var_Re", "T1")

save(results, file=paste("alpha_simu.RData"))

results %>%
  select(-Var_Re) %>%
  spread(alpha, T1)


