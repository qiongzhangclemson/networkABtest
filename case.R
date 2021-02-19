rm(list=ls())

source("evaluate.R")
source("opt.R")
require(expm)
require(gurobi)
require(Matrix)


rhos <- c(0.1, 0.5, 0.9)

collect.results <- NULL

for(n in c(1000, 2000, 5000)) {
  load("HR.Rdata")
  p <- 30
  Z <- X[1:n, 1:p]
  W <- as.matrix(W[1:n, 1:n])
  rho_design <- 0.5
  
  # optimal design without network
  x.opt0 <- opt.design(W, Z, rho=rho_design, network=FALSE)
  
  # modified optimal design  (T1 as constraint)
  x.opt <- opt.design.modif2(alpha=0.00001, W, Z, rho=rho_design)
  
  for(rho_true in rhos) {
    x.opt.true <- opt.design.modif2(alpha=0.00001, W, Z, rho=rho_true)
    var_re_opt0 <- evaluate(x.opt0, W, Z, rho=rho_true)
    var_re_opt <- evaluate(x.opt, W, Z, rho=rho_true)
    var_re_opt_true <- evaluate(x.opt.true, W, Z, rho=rho_true)
    
    # attach results
    re.vector <- cbind(n, p, rho_true, rho_design, 
                       rbind(var_re_opt0, var_re_opt))
    re.vector <- rbind(re.vector, c(n, p, rho_true, rho_true, var_re_opt_true))
    re.vector <- data.frame(re.vector)
    names(re.vector) <- c("n", "p", "network_cor", 
                          "design_network_cor", "Var_Re", "Obj", "Rand","T1", "T2")
    re.vector$method <- c("Optim_no_network","Optim","Optim_true")
    collect.results <- rbind(collect.results, re.vector)
  }
}

save(collect.results, file="sub.RData")


