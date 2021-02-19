rm(list=ls())
require(tidyverse)

require(xtable)

n <- 50
net.p <- 0.01
p <- 10
beta <- 1
mymodel <- 6

# mymodel: the true model
# 1: car model
# 2: linear model
# 3: sar model
# 4: car model with binary response
# 5: car model with gamma response
# 6: car model with varying rho

rho_design <- 0.5
alpha <- 0.001


load(file=paste("n=",n, "net.p=",net.p, "model=",mymodel,
                                 "alpha", alpha,"rho_design", rho_design, 
                                 "evaluate_simu_robust.RData"))

re <- collect.results %>%
  filter(method %in% c("Optim","Optim_No_Network","Optim_linear")) %>%
  select(-n, -network_density, -model, -rho_design, -alpha, -network_cor) %>%
  mutate(
    method=factor(method, labels=c("with CAR Network (Bayesian)", "with Linear Network", "w/o Network"))
  )

print(xtable(re[c(6,4, 1, 2, 3, 5)], digits=c(1, 1,0, 3,3,3,3)),
      include.rownames=FALSE)