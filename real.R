rm(list=ls())
require(Matrix)
require(tidyverse)
require(gurobi)

source("evaluate.R")
source("opt.R")


# Performance with Large Networks

load("HR.RData")
W0 <- W
microrep <- 10 # generate different dataset
rho_true_list <- c(0.1, 0.3, 0.5, 0.7, 0.9)
rho_design <- 0.5

cases <- list(
  nt = c(2000, 3000),
  p = c(5, 10, 20)
)

cases.run <- expand.grid(cases)

set.seed(1)

collect.results <- NULL
for(i in 1:nrow(cases.run)) {
  nt <- cases.run$nt[i]
  p <- cases.run$p[i]
  net.p <- cases.run$net.p[i]
  alpha <- cases.run$alpha[i]
    
  
  results <- NULL
  for(j in 1:microrep) {
    # generate data
    id <- sample(1:nrow(W0), nt)
    W <- as.matrix(W0[id, id])
    Z <- as.matrix(X[id,])
    id1 <- apply(W, 1, sum)!=0
    W <- W[id1, id1]
    Z <- Z[id1,]
    Z <- Z[ , apply(Z, 2, sd)>0]
    Z <- Z[,1:p]
    n <- nrow(W)
    net.p <- mean(W)
    #Z <- matrix(sample(c(-1, 1), n*p, replace=TRUE), n, p)
    
    # optimal design with network
    x.opt <- opt.design(Z, W, alpha=0.001, rho = rho_design)
    
    # optimal design without network
    x.opt0 <- opt.design(Z)
    
    for(rho_true in rho_true_list) {
      opt.re <- c(1,evaluate(x.opt, W, Z, rho=rho_true))
      opt.re0 <- c(0, evaluate(x.opt0, W, Z, rho=rho_true))
      results <- rbind(results, opt.re, opt.re0)
    }
  }
  collect.results <- rbind(collect.results, cbind(n, p, net.p, rho_design, results))
}



collect.results <- data.frame(collect.results)
names(collect.results) <- c("n", "p", "net.density", "rho_design","Network", "rho_true", 
                            "ImprovObj", "obj", "E(obj)", "T1", "T2", "E(T1)", "E(T2)", "Error")



save(collect.results, file="real_var.RData")

collect.results %>%
  mutate(
    n = factor(paste("n approx.", round(n, -3)), levels=c("n approx. 1000", "n approx. 2000")),
    p = factor(paste("p =", p), levels=c("p = 5", "p = 10", "p = 20")),
    Network = ifelse(Network==0, "No", "Yes")
  ) %>%
  ggplot(aes(x=as.factor(rho_true), y=ImprovObj*100, fill=Network, color=Network)) +
  geom_boxplot()+
  facet_grid(n~p)+
  xlab("True Correlation Parameter")+
  ylab("Improvement in Precision (%)")+
  theme_bw()+
  theme(legend.position = "bottom")

ggsave("real.pdf", width=9, height=7)


