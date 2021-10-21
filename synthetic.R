rm(list=ls())
require(Matrix)
require(tidyverse)
require(gurobi)

source("evaluate.R")
source("opt.R")


# Performance with Synthetic Networks
microrep <- 10
rho_true_list <- c(0.1, 0.3, 0.5, 0.7, 0.9)
rho_design <- 0.5

cases.small <- list(
  n = c(50, 100),
  p = c(5, 10),
  net.p = seq(0.02, 0.1, 0.02), # network density
  alpha = c(0.1, 0.01, 0.001, 0.0001)
)

cases.run <- expand.grid(cases.small)


cases.large <- list(
  n = c(500, 1000),
  p = c(10),
  net.p = 0.02, # network density
  alpha = 0.01
)

cases.run <- rbind(cases.run, expand.grid(cases.large))

set.seed(1)

collect.results <- NULL
for(i in 1:nrow(cases.run)) {
  n <- cases.run$n[i]
  p <- cases.run$p[i]
  net.p <- cases.run$net.p[i]
  alpha <- cases.run$alpha[i]
    
  
  results <- NULL
  for(j in 1:microrep) {
    # generate data
    Z <- matrix(sample(c(-1, 1), n*p, replace=TRUE), n, p)
    W <- sample(c(0, 1), n*n, replace=TRUE, prob=c(1-net.p/2,net.p/2))
    W <- matrix(W, n, n)
    W <- 1*((W+t(W))>0)
    nedge.total <- sum(W)
    diag(W) <- 0
    
    # optimal design with network
    x.opt <- opt.design(Z, W, alpha=alpha, rho = rho_design)
    
    # optimal design without network
    x.opt0 <- opt.design(Z)
    
    for(rho_true in rho_true_list) {
      x.opt.true <- opt.design(Z, W, alpha=alpha, rho = rho_true)
      opt.re <- c(1,evaluate(x.opt, W, Z, rho=rho_true))
      opt.re_true <- c(1,evaluate(x.opt.true, W, Z, rho=rho_true))
      opt.re0 <- c(0, evaluate(x.opt0, W, Z, rho=rho_true))
      results <- rbind(results, opt.re, opt.re_true, opt.re0)
    }
  }
  collect.results <- rbind(collect.results, cbind(n, p, net.p, alpha, rho_design, results))
}

collect.results$rho_design <- rep(c(0.5, 1, 0), 4100) # rho_design 0 referring to no network
# rho_design 1 referring to rho_design=rho_true.

collect.results <- data.frame(collect.results)
names(collect.results) <- c("n", "p", "net.density","alpha", "rho_design","Network", "rho_true", 
                            "ImprovObj", "obj", "E(obj)", "T1", "T2",  "E(T1)", "E(T2)","Error")



save(collect.results, file="syn.RData")


# robust to alpha

collect.results %>%
  filter(p == 10, Network == 1, net.density == 0.08, n %in% c(50, 100)) %>%
  mutate(
    n = factor(paste("n =", n), levels=c("n = 50", "n = 100")),    
    alpha=factor(paste("alpha =", alpha), levels=c("alpha = 0.1", "alpha = 0.01", "alpha = 0.001", "alpha = 1e-04"))
  ) %>%
  filter(rho_design == 0.5, net.density == 0.08) %>%
  ggplot(aes(x=as.factor(rho_true), y=ImprovObj*100)) +
  geom_boxplot()+
  facet_grid(n~alpha)+
  xlab("True Correlation Parameter")+
  ylab("Improvement in Precision (%)")+
  theme_bw()

ggsave("small_alpha.pdf", height=6, width=12)

# robust to rho

collect.results %>%
  filter(Network == 1, n %in% c(50, 100),net.density == 0.08, alpha == 0.001) %>%
  select(n, p, rho_design, rho_true, ImprovObj, `E(obj)`) %>%
  mutate(
    n = factor(paste("n =", n), levels=c("n = 50", "n = 100")),    
    p = factor(paste("p =", p), levels=c("p = 5", "p = 10"))   
  ) %>%
  group_by(n, p, `E(obj)`, rho_true) %>%
  summarize(
    Diff=first(ImprovObj)-last(ImprovObj)
  ) %>%
  ggplot(aes(x=as.factor(rho_true), y=Diff*100)) +
  geom_boxplot()+
  facet_grid(n~p)+
  xlab("True Correlation Parameter")+
  ylab("Difference of Improvement in Precision")+
  theme_bw()

ggsave("small_rho.pdf", height=6, width=6)


# robust to net.density

collect.results %>%
  filter(p == 10, n %in% c(50, 100), rho_design!=1, net.density>0.02) %>%
  mutate(
    n = factor(paste("n =", n), levels=c("n = 50", "n = 100")),
    Network = ifelse(Network==0, "No", "Yes")
  ) %>%
  ggplot() +
  geom_boxplot(aes(x=as.factor(rho_true), y=ImprovObj*100, fill=Network, color=Network))+
  facet_grid(n~paste("network density = ",net.density))+
  xlab("True Correlation Parameter")+
  ylab("Improvement in Precision (%)")+
  theme_bw()+
  theme(legend.position = "bottom")

ggsave("small_netp.pdf", height=6, width=12)

# robust to n
collect.results %>%
  filter(p == 10, rho_design!=1, net.density==0.02) %>%
  mutate(
    n = factor(paste("n =", n), levels=c("n = 50", "n = 100", "n = 500", "n = 1000")),
    Network = ifelse(Network==0, "No", "Yes")
  ) %>%
  ggplot() +
  geom_boxplot(aes(x=as.factor(rho_true), y=ImprovObj*100, fill=Network, color=Network))+
  facet_wrap(~n, scale="free_y")+
  xlab("True Correlation Parameter")+
  ylab("Improvement in Precision (%)")+
  theme_bw()+
  theme(legend.position = "bottom")

ggsave("n.pdf", height=3, width=3)


collect.results %>%
  filter(p == 10, rho_design!=1, net.density==0.02) %>%
  mutate(
    n = factor(paste("n =", n), levels=c("n = 50", "n = 100", "n = 500", "n = 1000")),
    Network = ifelse(Network==0, "No", "Yes")
  ) %>%
  ggplot() +
  geom_boxplot(aes(x=as.factor(rho_true), y=-T1, fill=Network, color=Network))+
  facet_wrap(~n, scale="free_y")+
  xlab("True Correlation Parameter")+
  ylab("Improvement in T1")+
  theme_bw()+
  theme(legend.position = "bottom")

ggsave("t1.pdf", height=3, width=3)


collect.results %>%
  filter(p == 10, rho_design!=1, net.density==0.02) %>%
  mutate(
    n = factor(paste("n =", n), levels=c("n = 50", "n = 100", "n = 500", "n = 1000")),
    Network = ifelse(Network==0, "No", "Yes")
  ) %>%
  ggplot() +
  geom_boxplot(aes(x=as.factor(rho_true), y=`E(T2)`-T2, fill=Network, color=Network))+
  facet_wrap(~n, scale="free_y")+
  xlab("True Correlation Parameter")+
  ylab("Improvement in T2")+
  theme_bw()+
  theme(legend.position = "bottom")

ggsave("t2.pdf", height=3, width=3)

