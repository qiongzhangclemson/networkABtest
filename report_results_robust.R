rm(list=ls())
require(tidyverse)
library(xtable)

nrep <- 100
cases <- list(
  p = rep(c(5, 20), each=25),  
  beta = 1,
  rho_design = 0.5,
  alpha = 0.001
)
cases.run <- expand.grid(cases)


results <- NULL
for(i in 1:nrow(cases.run)) {
  n <- cases.run$n[i]
  net.p <- cases.run$net.p[i]
  beta <- cases.run$beta[i]
  p <- cases.run$p[i]
  mymodel <- cases.run$mymodel[i]
  rho_design <- cases.run$rho_design[i]
  alpha <- cases.run$alpha[i]
  
  load(file=paste("case", i, "p=", p, "robust.RData"))
  
  
  results0 <- collect.results
  results0$case <- i-floor(i/25)*25+1
  results0$p <- p
  results <- rbind(results, results0)
} 


results %>%
  mutate(
    p = factor(paste("p =", p), levels=c("p = 5", "p = 20")),
    method = ifelse(method=="Optim", "Yes", "No")
  ) %>%
  ggplot(aes(x=method, y=quantile)) +
  geom_boxplot() +
  facet_grid(.~p) +
  ylab("Percentiles of MSE")+
  xlab("Network")+
  theme_bw()

ggsave("mse.pdf", height=3, width=6)