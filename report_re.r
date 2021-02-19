require(tidyverse)
require(xtable)

## Simulated network: section 6.2: small networks #######


cases <- list(
  n = c(50, 100),
  p = c(5, 10),
  net.p = c(0.01) # network density
)

cases.run <- expand.grid(cases)
re.all <- NULL
for(i in 1:nrow(cases.run)) {
  n <- cases.run$n[i]
  p <- cases.run$p[i]
  net.p <- cases.run$net.p[i]
  rho_design <- 0.5
  for(rho_true in c(0.1, 0.5, 0.9)) {
     load(paste("n=",n, "p=",p,"net.p=",net.p, "rho_true", rho_true, "simu.RData"))
     re.all <- rbind(re.all, re.vector)
  }
}


row.names(re.all) <- NULL
results <- re.all %>%
   select(-design_network_cor,  -Obj, -T1, -T2, -Rand) %>%
   arrange(method, n, p, network_cor) 


results %>%
  #filter(network_density == 1e-04) %>%
  mutate(
    method = factor(method, labels=c("with CAR Network (Bayesian)", "with Linear Network", "w/o Network",
                                     "with CAR Network (Exact)")),
    p=factor(paste("p =", p), levels=c("p = 5", "p = 10")),
    n=factor(paste("n =", n), levels=c("n = 50", "n = 100")),
  ) %>%
  ggplot(aes(network_cor, Var_Re, group=method, color=method, linetype=method))+
  geom_point()+
  geom_line()+
  facet_grid(n~p)+
  ylab("variance reduction")+
  xlab("network correlation of the true model")+
  theme_bw()+
  theme(legend.position = "bottom")
ggsave("small.pdf", width=7.5, height=5)

## Simulated network: section 6.2: large networks #######


cases <- list(
  n = c(1000, 2000),
  p = c(20, 50),
  net.p = c(0.0001) # network density
)


cases.run <- expand.grid(cases)
re.all <- NULL
for(i in 1:nrow(cases.run)) {
  n <- cases.run$n[i]
  p <- cases.run$p[i]
  net.p <- cases.run$net.p[i]
  rho_design <- 0.5
  for(rho_true in c(0.1, 0.5, 0.9)) {
    load(paste("n=",n, "p=",p,"net.p=",net.p, "rho_true", rho_true, "simu.RData"))
    re.all <- rbind(re.all, re.vector)
  }
}


row.names(re.all) <- NULL
results <- re.all %>%
  select(-design_network_cor,  -Obj, -T1, -T2, -Rand) %>%
  arrange(method, n, p, network_cor) 


results %>%
  mutate(
    method = factor(method, labels=c("with CAR Network (Bayesian)", "with Linear Network", "w/o Network",
                                     "with CAR Network (Exact)")),
    p=factor(paste("p =", p), levels=c("p = 20", "p = 50")),
    n=factor(paste("n =", n), levels=c("n = 1000", "n = 2000")),
  ) %>%
  ggplot(aes(network_cor, Var_Re, group=method, color=method, linetype=method))+
  geom_point()+
  geom_line()+
  facet_grid(n~p)+
  ylab("variance reduction")+
  xlab("network correlation of the true model")+
  theme_bw()+
  theme(legend.position = "bottom")

ggsave("large.pdf", width=7.5, height=5)






#### case study: Table########
load("sub.RData")
row.names(collect.results) <- NULL


results <- collect.results %>%
  select(-T1, -T2, -design_network_cor, -Obj, -Rand) %>%
  arrange(method, n, p, network_cor) 


report_re <- results %>%
  spread(method, Var_Re)

# copy paste the following for Table 2
print(xtable(report_re[c(1,2,3,5,4,6)], digits=c(1, 0, 0, 1, 3, 3, 3)),
      include.rownames=FALSE)



