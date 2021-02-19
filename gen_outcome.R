gen_network <- function(n, net.p, ego=FALSE) {
  W <- sample(c(0, 1), n*n, replace=TRUE, prob=c(1-net.p/2,net.p/2))
  W <- matrix(W, n, n)
  W <- 1*((W+t(W))>0)
  nedge.total <- sum(W)
  diag(W) <- 0
  for(j in 1:n) {
    if(sum(W[j,])==0) {
      id <- sample(c(1:(j-1),(j+1):n), 1)
      W[j,id] <- 1
      W[id, j] <- 1
    }
  }
  if(ego==TRUE) {
    W[1, -1] <- 1
    W[-1, 1] <- 1
  }
  return(W)
}

gen_outcome <- function(x, W, pars.true, Z=NULL, iter=100, model=1) {
  n <- nrow(W)
  rho <- pars.true[1]
  beta <- pars.true[2]
  z <- apply(Z*beta, 1, sum)
  betas <- rep(0, iter)
  if(model==1) {
    for(it in 1:iter) {
      delta <- CAR.simWmat1(rho, 1, W)
      # these lines are useful if you would like to generate it for different distributions
      #delta <- pnorm(delta)
      #delta <- ifelse(delta>0.999, 0.999, delta)
      y <- delta + x*beta+z+1
      data.splm <- data.frame(y=y, cbind(x, Z))
      listW <- mat2listw(W)
      fit.splm <- spautolm(y~., data = data.splm, listw=listW, family = "CAR")
      betas[it] <- fit.splm$fit$coefficients[2]
    }
  }
  if(model==2) {
    for(it in 1:iter) {
      y <- rnorm(length(x)) + x*beta+z+1+(row_sums(W)+W %*% x)*rho
      data.lm <- data.frame(y=y, cbind(x, Z, row_sums(W), W %*% x))
      fit.lm <-lm(y~., data = data.lm)
      betas[it] <- fit.lm$coefficients[2]
    }
  }
  if(model==3) {
    for(it in 1:iter) {
      delta <- SAR.simWmat1(rho, 1, W)
      # these lines are useful if you would like to generate it for different distributions
      #delta <- pnorm(delta)
      #delta <- ifelse(delta>0.999, 0.999, delta)
      y <- delta + x*beta+z+1
      data.splm <- data.frame(y=y, cbind(x, Z))
      listW <- mat2listw(W)
      fit.splm <- spautolm(y~., data = data.splm, listw=listW, family = "SAR")
      betas[it] <- fit.splm$fit$coefficients[2]
    }
  }
  if(model==4) {
    for(it in 1:iter) {
     pars <- c(rho, 1, beta, rep(beta, ncol(Z)))
     delta <- CAR.simGLM1(method = "binom", W, pars=pars, Xs = cbind(x, Z), n.trial = 1) 
     betas[it]=mean(delta$y[x==1])-mean(delta$y[x==-1])
    }
  }
  if(model==5) {
    for(it in 1:iter) {
      delta <- CAR.simWmat1(rho, 1, W)
      # these lines are useful if you would like to generate it for different distributions
      delta <- pnorm(delta)
      delta <- ifelse(delta>0.999, 0.999, delta)
      y <- qgamma(delta, 1) + x*beta+z+1
      data.splm <- data.frame(y=y, cbind(x, Z))
      listW <- mat2listw(W)
      fit.splm <- spautolm(y~., data = data.splm, listw=listW, family = "CAR")
      betas[it] <- fit.splm$fit$coefficients[2]
    }
  }
  if(model==6) {
    rhos <- matrix(runif(n*n, rho-0.4,rho+0.4),n,n)
    rhos <- (rhos+t(rhos))/2
    Wt <- W*rhos
    for(it in 1:iter) {
      delta <- CAR.simWmat1(0.99, 1, Wt)
      y <- delta + x*beta+z+1
      data.splm <- data.frame(y=y, cbind(x, Z))
      listW <- mat2listw(W)
      fit.splm <- spautolm(y~., data = data.splm, listw=listW, family = "CAR")
      betas[it] <- fit.splm$fit$coefficients[2]
    }
  }
  return(betas)
}







CAR.simWmat1 <- function (rho, prec, W) 
{
  n <- nrow(W)
  I <- diag.spam(1, n)
  Q <- as.spam(prec * (diag(apply(W, 2, sum)) - rho * W))
  cholR <- chol.spam(Q, pivot = "MMD", memory = list(nnzcolindices = 6.25 * 
                                                       n))
  X <- backsolve(cholR, rnorm(n))
  X
}

SAR.simWmat1 <- function (rho, prec, W) 
{
  n <- nrow(W)
  I <- diag.spam(1, n)
  Q <- as.spam(prec * (I - rho * W)%*%(I-rho*W))
  cholR <- chol.spam(Q, pivot = "MMD", memory = list(nnzcolindices = 6.25 * 
                                                       n))
  X <- backsolve(cholR, rnorm(n))
  X
}



CAR.simGLM1 <- function (method = c("binom", "poisson"), W, n = NULL, pars, 
          Xs = NULL, n.trial = 1) 
{
  rho <- pars[1]
  sigma <- pars[2]
  if (is.null(n)) {
    X <- CAR.simWmat1(rho, 1/sigma, W)
    Z.car <- list(X = X, W = W)
  }
  else {
    Z.car <- CAR.simTorus(n[1], n[2], rho, 1/sigma)
  }
  if (length(pars) > 2) {
    beta <- pars[-c(1, 2)]
    eta <- Xs %*% beta + Z.car$X
  }
  else {
    eta <- Z.car$X
  }
  if (method == "binom") {
    ps <- exp(eta)/(1 + exp(eta))
    Emean <- ps
    N <- length(ps)
    Y <- rbinom(N, n.trial, prob = ps)
  }
  else {
    lambda <- exp(eta)
    Emean <- lambda
    N <- length(lambda)
    Y <- rpois(N, lambda)
  }
  list(rho = rho, sigma = sigma, beta = beta, y = Y, covX = Xs, 
       W = Z.car$W, Z.true = Z.car$X, eta = eta, Emean = Emean, 
       n.trial = n.trial)
}