gen_outcome <- function(x, W, pars.true, Z=NULL, iter=100) {
  n <- nrow(W)
  rho <- pars.true[1]
  beta <- pars.true[2]
  z <- apply(Z*beta, 1, sum)
  betas <- rep(0, iter)
  rhos <- runif(n, 0.1, 0.9)
  for(it in 1:iter) {
      delta <- CAR.simWmat1(rhos, 1, W)
      y <- delta + x*beta+z+1
      data.splm <- data.frame(y=y, cbind(x, Z))
      listW <- mat2listw(W)
      fit.splm <- spautolm(y~., data = data.splm, listw=listW, family = "CAR")
      betas[it] <- fit.splm$fit$coefficients[2]
  }
  return(betas)
}


CAR.simWmat1 <- function (rho, prec, W) 
{
  n <- nrow(W)
  Q <- as.spam(prec * (diag(apply(W, 2, sum)) - diag(sqrt(rho)) %*% W %*% diag(sqrt(rho))))
  cholR <- chol.spam(Q, pivot = "MMD", memory = list(nnzcolindices = 6.25 * 
                                                       n))
  X <- backsolve(cholR, rnorm(n))
  X
}





