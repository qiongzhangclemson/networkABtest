evaluate <- function(x, W, Z, rho) {
  ##################################################################
  # Function to compute evaluation criterion for a design under car model assumption
  # Inputs:
  # x: design vector allocation (-1, 1)
  # W: network adjacency matrix 
  # Z: covariates
  # rho: network correlation parameter
  #
  # Outputs:
  # Precentage of Improvement to Precision and others
  #################################################################
  
  # compute statistics
  n <- nrow(W)
	nedge <- rowSums(W)
  nedge.total <- sum(W)
  p <- ncol(Z)
  Erand <- -999
  T1 <- -999
  T2 <- -999
  ImPre <- -999
  obj <- -999
  flag <- length(x)
  
  if(flag!=0) {
    Rmat <- Diagonal(n, x =nedge)-rho*W
    Fmat <- cbind(1, Z)
	  Fsig <- t(Fmat) %*% Rmat
	  xFsig <- Fsig %*% x
    Sigma <- Fsig %*% Fmat
          
    # compute the objective
    T1 <- rho * as.numeric(t(x)%*% W %*% x)
    T2 <- as.numeric(t(xFsig) %*% solve(Sigma, xFsig))
    obj <- nedge.total-T1-T2
          
    # compute the expected objective of random designs
    KF <- t(Fsig) %*% solve(Sigma, Fsig)
    K <- Rmat - KF
    Erand <- sum(diag(K))
    ET1 <- sum(diag(W))
    ET2 <- sum(diag(KF))
    
    # compute improvementment of precision 
    ImPre <- 1-Erand/obj
  }
  return(c(rho, ImPre, obj, Erand, T1, T2, ET1, ET2, 1*(flag==0)))
}





