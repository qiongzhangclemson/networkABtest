opt.design <- function(Z, W=NA, alpha=0.001, rho=0.5, timelimit=500) {
  ##################################################################
  # Function for optimal design with and without network
  # Inputs:
  # Z: covariate matrix, intercept not included
  # W: network adjacency matrix (NOT USED FOR network==F)
  # alpha: quantile parameter (NOT USED FOR network==F)
  # rho: network correlation parameter  (NOT USED FOR network==F)
  # timelimit (in sec): timelime in solving the optimization problem
  #
  # Outputs:
  # optimal allocation (-1, 1)
  #################################################################
  require(gurobi)
  require(Matrix)
  
  # common statistics
  Fmat <- cbind(1, Z)
  
  
  # setup global parameters
  params <- list(MIPGap=0.01, TimeLimit=timelimit)
  
  
  # setup optimization model
  model <- list()
  
  # setup linear constraint 
  model$A     <- matrix(1, 1, n)
  model$rhs   <- round(n/2) # to compile with the case that n is odd
  model$sense <- c("=")
  
  # set up decision variable type
  model$vtype <- rep('B', n)
  
  # with or without network
  
  network <- !is.na(sum(W))
  
  if(network == TRUE) {
    # network statistics
    n <- nrow(W)
    nedge <- rowSums(W)
    
    # set up quadratic objective
    Sigma.net <- Diagonal(n, x = nedge)-rho*W
    Fsig <- t(Fmat) %*% Sigma.net
    Sigma <- Fsig %*% Fmat		
    Q <- as.matrix(t(Fsig) %*% solve(Sigma, Fsig))
    model$Q     <- Q
    model$obj   <- -rowSums(Q)
    
    # set up quadratic constraint
    model$quadcon <- list()
    qc <- list()
    qc$Qc <- 4*W
    qc$q <- -4*rowSums(W)
    qc$rhs <- sqrt(sum(nedge))*qnorm(alpha)-sum(W)
    model$quadcon[[1]] <- qc
    
  } 

  if(network == FALSE) {
    # set up quadratic objective
    Finv <- solve(t(Fmat) %*% Fmat, t(Fmat))
    Q <- as.matrix(Fmat %*% Finv)
    model$Q     <- Q
    model$obj   <- -rowSums(Q)
  }
  
  # solve the problem with gurobi
  result <-gurobi(model, params)  
  return(2*result$x-1)
}






