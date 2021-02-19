opt.design <- function(W, Z, rho=0.5, network = TRUE, timelimit=500) {
  ##################################################################
  # Function for optimal design
  # Inputs:
  # W: adjacency matrix
  # Z: covariate matrix, intercept not included
  # rho: network correlation parameter
  #
  # Outputs:
  # optimal allocation (-1, 1)
  #################################################################
          require(gurobi)
          require(Matrix)
  
      if(network==FALSE) {
        Fmat <- cbind(1, Z)
        Finv <- solve(t(Fmat) %*% Fmat, t(Fmat))
        Q <- as.matrix(Fmat %*% Finv)
        
        
        model <- list()
        
        # set global parameters
        params <- list(MIPGap=0.01, TimeLimit=500)
        
        # set up linear constraint 
        model$A     <- matrix(1, 1, n)
        model$rhs   <- n/2
        model$sense <- c("=")
        
        # set up quadratic objective
        model$Q     <- Q
        model$obj   <- -rowSums(Q)
        
        # set up decision variable type
        model$vtype <- rep('B', n)
        
        # solve the problem and return the results
        result <-gurobi(model, params)
        return(2*result$x-1)
        
      } else {
			    n <- nrow(W)
			    nedge <- rowSums(W)
          Sigma.net <- Diagonal(n, x =nedge)-rho*W
          Fmat <- cbind(1, Z)
			    Fsig <- t(Fmat) %*% Sigma.net
          Sigma <- Fsig %*% Fmat
          Q <- as.matrix(rho*W + t(Fsig) %*% solve(Sigma, Fsig))
          
          
          model <- list()
          
          # set global parameters
          params <- list(MIPGap=0.01, TimeLimit=timelimit)
          
          # set up linear constraint 
          model$A     <- matrix(1, 1, n)
          model$rhs   <- n/2
          model$sense <- c("=")
          
          # set up quadratic objective
          model$Q     <- Q
          model$obj   <- -rowSums(Q)
          
          # set up decision variable type
          model$vtype <- rep('B', n)
          
          # solve the problem and return the results
          result <-gurobi(model, params)
          return(2*result$x-1)
      }
}




opt.design.modif <- function(alpha, W, Z, rho=0.5) {
  ##################################################################
  # Function for optimal design with constraints (T2 as constraint)
  # Inputs:
  # W: adjacency matrix
  # Z: covariate matrix, intercept not included
  # rho: network correlation parameter
  #
  # Outputs:
  # optimal allocation (-1, 1)
  #################################################################
  require(gurobi)
  require(Matrix)
  
  # network statistics
  n <- nrow(W)
  nedge <- rowSums(W)
  Sigma.net <- Diagonal(n, x =nedge)-rho*W
  Fmat <- cbind(1, Z)
  Fsig <- t(Fmat) %*% Sigma.net
  Sigma <- Fsig %*% Fmat		
  Sigma2 <- sqrtm(Fsig %*% t(Fsig))		
  T2mat <- Sigma2 %*% solve(Sigma, Sigma2)
  lambdas <- eigen(T2mat, only.values = TRUE)$values
  Qc <- as.matrix(t(Fsig) %*% solve(Sigma, Fsig))
  chi_up <- quantile_mix_chi(alpha, lambdas)
  
  
  
  model <- list()
  
  # set global parameters
  params <- list(MIPGap=0.01, TimeLimit=500)
  
  # set up linear constraint 
  model$A     <- matrix(1, 1, n)
  model$rhs   <- n/2
  model$sense <- c("=")
  
  # set up quadratic objective
  model$Q     <- W
  model$obj   <- -rowSums(W)
  
  # set up decision variable type
  model$vtype <- rep('B', n)
  
  # set up quadratic constraint
  model$quadcon <- list()
  qc <- list()
  qc$Qc <- 4*Qc
  qc$q <- -4*rowSums(Qc)
  qc$rhs <- chi_up-sum(Qc)
  model$quadcon[[1]] <- qc
  
  # solve the problem  
  result <-gurobi(model, params)
  return(2*result$x-1)
}

opt.design.modif2 <- function(alpha, W, Z, rho=0.5) {
  ##################################################################
  # Function for optimal design with constraints (T1 as constraint)
  # Inputs:
  # W: adjacency matrix
  # Z: covariate matrix, intercept not included
  # rho: network correlation parameter
  #
  # Outputs:
  # optimal allocation (-1, 1)
  #################################################################
  require(gurobi)
  require(Matrix)
  
  # network statistics
  n <- nrow(W)
  nedge <- rowSums(W)
  Sigma.net <- Diagonal(n, x =nedge)-rho*W
  Fmat <- cbind(1, Z)
  Fsig <- t(Fmat) %*% Sigma.net
  Sigma <- Fsig %*% Fmat		
  Sigma2 <- sqrtm(Fsig %*% t(Fsig))		
  T2mat <- Sigma2 %*% solve(Sigma, Sigma2)
  lambdas <- eigen(T2mat, only.values = TRUE)$values
  Q <- as.matrix(t(Fsig) %*% solve(Sigma, Fsig))
  z_alpha <- qnorm(alpha)
  
  
  model <- list()
  
  # set global parameters
  params <- list(MIPGap=0.01, TimeLimit=500)
  
  # set up linear constraint 
  model$A     <- matrix(1, 1, n)
  model$rhs   <- n/2
  model$sense <- c("=")
  
  # set up quadratic objective
  model$Q     <- Q
  model$obj   <- -rowSums(Q)
  
  # set up decision variable type
  model$vtype <- rep('B', n)
  
  # set up quadratic constraint
  model$quadcon <- list()
  qc <- list()
  qc$Qc <- 4*W
  qc$q <- -4*rowSums(W)
  qc$rhs <- sqrt(sum(nedge))*z_alpha-sum(W)
  model$quadcon[[1]] <- qc
  
  # solve the problem  
  result <-gurobi(model, params)
  return(2*result$x-1)
}




quantile_mix_chi <- function(alpha, lambdas, nrep = 10000) {
  # function to calculate empirical quantile of sum lambda chisq statistics
  p <- length(lambdas)
  chi_mat <- matrix(rchisq(nrep*p, 1), nrep, p)
  mix_alpha<- as.numeric(quantile(chi_mat %*% lambdas, prob=alpha))
  return(mix_alpha)
}

opt.design.linear <- function(W, Z, covariates = TRUE) {
  ##########################################################################
  # Function for optimal design with the linear model in Parker et al, 2017
  # Inputs:
  # W: adjacency matrix
  # Z: covariate matrix, intercept not included
  # rho: network correlation parameter
  #
  # Outputs:
  # optimal allocation (-1, 1)
  ##########################################################################
  n <- nrow(W)
  nedges <- apply(W, 1, sum)
  nrep <- 1000
  if(covariates==FALSE) {
    linear.obj <- function(x) {
      Fmat <- cbind(1, nedges, W%*%x)
      Fx <- t(Fmat) %*% x
      obj<-try(as.vector(t(Fx)%*%solve(t(Fmat)%*%Fmat, Fx)), TRUE)
      if(!is.numeric(obj)) obj <- NA
      return(obj)
    }
  }
  if(covariates==TRUE) {
    linear.obj <- function(x) {
      Fmat <- cbind(1, nedges, W%*%x, Z)
      Fx <- t(Fmat) %*% x
      obj<-try(as.vector(t(Fx)%*%solve(t(Fmat)%*%Fmat, Fx)), TRUE)
      if(!is.numeric(obj)) obj <- NA
      return(obj)
    }
  }
  designs <- replicate(nrep, sample(c(rep(1, n/2), rep(-1, n/2))))
  objs <- apply(designs, 2, linear.obj)
  x <- designs[,which.min(objs)]
  return(list(x=x, Erand=mean(objs)))
}



