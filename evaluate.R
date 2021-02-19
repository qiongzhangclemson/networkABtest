evaluate <- function(x, W, Z, rho, linear=NA) {
	        # function to compute the D-optimal efficiency given an allocation x
          var_re <- 0
          if(length(x)>0) {
            #car model evaluate
  			    n <- nrow(W)
	  		    nedge <- rowSums(W)
		  	    nedge.total <- sum(W)
			      p <- ncol(Z)
            Sigma.net <- Diagonal(n, x =nedge)-rho*W
            Fmat <- cbind(1, Z)
			      Fsig <- t(Fmat) %*% Sigma.net
			      xFsig <- Fsig %*% x
            Sigma <- Fsig %*% Fmat
            Sigma2 <- sqrtm(Fsig %*% t(Fsig))		
            T2mat <- Sigma2 %*% solve(Sigma, Sigma2)
            lambdas <- eigen(T2mat, only.values = TRUE)$values
            T1 <- rho * as.numeric(t(x)%*% W %*% x)
            T2 <- as.numeric(t(xFsig) %*% solve(Sigma, xFsig))
            obj <- nedge.total-T1-T2
            eff <- (obj/((1+rho)*sum(W)))
            var_re <- 1-(nedge.total-sum(lambdas))/obj
            Erand <- nedge.total-sum(lambdas)
            
            #linear model 
            linear.obj <- function(x) {
              Fmat <- cbind(1, nedge, W%*%x, Z)
              Fx <- t(Fmat) %*% x
              obj<-try(as.vector(t(Fx)%*%solve(t(Fmat)%*%Fmat, Fx)), TRUE)
              if(!is.numeric(obj)) obj <- NA
              return(obj)
            }
            var_re.linear<- (linear - linear.obj(x))/linear
          }
          return(c(var_re,obj, Erand, T1, T2))
}



evaluate_sar <- function(x, W, Z, rho) {
  # function to compute the D-optimal efficiency given an allocation x
  var_re <- 0
  if(length(x)>0) {
    n <- nrow(W)
    nedge <- rowSums(W)
    nedge.total <- sum(W)
    p <- ncol(Z)
    Sigma.net <- Diagonal(n, x =1)-rho*W
    Rx <- Sigma.net %*% x
    Sigma.net <- Sigma.net %*% Sigma.net
    Fmat <- cbind(1, Z)
    Fsig <- t(Fmat) %*% Sigma.net
    xFsig <- Fsig %*% x
    Sigma <- Fsig %*% Fmat
    obj <- as.numeric(sum(Rx^2)-t(xFsig) %*% solve(Sigma, xFsig))
    
    Sigma2 <- sqrtm(Fsig %*% t(Fsig))		
    T2mat <- Sigma2 %*% solve(Sigma, Sigma2)
    lambdas <- eigen(T2mat, only.values = TRUE)$values
    obj_rand <- n - sum(lambdas)+rho^2*sum(diag(W%*%W))
  }
  return(c(obj, obj_rand))
}
