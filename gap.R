
require(Matrix)

rm(list=ls())
T_value <- function(x, W, Z, rho) {
        # function to compute the T value given an allocation x
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
        T1 <- rho * as.numeric(t(x)%*% W %*% x)
        T2 <- as.numeric(t(xFsig) %*% solve(Sigma, xFsig))
        T_value <- nedge.total-T1-T2
        return(T_value)
}

#random graph of n=50 and possibility of edge=0.25

set.seed(3)
n <- 50
W <- matrix(0, nrow=n, ncol=n)
for (i in 1:(n-1)) {
    W[i,(i+1):n] <- rbinom(n-i, size=1, prob=0.25)        
}
W <- W+t(W)
D <- diag(rowSums(W))
m <- sum(rowSums(W))
Z <- rnorm(n, mean=0, sd=10)

rho_0 <-0.5

S <- D-rho_0*W
var_rho <- 1/12
Eigen_S <- eigen(S)
lamda_max <- Eigen_S$values[1]
lamda_min <- Eigen_S$values[n]
Eigen_W <- eigen(W)
lamda_spec <- max(abs(Eigen_W$values))

upper_bound <- min(c(n*lamda_max,(1+rho_0)*m))*(lamda_spec^2)*var_rho/(lamda_min^2)
upper_bound2 <- (m+2*sqrt(m))*(lamda_spec^2)*var_rho/(lamda_min^2)

K <- 400
diff <- numeric(K)
T_0s <- numeric(K)
for (k in 1:K) {
x <- rbinom(n, size=1, prob=0.5)
x[x==0] <- -1
T_rho0 <- T_value(x,W,Z,rho=rho_0)
T_0s[k] <- T_rho0

B <- 200
rhos <- runif(B, min=0, max=1)
Ts <- numeric(B)
for (i in 1:B) {
        Ts[i] <- T_value(x,W,Z,rho=rhos[i])
}
mean_T <- mean(Ts)

diff[k] <- T_rho0-mean_T
}

par(mfrow=c(1,2))
hist(T_0s, main=expression(paste("Histogram of T(", bold(x), ',', rho[0],")")), xlab=expression(paste("T(",bold(x),",",rho[0],")")))
hist(diff, main=expression(paste("Histogram of T(", bold(x), ',', rho[0],")-E(T(", bold(x), ',', rho,"))")), xlab=expression(paste("T(", bold(x), ',', rho[0],")-E(T(", bold(x), ',', rho,"))")))

par(mfrow=c(2,2))
x <- rbinom(n, size=1, prob=0.5)
x[x==0] <- -1
T_rho0 <- T_value(x,W,Z,rho=rho_0)

B <- 200
rhos <- runif(B, min=0, max=1)
Ts <- numeric(B)
for (i in 1:B) {
        Ts[i] <- T_value(x,W,Z,rho=rhos[i])
}
plot(rhos,Ts,type='l',xlab=expression(rho),ylab=expression(paste("T(",bold(x),",",rho,")")))
abline(h=T_rho0,col='red')
abline(h=mean(Ts),col='green')
legend("topleft",legend=c(expression(paste("T(",bold(x),",",rho[0],")")), "sample mean"),col=c('red','green'),lty='solid')