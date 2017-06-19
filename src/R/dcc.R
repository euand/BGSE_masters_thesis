library(fGarch)
library(KFAS)

##############
# Read in data
##############

data1 <- read.csv('/home/euan/documents/financial-econ/data/tickers/AAPL.csv')
data1 <- data1[seq(nrow(data1),1,-1),]
data1$Date <- as.Date(data1$Date)
price1 <- data1$Adj.Close
x <- ( log(price1[2:length(price1)]) - log(price1[1:(length(price1)-1)]) ) * 100

data2 <- read.csv('/home/euan/documents/financial-econ/data/tickers/BAC.csv')
data2 <- data2[seq(nrow(data2),1,-1),]
data2$Date <- as.Date(data2$Date)
price2 <- data2$Adj.Close
y <- ( log(price2[2:length(price2)]) - log(price2[1:(length(price2)-1)]) ) * 100

data3 <- read.csv('/home/euan/documents/financial-econ/data/tickers/CVX.csv')
data3 <- data3[seq(nrow(data3),1,-1),]
data3$Date <- as.Date(data3$Date)
price3 <- data3$Adj.Close
z <- ( log(price3[2:length(price3)]) - log(price3[1:(length(price3)-1)]) ) * 100

z <- z[data3$Date >= data1$Date[1]]
z <- z[!is.na(z)]
y <- y[data2$Date >= data1$Date[1]]
y <- y[!is.na(y)]
x <- x[1:length(y)]
z <- z[1:length(y)]

X <- rbind(x,y,z)

# Step 1: degarching the data
x_garch = garchFit(data = x); eps_x = x/sqrt(x_garch@h.t)
y_garch = garchFit(data = y); eps_y = y/sqrt(y_garch@h.t)
z_garch = garchFit(data = z); eps_z = z/sqrt(z_garch@h.t)
eps = rbind(eps_x, eps_y, eps_z)

############################################
# Calculating the log likelihood
###########################################

omega_gen = function(eps, alpha, beta){
  N <- nrow(eps)
  T <- ncol(eps)
  R <- (1/T)*(eps %*% t(eps))
  omega_hat <- (1 - alpha - beta) * R
  return(omega_hat)
}

rank_one_update <- function(L, x, a=1, b=1){
  n = length(x)
  work = sqrt(b/a) * x
  for(i in 1:(n-1)){
    r <- sqrt(L[i,i])
    c <- r / L[i,i]
    s <- work[i] / L[i,i]
    L[i,i] <- r
    for(j in (i+1):(n)){
      L[j,i]  <- (sqrt(a) * L[j,i] + s * work[j]) / c
      work[j] <- c * work[j] - s * L[j,i]
    }
  }
  r <- sqrt(L[n,n])
  c <- r / L[n,n]
  s <- work[i] / L[n,n]
  L[n,n] <- r
  return(sqrt(a) * L)
}

loglik_part2 <- function(L,eps){
  N <- nrow(L)
  A <- diag(1/sqrt(rowSums(L**2)))
  Q <- A %*% L
  y <- solve(Q %*% t(Q), eps)
  return( t(y) %*% (Q %*% t(Q)) %*% y )
}

linear_shrinkage <- function(eps, n, rho){
  # Linear shrinkage estimator of covariance function
  # eps: matrix of degarched returns
  # n: number of eigenvalues to take
  T             <- ncol(eps)
  N             <- nrow(eps)
  
  sample_covariance_matrix  <- (1/T)*(eps %*% t(eps))
  eigenvals                 <- eigen(sample_covariance_matrix)
  
  eigenvectors  <- eigenvals$vectors
  lambda        <- eigenvals$values
  
  lambda_bar    <- mean(lambda)
  C_bar         <- matrix(rep(0,(N**2)), N,N)
  shrink_coef   <- rep(0,n)
  for (i in 1:n){
    shrink_coef[i]  <- rho * lambda_bar + (1 - rho)*lambda[i]
    C_bar           <- C_bar + shrink_coef[i] * (eigenvectors[,i] %*% t(eigenvectors[,i]))
  }
  return(list(coefficients = shrink_coef, vectors = eigenvectors[,1:n]))
}

shrink_est_update <- function(L, shrink_est, alpha, beta){
  coefs   <- shrink_est$coefficients
  vectors <- shrink_est$vectors
  for( i in 1:length(coefs)){
    L     <- rank_one_update(L, vectors[,i], a = 1, b = (1 - alpha - beta) * coefs[i]) 
  }
  return(L)
}

log_likelihood_DCC <- function(params){
  alpha <- params[1]
  beta <- params[2]
  
  # Check parameters are in bounds, else return large value
  if((alpha+beta) >= (1-10E-6)){return(10E6)}
  
  # Calculate omega hat matrix and use to initialise log-likelihood calculation
  T           <- ncol(eps)
  omega       <- (1/T) * (1 - alpha - beta) * (eps %*% t(eps))
  shrink_est  <- linear_shrinkage(eps, 2, 0.5)
  
  # Cholesky decomposition of omega, and 
  L   <- t(chol(omega))
  l   <- sum(log(diag(L)))
  l2  <- loglik_part2(L, eps[,1]) 
  loglik2 <- rep(0,T)
  l   <- l + l2
  for(t in 2:T){
    L   <- rank_one_update(L,eps[,t], a = alpha, b = beta)
    L   <- shrink_est_update(L, shrink_est, alpha, beta)
    l   <- l + sum(log(diag(L)))
    l2  <- loglik_part2(L, eps[,t]) 
    loglik2[t] <- l2
    l   <- l + l2
  }
  return(l)
}

# Step 3: Estimate Parameters and Compute Numerically Hessian:
params <- c(0.6,0.2)
lowerBounds <- c(10E-6, 10E-6)
upperBounds <- 1 - lowerBounds
fit = nlminb(start = params, objective = log_likelihood_DCC,
             lower = lowerBounds, upper = upperBounds ,  control = list(trace=3))

log_likelihood_DCC(params)
