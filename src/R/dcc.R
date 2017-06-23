source('/home/euan/documents/BGSE_masters_thesis/src/R/data_prep.R')

############################################
# Calculating the log likelihood
###########################################

rank_one_update <- function(L, x, a=1, b=1){
  n = length(x)
  work = sqrt(a/b) * x
  for(i in 1:(n-1)){
    r <- sqrt(L[i,i] + work[i]*work[i])
    c <- r / L[i,i]
    s <- work[i] / L[i,i]
    L[i,i] <- r
    for(j in (i+1):n){
      L[j,i]  <- (L[j,i] + s * work[j]) / c
      work[j] <- c * work[j] - s * L[j,i]
    }
  }
  r <- sqrt(L[n,n])
  c <- r / L[n,n]
  s <- work[i] / L[n,n]
  L[n,n] <- r
  return(sqrt(b) * L)
}

linear_shrinkage <- function(eps, n, rho){
  # Linear shrinkage estimator of covariance function
  # eps: matrix of de-garched returns
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
    L     <- rank_one_update(L, vectors[,i], a = (1 - alpha - beta) * coefs[i], b = 1) 
  }
  return(L)
}

loglik_part2 <- function(L_tilde,eps){
  y <- solve(L_tilde %*% t(L_tilde), eps)
  return( eps %*% y )
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
  # Calculate cholesky decomposition of R_t 
  L_tilde <- diag(1/sqrt(rowSums(L**2))) %*% L
  l   <- sum(log(diag(L_tilde)))
  l2  <- loglik_part2(L_tilde, eps[,1]) 
  l   <- l + l2
  for(t in 2:T){
    L   <- rank_one_update(L,eps[,t], a = alpha, b = beta)
    L   <- shrink_est_update(L, shrink_est, alpha, beta)
    L_tilde <- diag(1/sqrt(rowSums(L**2))) %*% L
    y <- solve(L_tilde %*% t(L_tilde), eps[,t])
    l   <- l + sum(log(diag(L_tilde))) + eps[,t] %*% y
    }
  return(l)
}

###################
# Fitting the model
###################

params <- c(0.01,0.6)
lowerBounds <- c(10E-6, 10E-6)
upperBounds <- 1 - lowerBounds
fit = nlminb(start = params, objective = log_likelihood_DCC,
             lower = lowerBounds, upper = upperBounds ,  control = list(trace=3))

time0 = Sys.time()
log_likelihood_DCC(params)
time1 = Sys.time()
print(time1 - time0)



