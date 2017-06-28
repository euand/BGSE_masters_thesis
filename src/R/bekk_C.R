rm(list = ls())
source('/home/euan/documents/BGSE_masters_thesis/src/R/dcc_C.R')

linear_shrinkage <- function(eps, n){
  # Linear shrinkage estimator of covariance function
  # eps: matrix of de-garched returns
  # n: number of eigenvalues to take
  T             <- ncol(eps)
  N             <- nrow(eps)
  
  sample_covariance_matrix  <- (1/T)*(eps %*% t(eps))
  eigenvals                 <- eigen(sample_covariance_matrix)
  
  m    <- tr(sample_covariance_matrix %*% diag(rep(1,N)))/N
  d_sq <- norm(sample_covariance_matrix - diag(rep(m,N)), type = 'f')
  b_bar <- norm(eps[,1] %*% t(eps[,1]) - sample_covariance_matrix, type = 'f')^2
  for(t in 2:T){
    b_bar <- b_bar + norm(eps[,t] %*% t(eps[,t]) - sample_covariance_matrix, type = 'f')^2
  }
  b_bar <- b_bar/(T**2)
  b <- min(b_bar, d_sq)
  rho <- b/d_sq
  
  eigenvectors  <- eigenvals$vectors
  lambda        <- eigenvals$values
  
  lambda_bar    <- mean(lambda)
  C_bar         <- matrix(rep(0,(N**2)), N,N)
  shrink_coef   <- rep(0,n)
  for (i in 1:N){
    shrink_coef[i]  <- rho * lambda_bar + (1 - rho)*lambda[i]
    C_bar           <- C_bar + shrink_coef[i] * (eigenvectors[,i] %*% t(eigenvectors[,i]))
  }
  return(list(coefficients = shrink_coef[1:n], vectors = eigenvectors[,1:n], est = C_bar))
}

BEKK.filter <- function( shrinkage_coefs, shrinkage_vectors, eps , params , L){
  
  n   <- nrow(shrinkage_vectors)
  T   <- nrow(eps)
  N   <- ncol(eps)
  
  if( any(!is.finite(params)) ){ 
    filter = list( loglik=Inf , sigma2=rep(NA,T) )    
    return( filter ) 
  }  
  
  result <- .C( 'bekk_filter', 
                status = as.integer(0), 
                L = as.double(L) , 
                eps = as.double(eps) ,
                loglik = as.double(0) , 
                params = as.double(params) , 
                as.double(shrinkage_coefs),
                as.double(shrinkage_vectors),
                as.integer(T),
                as.integer(N),
                as.integer(n),
                PACKAGE="dynamo"
  )
  
  return(list('loglik' = result$loglik))
}

BEKK.fit <- function(X, n, Trace = 3){
  # Given matrix of returns, fits DCC model
  # Step 1: degarching the data
  N <- ncol(X)
  T <- nrow(X)
  eps <- matrix(rep(0, N*T), N,T)
  for(i in 1:N){
    eps[i,] <- X[,i]/sqrt(garchFit(data = X[,i], trace = F)@h.t)
  }
  # Step 2: Optimising the DCC model parameters
  shrink_est  <- linear_shrinkage(t(X), n)
  omega <- shrink_est$est
  omega <- diag(diag(omega**(-1/2))) %*% omega %*% diag(diag(omega**(-1/2)))
  L   <- t(chol(omega))
  
  coefs <- shrink_est$coefficients; vectors = shrink_est$vectors
  llh <- function(params){
    return(as.matrix(BEKK.filter(coefs, t(vectors), t(eps), params, t(L))$loglik))
  }
  params <- c(0.07,0.9)
  
  lowerBounds <- c(10E-6, 10E-6)
  upperBounds <- 1 - lowerBounds
  if(Trace > 0){fit = nlminb(start = params, objective = llh,
                             lower = lowerBounds, upper = upperBounds, control =  list(trace=Trace))}
  else{  fit = nlminb(start = params, objective = llh,
                      lower = lowerBounds, upper = upperBounds)}
  return(list(par = fit$par, loglik = -0.5*fit$objective))
}
