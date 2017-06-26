dyn.load('/home/euan/documents/BGSE_masters_thesis/R-Package-dynamo/src/dynamo.so')
source('/home/euan/documents/BGSE_masters_thesis/src/R/data_prep.R')

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

dcc.filter <- function( shrinkage_coefs, shrinkage_vectors, eps , params , L){
  
  n   <- nrow(shrinkage_vectors)
  T   <- nrow(eps)
  N   <- ncol(eps)
  
  if( any(!is.finite(params)) ){ 
    filter = list( loglik=Inf , sigma2=rep(NA,T) )    
    return( filter ) 
  }  
  
  result <- .C( 'dcc_filter', 
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

DCC.fit <- function(X, n, rho, Trace = 3){
  # Given matrix of returns, fits DCC model
  # Step 1: degarching the data
  N <- ncol(X)
  T <- nrow(X)
  eps <- matrix(rep(0, N*T), N,T)
  for(i in 1:N){
    eps[i,] <- X[,i]/sqrt(garchFit(data = X[,i], trace = F)@h.t)
  }
  dim(eps)
  # Step 2: Optimising the DCC model parameters
  omega  <- (1/T)*(eps %*% t(eps))
  L   <- t(chol(omega))
  
  shrink_est  <- linear_shrinkage(eps, n, rho)
  
  coefs <- shrink_est$coefficients; vectors = shrink_est$vectors
  
  llh <- function(params){
    return(as.matrix(dcc.filter(coefs, t(vectors), t(eps), params, L)$loglik))
  }
  
  params <- c(0.05,0.9)
  lowerBounds <- c(10E-6, 10E-6)
  upperBounds <- 1 - lowerBounds
  if(Trace > 0){fit = nlminb(start = params, objective = llh,
                             lower = lowerBounds, upper = upperBounds, control =  list(trace=Trace))}
  else{  fit = nlminb(start = params, objective = llh,
                      lower = lowerBounds, upper = upperBounds)}
  return(list(par = fit$par, loglik = -0.5*fit$objective))
}

DCC.fit.eps <- function(eps, n, rho, Trace = 3){
  # Given matrix of epsilons, fits DCC model
  
  N <- nrow(eps)
  T <- ncol(eps)
  
  # Step 2: Optimising the DCC model parameters
  omega  <- (1/T)*(eps %*% t(eps))
  L   <- t(chol(omega))
  
  shrink_est  <- linear_shrinkage(eps, n, rho)
  
  coefs <- shrink_est$coefficients; vectors = shrink_est$vectors
  
  llh <- function(params){
    return(as.matrix(dcc.filter(coefs, t(vectors), t(eps), params, L)$loglik))
  }
  
  params <- c(0.05,0.9)
  lowerBounds <- c(10E-6, 10E-6)
  upperBounds <- 1 - lowerBounds
  if(Trace > 0){fit = nlminb(start = params, objective = llh,
                             lower = lowerBounds, upper = upperBounds, control =  list(trace=Trace))}
  else{  fit = nlminb(start = params, objective = llh,
                      lower = lowerBounds, upper = upperBounds)}
  return(list(par = fit$par, loglik = -0.5*fit$objective))
}