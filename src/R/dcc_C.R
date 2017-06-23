dyn.load('/home/euan/documents/BGSE_masters_thesis/R-Package-dynamo/src/dynamo.so')

source('/home/euan/documents/BGSE_masters_thesis/src/R/data_prep.R')
#dcc_filter(int *status, double **omega, double **eps,
# double *loglik, double *param, int *T, int *N)

dcc.filter <- function( eps , params , L){
  
  T   <- nrow(eps)
  N   <- ncol(eps)
  
  if( any(!is.finite(params)) ){ 
    filter = list( loglik=-Inf , sigma2=rep(NA,T) )    
    return( filter ) 
  }  
  
  result <- .C( 'dcc_filter', 
                status = as.integer(0), 
                L = as.double(L) , 
                eps = as.double(eps) ,
                loglik = as.double(0) , 
                params = as.double(params) , 
                as.integer(T),
                as.integer(N),
                PACKAGE="dynamo"
  )
  
  return(result)
}

params <- c(0.01,0.6)
alpha = params[1]; beta = params[2]
T <- ncol(eps)
omega  <- (1/T) * (1 - alpha - beta) * (eps %*% t(eps))
L   <- t(chol(omega))

L[1,1] + L[2,2]

params = c(0.05,0.9)
eps = eps[c(1,2),]
dcc.filter(t(eps), params, L)

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
  r <- sqrt(L[n,n] + work[n]*work[n])
  c <- r / L[n,n]
  s <- work[i] / L[n,n]
  L[n,n] <- r
  return(sqrt(b) * L)
}

testing_for_C <- function(params, eps, L){
  # Cholesky decomposition of omega, and 
  # Calculate cholesky decomposition of R_t 
  alpha = params[1]
  beta  = params[2]
  L_tilde <- diag(1/sqrt(rowSums(L**2))) %*% L
  l   <- sum(log(diag(L_tilde)))
  t = 1
  for(t in 2:T){
    L   <- rank_one_update(L,eps[,t], a = alpha, b = beta)
    L_tilde <- diag(1/sqrt(rowSums(L**2))) %*% L
    
    l   <- l + sum(log(diag(L_tilde)))
  }
  return(l)
}

params = c(0.05,0.9)
# Calculate omega hat matrix and use to initialise log-likelihood calculation
T           <- ncol(eps)
omega       <- (1/T) * (1 - alpha - beta) * (eps %*% t(eps))
L           <- t(chol(omega))
testing_for_C(params, eps, L)
dcc.filter(t(eps), params, L)$loglik

