rm(list=ls())
library(fGarch)

DCC.fit <- function(X){
  
  N <- nrow(X)
  T <- ncol(X)
  
  eps <- matrix(0,N,T)
  for(i in 1:N){
    eps[i,] <- X[i,]/sqrt(garchFit(data = X[i,])@h.t)
  }  
  
  omega_gen = function(eps, alpha, beta){
    R = matrix(rep(0,nrow(eps)), nrow(eps), nrow(eps))
    for(i in 1:nrow(eps)){
      R = R + eps[,i] %*% t(eps[,i])
    }
    R = R / nrow(eps)
    omega_hat = (1 - alpha - beta) * R
    return(omega_hat)
  }
  
  rank_one_update <- function(L, x, a=1, b=1){
    n = length(x)
    work = x
    for(i in 1:n){
      r = sqrt(a*L[i,i])
      c = r/(sqrt(a)*L[i,i])
      s = (work[i])/(sqrt(a)*L[i,i])
      L[i,i] = r
      if(i == n){return(L)}
      for(j in (i+1):n){
        L[j,i] = (sqrt(a) * L[j,i] + s*work[j]) / c
        work[j] = c*work[j] - s*L[j,i]
      }
    }
  }
  
  rank_r_update <- function(L, W){
    # L is upper triangular
    # Calculate LDL decomposition
    dd <- diag(L) 
    ch <- L/dd
    DD <- dd^2 
    
    work = W
    r = ncol(W)
    n = nrow(L)
    alpha = rep(0,r); gamma = rep(0,r)
    for(i in 1:r){
      alpha[i] = 1; gamma[i] = 1
    }
    for(j in 1:n){
      for(i in 1:r){
        alpha_bar = alpha[i] + work[j,i]*work[j,i]/DD[j]
        DD[j] = DD[j] * alpha_bar
        gamma[i] = work[j,i]/DD[j]
        DD[j] = DD[j]/alpha[i]
        alpha[i] = alpha_bar
      }
      for(p in j:n){
        for(i in 1:r){
          work[p,i] = work[p,i] - work[j,i] * ch[p,j];
          ch[p,j] = ch[p,j] + gamma[i] * work[p,i];
        }
      }
    }
    return(ch*sqrt(DD))
  }
  
  frob_dist <- function(U){
    U = matrix(U, N, r)
    dist = norm(omega_hat - U %*% t(U), type = 'f')
    return(dist)
  }

  log_likelihood_full_updates <- function(params){
    alpha <- params[1]
    beta <- params[2]
    
    if((alpha+beta) >= (1-10E-6)){return(10E6)}
    
    T <- ncol(eps)
    omega = omega_gen(eps, alpha, beta)
    
    loglik = rep(0,T)
    
    L <- t(chol(omega))
    l = sum(log(diag(L)))
    y = solve(diag(diag(omega**(-1/2))) %*% omega %*% diag(diag(omega**(-1/2))), eps[,1])
    l2 = y %*% diag(diag(omega**(-1/2))) %*% omega %*% diag(diag(omega**(-1/2))) %*% t(t(y))
    l = l + l2
    Q_t_minus = omega
    for(t in 2:T){
      Qt = omega + alpha * t(t(eps[,t])) %*% eps[,t] + beta * Q_t_minus
      L = t(chol(Qt))
      y = solve(diag(diag(Qt**(-1/2))) %*% Qt %*% diag(diag(Qt**(-1/2))), eps[,t])
      l = l + y %*% diag(diag(omega**(-1/2))) %*% omega %*% diag(diag(omega**(-1/2))) %*% t(t(y))
      Q_t_minus = Qt
    }
    return(-l)
  }
  
  lower = c(10E-6,10E-6)
  upper = c(1-10E-6, 1-10E-6)
  
  params = c(0.3,0.3)
  
  fit = nlminb(start = params, objective = log_likelihood_full_updates,
               lower = lower, upper = upper,  control = list(trace=3))
  
  return(fit$par)
}

DCC.fit(X)
