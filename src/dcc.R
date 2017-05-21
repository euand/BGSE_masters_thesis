rm(list=ls())
library(fGarch)
library(KFAS)

data1 <- read.csv('/home/euan/documents/financial-econ/data/tickers/DBB.csv')
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

# Step 1: degarching the data

x_garch = garchFit(data = x)
eps_x = x/sqrt(x_garch@h.t)

y_garch = garchFit(data = y)
eps_y = y/sqrt(y_garch@h.t)

z_garch = garchFit(data = z)
eps_z = z/sqrt(z_garch@h.t)

eps = rbind(eps_y, eps_x, eps_z)

r_bar = function(eps){
  R = matrix(rep(0,nrow(eps)), nrow(eps), nrow(eps))
  for(i in 1:nrow(eps)){
    R = R + eps[,i] %*% t(eps[,i])
  }
  R = R / nrow(eps)
  return(R)
}

omega_gen = function(eps, alpha, beta){
  R = r_bar(eps)
  omega_hat = (1 - alpha - beta) * R
  return(omega_hat)
}

get_W <- function(omega,k){
  s = svd(omega)
  if(k == 1){
    A = s$d[1] * t(t(s$u[,1])) %*% t(s$v)[1,]
  }
  else{
    A = s$u[,1:k] %*% diag(s$d[1:k]) %*% t(s$v)[1:k,]
  }
  W = chol(A)
}

rank_one_update <- function(L, x, a, b){
  n = length(x)
  work = x
  L_new = L
  for(i in 1:n){
    r = sqrt(a*L[i,i])
    c = r/(sqrt(a)*L[i,i])
    s = (work[i])/(sqrt(a)*L[i,i])
    L_new[i,i] = r
    if(i == n){return(L_new)}
    for(j in (i+1):n){
      L_new[j,i] = (sqrt(a) * L[j,i] + s*work[j]) / c
      work[j] = c*work[j] - s*L_new[j,i]
    }
  }
}

rank_r_update <- function(L, W){
  # L is upper triangular
  # Calculate LDL decomposition
  dd <- diag(L) 
  ch <- t(L/dd) 
  DD <- dd^2 

  work = W
  r = ncol(W)
  n = nrow(L)
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

rank_r_update(L,W)
t(chol(L%*%t(L) + W%*%t(W)))

log_likelihood_full_updates <- function(params){
  alpha <- params[1]
  beta <- params[2]
  
  if((alpha+beta) >= (1-10E-6)){return(10E6)}
  
  T <- ncol(eps)
  omega = omega_gen(eps, alpha, beta)
  
  loglik = rep(0,T)
  
  L <- chol(omega)
  l = sum(log(diag(L)))
  y = solve(diag(diag(omega**(-1/2))) %*% omega %*% diag(diag(omega**(-1/2))), eps[,1])
  l2 = y %*% diag(diag(omega**(-1/2))) %*% omega %*% diag(diag(omega**(-1/2))) %*% t(t(y))
  l = l + l2
  Q_t_minus = omega
  for(t in 2:T){
    #Qt = omega + alpha * t(t(eps[,t])) %*% eps[,t] + beta * Q_t_minus
    #L = t(chol(Qt))
    L = rank_one_update(L,eps[,t])
    l = l + sum(log(diag(L)))
    y = solve(diag(diag(Qt**(-1/2))) %*% Qt %*% diag(diag(Qt**(-1/2))), eps[,t])
    l = l + y %*% diag(diag(omega**(-1/2))) %*% omega %*% diag(diag(omega**(-1/2))) %*% t(t(y))
    Q_t_minus = Qt
  }
  return(l)
}
lower = c(10E-6,10E-6)
upper = c(1-10E-6, 1-10E-6)

params = c(0.3,0.3)

fit = nlminb(start = params, objective = log_likelihood_full_updates,
             lower = lower, upper = upper,  control = list(trace=3))
