rm(list = ls())
source('/home/euan/documents/BGSE_masters_thesis/src/R/dcc_C.R')
library(ccgarch)
set.seed(42)

nobs <- 1000; cut <- 1000; nu <- 8
N <- 50
a <- runif(N,0.01,0.05)
A <- diag(runif(N,0.1,0.2))
B <- diag(runif(N,0.6,0.8))
x <- (N:1)/N
uncR <- toeplitz(x)
dcc.para <- c(0.05,0.9)

# Simulating data from the original DCC-GARCH(1,1) process, storing bias
bias2 <- c(0,0)
for(i in 1:25){
  print(i)
  X <- dcc.sim(nobs, a, A, B, uncR, dcc.para, model="diagonal")
  est <- DCC.fit.eps(X$eps, 5, 0.5)
  bias2 <- bias2 - dcc.para + est$par
}
bias2 = bias2/25
print(bias2)

n = c(5,10,20,30,40,50)
biases = matrix(rep(0,12),2,6)
for(j in 1:length(n)){
  for(i in 1:25){
    print(i)
    X <- dcc.sim(nobs, a, A, B, uncR, dcc.para, model="diagonal")
    est <- DCC.fit(X$eps, n[j], 0.5)
    bias2 <- bias2 - dcc.para + est$par
  }
  bias2 = bias2/25
  biases[,j] = bias2
}
