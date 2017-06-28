rm(list = ls())
source('/home/euan/documents/BGSE_masters_thesis/src/R/dcc_C.R')
library(ccgarch)
set.seed(42)

nobs <- 1000; cut <- 1000; nu <- 8
N <- 10
a <- rep(0.01, N)
A <- diag(rep(0.05, N))
B <- diag(rep(0.9,N))
x <- (N:1)/N
uncR <- toeplitz(x)
d.f <- Inf
dcc.para <- c(0.05,0.93)

# Simulating data from the original DCC-GARCH(1,1) process, storing bias
n <- c(1,2,3,4,5,6,7,8,9,10)
biases = matrix(rep(0,length(n)*2),2,length(n))
for(j in 1:length(n)){
  bias2 <- 0
  for(i in 1:25){
    print(i)
    X <- dcc.sim(nobs, a, A, B, uncR, dcc.para, model="diagonal")
    est <- DCC.fit(X$eps, n[j], Trace = 3)
    bias2 <- bias2 - dcc.para + est$par
  }
  bias2 <- bias2/25
  biases[,j] <- bias2
}
bias2 = bias2/25
print(bias2*25)
print(biases)

#plot(abs(biases[1,]), type = 'l', ylim = c(0,0.1))

X <- X$eps
