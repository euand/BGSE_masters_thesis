rm(list = ls())
source('/home/euan/documents/BGSE_masters_thesis/src/R/dcc_C.R')

time0 <- Sys.time()
n <- ceiling(nrow(eps)/2)
DCC.fit(X, 2, 0.5)
time1 <- Sys.time()
print(time1 - time0)