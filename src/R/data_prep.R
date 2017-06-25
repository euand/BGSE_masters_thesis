#############################################################
# Loading and processing data for use in dcc model estimation
#############################################################
library(fGarch)
library(KFAS)

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

X <- cbind(x,y,z)
dim(X)

# Step 1: degarching the data
x_garch = garchFit(data = x, trace = F); eps_x = x/sqrt(x_garch@h.t)
y_garch = garchFit(data = y, trace = F); eps_y = y/sqrt(y_garch@h.t)
z_garch = garchFit(data = z, trace = F); eps_z = z/sqrt(z_garch@h.t)
eps = rbind(eps_x, eps_y, eps_z)
