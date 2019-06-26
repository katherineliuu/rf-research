# Simulation studies to investigate RF and LM prediction intervals
library(randomForestSRC)
###################################################################
# Simulation 1: Linear Model (one variable)

set.seed(06232019)#set seed
iter<-2
MSE.lr.1 <- c(rep(NA,iter))
MSE.rf.1 <- c(rep(NA,iter))
meanWidth.lr.1 <- c(rep(NA,iter))
meanWidth.rf.1 <- c(rep(NA,iter))
percntPI.lr.1 <- c(rep(NA,iter))
percntPI.rf.1 <- c(rep(NA,iter))

for (i in 1:iter){
n <- 2000  #number of observations (trainig and test set combined)
sig <- 5  #error standard deviation
x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
e <- rnorm(n, 0, sig) #generate error term 
y <- 2*x+3 + e  #generate response values 
data <- data.frame(x, y) #combine explanatory and response varables into data frame
train <- data[1:(n/2), ] #use first half of observations as training data
test <- data[(n/2+1):n, ] #use second half of observations as test data
#-----------------------------------------------------------------
#LR
M <- lm(data=train, y~x)
yhat <- predict(M, newdata=test) # predictions
PI<- predict(M, newdata=test, interval="prediction", level=0.95) # prediction intervals with 95%
#head(PI)

MSElr1 <- mean((yhat - test$y)^2)
MSE.lr.1[i] <- MSElr1
################################# look into width
meanWidthlr1 <- mean(PI[,3]-PI[,2])
meanWidth.lr.1[i] <- meanWidthlr1
################################# percentage in PI
Rlr1 <- test$y>=PI[,2] & test$y<=PI[,3] 
#sum(Rlr1)
percntPI.lr.1[i] <- sum(Rlr1)/nrow(test)

#RF
o <- quantreg(y ~ ., train, splitrule = "mse", nodesize = 5)  ## quantile regression with mse splitting
o.test <- quantreg(object=o, newdata=test)   # test call
#head(o.test$predicted)  #predictions

quant.dat <- get.quantile(o.test, c(.025, .975))
#head(quant.dat)    #prediction intervals
MSErf1 <- o.test$err.rate[n/2]
MSE.rf.1[i] <- MSErf1
################################# look into width
meanWidthrf1 <- mean(quant.dat[,2]-quant.dat[,1])
meanWidth.rf.1[i] <- meanWidthrf1
################################# percentage in PI
Rrf1 <- test$y>=quant.dat[,1] & test$y<=quant.dat[,2] 
#sum(Rrf1)
percntPI.rf.1[i] <- sum(Rrf1)/nrow(test)
} #for loop

MSE.lr.1
MSE.rf.1 
meanWidth.lr.1 
meanWidth.rf.1 
percntPI.lr.1 
percntPI.rf.1 


save.image(file="filename.Rdata")


###################################################################
# Simulation 2: Nonlinear Model (one variable)

set.seed(06232019) #set seed
MSE.lr.2 <- c(rep(NA,iter))
MSE.rf.2 <- c(rep(NA,iter))
meanWidth.lr.2 <- c(rep(NA,iter))
meanWidth.rf.2 <- c(rep(NA,iter))
percntPI.lr.2 <- c(rep(NA,iter))
percntPI.rf.2 <- c(rep(NA,iter))

for (i in 1:iter){
n <- 2000  #number of observations (trainig and test set combined)
sig <- 5  #error standard deviation
x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
e <- rnorm(n, 0, sig) #generate error term 
y <- 0.1*(x-7)^2-3*cos(x) + 5*log(abs(x)) + 3 + e  #generate response values 
data <- data.frame(x, y) #combine explanatory and response varables into data frame
train <- data[1:(n/2), ] #use first half of observations as training data
test <- data[(n/2+1):n, ] #use second half of observations as test data
-----------------------------------------------------------------
#LR

M <- lm(data=train, y~x)
yhat <- predict(M, newdata=test) # predictions
PI<- predict(M, newdata=test, interval="prediction", level=0.95) # prediction intervals with 95%
head(PI)

MSElr2 <- mean((yhat - test$y)^2)
MSE.lr.2[i] <- MSElr2
################################# look into width
meanWidthlr2 <- mean(PI[,3]-PI[,2])
meanWidth.lr.2[i] <- meanWidthlr2
################################# percentage in PI
Rlr2 <- yhat>=PI[,2] & yhat<=PI[,3] 
percntPI.lr.2[i] <- sum(Rlr2)/nrow(test)

#RF
o <- quantreg(y ~ ., train, splitrule = "mse", nodesize = 5)  ## quantile regression with mse splitting
o.test <- quantreg(object=o, newdata=test)   # test call
#head(o.test$predicted)  #predictions

quant.dat <- get.quantile(o.test, c(.025, .975))
#head(quant.dat)   #prediction intervals
MSErf2 <- o.test$err.rate[n/2]
MSE.rf.2[i] <- MSErf2
################################# look into width
meanWidthrf2 <- mean(quant.dat[,2]-quant.dat[,1])
meanWidth.rf.2[i] <- meanWidthrf2
################################# percentage in PI
Rrf2 <- o.test$predicted>=quant.dat[,1] & o.test$predicted<=quant.dat[,2] 
percntPI.rf.2[i] <- sum(Rrf2)/nrow(test)
} #for loop
MSE.lr.2
MSE.rf.2
meanWidth.lr.2
meanWidth.rf.2 
percntPI.lr.2 
percntPI.rf.2 

###################################################################
# Simulation 3: Linear Model (multivariate)

set.seed(06232019) #set seed
MSE.lr.3 <- c(rep(NA,iter))
MSE.rf.3 <- c(rep(NA,iter))
meanWidth.lr.3 <- c(rep(NA,iter))
meanWidth.rf.3 <- c(rep(NA,iter))
percntPI.lr.3 <- c(rep(NA,iter))
percntPI.rf.3 <- c(rep(NA,iter))

for (i in 1:iter){
n <- 2000  #number of observations (trainig and test set combined)
xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
sig <- 5  #error standard deviation
e <- rnorm(n, 0, sig) #generate error term 
y <- 2*X[,1]+3*X[,4]+4*X[,6]-3*X[,7]+ X[,9] + e  #generate response values 
data <- data.frame(X, y) #combine explanatory and response varables into data frame
train <- data[1:(n/2), ] #use first half of observations as training data
test <- data[(n/2+1):n, ] #use second half of observations as test data
-----------------------------------------------------------------
  #LR
M <- lm(data=train, y~.)
yhat <- predict(M, newdata=test) # predictions
PI<- predict(M, newdata=test, interval="prediction", level=0.95) # prediction intervals with 95%
#head(PI)

MSElr3 <- mean((yhat - test$y)^2)
MSE.lr.3[i] <- MSElr3
###################     look into width
meanWidthlr3 <- mean(PI[,3]-PI[,2])
meanWidth.lr.3[i] <- meanWidthlr3
################################# percentage in PI
Rlr3 <- yhat>=PI[,2] & yhat<=PI[,3] 
percntPI.lr.3[i] <- sum(Rlr3)/nrow(test)


#RF
o <- quantreg(y ~ ., train, splitrule = "mse", nodesize = 5)  ## quantile regression with mse splitting
o.test <- quantreg(object=o, newdata=test)   # test call
#head(o.test$predicted)  #predictions

quant.dat <- get.quantile(o.test, c(.025, .975))
#head(quant.dat)   #prediction intervals
MSErf3 <- o.test$err.rate[n/2]
MSE.rf.3[i] <- MSErf3
################################# look into width
meanWidthrf3 <- mean(quant.dat[,2]-quant.dat[,1])
meanWidth.rf.3[i] <- meanWidthrf3
################################# percentage in PI
Rrf3 <- o.test$predicted>=quant.dat[,1] & o.test$predicted<=quant.dat[,2] 
percntPI.rf.3[i] <- sum(Rrf3)/nrow(test)
} #for loop
MSE.lr.3
MSE.rf.3
meanWidth.lr.3
meanWidth.rf.3
percntPI.lr.3
percntPI.rf.3


###################################################################
# Simulation 4: Nonlinear Model (multivariate)

set.seed(06232019) #set seed
set.seed(06232019) #set seed
MSE.lr.4 <- c(rep(NA,iter))
MSE.rf.4 <- c(rep(NA,iter))
meanWidth.lr.4 <- c(rep(NA,iter))
meanWidth.rf.4 <- c(rep(NA,iter))
percntPI.lr.4 <- c(rep(NA,iter))
percntPI.rf.4 <- c(rep(NA,iter))

for (i in 1:iter){
n <- 2000  #number of observations (trainig and test set combined)
xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
sig <- 5  #error standard deviation
e <- rnorm(n, 0, sig) #generate error term 
y <- (X[,1]-6)^2 + 12*cos(X[,3]) + (X[,7]-5)*(X[,8]-3) + 0.02*(X[,10]-5)^5+ e  #generate response values 
data <- data.frame(X, y) #combine explanatory and response varables into data frame
train <- data[1:(n/2), ] #use first half of observations as training data
test <- data[(n/2+1):n, ] #use second half of observations as test data
-----------------------------------------------------------------
  #LR
M <- lm(data=train, y~.)
yhat <- predict(M, newdata=test) # predictions
PI<- predict(M, newdata=test, interval="prediction", level=0.95) # prediction intervals with 95%
#head(PI)

MSElr4 <- mean((yhat - test$y)^2)
MSE.lr.4[i] <- MSElr4
###################     look into width
meanWidthlr4 <- mean(PI[,3]-PI[,2])
meanWidth.lr.4[i] <- meanWidthlr4
################################# percentage in PI
Rlr4 <- yhat>=PI[,2] & yhat<=PI[,3] 
percntPI.lr.4[i] <- sum(Rlr4)/nrow(test)


#RF
o <- quantreg(y ~ ., train, splitrule = "mse", nodesize = 5)  ## quantile regression with mse splitting
o.test <- quantreg(object=o, newdata=test)   # test call
#head(o.test$predicted)  #predictions

quant.dat <- get.quantile(o.test, c(.025, .975))
#head(quant.dat)   #prediction intervals
MSErf4 <- o.test$err.rate[n/2]
MSE.rf.4[i] <- MSErf4
################################# look into width
meanWidthrf4 <- mean(quant.dat[,2]-quant.dat[,1])
meanWidth.rf.4[i] <- meanWidthrf4
################################# percentage in PI
Rrf4 <- o.test$predicted>=quant.dat[,1] & o.test$predicted<=quant.dat[,2] 
percntPI.rf.4[i] <- sum(Rrf4)/nrow(test)
} #for loop
MSE.lr.4
MSE.rf.4
meanWidth.lr.4
meanWidth.rf.4
percntPI.lr.4
percntPI.rf.4







