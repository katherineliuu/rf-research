# Simulation studies to investigate how QRF and RFOOB perform. 
library(randomForestSRC)
library(ggformula)
iter<-10
itertune<-10
getResults <- function(train, test, nodesize, mtry) {
  results <-  array(NA, dim=c(4,3)) #MSE, meanWIdth, coverage rate
  #LR
  M <- lm(data=train, y~.)
  yhat <- predict(M, newdata=test) # predictions
  PI<- predict(M, newdata=test, interval="prediction", level=0.95) # prediction intervals with 95%
  results[1,1]  <- mean((yhat - test$y)^2) #MSE
  results[1,2] <- mean(PI[,3]-PI[,2]) #meanWidth
  Rlr <- test$y>=PI[,2] & test$y<=PI[,3] 
  results[1,3] <- sum(Rlr)/nrow(test) #coverage rate
  #RF-QRF
  o <- quantreg(y ~ ., train, splitrule = "mse", nodesize = nodesize, mtry = mtry)  ## quantile regression with mse splitting
  o.test <- quantreg(object=o, newdata=test)   # test call
  quant.dat <- get.quantile(o.test, c(.025, .975))
  results[2,1] <- o.test$err.rate[length(o.test$err.rate)] #MSE
  results[2,2] <- mean(quant.dat[,2]-quant.dat[,1]) #meanWidth
  Rrf <- test$y>=quant.dat[,1] & test$y<=quant.dat[,2] 
  results[2,3] <- sum(Rrf)/nrow(test) #coverage rate
  #RFOOB
  ## symmetry = TRUE leads to the symmetric OOB Intervals
  ## symmetry = FALSE leads to the standard OOB Intervals
  oob_error <- (o$predicted.oob-train$y)
  oob_abs_error = sort(abs(o$predicted.oob-train$y))
  
  Dl = o.test$predicted - quantile(oob_abs_error, 0.95) # alpha=0.05, 1-alpha
  Du = o.test$predicted + quantile(oob_abs_error, 0.95)
  results[4,1]<-results[2,1]
  results[4,2] <- 2*quantile(oob_abs_error, 0.95) #meanWidth
  Rrf.oob.symm <- test$y>=Dl & test$y<=Du 
  results[4,3] <- sum(Rrf.oob.symm)/nrow(test) #coverage rate
  
  Dl = o.test$predicted + quantile(oob_error, 0.025) # alpha/2
  Du = o.test$predicted + quantile(oob_error, 0.975) # 1-alpha/2
  results[3,1]<-results[2,1]
  results[3,2] <- mean(Du-Dl) #meanWidth
  Rrf.oob <- test$y>=Dl & test$y<=Du 
  results[3,3] <- sum(Rrf.oob)/nrow(test) #coverage rate
  
  return(results)
}  ##  function
avgResults <- array(NA, dim=c(40,3)) # average MSE, meanWIdth, coverage rate for 4 simulations
allResults <- array(NA, dim=c(4,3, iter))
optimal <- array(NA, dim=c(2,itertune))


###################################################################
# Simulation 1: heteroscedastic-symmetric (normal)

set.seed(06232019)#set seed

for (i in 1:itertune){
  n <- 2000  #number of observations (trainig and test set combined)
  xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
  X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
  e <- rnorm(n, mean = 0, sd =0.2*abs(X[,1]^2+X[,3])) #generate heterosedastic error term 
  y <- (X[,1]-6)^2 + 12*cos(X[,3]) + e  #generate response values  (nonlinear)
  data <- data.frame(X, y) #combine explanatory and response varables into data frame
  train <- data[1:(n/2), ] #use first half of observations as training data
  test <- data[(n/2+1):n, ] #use second half of observations as test data
  ## default tuning call
  #ot <- tune(y ~ ., train)
  ## use optimized nodesize
  #optimal[,i] <- ot$optimal
}
optimal

plot(X[,1], y)
plot(X[,3], y)
hist(e)
j=0
nodesize <- c(1,5,10,15,20,25,30,35,40,45)

for(nd in nodesize){ #loop through a range of nodesizes
  set.seed(062)
  for (i in 1:10){
    n <- 2000  #number of observations (trainig and test set combined)
    sig <- 5  #error standard deviation
    x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
    e <- rnorm(n, 0, sig) #generate error term 
    y <- 2*x+3 + e  #generate response values 
    data <- data.frame(x, y) #combine explanatory and response varables into data frame
    train <- data[1:(n/2), ] #use first half of observations as training data
    test <- data[(n/2+1):n, ] #use second half of observations as test data
    allResults[,,i] <- getResults(train, test, nodesize = nd, mtry = 10)
  }
  j = j+1
  avgResults[j,1]<- mean(allResults[1,1,])
  avgResults[j,2]<- mean(allResults[1,2,])
  avgResults[j,3]<- mean(allResults[1,3,])
  avgResults[j+10,1 ] <- mean(allResults[2,1,])
  avgResults[j+10,2 ] <- mean(allResults[2,2,])
  avgResults[j+10,3 ] <- mean(allResults[2,3,])
  avgResults[j+20,1 ] <- mean(allResults[3,1,])
  avgResults[j+20,2 ] <- mean(allResults[3,2,])
  avgResults[j+20,3 ] <- mean(allResults[3,3,])
  avgResults[j+30,1 ] <- mean(allResults[4,1,])
  avgResults[j+30,2 ] <- mean(allResults[4,2,])
  avgResults[j+30,3 ] <- mean(allResults[4,3,])
}
avgResults
save("MoreSimStudy1.Rdata")


# Simulation 2: heteroscedastic-symmetric (non-normal, T distribution)

set.seed(06232019)#set seed

for (i in 1:itertune){
  n <- 2000  #number of observations (trainig and test set combined)
  df=5
  xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
  X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
  e <- rt(n, df)*sqrt(0.3*(abs(X[,1]^3+X[,3])) * (df-2)/df)
  
  #generate heterosedastic error term 
  y <- (X[,1]-6)^2 + 12*cos(X[,3]) + e  #generate response values  (nonlinear)
  data <- data.frame(X, y) #combine explanatory and response varables into data frame
  train <- data[1:(n/2), ] #use first half of observations as training data
  test <- data[(n/2+1):n, ] #use second half of observations as test data #
  # default tuning call
  #ot <- tune(y ~ ., train)
  ## use optimized nodesize
  #optimal[,i] <- ot$optimal
}
optimal

#plot(X[,1], y)
#plot(X[,3], y)
#hist(e)
j=0
nodesize <- c(1,2,3,4,5,6,7,8,9,10)

for(nd in nodesize){ #loop through a range of nodesizes
  set.seed(062)
  for (i in 1:10){
    n <- 2000  #number of observations (trainig and test set combined)
    sig <- 5  #error standard deviation
    x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
    e <- rnorm(n, 0, sig) #generate error term 
    y <- 2*x+3 + e  #generate response values 
    data <- data.frame(x, y) #combine explanatory and response varables into data frame
    train <- data[1:(n/2), ] #use first half of observations as training data
    test <- data[(n/2+1):n, ] #use second half of observations as test data
    allResults[,,i] <- getResults(train, test, nodesize = nd, mtry = 10)
  }
  j = j+1
  avgResults[j,1]<- mean(allResults[1,1,])
  avgResults[j,2]<- mean(allResults[1,2,])
  avgResults[j,3]<- mean(allResults[1,3,])
  avgResults[j+10,1 ] <- mean(allResults[2,1,])
  avgResults[j+10,2 ] <- mean(allResults[2,2,])
  avgResults[j+10,3 ] <- mean(allResults[2,3,])
  avgResults[j+20,1 ] <- mean(allResults[3,1,])
  avgResults[j+20,2 ] <- mean(allResults[3,2,])
  avgResults[j+20,3 ] <- mean(allResults[3,3,])
  avgResults[j+30,1 ] <- mean(allResults[4,1,])
  avgResults[j+30,2 ] <- mean(allResults[4,2,])
  avgResults[j+30,3 ] <- mean(allResults[4,3,])
}
avgResults
save("MoreSimStudy2.Rdata")


###################################################################
# Simulation 3: heteroscedastic-asymmetric (with mean 0)
set.seed(06232019)#set seed

for (i in 1:itertune){
  n <- 2000  #number of observations (trainig and test set combined)
  df=5
  xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
  X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
  e <- (rexp(n, rate = 2)-1/2)*0.1*abs(X[,1]^3)
  
  #generate heterosedastic error term 
  y <- (X[,1]-6)^2 + 12*cos(X[,3]) + e  #generate response values  (nonlinear)
  data <- data.frame(X, y) #combine explanatory and response varables into data frame
  train <- data[1:(n/2), ] #use first half of observations as training data
  test <- data[(n/2+1):n, ] #use second half of observations as test data
  # default tuning call
  #ot <- tune(y ~ ., train)
  ## use optimized nodesize
  #optimal[,i] <- ot$optimal
}
optimal

plot(X[,1], y)
plot(X[,3], y)
hist(e)

j=0
nodesize <- c(1,2,3,4,5,6,7,8,9,10)

for(nd in nodesize){ #loop through a range of nodesizes
  set.seed(062)
  for (i in 1:10){
    n <- 2000  #number of observations (trainig and test set combined)
    sig <- 5  #error standard deviation
    x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
    e <- rnorm(n, 0, sig) #generate error term 
    y <- 2*x+3 + e  #generate response values 
    data <- data.frame(x, y) #combine explanatory and response varables into data frame
    train <- data[1:(n/2), ] #use first half of observations as training data
    test <- data[(n/2+1):n, ] #use second half of observations as test data
    allResults[,,i] <- getResults(train, test, nodesize = nd, mtry = 10)
  }
  j = j+1
  avgResults[j,1]<- mean(allResults[1,1,])
  avgResults[j,2]<- mean(allResults[1,2,])
  avgResults[j,3]<- mean(allResults[1,3,])
  avgResults[j+10,1 ] <- mean(allResults[2,1,])
  avgResults[j+10,2 ] <- mean(allResults[2,2,])
  avgResults[j+10,3 ] <- mean(allResults[2,3,])
  avgResults[j+20,1 ] <- mean(allResults[3,1,])
  avgResults[j+20,2 ] <- mean(allResults[3,2,])
  avgResults[j+20,3 ] <- mean(allResults[3,3,])
  avgResults[j+30,1 ] <- mean(allResults[4,1,])
  avgResults[j+30,2 ] <- mean(allResults[4,2,])
  avgResults[j+30,3 ] <- mean(allResults[4,3,])
}
avgResults
save("MoreSimStudy3.Rdata")


###################################################################
# Simulation 4: homoscedastic-symmetric(Normal)
set.seed(06232019)#set seed

for (i in 1:itertune){
  n <- 2000  #number of observations (trainig and test set combined)
  xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
  X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
  e <- rnorm(n, mean=0, 5) #generate homosedastic error term that follows T distribution
  y <- (X[,1]-6)^2 + 12*cos(X[,3]) + e  #generate response values  (nonlinear)
  data <- data.frame(X, y) #combine explanatory and response varables into data frame
  train <- data[1:(n/2), ] #use first half of observations as training data
  test <- data[(n/2+1):n, ] #use second half of observations as test data
  # default tuning call
  ot <- tune(y ~ ., train)
  ## use optimized nodesize
  optimal[,i] <- ot$optimal
}
optimal
plot(X[,1], y)
plot(X[,3], y)
hist(e)

j=0
nodesize <- c(1,2,3,4,5,6,7,8,9,10)

for(nd in nodesize){ #loop through a range of nodesizes
  set.seed(062)
  for (i in 1:10){
    n <- 2000  #number of observations (trainig and test set combined)
    sig <- 5  #error standard deviation
    x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
    e <- rnorm(n, 0, sig) #generate error term 
    y <- 2*x+3 + e  #generate response values 
    data <- data.frame(x, y) #combine explanatory and response varables into data frame
    train <- data[1:(n/2), ] #use first half of observations as training data
    test <- data[(n/2+1):n, ] #use second half of observations as test data
    allResults[,,i] <- getResults(train, test, nodesize = nd, mtry = 10)
  }
  j = j+1
  avgResults[j,1]<- mean(allResults[1,1,])
  avgResults[j,2]<- mean(allResults[1,2,])
  avgResults[j,3]<- mean(allResults[1,3,])
  avgResults[j+10,1 ] <- mean(allResults[2,1,])
  avgResults[j+10,2 ] <- mean(allResults[2,2,])
  avgResults[j+10,3 ] <- mean(allResults[2,3,])
  avgResults[j+20,1 ] <- mean(allResults[3,1,])
  avgResults[j+20,2 ] <- mean(allResults[3,2,])
  avgResults[j+20,3 ] <- mean(allResults[3,3,])
  avgResults[j+30,1 ] <- mean(allResults[4,1,])
  avgResults[j+30,2 ] <- mean(allResults[4,2,])
  avgResults[j+30,3 ] <- mean(allResults[4,3,])
}
avgResults
save("MoreSimStudy4.Rdata")

###################################################################
# Simulation 5: homoscedastic-symmetric(T)
set.seed(06232019)#set seed

for (i in 1:itertune){
  n <- 2000  #number of observations (trainig and test set combined)
  xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
  X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
  e <- rt(n, df=df) #generate homosedastic error term that follows T distribution
  y <- (X[,1]-6)^2 + 12*cos(X[,3]) + e  #generate response values  (nonlinear)
  data <- data.frame(X, y) #combine explanatory and response varables into data frame
  train <- data[1:(n/2), ] #use first half of observations as training data
  test <- data[(n/2+1):n, ] #use second half of observations as test data
  # default tuning call
  ot <- tune(y ~ ., train)
  ## use optimized nodesize
  optimal[,i] <- ot$optimal
}
optimal

plot(X[,1], y)
plot(X[,3], y)
hist(e)

j=0
nodesize <- c(1,2,3,4,5,6,7,8,9,10)

for(nd in nodesize){ #loop through a range of nodesizes
  set.seed(062)
  for (i in 1:10){
    n <- 2000  #number of observations (trainig and test set combined)
    sig <- 5  #error standard deviation
    x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
    e <- rnorm(n, 0, sig) #generate error term 
    y <- 2*x+3 + e  #generate response values 
    data <- data.frame(x, y) #combine explanatory and response varables into data frame
    train <- data[1:(n/2), ] #use first half of observations as training data
    test <- data[(n/2+1):n, ] #use second half of observations as test data
    allResults[,,i] <- getResults(train, test, nodesize = nd, mtry = 10)
  }
  j = j+1
  avgResults[j,1]<- mean(allResults[1,1,])
  avgResults[j,2]<- mean(allResults[1,2,])
  avgResults[j,3]<- mean(allResults[1,3,])
  avgResults[j+10,1 ] <- mean(allResults[2,1,])
  avgResults[j+10,2 ] <- mean(allResults[2,2,])
  avgResults[j+10,3 ] <- mean(allResults[2,3,])
  avgResults[j+20,1 ] <- mean(allResults[3,1,])
  avgResults[j+20,2 ] <- mean(allResults[3,2,])
  avgResults[j+20,3 ] <- mean(allResults[3,3,])
  avgResults[j+30,1 ] <- mean(allResults[4,1,])
  avgResults[j+30,2 ] <- mean(allResults[4,2,])
  avgResults[j+30,3 ] <- mean(allResults[4,3,])
}
avgResults
save("MoreSimStudy5.Rdata")

###################################################################
# Simulation 6: homoscedastic-asymmetric

set.seed(06232019)#set seed

for (i in 1:itertune){
  n <- 2000  #number of observations (trainig and test set combined)
  df=5
  xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
  X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
  e <- rexp(n, rate = 2)-1/2
  
  #generate heterosedastic error term 
  y <- (X[,1]-6)^2 + 12*cos(X[,3]) + e  #generate response values  (nonlinear)
  data <- data.frame(X, y) #combine explanatory and response varables into data frame
  train <- data[1:(n/2), ] #use first half of observations as training data
  test <- data[(n/2+1):n, ] #use second half of observations as test data
  # default tuning call
  ot <- tune(y ~ ., train)
  ## use optimized nodesize
  optimal[,i] <- ot$optimal
}
optimal

plot(X[,1], y)
plot(X[,3], y)
hist(e)

j=0
nodesize <- c(1,2,3,4,5,6,7,8,9,10)

for(nd in nodesize){ #loop through a range of nodesizes
  set.seed(062)
  for (i in 1:10){
    n <- 2000  #number of observations (trainig and test set combined)
    sig <- 5  #error standard deviation
    x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
    e <- rnorm(n, 0, sig) #generate error term 
    y <- 2*x+3 + e  #generate response values 
    data <- data.frame(x, y) #combine explanatory and response varables into data frame
    train <- data[1:(n/2), ] #use first half of observations as training data
    test <- data[(n/2+1):n, ] #use second half of observations as test data
    allResults[,,i] <- getResults(train, test, nodesize = nd, mtry = 10)
  }
  j = j+1
  avgResults[j,1]<- mean(allResults[1,1,])
  avgResults[j,2]<- mean(allResults[1,2,])
  avgResults[j,3]<- mean(allResults[1,3,])
  avgResults[j+10,1 ] <- mean(allResults[2,1,])
  avgResults[j+10,2 ] <- mean(allResults[2,2,])
  avgResults[j+10,3 ] <- mean(allResults[2,3,])
  avgResults[j+20,1 ] <- mean(allResults[3,1,])
  avgResults[j+20,2 ] <- mean(allResults[3,2,])
  avgResults[j+20,3 ] <- mean(allResults[3,3,])
  avgResults[j+30,1 ] <- mean(allResults[4,1,])
  avgResults[j+30,2 ] <- mean(allResults[4,2,])
  avgResults[j+30,3 ] <- mean(allResults[4,3,])
}
avgResults
save("MoreSimStudy6.Rdata")
