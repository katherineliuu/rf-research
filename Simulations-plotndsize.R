# Simulation studies to investigate MSPE, width and coverage rate against nodesize
library(randomForestSRC)
library(ggformula)
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

iter <- 10
itertune<- 10
avgResults <- array(NA, dim=c(40,3)) # average MSE, meanWIdth, coverage rate for 4 simulations
allResults <- array(NA, dim=c(4,3, iter))
optimal <- array(NA, dim=c(2,itertune))

###################################################################
# Simulation 1: Linear Model (one variable)
set.seed(06232016)#set seed
for(i in 1:itertune){
  n <- 2000  #number of observations (trainig and test set combined)
  sig <- 5  #error standard deviation
  x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
  e <- rnorm(n, 0, sig) #generate error term 
  y <- 2*x+3 + e  #generate response values 
  data <- data.frame(x, y) #combine explanatory and response varables into data frame
  train <- data[1:(n/2), ] #use first half of observations as training data
  test <- data[(n/2+1):n, ] #use second half of observations as test data
  ## default tuning call
  ot <- tune(y ~ ., train)
  ## use optimized nodesize
  optimal[,i] <- ot$optimal
} # get optimal parameters
optimal ##15,20,25

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
    allResults[,,i] <- getResults(train, test, nodesize = nd, mtry = 1)
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
save("Simulations-plotndsize.Rdata")

#plotting
ndsize<- c(rep(nodesize,4))
method <- c(rep('LR',10), rep('RF-QRF',10), rep('RF-OOB',10),rep('RF-OOBsym',10))
MSPE <- avgResults[,1]
PIWidth <- avgResults[,2]
Cov.Rate <- avgResults[,3]
dataf<- data.frame(ndsize,method,MSPE,PIWidth,Cov.Rate)
p <- ggplot(dataf, aes(x=ndsize,y=MSPE, color=method))
p + geom_point() +geom_line()
p <- ggplot(dataf, aes(x=ndsize,y=PIWidth, color=method))
p + geom_point() +geom_line()
p <- ggplot(dataf, aes(x=ndsize,y=Cov.Rate, color=method))
p + geom_point() +geom_line() + geom_hline(yintercept=0.95)
###################################################################
# Simulation 2: Nonlinear Model (one variable)

set.seed(06232021) #set seed

for(i in 1:itertune){
  n <- 2000  #number of observations (trainig and test set combined)
  sig <- 5  #error standard deviation
  x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
  e <- rnorm(n, 0, sig) #generate error term 
  y <- 0.1*(x-7)^2-3*cos(x) + 5*log(abs(x)) + 3 + e  #generate response values 
  data <- data.frame(x, y) #combine explanatory and response varables into data frame
  train <- data[1:(n/2), ] #use first half of observations as training data
  test <- data[(n/2+1):n, ] #use second half of observations as test data
  ## default tuning call
  ot <- tune(y ~ ., train)
  optimal[,i] <- ot$optimal
}
optimal ##7,8,9,10

j=0
nodesize <- c(1,3,5,7,9,11,13,15,17,19)

for(nd in nodesize){ #loop through a range of nodesizes
  set.seed(0621)
  for (i in 1:10){
    n <- 2000  #number of observations (trainig and test set combined)
    sig <- 5  #error standard deviation
    x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
    e <- rnorm(n, 0, sig) #generate error term 
    y <- 0.1*(x-7)^2-3*cos(x) + 5*log(abs(x)) + 3 + e  #generate response values 
    data <- data.frame(x, y) #combine explanatory and response varables into data frame
    train <- data[1:(n/2), ] #use first half of observations as training data
    test <- data[(n/2+1):n, ] #use second half of observations as test data
    allResults[,,i] <- getResults(train, test, nodesize = nd, mtry = 1)
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

#plotting
ndsize<- c(rep(nodesize,4))
method <- c(rep('LR',10), rep('RF-QRF',10), rep('RF-OOB',10),rep('RF-OOBsym',10))
MSPE <- avgResults[,1]
PIWidth <- avgResults[,2]
Cov.Rate <- avgResults[,3]
dataf<- data.frame(ndsize,method,MSPE,PIWidth,Cov.Rate)
p <- ggplot(dataf, aes(x=ndsize,y=MSPE, color=method))
p + geom_point() +geom_line()
p <- ggplot(dataf, aes(x=ndsize,y=PIWidth, color=method))
p + geom_point() +geom_line()
p <- ggplot(dataf, aes(x=ndsize,y=Cov.Rate, color=method))
p + geom_point() +geom_line() + geom_hline(yintercept=0.95)

#(data = filter(dataf, method%in%c("LR","RF-QRF","RF-OOB","RF-OOBsym")), aes(x=x, y=pred, colour=factor(method, labels=c("LR","RF-QRF","RF-OOB","RF-OOBsym"))), size = 1.2)+
 # scale_colour_manual(breaks = c("LR","RF-QRF","RF-OOB","RF-OOBsym"),values=c("red", "blue", "green","yellow"))+
  #labs(color = "Curve")+ theme(legend.text=element_text(size=12))+ theme(legend.title=element_text(size=12))+theme(legend.position = "bottom")+theme_bw()


save("Simulations-plotndsize.Rdata")
###################################################################
# Simulation 3: Linear Model (multivariate)

set.seed(06232012) #set seed
for(i in 1:itertune){
  n <- 2000 #number of observations (trainig and test set combined)
  xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
  X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
  sig <- 5  #error standard deviation
  e <- rnorm(n, 0, sig) #generate error term 
  y <- 2*X[,1]+3*X[,4]+4*X[,6]-3*X[,7]+ X[,9] + e  #generate response values 
  data <- data.frame(X, y) #combine explanatory and response varables into data frame
  train <- data[1:(n/2), ] #use first half of observations as training data
  test <- data[(n/2+1):n, ] #use second half of observations as test data
  ## default tuning call
  ot <- tune(y ~ ., train)
  optimal[,i] <- ot$optimal
}
optimal #1,222,33,4


j=0
nodesize <- c(1,2,3,4,5,6,7,8,9,10)

for(nd in nodesize){ #loop through a range of nodesizes
  set.seed(0621)
  for (i in 1:10){
    n <- 2000  
    xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
    X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
    sig <- 5  #error standard deviation
    e <- rnorm(n, 0, sig) #generate error term 
    y <- 2*X[,1]+3*X[,4]+4*X[,6]-3*X[,7]+ X[,9] + e  #generate response values 
    data <- data.frame(X, y) #combine explanatory and response varables into data frame
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
save("Simulations-plotndsize.Rdata")

#plotting
ndsize<- c(rep(nodesize,4))
method <- c(rep('LR',10), rep('RF-QRF',10), rep('RF-OOB',10),rep('RF-OOBsym',10))
MSPE <- avgResults[,1]
PIWidth <- avgResults[,2]
Cov.Rate <- avgResults[,3]
dataf<- data.frame(ndsize,method,MSPE,PIWidth,Cov.Rate)
p <- ggplot(dataf, aes(x=ndsize,y=MSPE, color=method))
p + geom_point() +geom_line()
p <- ggplot(dataf, aes(x=ndsize,y=PIWidth, color=method))
p + geom_point() +geom_line()
p <- ggplot(dataf, aes(x=ndsize,y=Cov.Rate, color=method))
p + geom_point() +geom_line() + geom_hline(yintercept=0.95)

###################################################################
# Simulation 4: Nonlinear Model (multivariate)

set.seed(06231) #set seed

for(i in 1:itertune){
  n <- 2000  #number of observations (trainig and test set combined)
  xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
  X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
  sig <- 5  #error standard deviation
  e <- rnorm(n, 0, sig) #generate error term 
  y <- (X[,1]-6)^2 + 12*cos(X[,3]) + (X[,7]-5)*(X[,8]-3) + 0.02*(X[,10]-5)^5+ e  #generate response values 
  data <- data.frame(X, y) #combine explanatory and response varables into data frame
  train <- data[1:(n/2), ] #use first half of observations as training data
  test <- data[(n/2+1):n, ] #use second half of observations as test data
  ## default tuning call
  ot <- tune(y ~ ., train)
  optimal[,i] <- ot$optimal
}
optimal


j=0
nodesize <- c(1,2,3,4,5,6,7,8,9,10)

for(nd in nodesize){ #loop through a range of nodesizes
  set.seed(0621)
  for (i in 1:10){
    n <- 2000  
    xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
    X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
    sig <- 5  #error standard deviation
    e <- rnorm(n, 0, sig) #generate error term 
    y <- (X[,1]-6)^2 + 12*cos(X[,3]) + (X[,7]-5)*(X[,8]-3) + 0.02*(X[,10]-5)^5+ e  #generate response values 
    data <- data.frame(X, y) #combine explanatory and response varables into data frame
    train <- data[1:(n/2), ] #use first half of observations as training data
    test <- data[(n/2+1):n, ] #use second half of observations as test data
    allResults[,,i] <- getResults(train, test, nodesize = nd, mtry = 8)
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
save("Simulations-plotndsize.Rdata")

#plotting
ndsize<- c(rep(nodesize,4))
method <- c(rep('LR',10), rep('RF-QRF',10), rep('RF-OOB',10),rep('RF-OOBsym',10))
MSPE <- avgResults[,1]
PIWidth <- avgResults[,2]
Cov.Rate <- avgResults[,3]
dataf<- data.frame(ndsize,method,MSPE,PIWidth,Cov.Rate)
p <- ggplot(dataf, aes(x=ndsize,y=MSPE, color=method))
p + geom_point() +geom_line()
p <- ggplot(dataf, aes(x=ndsize,y=PIWidth, color=method))
p + geom_point() +geom_line()
p <- ggplot(dataf, aes(x=ndsize,y=Cov.Rate, color=method))
p + geom_point() +geom_line() + geom_hline(yintercept=0.95)
