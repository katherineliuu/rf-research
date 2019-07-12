# Simulation studies to investigate RF and LM prediction intervals
library(randomForestSRC)
avgResults <-  array(NA, dim=c(8,3)) # average MSE, meanWIdth, coverage rate for 4 simulations

RFOOBInterval <- function(x, 
                          y, 
                          x0, 
                          ntree = 1000, 
                          alpha = 0.10, 
                          symmetry = TRUE,
                          mtry = if (!is.null(y) && !is.factor(y))
                            max(floor(ncol(x)/3), 1) else floor(sqrt(ncol(x))), 
                          nodesize = if (!is.null(y) && !is.factor(y)) 5 else 1){
  x <- as.matrix(x)
  x0 <- as.matrix(x0)
  colnames(x) <- 1:ncol(x)
  rownames(x) <- 1:nrow(x)
  colnames(x0) <- 1:ncol(x0)
  rownames(x0) <- 1:nrow(x0)
  
  n = nrow(x)
  ntest = nrow(x0)
  
  rf = randomForest(x=x, y=y, ntree=ntree, mtry = mtry, nodesize = nodesize,
                    keep.forest=TRUE, keep.inbag=TRUE)
  
  test_pred <- predict(rf, x0)
  oob_abs_error = sort(abs(y - rf$predicted))
  oob_error = sort(y - rf$predicted)
  
  upper_pred = rep(NA, ntest)
  lower_pred = rep(NA, ntest)
  
  ## symmetry = TRUE leads to the symmetric OOB Intervals
  ## symmetry = FALSE leads to the standard OOB Intervals
  if(symmetry){
    for (i in 1:ntest){
      upper_pred[i] = test_pred[i] + quantile(oob_abs_error,1-alpha)
      lower_pred[i] = test_pred[i] - quantile(oob_abs_error,1-alpha)
    }
  }
  else{
    for (i in 1:ntest){
      upper_pred[i] = test_pred[i] + quantile(oob_error, 1-alpha/2)
      lower_pred[i] = test_pred[i] + quantile(oob_error, alpha/2)
    }
    
  }
  
  return(list(pred = matrix(test_pred,ntest,1),
              lo = matrix(lower_pred,ntest,1), 
              up = matrix(upper_pred,ntest,1), 
              fit = matrix(predict(rf,x),n,1)))
} #### function OOB prediction interval 




getResults <- function(train, test, nodesize) {
  results <-  array(NA, dim=c(2,3)) #MSE, meanWIdth, coverage rate
  #LR
  M <- lm(data=train, y~.)
  yhat <- predict(M, newdata=test) # predictions
  PI<- predict(M, newdata=test, interval="prediction", level=0.95) # prediction intervals with 95%
  results[1,1]  <- mean((yhat - test$y)^2) #MSE
  results[1,2] <- mean(PI[,3]-PI[,2]) #meanWidth
  Rlr <- test$y>=PI[,2] & test$y<=PI[,3] 
  results[1,3] <- sum(Rlr)/nrow(test) #percentPI
  #RF
  o <- quantreg(y ~ ., train, splitrule = "mse", nodesize = nodesize)  ## quantile regression with mse splitting
  o.test <- quantreg(object=o, newdata=test)   # test call
  quant.dat <- get.quantile(o.test, c(.025, .975))
  results[2,1] <- o.test$err.rate[n/2] #MSE
  results[2,2] <- mean(quant.dat[,2]-quant.dat[,1]) #meanWidth
  Rrf <- test$y>=quant.dat[,1] & test$y<=quant.dat[,2] 
  results[2,3] <- sum(Rrf)/nrow(test) #percentPI
  return(results)
}  ##function

iter <- 10
allResults <- array(NA, dim=c(2,3, iter))
OOBResults <- array(NA, dim=c(2,3,iter)) 
###################################################################
# Simulation 1: Linear Model (one variable)
set.seed(06232016)#set seed

for(i in 1:iter){
  n <- 100 #number of observations (trainig and test set combined)
  sig <- 5  #error standard deviation
  x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
  e <- rnorm(n, 0, sig) #generate error term 
  y <- 2*x+3 + e  #generate response values 
  data <- data.frame(x, y) #combine explanatory and response varables into data frame
  train <- data[1:(n/2), ] #use first half of observations as training data
  test <- data[(n/2+1):n, ] #use second half of observations as test data
  allResults[,,i] <- getResults(train, test, 10)
  OOBResults[,,i] <- RFOOBInterval(x= train$x, y=train$y, x0 = test$x,  symmetry = TRUE)
}
  allResults1
  OOBResults
  avgResults[1,1] <- mean(allResults1[1,1,])
  avgResults[1,2] <- mean(allResults1[1,2,])
  avgResults[1,3] <- mean(allResults1[1,3,])
  avgResults[2,1] <- mean(allResults1[2,1,])
  avgResults[2,2] <- mean(allResults1[2,2,])
  avgResults[2,3] <- mean(allResults1[2,3,])
  avgResults
  
  
  
  ###################################################################
  # Simulation 2: Nonlinear Model (one variable)
  
    set.seed(06232021) #set seed
  
  for(i in 1:iter){
    n <- 2000  #number of observations (trainig and test set combined)
    sig <- 5  #error standard deviation
    x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
    e <- rnorm(n, 0, sig) #generate error term 
    y <- 0.1*(x-7)^2-3*cos(x) + 5*log(abs(x)) + 3 + e  #generate response values 
    data <- data.frame(x, y) #combine explanatory and response varables into data frame
    train <- data[1:(n/2), ] #use first half of observations as training data
    test <- data[(n/2+1):n, ] #use second half of observations as test data
    allResults[,,iter] <- getResults(train, test, 10)
  }
  allResults2
  avgResults[3,1] <- mean(allResults2[1,1,])
  avgResults[3,2] <- mean(allResults2[1,2,])
  avgResults[3,3] <- mean(allResults2[1,3,])
  avgResults[4,1] <- mean(allResults2[2,1,])
  avgResults[4,2] <- mean(allResults2[2,2,])
  avgResults[4,3] <- mean(allResults2[2,3,])
  avgResults
    
    ###################################################################
    # Simulation 3: Linear Model (multivariate)
    
    set.seed(06232012) #set seed
  for(i in 1:iter){
      n <- 2000  #number of observations (trainig and test set combined)
      xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
      X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
      sig <- 5  #error standard deviation
      e <- rnorm(n, 0, sig) #generate error term 
      y <- 2*X[,1]+3*X[,4]+4*X[,6]-3*X[,7]+ X[,9] + e  #generate response values 
      data <- data.frame(X, y) #combine explanatory and response varables into data frame
      train <- data[1:(n/2), ] #use first half of observations as training data
      test <- data[(n/2+1):n, ] #use second half of observations as test data
      allResults[,,iter] <- getResults(train, test, 1)
      }
    allResults3
    avgResults[5,1] <- mean(allResults3[1,1,])
    avgResults[5,2] <- mean(allResults3[1,2,])
    avgResults[5,3] <- mean(allResults3[1,3,])
    avgResults[6,1] <- mean(allResults3[2,1,])
    avgResults[6,2] <- mean(allResults3[2,2,])
    avgResults[6,3] <- mean(allResults3[2,3,])
    avgResults
      
      ###################################################################
      # Simulation 4: Nonlinear Model (multivariate)
      
      set.seed(06231) #set seed

for(i in 1:iter){
        n <- 2000  #number of observations (trainig and test set combined)
        xvec <- runif(n*10, 0, 10) #generate observations for 10 explanatory variables
        X <- matrix(xvec, ncol=10) #arrange explanatory variables into matrix
        sig <- 5  #error standard deviation
        e <- rnorm(n, 0, sig) #generate error term 
        y <- (X[,1]-6)^2 + 12*cos(X[,3]) + (X[,7]-5)*(X[,8]-3) + 0.02*(X[,10]-5)^5+ e  #generate response values 
        data <- data.frame(X, y) #combine explanatory and response varables into data frame
        train <- data[1:(n/2), ] #use first half of observations as training data
        test <- data[(n/2+1):n, ] #use second half of observations as test data
        allResults4[,,iter] <- getResults(train, test, 1)
        }
  allResults4
  avgResults[7,1] <- mean(allResults4[1,1,])
  avgResults[7,2] <- mean(allResults4[1,2,])
  avgResults[7,3] <- mean(allResults4[1,3,])
  avgResults[8,1] <- mean(allResults4[2,1,])
  avgResults[8,2] <- mean(allResults4[2,2,])
  avgResults[8,3] <- mean(allResults4[2,3,])
  
  avgResults
  colnames(avgResults) <- c("MSPE", "PIWidth", "CoverageRate")
  rownames(avgResults) <- c("Sim1.LR", "Sim1.RF", 
                            "Sim2.LR", "Sim2.RF",
                            "Sim3.LR", "Sim3.RF",
                            "Sim4.LR", "Sim4.RF" )
  
  #plot for coverage rate
  boxplot(len~supp*dose, notch=TRUE,
          main = "Multiple boxplots for comparision",
          at = c(1,2,3,4),
          names = c("1", "2", "3", "4"),
          col=(c("navy","maroon")))