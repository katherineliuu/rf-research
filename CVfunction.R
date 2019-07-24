library(randomForestSRC)
library(ggformula)
getResults <- function(train, test) {
  results <-  array(NA, dim=c(4,3)) #MSE, meanWIdth, coverage rate
  #LR
  M <- lm(data=train, y~.)
  yhat <- predict(M, newdata=test) # predictions
  PI<- predict(M, newdata=test, interval="prediction", level=0.95) # prediction intervals with 95%
  results[1,1]  <- mean((yhat - test$y)^2) #MSE
  results[1,2] <- mean(PI[,3]-PI[,2]) #meanWidth
  Rlr <- test$y>=PI[,2] & test$y<=PI[,3] 
  results[1,3] <- sum(Rlr)/nrow(test) #coverage rate
  #tune
  ot <- tune(y ~ ., train)
  #RF-QRF
  o <- quantreg(y ~ ., train, splitrule = "mse", nodesize = ot$optimal[1], mtry = ot$optimal[2])  ## quantile regression with mse splitting
  o <- quantreg(y ~ ., train, splitrule = "mse", nodesize = 10, mtry = 10)  ## quantile regression with mse splitting
  
    o.test <- quantreg(object=o, newdata=test)   # test call
  #o.test <- predict(object=o, newdata=test)   # test call
  quant.dat <- get.quantile(o.test, c(.025, .975))
  results[2,1] <- o.test$err.rate[length(o.test$err.rate)] #MSE
  results[2,2] <- mean(quant.dat[,2]-quant.dat[,1]) #meanWidth
  Rrf <- test$y>=quant.dat[,1] & test$y<=quant.dat[,2] 
  results[2,3] <- sum(Rrf)/nrow(test) #coverage rate
  #RFOOB
  ## symmetry = TRUE leads to the symmetric OOB Intervals
  ## symmetry = FALSE leads to the standard OOB Intervals
  oob_error <- o$predicted.oob-train$y
  oob_abs_error <- abs(o$predicted.oob-train$y)
  
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


CV <- function(data, nfolds){
  v <- sample(1:nrow(data)) # Randomly assign cases to folds
  Results <-  array(NA, dim=c(4,3,nfolds))
  for(fold in 1:nfolds){
    foldsize <- nrow(data)/nfolds
    test<- data[v[((fold-1)*foldsize+1):(foldsize*fold)], ]
    train <- data[-v[((fold-1)*foldsize+1):(foldsize*fold)], ]
    Results[,,fold]<- getResults(train, test)
  }
  return(Results)
}