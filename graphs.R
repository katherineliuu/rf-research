# Simulation studies to investigate RF and LM prediction intervals
library(randomForestSRC)
library(ggformula)
###################################################################
# Simulation 1: Linear Model (one variable)

set.seed(06232019)#set seed


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
  yhat.train <- predict(M, newdata=train)
  PI<- predict(M, newdata=test, interval="prediction", level=0.95) # prediction intervals with 95%
  #head(PI)
  
  summary(M)
  
  #RF
  o <- quantreg(y ~ ., train, splitrule = "mse", nodesize = 10)  ## quantile regression with mse splitting
  o.test <- quantreg(object=o, newdata=test)   # test call
  o.train <- quantreg(object=o, newdata=train) 
  #head(o.test$predicted)  #predictions
  
  quant.dat <- get.quantile(o.test, c(.025, .975))
  #head(quant.dat)    #prediction intervals
  o.test1 <- quantreg(object=o, newdata=test)

  fun.1<- function(x) 2*x+3

#gf_point(data=data, y~x) %>% gf_lm() %>% gf_labs("", x="") %>% gf_function(fun = fun.1, col="green")
#plot(test$x, o.test$predicted, xlim=range(data$x), ylim=range(o.test$predicted), xlab="x", ylab="y", main = "",pch=16)
#lines( test$x[order(test$x)],  o.test$predicted[order(test$x)], xlim=range(test$x), ylim=range(o.test$predicted), pch=16)
  x<- rep(train$x, 3)
  y<-rep(train$y, 3)
  method <- c(rep('Truef',1000), rep('RF',1000), rep('LR',1000))
  pred <- c(fun.1(train$x), o$predicted.oob, yhat.train)
  dataf<- data.frame(x,y,method,pred)

  p <- ggplot(train, aes(x=x,y=y))
  p + geom_point() +geom_line(data = filter(dataf, method%in%c("Truef")), aes(x=x, y=pred, colour=factor(method, labels=c("Truef"))), size = 1.2)+
    scale_colour_manual(breaks = c("Truef"),values=c("green"))+
    labs(color = "Curve")+ theme(legend.text=element_text(size=12))+ theme(legend.title=element_text(size=12))+theme(legend.position = "bottom")+theme_bw()
  
  p + geom_point() +geom_line(data = filter(dataf, method%in%c("LR","Truef")), aes(x=x, y=pred, colour=factor(method, labels=c("LR","Truef"))), size = 1.2)+
    scale_colour_manual(breaks = c("LR","Truef"),values=c("red", "green"))+
    labs(color = "Curve")+ theme(legend.text=element_text(size=12))+ theme(legend.title=element_text(size=12))+theme(legend.position = "bottom")+theme_bw()
  
  p + geom_point() +geom_line(data = filter(dataf, method%in%c("LR","RF","Truef")), aes(x=x, y=pred, colour=factor(method, labels=c("LR","RF","Truef"))), size = 1.2)+
    scale_colour_manual(breaks = c("LR","RF","Truef"),values=c("red", "blue", "green"))+
    labs(color = "Curve")+ theme(legend.text=element_text(size=12))+ theme(legend.title=element_text(size=12))+theme(legend.position = "bottom")+theme_bw()

###################################################################
# Simulation 2: Nonlinear Model (one variable)
set.seed(0623201) #set seed

  n <- 2000  #number of observations (trainig and test set combined)
  sig <- 5  #error standard deviation
  x <- runif(n, 0, 10) #generate x variable Unif(0,10) distribution
  e <- rnorm(n, 0, sig) #generate error term 
  y <- 0.1*(x-7)^2-3*cos(x) + 5*log(abs(x)) + 3 + e  #generate response values 
  data <- data.frame(x, y) #combine explanatory and response varables into data frame
  train <- data[1:(n/2), ] #use first half of observations as training data
  test <- data[(n/2+1):n, ] #use second half of observations as test data

    #LR
    
    M <- lm(data=train, y~x)
  yhat <- predict(M, newdata=test) # predictions
  yhat.train <- predict(M, newdata=train)
  PI<- predict(M, newdata=test, interval="prediction", level=0.95) # prediction intervals with 95%
  
  
  #RF
  o <- quantreg(y ~ ., train, splitrule = "mse", nodesize = 10)  ## quantile regression with mse splitting
  o.test <- quantreg(object=o, newdata=test)
  o.train <- quantreg(object=o, newdata=train) 

  # test call
  #head(o.test$predicted)  #predictions
  


fun.2 <- function(x) 0.1*(x-7)^2-3*cos(x) + 5*log(abs(x)) + 3 

  
#gf_point(data=data, y~x) %>% gf_lm() %>% gf_labs("", x="") %>% gf_function(fun = fun.2, col="green") 


x<- rep(train$x, 3)
y<-rep(train$y, 3)
method <- c(rep('Truef',1000), rep('RF',1000), rep('LR',1000))
pred <- c(fun.2(train$x), o$predicted.oob, yhat.train)
dataf<- data.frame(x,y,method,pred)

p <- ggplot(train, aes(x=x,y=y))
p + geom_point() +geom_line(data = filter(dataf, method%in%c("Truef")), aes(x=x, y=pred, colour=factor(method, labels=c("Truef"))), size = 1.2)+
  scale_colour_manual(breaks = c("Truef"),values=c("green"))+
  labs(color = "Curve")+ theme(legend.text=element_text(size=12))+ theme(legend.title=element_text(size=12))+theme(legend.position = "bottom")+theme_bw()

p + geom_point() +geom_line(data = filter(dataf, method%in%c("LR","Truef")), aes(x=x, y=pred, colour=factor(method, labels=c("LR","Truef"))), size = 1.2)+
  scale_colour_manual(breaks = c("LR","Truef"),values=c("red", "green"))+
  labs(color = "Curve")+ theme(legend.text=element_text(size=12))+ theme(legend.title=element_text(size=12))+theme(legend.position = "bottom")+theme_bw()

p + geom_point() +geom_line(data = filter(dataf, method%in%c("LR","RF","Truef")), aes(x=x, y=pred, colour=factor(method, labels=c("LR","RF","Truef"))), size = 1.2)+
  scale_colour_manual(breaks = c("LR","RF","Truef"),values=c("red", "blue", "green"))+
  labs(color = "Curve")+ theme(legend.text=element_text(size=12))+ theme(legend.title=element_text(size=12))+theme(legend.position = "bottom")+theme_bw()