  ########libraries########
  library(readr)
  library(car)
  library(MASS)
  library(fitdistrplus)
  #########################
  setwd("C:/Users/Christian/Desktop/5. Semester/StatMod/Assignments/Assignment Done")
  
  tuno <- read.table("tuno.txt", header=TRUE, sep=" ")
  
  #######################DESCRIPTIVE STATISTICS######################
  #Summary statistics
  summary(tuno)

  #Graphical presentation
  scatterplotMatrix(~ r.day + month + day  +pow.obs+ ws30 +wd30 , diagonal="boxplot", data = tuno)
  
  #take a closer look at the relevant correlation
  scatterplotMatrix(~pow.obs+ ws30 , diagonal="boxplot", data = tuno)
  
  ##We see correlation between wind speed at ground and the average daily wind power production
  #this makes sense.
  cor(tuno$ws30, tuno$pow.obs)
  
  ##Normalize power production to be STRICTLY between 0 and 1:
  norm.pow <- tuno$pow.obs/5000
  
  #histogram of data w/o normalization and with:
  par(mfrow=c(1,2))
  hist(tuno$pow.obs)
  hist(norm.pow)
  par(mfrow=c(1,1))
  
  #replace pow.obs with normalization version
  tuno$pow.obs<-norm.pow
  ###################################################################
  
  
  ###########################SIMPLE MODELS###########################
  pow<-1
  if(pow == 1){
    y.test<-tuno$pow.obs
  }else{
    y.test<-tuno$ws30
  }
  
  ##Gamma
  l.gamma <- function(y, pars){
    -sum(dgamma(y, shape=pars[1], rate=pars[2], log=TRUE))
  }
  opt.gamma<- nlminb(c(1,1), l.gamma, lower=c(0,0), y=y.test)
  y.gamma <-dgamma(y.test, shape=opt.gamma$par[1],rate=opt.gamma$par[2])
  hist(y.test,prob=TRUE,ylab="Density", xlab="Predicted wind speed 30 meters above ground level, m/s", main="Different fits")
  points(y.test, y.gamma , pch=19, col=2)
  
  
  ##BETA HAS TO BE BETWEEN 0 and 1
  l.beta <- function(y, pars){
    -sum(dbeta(y, shape1=pars[1], shape2=pars[2], log=TRUE))
  }
  opt.beta<- nlminb(c(1,1), l.beta, lower=c(0,0), y=y.test)
  y.beta<-dbeta(y.test, shape1=opt.beta$par[1],shape2=opt.beta$par[2])
  points(y.test, y.beta, pch=19, col=3)
  
  ##log normal
  l.normal <- function(y, pars){
    -sum(dnorm(y, mean=pars[1], sd=pars[2], log=TRUE))
  }
  opt.normal<- nlminb(c(1,1), l.normal, lower=c(0,0), y=y.test)
  y.normal<-dnorm(y.test, mean=opt.normal$par[1],sd=opt.normal$par[2])

  points(y.test, y.normal, pch=19, col=4)

  
  l.weibull <- function(y, pars){
    -sum(dweibull(y, shape=pars[1], scale=pars[2], log=TRUE))
  }
  opt.weibull<- nlminb(c(1,1), l.weibull, lower=c(0,0), y=y.test)
  y.weibull<-dweibull(y.test, shape=opt.weibull$par[1],scale=opt.weibull$par[2])
  points(y.test, y.weibull, pch=19, col=5)
  
  if(pow == 1){
    legend(0.8,3.5,c("Gamma","Beta","Normal", "Weibull"),lty=c(1,1,1,1),lwd=c(2.5,2.5,2.5,2.5),col=c(2,3,4,5))
  }else{
    legend(21, 0.1,c("Gamma", "Normal", "Weibull"),lty=c(1,1,1),lwd=c(2.5,2.5,2.5),col=c(2,4,5))
  }
  
  
  #Gamma is best for ws30 AIC is lowest and normality seems to be best.
  2*opt.gamma$objective+4
  
  #bedste for pow.obs
  2*opt.beta$objective+4
  
  2*opt.normal$objective+4
  
  2*opt.weibull$objective+4
  
  ###################TRANSFORMATIONS FOR POW.OBS############
  ##ASYMMETRICAL TRANSFORMATION Lambda >0
  bc.trans <- function(lambda,y){
    y.l <- (1/lambda)*log((y^lambda)/(1-y^lambda))
    y.diff <- (1/(y*(1-y^lambda)))
    return(list(y=y.l, dy=y.diff))
  }
  
  lp.lambda <- function(lambda, y){
    n <- length(y)
    y.l <- bc.trans(lambda, y)
    ## The variance estimate which max the Lik.
    sigmasq <- 1/n * sum((y.l$y-mean(y.l$y))^2)
    ## The profile log-likelihood (change of variable principle, that's why y.diff)
    return(- (n/2)*log(sigmasq) - (n/2)*log(2*pi)-n/2+sum(log(y.l$dy)))
  }
  
  lambda <- seq(0.01,1,by=0.01)
  lp <- sapply(lambda, lp.lambda, y=tuno$pow.obs)
  plot(lambda, lp-max(lp), type="l")
  lines(range(lambda), -qchisq(0.95,df=1)/2*c(1,1), lty=2, col=2)
  lambda.as1<-optimize(lp.lambda, c(0,1),y=tuno$pow.obs, maximum=TRUE)
  
  #transform data with optimal lambda
  ty.as1<-bc.trans(lambda.as1$maximum, tuno$pow.obs)
  
  
  ##ASYMMETRICAL TRANSFORMATION Lambda in (0,1)
  bc.trans2 <- function(lambda,y){
    y.l <- 2*log((y^lambda)/((1-y)^(1-lambda)))
    return(y.l)
  }
  
  lp.lambda <- function(lambda, y){
    n <- length(y)
    y.l <- bc.trans2(lambda, y)
    ## The variance estimate which max the Lik.
    sigmasq <- 1/n * sum((y.l-mean(y.l))^2)
    y.diff <- lambda/y + (1-lambda)/(1-y)
    ## The profile log-likelihood (change of variable principle, that's why y.diff)
    return(-n/2 * log(sigmasq) - n/2*log(2*pi)-n/2+sum(log(y.diff)))
  }
  
  lambda <- seq(0.01,1,by=0.01)
  lp <- sapply(lambda, lp.lambda, y=tuno$pow.obs)
  plot(lambda, lp-max(lp), type="l")
  lines(range(lambda), -qchisq(0.95,df=1)/2*c(1,1), lty=2, col=2)
  lambda.as2<-optimize(lp.lambda, c(0,1),y=tuno$pow.obs, maximum=TRUE)
  
  ty.as2<-bc.trans2(lambda.as2$maximum, tuno$pow.obs)
  
  ##BOX COX TRANSFORMATION
  bc.trans3 <- function(lambda,y){
    y.l <- (y^lambda-1)/lambda
    if(lambda == 0){
      y.l<-log(y)
    }
    return(y.l)
  }
  
  lp.lambda <- function(lambda, y){
    n <- length(y)
    y.l <- bc.trans3(lambda, y)
    ## The variance estimate which max the Lik.
    sigmasq <- 1/n * sum((y.l-mean(y.l))^2)
    y.diff <- (y^lambda)/y
    ## The profile log-likelihood
    return(-n/2 * log(sigmasq) - n/2*log(2*pi)-n/2+sum(log(y.diff)))
  }
  
  lambda <- seq(0.01,1,by=0.01)
  lp <- sapply(lambda, lp.lambda, y=tuno$pow.obs)
  plot(lambda, lp-max(lp), type="l")
  lines(range(lambda), -qchisq(0.95,df=1)/2*c(1,1), lty=2, col=2)
  lambda.bc<-optimize(lp.lambda, c(0,1),y=tuno$pow.obs, maximum=TRUE)
  
  ty.bc<-bc.trans3(lambda.bc$maximum, tuno$pow.obs)

  #or using R
  #BC<-boxcox(lm(tuno$pow.obs~1), lambda=lambda)
  #BC$x[which(BC$y==max(BC$y))]
  
  ###################################################################
  
  ##########################REGRESSION MODELS########################
  
  #by the scatterplot we can see ws30 looks a bit curvy, polynomial expansion
  #wd30 would make sense that the wind moves a bit more periodic, it's difficult to observe
  #however we will use a fourier expansion for wd30
  scatterplotMatrix(~ pow.obs + ws30 + wd30,diagonal="boxplot", data=tuno)
  #We don't include time such as month and days due to lack of information about them
  
  ##linear regression without transformed or expanded variables.
  y.reg<-lm(tuno$pow.obs~I(tuno$ws30)+I(tuno$wd30))
  par(mfrow=c(2,2))
  plot(y.reg, which=1:4)
  par(mfrow=(c(1,1)))
  residualPlots(y.reg)
  ##from the normal QQ plot we see we have a few outliers, errors seems to be fairly normal distributed
  ##in cooks distance which tells us about leverage and points with large residuals, we see we have some extremes
  #however we can't remove the extremes, for example 25 is the maximum in ws30, which could have been a storm.
  ##residuals seems to follow a pattern and is not independent
  ##variance in scale location doesn't follow homogeneity, since we see a small bump in the graph and a bunch of points collected
  ##from the residual plots we see that there might be some non-linearity going on with the explanatory variables, we need to transform
  #the variables.
  
  #to use a generalized linear model or polynomial regression we make expansions of the data
  #ws30 is expanded as a taylor polynomial and wd30 as a simple fourier series
  
  y.glm <- lm(pow.obs ~ I(ws30) + I(ws30^2)+I(ws30^3)+I(sin(wd30))+I(cos(wd30)), data=tuno)
  par(mfrow=c(2,2))
  plot(y.glm, which=1:4)
  par(mfrow=(c(1,1)))
  #seems to be a slight improvement from the initial linear regression model
  
  residualPlots(y.glm)
  #huge improvement, the prediction variables seems to follow a linear pattern and the residuals seems to be spread fairly equally.
  #the clustering in ws30 is due to outliers.
  
  #by transforming variables we might get more evenly spread out residuals
  y.glm$ty1 <- lm(ty.as1$y ~ I(ws30) + I(ws30^2)+I(ws30^3)+I(sin(wd30))+I(cos(wd30)), data=tuno)
  y.glm$ty2 <- lm(ty.as2 ~ I(ws30) + I(ws30^2)+I(ws30^3)+I(sin(wd30))+I(cos(wd30)), data=tuno)
  y.glm$ty3 <- lm(ty.bc ~ I(ws30) + I(ws30^2)+I(ws30^3)+I(sin(wd30))+I(cos(wd30)), data=tuno)
  AIC(y.glm$ty1)
  AIC(y.glm$ty2)
  AIC(y.glm$ty3)
  
  par(mfrow=c(2,2))
  plot(y.glm$ty1, which=1:4)
  par(mfrow=(c(1,1)))
  
  par(mfrow=c(2,2))
  plot(y.glm$ty2, which=1:4)
  par(mfrow=(c(1,1)))
  
  par(mfrow=c(2,2))
  plot(y.glm$ty3, which=1:4)
  par(mfrow=(c(1,1)))
  
  residualPlots((y.glm$ty3))
  #it can be seen that box cox seems to be the most reasonable transformation, the variance seems to be most homogeneous
  #at BC and the residuals are independent and normal. Cooks distance looks fine as well. The residual plots is also very good.
  
  #############BACKWARD AND FORWARD SELECTION###############
  y.fit<-y.glm$ty3
  
  #bruger f og t test, se wiki
  summary(y.fit)
  drop1(y.fit, test='F')
  #0.05 signifikans niveau, dropper I(ws30^3) størst P-værdi, backward selection
  
  y.fit <- update(y.fit, .~. -I(ws30^3))
 
  summary(y.fit)
  drop1(y.fit, test='F')
  
  y.fit<-update(y.fit, .~. -I(sin(wd30)))
  summary(y.fit)
  drop1(y.fit, test='F')
  
  
  ####Forward
  y.fitf<-lm(ty.bc ~ 1, data=tuno)

  add1(y.fitf, scope=~ I(ws30) + I(ws30^2)+I(ws30^3)+I(sin(wd30))+I(cos(wd30)), data=tuno,test='F')
    #aic is smallest for ws30 and smallest p value, most significant, so we add it to the model
  y.fitf<-update(y.fitf, .~. +I(ws30))
  
  add1(y.fitf, scope=~ I(ws30) + I(ws30^2)+I(ws30^3)+I(sin(wd30))+I(cos(wd30)), data=tuno,test='F')
  
  y.fitf<-update(y.fitf, .~. +I(ws30^2))
  #same procedure
  add1( y.fitf, scope=~ I(ws30) + I(ws30^2)+I(ws30^3)+I(sin(wd30))+I(cos(wd30)), data=tuno,test='F')
  
  y.fitf<-update(y.fitf, .~. +I(cos(wd30)))
  add1( y.fitf, scope=~ I(ws30) + I(ws30^2)+I(ws30^3)+I(sin(wd30))+I(cos(wd30)), data=tuno,test='F')
  #no more significant parameters to add. we get the same as with backward selection
  
  
  par(mfrow=c(2,2))
  plot(y.fit, which=1:4)
  par(mfrow=(c(1,1)))
  
  residualPlots(y.fit)
  confint(y.fit)
  summary(y.fit)
  #still a small down trend in the variance homogeneity, everything else looks fine
  #y.fit is the final model
  
  ##transforming data back
  orig.trans <- function(lambda,ytrans){
    y <- exp(log(lambda*ytrans + 1)/lambda)
    if(lambda == 0){
      y<-exp(ytrans)
    }
    return(y)
  }

  #sst og sse
  #sum(tuno$pow.obs-y.pred(tuno$ws30,tuno$wd30))^2
  #sum(y.pred(tuno$ws30,tuno$wd30)-mean(tuno$pow.obs))^2
  
  newwd<-seq(0, 6, length.out=288)
  newws <- seq(0,24, length.out=288)
  xnew <- data.frame(ws30 = 10, wd30=newwd)
  xnew2<- data.frame(ws30=newws, wd30=6)
  
  y.val<-predict(y.fit,newdata=xnew,interval="confidence", level=0.95)
  y.val2<-predict(y.fit,newdata=xnew,interval="prediction", level=0.95)
  
  y.val3<-predict(y.fit,newdata=xnew2,interval="confidence", level=0.95)
  y.val4<-predict(y.fit,newdata=xnew2,interval="prediction", level=0.95)
  
  y<-orig.trans(lambda.bc$maximum, y.val[,1])
  y.conf1<-orig.trans(lambda.bc$maximum, y.val[,2])
  y.conf2<-orig.trans(lambda.bc$maximum, y.val[,3])
  y.pred1<-orig.trans(lambda.bc$maximum, y.val2[,2])
  y.pred2<-orig.trans(lambda.bc$maximum, y.val2[,3])
  plot(tuno$wd30,tuno$pow.obs*5000,ylab="Avg daily wind power prod. kW", xlab="Wind direction, rad 30 meters",main="Prediction Results, ws30 constant 10 m/s", xlim=c(0,5.8))
  lines(newwd,y*5000)
  lines(newwd,y.conf1*5000,col="blue")
  lines(newwd,y.conf2*5000, col="blue")
  lines(newwd,y.pred1*5000, col="red")
  lines(newwd,y.pred2*5000, col="red")
  legend(0,4800,c("Prediction Interval","Confidence interval","Model Prediction"),lty=c(1,1,1),lwd=c(2.5,2.5,2.5),col=c("red","blue","black"))
  
  
  #predict for constant wind direction
  y.2<-orig.trans(lambda.bc$maximum, y.val3[,1])
  y.conf3<-orig.trans(lambda.bc$maximum, y.val3[,2])
  y.conf4<-orig.trans(lambda.bc$maximum, y.val3[,3])
  y.pred3<-orig.trans(lambda.bc$maximum, y.val4[,2])
  y.pred4<-orig.trans(lambda.bc$maximum, y.val4[,3])
  plot(tuno$ws30,tuno$pow.obs*5000, ylab="Avg daily wind power prod. kW", xlab="Wind Speed, m/s 30 meters", main="Prediction Results, wd30 constant 6 rad"
       , xlim=c(0,23))
  lines(newws,y.2*5000)
  lines(newws,y.conf3*5000,col="blue")
  lines(newws,y.conf4*5000, col="blue")
  lines(newws,y.pred3*5000, col="red")
  lines(newws,y.pred4*5000, col="red")
  legend(0,4500,c("Prediction Interval","Confidence interval","Model Prediction"),lty=c(1,1,1),lwd=c(2.5,2.5,2.5),col=c("red","blue","black"))
  
  
  ###################################################################

  
  