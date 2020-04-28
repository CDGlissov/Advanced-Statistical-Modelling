rm(list=ls())
#data
library(gclus)
#Tables
library(xtable);
#Statistics
library(mgcv);
library(car);
library(grid);
library(gridExtra);
library(repr);
library(nlme);
library(glmmTMB)
#Plots
library(ggplot2, quietly = TRUE)
library(ggpubr); 
library(ggfortify);
library(GGally);
library(reshape2)
ggplot2::theme_set(ggplot2::theme_grey())
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
path="C:/Users/Christian/Dropbox/10. Semester/Dataanalyse og Statistisk Modellering/Assignments/Assignment 3"
setwd(path)

concrete <- read.table("concrets.csv", quote="\"",header = TRUE, comment.char="")
concrete$batch<-factor(concrete$batch)
concrete$date<-as.Date(concrete$date)

#### Problem A ####
## 1: Plot the observations of the 7-days and 28-days strength as a function of time.
ggplot(data = concrete, aes(x=date)) + geom_point(aes(y=y7,color = batch,shape = "y7")) + geom_point(aes(y=y28,color=batch,shape = "y28"))  + 
  scale_x_date(date_minor_breaks = "1 month") + 
  labs(shape = "Time",color="Batch") + 
  xlab("Date")+ ylab("Strength")

ggplot(data = concrete,aes(x =air.temp)) + 
  geom_point(aes(y=y7,color = batch,shape = "y7")) + 
  geom_point(aes(y=y28,color=batch,shape = "y28")) + 
  labs(shape = "Time",color="Batch") + xlab("Air Temperature")+ ylab("Strength")

ggplot(data = concrete,aes(x=y7,y=y28)) + geom_point(aes(color=batch))

dfmelt <- melt(concrete, measure.vars=2:3)
ggplot(data=dfmelt, aes(x=batch,y=value, na.rm=TRUE))+ geom_boxplot(aes(fill=variable))+
  geom_point(aes(y=value, group=variable), position = position_dodge(width=0.75))+ 
  xlab("batch") + ylab("Strength")

## 2: Estimate the mean concrete strength within batches.
dfMean<-data.frame(y7Mean = rep(0,5), y28Mean=rep(0,5))
for (i in 1:5) {
  dfMean$y7Mean[i]<-mean(concrete$y7[concrete$batch==i])
  dfMean$y28Mean[i]<-mean(concrete$y28[concrete$batch==i])
}
dfMean

## 3: Formulate a mixed effect model for the 28-days strength.
concrete.lme1 <- lme(y28~y7*air.temp, 
                    random=~1|batch, data=concrete, method = "ML")
concrete.lm1 <- lm(y28~y7*air.temp, data=concrete)
#Fixed part
anova(concrete.lme1,type = "marginal")
#check for random part, lower AIC -> better with random part
anova(concrete.lme1,concrete.lm1)

concrete.lme2<-update(concrete.lme1,~.-y7:air.temp)
anova(concrete.lme2,type = "marginal")

concrete.lme3<-update(concrete.lme2,~.-air.temp)
anova(concrete.lme3,type = "marginal")
concrete.lme3

#fixed.effects(concrete.lme3)
#random.effects(concrete.lme3) 
#VarCorr(concrete.lme3)

#Diag plots
plot(concrete.lme3, y28 ~ fitted(.) | batch, abline = c(0,1))
plot(concrete.lme3, col=concrete$batch)
qqnorm(resid(concrete.lme3,type="p"), col=concrete$batch)
qqline(resid(concrete.lme3,type="p"), col="blue")

## 4: Test for dependency on the air temperature.
#as air temp got removed from the model, there is no statistical dependence

## Problem B

#group = 1...k, i
# j = samples

n_batch = length(unique(concrete$batch))
#ni, length in each batch
ni = matrix(0,ncol=1,nrow=n_batch)
for (i in 1:n_batch){
  ni[i] = length(concrete$batch[concrete$batch == i])
}


# mu (eq. 5.85) group mean
xi_p = matrix(0,ncol=2,nrow=5)
for(i in 1:n_batch){
  for(j in 1:ni[i]){
    temp1=concrete$y7[concrete$batch == i][j]/ni[i]
    temp2=concrete$y28[concrete$batch == i][j]/ni[i]
    xi_p[i,1]=xi_p[i,1]+temp1
    xi_p[i,2]=xi_p[i,2]+temp2
  }
}


#SSE (eq (5.87)), mu_i (eq (5.86))
SSE <- matrix(0,ncol=2,nrow=2)
mui <- matrix(0,ncol=5,nrow=2)
batch <-1:5
for (i in 1:5){
  mui[,i] <- c(mean(concrete$y7[concrete$batch==batch[i]]),
               mean(concrete$y28[concrete$batch==batch[i]]))
  
  SSE <- SSE + (ni[i]-1)*var(cbind(concrete$y7[concrete$batch==batch[i]],
                                   concrete$y28[concrete$batch==batch[i]]))
  
}
SSE

#SST (eq. (5.89))
SST<-(sum(ni)-1) * var(cbind(concrete$y7,
                    concrete$y28))


SSB = matrix(0,ncol = 2,nrow = 2)
for (i in 1:5){
  SSB<- SSB + ni[i] *(mui[,i]-mu)%*%t(mui[,i]-mu)
}
SSB


N <- sum(ni)
k = n_batch
n0 <- (N-sum(ni^2)/N)/(k-1)
Sigma.m <- SSE/(N-k)
Sigma0.m <- (SSB/(k-1) - Sigma.m)/n0

#Estimate the correlation between the levels of early and 28-day concrete
#strength, including a confidence interval. HOW?

ev <- eigen(Sigma.m)
# extract components
(values <- ev$values)


ev <- eigen(Sigma0.m)
# extract components
(values <- ev$values)

##### Part 2
## Problem A

clo <- read.csv("dat_count3.csv", sep=";")
clo$subjId = as.factor(clo$subjId)
clo$day=as.factor(clo$day)
clo$sex=as.factor(clo$sex)
summary(clo)
str(clo)


fit.glm.bin1=glmmTMB(cbind(clo,nobs-clo)~time+sex+tOut+tInOp+(1|subjId), 
                    family=binomial(link = "logit"), data=clo)

fit.glm.bin2=glmmTMB(cbind(clo,nobs-clo)~time+sex+tOut+tInOp+ar1(as.factor(day)-1|subjId), 
        family=binomial(link = "logit"), data=clo)
fit.glm.bin1
fit.glm.bin2

fit.glm.bin2=glmmTMB(cbind(clo,nobs-clo)~time+sex+tOut+tInOp+ar1(as.factor(day)-1|subjId), 
                     family=binomial(link = "logit"), data=clo)

#t?nk over offset(log(time))
fit.glm.pois<-glmmTMB(clo~time+nobs+sex+tOut+tInOp,family=poisson(link="log"),
                 data=clo)
predict(fit.glm.pois)
#lav plots

### Part B ###
library(numDeriv)

X <- matrix(0,ncol=6,nrow=dim(clo)[1])
X[,1] = 1
X[,2] = clo$time
X[,3] = clo$nobs
X[,4] = clo$sex
X[,5] = clo$tOut
X[,6] = clo$tInOp

nll <- function(u,beta,sigma.u,X){
  # One u for every subject that are the same for each
  U <- u[clo$subjId]
  mu <- exp(X%*%beta)
  -sum(dpois(clo$clo,lambda=mu*U,log=TRUE)
         + dgamma(u, shape=sigma.u, scale=1/sigma.u, log=TRUE))
}

nll.LA <- function(theta,X){
  beta <- theta[1:dim(X)[2]]
  sigma.u <- exp(theta[dim(X)[2]+1])
  est <- nlminb(exp(rep(0,47)), objective = nll,
                beta=beta, sigma.u=sigma.u,X=X,lower = rep(0.0000001,47))
  u <- est$par
  l.u <- est$objective
  H <- hessian(func = nll, x = u, beta = beta, sigma.u = sigma.u,
               X=X)
  l.u + 0.5 * log(det(H/(2*pi)))
}

system.time((fit1 <- nlminb(c(rep(0,6),1),nll.LA,X=X)))
fit1 

nll.LA2 <- function(theta,X){
  beta <- theta[1:dim(X)[2]]
  sigma.u <- exp(theta[dim(X)[2]+1])
  fun.tmp <- function(ui,u,beta,sigma.u,X,i){
    u <- exp(u*0)
    u[i]<-ui
    nll(u,beta,sigma.u,X)
  }
  u <- numeric(47)
  ## Use grouping structure
  for(i in 1:length(u)){
    u[i] <- nlminb(1,objective = fun.tmp, u=u,beta=beta,
                   sigma.u=sigma.u,X=X,i=i, lower=rep(0,47))$par
  }
  l.u <- nll(u,beta,sigma.u,X)
  H <- numeric(length(u))
  for(i in 1:length(u)){
    H[i] <- hessian(func = fun.tmp, x = u[i],u=u, beta = beta, 
                    sigma.u = sigma.u, X=X,i=i)}
  l.u + 0.5 * log(prod(H/(2*pi)))
}

system.time((fit2 <- nlminb(c(rep(0,6),1),nll.LA2,X=X)))
fit2

## Check accuracy by importance sampling

nll.LA4Rexsamp <- function(theta,X,k,seed){
  set.seed(seed)
  beta <- theta[1:dim(X)[2]]
  sigma.u <- exp(theta[dim(X)[2]+1])
  est <- nlminb(rep(1,47),objective = nll,
                beta=beta, sigma.u=sigma.u,X=X, lower=0)
  u <- est$par
  l.u <- est$objective
  H <- diag(hessian(func = nll, x = u, beta = beta, sigma.u = sigma.u,
                    X=X))
  L <- numeric(k)
  s <- sqrt(1/H)
  for(i in 1:k){## Simulation part
    u.sim <- rgamma(length(u), shape=sigma.u, scale=1/sigma.u)
    L[i] <- exp(-nll(u.sim,beta,sigma.u,X))/
      prod(dnorm(u.sim,mean=u,sd=s))
  }
  -log(mean(L))
}

k <- 100000
L <- nll.LA4Rexsamp(fit1$par,X,k=k,seed=16345)
c(L,fit1$objective, fit2$objective)

#Plotting the confidence. This can take a long time
n<-1000
L.sim <- numeric(n)
for(i in 1:n){
  L.sim[i] <- exp(-nll.LA4Rexsamp(fit1$par,X,k=1,seed=i))
}

## For plotting
cs <- cumsum(L.sim)/(1:n)
css <- cumsum(L.sim^2)/1:n
csse <- css/(1:n)-(cs/1:n)^2


#set y_lim for confidence intervals, but very small values
y_lim = c(min(-log((cs+2*sqrt(csse))[-(1:2)])),max(-log((cs-2*sqrt(csse))[-(1:2)])))

plot(log10(1:n),-log(cs),type="l")# ylim=y_lim
lines(log10(3:n),-log((cs+2*sqrt(csse))[-(1:2)]),type="l",col=2)#
lines(log10(3:n),-log((cs-2*sqrt(csse))[-(1:2)]),type="l",col=2)#
lines(log10(c(1,n)),fit1$objective*c(1,1),col=3)

### Problem C







