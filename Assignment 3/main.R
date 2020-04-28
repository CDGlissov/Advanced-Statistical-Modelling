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
setwd("~/Google Drev/UNI/Semester 10/Statmod 2/Advanced-Statistical-Modelling/Assignment 3")

###############################################################################################################################
##### PART 1 ##################################################################################################################
###############################################################################################################################

##### A #######################################################################################################################
concrete <- read.csv("concrets.csv", sep = " ")
concrete$date <- as.Date(concrete$date)
concrete$batch <- as.factor(concrete$batch)

par(mfrow=c(1,1))
plot(y28 ~ date, data = concrete, ylim = c(0,30), col = batch)
points(y7 ~ date, data = concrete, col = batch)

plot(y28 ~ y7, data = concrete, col = batch)

# Estimating the mean concrete strengt within the batches
par(mfrow=c(1,3))
plot.design(concrete)

fit <- lm(y28 ~ batch -1, data = concrete)
summary(fit)

# Date does not seem to show any dependence on the strength, only the batches and that
# information we already have in 'batch'
# Trying different mixed effects models
fit1 <- lme(y28 ~ air.temp, random = ~1 | batch, data = concrete, method = "ML")
fit2 <- lme(y28 ~ y7, random = ~1 | batch, data = concrete, method = "ML")
fit3 <- lme(y28 ~ air.temp + y7, random = ~1 | batch, data = concrete, method = "ML")
fit4 <- lme(y28 ~ air.temp*y7, random = ~1 | batch, data = concrete, method = "ML")
fit5 <- lme(y28 ~ air.temp, random = ~1 | y7, data = concrete, method = "ML")

round(c(summary(fit1)$AIC,summary(fit2)$AIC,summary(fit3)$AIC,summary(fit4)$AIC,summary(fit5)$AIC))
# The AIC is lowest for fit 2,3 and 4 which include y7, y7 + temp and y7 + temp + y7:temp
summary(fit2)
# So is the air temperature significant?
anova.lme(fit2,fit3,type = "marginal")
fitFinal <- lme(y28 ~ y7, random = ~1 | batch, data = concrete)
summary(fitFinal)
# Air temperature does not seem to be significant

##### B #######################################################################################################################
# moment estimates

nBatch <- length(levels(concrete$batch))
ni <- summary(concrete$batch)
N = length(concrete$y28)

mu <- with(concrete,c(mean(y7),mean(y28)))

# Calculating the different errors
mui <- matrix(0,nrow = 2, ncol = nBatch)
SSE <- matrix(0,nrow = 2, ncol = 2)
SSB <- matrix(0,nrow = 2, ncol = 2)

for (i in 1:nBatch) {
  mui[,i] <- cbind(mean(concrete$y7[concrete$batch == i]),mean(concrete$y28[concrete$batch == i]))
  SSE <- SSE + (ni[i]-1)*var(cbind(concrete$y7[concrete$batch == i],concrete$y28[concrete$batch == i]))
  XP <- mui[,i]-(mu)
  SSB <- SSB + ni[i]*XP%*%t(XP)
}
mui
SSE
SSB
SST <- (length(concrete$y28)-1)*var(cbind(concrete$y7,concrete$y28))
SST

# The moment estimates:
n0 <- (N-sum(ni^2)/N)/(nBatch-1)
Sigma.m <- SSE/(N-nBatch)
Sigma0.m <- 1/n0*(SSB/(nBatch-1)-Sigma.m)

cov2cor(Sigma.m)
cov2cor(Sigma0.m)

## Within meter covariance matrix
SigBlock <- kronecker(matrix(1,ncol=3,nrow=3),Sigma0.m) +
  kronecker(diag(3),Sigma.m)
Sig <- kronecker(diag(6),SigBlock)

SigBlock

## Within batch corelation matrix
round(cov2cor(SigBlock),digits=TRUE)

#########################################################################################################################
###  PART 2  ############################################################################################################
#########################################################################################################################

# b #####################################################################################################################
# Reading the data
clo <- read.csv('dat_count3.csv', sep =';')
clo$day <- as.factor(clo$day)
clo$subjId <- as.factor(clo$subjId)
clo$sex <- factor(clo$sex)
library(numDeriv)

# Setting up the design matrix
X <- matrix(0,ncol=6,nrow=dim(clo)[1])
X[,1] = 1
X[,2] = clo$time
X[,3] = clo$nobs
X[,4] = clo$sex
X[,5] = clo$tOut
X[,6] = clo$tInOp

# Joint probability
nll <- function(u,beta,sigma.u,X){
  # One u for every observation, i.e. will have more of the same
  U <- u[clo$subjId]
  eta <- X%*%beta
  -sum(dpois(clo$clo,lambda=exp(eta + log(U)),log=TRUE) +
         dgamma(U, shape=sigma.u, rate=sigma.u,log=TRUE))
}

# Laplacian log likehood
nll.LA <- function(theta,X){
  beta <- theta[1:dim(X)[2]]
  sigma.u <- exp(theta[dim(X)[2]+1])
  est <- nlminb(rep(1,47), objective = nll,
                beta=beta, sigma.u=sigma.u,X=X, lower = rep(0,47),control = list(
                  xf.tol = .Machine$double.eps))
  u <- est$par
  l.u <- est$objective
  H <- hessian(func = nll, x = u, beta = beta, sigma.u = sigma.u,
               X=X)
  l.u + 0.5 * log(det(H/(2*pi)))
}

# optimizing.. Takes 3-4 minutes
system.time((fit1 <- nlminb(c(rep(0,6),1),nll.LA,X=X, control = list(trace = 5))))
fit1
# Getting the resulting random effects
opt <- nlminb(rep(1,47),objective = nll, beta = fit1$par[1:6], sigma.u = fit1$par[7], X = X, lower = rep(0,47))

#### Trying to optimize one dimension at a time, since u's are independent
# Wont work for some reason.. Gives loglikelihood of -228
nll.LA2 <- function(theta,X){
  beta <- theta[1:dim(X)[2]]
  sigma.u <- exp(theta[dim(X)[2]+1])
  fun.tmp <- function(ui,u,beta,sigma.u,X,i){
    ut <- exp(u*0)
    ut[i]<-ui
    nll(ut,beta,sigma.u,X)
  }
  u <- numeric(47)
  ## Use grouping structure
  for(i in 1:length(u)){
    u[i] <- nlminb(1,objective = fun.tmp, u=u,beta=beta,
                   sigma.u=sigma.u,X=X,i=i, lower=0)$par
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

# Checking the accuracy by importance sampling
nll.LA4Rexsamp <- function(theta,X,k,seed){
  set.seed(seed)
  beta <- theta[1:dim(X)[2]]
  sigma.u <- exp(theta[dim(X)[2]+1])
  est <- nlminb(rep(1,47),objective = nll,
                beta=beta, sigma.u=sigma.u,X=X, lower=rep(0,47))
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

# Plotting the confidence
n<-1000
L.sim <- numeric(n)
for(i in 1:n){
  L.sim[i] <- exp(-nll.LA4Rexsamp(fit1$par,X,k=1,seed=i))
  if (mod(i,10)==0) {
    print(i)
  }
}

## For plotting
cs <- cumsum(L.sim)/(1:n)
css <- cumsum(L.sim^2)/1:n
csse <- css/(1:n)-(cs/1:n)^2


#set y_lim for confidence intervals, but very small values
y_lim = c(min(-log((cs+2*sqrt(csse))[-(1:2)])),max(-log((cs-2*sqrt(csse))[-(1:2)])))

plot(log10(1:n),-log(cs),type="l", ylim = y_lim)# ylim=y_lim
lines(log10(3:n),-log((cs+2*sqrt(csse))[-(1:2)]),type="l",col=2)#
lines(log10(3:n),-log((cs-2*sqrt(csse))[-(1:2)]),type="l",col=2)#
lines(log10(c(1,n)),fit1$objective*c(1,1),col=3)

