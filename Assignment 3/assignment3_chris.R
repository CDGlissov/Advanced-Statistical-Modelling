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
library(DHARMa)
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
#clo$subjId = as.factor(clo$subjId)
clo$day=as.factor(clo$day)
clo$sex=as.factor(clo$sex)
summary(clo)
str(clo)

################## Binomial#################
#let's try with a simple model first
fit.glm.bin1=glmmTMB(cbind(clo,nobs-clo)~time+day+sex+tOut+tInOp+(1|subjId), 
                    family=binomial(link = "logit"), data=clo, REML=FALSE)

# Reducing
Anova(fit.glm.bin1, type="II")
fit.glm.bin2<-update(fit.glm.bin1,~.-tOut)
Anova(fit.glm.bin2, type="II")
fit.glm.bin3<-update(fit.glm.bin2,~.-time)
Anova(fit.glm.bin3, type="II")
fit.glm.bin4<-update(fit.glm.bin3,~.-tInOp)
Anova(fit.glm.bin4, type="II")
fit.glm.bin5<-update(fit.glm.bin4,~.-day)
Anova(fit.glm.bin5, type="II")


#FINAL MODEL
fit.glm.bin.ml=glmmTMB(cbind(clo,nobs-clo)~sex+(1|subjId), 
                       family=binomial(link = "logit"), data=clo, REML=FALSE)
fit.glm.bin=glmmTMB(cbind(clo,nobs-clo)~sex+(1|subjId), 
                     family=binomial(link = "logit"), data=clo, REML=TRUE)



## POISSON, https://stats.stackexchange.com/questions/11182/when-to-use-an-offset-in-a-poisson-regression


fit.glm.pois1<-glmmTMB(clo~day+sex+tOut+tInOp+(1|subjId),offset = log(time),family=poisson(link="log"),
                      data=clo, REML=FALSE)


Anova(fit.glm.pois1, type="II")
fit.glm.pois2<-update(fit.glm.pois1,~.-nobs)
Anova(fit.glm.pois2, type="II")
fit.glm.pois3<-update(fit.glm.pois2,~.-tOut)
Anova(fit.glm.pois3, type="II")
fit.glm.pois4<-update(fit.glm.pois3,~.-tInOp)
Anova(fit.glm.pois4, type="II")
fit.glm.pois5<-update(fit.glm.pois4,~.-day)

summary(fit.glm.pois5)

fit.glm.pois.ml<-glmmTMB(clo~sex+(1|subjId), offset=log(time), family=poisson(link="log"),data=clo, REML=FALSE)
fit.glm.pois<-glmmTMB(clo~sex+(1|subjId), offset=log(time), family=poisson(link="log"),data=clo, REML=TRUE)

###COMPARE
summary(fit.glm.pois.ml)
summary(fit.glm.bin.ml)

#reml estimates
c(exp(-2.1378), exp(-1.0834))
c(exp(-3.588)/(1+exp(-3.588)), exp(-1.235)/(1+exp(-1.235)))

#goodness of fit
(1-pchisq(259.5,133)) #poisson
(1-pchisq(257.4,133)) #binomial


##### DIAGNOSTICS
#random effect
ran.ef.pois = exp(as.vector(unlist(ranef(fit.glm.pois))))
data1 = data.frame(x=1:47, y=ran.ef.pois)
ran.ef.bin = exp(as.vector(unlist(ranef(fit.glm.bin))))
data2 = data.frame(x=1:47, y=ran.ef.bin)


#plot random effects
plot1=ggplot(data1, aes(x=factor(x), y=ran.ef.pois)) + geom_point() + xlab("subject ID") +
  ylab("Random effect") +ylim(c(0,4.2))+ geom_linerange(aes(x=x, ymax=ran.ef.pois, ymin=0))+
  geom_text(aes(label=round(ran.ef.pois,2)),hjust=-0.2, vjust=0.3,angle = 90)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))

plot2=ggplot(data2, aes(x=factor(x), y=ran.ef.bin)) + geom_point() + xlab("subject ID") +
  ylab("Random effect") +ylim(c(0,4.2))+ geom_linerange(aes(x=x, ymax=ran.ef.bin, ymin=0))+
  geom_text(aes(label=round(ran.ef.bin,2)),hjust=-0.2, vjust=0.3,angle = 90)+
  theme(axis.text.x = element_text(angle = 90, vjust = 0.5))
grid.arrange(plot1, plot2, ncol=2)

#qqplots
par(mfrow=c(2,2))
qqnorm(ran.ef.pois, col = "blue",main="Poisson: Random Effects")
qqline(ran.ef.pois, col= "blue")
qqnorm(ran.ef.bin, col = "blue",main="Binomial: Random Effects")
qqline(ran.ef.bin, col= "blue")
qqnorm(resid(fit.glm.pois), col =gg_color_hue(2), main="Poisson: Fixed Effects")
qqline(resid(fit.glm.pois), col="blue")
qqnorm(resid(fit.glm.bin), col =gg_color_hue(2), main="Binomial: Fixed Effects")
qqline(resid(fit.glm.bin), col="blue")

#residuals
par(mfrow=c(1,2))
plot(fitted(fit.glm.pois),residuals(fit.glm.pois, "pearson") , main="Poisson: Residuals against fitted",
     col=gg_color_hue(2), ylab="Pearson Residual", xlab="Fitted")
abline(0,0)
plot(fitted(fit.glm.bin),residuals(fit.glm.bin, "pearson") , main="Binomial: Residuals against fitted",
     col=gg_color_hue(2), ylab="Pearson Residual", xlab="Fitted")
abline(0,0)

##### PREDICTING POPULATION #POISSON
maxT<-max(clo$time)
minT<-min(clo$time)
newdata<-data.frame(time=rep(seq(minT,maxT,length.out=68),2), subjId=NA)
nPois<-length(newdata$time)
#half male and half female
newdata$sex<-c(rep("female",nPois/2),rep("male",nPois/2))
predPois <- predict(fit.glm.pois,newdata,se.fit=TRUE)

newdata$Pred<-exp(predPois$fit)
newdata$UpperConf<-exp(predPois$fit + 1.96 *predPois$se.fit)
newdata$LowerConf<-exp(predPois$fit - 1.96 *predPois$se.fit)

# Plotting
p1<-ggplot(newdata)+geom_ribbon(aes(x = time,ymin = LowerConf,ymax=UpperConf,color = sex),alpha = 0.2)+
  geom_line(aes(x= time,y=Pred,color = sex))+ geom_point(data=clo,aes(x=time,y=clo,colour= sex))
p1

##### PREDICTING POPULATION #BINOMIAL
maxT<-max(clo$time)
minT<-min(clo$time)
newdata<-data.frame(time=rep(seq(minT,maxT,length.out=68),2), subjId=NA)
nPois<-length(newdata$time)
#half male and half female
newdata$sex<-c(rep("female",nPois/2),rep("male",nPois/2))
predPois <- predict(fit.glm.bin,newdata,se.fit=TRUE)

invlogit = function(x){
  return(exp(x)/(1+exp(x)))
}

newdata$Pred<-invlogit(predPois$fit)
newdata$UpperConf<-invlogit(predPois$fit + 1.96 *predPois$se.fit)
newdata$LowerConf<-invlogit(predPois$fit - 1.96 *predPois$se.fit)

# Plotting
p2<-ggplot(newdata)+geom_ribbon(aes(x = time,ymin = LowerConf,ymax=UpperConf,color = sex),alpha = 0.2)+
  geom_line(aes(x= time,y=Pred,color = sex))+ geom_point(data=clo,aes(x=time,y=clo,colour= sex))
p2

##### PREDICTING BASED ON SUBJECTS NEED TO MAKE FOR MORE SUBJECTS AND BINOMIAL/POISSON
maxT<-max(clo$time)
minT<-min(clo$time)

newdata<-data.frame(time=rep(seq(minT,maxT,length.out=68),2), subjId=5)
nPois<-length(newdata$time)
newdata$sex<-c(rep("female",nPois/2),rep("male",nPois/2))
predPois <- predict(fit.glm.pois,newdata,se.fit=TRUE, allow.new.levels=FALSE)
newdata$Pred<-exp(predPois$fit)
newdata$UpperConf<-exp(predPois$fit + 1.96 *predPois$se.fit)
newdata$LowerConf<-exp(predPois$fit - 1.96 *predPois$se.fit)

# Plotting
p<-ggplot(newdata)+geom_ribbon(aes(x = time,ymin = LowerConf,ymax=UpperConf,color = sex),alpha = 0.2)+
  geom_line(aes(x= time,y=Pred,color = sex))+ geom_point(data=clo,aes(x=time,y=clo,colour= sex))
p



###### Part B #######
library(numDeriv)

X <- matrix(0,ncol=2,nrow=dim(clo)[1])
X[,1] = 1
X[,2][clo$sex=="male"] = 1 

# Joint probability
nll <- function(u,beta,sigma.u,X){
  # One u for every observation, i.e. will have more of the same
  U <- u[clo$subjId]
  eta <- X%*%beta+log(clo$time)
  -sum(dpois(clo$clo,lambda=exp((eta + log(U))),log=TRUE)) -
         sum(dgamma(u, shape=sigma.u, scale=1/sigma.u,log=TRUE))
}

nll.LA2 <- function(theta,X){
  beta <- theta[1:dim(X)[2]]
  sigma.u <- exp(theta[dim(X)[2]+1])
  fun.tmp <- function(ui,u,beta,sigma.u,X,i){
    u <- exp(u*0)
    u[i]<-ui
    nll(u,beta,sigma.u,X)
  }
  u <- rep(1,47)
  ## Use grouping structure
  for(i in 1:length(u)){
    u[i] <- nlminb(1,objective = fun.tmp, u=u,beta=beta,
                   sigma.u=sigma.u,X=X,i=i, lower=rep(0.01,47))$par
  }
  l.u <- nll(u,beta,sigma.u,X)
  H <- numeric(length(u))
  for(i in 1:length(u)){
    H[i] <- hessian(func = fun.tmp, x = u[i], u=u, beta = beta, 
                    sigma.u = sigma.u, X=X,i=i)}
  l.u + 0.5 * log(prod(H/(2*pi)))
}

system.time((fit1 <- nlminb(c(-2,-1,1),nll.LA2,X=X)))
fit1

opt <- nlminb(rep(1,47),objective = nll, beta = fit1$par[1:2], sigma.u = fit1$par[3], 
              X = X, lower = rep(0,47))
hist(opt$par)

## Check accuracy by importance sampling
nll.LA.ImportWeight <- function(theta,X,k,seed){
  set.seed(seed)
  beta <- theta[1:dim(X)[2]]
  sigma.u <- exp(theta[dim(X)[2]+1])
  est <- nlminb(rep(1,47),objective = nll,
                beta=beta, sigma.u=sigma.u,X=X, lower=rep(0.001,47))
  u <- est$par
  l.u <- est$objective
  H <- diag(hessian(func = nll, x = u, beta = beta, sigma.u = sigma.u,
                    X=X))
  L <- numeric(k)
  s <- sqrt(1/H)
  for(i in 1:k){## Simulation part
    u.sim <- rgamma(length(u), shape=sigma.u, scale=1/sigma.u)
    L[i] <- exp(-nll(u.sim,beta,sigma.u,X))/
      prod(dgamma(u.sim, shape=sigma.u, scale=1/sigma.u))
  }
  -log(mean(L))
}

k <- 1000
L <- nll.LA.ImportWeight(fit1$par,X,k=k,seed=16345)
c(-L,logLik(fit.glm.pois), logLik(fit.glm.bin),-fit1$objective)





#pois_fitted_value = predict(fit.glm.pois,newdata = clo[order(clo$time),], type="response", se.fit=TRUE)
#plot(clo$time, clo$clo)
#points(clo$time[order(clo$time)],pois_fitted_value$fit, col=gg_color_hue(2))

