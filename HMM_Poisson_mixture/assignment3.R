#####PART 1######

setwd("C:/Users/Christian/Desktop/5. Semester/StatMod/Assignment 3")
soap<-read.table("soap.txt", header=F)

#weekly sales of soap product
names(soap)="x"
#we see that the mean and sd is not equal to each other, might indicate
#that the data is not poisson and that underdispersion might take place
summary(soap$x)



## looking into overdispersion
summary(fit1 <- glm(soap$x~1,data=dat,family=poisson))
summary(glm(soap$x~1,data=dat,family=quasipoisson))

summary(fit1)
## Goodness of fit test
1-pchisq(630,df=241)
## Hence we need either overdispersion or other model.

#sort into frequencies
freq <- vector(mode="integer", length=23)
for(i in 1:length(soap$x)){
  freq[soap$x[i]+1]=freq[soap$x[i]+1]+1
}
freq<-as.integer(freq)
n<-0:22
dat = data.frame(n,freq)
dat

#look at plot
plot(dat$n,dat$freq)
hist(dat$freq)

#we try to fit the data
nlp <- function(theta,x){
  -sum(dpois(x,lambda=theta, log=TRUE))
}
opt.nlp<-optimize(nlp, c(0.01,22), x=soap$x)
plot(0:22, dat$freq/sum(dat$freq),type="h")
points(0:22,dpois(0:22,lambda=opt.nlp$minimum),pch=19,
       col=2)
#the fit looks bad.

#lav goodness of fit
#expected values:
e <- dpois(0:22,lambda = opt.nlp$minimum)*sum(dat$freq)
e

obs <- dat$freq

#her ses stor under-dispersion og en smule over-dispersion omkring mean.
e-obs

#finder residuals
r<-(obs-e)^2/e

chisq.stat <- sum(r)
chisq.stat

df<-length(obs)-1
df

1-pchisq(chisq.stat,df)
#meget signifikant og vi afviser nulhypotesen (The data are consistent with the specified distribution), 
#modellen er derfor ikke en poisson

############ PART 2 ##############


############ TRANSFORM ##############
#To make mixture, we first have to make the variables unconstrained to be
#able to maximize, slide 13 uge 10
p2w <- function(m, lambda, delta){
  if(sum(delta) >= 1){return("sum(delta) should be < 1")}
  if(length(lambda) != m){return("length(lambda) should be m")}
  if(length(delta) != (m-1)){return("length(delta) should be m-1")}
  eta<-log(lambda)
  tau <- log(delta/(1-sum(delta)))
  return(list(eta=eta,tau=tau))
}

p2n <- function(m, eta, tau){
  if(m==1){return(exp(eta))}
  if(length(eta) != m){return("length(lambda) should be m")}
  if(length(tau) != (m-1)){return("length(delta) should be m-1")}
  lambda <- exp(eta)
  delta <- exp(tau)/(1+sum(exp(tau)))
  #delta1 = 1-sum(delta)
  delta <- c(1-sum(delta),delta)
  return(list(lambda=lambda,delta=delta))
}

#Let's define the nll:

nll <- function(theta, state, y){
  if(state==1){
    return(-sum(dpois(y,lambda=exp(theta), log = TRUE)))
  }
  #this is our parameters delta and lambda transformed
  eta <- theta[1:state]
  tau <- theta[(state+1):(2*state-1)]
  #make to natural parameters for natural likelihood
  n.pars <- p2n(state,eta,tau)
  n <- length(y)
  #likelihood
  nll <- 0
  for(i in 1:n){
    nll <- nll - log(sum(n.pars$delta * dpois(y[i],lambda=n.pars$lambda)))
  }
  return(nll)
}

y <- soap$x

#fit 1 mixture:
m <- 1; lambda <- mean(y); delta <- c()
#make to working parameters so we can optimize
wpars <- p2w(m,lambda,delta)
theta <- c(wpars$eta,wpars$tau)

#using nlminb to optimize the likelihood and parameters
opt1 <- nlminb(theta,nll,state=m,y=y)
#makes working parameters to natural, which gives us the natural parameters
pars1 <- p2n(m,opt1$par[1:m],opt1$par[(m+1):(2*m-1)])

#let's do it for m=2 mixtures, with C=delta equal probability
#lambda is split for each mean.
m <- 2; lambda <- c(1/2,3/2)*mean(y); delta <- c(0.5)
wpars <- p2w(m,lambda,delta)
theta <- c(wpars$eta,wpars$tau)
opt2 <- nlminb(theta,nll,state=m,y=y)
pars2 <- p2n(m,opt2$par[1:m],opt2$par[(m+1):(2*m-1)])

#for 3:
m <- 3; lambda <- c(1/2,1,3/2)*mean(y); delta <- c(1,1)/3
wpars <- p2w(m,lambda,delta)
theta <- c(wpars$eta,wpars$tau)
opt3 <- nlminb(theta,nll,state=m,y=y)
pars3 <- p2n(m,opt3$par[1:m],opt3$par[(m+1):(2*m-1)])

#for 4
m <- 4; lambda <- c(1/4,3/4,5/4,7/4)*mean(y); delta <- c(1,1,1)/4
wpars <- p2w(m,lambda,delta)
theta <- c(wpars$eta,wpars$tau)
opt4 <- nlminb(theta,nll,state=m,y=y)
pars4 <- p2n(m,opt4$par[1:m],opt4$par[(m+1):(2*m-1)])

#plot it all
plot(0:22, dat$freq/sum(dat$freq),type="h",ylim=c(0,0.2))
lines(0:22, dpois(0:22, lambda = pars1), col=2, lwd=2)

lines(0:22,pars2$delta[1]*dpois(0:22,lambda=pars2$lambda[1]) +
         pars2$delta[2]*dpois(0:22,lambda=pars2$lambda[2]),
       col=3, lwd=2)
lines(0:22,pars3$delta[1]*dpois(0:22,lambda=pars3$lambda[1]) +
         pars3$delta[2]*dpois(0:22,lambda=pars3$lambda[2]) +
         pars3$delta[3]*dpois(0:22,lambda=pars3$lambda[3]),
       col=4, lwd=2)
lines(0:22,pars4$delta[1]*dpois(0:22,lambda=pars4$lambda[1]) +
        pars4$delta[2]*dpois(0:22,lambda=pars4$lambda[2]) +
        pars4$delta[3]*dpois(0:22,lambda=pars4$lambda[3])+
        pars4$delta[4]*dpois(0:22,lambda=pars4$lambda[4]),
      col=5, lwd=2)
legend("topright", c("m=1", "m=2", "m=3", "m=4"), lty=c(1,1,1,1),lwd=c(2,2,2,2), col=2:5)

#choosing a model:
AIC <- 2 * c(opt1$objective, opt2$objective,
              opt3$objective, opt4$objective) +
    2 * c(length(opt1$par), length(opt2$par),
          length(opt3$par), length(opt4$par))
AIC
#we see for m=3 is the best model
1-pchisq(2*(opt1$objective-opt2$objective),
         df=length(opt2$par)-length(opt1$par))

1-pchisq(2*(opt2$objective-opt3$objective),
         df=length(opt3$par)-length(opt2$par))

1-pchisq(2*(opt3$objective-opt4$objective),
         df=length(opt4$par)-length(opt3$par))


#wald test for working parameters
library(numDeriv)

#find hessian
H <- hessian(nll, opt3$par, state=3, y=soap$x)
V.pars <- solve(H)
se <- sqrt(diag(V.pars))
#for the first parameter
exp(opt3$par[1]+qnorm(c(0.025,0.975))*se[1])
exp(opt3$par[2]+qnorm(c(0.025,0.975))*se[2])
exp(opt3$par[3]+qnorm(c(0.025,0.975))*se[3])
exp(opt3$par[4]+qnorm(c(0.025,0.975))*se[4])
exp(opt3$par[5]+qnorm(c(0.025,0.975))*se[5])
###shouldn't we have a sixth working parameter parameter???###
#exp((1-sum(opt3$par[4:5]))+qnorm(c(0.025,0.975))*se[6])

#we choose the profile likelihood for the first parameter
nllp <- function(eta1, m, y){
  f.tmp<-function(tau, eta1, y, m){
    theta <- c(eta1, tau) #vi vil optimere tau, det er vores nuissance
    nll(theta, m, y)
  }
  m <- 2; lambda <- c(1/2,3/2)*mean(y); delta <- c(0.5)
  wpars <- p2w(m,lambda,delta) #transform to working parameters
  theta <- c(wpars$eta[2],wpars$tau) #excluding eta[1]
  nlminb(theta, f.tmp, eta1=eta1, y=y, m=m)$objective
}

eta1 <- seq(0.5, 3, by=0.01 )
profile<-sapply(eta1, nllp, m=2, y=y)
plot(eta1, exp(-profile-max(-profile)))

###REPARAMATIZING so only 1 maxima:

p2w2 <- function(m, lambda, delta){
  if(sum(delta) >= 1){return("sum(delta) should be < 1")}
  if(length(lambda) != m){return("length(lambda) should be m")}
  if(length(delta) != (m-1)){return("length(delta) should be m-1")}
  eta<-log(lambda-c(0,lambda[-m]))
  tau <- log(delta/(1-sum(delta)))
  return(list(eta=eta,tau=tau))
}

p2n2 <- function(m, eta, tau){
  if(m==1){return(exp(eta))}
  if(length(eta) != m){return("length(lambda) should be m")}
  if(length(tau) != (m-1)){return("length(delta) should be m-1")}
  lambda <- cumsum(exp(eta))
  delta <- exp(tau)/(1+sum(exp(tau)))
  #delta1 = 1-sum(delta)
  delta <- c(1-sum(delta),delta)
  return(list(lambda=lambda,delta=delta))
}

nll2 <- function(theta, state, y){
  if(state==1){
    return(-sum(dpois(y,lambda=exp(theta), log = TRUE)))
  }
  #this is our parameters delta and lambda transformed
  eta <- theta[1:state]
  tau <- theta[(state+1):(2*state-1)]
  #make to natural parameters for natural likelihood
  n.pars <- p2n2(state,eta,tau)
  n <- length(y)
  #likelihood
  nll2 <- 0
  for(i in 1:n){
    nll2 <- nll2 - log(sum(n.pars$delta * dpois(y[i],lambda=n.pars$lambda)))
  }
  return(nll2)
}

nllp2 <- function(eta1, m, y){
  f.tmp<-function(tau, eta1, y, m){
    theta <- c(eta1, tau) #vi vil optimere tau, det er vores nuissance
    nll2(theta, m, y)
  }
  m <- 2; lambda <- c(1/2,3/2)*mean(y); delta <- c(0.5)
  wpars <- p2w2(m,lambda,delta) #transform to working parameters
  theta <- c(wpars$eta[2],wpars$tau) #excluding eta[1]
  nlminb(theta, f.tmp, eta1=eta1, y=y, m=m)$objective
}

eta1 <- seq(0, 3, by=0.01 )
profile<-sapply(eta1, nllp2, m=2, y=y)
plot(eta1, exp(-profile-max(-profile)))
a<--(profile-max(profile))
eta1[which(a==max(a))]

###############UNCERTAINTY OF NATURAL PARAMETERS
##SIMULATION:
## Confidence intervals for natural pars
## by simulation from normal dist.
library(mvtnorm)
k <- 100000
PARS <- rmvnorm(k,mean=opt3$par,sigma=V.pars)
dim(PARS)
## Simulated (lambda1)
(CIlambda1 <- quantile(exp(PARS[ ,1]),probs=c(0.025,0.975)))
## Wald based
exp(opt3$par[1]+qnorm(c(0.025,0.975))*se[1])
#we can see that they are close to each other

## Simulated (lambda2)
CIlambda1 <- rbind(CIlambda1,quantile(exp(PARS[ ,2]),probs=c(0.025,0.975)))
CIlambda1[2, ]
## Wald based
exp(opt3$par[2]+qnorm(c(0.025,0.975))*se[2])

## Simulated (lambda3)
CIlambda1 <- rbind(CIlambda1,quantile(exp(PARS[ ,3]),probs=c(0.025,0.975)))
CIlambda1[3, ]
## Wald based
exp(opt3$par[3]+qnorm(c(0.025,0.975))*se[3])

## Simulated delta1-delta3
delta2 <- exp(PARS[ ,4])/(1+rowSums(exp(PARS[ ,4:5])))
delta3 <- exp(PARS[ ,5])/(1+rowSums(exp(PARS[ ,4:5])))
delta1 <- 1-delta2-delta3

## Estimated values of delta 1-3
delta2.hat <- exp(opt3$par[4]) /
  (1+sum(exp(opt3$par[4:5])))
delta3.hat <- exp(opt3$par[5]) /
  (1+sum(exp(opt3$par[4:5])))
delta1.hat <- 1-delta2.hat-delta3.hat

## CI for delta
CIdelta1 <- c(delta1.hat,quantile(delta1,probs=c(0.025,0.5,0.975)))
CIdelta1 <- rbind(CIdelta1,
                  c(delta2.hat,quantile(delta2,
                                        probs=c(0.025,0.5,0.975))))
CIdelta1 <- rbind(CIdelta1,
                  c(delta3.hat,quantile(delta3,
                                        probs=c(0.025,0.5,0.975))))
CIdelta1
############ PART 3 ##############

## book scripts
source("A1.R")
y <- soap$x
length(y)

acf(soap$x)

## 1 - state 
## Initial values
m <- 1
lambda0 <- mean(y)
gamma0 <- 0

## optimize
fit1 <- pois.HMM.mle(y,m,lambda0,gamma0)
fit1

## 2 - state 
## Initial values
m <- 2
lambda0 <- quantile(y,c(0.25,0.75))
gamma0 <- matrix(0.05,ncol=m,nrow=m)
diag(gamma0) <- 1-(m-1)*gamma0[1,1]

## optimize
fit2 <- pois.HMM.mle(y,m,lambda0,gamma0)
fit2

## 3 - state 
## Initial values
m <- 3
lambda0 <- quantile(y,c(0.25,0.5,0.75))
gamma0 <- matrix(0.05,ncol=m,nrow=m)
diag(gamma0) <- 1-(m-1)*gamma0[1,1]

## optimize
fit3 <- pois.HMM.mle(y,m,lambda0,gamma0)
fit3


## 4 - state 
## Initial values
m <- 4
lambda0 <- quantile(y,c(0.2,0.4,0.6,0.8))
gamma0 <- matrix(0.025,ncol=m,nrow=m)
diag(gamma0) <- 1-(m-1)*gamma0[1,1]

## optimize
fit4 <- pois.HMM.mle(y,m,lambda0,gamma0)
fit4

AIC <- c(fit1$AIC,fit2$AIC,fit3$AIC,fit4$AIC)
ll <-  -c(fit1$mllk,fit2$mllk,fit3$mllk,fit4$mllk)
AIC
m <- c(1,2,3,4)
df <- m + (m^2-m) ## lambda + gamma
## What should we report
cbind(df, AIC, ll)

###############Confidence intervals for working Parameters#########
## Finding the se for model 3, (the best model)
m <- 3
parvect  <- pois.HMM.pn2pw(m,fit3$lambda,fit3$gamma)
mod <- nlm(pois.HMM.mllk,parvect,x=y,m=m,
           hessian=TRUE)  
mod

parvect <- mod$estimate
names(parvect) <- c("lambda1","lambda2","lambda3","tau21",
                    "tau31","tau12","tau32","tau13","tau23")


se <- sqrt(diag(solve(mod$hessian)))

## Working pars + standard error
round(cbind(parvect,se),digits=2)
fit3$gamma

###confidence intervals
confwp<-exp(parvect[i]+qnorm(c(0.025,0.975))*se[i])
for(i in 1:length(parvect)){
  confwp<-parvect[i]+qnorm(c(0.025,0.975))*se[i]
  print(confwp)
}

#########confidence interval of natural parameters
## Profile likelihood for lambda 1
PL.lambda1 <- function(lambda1,m,y,lambda0,gamma0){
  fun.tmp <- function(pars,lambda1,y,m){
    parvect <- c(log(lambda1),pars)
    pois.HMM.mllk(parvect,y,m)
  }
  lambda0 <- c(lambda1,lambda0)
  parvect0 <- pois.HMM.pn2pw(m, lambda0, gamma0)
  parvect0 <- parvect0[-1]
  np <- length(parvect0)
  lower    <- rep(-10,np)
  upper    <- c(rep(max(y),m-1),rep(10,np+1-m))
  nlminb(parvect0,fun.tmp, lambda1=lambda1,
         y=y, m=m, lower=lower,
         upper=upper)$objective    
}

## Initial values for estimation
lambda0 <- quantile(y,probs=c(1/3,2/3))
#PL.lambda1(lambda1=29,m=m,y=y,lambda0,gamma0)

## Which lamdas should we look at=
lambda1 <- seq(min(y),max(y),length=100)

## The profile liklielihood 
llp.lambda1 <- sapply(lambda1,PL.lambda1,m=m,y=y,
                      lambda0=lambda0,gamma0=gamma0)


par(mfrow=c(1,1))
plot(lambda1,exp(-(llp.lambda1-fit3$mllk)),
     type="l")
lines(range(lambda1),c(1,1)*exp(-qchisq(0.95,df=1)/2),col=2,lty=2,lwd=2)
rug(fit3$lambda,col=3,lwd=2)

## Wald statistic
cbind(mod$estimate,se)

## Quadratic (local) approximation
cbind(exp(mod$estimate-1.96*se),exp(mod$estimate+1.96*se))[1:3, ]
lines(exp(mod$estimate[1]-1.96*se[1]*c(-1,1)),
      c(1,1)*exp(-qchisq(0.95,df=1)/2),
      col=4,lwd=2)
lines(exp(mod$estimate[2]-1.96*se[2]*c(-1,1)),
      c(1,1)*exp(-qchisq(0.95,df=1)/2),
      col=4,lwd=2)
lines(exp(mod$estimate[3]-1.96*se[3]*c(-1,1)),
      c(1,1)*exp(-qchisq(0.95,df=1)/2),
      col=4,lwd=2)

######
