library(lmtest)
setwd("C:/Users/Christian/Desktop/5. Semester/StatMod/Assignment 2")

logistic <- read.delim("Logistic.txt", header=TRUE, sep="")

#n = total number og ikke no, explanatory er logistic$AZT
fitlog <- glm(cbind(logistic$AIDS_yes, n-logistic$AIDS_yes)~logistic$AZT, family=binomial("logit"), data=logistic)
summary(fitlog)

########odds ratio, it's logit so e^(b1) og e^(b0)###########
#When patient is not treated with AZT, odds ratio is 0.355
exp(fitlog$coefficients[1])
#When patient is treated with AZT, odds ratio is 0.486, this is how much the odds 
#of AIDS decreases
exp(fitlog$coefficients[2])

#confidence interval and odds ratio, interval doesn't cover 0, indicates significance
round(exp(cbind(OR = coef(fitlog), confint.default(fitlog, level=0.95))),5)

######################MODEL COMPARISON########################
#Log ratio test, givet ved deviance, page 168, model comparison, IAL
#null model deviance and the deviance of the model with AZT
LRT<- fitlog$null.deviance-fitlog$deviance
#pvalue
1-pchisq(LRT, df=1)

#this is less than 5%, we reject H0: AZT have no effect on AIDS.
#LRT sees if there is a significant difference between the null model
#and the model with AZTyes

#wald test, works best for large data sets. Tests if B1 is significant
#it does this by looking at if it improves the MLE or not
#H0: B1=0, if the nulhypothesis is accepted, B1 is insignificant
Z<-(fitlog$coefficients[2]/summary(fitlog)$coefficients[, 2][2])^2
1-pchisq(Z, df=1)
#we see it's significant, H0 is rejected

#confidence interval for the coefficients, they don't span over 0. They are significant
confint(fitlog)

#Score test, what's wrong? regn den ud selv
s <- function(psi,x,N){
    theta <- exp(psi)/(1+exp(psi))
    n *(mean(x) - N * theta) /(theta*(1-theta))
}

I <- function(psi, x, N){
  n*N*exp(psi)/(1+exp(psi))^2
}

#score test med ANOVA, den er signifikant
anova(fitlog, test="Rao")



###################### PART 3 SURVIVAL ANALYSIS ##########################
library(survival)
actg<-read.delim("actg320.txt", header=TRUE, sep="")

# Number of events indicate AIDS or death, 1=dead, 63 deaths in control treatment
#33 in new treatment
by(actg$event, actg$tx, sum)

# Total follow up time in days
by(actg$time, actg$tx, sum)
#average of 236 days
by(actg$time, actg$tx, summary)

#survival fit
Surv.tx <- survfit(Surv(time, event == 1)~tx, conf.type = "log-log", data = actg)
Surv.tx

#clearly the new treatment group are doing better, confidence interval doesn't cross each other
#and we see that overall the new treatment group has a higher survival probability
plot(Surv.tx, conf.int = TRUE, las = 1, xlab = "Days since test start",
     ylab = "Estimated Survival Probability", col=1:2, lwd = 1)
legend("bottomleft", col=1:2, c("Control treatment","New treatment"), lwd=2)

#Cumulative 
plot(Surv.tx, fun=function(x) {1- x}, conf.int = TRUE, las = 1,
     xlab = "Days since test start.", ylab = "Estimated Failure Probability", col = 1:2)
legend("topleft", col=1:2, c("Control Treatment", "New treatment"), lwd=2)

#log rank test
survdiff(Surv(time, event == 1) ~ tx, data = actg)
#low p-value, reject that they are the same.
1-pchisq(10.5, df=1)
###########################PARAMETRIC SURVIVAL MODELS#########################

par(mfrow = c(2,2))

# The exponential model
modexp <- survreg(Surv(time, event) ~ tx + cd4, data = actg, dist = "exponential")
summary(modexp)

# The Weibull model
modWeibull <- survreg(Surv(time, event) ~ tx + cd4, data = actg, dist = "weibull")
summary(modWeibull)

# The log-logistic model
modLog <- survreg(Surv(time, event) ~ tx + cd4, data = actg, dist = "loglogistic")
summary(modLog)


par(mfrow=c(1,3))
# Model check (Coxsnell) exponential
actg$CoxSnell <- actg$time*exp(-modexp$linear.predictors)
survexp <- survfit(Surv(CoxSnell, event == 1) ~ 1 , data = actg)
plot(survexp$time, -log(survexp$surv))
abline(a=0, b=1, col = 2)


# Model check (Coxsnell) weibull
actg$CoxSnell2 <- exp((log(actg$time)-modWeibull$linear.predictors)/modWeibull$scale)
survWeibull <- survfit(Surv(CoxSnell2, event==1)~1 , data = actg)
plot(survWeibull$time, -log(survWeibull$surv))
abline(a=0, b=1, col = 2)


#model check (coxsnell) log logistic
actg$CoxSnell3 <- (log(actg$time)-modLog$linear.predictors)/modLog$scale
actg$CS3 <- log(1+exp(actg$CoxSnell3))
survLog <- survfit(Surv(CS3, event==1)~1 , data = actg)
plot(survLog$time, -log(survLog$surv))
abline(a=0, b=1, col = 2)

par(mfrow = c(1,1))
#Looking at the AIC, log-logistic model seems to be the best.
AICexp <- -2*modexp$loglik[2]+ 2*modexp$df; AICexp
AICWeibull <- -2*modWeibull$loglik[2]+ 2*modWeibull$df; AICWeibull
AICLog <- -2*modLog$loglik[2]+ 2*modLog$df; AICLog

confint(modLog)
summary(modLog)

#Time Ratio, it is seen that 1 is not included, this means it is significant
exp(cbind(coef(modLog),confint(modLog)))
exp(cbind(coef(modLog),confint(modLog))*50)

xrange <- range(actg$time); xrange
t <-seq(xrange[1],xrange[2],length=100)

# Log-Logistic
coef4 <- modLog$coefficients
# The model for the new treatment (x = 1)
z411 <- (log(t)-(coef4[1]+coef4[2]+coef4[3]))/modLog$scale
# The model for the placebo treatment (x = 0)
z400 <- (log(t)-(coef4[1]+coef4[3]))/modLog$scale
S411 <- (1+exp(z411))^-1
S400 <- (1+exp(z400))^-1

# Count = 50
# The model for the new treatment w. CD4=80 (x = 1)
z41 <- (log(t)-(coef4[1]+coef4[2]+coef4[3]*50))/modLog$scale
# The model for the placebo treatment w. CD4=50 (x = 0)
z40 <- (log(t)-(coef4[1]+coef4[3]*50))/modLog$scale
S41 <- (1+exp(z41))^-1
S40 <- (1+exp(z40))^-1

# get the range for the y axis
yrange <- range(seq(0,0.1,0.001)); yrange
# set up the plot
plot(xrange, yrange, type="n",xlab= "Days since test start", ylab=
       "Estimated probability of failure",las= 1)
# Log logistic
lines(t, 1-S411, type="l", col=1, lwd=2)
lines(t, 1-S400, type="l", col=2, lwd=2)
lines(t, 1-S41, type="l", col=1, lty = 2, lwd=2)
lines(t, 1-S40, type="l", col=2, lty = 2, lwd=2)
legend(x = "bottomright", lty = c(1,1,2,2), lwd = 2, col = 1:2, legend=
         c("New Treatment","Controlloed Treatment", "New treatment CD4*50", "Controlled treament CD4*50"))

