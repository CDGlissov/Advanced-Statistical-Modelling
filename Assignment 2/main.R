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
#Plots
library(ggplot2, quietly = TRUE)
library(ggpubr); 
library(ggfortify);
library(GGally);
ggplot2::theme_set(ggplot2::theme_grey())
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}
# Load path
path="C:/Users/Christian/Dropbox/10. Semester/Dataanalyse og Statistisk Modellering/Assignments/Assignment 2"
setwd(path)
# Load sources
source("functions.R")
# Load data
dat_count = read.csv("dat_count.csv", header=TRUE, sep=";")
data(ozone)
dat_count[,c(1,2,3,5,6)] <- lapply(dat_count[,c(1,2,3,5,6)], factor)
head(ozone)
head(dat_count)
n=length(ozone$Ozone)

# SUMMARY STATISTICS 1.1 ###################################################################
#plot ozone

plot_fun <- function(data, mapping, pts=list(), smt=list(), ...){
  ggplot(data = data, mapping = mapping, ...) + 
    do.call(geom_point, pts) +
    do.call(geom_smooth, smt)+
    theme(axis.text.x = element_text(angle = 90, hjust = 1))
}

p=ggpairs(ozone,
          upper = list(continuous="cor"),
          diag =list(continuous=wrap("densityDiag")),
          lower=list(continuous =wrap(plot_fun, pts=list(size=0.4, colour="black"), 
                     smt=list(method="loess", se=T, size=0.2, colour="red"))) )
p

qplot(ozone$Ozone,geom="histogram",binwidth = 1.9,xlab = "Ozone",col=I("black"),ylab = "Count",alpha=I(.4))

#remove invtmp and Hgt due to high correlation with other variables
drops <- c("InvTmp", "Hgt")
ozone=ozone[ , !(names(ozone) %in% drops)]

# Autocorrelation plot
autoCorr <- with(acf(ozone$Ozone, lag.max = 50, plot = F),data.frame(lag,acf))
ggplot(autoCorr,aes(x=lag,y=acf)) + 
  geom_hline(aes(yintercept = 0)) +
  geom_segment(mapping = aes(xend = lag, yend = 0)) +
  geom_hline(aes(yintercept = qnorm(0.975)/sqrt(length(ozone$Ozone))),lty = 2, col = 4, lwd = 0.7) +
  labs(x = "Lag", y = "ACF", title = "Autocorrelation of Ozone measurements")


# TASK 1.2-1.3 #############################################################################

LMInit <- lm(Ozone ~ ., data = ozone)
Anova(LMInit, type = 2)
LM1=model.select(LMInit, verbose = T)
summary(LM1)
autoplot(LM1, which=c(1:3,5), nrow=2, ncol=2) + theme(legend.position="none")

# building simple models
par(mfrow=c(1,1))
l = data.frame(boxCox(LM1, lambda = seq(0,1,0.001)))
l_opt=l$x[l$y==max(l$y)]
ggplot(l,aes(x = x,y = y)) + 
  geom_line() + 
  geom_point(aes(x=l_opt,y=max(y)), col = 2) +
  geom_line(aes(x = rep(l_opt,100),y=seq(max(y),-1460,length.out = 100)), lty = 2) +
  geom_abline(slope = 0,intercept = max(l$y)-qchisq(.95,1)/2, lty = 2, col = 2) +
  xlim(0,0.85) +
  annotate("text", x=0.6,y=-1385,label = "95 %", col = 2) + 
  labs(title = expression(paste("Profile log-likelihood of ",lambda)),x = expression(lambda),y="log-likelihood")
# Seems 1/3 i.e the cubic-root transformation

# Making a new data set with the transformed response
ozoneT <- ozone
ozoneT$Ozone <- (ozone$Ozone^(l_opt)-1)/l_opt

gaus_bc <- lm(Ozone ~ Temp+InvHt+Hum, data = ozoneT)
drop1(gaus_bc, test="Chisq")
summary(gaus_bc)
Anova(gaus_bc, type = 2)
autoplot(gaus_bc, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")
par(mfrow=c(2,1))

acf1=acf(gaus_bc$residuals, plot=FALSE, lag.max=50)
pacf1=pacf(gaus_bc$residuals, plot=FALSE, lag.max=50)
plot(acf1, main="ACF of residuals")
plot(pacf1, main="PACF of residuals")
#sd(ozone$Ozone)/mean(ozone$Ozone)
# GENERALIZED MODEL 1.4-1.5 ################################################################

#inverse/log/identiy, identity doesn't work
gam_inv_init <- glm(Ozone ~ . , family = Gamma(link="inverse"), data = ozone) #identity
gam_log_init <- glm(Ozone ~ . , family = Gamma(link="log"), data = ozone)

#test for fit
1-pchisq(summary(gam_log_init)$deviance, summary(gam_log_init)$df.residual)
1-pchisq(summary(gam_inv_init)$deviance, summary(gam_inv_init)$df.residual)

gam_inv=model.select(gam_inv_init)
gam_log=model.select(gam_log_init)
summary(gam_inv)
summary(gam_log)
anova(gam_log,test = "F")
autoplot(gam_inv, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")
autoplot(gam_log, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

#1/mu^2, inverse, identity and log. 1/mu^2 and identity doesn't work.
igaus_inv_init <- glm(Ozone ~ ., family = inverse.gaussian(link="inverse"), data = ozone)
igaus_log_init <- glm(Ozone ~ ., family = inverse.gaussian(link="log"), data = ozone)

#test for fit
1-pchisq(summary(igaus_inv_init)$deviance, summary(igaus_inv_init)$df.residual)
1-pchisq(summary(igaus_log_init)$deviance, summary(igaus_log_init)$df.residual)

igaus_inv=model.select(igaus_inv_init)
anova(igaus_inv_init, test="F") #type 2/3 testing doesn't work? So use type 1
#igaus_inv=update(igaus_inv_init, ~.-InvTmp)    - Cannot remove since dropped
#anova(igaus_inv, test="F")
igaus_inv=update(igaus_inv_init, ~.-Wind)
anova(igaus_inv, test="F")
#igaus_inv=update(igaus_inv, ~.-Hgt)            - Cannot remove since dropped
#anova(igaus_inv, test="F")

igaus_log=model.select(igaus_log_init)
summary(igaus_log)
summary(igaus_inv)
autoplot(igaus_log, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")
autoplot(igaus_inv, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

AIC(igaus_log)
AIC(igaus_inv)
AIC(gam_log)
AIC(gam_inv)
AIC(LM1)
AIC_transformed(gaus_bc, l_opt, ozone$Ozone)

# Q1.6 weighted matrix as a function of mu #################################################

#FOR GAUSSIAN
X = model.matrix(gaus_bc)
fitted_values=as.vector((X %*% coef(gaus_bc)))
deviance_residuals=sum(((ozone$Ozone^l_opt-1)/l_opt-fitted_values)^2)
df = length(ozone$Ozone)-(length(gaus_bc$coefficients))
dispersion=deviance_residuals/df
scaled_cov=solve(t(X) %*% X)*dispersion
scaled_cov
(summary(gaus_bc)$cov.scaled)

#FOR GAMMA
X = model.matrix(gam_log)
fitted_values=as.vector(exp(X %*% coef(gam_log)))
y=ozone$Ozone
deviance_residuals=sum(2*(y/fitted_values - log(y/fitted_values)-1))
df = length(ozone$Ozone)-(length(gam_log$coefficients))
chi_2=sum((y-fitted_values)^2/fitted_values^2)
dispersion=chi_2/df
scaled_cov=solve(t(X) %*% X) * dispersion
scaled_cov
(summary(gam_log)$cov.scaled)

# Q2.1 and 2.2 develop, present the model ##################################################

# Choosing the transformed gaussian as seen from the residuals a second order term seems relevant
# Also looking at the GAM indicates higher order terms
pgam=gam(Ozone ~ s(Temp)+s(InvHt)+s(Pres)+s(Vis)+s(Hum)+s(Wind),data=ozone)
par(mfrow=c(2,4))
plot(pgam)
par(mfrow=c(1,1))

degree2=paste0("I(",names(ozoneT)[2:7],"^2)",collapse=" + ")
form=as.formula(paste0("Ozone ~ . * . +", degree2))
# Estimating the lambda parameter again for the full model
glm_init <- lm(form, data = ozone)
l <- data.frame(boxCox(glm_init, lambda = seq(0.07,0.5,0.001)))
l_final <- l$x[l$y==max(l$y)]
# Plotting the log-likelihood again
ggplot(l,aes(x = x,y = y)) + 
  geom_line() + 
  geom_point(aes(x=l_final,y=max(y)), col = 2) +
  geom_line(aes(x = rep(l_final,100),y=seq(max(y),-1338,length.out = 100)), lty = 2) +
  geom_abline(slope = 0,intercept = max(l$y)-qchisq(.95,1)/2, lty = 2, col = 2) +
  annotate("text", x=0.45,y=-1332.2,label = "95 %", col = 2) + 
  labs(title = expression(paste("Profile log-likelihood of ",lambda)),x = expression(lambda),y="log-likelihood")

ozoneT$Ozone <- (ozone$Ozone^l_final-1)/l_final
glm_initBC <- lm(form, data = ozoneT)
# The final model
model_final=model.select(glm_initBC)
model_final
Anova(model_final,type = 2)
BIC(model_final)
AIC_transformed(model_final, l_final, ozone$Ozone)
BIC_transformed(model_final, l_final, ozone$Ozone[2:n])
par(mfrow=c(2,1))
acf(model_final$residuals, main = "ACF of residuals")
pacf(model_final$residuals, main = "PACF of residuals")
par(mfrow=c(1,1))
autoplot(model_final, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

Res_glm1<-residuals(model_final)
p_glm = list()
p_glm[[1]]<-ggplot(ozone,aes(x =Temp,y=Res_glm1))+geom_point()+labs(y="Resisduals")
p_glm[[2]]<-ggplot(ozone,aes(x =InvHt,y=Res_glm1))+geom_point()+labs(y="Residuals")
p_glm[[3]]<-ggplot(ozone,aes(x =Pres,y=Res_glm1))+geom_point()+labs(y="Residuals")
p_glm[[4]]<-ggplot(ozone,aes(x =Hum,y=Res_glm1))+geom_point()+labs(y="Residuals")
ggarrange(plotlist=p_glm,ncol=1,nrow=4,labels=c("A","B","C","D"))

###### MODEL WITH AR(1) ######
degree2=paste0("I(",names(ozone)[2:9],"^2)",collapse=" + ")
form=as.formula(paste0("((ozone$Ozone^(l_opt)-1)/l_opt)[2:n] ~ ((ozone$Ozone^(l_opt)-1)/l_opt)[1:(n-1)]+. * .+", degree2))
lag0=ozone[2:n,]
glm_lag_init = glm(form, family=gaussian(link="identity"), data = lag0)
model_lag_final=model.select(glm_lag_init)
BIC_transformed(model_lag_final, l_opt, ozone$Ozone[2:n])
AIC_transformed(model_lag_final, l_opt, ozone$Ozone[2:n])
par(mfrow=c(2,1))
acf(model_lag_final$residuals)
pacf(model_lag_final$residuals)
par(mfrow=c(1,1))
autoplot(model_lag_final, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

# PRED INTERVALS ###########################################################################
glm_final = model_final
n_new = 70
colors = c("Fit" = "black","Confidence interval"="red","Prediction interval" = "blue")

# WIND ##########
wind = with(ozone,seq(min(Wind),max(Wind),length.out=n_new))
newdat = with(ozone,data.frame(Temp = mean(Temp), InvHt = mean(InvHt), Pres = mean(Pres),
                               Vis = mean(Vis), Hum = mean(Hum), Wind = wind))
preddataP = data.frame(predict(glm_final, newdata = newdat, interval = "p"))
preddataC = predict(glm_final, newdata = newdat, interval = "c")

ggplot(ozone, aes(x = Wind, y = Ozone)) + geom_point() +
  geom_line(data=preddataP, aes(x=wind, y = boxcox_inv(lwr,l_final), colour = "Prediction interval"),lwd = 0.5, lty = 2) +
  geom_line(data=preddataP, aes(x=wind, y = boxcox_inv(upr,l_final), colour = "Prediction interval"),lwd = 0.5, lty = 2) +
  geom_ribbon(data=preddataC,aes(x=wind,ymin=boxcox_inv(lwr,l_final),
                     ymax=boxcox_inv(upr,l_final), fill = "Confidence interval"),lwd = 0.3,alpha = 0.5,inherit.aes=FALSE) +
  geom_line(data=preddataP, aes(x=wind, y = boxcox_inv(fit,l_final),colour = "Fit"),lwd = 1, lty = 1) +
  labs(title = "Predictions given Wind", fill = "Legend") +
  scale_color_manual("",breaks = c("Fit","Prediction interval","Confidence interval"),values = colors)

  
# TEMP ##########
temp = with(ozone,seq(min(Temp),max(Temp),length.out=n_new))
newdat1 = data.frame(Temp = temp, InvHt = mean(ozone$InvHt), Pres = mean(ozone$Pres),
                     Vis = mean(ozone$Vis), Hum = mean(ozone$Hum), Wind= mean(ozone$Wind))
preddata1P = predict(glm_final, newdata = newdat1, interval = "p")
preddata1C = predict(glm_final, newdata = newdat1, interval = "c")

ggplot(ozone, aes(Temp, Ozone)) + geom_point() +
  geom_line(data=preddata1P, aes(x=temp, y = boxcox_inv(lwr,l_final), colour = "Prediction interval"),lwd = 0.5, lty = 2) +
  geom_line(data=preddata1P, aes(x=temp, y = boxcox_inv(upr,l_final), colour = "Prediction interval"),lwd = 0.5, lty = 2) +
  geom_ribbon(data=preddata1C,aes(x=temp,ymin=boxcox_inv(lwr,l_final),
                                 ymax=boxcox_inv(upr,l_final), fill = "Confidence interval"),lwd = 0.3,alpha = 0.5,inherit.aes=FALSE) +
  geom_line(data=preddata1P, aes(x=temp, y = boxcox_inv(fit,l_final),colour = "Fit"),lwd = 1, lty = 1) +
  labs(title = "Predictions given Temperature", fill = "Legend") +
  scale_color_manual("",breaks = c("Fit","Prediction interval","Confidence interval"),values = colors)

# InvHt ##########
invht= with(ozone,seq(min(InvHt),max(InvHt),length.out=n_new))
newdat2 = data.frame(Temp = mean(ozone$Temp), InvHt = invht, 
                     Pres = mean(ozone$Pres),
                     Vis = mean(ozone$Vis), Hum = mean(ozone$Hum), 
                     Wind = mean(ozone$Wind))
preddata2P = predict(glm_final, newdata = newdat2, interval = "p")
preddata2C = predict(glm_final, newdata = newdat2, interval = "c")

ggplot(ozone, aes(InvHt, Ozone)) + geom_point() +
  geom_line(data=preddata2P, aes(x=invht, y = boxcox_inv(lwr,l_final), colour = "Prediction interval"),lwd = 0.5, lty = 2) +
  geom_line(data=preddata2P, aes(x=invht, y = boxcox_inv(upr,l_final), colour = "Prediction interval"),lwd = 0.5, lty = 2) +
  geom_ribbon(data=preddata2C,aes(x=invht,ymin=boxcox_inv(lwr,l_final),
                                 ymax=boxcox_inv(upr,l_final), fill = "Confidence interval"),lwd = 0.3,alpha = 0.5,inherit.aes=FALSE) +
  geom_line(data=preddata2P, aes(x=invht, y = boxcox_inv(fit,l_final),colour = "Fit"),lwd = 1, lty = 1) +
  labs(title = "Predictions given Inversion Base Height", fill = "Legend") +
  scale_color_manual("",breaks = c("Fit","Prediction interval","Confidence interval"),values = colors)

# PRES ##########
pres = with(ozone,seq(min(Pres),max(Pres),length.out=n_new))
newdat3 = with(ozone,data.frame(Temp = mean(Temp), InvHt = mean(InvHt), Pres = pres ,
                                Vis = mean(Vis), Hum = mean(Hum), Wind =mean(Wind)))

preddata3P = predict(glm_final, newdata = newdat3, interval = "p")
preddata3C = predict(glm_final, newdata = newdat3, interval = "c")

ggplot(ozone, aes(x = Pres, y=Ozone)) + geom_point() + 
  geom_line(data=preddata3P, aes(x=pres, y = boxcox_inv(lwr,l_final), colour = "Prediction interval"),lwd = 0.5, lty = 2) +
  geom_line(data=preddata3P, aes(x=pres, y = boxcox_inv(upr,l_final), colour = "Prediction interval"),lwd = 0.5, lty = 2) +
  geom_ribbon(data=preddata3C,aes(x=pres,ymin=boxcox_inv(lwr,l_final),
                                 ymax=boxcox_inv(upr,l_final), fill = "Confidence interval"),lwd = 0.3,alpha = 0.5,inherit.aes=FALSE) +
  geom_line(data=preddata3P, aes(x=pres, y = boxcox_inv(fit,l_final),colour = "Fit"),lwd = 1, lty = 1) +
  labs(title = "Predictions given Pressure", fill = "Legend") +
  scale_color_manual("",breaks = c("Fit","Prediction interval","Confidence interval"),values = colors)

# HUM ##########
hum = with(ozone,seq(min(Hum),max(Hum),length.out=n_new))
newdat4 = with(ozone,data.frame(Temp = mean(Temp), InvHt = mean(InvHt), Pres = mean(Pres) ,
                                Vis = mean(Vis), Hum = hum, Wind =mean(Wind)))

preddata4P = predict(glm_final, newdata = newdat4, interval = "p")
preddata4C = predict(glm_final, newdata = newdat4, interval = "c")

ggplot(ozone, aes(x = Hum, y=Ozone)) + geom_point() + 
  geom_line(data=preddata4P, aes(x=hum, y = boxcox_inv(lwr,l_final), colour = "Prediction interval"),lwd = 0.5, lty = 2) +
  geom_line(data=preddata4P, aes(x=hum, y = boxcox_inv(upr,l_final), colour = "Prediction interval"),lwd = 0.5, lty = 2) +
  geom_ribbon(data=preddata4C,aes(x=hum,ymin=boxcox_inv(lwr,l_final),
                                 ymax=boxcox_inv(upr,l_final), fill = "Confidence interval"),lwd = 0.3,alpha = 0.5,inherit.aes=FALSE) +
  geom_line(data=preddata4P, aes(x=hum, y = boxcox_inv(fit,l_final),colour = "Fit"),lwd = 1, lty = 1) +
  labs(title = "Predictions given Humidity", fill = "Legend") +
  scale_color_manual("",breaks = c("Fit","Prediction interval","Confidence interval"),values = colors)

# VIS ##########
vis = with(ozone,seq(min(Vis),max(Vis),length.out=n_new))
newdat5 = with(ozone,data.frame(Temp = mean(Temp), InvHt = mean(InvHt), Pres = mean(Pres) ,
                                Vis = vis, Hum = mean(Hum), Wind =mean(Wind)))

preddata5P = predict(glm_final, newdata = newdat5, interval = "p")
preddata5C = predict(glm_final, newdata = newdat5, interval = "c")

ggplot(ozone, aes(x = Vis, y=Ozone)) + geom_point() + 
  geom_line(data=preddata5P, aes(x=vis, y = boxcox_inv(lwr,l_final), colour = "Prediction interval"),lwd = 0.5, lty = 2) +
  geom_line(data=preddata5P, aes(x=vis, y = boxcox_inv(upr,l_final), colour = "Prediction interval"),lwd = 0.5, lty = 2) +
  geom_ribbon(data=preddata5C,aes(x=vis,ymin=boxcox_inv(lwr,l_final),
                                 ymax=boxcox_inv(upr,l_final), fill = "Confidence interval"),lwd = 0.3,alpha = 0.5,inherit.aes=FALSE) +
  geom_line(data=preddata5P, aes(x=vis, y = boxcox_inv(fit,l_final),colour = "Fit"),lwd = 1, lty = 1) +
  labs(title = "Predictions given Visibility", fill = "Legend") +
  scale_color_manual("",breaks = c("Fit","Prediction interval","Confidence interval"),values = colors)


# QUESTIONS AND NOTES ######################################################################
# Data is not in metric units. Should we convert units to metric, for inference?

# Ozone is positive, not a many values are 0, inverse gaussian or gamma distribution?
# https://cran.r-project.org/web/packages/GlmSimulatoR/vignettes/dealing_with_right_skewed_data.html

# Should look into correlation between some of the variables.

# if ggfortify or gpubr doesnt work, do remotes::update_packages("rlang")

# 



