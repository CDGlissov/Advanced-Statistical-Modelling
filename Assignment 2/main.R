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



# TASK 1.2-1.3 #############################################################################

LMInit <- lm(Ozone ~ ., data = ozone)
LM1=model.select(LMInit)
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
gaus_bc <- glm((Ozone^(l_opt)-1)/l_opt ~ Temp+InvHt+Hum, family=gaussian, data = ozone)
drop1(gaus_bc, test="Chisq")
summary(gaus_bc)
autoplot(gaus_bc, which=c(1:3,5), nrow=2,ncol=2)+ theme(legend.position="none")

# GENERALIZED MODEL 1.4-1.5 ################################################################

#inverse/log/identiy, identity doesn't work
gam_inv_init <- glm(Ozone ~ ., family = Gamma(link="inverse"), data = ozone) #identity
gam_log_init <- glm(Ozone ~ ., family = Gamma(link="log"), data = ozone)
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

anova(igaus_inv_init, test="F") #type 2/3 testing doesn't work? So use type 1
igaus_inv=update(igaus_inv_init, ~.-InvTmp)
anova(igaus_inv, test="F")
igaus_inv=update(igaus_inv, ~.-Wind)
anova(igaus_inv, test="F")
igaus_inv=update(igaus_inv, ~.-Hgt)
anova(igaus_inv, test="F")

igaus_log=model.select(igaus_log_init)
summary(igaus_log)
summary(igaus_inv)
autoplot(igaus_log, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")
autoplot(igaus_inv, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

AIC(igaus_log)
AIC(igaus_inv)
AIC(gam_log)
AIC(gam_inv)
AIC_transformed(gaus_bc, l_opt, ozone$Ozone)

# Q1.6 weighted matrix as a function of mu #################################################

#FOR GAUSSIAN
fitted_values=as.vector((X %*% coef(gaus_bc)))
deviance_residuals=sum((ozone$Ozone^l_opt-fitted_values)^2)
df = length(ozone$Ozone)-(length(gaus_bc$coefficients))
dispersion=deviance_residuals/df
X = model.matrix(gaus_bc)
scaled_cov=solve(t(X) %*% X)*dispersion
scaled_cov
(summary(gaus_bc)$cov.scaled)

#FOR GAMMA
X <- cbind(rep(1,330),as.matrix(gam_log$model)[,c(2,3,4)])
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
pgam=gam(Ozone ~ s(Temp)+s(InvHt)+s(Pres)+s(Vis)+s(Hgt)+s(Hum)+s(InvTmp)+s(Wind),data=ozone)
par(mfrow=c(2,4))
plot(pgam)
par(mfrow=c(1,1))

glm_init <- glm(Ozone ~ .*.*I(Pres^2)+, family=gaussian(link="identity"), data = ozone)

glm1 = model.select(glm_init)
l = boxCox(glm1, lambda = seq(0,1,0.01))
l_final=l$x[l$y==max(l$y)]

glm1_bc = glm(Ozone^l_final ~ Temp + InvHt + Pres + Hum + InvTmp + I(Pres^2) + Temp:InvTmp + 
                Hum:InvTmp, family=gaussian(link="identity"), data = ozone)
drop1(glm1_bc, test="F")
glm1_bc = update(glm1_bc, ~.-Temp:InvTmp)
summary(glm1_bc)
AIC_transformed(glm1_bc, l_final, ozone$Ozone)

#residuals look better but aic is worse
autoplot(gaus_bc, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

Res_glm1<-residuals(glm1_bc)
p_glm = list()
p_glm[[1]]<-ggplot(ozone,aes(x =Temp,y=Res_glm1))+geom_point()+labs(y="Resisduals")
p_glm[[2]]<-ggplot(ozone,aes(x =InvHt,y=Res_glm1))+geom_point()+labs(y="Residuals")
p_glm[[3]]<-ggplot(ozone,aes(x =Pres,y=Res_glm1))+geom_point()+labs(y="Residuals")
p_glm[[4]]<-ggplot(ozone,aes(x =Hum,y=Res_glm1))+geom_point()+labs(y="Residuals")
ggarrange(plotlist=p_glm,ncol=1,nrow=4,labels=c("A","B","C","D"))




# QUESTIONS AND NOTES ######################################################################
# Data is not in metric units. Should we convert units to metric, for inference?

# Ozone is positive, not a many values are 0, inverse gaussian or gamma distribution?
# https://cran.r-project.org/web/packages/GlmSimulatoR/vignettes/dealing_with_right_skewed_data.html

# Should look into correlation between some of the variables.

# if ggfortify or gpubr doesnt work, do remotes::update_packages("rlang")

# 

