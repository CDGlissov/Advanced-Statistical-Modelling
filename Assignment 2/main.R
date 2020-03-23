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
source("reduce.R")
# Load data
dat_count = read.csv("dat_count.csv", header=TRUE, sep=";")
data(ozone)
dat_count[,c(1,2,3,5,6)] <- lapply(dat_count[,c(1,2,3,5,6)], factor)
head(ozone)
head(dat_count)

# SUMMARY STATISTICS 1.1 ###################################################################
summary(ozone)
summary(dat_count)

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
# TASK 1.2-1.3 #############################################################################

# Initial model; 2 way interactions with Pres squared
LMInit <- lm(Ozone ~ . , data = ozone)
anova(LMInit)

#Reduce the model
LM1=model.select(LMInit)
summary(LM1)
# In the end only temperature, inversion base height, and humidity are significant

# Diagnostics
autoplot(LM1, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

# Trying box-cox:
par(mfrow=c(1,1))
l = boxCox(LM1, lambda = seq(0,1,0.01))
l_opt=l$x[l$y==max(l$y)]
# Seems around 1/3 i.e the cubic-root transformation
LM1_bc <- glm(Ozone^(l_opt) ~ Temp + InvHt + Hum, family=gaussian, data = ozone)

# Checking if everything is still significant
drop1(LM1_bc, test = "F")

# Residuals have improved.
summary(LM1_bc)
autoplot(LM1_bc, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")


# GENERALIZED MODEL 1.4-1.5 ################################################################

## Gamma family
# Identity link function not working - trying log and inverse link functions
# Trying the log link function
GLMGam_init_log <- glm(Ozone ~ . , family = Gamma(link="log"), data = ozone)
GLMGam_log=model.select(GLMGam_init_log)
summary(GLMGam_log)
autoplot(GLMGam_init_log, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

# Trying the inverse link function
GLMGam_init_inv <- glm(Ozone ~ . , family = Gamma(link="inverse"), data = ozone)
GLMGam_inv=model.select(GLMGam_init_inv)
summary(GLMGam_inv)
autoplot(GLMGam_init_inv, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

## Inverse gaussian family
# Only the log link function seems to be working
# 1/mu^2, inverse, identity and log.
InvGaus_init <- glm(Ozone ~ ., family = inverse.gaussian("log"), data = ozone)
InvGaus=model.select(InvGaus_init)
summary(InvGaus)
autoplot(InvGaus, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

# COMPARING MODELS 1.5 #####################################################################

# The classical linear model
LM1_bcJacobian <- sum(log(BoxCoxData(ozone$Ozone,l_opt)$jacobian))
AICRegModel(unclass(logLik(LM1_bc))[1],LM1_bcJacobian,5)

# The generalized models
AIC(GLMGam_init_log)
AIC(GLMGam_init_inv)
AIC(InvGaus)
# Looks like we are going to use the classical glm.


# QUESTIONS AND NOTES ######################################################################
# Data is not in metric units. Should we convert units to metric, for inference?

# Ozone is positive, not a many values are 0, inverse gaussian or gamma distribution?
# https://cran.r-project.org/web/packages/GlmSimulatoR/vignettes/dealing_with_right_skewed_data.html

# Should look into correlation between some of the variables.

# if ggfortify or gpubr doesnt work, do remotes::update_packages("rlang")



par(mfrow=c(3,3))
# Possibly 2. order for Pres variable.
plot(gam(Ozone ~ s(Temp)+s(InvHt)+s(Pres)+s(Vis)+s(Hgt)+s(Hum)+s(InvTmp)+s(Wind),data=ozone))

# Initial model; 2 way interactions with Pres squared
LMInit <- lm(Ozone ~ .*.*I(Pres^2), data = ozone)

# GENERALIZED MODEL 1.4-1.5 ################################################################

GLMGam_init <- glm(Ozone ~ .*.*I(Pres^2), family = Gamma(link="identity"), data = ozone)
#inverse/log/identiy
GLMGam=model.select(GLMGam_init)
summary(GLMGam)
autoplot(GLMGam, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

InvGaus_init <- glm(Ozone ~ .*.*I(Pres^2), family = inverse.gaussian(link="inverse"), data = ozone)
#1/mu^2, inverse, identity and log.
InvGaus=model.select(InvGaus_init)
summary(InvGaus)
autoplot(InvGaus, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")
