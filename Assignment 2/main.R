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
par(mfrow=c(3,3))
# Possibly 2. order for Pres variable.
gam(Ozone ~ s(Temp)+s(InvHt)+s(Pres)+s(Vis)+s(Hgt)+s(Hum)+s(InvTmp)+s(Wind),data=ozone)

# Initial model; 2 way interactions with Pres squared
LMInit <- lm(Ozone ~ .*.*I(Pres^2), data = ozone)

#Reduce the model
LM1=model.select(LMInit)
summary(LM1)

#Residuals
autoplot(LM1, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

# Trying box-cox:
par(mfrow=c(1,1))
l = boxCox(LM1, lambda = seq(0,1,0.01))
l_opt=l$x[l$y==max(l$y)]
# Seems 1/3 i.e the cubic-root transformation
LM1_bc <- glm(Ozone^(l_opt) ~ Temp + InvHt + Pres + Hum + InvTmp + I(Pres^2) + 
               Temp:InvTmp + Hum:InvTmp, family=gaussian, data = ozone)

# Now Temp:InvTmp becomes insignificant
drop1(LM1_bc, test = "F")
LM1_bc <- update(LM1_bc, .~. - Temp:InvTmp)

#Residuals have improved.
summary(LM1_bc)
autoplot(LM1_bc, which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

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


##### Clothing Insulation Level: Count data #####

# Loading data
clo_count <- read.csv("~/Library/Mobile Documents/com~apple~CloudDocs/DTU/10.semester/02424VideregåendeDataanalyseStatistiskModellering/Advanced-Statistical-Modelling/Assignment 2/dat_count.csv", sep=";")
# clo is now the number of times a person changes clothes 

#### 1: Develop a Generalized Model based on the Binomial Distribution #####
## Model 1 ##
# Ignoring subject ID and day
model<-cbind(clo,nobs-clo)~time+sex+tOut+tInOp
Binom.glm<-glm(model,family=binomial,
               data=clo_count)

# Reducing the model
anova(Binom.glm,test="Chisq")

# Remove tOut
Binom.glm<-update(Binom.glm,~. -tOut)
anova(Binom.glm,test="Chisq")

# Removing tInOp
Binom.glm<-update(Binom.glm,~. -tInOp)
anova(Binom.glm,test="Chisq")

# Removing time 
Binom.glm<-update(Binom.glm,~. -time)
anova(Binom.glm,test="Chisq")
summary(Binom.glm)
pval<-1- pchisq(172.57,134)
pval
# This p-value is to small

# Checking the residuals 
plot(Binom.glm)

## Model 2: Including over dispersion ##
Quasibinom.glm<-glm(model,family=quasibinomial,
               data=clo_count)

summary(Quasibinom.glm)

anova(Quasibinom.glm,test="Chisq")

# Remove tOut
Quasibinom.glm<-update(Quasibinom.glm,~. -tOut)
anova(Quasibinom.glm,test="Chisq")

# Remove tInOp
Quasibinom.glm<-update(Quasibinom.glm,~. -tInOp)
anova(Quasibinom.glm,test="Chisq")

# Remove time 
Quasibinom.glm<-update(Quasibinom.glm,~. -time)
anova(Quasibinom.glm,test="Chisq")

plot(Quasibinom.glm)


## Model 3: Using another link-function ##
# Check which other link function might be appropriate
anova(glm(model,family = binomial,data=clo_count))
anova(glm(model,family=binomial(probit),data=clo_count))
anova(glm(model,family=binomial(cauchit),data=clo_count))
anova(glm(model,family=binomial(cloglog),data=clo_count))

# Try the log-log (aloth)
cloglog.glm<-glm(model,data=clo_count,
                 family=binomial(cloglog))

# Reduce the model
anova(cloglog.glm,test="Chisq")

# Remove tOut
cloglog.glm<-update(cloglog.glm,~. -tOut)
anova(cloglog.glm,test="Chisq")

# Remove tInOp
cloglog.glm<-update(cloglog.glm,~. -tInOp)
anova(cloglog.glm,test="Chisq")

# Remove time 
cloglog.glm<-update(cloglog.glm,~. -time)
anova(cloglog.glm,test="Chisq")

# Checking for sufficiency
summary(cloglog.glm)
1- pchisq(172.57,134)
# Exactly the same problem as before... 

#### 2: Develop a generalized linear model with the Poisson model ####

# QUESTIONS AND NOTES ######################################################################
# Data is not in metric units. Should we convert units to metric, for inference?

# Ozone is positive, not a many values are 0, inverse gaussian or gamma distribution?
# https://cran.r-project.org/web/packages/GlmSimulatoR/vignettes/dealing_with_right_skewed_data.html

# Should look into correlation between some of the variables.

# if ggfortify or gpubr doesnt work, do remotes::update_packages("rlang")




