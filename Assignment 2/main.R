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
#

path="C:/Users/Christian/Dropbox/10. Semester/Dataanalyse og Statistisk Modellering/Assignments/Assignment 2"
setwd(path)
dat_count = read.csv("dat_count.csv", header=TRUE, sep=";")
data(ozone)

dat_count[,c(1,2,3,5,6)] <- lapply(dat_count[,c(1,2,3,5,6)], factor)

head(ozone)
head(dat_count)
# SUMMARY STATISTICS #######################################################################
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
# MODELLING ################################################################################

#### General linear model ####
# Checking for higher order terms
par(mfrow=c(3,3))
Model1GAM<-gam(Ozone ~ s(Temp)+s(InvHt)+s(Pres)+s(Vis)+s(Hgt)+s(Hum)+s(InvTmp)+s(Wind),data=ozone)
plot(Model1GAM)
# Possibly 2 order for Pres variable.

# Initial model; 2 way interactions with Pres squared
model <- lm(Ozone ~ .*.*I(Pres^2), data = ozone)
summary(model)

# Reducing model first using "step"
model <- step(model, trace = 1)

# The rest manually
drop1(model, test = "F")
model <- update(model, .~. - InvHt:Hum:I(Pres^2))

drop1(model, test = "F")
model <- update(model, .~. - InvHt:Hum)

drop1(model, test = "F")
model <- update(model, .~. - Hgt:InvTmp:I(Pres^2))

drop1(model, test = "F")
model <- update(model, .~. - Hgt:InvTmp)

drop1(model, test = "F")
model <- update(model, .~. - InvHt:Hgt:I(Pres^2))

drop1(model, test = "F")
model <- update(model, .~. - InvHt:Hgt)

drop1(model, test = "F")
model <- update(model, .~. - Hgt:I(Pres^2))

drop1(model, test = "F")
model <- update(model, .~. - Hgt)

drop1(model, test = "F")
model <- update(model, .~. - Hum:InvTmp:I(Pres^2))

drop1(model, test = "F")
model <- update(model, .~. - Temp:InvHt:I(Pres^2))

drop1(model, test = "F")
model <- update(model, .~. - InvHt:I(Pres^2))

drop1(model, test = "F")
model <- update(model, .~. - Temp:InvTmp:I(Pres^2))

drop1(model, test = "F")
model <- update(model, .~. - Temp:InvHt)

drop1(model, test = "F")
model <- update(model, .~. - Temp:I(Pres^2))

drop1(model, test = "F")
model <- update(model, .~. - Vis)

drop1(model, test = "F")
model <- update(model, .~. - InvTmp:Wind:I(Pres^2))

drop1(model, test = "F")
model <- update(model, .~. - InvTmp:Wind)

drop1(model, test = "F")
model <- update(model, .~. - Wind:I(Pres^2))

drop1(model, test = "F")
model <- update(model, .~. - Wind)

drop1(model, test = "F")
model <- update(model, .~. - InvHt:Pres)

drop1(model, test = "F")
model <- update(model, .~. - InvTmp:I(Pres^2))

drop1(model, test = "F")
model <- update(model, .~. - Hum:I(Pres^2))

drop1(model, test = "F")
summary(model)

par(mfrow = c(2,2))
plot(model)
# Aaaand what do you know? We need a transformation! Trying box-cox:
par(mfrow=c(1,1))
boxCox(model, lambda = seq(0,1,0.01))
# Seems 1/3 would do the job i.e. the cubic-root transformation

modelT <- lm(formula = I(Ozone^(1/3)) ~ Temp + InvHt + Pres + Hum + InvTmp + I(Pres^2) + 
               Temp:InvTmp + Hum:InvTmp, data = ozone)
drop1(modelT, test = "F") # Now Temp:InvTmp becomes insignificant

modelT <- update(modelT, .~. - Temp:InvTmp)
drop1(modelT, test = "F")
summary(modelT)

# Diagnostics plot again
par(mfrow=c(2,2))
plot(modelT)
# ERROR - but residuals look better.




# QUESTIONS AND NOTES ######################################################################
# Data is not in metric units. Should we convert units to metric, for inference?
# Ozone is positive, not a many values are 0, inverse gaussian or gamma distribution?
# Should look into correlation between some of the variables.


