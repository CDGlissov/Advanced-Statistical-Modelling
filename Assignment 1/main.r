##### 02424 Assignment 1 #####
rm(list=ls())
# Loading packages
library(ggplot2, quietly = TRUE)
library(xtable);
library(ggpubr);
library(mgcv);
library(car);
library(grid);
library(gridExtra);
library(repr);
library(ggfortify);
library(GGally);
ggplot2::theme_set(ggplot2::theme_grey())
gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

##### Loading data #####
# Setting working directory
path="C:/Users/Christian/Dropbox/10. Semester/Dataanalyse og Statistisk Modellering/Afleveringer @ Statmod 2.0"
setwd(path)

# Loading the data
clothingSum <- read.csv("clothingSum.csv")
clothingSum$sex <-factor(clothingSum$sex)
clothingSum$subjId <-factor(clothingSum$subjId)

clothingFull <- read.csv("clothingFull.csv")
clothingFull$sex <- factor(clothingFull$sex)
clothingFull$subjId <- factor(clothingFull$subjId)
clothingFull$day <- factor(clothingFull$day)
clothingFull$obs.no <- factor(clothingFull$obs.no)

##### Explorative Analysis #####
# Plotting clothes as function of both indoor and outdoor temperature.
# Colours indicate sex
plotOut<-ggplot(clothingSum,aes(x=tOut,y=clo,col=sex))+geom_point()
plotIn<-ggplot(clothingSum,aes(x=tInOp,y=clo,col=sex))+geom_point()

ggarrange(plotIn,plotOut,ncol=1,nrow=2)

# Plotting boxplots
b.1=ggplot(clothingSum, aes(tOut,x=sex)) +
  geom_boxplot(aes(fill=sex))  +guides(fill=FALSE)+theme(axis.text=element_text(size=12),
                                                         axis.title=element_text(size=14,face="bold"))

b.2=ggplot(clothingSum, aes(tInOp,x=sex)) +
  geom_boxplot(aes(fill=sex))  +guides(fill=FALSE)+theme(axis.text=element_text(size=12),
                                                         axis.title=element_text(size=14,face="bold"))

b.3=ggplot(clothingSum, aes(x=sex, clo)) +
  geom_boxplot(aes(fill=sex)) +theme(axis.text=element_text(size=12),
                                     axis.title=element_text(size=14,face="bold"), 
                                     legend.box.just = c("top"), 
                                     legend.background = element_rect(fill=alpha('white', 0.4)),
                                     legend.position = c(0.87, 0.93))
b.3.plot=grid.arrange(b.1,b.2,b.3,ncol=3)
#ggsave("box_plot_clothing.png", height=5, width=12, plot=b.3.plot)

# density plots
#options(repr.plot.width=10, repr.plot.height=6)
p=ggpairs(clothingSum[,c(3,4,5,6)],
          upper = list(continuous=wrap("cor", size = 5,alignPercent=0.9)), 
          diag =list(discrete="barDiag", continuous = wrap("densityDiag", alpha=0.5 )),
          mapping = ggplot2::aes(colour=sex),
          lower=list(continuous ="smooth"))

dens.p1 = getPlot(p,1,1)+theme(legend.position="none")+ 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold"))
dens.p2 = getPlot(p,2,2)+theme(legend.position="none") + 
  theme(axis.text=element_text(size=12),
        axis.title=element_text(size=14,face="bold")) + ylab("")
dens.p3 = getPlot(p,3,3) + ylab("")+theme(axis.text=element_text(size=12),
                                          axis.title=element_text(size=14,face="bold"), 
                                          legend.background = element_rect(fill=alpha('white', 0.4)),
                                          legend.position = c(0.87, 0.93)) 

density.p=grid.arrange(dens.p1,dens.p2,dens.p3,ncol=3)

# Chekcing summmary statistics
summary(clothingSum)

##### General Linear Model #####
# Fitting a general linear model
par(mfrow=c(1,2))
Model1GAM<-gam(clo~s(tInOp) +s(tOut),data=clothingSum)
plot(Model1GAM)

# Should include tOut^2, tInOp^3 didn't make sense to include
fitM<-lm(clo~sex*I(tOut-mean(tOut))*I(tInOp-mean(tInOp))+I((tOut-mean(tOut))^2)*sex + 
              I((tInOp-mean(tInOp))^2)*sex, data=clothingSum)
Anova(fitM,type=3)
xtable(anova(fitM))

# Removing I((tInOp-mean(tInOp))^2):sex
fitM<-update(fitM,~. -I((tInOp-mean(tInOp))^2):sex)
Anova(fitM,type=3)

# Removing I((tOut-mean(tOut))^2):sex
fitM<-update(fitM,~. -I((tOut-mean(tOut))^2):sex)
Anova(fitM,type=3)

# Removing I((tInOp-mean(tInOp))^2)
fitM<-update(fitM,~. -I((tInOp-mean(tInOp))^2))
Anova(fitM,type=3)

# Removing sex:I(tOut - mean(tOut)):I(tInOp - mean(tInOp))
fitM<-update(fitM,~. -sex:I(tOut - mean(tOut)):I(tInOp - mean(tInOp)))
Anova(fitM,type=3)

# Removing I(tOut - mean(tOut)):I(tInOp - mean(tInOp))
fitM<-update(fitM,~. -I(tOut - mean(tOut)):I(tInOp - mean(tInOp)))
Anova(fitM,type=3)

# Removing sex:I(tOut - mean(tOut)) 
fitM<-update(fitM,~. -sex:I(tOut - mean(tOut)))
Anova(fitM,type=3)

# Removing I((tOut - mean(tOut))^2) 
fitM<-update(fitM,~. -I((tOut - mean(tOut))^2) )
Anova(fitM,type=3)

ModelFinalLM<-fitM

#to read the model without intercepts.
ModelFinalRead<-lm(clo~-1+sex+I(tOut-mean(tOut))+sex:I(tInOp - mean(tInOp)),data=clothingSum)

# We cannot reduce any more
summary(ModelFinalRead)
xtable(anova(ModelFinalLM))

# Making a new data-set
mintInOp<-min(clothingSum$tInOp)
maxtInOp<-max(clothingSum$tInOp)

male_tInOp<-data.frame(tInOp=seq(mintInOp,maxtInOp,length.out=length(clothingSum$clo)))
nMaleTinOp<-length(male_tInOp)
male_tInOp$tOut<-rep(mean(clothingSum$tOut),nMaleTinOp)
male_tInOp$sex<-rep("male",nMaleTinOp)
male_tInOp$pred<-predict(ModelFinalLM,male_tInOp)
male_tInOp$lower<-predict(ModelFinalLM, newdata = male_tInOp, interval = "confidence")[,2]
male_tInOp$upper<-predict(ModelFinalLM, newdata = male_tInOp, interval = "confidence")[,3]
male_tInOp$lowerPred<-predict(ModelFinalLM, newdata = male_tInOp, interval = "prediction")[,2]
male_tInOp$upperPred<-predict(ModelFinalLM, newdata = male_tInOp, interval = "prediction")[,3]

female_tInOp<-data.frame(tInOp=seq(mintInOp,maxtInOp,length.out=length(clothingSum$clo)))
nFemaletinOp<-length(female_tInOp)
female_tInOp$tOut<-rep(mean(clothingSum$tOut),nFemaletinOp)
female_tInOp$sex<-rep("female",nFemaletinOp)
female_tInOp$pred<-predict(ModelFinalLM,female_tInOp)
female_tInOp$lower<-predict(ModelFinalLM, newdata = female_tInOp, interval = "confidence")[,2]
female_tInOp$upper<-predict(ModelFinalLM, newdata = female_tInOp, interval = "confidence")[,3]
female_tInOp$lowerPred<-predict(ModelFinalLM, newdata = female_tInOp, interval = "prediction")[,2]
female_tInOp$upperPred<-predict(ModelFinalLM, newdata = female_tInOp, interval = "prediction")[,3]

maxtOut<-max(clothingSum$tOut)
mintOut<-min(clothingSum$tOut)

male_tOut<-data.frame(tOut=seq(mintOut,maxtOut,length.out=length(clothingSum$clo)))
nMaletOut<-length(male_tOut)
male_tOut$tInOp<-rep(mean(clothingSum$tInOp),nMaletOut)
male_tOut$sex<-rep("male",nMaletOut);
male_tOut$pred<-predict(ModelFinalLM,male_tOut)
male_tOut$lower<-predict(ModelFinalLM, newdata = male_tOut, interval = "confidence")[,2]
male_tOut$upper<-predict(ModelFinalLM, newdata = male_tOut, interval = "confidence")[,3]
male_tOut$lowerPred<-predict(ModelFinalLM, newdata = male_tOut, interval = "prediction")[,2]
male_tOut$upperPred<-predict(ModelFinalLM, newdata = male_tOut, interval = "prediction")[,3]

female_tOut<-data.frame(tOut=seq(mintOut,maxtOut,length.out=length(clothingSum$clo)))
nFemaletOut<-length(female_tOut)
female_tOut$tInOp<-rep(mean(clothingSum$tInOp),nFemaletOut)
female_tOut$sex<-rep("female",nFemaletOut);
female_tOut$pred<-predict(ModelFinalLM,female_tOut)
female_tOut$lower<-predict(ModelFinalLM, newdata = female_tOut, interval = "confidence")[,2]
female_tOut$upper<-predict(ModelFinalLM, newdata = female_tOut, interval = "confidence")[,3]
female_tOut$lowerPred<-predict(ModelFinalLM, newdata = female_tOut, interval = "prediction")[,2]
female_tOut$upperPred<-predict(ModelFinalLM, newdata = female_tOut, interval = "prediction")[,3]

# tInOp plot
# Adding lines
plottInOpLM<-ggplot(male_tInOp,aes(x=tInOp,y=pred))+
  geom_line(data=male_tInOp,aes(x=tInOp,y=pred))+
  geom_line(data=male_tInOp,aes(x=tInOp,y=lowerPred,col='male'),linetype=2)+
  geom_line(data=male_tInOp,aes(x=tInOp,y=upperPred,col='male'),linetype=2,show.legend=FALSE)+
  geom_ribbon(data=male_tInOp,aes(x=tInOp,ymin=lower,ymax=upper,alpha=0.1,fill='male'),show.legend=FALSE)
plottInOpLM<-plottInOpLM+
  geom_line(data=female_tInOp,aes(x=tInOp,y=pred))+
  geom_line(data=female_tInOp,aes(x=tInOp,y=lowerPred,col='female'),linetype=2)+
  geom_line(data=female_tInOp,aes(x=tInOp,y=upperPred,col='female'),linetype=2,show.legend=FALSE)+
  geom_ribbon(data=female_tInOp,aes(x=tInOp,ymin=lower,ymax=upper,alpha=0.1,fill='female'),show.legend=FALSE)
# Adding Points
plottInOpLM<-plottInOpLM+
  geom_point(data=clothingSum,aes(x=tInOp,y=clo,col=sex))+
  labs(y="Level of Clothing",colour="Sex")

#tOut plot
plottOutLM<-ggplot(male_tOut,aes(x=tOut,y=pred))+geom_line(data=male_tOut,aes(x=tOut,y=pred))+
  labs(y="Level of Clothing",colour="Sex")+
  geom_line(data=male_tOut,aes(x=tOut,y=lowerPred,col='male'),linetype=2)+
  geom_line(data=male_tOut,aes(x=tOut,y=upperPred,col='male'),linetype=2,show.legend=FALSE)+geom_ribbon(data=male_tOut,aes(x=tOut,ymin=lower,ymax=upper,alpha=0.1,fill='male'),show.legend=FALSE)
plottOutLM<-plottOutLM+geom_line(data=female_tOut,aes(x=tOut,y=pred))+geom_line(data=female_tOut,aes(x=tOut,y=lowerPred,col='female'),linetype=2)+geom_line(data=female_tOut,aes(x=tOut,y=upperPred,col='female'),linetype=2,show.legend=FALSE)+geom_ribbon(data=female_tOut,aes(x=tOut,ymin=lower,ymax=upper,alpha=0.1,fill='female'),show.legend=FALSE)
plottOutLM<-plottOutLM+geom_point(data=clothingSum,aes(x=tOut,y=clo,col=sex))+labs(y="Level of Clothing",colour="Sex")

ggarrange(plottInOpLM,plottOutLM,ncol=1,nrow=2,labels=c('A','B'))

# Finding the residuals
clothingSum$Res<-residuals(ModelFinalLM)
# Plotting residuals vs. sex
ggplot(clothingSum,aes(x =sex,y=Res,fill=sex))+geom_boxplot()+labs(y="Resisduals")

# Plotting Diagnoistics
autoplot(ModelFinalLM, colour="sex", which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")


## Weighted Model ##
fullModel<-lm(clo~sex*I(tOut-mean(tOut))*I(tInOp-mean(tInOp))+I((tOut-mean(tOut))^2)*sex + 
                       I((tInOp-mean(tInOp))^2)*sex, data=clothingSum)
X = model.matrix(fullModel)
n<-length(clothingSum$sex)

# ll of model
llWeight<-function(c,X){
  # Defining the sigma-matrix
  n<-length(clothingSum$sex)
  weightsVec=rep(1,n)
  weightsVec[clothingSum$sex=="female"]=c
  #tempFit<-lm(clo~sex*tOut*tInOp+tInOp2*sex+tOut2*sex+tInOp3*sex,data=clothingSum,weight=weightsVec)
  sigmaMat<-diag(weightsVec)
  beta<-solve(t(X)%*%solve(sigmaMat)%*%X)%*%t(X)%*%solve(sigmaMat)%*%clothingSum$clo
  # Estimating sigma
  m1<-dim(beta)[1]
  sigma<-as.numeric((t(clothingSum$clo-X%*%beta)%*%solve(sigmaMat)%*%(clothingSum$clo-X%*%beta))/(n-m1))
  mu <- X%*%beta
  ll<--(-n/2*log(2*pi)-(1/2)*log(det(sigma*sigmaMat))-1/(2*sigma)*t(clothingSum$clo-mu)%*%solve(sigmaMat)%*%(clothingSum$clo-mu))
  return(ll)
}


# Weight Statistics
getWeight<-function(X,plotWeights=FALSE){
  w = 1.5
  fun=llWeight
  n<-length(clothingSum$sex)
  optWeights<-optim(par = w,X=X, fn=fun,method='L-BFGS-B',lower=c(0.1),upper=c(10))
  weights=seq(0.1,5,0.01)
  dfWeights<-data.frame(weights=weights)
  dfWeights$logLik<-sapply(weights,fun,X)
  dfWeights$logLik<--(dfWeights$logLik-min(dfWeights$logLik))
  dfWeights$Critical<--qchisq(0.95,df=1)/2
  if(plotWeights == TRUE){
    ggplot(dfWeights,aes(x=weights,y=logLik))+
      geom_line()+geom_point(aes(x=optWeights$par,y=(0),col="Optimal"))+
      geom_line(aes(x=weights,y=Critical),linetype=2)+labs(x="c",y="Log-likelihood")+
      ylim(-10,0)
    
  }
  chiObs<- -2*(-fun(1,X)+fun(optWeights$par,X))
  p_val = (1-pchisq(chiObs,1))/2
  df_stat = data.frame(chiObs = chiObs, p_val = p_val, optWeight = optWeights$par)
  
  weightsVec=rep(1,n)
  #lm inverts the weights
  weightsVec[clothingSum$sex=="female"]=1/optWeights$par
  
  data_list = list("dfWeights" = dfWeights , "df_stat"=df_stat, "optWeights" = weightsVec)
  return(data_list)
}

# Fitting models with weights
df_weight=getWeight(X)
df_weight$df_stat

ModelWeight<-lm(clo~sex*I(tOut-mean(tOut))*I(tInOp-mean(tInOp))+I((tOut-mean(tOut))^2)*sex + 
                I((tInOp-mean(tInOp))^2)*sex, weights=df_weight$optWeights, data=clothingSum)
Anova(ModelWeight,type=3)
xtable(Anova(ModelWeight,type=3))

# Removing sex:I((tOut - mean(tOut))^2)
# 
ModelWeight<-update(ModelWeight,~. -sex:I((tOut - mean(tOut))^2))
X = model.matrix(ModelWeight)
df_weight=getWeight(X)
df_weight$df_stat
ModelWeight<-update(ModelWeight,~.,weights=df_weight$optWeights)
Anova(ModelWeight,type=3)


# Removing sex:I((tInOp - mean(tInOp))^2)
ModelWeight<-update(ModelWeight,~. -sex:I((tInOp - mean(tInOp))^2))
X = model.matrix(ModelWeight)
df_weight=getWeight(X)
df_weight$df_stat
ModelWeight<-update(ModelWeight,~.,weights=df_weight$optWeights)
Anova(ModelWeight,type=3)


# Removing I((tInOp - mean(tInOp))^2) 
ModelWeight<-update(ModelWeight,~. -I((tInOp - mean(tInOp))^2))
X = model.matrix(ModelWeight)
df_weight=getWeight(X)
df_weight$df_stat
ModelWeight<-update(ModelWeight,~.,weights=df_weight$optWeights)
Anova(ModelWeight,type=3)

# Removing sex:I(tOut - mean(tOut)):I(tInOp - mean(tInOp)) 
ModelWeight<-update(ModelWeight,~. -sex:I(tOut - mean(tOut)):I(tInOp - mean(tInOp)))
X = model.matrix(ModelWeight)
df_weight=getWeight(X)
df_weight$df_stat
ModelWeight<-update(ModelWeight,~.,weights=df_weight$optWeights)
Anova(ModelWeight,type=3)

# Removing I(tOut - mean(tOut)):I(tInOp - mean(tInOp))
ModelWeight<-update(ModelWeight,~. -I(tOut - mean(tOut)):I(tInOp - mean(tInOp)))
X = model.matrix(ModelWeight)
df_weight=getWeight(X)
df_weight$df_stat
ModelWeight<-update(ModelWeight,~.,weights=df_weight$optWeights)
Anova(ModelWeight,type=3)

# Removing sex:I(tOut - mean(tOut)) 
ModelWeight<-update(ModelWeight,~. -sex:I(tOut - mean(tOut)))
X = model.matrix(ModelWeight)
df_weight=getWeight(X)
df_weight$df_stat
ModelWeight<-update(ModelWeight,~.,weights=df_weight$optWeights)
Anova(ModelWeight,type=3)

#Plotting
autoplot(ModelWeight, colour="sex", which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

# Finding the model residuals
clothingSum$ResW<-residuals(ModelWeight)
# Plotting residuals vs. sex
resSexW<-ggplot(clothingSum,aes(x =sex,y=ResW,fill=sex))+geom_boxplot()+labs(y="Resisduals")
# Plotting Residuals vs. tInOp and tOut
restInOpW<-ggplot(clothingSum,aes(x =tInOp,y=ResW,col=sex))+geom_point()+labs(y="Residuals")
restOutW<-ggplot(clothingSum,aes(x =tOut,y=ResW,col=sex))+geom_point()+labs(y="Residuals")
ggarrange(resSexW,restInOpW,restOutW,ncol=1,nrow=3,labels=c("A","B","C"))

# Plots
weights_male = rep(1,n)
weights_female=1/(rep(1,n)*df_weight$df_stat$optWeight)

male_tInOp$predW<-predict(ModelWeight,male_tInOp)
male_tInOp$lowerW<-predict(ModelWeight, newdata = male_tInOp, interval = "confidence", weights=weights_male)[,2]
male_tInOp$upperW<-predict(ModelWeight, newdata = male_tInOp, interval = "confidence", weights=weights_male)[,3]
male_tInOp$lowerPredW<-predict(ModelWeight, newdata = male_tInOp, interval = "prediction", weights=weights_male)[,2]
male_tInOp$upperPredW<-predict(ModelWeight, newdata = male_tInOp, interval = "prediction", weights=weights_male)[,3]

female_tInOp$predW<-predict(ModelWeight,female_tInOp)
female_tInOp$lowerW<-predict(ModelWeight, newdata = female_tInOp, interval = "confidence",weights=weights_female)[,2]
female_tInOp$upperW<-predict(ModelWeight, newdata = female_tInOp, interval = "confidence",weights=weights_female)[,3]
female_tInOp$lowerPredW<-predict(ModelWeight, newdata = female_tInOp, interval = "prediction",weights=weights_female)[,2]
female_tInOp$upperPredW<-predict(ModelWeight, newdata = female_tInOp, interval = "prediction",weights=weights_female)[,3]

male_tOut$predW<-predict(ModelWeight,male_tOut)
male_tOut$lowerW<-predict(ModelWeight, newdata = male_tOut, interval = "confidence", weights=weights_male)[,2]
male_tOut$upperW<-predict(ModelWeight, newdata = male_tOut, interval = "confidence", weights=weights_male)[,3]
male_tOut$lowerPredW<-predict(ModelWeight, newdata = male_tOut, interval = "prediction", weights=weights_male)[,2]
male_tOut$upperPredW<-predict(ModelWeight, newdata = male_tOut, interval = "prediction", weights=weights_male)[,3]

female_tOut$predW<-predict(ModelWeight,female_tOut)
female_tOut$lowerW<-predict(ModelWeight, newdata = female_tOut, interval = "confidence",weights=weights_female)[,2]
female_tOut$upperW<-predict(ModelWeight, newdata = female_tOut, interval = "confidence",weights=weights_female)[,3]
female_tOut$lowerPredW<-predict(ModelWeight, newdata = female_tOut, interval = "prediction",weights=weights_female)[,2]
female_tOut$upperPredW<-predict(ModelWeight, newdata = female_tOut, interval = "prediction",weights=weights_female)[,3]

# tInOp plot, Adding lines
plottInOpLM<-ggplot(male_tInOp,aes(x=tInOp,y=predW))+
  geom_line(data=male_tInOp,aes(x=tInOp,y=predW))+
  geom_line(data=male_tInOp,aes(x=tInOp,y=lowerPredW,col='male'),linetype=2)+
  geom_line(data=male_tInOp,aes(x=tInOp,y=upperPredW,col='male'),linetype=2,show.legend=FALSE)+
  geom_ribbon(data=male_tInOp,aes(x=tInOp,ymin=lowerW,ymax=upperW,alpha=0.1,fill='male'),show.legend=FALSE)
plottInOpLM<-plottInOpLM+geom_line(data=female_tInOp,aes(x=tInOp,y=predW))+
  geom_line(data=female_tInOp,aes(x=tInOp,y=lowerPredW,col='female'),linetype=2)+
  geom_line(data=female_tInOp,aes(x=tInOp,y=upperPredW,col='female'),linetype=2,show.legend=FALSE)+
  geom_ribbon(data=female_tInOp,aes(x=tInOp,ymin=lowerW,ymax=upperW,alpha=0.1,fill='female'),show.legend=FALSE)
plottInOpLM<-plottInOpLM+geom_point(data=clothingSum,aes(x=tInOp,y=clo,col=sex))+
  labs(y="Level of Clothing",colour="Sex")

#tOut plot
plottOutLM<-ggplot(male_tOut,aes(x=tOut,y=predW))+
  geom_line(data=male_tOut,aes(x=tOut,y=predW))+
  geom_line(data=male_tOut,aes(x=tOut,y=lowerPredW,col='male'),linetype=2)+
  geom_line(data=male_tOut,aes(x=tOut,y=upperPredW,col='male'),linetype=2,show.legend=FALSE)+
  geom_ribbon(data=male_tOut,aes(x=tOut,ymin=lowerW,ymax=upperW,alpha=0.1,fill='male'),show.legend=FALSE)
plottOutLM<-plottOutLM+geom_line(data=female_tOut,aes(x=tOut,y=predW))+
  geom_line(data=female_tOut,aes(x=tOut,y=lowerPredW,col='female'),linetype=2)+
  geom_line(data=female_tOut,aes(x=tOut,y=upperPredW,col='female'),linetype=2,show.legend=FALSE)+
  geom_ribbon(data=female_tOut,aes(x=tOut,ymin=lowerW,ymax=upperW,alpha=0.1,fill='female'),show.legend=FALSE)
plottOutLM<-plottOutLM+geom_point(data=clothingSum,aes(x=tOut,y=clo,col=sex))+
  labs(y="Level of Clothing",colour="Sex")

ggarrange(plottInOpLM,plottOutLM,ncol=1,nrow=2,labels=c('A','B'))



######## Including subject ID #########
par(mfrow=c(1,1))
ggplot(clothingSum,aes(x=subjId,y=residuals(ModelWeight)))+geom_boxplot(aes(fill=sex))+
  theme(axis.text.x = element_text(angle = 90, hjust = 1))+labs(y="Residuals of weighted model", x="Subject ID")

modelSubject<-lm(clo ~ -1 + subjId + sex:I(tInOp-mean(tInOp)) + sex:I(tOut-mean(tOut)) + I(tInOp-mean(tInOp)) + 
                   I(tOut-mean(tOut)) + I((tInOp-mean(tInOp))^2) + I((tOut-mean(tOut))^2),data = clothingSum)
X = model.matrix(modelSubject)
df_weight=getWeight(X)
df_weight$df_stat
modelSubject<-update(modelSubject,~.,weights=df_weight$optWeights)
Anova(modelSubject,type=3)


# Removing I(tOut - mean_tOut):sex
modelSubject<-update(modelSubject,~. -I(tOut - mean(tOut)):sex)
X = model.matrix(modelSubject)
df_weight=getWeight(X)
df_weight$df_stat
modelSubject<-update(modelSubject,~.,weights=df_weight$optWeights)
Anova(modelSubject,type=3)

# Removing I((tInOp - mean_tInOp)^2)
modelSubject<-update(modelSubject,~. -I((tInOp - mean(tInOp))^2))
X = model.matrix(modelSubject)
df_weight=getWeight(X)
df_weight$df_stat
modelSubject<-update(modelSubject,~.,weights=df_weight$optWeights)
Anova(modelSubject,type=3)

# Removing I(tInOp - mean_tInOp):sex
modelSubject<-update(modelSubject,~. -I(tInOp - mean(tInOp)):sex)
X = model.matrix(modelSubject)
df_weight=getWeight(X)
df_weight$df_stat
modelSubject<-update(modelSubject,~.,weights=df_weight$optWeights)
Anova(modelSubject,type=3)

# Removing I(tInOp - mean_tInOp)
modelSubject<-update(modelSubject,~. -I(tInOp - mean(tInOp)))
X = model.matrix(modelSubject)
df_weight=getWeight(X)
df_weight$df_stat
modelSubject<-update(modelSubject,~.,weights=df_weight$optWeights)
Anova(modelSubject,type=3)

autoplot(modelSubject, colour=gg_color_hue(2)[clothingSum$sex], which=c(1:3,5), nrow=2,ncol=2)+theme(legend.position="none")

# Making some plots of the intercepts
# SexByID has the sex of the IDs
sexByID <- rep(0,length(levels(clothingSum$subjId)))
i <- 0
for (ID in levels(clothingSum$subjId)) {
  i = i + 1
  if (clothingSum$sex[clothingSum$subjId==ID][1] == "female") {
    sexByID[i] <- "female"
  } else {
    sexByID[i] <- "male"
  }
}

# MF has the individual intercepts and corresponding sex
SIbeta0s <- coefficients(modelSubject)[1:length(levels(clothingSum$subjId))]
MF <- data.frame(IS = SIbeta0s, sex = sexByID)

# Boxplot
ggplot(MF,aes(x=sex,y=SIbeta0s,fill = sex)) +
  geom_boxplot() + 
  labs(title = "Intercepts by sex")

# Density
ggplot(MF,aes(x = IS, fill = sex)) + 
  theme_bw() + 
  geom_density(alpha = 0.5) + 
  labs(title = "Densities of the intercepts by sex",x="Intercept",y="Density")

# Predictions using mean of male and females respectively
maleMean <- mean(MF$IS[MF$sex == "male"])
femaleMean <- mean(MF$IS[MF$sex == "female"])
# beta1 and beta2
SIbeta2 <- modelSubject$coefficients[length(modelSubject$coefficients)]
SIbeta1 <- modelSubject$coefficients[length(modelSubject$coefficients)-1]

predsMale <- data.frame(tOut = clothingSum$tOut, line = maleMean + SIbeta1*(clothingSum$tOut-mean(clothingSum$tOut)) + SIbeta2*(clothingSum$tOut-mean(clothingSum$tOut))^2)
predsFemale <- data.frame(tOut = clothingSum$tOut, line = femaleMean + SIbeta1*(clothingSum$tOut-mean(clothingSum$tOut)) + SIbeta2*(clothingSum$tOut-mean(clothingSum$tOut))^2)
ggplot(clothingSum,aes(x = tOut, y = clo, col = sex)) +
  theme_bw() +
  geom_jitter() +
  geom_line(data = predsMale, aes(x=tOut,y=line), col = gg_color_hue(2)[2]) + 
  geom_line(data = predsFemale, aes(x=tOut,y=line), col = gg_color_hue(2)[1]) + 
  labs(title = "Predictions using mean intercepts for males and females")

# Looking at the residuls
clothingSum$ResSI<-residuals(modelSubject)

# Plotting residuals vs. sex
resSexSI<-ggplot(clothingSum,aes(x =sex,y=ResSI,fill=sex))+geom_boxplot()+labs(y="Residuals")

# Plotting Residuals vs. tInOp and tOut
restInOpSI<-ggplot(clothingSum,aes(x =tInOp,y=ResSI,col=sex))+geom_point()+labs(y="Residuals")
restOutSI<-ggplot(clothingSum,aes(x =tOut,y=ResSI,col=sex))+geom_point()+labs(y="Residuals")

ggarrange(resSexW,restInOpSI,restOutSI,ncol=1,nrow=3,labels=c("A","B","C"))

####### FULL DATA SET ########

#PLOTTING
plotOut<-ggplot(clothingFull,aes(x=tOut,y=clo,col=sex))+geom_point()+theme(legend.position="none")
plotIn<-ggplot(clothingFull,aes(x=tInOp,y=clo,col=sex))+geom_point()+theme(legend.position="none")
plotInOut<-ggplot(clothingFull,aes(x=tInOp,y=tOut,col=sex))+geom_point()+theme(legend.position="none")
# we should still weight it
grid.arrange(plotIn,plotOut,plotInOut,ncol=3,nrow=1)

# Plotting boxplots
boxPlotIn<- ggplot(clothingFull,aes(x =sex,y=tInOp))+
  theme(legend.position="none") + geom_boxplot(aes(fill=sex))
boxPlotOut<-ggplot(clothingFull,aes(x=sex,y=tOut))+
  theme(legend.position="none")+ geom_boxplot(aes(fill=sex))
boxPlotClo<-ggplot(clothingFull,aes(x=sex,y=clo))+
  theme(legend.position="none") + geom_boxplot(aes(fill=sex))
boxPlotDay<- ggplot(clothingFull,aes(x =day,y=clo))+
  theme(legend.position="none")+ geom_boxplot(aes(fill=sex))
boxPlotDay1<- ggplot(clothingFull,aes(x =day,y=tOut))+
  theme(legend.position="none")+ geom_boxplot(aes(fill=sex))
boxPlotDay2<- ggplot(clothingFull,aes(x =day,y=tInOp))+ 
  geom_boxplot(aes(fill=sex))
boxPlotsub1<- ggplot(clothingFull,aes(x =subjId,y=clo))+
  theme(legend.position="none")+ theme(axis.text.x = element_text(angle = 90))+ 
  geom_boxplot(aes(fill=sex))
boxPlotsub2<- ggplot(clothingFull,aes(x =subjId,y=tOut))+
  theme(legend.position="none")+ theme(axis.text.x = element_text(angle = 90))+ 
  geom_boxplot(aes(fill=sex))
boxPlotsub3<- ggplot(clothingFull,aes(x =subjId,y=tInOp))+geom_boxplot()+
  theme(legend.position="none")+ theme(axis.text.x = element_text(angle = 90))+ 
  geom_boxplot(aes(fill=sex))
grid.arrange(boxPlotIn, boxPlotOut,boxPlotClo,boxPlotDay,boxPlotDay1,boxPlotDay2,boxPlotsub1,boxPlotsub2,boxPlotsub3,ncol=3,nrow=3)

boxPlotobs1<- ggplot(clothingFull,aes(x =obs.no,y=clo))+theme(legend.position="none")+ 
  geom_boxplot(aes(fill=sex))
boxPlotobs2<- ggplot(clothingFull,aes(x =obs.no,y=tOut))+theme(legend.position="none")+ 
  geom_boxplot(aes(fill=sex))
boxPlotobs3<- ggplot(clothingFull,aes(x =obs.no,y=tInOp))+theme(legend.position="none")+ 
  geom_boxplot(aes(fill=sex))
grid.arrange(boxPlotobs1, boxPlotobs2,boxPlotobs3,ncol=3,nrow=1)

min.ll = function(wt, fit){
  n=length(clothingFull$clo)
  wopt<-rep(1,n)
  wopt[clothingFull$sex=="female"]=wt
  form=as.formula(paste("clo", "~", as.character(fit$call$formula))[3])
  tempFit<-lm(form, data=clothingFull, weights=wopt)
  ll= -logLik(tempFit)
  return(ll)
}

getWeightFull <- function(fit,plotWeights=FALSE){
  n<-length(clothingFull$sex)
  wt=rep(1,n)
  optWeights<-optim(par = c(1.5),fit=fit, fn = min.ll, method='L-BFGS-B', lower=c(0.1), upper=c(200))
  
  weights=seq(0.1,5,0.01)
  dfWeights<-data.frame(weights=weights)
  dfWeights$logLik<-sapply(weights,min.ll, fit)
  dfWeights$logLik<--(dfWeights$logLik-min(dfWeights$logLik))
  dfWeights$Critical<--qchisq(0.95,df=1)/2
  if(plotWeights == TRUE){
    ggplot(dfWeights,aes(x=weights,y=logLik))+
      geom_line()+geom_point(aes(x=optWeights$par,y=(0),col="Optimal"))+
      geom_line(aes(x=weights,y=Critical),linetype=2)+labs(x="c",y="Log-likelihood")+
      ylim(-10,0)
  }
  chiObs<- -2*(-min.ll(1,fit)+min.ll(optWeights$par,fit))
  p_val = (1-pchisq(chiObs,1))/2
  df_stat = data.frame(chiObs = chiObs, p_val = p_val, optWeight = optWeights$par)
  wt[clothingFull$sex=="female"]=1/optWeights$par
  data_list = list("dfWeights" = dfWeights , "df_stat"=df_stat, "optWeights" = wt)
  return(data_list)
}

#### Subject ID full ####
wt <- rep(1,length(clothingFull$tOut))
FullModelSubject <- lm(clo ~ -1 + subjId + I(tOut-mean(tOut)) + I((tOut-mean(tOut))^2), data = clothingFull, weights = wt)
wt <- getWeightFull(FullModelSubject)
FullModelSubject <- update(FullModelSubject,weights = wt$optWeights)
par(mfrow=c(2,2))

autoplot(FullModelSubject, colour=gg_color_hue(2)[clothingFull$sex], which=c(1:3,4), nrow=2,ncol=2)+theme(legend.position="none")


FullModelWeight<-lm(clo ~ sex + I(tOut - mean(tOut)) + I(tInOp - mean(tInOp)) + I((tOut -mean(tOut))^2) + sex:I(tInOp - mean(tInOp)), data=clothingFull)
wt <- getWeightFull(FullModelWeight)
FullModelWeight<-update(FullModelWeight,~.,weights=wt$optWeights)
Anova(ModelWeight,type=3)

autoplot(FullModelSubject, colour=gg_color_hue(2)[clothingFull$sex], which=c(1:3,4), nrow=2,ncol=2)+theme(legend.position="none")



