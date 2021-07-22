library(fda)
library(tidyverse)
set.seed(01327001)

data("CanadianWeather")

str(CanadianWeather)

subs<-c(1,15,32,46,60,74,91,105,121,135,152,166,
        182,196,213,227,244,258,274,288,305,319,335,349)


avgtemp<-t(CanadianWeather$dailyAv[,,1])[,subs]

matplot(t(avgtemp),type="l",col=(CanadianWeather$place%in%c("Dawson","Scheffervll"))+1)

y<-log(colSums(CanadianWeather$dailyAv[,,2]))
y<-y-mean(y)

#X<-vector("list",2)
#names(X)<-c("x","pp")

#X$x<-split(avgtemp,row(avgtemp))
#X$pp<-split(rep(subs,each=35),row(avgtemp))

#irregular
sparseObs<-sample(16:24,35,replace=T)
X<-fdapace::Sparsify(avgtemp,subs,sparseObs)
names(X)<-c("pp","x")

trainInd<-1:35
#testInd<-(1:35)[-trainInd]
testInd<-trainInd
Xtrain<-X
Xtrain$x<-X$x[trainInd]
Xtrain$pp<-X$pp[trainInd]
ytrain<-y[trainInd]

Xtest<-X
Xtest$x<-X$x[testInd]
Xtest$pp<-X$pp[testInd]
ytest<-y[testInd]

source("./functions.R")
K<-8

resFPCR<-doFPCReg(Xtrain,ytrain,Xtest,ytest,cl=subs,K=K,method="Irreg")

resRFPCR<-doRFPCReg(Xtrain,ytrain,Xtest,ytest,cl=subs,K=K,method="Irreg")

hs<-c(15,60,180)

resSLS<-doSparseFPCAReg(Xtrain,ytrain,Xtest,ytest,K=K,hs.mu = hs,hs.cov = hs,ncov=50,method="LS")

resSRob<-doSparseFPCAReg(Xtrain,ytrain,Xtest,ytest,K=K,hs.mu = hs,hs.cov = hs,ncov=50,method="Rob")
  geom_point()+
  geom_line()+
  labs(x="Number of Components",y="LOOCV Error")+
  theme(text = element_text(size = 20))
dev.off()
