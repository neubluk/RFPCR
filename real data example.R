library(fda)
library(tidyverse)
set.seed(01327001)

data("CanadianWeather")

str(CanadianWeather)

subs<-c(1,15,32,46,60,74,91,105,121,135,152,166,
        182,196,213,227,244,258,274,288,305,319,335,349)

#columns 11-16
#rows 25,28,29

avgtemp<-t(CanadianWeather$dailyAv[,,1])[,subs]
# 
# 
 matplot(t(avgtemp),type="l",col=(CanadianWeather$place%in%c("Dawson","Scheffervll"))+1)
# 
 y<-log(colSums(CanadianWeather$dailyAv[,,2]))
 y<-y-mean(y)
hist(y,freq=F)
#abline(v=y[(CanadianWeather$place%in%c("Pr. Rupert","Dawson","Scheffervll"))])
phist<-ggplot(data.frame(y=y),aes(x=y))+
  geom_histogram(aes(y=..density..),bins = 10,color="black",fill="white")+
  labs(x="log(Annual precipitation) centered",y="Density")+
  theme(text = element_text(size = 20))

X<-vector("list",2)
names(X)<-c("x","pp")

X$x<-split(avgtemp,row(avgtemp))
X$pp<-split(rep(subs,each=35),row(avgtemp))

#irregular
# sparseObs<-sample(16:24,35,replace=T)
# X<-fdapace::Sparsify(avgtemp,subs,sparseObs)
# names(X)<-c("pp","x")
# 

# save(CanadianWeather,X,y,file="canWeather.Rdata")

load("canWeather.Rdata")
# 
df<-data.frame(x=numeric(0),
               y=numeric(0),
               obs=numeric(0),
               color=numeric(0))

a<-c(35,26,27,29)
for (i in 1:35){
  df<-rbind(df,cbind(x=X$pp[[i]],
                     y=X$x[[i]],
                     obs=i,
                     color=(i%in%a)+1))
}

pcurves<-ggplot(df,aes(x=x,y=y,group=obs,color=color))+
  geom_line()+
  geom_point()+
  theme(legend.position="none")+
  labs(x="Day",y="Temperature")+
  theme(text = element_text(size = 20))
# 
pdf("canweather.pdf",width=10)

gridExtra::grid.arrange(pcurves,phist,nrow=1)
dev.off()

# X$x[[28]]<-X$x[[28]][-(11:16)]
# X$pp[[28]]<-X$pp[[28]][-(11:16)]
# 
# X$x[[29]]<-X$x[[29]][-(11:16)]
# X$pp[[29]]<-X$pp[[29]][-(11:16)]
# 
# X$x[[25]]<-X$x[[25]][-(11:16)]
# X$pp[[25]]<-X$pp[[25]][-(11:16)]
  

trainInd<-1:35#sample(1:35,0.8*35)#

testInd<-(1:35)[-trainInd]
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

plot(resFPCR$grid,resFPCR$betahat[[2]],type="l")
qqnorm(ytest-resFPCR$yhatTest[[2]])
qqline(ytest-resFPCR$yhatTest[[2]])
plot(resFPCR$yhatTest[[2]]~ytest)
abline(a=0,b=1)

resRFPCR<-doRFPCReg(Xtrain,ytrain,Xtest,ytest,cl=subs,K=K,method="Irreg")

plot(resRFPCR$grid,resRFPCR$betahat[[2]],type="l")
qqnorm(ytest-resRFPCR$yhatTest[[2]])
qqline(ytest-resRFPCR$yhatTest[[2]])
plot(resRFPCR$yhatTest[[2]]~ytest)
abline(a=0,b=1)

hs<-c(15,60,180)

resSLS<-doSparseFPCAReg(Xtrain,ytrain,Xtest,ytest,K=K,hs.mu = hs,hs.cov = hs,ncov=50,method="LS",preSmooth = F)

plot(resSLS$grid,resSLS$betahat[[1]],type="l")
qqnorm(ytest-resSLS$yhatTest[[1]])
qqline(ytest-resSLS$yhatTest[[1]])
plot(resSLS$yhatTest[[1]]~ytest)
abline(a=0,b=1)


plot(resSLS$grid,resSLS$betahat[[2]],type="l")
qqnorm(ytest-resSLS$yhatTest[[2]])
qqline(ytest-resSLS$yhatTest[[2]])
plot(resSLS$yhatTest[[2]]~ytest)
abline(a=0,b=1)


resSRob<-doSparseFPCAReg(Xtrain,ytrain,Xtest,ytest,K=K,hs.mu = hs,hs.cov = hs,ncov=50,method="Rob",preSmooth=F)

plot(resSRob$grid,resSRob$betahat[[1]],type="l")
qqnorm(ytest-resSRob$yhatTest[[1]])
qqline(ytest-resSRob$yhatTest[[1]])
plot(resSRob$yhatTest[[1]]~ytest)
abline(a=0,b=1)

plot(resSRob$grid,resSRob$betahat[[2]],type="l")
qqnorm(ytest-resSRob$yhatTest[[2]])
qqline(ytest-resSRob$yhatTest[[2]])
plot(resSRob$yhatTest[[2]]~ytest)
abline(a=0,b=1)

##ggplots

df<-data.frame(grid=resFPCR$grid,
               betahat=resFPCR$betahat[[2]],
               Method="FPCPR",
               type="(R)FPCPR")

df<-rbind(df,data.frame(grid=resRFPCR$grid,
                        betahat=resRFPCR$betahat[[2]],
                        Method="RFPCPR",
                        type="(R)FPCPR"))

df<-rbind(df,data.frame(grid=resSLS$grid,
                        betahat=resSLS$betahat[[1]],
                        Method="S-LS-LS",
                        type="S-Methods"))

df<-rbind(df,data.frame(grid=resSLS$grid,
                        betahat=resSLS$betahat[[2]],
                        Method="S-LS-ROB",
                        type="S-Methods"))

df<-rbind(df,data.frame(grid=resSRob$grid,
                        betahat=resSRob$betahat[[1]],
                        Method="S-ROB-LS",
                        type="S-Methods"))

df<-rbind(df,data.frame(grid=resSRob$grid,
                        betahat=resSRob$betahat[[2]],
                        Method="S-ROB-ROB",
                        type="S-Methods"))

pdf("can weather beta.pdf",width=10)
ggplot(df,aes(x=grid,y=betahat,color=Method))+
  geom_line(aes(linetype=Method))+
  facet_grid(rows=vars(type),scales="free_y")+
  labs(x="t",y=bquote(hat(beta)(t)))+
  theme(text = element_text(size = 20))

dev.off()
###################


threshold<-2.5*sd(ytest-resFPCR$yhatTest[[2]])#sort(abs(ytest-resFPCR$yhatTest[[2.5]]),decreasing = T)[3]
df<-data.frame(y=ytest,
               yhat=resFPCR$yhatTest[[2]],
               res=ytest-resFPCR$yhatTest[[2]],
               Method="FPCPR",
               type="(R)FPCPR",
               robust="non-robust",
               place=CanadianWeather$place)
df<-df%>%mutate(outlier=abs(res) >=threshold)


threshold<-2.5*sd(ytest-resRFPCR$yhatTest[[2]])#sort(abs(ytest-resRFPCR$yhatTest[[2.5]]),decreasing = T)[3]
df<-rbind(df,data.frame(y=ytest,
                        yhat=resRFPCR$yhatTest[[2]],
                        res=ytest-resRFPCR$yhatTest[[2]],
                        Method="RFPCPR",
                        type="(R)FPCPR",
                        robust="robust",
                        outlier=abs(ytest-resRFPCR$yhatTest[[2]])>=threshold,
                        place=CanadianWeather$place))

threshold<-2.5*sd(ytest-resSLS$yhatTest[[1]])#sort(abs(ytest-resSLS$yhatTest[[1]]),decreasing = T)[3]
df<-rbind(df,data.frame(y=ytest,
                        yhat=resSLS$yhatTest[[1]],
                        res=ytest-resSLS$yhatTest[[1]],
                        Method="S-LS-LS",
                        type="S-LS-",
                        robust="non-robust",
                        outlier=abs(ytest-resSLS$yhatTest[[1]])>=threshold,
                        place=CanadianWeather$place))

threshold<-2.5*sd(ytest-resSLS$yhatTest[[2]])#sort(abs(ytest-resSLS$yhatTest[[2.5]]),decreasing = T)[3]
df<-rbind(df,data.frame(y=ytest,
                        yhat=resSLS$yhatTest[[2]],
                        res=ytest-resSLS$yhatTest[[2]],
                        Method="S-LS-ROB",
                        type="S-LS-",
                        robust="robust",
                        outlier=abs(ytest-resSLS$yhatTest[[2]])>=threshold,
                        place=CanadianWeather$place))

threshold<-2.5*sd(ytest-resSRob$yhatTest[[1]])#sort(abs(ytest-resSRob$yhatTest[[1]]),decreasing = T)[3]
df<-rbind(df,data.frame(y=ytest,
                        yhat=resSRob$yhatTest[[1]],
                        res=ytest-resSRob$yhatTest[[1]],
                        Method="S-ROB-LS",
                        type="S-ROB-",
                        robust="non-robust",
                        outlier=abs(ytest-resSRob$yhatTest[[1]])>=threshold,
                        place=CanadianWeather$place))

threshold<-2.5*sd(ytest-resSRob$yhatTest[[2]])#sort(abs(ytest-resSRob$yhatTest[[2.5]]),decreasing = T)[3]
df<-rbind(df,data.frame(y=ytest,
                        yhat=resSRob$yhatTest[[2]],
                        res=ytest-resSRob$yhatTest[[2]],
                        Method="S-ROB-ROB",
                        type="S-ROB-",
                        robust="robust",
                        outlier=abs(ytest-resSRob$yhatTest[[2]])>=threshold,
                        place=CanadianWeather$place))

df<-df%>%group_by(Method)%>%arrange((res),.by_group=T)

# pdf("can weather fitted.pdf",width=10)
ggplot(df,aes(x=y,y=yhat,color=Method))+
  geom_point()+
  facet_grid(robust~type,scales="free")+
  geom_abline(intercept=0,slope=1)+
  labs(x="y",y=bquote(hat(y)))+
  theme(text = element_text(size = 20))
#   
# dev.off()
  #geom_text(label=ifelse(df$outlier==1,as.character(df$place),""),hjust=0,vjust=0)

# df<-df%>%group_by(Method)%>%mutate(sd=sd(res))
# ggplot(df,aes(x=place,y=res,color=Method))+
#   geom_point()+
#   facet_grid(robust~type,scales="free")+
#   geom_hline(aes(yintercept=sd*c(-2.5,2.5)))+
#   geom_text(label=ifelse(abs(df$res)>2.5*df$sd,as.character(df$place),""),hjust=0,vjust=0)

tmp<-ifelse(df$Method=="FPCPR" & df$place=="Quebec",1.5,
            ifelse(df$Method=="FPCPR" & df$place=="Scheffervll",1,0))

pdf("can weather qq.pdf",width=10)
ggplot(df,aes(sample=res,color=Method))+
  geom_point(stat="qq")+#stat_qq()+
  stat_qq_line(fullrange=T)+
  facet_grid(robust~type,scales="free")+
  geom_text(label=ifelse(df$outlier==1,as.character(df$place),""),
            hjust=ifelse(df$Method=="FPCPR" & df$place=="Scheffervll",0,1),
            vjust=ifelse(df$Method=="FPCPR" & df$place=="Scheffervll",0.8,0),stat="qq",
            size=4,)+
  scale_x_continuous(limits = c(-3,3))+
  theme(text = element_text(size = 20))+
  labs(x="Theoretical Quantiles",y="Sample Quantiles")
dev.off()

#############

cverr<-data.frame(err=numeric(0),
                  Method=character(0),
                  k=numeric(0))

cverr<-rbind(cverr,data.frame(err=resFPCR$cvErr,
                              Method="FPCPR",
                              k=1:K))
cverr<-rbind(cverr,data.frame(err=resRFPCR$cvErr,
                              Method="RFPCPR",
                              k=1:K))
cverr<-rbind(cverr,data.frame(err=resSLS$cvErr[[1]],
                              Method="S-LS-LS",
                              k=1:K))
cverr<-rbind(cverr,data.frame(err=resSLS$cvErr[[2]],
                              Method="S-LS-ROB",
                              k=1:K))
cverr<-rbind(cverr,data.frame(err=resSRob$cvErr[[1]],
                              Method="S-ROB-LS",
                              k=1:K))
cverr<-rbind(cverr,data.frame(err=resSRob$cvErr[[2]],
                              Method="S-ROB-ROB",
                              k=1:K))

pdf("can weather cv.pdf", width=10)
ggplot(cverr,aes(y=err,x=k,group=Method,color=Method,linetype=Method))+
  geom_point()+
  geom_line()+
  labs(x="Number of Components",y="LOOCV Error")+
  theme(text = element_text(size = 20))
dev.off()
