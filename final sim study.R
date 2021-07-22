library(parallel)

n.cores<-detectCores()
cl<-makeCluster(10) #check number of cores
clusterEvalQ(cl,{
  source("./functions.R")
  
  N<-100 #number of curves
  
  sigma<-0.1 #noise standard deviation
  
  n<-100 #number of observations per curve
  
  #betaCoef<-function(t)log(1.5*t^2+10)+cos(4*pi*t) #does notwork for sigma 0.5
  betaCoef<-function(t)2*sin(0.5*pi*t)+4*sin(1.5*pi*t)+5*sin(2.5*pi*t)
  
  K<-8
  
  fixK<-F
  
  sparseAgg<-T
  sparseObs<-3:5
  
})
clusterSetRNGStream(cl,01327001)


m<-30

##REGULAR
system.time({
  res1<-parSapply(cl,1:m,function(mi){#replicate(m,{
    doSimulation(N,n,K,method = "Reg",betaCoef = betaCoef,fixK=fixK)
  },simplify=F)
})

system.time({
  res2<-parSapply(cl,1:m,function(mi){#replicate(m,{
    doSimulation(N,n,K,method = "Reg",betaCoef = betaCoef,fixK=fixK,eps=0.1)
  },simplify=F)
})

system.time({
  res3<-parSapply(cl,1:m,function(mi){#replicate(m,{
    doSimulation(N,n,K,method = "Reg",betaCoef = betaCoef,fixK=fixK,eps=0.2)
  },simplify=F)
})

stopCluster(cl)

save(res1,res2,res3,
     file="results2/sigma01/K8CVn100Reg.Rdata")

# #IRREGULAR
# system.time({
#   res1<-parSapply(cl,1:m,function(mi){#replicate(m,
#     doSimulation(N,n,K,method = "Irreg",betaCoef = betaCoef,fixK=fixK,sparseObs=sparseObs,aggressive=sparseAgg,
#                  contFactor=contFactor)
#   },simplify = F)
# })
# 
# system.time({
#   res2<-parSapply(cl,1:m,function(mi){#replicate(m,
#     doSimulation(N,n,K,method = "Irreg",betaCoef = betaCoef,fixK=fixK,sparseObs=sparseObs,eps=0.1,aggressive=sparseAgg,
#                  contFactor=contFactor)
#   },simplify = F)
# })
# 
# system.time({
#   res3<-parSapply(cl,1:m,function(mi){#replicate(m,
#     doSimulation(N,n,K,method = "Irreg",betaCoef = betaCoef,fixK=fixK,sparseObs=sparseObs,eps=0.2,aggressive=sparseAgg,
#                  contFactor=contFactor)
#   },simplify = F)
# })
# 
# stopCluster(cl)
# 
# save(res1,res2,res3,
#      file="results/sigma01/K8CVn35Cont10AggIrreg.Rdata")

