#source("rfpcr_vanaelst/PP.R") #consult Kalogridis or Van Aelst for corresponding file
library(sparseFPCA) #to install: devtools::install_github('msalibian/sparseFPCA', ref = "master")
library(fdapace)

doSparseFPCAReg <-
  function(Xtrain,
           ytrain,
           Xtest,
           ytest,
           hs.mu,
           hs.cov,
           K,
           ncov,
           method,
           methodReg="robust",
           preSmooth = FALSE,
           fixK=FALSE,
           rho=1e-3,
           ...) {
    
    #presmooth
    if (preSmooth) {
      Xs <- Xtrain
      for (i in 1:length(X$x)) {
        ys <- Xtrain$x[[i]]
        xs <- Xtrain$pp[[i]]
        if (length(xs) >= 4) {
          Xs$x[[i]] <- smooth.spline(xs, ys, spar = NULL, all.knots = T)$y
        }
      }
      
      Xtrain <- Xs
    }
    tryCatch({
      if (method == "LS") {
        sfpca <- sparseFPCA::lsfpca(
          X = Xtrain,
          ncpus = 4,
          hs.mu = hs.mu,
          hs.cov = hs.cov,
          k = K,
          s = K,
          trace = FALSE,
          seed = 01327001,
          k.cv = 5,
          ncov = ncov,
          rho.param = rho
        )
      }
      else if (method == "Rob") {
        sfpca <- sparseFPCA::efpca(
          X = Xtrain,
          ncpus = 1,
          hs.mu = hs.mu,
          hs.cov = hs.cov,
          alpha = 0.2,
          k = K,
          s = K,
          trace = F,
          seed = 01327001,
          k.cv = 5,
          ncov = ncov,
          rho.param = rho
        )
      }
      
      sfpca.eigen <- eigen(sfpca$cov.fun)
      sfpca.eigen$values <- sfpca.eigen$values[1:K]
      sfpca.eigen$values[sfpca.eigen$values < 0] <- 0
      sfpca.eigen$vectors <- sfpca.eigen$vectors[, 1:K]
      
      nonZeroEV <- TRUE
      lambda <- sfpca.eigen$values[nonZeroEV] / (ncov - 1)
      
      sfpca.eigen$values <- sfpca.eigen$values[nonZeroEV]
      sfpca.eigen$vectors <- sfpca.eigen$vectors[, nonZeroEV]
      
      ct <-
        lmrob.control(
          tuning.psi = 4.685,
          k.max = 5000, max.it = 5000, maxit.scale = 5000, nResample = 2000
        )
      
      comp<-K
      comp.r<-K
      if (fixK==FALSE){
        Mf.r <- rep(NA, min(K, ncol(sfpca.eigen$vectors)))
        Mf<-Mf.r
        for (f in 1:length(Mf)) {
              fit.r1 <- lmrob(ytrain ~ sfpca$xis[, nonZeroEV][, 1:f], control = ct)
              hat.matrix <-
                cbind(rep(1, length(ytrain)), sfpca$xis[, nonZeroEV][, 1:f]) %*% solve(t(cbind(rep(
                  1, length(ytrain)
                ) , sfpca$xis[, nonZeroEV][, 1:f])) %*% diag(fit.r1$rweights) %*% cbind(rep(1, length(ytrain)) , sfpca$xis[, nonZeroEV][, 1:f])) %*%
                t(cbind(rep(1, length(ytrain)) , sfpca$xis[, nonZeroEV][, 1:f])) %*% diag(fit.r1$rweights)
              press.res.r <- fit.r1$residuals / (1 - diag(hat.matrix))
              Mf.r[f] <- scaleTau2(press.res.r, c2 = 5) ^ 2

              fit.1 <- lm(ytrain ~ sfpca$xis[,nonZeroEV][,1:f])
              press.res <- fit.1$residuals / (1-hatvalues(fit.1))
              Mf[f] <- mean(press.res^2)
            
        }

        se<-sd(Mf.r)/sqrt(length(Mf.r))
        comp.r <- min(which(abs(Mf.r-min(Mf.r))<se))
        
        se<-sd(Mf)/sqrt(length(Mf))
        comp <- min(which(abs(Mf-min(Mf))<se))
         
      }

      fit.r<-lmrob(ytrain ~ sfpca$xis[, nonZeroEV][, 1:comp.r], control = ct)
      
      fit<-lm(ytrain ~ sfpca$xis[, nonZeroEV][, 1:comp])
      

      norms <- apply(sfpca.eigen$vectors,
                     2, L2.norma.mesh, mesh = sfpca$tt)
      sfpca.ef <- scale(sfpca.eigen$vectors
                        , center = FALSE, scale = norms)

      
      betahat.r <- as.numeric(sfpca.ef[, 1:comp.r] %*% as.matrix(fit.r$coef[-1]))
      betahat <- as.numeric(sfpca.ef[, 1:comp] %*% as.matrix(fit$coef[-1]))
      
      yhatTrain.r<-fit.r$fitted.values
      yhatTrain<-fit$fitted.values
      
      
      tts <- unlist(Xtrain$pp)
      mus <- unlist(sfpca$muh)
      mu.fn <- approxfun(x=tts, y=mus)
      mu.fn.r <- lapply(Xtest$pp,mu.fn)

      if (comp.r==comp){
        testScores<-sparseFPCA::pred.scores(Xtrain,sfpca$muh,Xtest,mu.fn.r,sfpca$cov.fun,sfpca$tt,
                                           k=comp,s=comp,rho=rho)
        testScores.r<-testScores
      }
      else{
        testScores.r<-sparseFPCA::pred.scores(Xtrain,sfpca$muh,Xtest,mu.fn.r,sfpca$cov.fun,sfpca$tt,
                                            k=comp.r,s=comp.r,rho=rho)
        testScores<-sparseFPCA::pred.scores(Xtrain,sfpca$muh,Xtest,mu.fn.r,sfpca$cov.fun,sfpca$tt,
                                            k=comp,s=comp,rho=rho)
      }

       yhatTest.r<-as.numeric(testScores.r%*%fit.r$coef[-1]+fit.r$coef[1])
       yhatTest<-as.numeric(testScores%*%fit$coef[-1]+fit$coef[1])#predict(fit,newdata=data.frame(testScores))
      
      res <- list(betahat=list(betahat,betahat.r),
                  grid=sfpca$tt,
                  yhatTest=list(yhatTest,yhatTest.r),
                  comp=list(comp,comp.r),
                  cvErr=list(Mf,Mf.r))
      return(res)
    },
    error = function(cond) {
      print(cond)      
      return(NULL)
    })
  }


doRFPCReg <- function(Xtrain, ytrain, Xtest, ytest, K, cl, method="Reg",fixK=FALSE,...) {
  tryCatch({
    if (method=="Reg"){
      Xmat <- matrix(unlist(Xtrain$x),
                     nrow = length(Xtrain$x),
                     byrow = TRUE)
      XmatTest <- matrix(unlist(Xtest$x),
                     nrow = length(Xtest$x),
                     byrow = TRUE)
      
      rfpcrReg <- flm.rob(X = Xmat, y = ytrain,cl,Xtest=XmatTest,ytest=ytest,
                          k = K,method=method,fixK=fixK,...)
    }
    else if (method=="Irreg"){
      tmp<-imputeIrregData(Xtrain)
      XImp<-tmp$X
      cl<-tmp$cl
      Xmat <- matrix(unlist(XImp$x),
                     nrow = length(XImp$x),
                     byrow = TRUE)
      
      rfpcrReg <- flm.rob(X = Xmat, y = ytrain,cl,Xtest=Xtest,ytest=ytest,
                          k = K,method=method,fixK=fixK)
    }
    
  
    yhatTrain<-rfpcrReg$yhatTrain
    yhatTest<-rfpcrReg$yhatTest
    
    yhatTrain1<-rfpcrReg$yhatTrain1
    yhatTest1<-rfpcrReg$yhatTest1
    
    res<-list(betahat=list(rfpcrReg$beta,rfpcrReg$beta1),
              grid=cl,
              yhatTest=list(yhatTest,yhatTest1),
              comp=rfpcrReg$comp.ret,
              cvErr=rfpcrReg$cvErr)
    
    return(res)
  },
  error = function(cond) {
    print(cond)
    return(NULL)
  })
}

doFPCReg <- function(Xtrain, ytrain, Xtest, ytest, K, cl, method="Reg",fixK=FALSE,...) {
  tryCatch({
    if (method=="Reg"){
      Xmat <- matrix(unlist(Xtrain$x),
                     nrow = length(Xtrain$x),
                     byrow = TRUE)
      XmatTest <- matrix(unlist(Xtest$x),
                         nrow = length(Xtest$x),
                         byrow = TRUE)
      
      rfpcrReg <- flm.nonrob(X = Xmat, y = ytrain,cl,Xtest=XmatTest,ytest=ytest,
                          k = K,method=method,fixK=fixK,...)
    }
    else if (method=="Irreg"){
      tmp<-imputeIrregData(Xtrain)
      XImp<-tmp$X
      cl<-tmp$cl
      Xmat <- matrix(unlist(XImp$x),
                     nrow = length(XImp$x),
                     byrow = TRUE)
      
      rfpcrReg <- flm.nonrob(X = Xmat, y = ytrain,cl,Xtest=Xtest,ytest=ytest,
                          k = K,method=method,fixK=fixK)
    }
   
    yhatTrain<-rfpcrReg$yhatTrain
    yhatTest<-rfpcrReg$yhatTest
    
    yhatTrain1<-rfpcrReg$yhatTrain1
    yhatTest1<-rfpcrReg$yhatTest1
    
    res<-list(betahat=list(rfpcrReg$beta,rfpcrReg$beta1),
              grid=cl,
              yhatTest=list(yhatTest,yhatTest1),
              comp=rfpcrReg$comp.ret,
              cvErr=rfpcrReg$cvErr)
    
    return(res)
  },
  error = function(cond) {
    print(cond)    
    return(NULL)
  })
}

imputeIrregData <- function(X) {
  pps <- sort(unique(unlist(X$pp)))
  XImp <- vector("list", 2)
  names(XImp) <- c("x", "pp")
  N<-length(X$x)
  XImp$pp <- replicate(N, pps, simplify = FALSE)
  
  
  for (i in 1:N) {
    tmp <- approxfun(X$pp[[i]], X$x[[i]])
    xtemp <- numeric(length(pps))
    for (j in 1:length(pps)) {
      if (max(X$pp[[i]]) >= pps[j] & pps[j] >= min(X$pp[[i]])) {
        xtemp[j] <- tmp(pps[j])
      }
      else if (pps[j] < min(X$pp[[i]])) {
        xtemp[j] <- tmp(min(X$pp[[i]]))
      }
      else{
        xtemp[j] <- tmp(max(X$pp[[i]]))
      }
    }
    XImp$x[[i]] <- xtemp#smooth.spline(pps, xtemp, all.knots = T, spar = NULL)$y
  }
  
  return(list(X=XImp,cl=pps))
  #return(XImp)
  
}

doSimulation<-function(N,n,K,method,sparseObs=3:5,fixK=FALSE,betaCoef=NULL,eps=NULL,contFactorX=5,contFactorY=10,sigma=0.1,...){
  ncov <- 50 
  k.cv <- 5 
  hs.mu <- c(0.1,0.2) 
  hs.cov <- hs.mu
  kmax<-50
  cl<-seq(0,1,length.out=n)

  #data
  if (method=="Reg"){
    XReg<-Wiener(2*N,pts=cl,K=kmax,sparsify = rep(n,N))
    names(XReg)<-c("pp","x")
    errors<-rnorm(2*N)
    y0<-sapply(1:(2*N),function(i){
      mean(XReg$x[[i]]*betaCoef(XReg$pp[[i]]))
    })
    yReg<-sapply(1:(2*N),function(i){
      y0[i]+sd(y0)*sigma*errors[i]
    })
    
    XCont<-XReg
    XCont$x<-XReg$x[1:N]
    XCont$pp<-XReg$pp[1:N]
    yCont<-yReg[1:N]
    
    Xtest<-XReg
    Xtest$x<-XReg$x[-(1:N)]
    Xtest$pp<-XReg$pp[-(1:N)]
    y0test<-y0[-(1:N)]
  }
  else if (method=="Irreg"){
    XReg<-Wiener(2*N,pts=cl,K=kmax)
    sparseSample<-sample(sparseObs,2*N,replace=T)
    XIrreg<-Sparsify(XReg,cl,sparsity = sparseSample,...)
    names(XIrreg)<-c("pp","x")
    errors<-rnorm(2*N)
    y0<-as.numeric(XReg%*%betaCoef(cl)/ncol(XReg))
    yIrreg<-sapply(1:(2*N),function(i){
      y0[i]+sigma*errors[i]
    })
    
    XCont<-XIrreg
    XCont$x<-XIrreg$x[1:N]
    XCont$pp<-XIrreg$pp[1:N]
    yCont<-yIrreg[1:N]
    
    Xtest<-XIrreg
    Xtest$x<-XIrreg$x[-(1:N)]
    Xtest$pp<-XIrreg$pp[-(1:N)]
    
    y0test<-y0[-(1:N)]
  }
  
  #contamination
  if (!is.null(eps)){
    contInd<-1:(N*eps)
    for (i in contInd){
      XCont$x[[i]]<-XCont$x[[i]]*contFactorX
      yCont[i]<-contFactorY*y0[i]+sigma*errors[i]
    }
  }
  

  resErr<-data.frame(estError=numeric(0),
                     predErrorTest=numeric(0),
                     comp=numeric(0))
  resSim<-list()
  

  resFPCR<-doFPCReg(XCont,yCont,Xtest,y0test,K = K,cl=cl,method=method,fixK=fixK)
  
  if (!is.null(resFPCR)){
    resErr<-rbind(resErr,
                  data.frame(estError=mean( (resFPCR$betahat[[2]]-betaCoef(resFPCR$grid) )^2/
                                              mean(betaCoef(resFPCR$grid)^2)),
                             predError=mean( (resFPCR$yhatTest[[2]]-y0test)^2/
                                               y0test^2),
                             comp=resFPCR$comp))
  }
  else{
    resErr<-rbind(resErr,data.frame(estError=c(NA,NA),
                                    predError=c(NA,NA),
                                    comp=c(NA,NA)))
  }
  
  resSim<-c(resSim,list(fpcr=resFPCR))
  
  resRFPCR<-doRFPCReg(XCont,yCont,Xtest,y0test,K = K,cl=cl,method=method,fixK=fixK)
  if (!is.null(resRFPCR)){
    resErr<-rbind(resErr,
                  data.frame(estError=mean( (resRFPCR$betahat[[2]]-betaCoef(resRFPCR$grid) )^2/
                                              mean(betaCoef(resRFPCR$grid)^2)),
                             predError=mean( (resRFPCR$yhatTest[[2]]-y0test)^2/
                                               y0test^2),
                             comp=resRFPCR$comp))
  }
  else{
    resErr<-rbind(resErr,data.frame(estError=c(NA,NA),
                                    predError=c(NA,NA),
                                    comp=c(NA,NA)))
  }
  
  resSim<-c(resSim,list(rfpcr=resRFPCR))
  
  
  #sfpca
  resSFPCA<-doSparseFPCAReg(XCont,yCont,Xtest,y0test,hs.mu=hs.mu,hs.cov=hs.cov,K=K,ncov=ncov,method="LS",fixK=fixK)
  if (!is.null(resSFPCA)){
    resErr<-rbind(resErr,
                  data.frame(estError=mean( (resSFPCA$betahat[[1]]-betaCoef(resSFPCA$grid) )^2/
                                              mean(betaCoef(resSFPCA$grid)^2)),
                             predError=mean( (resSFPCA$yhatTest[[1]]-y0test)^2/
                                               y0test^2),
                             comp=resSFPCA$comp[[1]]),
                  data.frame(estError=mean( (resSFPCA$betahat[[2]]-betaCoef(resSFPCA$grid) )^2/
                                              mean(betaCoef(resSFPCA$grid)^2)),
                             predError=mean( (resSFPCA$yhatTest[[2]]-y0test)^2/
                                               y0test^2),
                             comp=resSFPCA$comp[[2]]))
  }
  else{
    resErr<-rbind(resErr,data.frame(estError=c(NA,NA),
                                    predError=c(NA,NA),
                                    comp=c(NA,NA)))
  }
  resSim<-c(resSim,list(sfpcaLS=resSFPCA))
  
  resSFPCArob<-doSparseFPCAReg(XCont,yCont,Xtest,y0test,hs.mu=hs.mu,hs.cov=hs.cov,K=K,ncov=ncov,method="Rob",fixK=fixK)
  
  if (!is.null(resSFPCArob)){
  resErr<-rbind(resErr,
                data.frame(estError=mean( (resSFPCArob$betahat[[1]]-betaCoef(resSFPCArob$grid) )^2/
                                            mean(betaCoef(resSFPCArob$grid)^2)),
                           predError=mean( (resSFPCArob$yhatTest[[1]]-y0test)^2/
                                             y0test^2),
                           comp=resSFPCArob$comp[[1]]),
                data.frame(estError=mean( (resSFPCArob$betahat[[2]]-betaCoef(resSFPCArob$grid) )^2/
                                            mean(betaCoef(resSFPCArob$grid)^2)),
                           predError=mean( (resSFPCArob$yhatTest[[2]]-y0test)^2/
                                             y0test^2),
                           comp=resSFPCArob$comp[[2]]))
  }
  else{
    resErr<-rbind(resErr,data.frame(estError=c(NA,NA),
                                    predError=c(NA,NA),
                                    comp=c(NA,NA)))
  }
  resSim<-c(resSim,list(sfpcaRobrob=resSFPCArob))
  
  resSim<-c(resSim,list(Xtrain=XCont,ytrain=yCont,Xtest=Xtest,ytest=y0test))
  
  rownames(resErr)<-c("FPCPR","RFPCPR","SLSLS","SLSROB","SROBLS","SROBROB")
  
  return(list(err=resErr,data=resSim))
}


