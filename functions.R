source("rfpcr_vanaelst/PP.R")
library(sparseFPCA)
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
      
      #print(sfpca$muh)
      
      sfpca.eigen <- eigen(sfpca$cov.fun)
      sfpca.eigen$values <- sfpca.eigen$values[1:K]
      sfpca.eigen$values[sfpca.eigen$values < 0] <- 0
      sfpca.eigen$vectors <- sfpca.eigen$vectors[, 1:K]
      
      nonZeroEV <- TRUE#sfpca.eigen$values>1e-5
      lambda <- sfpca.eigen$values[nonZeroEV] / (ncov - 1)
      
      #comp.exp<-min(which(cumsum(lambda)/sum(lambda)>=0.99))
      
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
        # comp.fit <-
        #   which.min(Mf) # Cross-validation for the selection of K
        
        se<-sd(Mf.r)/sqrt(length(Mf.r))
        comp.r <- min(which(abs(Mf.r-min(Mf.r))<se))#min(comp.fit,comp.exp)
        
        se<-sd(Mf)/sqrt(length(Mf))
        comp <- min(which(abs(Mf-min(Mf))<se))
        
        
      }

      fit.r<-lmrob(ytrain ~ sfpca$xis[, nonZeroEV][, 1:comp.r], control = ct)
      
      fit<-lm(ytrain ~ sfpca$xis[, nonZeroEV][, 1:comp])
      

      norms <- apply(sfpca.eigen$vectors,
                     2, L2.norma.mesh, mesh = sfpca$tt)
      sfpca.ef <- scale(sfpca.eigen$vectors
                        , center = FALSE, scale = norms)
      
      #matplot(sfpca$tt,sfpca.ef[,1:2],...)
      
      betahat.r <- as.numeric(sfpca.ef[, 1:comp.r] %*% as.matrix(fit.r$coef[-1]))
      betahat <- as.numeric(sfpca.ef[, 1:comp] %*% as.matrix(fit$coef[-1]))
      #lines(sfpca$tt,betahat,lty=3+2*(method=="Rob"),col=3+2*(method=="Rob"))
      #lines(sfpca$tt,betahat.r,lty=4+2*(method=="Rob"),col=4+2*(method=="Rob"))
      
      yhatTrain.r<-fit.r$fitted.values
      yhatTrain<-fit$fitted.values
      
      
      tts <- unlist(Xtrain$pp)
      mus <- unlist(sfpca$muh)
      mu.fn <- approxfun(x=tts, y=mus)
      mu.fn.r <- lapply(Xtest$pp,mu.fn)
      #browser()
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
      #matplot(testScores,...)
      # 
       yhatTest.r<-as.numeric(testScores.r%*%fit.r$coef[-1]+fit.r$coef[1])
       yhatTest<-as.numeric(testScores%*%fit$coef[-1]+fit$coef[1])#predict(fit,newdata=data.frame(testScores))
      
      # betahatfn<-splinefun(sfpca$tt,betahat)
      # #betahatfn <- approxfun(sfpca$tt, betahat)
      # # yhatTrain <-
      # #   sapply(1:length(Xtrain$x), function(i) {
      # #     #predict using smoothed or unsmoothed X?
      # #     mean(Xtrain$x[[i]] * betahatfn(Xtrain$pp[[i]]))
      # #   })
      # yhatTest <- sapply(1:length(Xtest$x), function(i) {
      #   mean(Xtest$x[[i]] * betahatfn(Xtest$pp[[i]]))
      # })
      
      # estError<-ifelse(!is.null(betaCoef),integrate( function(x) (betaCoef(x)-betahatfn(x))^2,
      #                                            lower=0,upper=1)$value,NA)
      # estError <-
      #   ifelse(!is.null(betaCoef), mean((betahat - betaCoef(sfpca$tt)) ^ 2), NA)
      # 
      # res <-  data.frame(
      #     estError = estError,
      #     predErrorTrain = mean((ytrain - yhatTrain) ^
      #                             2),
      #     predErrorTest = mean((ytest - yhatTest) ^
      #                            2),
      #     comp.s = comp.s
      #   )
      
      res <- list(betahat=list(betahat,betahat.r),
                  grid=sfpca$tt,
                  yhatTest=list(yhatTest,yhatTest.r),
                  comp=list(comp,comp.r),
                  cvErr=list(Mf,Mf.r))
      return(res)
    },
    error = function(cond) {
      print(cond)
      # res <-   data.frame(
      #                estError = NA,
      #                predErrorTrain = NA,
      #                predErrorTest = NA,
      #                comp.s = NA
      #              )
      
      return(NULL)
    })
    #return(res)
    #print(res)
    #return(res[which.min(res$predErrorTest),])
  }

doPaceReg <-
  function(Xtrain,
           ytrain,
           Xtest,
           ytest,
           K,
           ncov,
           preSmooth=FALSE,
           fixK=FALSE,
           ...) {
    
    myop<-list(error=T,
               verbose=F,
               nRegGrid=ncov,
               maxK=K,
               methodXi="CE")
    
    #presmooth
    if (preSmooth){
      Xs <- Xtrain
      for (i in 1:N) {
        ys <- Xtrain$x[[i]]
        xs <- Xtrain$pp[[i]]
        if (length(xs) >= 4) {
          Xs$x[[i]] <- smooth.spline(xs, ys, spar = NULL, all.knots = T)$y
        }
      }
      
      Xtrain <- Xs
    }
    
    tryCatch({
      paceIrreg <- FPCA(Ly = Xtrain$x,
                        Lt = Xtrain$pp,
                        optns = myop)
      
      #print(paceIrreg$optns)
      #comp.exp<-min(which(cumsum(paceIrreg$lambda)/sum(paceIreg$lambda)>=0.99))
      
      ct <-
        lmrob.control(
          tuning.psi = 4.685,
          k.max = 5000, max.it = 5000, maxit.scale = 5000, nResample = 2000
        )
      comp.s<-K
      if (fixK==FALSE){
        Mf <- rep(NA, min(K, ncol(paceIrreg$phi)))
        for (f in 1:length(Mf)) {
          fit.r1 <- lmrob(ytrain ~ paceIrreg$xiEst[, 1:f], control = ct)
          hat.matrix <-
            cbind(rep(1, length(ytrain)), paceIrreg$xiEst[, 1:f]) %*% solve(t(cbind(rep(
              1, length(ytrain)
            ) , paceIrreg$xiEst[, 1:f])) %*% diag(fit.r1$rweights) %*% cbind(rep(1, length(ytrain)) , paceIrreg$xiEst[, 1:f])) %*%
            t(cbind(rep(1, length(ytrain)) , paceIrreg$xiEst[, 1:f])) %*% diag(fit.r1$rweights)
          press.res <- fit.r1$residuals / (1 - diag(hat.matrix))
          #Mf[f] <- CVcrit(press.res,...)^2
          Mf[f] <- scaleTau2(press.res, c2 = 5) ^ 2
        }
        comp.fit <-
          which.min(Mf) # Cross-validation for the selection of K
        
        #comp.fit<-(sum(abs(diff(Mf))>1e-2)+1)
        #comp.fit<-which.max(Mf[3:length(Mf)]+Mf[1:(length(Mf)-2)]-2*Mf[2:(length(Mf)-1)])+1
        #plot(Mf[3:length(Mf)]+Mf[1:(length(Mf)-2)]-2*Mf[2:(length(Mf)-1)])
        #plot(Mf)
        
        #Mf2max<-which.max(Mf[3:length(Mf)]+Mf[1:(length(Mf)-2)]-2*Mf[2:(length(Mf)-1)])+1
        se<-sd(Mf)/sqrt(length(Mf))
        Mfmin<-min(which(abs(Mf-min(Mf))<se))
        #plot(Mf)
        
        #print(which(abs(Mf-min(Mf))<se))
        
        # if (abs(Mf2max-Mfmin)<1e-3)comp.fit<-Mf2max
        # else comp.fit<-Mfmin
        comp.fit<-Mfmin
        comp.s <- comp.fit#min(comp.exp,comp.fit)
      }
      
      fit<-lmrob(ytrain ~ paceIrreg$xiEst[, 1:comp.s], control = ct)
      coefBeta <-
        fit$coef[-1]
      
      #print(fit$coef)
      
      #matplot(paceIrreg$workGrid,paceIrreg$phi[, 1:comp.s],...)
      
      betahat <- as.numeric(paceIrreg$phi[, 1:comp.s] %*% as.matrix(coefBeta))
      #plot(paceIrreg$workGrid,betahat,...)
      # betahatfn<-splinefun(paceIrreg$workGrid,betahat)
      # #betahatfn <- approxfun(paceIrreg$workGrid, betahat)
      # yhatTrain <-
      #   sapply(1:length(Xtrain$x), function(i) {
      #     #predict using smoothed or unsmoothed X?
      #     mean(Xtrain$x[[i]] * betahatfn(Xtrain$pp[[i]]))
      #   })
      # yhatTest <- sapply(1:length(Xtest$x), function(i) {
      #   mean(Xtest$x[[i]] * betahatfn(Xtest$pp[[i]]))
      # })
      # estError<-ifelse(!is.null(betaCoef),integrate( function(x) (betaCoef(x)-betahatfn(x))^2,
      #                                            lower=0,upper=1)$value,NA)
      
      #yhatTrain<-paceIrreg$xiEst[, 1:comp.s]%*%coefBeta+fit$coef[1]#fit$fitted.values
      yhatTrain<-fit$fitted.values
      
      
      testPred<-predict(paceIrreg,Xtest$x,Xtest$pp,K=comp.s,xiMethod="CE")
      #print(testPred)
      
      testScores<-testPred$scores
      
      #str(testScores)
      #matplot(testScores,...)
      
      # betahatfn<-splinefun(paceIrreg$workGrid,betahat)
      # yhatTest <- sapply(1:length(Xtest$x), function(i) {
      #   mean(Xtest$x[[i]] * betahatfn(Xtest$pp[[i]]))
      # })
      yhatTest<-as.numeric(testScores%*%coefBeta+fit$coef[1])
      
      #plot(yhatTest,predict(fit,newdata=data.frame(testScores),type="response"))
      #yhatTest<-predict(fit,newdata=data.frame(testScores),type="response")
      #plot(ytest~yhatTest)
      
      # estError <-
      #   ifelse(!is.null(betaCoef), mean((betahat - betaCoef(
      #     paceIrreg$workGrid
      #   )) ^ 2), NA)
      # 
      # res <- data.frame(
      #     estError = estError,
      #     predErrorTrain = mean((ytrain - yhatTrain)^2),
      #     predErrorTest = mean((ytest - yhatTest)^2),
      #     comp.s = comp.s
      #   )
      
      res<-list(betahat=betahat,
                grid=paceIrreg$workGrid,
                yhatTest=yhatTest,
                comp.s=comp.s)
      
      
      return(res)
    },
    error = function(cond) {
      print(cond)
      # res <-   data.frame(
      #                estError = NA,
      #                predErrorTrain = NA,
      #                predErrorTest = NA,
      #                comp.s = NA
      #              )
      res <- list(betahat=NA,
                  grid=NA,
                  yhatTest=NA,
                  comp.s=comp.s)
      return(res)
    })
    
    #return(res)
    #return(res[which.min(res$predErrorTest),])
    
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
    
    
    # betahat<-splinefun(cl,rfpcrReg$betaCoef)
    # #betahat <- approxfun(cl, rfpcrReg$betaCoef)
    # yhatTrain <- sapply(1:length(Xtrain$x), function(i) {
    #   mean(Xtrain$x[[i]] * betahat(Xtrain$pp[[i]]))
    # })
    # yhatTest <- sapply(1:length(Xtest$x), function(i) {
    #   mean(Xtest$x[[i]] * betahat(Xtest$pp[[i]]))
    # })
    #lines(cl,rfpcrReg$beta1,lty=2,col=2)
    yhatTrain<-rfpcrReg$yhatTrain
    yhatTest<-rfpcrReg$yhatTest
    
    # estError=ifelse(!is.null(betaCoef),integrate( function(x) (betaCoef(x)-betahat(x))^2,
    #                                           lower=0,upper=1)$value,NA)
    # res <-
    #   data.frame(
    #     estError = ifelse(!is.null(betaCoef), mean((
    #       rfpcrReg$betaCoef - betaCoef(cl)
    #     ) ^ 2), NA),
    #     predErrorTrain = mean((ytrain - yhatTrain) ^ 2),
    #     predErrorTest = mean((ytest - yhatTest) ^ 2),
    #     comp.s = rfpcrReg$comp.ret
    #   )
    
    # betahat <- splinefun(cl, rfpcrReg$beta1)
    # yhatTrain <- sapply(1:length(Xtrain$x), function(i) {
    #   mean(Xtrain$x[[i]] * betahat(Xtrain$pp[[i]]))
    # })
    # yhatTest <- sapply(1:length(Xtest$x), function(i) {
    #   mean(Xtest$x[[i]] * betahat(Xtest$pp[[i]]))
    # })
    
    yhatTrain1<-rfpcrReg$yhatTrain1
    yhatTest1<-rfpcrReg$yhatTest1
    # plot(yhatTest~ytest)
    # abline(0,1)
    # estError=ifelse(!is.null(betaCoef),integrate( function(x) (betaCoef(x)-betahat(x))^2,
    #                                           lower=0,upper=1)$value,NA)
    
    # res <-
    #   rbind(
    #     res,
    #     data.frame(
    #       estError = ifelse(!is.null(betaCoef), mean((
    #         rfpcrReg$beta1 - betaCoef(cl)
    #       ) ^ 2), NA),
    #       predErrorTrain = mean((ytrain - yhatTrain) ^
    #                               2),
    #       predErrorTest = mean((ytest - yhatTest) ^ 2),
    #       comp.s = rfpcrReg$comp.ret
    #     )
    #   )
    res<-list(betahat=list(rfpcrReg$beta,rfpcrReg$beta1),
              grid=cl,
              yhatTest=list(yhatTest,yhatTest1),
              comp=rfpcrReg$comp.ret,
              cvErr=rfpcrReg$cvErr)
    
    return(res)
  },
  error = function(cond) {
    print(cond)
    # return(data.frame(
    #   estError = c(NA, NA),
    #   predErrorTrain = c(NA, NA),
    #   predErrorTest = c(NA, NA),
    #   comp.s = c(NA, NA)
    # ))
    
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
    
    
    # betahat<-splinefun(cl,rfpcrReg$betaCoef)
    # #betahat <- approxfun(cl, rfpcrReg$betaCoef)
    # yhatTrain <- sapply(1:length(Xtrain$x), function(i) {
    #   mean(Xtrain$x[[i]] * betahat(Xtrain$pp[[i]]))
    # })
    # yhatTest <- sapply(1:length(Xtest$x), function(i) {
    #   mean(Xtest$x[[i]] * betahat(Xtest$pp[[i]]))
    # })
    #plot(cl,rfpcrReg$beta1,type="l")
    yhatTrain<-rfpcrReg$yhatTrain
    yhatTest<-rfpcrReg$yhatTest
    
    # estError=ifelse(!is.null(betaCoef),integrate( function(x) (betaCoef(x)-betahat(x))^2,
    #                                           lower=0,upper=1)$value,NA)
    # res <-
    #   data.frame(
    #     estError = ifelse(!is.null(betaCoef), mean((
    #       rfpcrReg$betaCoef - betaCoef(cl)
    #     ) ^ 2), NA),
    #     predErrorTrain = mean((ytrain - yhatTrain) ^ 2),
    #     predErrorTest = mean((ytest - yhatTest) ^ 2),
    #     comp.s = rfpcrReg$comp.ret
    #   )
    
    # betahat <- splinefun(cl, rfpcrReg$beta1)
    # yhatTrain <- sapply(1:length(Xtrain$x), function(i) {
    #   mean(Xtrain$x[[i]] * betahat(Xtrain$pp[[i]]))
    # })
    # yhatTest <- sapply(1:length(Xtest$x), function(i) {
    #   mean(Xtest$x[[i]] * betahat(Xtest$pp[[i]]))
    # })
    
    yhatTrain1<-rfpcrReg$yhatTrain1
    yhatTest1<-rfpcrReg$yhatTest1
    # plot(yhatTest~ytest)
    # abline(0,1)
    # estError=ifelse(!is.null(betaCoef),integrate( function(x) (betaCoef(x)-betahat(x))^2,
    #                                           lower=0,upper=1)$value,NA)
    
    # res <-
    #   rbind(
    #     res,
    #     data.frame(
    #       estError = ifelse(!is.null(betaCoef), mean((
    #         rfpcrReg$beta1 - betaCoef(cl)
    #       ) ^ 2), NA),
    #       predErrorTrain = mean((ytrain - yhatTrain) ^
    #                               2),
    #       predErrorTest = mean((ytest - yhatTest) ^ 2),
    #       comp.s = rfpcrReg$comp.ret
    #     )
    #   )
    res<-list(betahat=list(rfpcrReg$beta,rfpcrReg$beta1),
              grid=cl,
              yhatTest=list(yhatTest,yhatTest1),
              comp=rfpcrReg$comp.ret,
              cvErr=rfpcrReg$cvErr)
    
    return(res)
  },
  error = function(cond) {
    print(cond)
    # return(data.frame(
    #   estError = c(NA, NA),
    #   predErrorTrain = c(NA, NA),
    #   predErrorTest = c(NA, NA),
    #   comp.s = c(NA, NA)
    # ))
    
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

makeIrregMatrix <- function(X){
  tt<-sort(unique(unlist(X$pp)))
  N<-length(X$x)
  
  Xmat<-matrix(c(0),nrow=N,ncol=length(tt))
  
  for (i in 1:N){
    for (j in 1:length(tt)){
      if (tt[j]%in%X$pp[[i]]){
        Xmat[i,j]<-X$x[[i]][which(X$pp[[i]]==tt[j])]
      }
      else{
        Xmat[i,j]<-NA
      }
    }
  }
  return(list(Xmat=Xmat,cl=tt))
}

doSimulation<-function(N,n,K,method,sparseObs=3:5,fixK=FALSE,betaCoef=NULL,eps=NULL,contFactorX=5,contFactorY=10,sigma=0.1,...){
  ncov <- 50 # size of grid to estimate cov_X
  k.cv <- 5 # number of cv folds for selecting optimal bandwidths h
  hs.mu <- c(0.1,0.2) # for step 1 kernel of estimating mean function
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
  
  #fpcr
  resErr<-data.frame(estError=numeric(0),
                     predErrorTest=numeric(0),
                     comp=numeric(0))
  resSim<-list()
  #van aelst
  resFPCR<-doFPCReg(XCont,yCont,Xtest,y0test,K = K,cl=cl,method=method,fixK=fixK)
  
  if (!is.null(resFPCR)){
    resErr<-rbind(resErr,
                  # data.frame(estError=mean( (resRFPCR$betahat[[1]]-betaCoef(resRFPCR$grid) )^2/
                  #                                           mean(betaCoef(resRFPCR$grid)^2)),
                  #            predError=mean( (resRFPCR$yhatTest[[1]]-y0test)^2/
                  #                              y0test^2),
                  #            comp=resRFPCR$comp),
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
  #browser()
  if (!is.null(resRFPCR)){
    resErr<-rbind(resErr,
                  # data.frame(estError=mean( (resRFPCR$betahat[[1]]-betaCoef(resRFPCR$grid) )^2/
                  #                                           mean(betaCoef(resRFPCR$grid)^2)),
                  #            predError=mean( (resRFPCR$yhatTest[[1]]-y0test)^2/
                  #                              y0test^2),
                  #            comp=resRFPCR$comp),
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
  
  resErr
  resSim<-c(resSim,list(rfpcr=resRFPCR))
  resSim
  #browser()
  #pace
  # resPACE<-doPaceReg(XCont,yCont,Xtest,y0test,ncov=ncov,K=K,fixK=fixK)
  # resErr<-rbind(resErr,
  #               data.frame(estError=mean( (resPACE$betahat-betaCoef(resPACE$grid) )^2/
  #                                           mean(betaCoef(resPACE$grid)^2)),
  #                          predError=mean( (resPACE$yhatTest-y0test)^2/
  #                                            y0test^2),
  #                          comp=resPACE$comp.s))
  # resErr
  # resSim<-c(resSim,list(pace=resPACE))
  # resSim
  
  #sfpca
  resSFPCA<-doSparseFPCAReg(XCont,yCont,Xtest,y0test,hs.mu=hs.mu,hs.cov=hs.cov,K=K,ncov=ncov,method="LS",fixK=fixK)
  #browser()
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


