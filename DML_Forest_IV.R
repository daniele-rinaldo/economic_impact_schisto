DML_Forest_IV <- function(y,d,x,instr,tune,data,boot=FALSE,id,cl){
  
  cat('2-Fold cross-validated random forest 2SLS Double ML Estimation','\n')
  
  
  xsm    <- data[,x]
  ysm    <- data[,y]
  dsm    <- data[,d]
  ism    <- data[,instr]
  ni     <- length(instr)
  train  <- sample(1:nrow(data),floor(nrow(data)/2))
  # change for number of trees and cross-validation splits   
  ntree  <- 1000
  K      <- 2
  id_n   <- data[,id]
  

  if(tune==TRUE){
    
    cat('Tuned Random Forest','\n')
    cat('Calculating residuals on',y,'\n')
    
    forest_t1        <- tuneRF(xsm[-train,], ysm[-train], mtryStart=floor(sqrt(ncol(xsm[-train,]))), 
                               stepFactor=1.5, improve=0.05, nodesize=5, doBest=TRUE, plot=FALSE, trace=FALSE)
    min1             <- forest_t1$mtry
    forest_t2        <- tuneRF(xsm[train,], ysm[train], mtryStart=floor(sqrt(ncol(xsm[train,]))), 
                               stepFactor=1.5, improve=0.05, nodesize=5, doBest=TRUE, plot=FALSE, trace=FALSE)
    min2             <- forest_t2$mtry
    forest1.y        <- randomForest(xsm[-train,],ysm[-train],mtry=min1, ntree=ntree,maxnodes = 10,data=data);
    forest2.y        <- randomForest(xsm[train,],ysm[train], mtry=min2,ntree=ntree, maxnodes = 10,data=data)
  }
  
  if(tune==FALSE){
    
    cat('Untuned Random Forest','\n')
    cat('Calculating residuals on',y,'\n')
    
    forest1.y      <- randomForest(xsm[-train,],ysm[-train],ntree=ntree, maxnodes = 10,data=data);
    forest2.y      <- randomForest(xsm[train,],ysm[train], ntree=ntree, maxnodes = 10,data=data)
  }
  
  yhat1              <- predict(forest1.y,newdata=xsm[train,])
  yhat2              <- predict(forest2.y,newdata=xsm[-train,])
  
  # Compute projections on both samples = #
  
  cat('Calculating residuals on',d,'\n')
  if(tune==TRUE){
    forest_t1        <- tuneRF(xsm[-train,], dsm[-train], mtryStart=floor(sqrt(ncol(xsm[-train,]))), 
                               stepFactor=1.5, improve=0.05, nodesize=5, doBest=TRUE, plot=FALSE, trace=FALSE)
    tune_d1          <- forest_t1$mtry
    forest_t2        <- tuneRF(xsm[train,], dsm[train], mtryStart=floor(sqrt(ncol(xsm[train,]))), 
                                 stepFactor=1.5, improve=0.05, nodesize=5, doBest=TRUE, plot=FALSE, trace=FALSE)
    tune_d2          <- forest_t2$mtry
    forest1.d        <- randomForest(xsm[-train,],mtry=tune_d1,dsm[-train],maxnodes = 10)
    forest2.d        <- randomForest(xsm[train,],mtry=tune_d2,dsm[train],maxnodes = 10)
  }
  
  if(tune==FALSE){
    forest1.d        <- randomForest(xsm[-train,],dsm[-train],ntree=ntree,maxnodes = 10)
    forest2.d        <- randomForest(xsm[train,],dsm[train],ntree=ntree,maxnodes = 10)
  }
  
  
  dhat1            <- predict(forest1.d,newdata=xsm[train,])
  dhat2            <- predict(forest2.d,newdata=xsm[-train,])
  dv1              <- dsm[train]-dhat1
  dv2              <- dsm[-train]-dhat2
  
  
  if(tune==TRUE){
    forest_t1        <- tuneRF(xsm[-train,], dsm[-train], mtryStart=floor(sqrt(ncol(xsm[-train,]))), 
                               stepFactor=1.5, improve=0.05, nodesize=5, doBest=TRUE, plot=FALSE, trace=FALSE)
    min1               <- forest_t1$mtry
    forest_t2          <- tuneRF(xsm[train,], dsm[train], mtryStart=floor(sqrt(ncol(xsm[train,]))), 
                                 stepFactor=1.5, improve=0.05, nodesize=5, doBest=TRUE, plot=FALSE, trace=FALSE)
    min2             <- forest_t2$mtry
    forest1.d        <- randomForest(xsm[-train,],mtry=min1,dsm[-train],maxnodes = 10)
    forest2.d        <- randomForest(xsm[train,],mtry=min2,dsm[train],maxnodes = 10)
    
  
    for (i in 1:ni){
      cat('Calculating residuals on instruments: ',colnames(ism[i]),'\n')
      ismuse          <- ism[,i]
      forest_t_i1     <- tuneRF(xsm[-train,], ismuse[-train], mtryStart=floor(sqrt(ncol(xsm[-train,]))), 
                          stepFactor=1.5, improve=0.05, nodesize=5, doBest=TRUE, plot=FALSE, trace=FALSE)
      tune_i1          <- forest_t_i1$mtry
      forest_t_i2      <- tuneRF(xsm[train,], ismuse[train], mtryStart=floor(sqrt(ncol(xsm[train,]))), 
                                   stepFactor=1.5, improve=0.05, nodesize=5, doBest=TRUE, plot=FALSE, trace=FALSE)
      tune_i2          <- forest_t_i2$mtry
      forest1.i        <- randomForest(xsm[-train,],mtry=tune_i1,ismuse[-train],ntree=ntree,maxnodes = 10)
      forest2.i        <- randomForest(xsm[train,],mtry=tune_i2,ismuse[train],ntree=ntree,maxnodes = 10)
      ihat1            <- predict(forest1.i,newdata=xsm[train,])
      ihat2            <- predict(forest2.i,newdata=xsm[-train,])
      iv1              <- ismuse[train]-ihat1
      iv2              <- ismuse[-train]-ihat2
      assign(paste("iv1",i,sep="_"),iv1)
      assign(paste("iv2",i,sep="_"),iv2)
    }
  }
  
  if(tune==FALSE){
    
    for (i in 1:ni){
    cat('Calculating residuals on instruments: ',colnames(ism[i]),'\n')
    ismuse          <- ism[,i]
    forest1.i       <- randomForest(xsm[-train,],ismuse[-train],ntree=ntree,maxnodes = 10)
    forest2.i       <- randomForest(xsm[train,],ismuse[train],ntree=ntree,maxnodes = 10)
    ihat1            <- predict(forest1.i,newdata=xsm[train,])
    ihat2            <- predict(forest2.i,newdata=xsm[-train,])
    iv1              <- ismuse[train]-ihat1
    iv2              <- ismuse[-train]-ihat2
    assign(paste("iv1",i,sep="_"),iv1)
    assign(paste("iv2",i,sep="_"),iv2)
  }
  }

  
  
  cat('Computing the estimator','\n')
  # Compute 2-fold Cross-Validated DML 2sls estimator
  
  cat('part1','\n')
  datatrain          <- data[train,]
  datanotrain        <- data[-train,]

  datatrain$ysm      <- ysm[train]
  datanotrain$ysm    <- ysm[-train]

  datatrain$yhat1      <- yhat1
  datanotrain$yhat2    <- yhat2

  
  # # uncomment first two if auxiliary datasets above are not used
  # form_1            <- as.formula(paste("(ysm[train] - yhat1) ~ -1+dv1  | ", paste(paste("iv1_",seq(1,ni),sep=""),collapse="+"),sep=""))
  # form_2            <- as.formula(paste("(ysm[-train] - yhat2) ~ -1+dv2  | ", paste(paste("iv2_",seq(1,ni),sep=""),collapse="+"),sep=""))
  form_1            <- as.formula(paste("(ysm- yhat1) ~ -1+dv1  | ", paste(paste("iv1_",seq(1,ni),sep=""),collapse="+"),sep=""))
  form_2            <- as.formula(paste("(ysm -yhat2) ~ -1+dv2  | ", paste(paste("iv2_",seq(1,ni),sep=""),collapse="+"),sep=""))

  # DML_1_tsls         <- ivreg(form_1,data=datatrain)
  # DML_2_tsls         <- ivreg(form_2,data=datanotrain)
  # DML_1_tsls         <- ivreg(form_1,data=data[train])
  # DML_2_tsls         <- ivreg(form_2,data=data[-train])
  cat('part2','\n')
  form_1_first            <- as.formula(paste("dv1  ~ ", paste(paste("-1+iv1_",seq(1,ni),sep=""),collapse="+"),sep=""))
  form_2_first            <- as.formula(paste("dv2  ~ ", paste(paste("-1+iv2_",seq(1,ni),sep=""),collapse="+"),sep=""))

  datatrain$first1           <- fitted(lm(form_1_first))
  datanotrain$first2         <- fitted(lm(form_2_first))

  DML_1_tsls           <- lm(ysm-yhat1~-1+first1,data=datatrain)
  DML_2_tsls           <- lm(ysm-yhat2~-1+first2,data=datanotrain)
  #
   # DML_1_tsls           <- ivreg(form_1,data=data[train])
   # DML_2_tsls           <- lm(form_2,data=data[-train])


  theta1tsls         <- summary(DML_1_tsls)$coefficients[1,1]
  theta1tsls_se      <- coeftest(DML_1_tsls,vcov=vcovCL,cluster = as.formula(paste("~ ",cl)))[1,2]
  theta2tsls         <- summary(DML_2_tsls)$coefficients[1,1]
  theta2tsls_se      <- coeftest(DML_2_tsls,vcov=vcovCL,cluster = as.formula(paste("~ ",cl)))[1,2]
  # cat('part3','\n')
  theta_cv_DML_IV      <- mean(c(theta1tsls,theta2tsls))
  theta_cv_DML_IV_se   <- sqrt(1/K^2*(theta1tsls_se^2+theta2tsls_se^2))
  # theta_cv_DML_IV_se   <- sqrt(mean((summary(DML_1_tsls)$coef[1,2])^2,(summary(DML_2_tsls)$coef[1,2])^2))
  

  if(boot==TRUE){    
    cat('Bootstrapping standard errors','\n')
    boot1             <- coeftest(DML_1_tsls,vcov=vcovBS,cluster = as.formula(paste("~ ",id,cl,sep="+")), type="wild",R=3000,cores=4)
    boot2             <- coeftest(DML_2_tsls,vcov=vcovBS,cluster = as.formula(paste("~ ",id,cl,sep="+")), type="wild",R=3000,cores=4)
    
    theta_cv_DML_IV_se   <- sqrt(1/K^2*(boot1[1,2]^2+boot2[1,2]^2))
    
  }   
  
  theta_cv_DML_IV_pval <- 2*pt(-abs(theta_cv_DML_IV/theta_cv_DML_IV_se), df = nrow(data)-2 - length(unique(id_n)))
  
  
  cat('Calculating MSE','\n')
  # Mean Squared Errors in predicting y and d
  err_y1        <- sqrt(mean((yhat1-fitted(DML_1_tsls))^2))
  err_y2        <- sqrt(mean((yhat2-fitted(DML_2_tsls))^2))
  mse.y         <- mean(err_y1,err_y2)
  err_d1        <- sqrt(mean((dhat1-dsm[train])^2))
  err_d2        <- sqrt(mean((dhat2-dsm[-train])^2))
  mse.d         <- mean(err_d1,err_d2)
  
  
  return(list(theta_2sls=theta_cv_DML_IV,theta_2sls.se=theta_cv_DML_IV_se,theta_2sls.p=theta_cv_DML_IV_pval,
              mod1 = theta1tsls, mod2 = theta2tsls,
              mse.y=mse.y,mse.d=mse.d
  ))
  
}



