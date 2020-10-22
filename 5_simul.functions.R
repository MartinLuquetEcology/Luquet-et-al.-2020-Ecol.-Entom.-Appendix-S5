## Functions applied on simulated datasets


### Comparing classification and quantification performances

  ### "HPLC" datasets -> model.comp.HPLC

model.comp.HPLC <- function(train.sets,test.set)
  
{
  
  acc <- data.frame(set=1:length(train.sets),RF=NA,GMM=NA,DA=NA,Thresh=NA)
  kld.unadj <- data.frame(set=1:length(train.sets),RF=NA,GMM=NA,DA=NA,Thresh=NA)
  kld.adj <- data.frame(set=1:length(train.sets),RF=NA,GMM=NA,DA=NA,Thresh=NA)
  bc.unadj <- kld.unadj 
  bc.adj <- kld.adj
  
  #For the tables we need the same treatments in both datasets
  levels(test.set$Tmt.f) <- levels(train.sets[[1]]$Tmt.f)
  prev.Real <- summary(test.set$Tmt.f)
  
  for (k in 1:length(train.sets))
  {
    train.set <- na.omit(train.sets[[k]])
    
    #Sometimes we have only 0 values for erlose -> scaling is impossible
    if(sum(train.set$Erlose)==0){train.set$Erlose[1]<-0.00001}
    
    #DA
    DA <- try(MVA.cmv(scale(train.set[,c(4:11,13)]),train.set$Tmt.f,model="PPLS-DA/LDA",crit.inn="NMC",
                  repet=10,kout=6,kinn=5,ncomp=7),silent=T)
    
    if(class(DA)[1] != 'try-error'){
    
    #Accuracy
    predict.DA <- try(predict(DA,stand(test.set[,c(4:11,13)],scale(train.set[,c(4:11,13)])))$Group,
                      silent=T)
    
    if(class(predict.DA) != 'try-error'){
    
    table.DA <- table(test.set$Tmt.f,predict.DA)
    acc$DA[k] <- sum(diag(table.DA))/sum(table.DA)
    
    #KLD
    unadj.prev.DA <- colSums(table.DA)
    pred.train.DA <- levels(train.set$Tmt.f)[apply(DA$pred.prob,1,which.max)]
    table.train.DA <- table(train.set$Tmt.f,pred.train.DA )
    adj.prev.DA <- adjust.prev(table.train.DA,table.DA)
    
    kld.unadj$DA[k] <- kld(prev.Real,unadj.prev.DA)
    bc.unadj$DA[k] <- bc(prev.Real,unadj.prev.DA)
    
    if(!is.na(adj.prev.DA[1])){
      kld.adj$DA[k] <- kld(prev.Real,adj.prev.DA)
      bc.adj$DA[k] <- bc(prev.Real,adj.prev.DA)
    }else{
      kld.adj$DA[k] <- NA 
      bc.adj$DA[k] <-NA
    }
    
    }
    }
  
    #Random Forest
    rf <- randomForest(formula_RF_HPLC , data=train.set, ntree=1000) 
    
    #Accuracy
    table.RF <- table(test.set$Tmt.f,predict(rf,newdata=test.set))
    acc$RF[k] <- sum(diag(table.RF))/sum(table.RF)
    
    #KLD
    unadj.prev.RF <-  colSums(table.RF)
    table.train.RF <- rf$confusion[,-ncol(rf$confusion)]
    adj.prev.RF <- adjust.prev(table.train.RF,table.RF)
    
    kld.unadj$RF[k] <- kld(prev.Real,unadj.prev.RF)
    bc.unadj$RF[k] <- bc(prev.Real,unadj.prev.RF)
    
    if(!is.na(adj.prev.RF[1])){
      kld.adj$RF[k] <- kld(prev.Real,adj.prev.RF)
      bc.adj$RF[k] <- bc(prev.Real,adj.prev.RF)
    }else{
      kld.adj$RF[k] <- NA 
      bc.adj$RF[k] <-NA
    }
    
    #GMM
    GMM <- try(MclustDA(train.set[,c(4:11,13)],train.set$Tmt.f,verbose=F),silent=T)
    
    if(class(GMM) != 'try-error'){
    #Accuracy
    predict.GMM <- try(predict(GMM,newdata=test.set[,names(train.set[,c(4:11,13)])])$classification,silent=T)
    if(class(predict.GMM) != 'try-error'){
    
    table.GMM <- table(test.set$Tmt.f,predict.GMM)
    acc$GMM[k] <-  sum(diag(table.GMM))/sum(table.GMM)
    
    #KLD
    unadj.prev.GMM <- colSums(table.GMM)
    table.train.GMM <- table(train.set$Tmt.f,predict(GMM)$classification)
    adj.prev.GMM <- adjust.prev(table.train.GMM,table.GMM)
    
    kld.unadj$GMM[k] <- kld(prev.Real,unadj.prev.GMM)
    bc.unadj$GMM[k] <- bc(prev.Real,unadj.prev.GMM)
    
    if(!is.na(adj.prev.GMM[1])){
      kld.adj$GMM[k] <- kld(prev.Real,adj.prev.GMM)
      bc.adj$GMM[k] <- bc(prev.Real,adj.prev.GMM)
    }else{
      kld.adj$GMM[k] <- NA 
      bc.adj$GMM[k] <-NA
    }
    }
    }
    #Threshold
    
    #Accuracy
    GF.Unfed <- train.set$GF_Ratio[train.set$Tmt.f=="Unfed"]
    TS.Unfed <- train.set$TotSug[train.set$Tmt.f=="Unfed"]
    limGF <- round(0.05*length(GF.Unfed))
    if(limGF==0){limGF <- 1}
    GF.thresh <- sort(GF.Unfed)[limGF]
    TS.thresh <- sort(TS.Unfed)[0.95*length(TS.Unfed)]
    HR.notHoneydew <- train.set$H_Ratio[train.set$Tmt.f=="EFN"|train.set$Tmt.f=="Unfed"]
    HR.thresh <- sort(HR.notHoneydew)[0.95*length(HR.notHoneydew)]
    
    pred.Thresh <- factor(rep(NA,nrow(test.set)),levels=levels(train.set$Tmt.f))
    pred.Thresh[which(test.set$GF_Ratio>GF.thresh  | test.set$TotSug<TS.thresh)] <- "Unfed"
    pred.Thresh[is.na( pred.Thresh) & test.set$H_Ratio > HR.thresh ] <- "ApH"
    pred.Thresh[is.na( pred.Thresh)] <- "EFN"
    
    table.thresh <- table(test.set$Tmt.f,pred.Thresh)
    acc$Thresh[k] <- sum(diag(table.thresh))/sum(table.thresh)
    
    #KLD
    unadj.prev.Thresh <- colSums(table.thresh)
    
    pred.Thresh.train <- factor(rep(NA,nrow(train.set)),levels=levels(train.set$Tmt.f))
    pred.Thresh.train[which(train.set$GF_Ratio>GF.thresh  | train.set$TotSug<TS.thresh)] <- "Unfed"
    pred.Thresh.train[is.na( pred.Thresh.train) & train.set$H_Ratio > HR.thresh ] <- "ApH"
    pred.Thresh.train[is.na( pred.Thresh.train)] <- "EFN"
    
    table.thresh.train <- table(train.set$Tmt.f,pred.Thresh.train)
    adj.prev.Thresh <- adjust.prev(table.thresh.train,table.thresh)
    
    kld.unadj$Thresh[k] <- kld(prev.Real,unadj.prev.Thresh)
    bc.unadj$Thresh[k] <- bc(prev.Real,unadj.prev.Thresh)
    
    if(!is.na(adj.prev.Thresh[1])){
      kld.adj$Thresh[k] <- kld(prev.Real,adj.prev.Thresh)
      bc.adj$Thresh[k] <- bc(prev.Real,adj.prev.Thresh)
    }else{
      kld.adj$Thresh[k] <- NA 
      bc.adj$Thresh[k] <-NA
    }
    
  }
  
  metrics <- list(Accuracy=acc,Kld.unadj=kld.unadj,Kld.adj=kld.adj,BC.unadj=bc.unadj,BC.adj=bc.adj)
  return(metrics)
  
}


  ### "Anthrone" datasets -> model.comp.Anth

model.comp.Anth <- function(train.sets,test.set)
  
{
  
  acc <- data.frame(set=1:length(train.sets),RF=NA,GMM=NA,DA=NA,Thresh=NA)
  kld.unadj <- data.frame(set=1:length(train.sets),RF=NA,GMM=NA,DA=NA,Thresh=NA)
  kld.adj <- data.frame(set=1:length(train.sets),RF=NA,GMM=NA,DA=NA,Thresh=NA)
  bc.unadj <- kld.unadj 
  bc.adj <- kld.adj
  
  prev.Real <- summary(test.set$Tmt.f)
  
  for (k in 1:length(train.sets))
  {
    
    train.set <- train.sets[[k]]
    train.set <- na.omit(train.sets[[k]])
    
    #Sometimes we have only 0 values for erlose -> scaling is impossible
    if(sum(train.set$Erlose)==0){train.set$Erlose[1]<-0.00001}
    
    #DA
    DA <- try(MVA.cmv(scale(train.set[,c(4,14,15)]),train.set$Tmt.f,model="PPLS-DA/LDA",crit.inn="NMC",
                  repet=10,kout=6,kinn=5,ncomp=7),silent=T)
    
    if(class(DA)[1] != 'try-error'){
    
    #Accuracy
    predict.DA <- try(predict(DA,stand(test.set[,c(4,14,15)],scale(train.set[,c(4,14,15)])))$Group,silent=T)
    
    if(class(predict.DA) != 'try-error'){
      
    table.DA <- table(test.set$Tmt.f,predict.DA)
    acc$DA[k] <- sum(diag(table.DA))/sum(table.DA)
    
    #KLD
    unadj.prev.DA <- colSums(table.DA)
    pred.train.DA <- levels(test.set$Tmt.f)[apply(DA$pred.prob,1,which.max)]
    table.train.DA <- table(train.set$Tmt.f,pred.train.DA )
    adj.prev.DA <- adjust.prev(table.train.DA,table.DA)
    
    kld.unadj$DA[k] <- kld(prev.Real,unadj.prev.DA)
    bc.unadj$DA[k] <- bc(prev.Real,unadj.prev.DA)
    
    if(!is.na(adj.prev.DA[1])){
      kld.adj$DA[k] <- kld(prev.Real,adj.prev.DA)
      bc.adj$DA[k] <- bc(prev.Real,adj.prev.DA)
    }else{
      kld.adj$DA[k] <- NA 
      bc.adj$DA[k] <-NA
    }
    }
    }
    
    #Random Forest
    rf <- randomForest(formula_RF_Anth , data=train.set, ntree=1000) 
    
    #Accuracy
    table.RF <- table(test.set$Tmt.f,predict(rf,newdata=test.set))
    acc$RF[k] <- sum(diag(table.RF))/sum(table.RF)
    
    
    #KLD
    unadj.prev.RF <-  colSums(table.RF)
    table.train.RF <- rf$confusion[,-ncol(rf$confusion)]
    adj.prev.RF <- adjust.prev(table.train.RF,table.RF)
    
    kld.unadj$RF[k] <- kld(prev.Real,unadj.prev.RF)
    bc.unadj$RF[k] <- bc(prev.Real,unadj.prev.RF)
    
    if(!is.na(adj.prev.RF[1])){
      kld.adj$RF[k] <- kld(prev.Real,adj.prev.RF)
      bc.adj$RF[k] <- bc(prev.Real,adj.prev.RF)
    }else{
      kld.adj$RF[k] <- NA 
      bc.adj$RF[k] <-NA
    }
    
    #GMM
    GMM <- try(MclustDA(train.set[,c(4,14,15)],train.set$Tmt.f,verbose=F),silent=T)
    
    if(class(GMM) != 'try-error'){
    
    #Accuracy
    predict.GMM <- try(predict(GMM,newdata=test.set[,names(train.set[,c(4,14,15)])])$classification,silent=T)
    
    if(class(predict.GMM) != 'try-error'){
      
    table.GMM <- table(test.set$Tmt.f,predict.GMM)
    acc$GMM[k] <-  sum(diag(table.GMM))/sum(table.GMM)
    
    #KLD
    unadj.prev.GMM <- colSums(table.GMM)
    table.train.GMM <- table(train.set$Tmt.f,predict(GMM)$classification)
    adj.prev.GMM <- adjust.prev(table.train.GMM,table.GMM)
    
    kld.unadj$GMM[k] <- kld(prev.Real,unadj.prev.GMM)
    bc.unadj$GMM[k] <- bc(prev.Real,unadj.prev.GMM)
    
    if(!is.na(adj.prev.GMM[1])){
      kld.adj$GMM[k] <- kld(prev.Real,adj.prev.GMM)
      bc.adj$GMM[k] <- bc(prev.Real,adj.prev.GMM)
    }else{
      kld.adj$GMM[k] <- NA 
      bc.adj$GMM[k] <-NA
    }
    }
    }
    
    #Threshold
    RS.Unfed <- train.set$RS_Ratio[train.set$Tmt.f=="Unfed"]
    TS.Unfed <- train.set$TotSug[train.set$Tmt.f=="Unfed"]
    
    limRS <- round(0.05*length(RS.Unfed))
    if(limRS==0){limRS <- 1}
    RS.thresh <- sort(RS.Unfed)[limRS]
    TS.thresh <- sort(TS.Unfed)[0.95*length(TS.Unfed)]
    FR.notNect <- train.set$F_Ratio[train.set$Tmt.f=="ApH"|train.set$Tmt.f=="Unfed"]
    FR.thresh <- sort(FR.notNect)[0.95*length(FR.notNect)]
    
    #Accuracy
    pred.Thresh <- factor(rep(NA,nrow(test.set)),levels=levels(train.set$Tmt.f))
    pred.Thresh[which(test.set$RS_Ratio>RS.thresh  | test.set$TotSug<TS.thresh)] <- "Unfed"
    pred.Thresh[is.na( pred.Thresh) & test.set$F_Ratio <= FR.thresh ] <- "ApH"
    pred.Thresh[is.na( pred.Thresh)] <- "EFN"
    
    table.thresh <- table(test.set$Tmt.f,pred.Thresh)
    acc$Thresh[k] <- sum(diag(table.thresh))/sum(table.thresh)
    
    #KLD
    unadj.prev.Thresh <- colSums(table.thresh)
    
    pred.Thresh.train <- factor(rep(NA,nrow(train.set)),levels=levels(train.set$Tmt.f))
    pred.Thresh.train[which(train.set$RS_Ratio>RS.thresh  | train.set$TotSug<TS.thresh)] <- "Unfed"
    pred.Thresh.train[is.na( pred.Thresh.train) & train.set$F_Ratio <= FR.thresh ] <- "ApH"
    pred.Thresh.train[is.na( pred.Thresh.train)] <- "EFN"
    
    table.thresh.train <- table(train.set$Tmt.f,pred.Thresh.train)
    adj.prev.Thresh <- adjust.prev(table.thresh.train,table.thresh)
    
    kld.unadj$Thresh[k] <- kld(prev.Real,unadj.prev.Thresh)
    bc.unadj$Thresh[k] <- bc(prev.Real,unadj.prev.Thresh)
    
    if(!is.na(adj.prev.Thresh[1])){
      kld.adj$Thresh[k] <- kld(prev.Real,adj.prev.Thresh)
      bc.adj$Thresh[k] <- bc(prev.Real,adj.prev.Thresh)
    }else{
      kld.adj$Thresh[k] <- NA 
      bc.adj$Thresh[k] <-NA
    }
    
  }
  
  metrics <- list(Accuracy=acc,Kld.unadj=kld.unadj,Kld.adj=kld.adj,BC.unadj=bc.unadj,BC.adj=bc.adj)
  return(metrics)
  
}



#### Estimating prevalence via adjusted counting

adjust.prev <- function(train,test){
  
  options(error = NULL)  
  
  #If one class has a true positive rate of 0
  if (sum(diag(train) == 0) == 0) {
    
    #Error rates matrix
    rates <- matrix(0,nrow(train),ncol(train))
    
    for (i in 1:nrow(train))
    {
      falserates <- train[i,-i]/sum(train[i,])
      rates [-i,i] <- falserates
      rates [i,i] <- 1-sum(falserates)
      
    }
    
    #Unadjusted prevalences matrix
    nb.Ind <- sum(test)
    unadj <- as.matrix( colSums(test) / nb.Ind )
    
    if(class(try(solve(rates))) == "matrix"){
      
      #Solving
      est.vector <- solve(rates,unadj)
      
      if (sum(est.vector>0) != length(est.vector))
      {
        est.vector[est.vector<0] <- 0
        est.vector[est.vector>0] <- est.vector[est.vector>0] / sum(est.vector[est.vector>0])
      }
        
        
        
        #Adjusted prevalence values are stored here
        adj.prev <- as.numeric(round(est.vector*nb.Ind))
        names(adj.prev) <- rownames(train)
        
        if(all.equal(sum(adj.prev),nb.Ind) != TRUE) {  
          max.deci <- which.max(est.vector*nb.Ind-floor(est.vector*nb.Ind))
          adj.prev[max.deci] <- adj.prev[max.deci]+1
        } 
        
      
    }else{adj.prev <- NA}
  }else{adj.prev <- NA}
  
  return(adj.prev)
}

### Metrics calculation

## KLD

kld <- function(real,est) {
  
  tot <- sum(real)
  
  if(any(real==0)==T){real[real==0] <- 0.5}
  if(any(est==0)==T){est[est==0] <- 0.5}
  if(any(est==tot)==T){est[est==tot] <- tot-0.5}
  
  kld <- sum((real/tot)*log(real/est))
  return(kld)
}

## Bray-Curtis distance

bc <- function(real,est) {
  
  tab <- rbind(real,est)
  
  num <- 2 * sum(apply(tab,2,min))
  denom <- sum(tab)
  
  bc <- 1 - (num/denom)
  return(bc)
}