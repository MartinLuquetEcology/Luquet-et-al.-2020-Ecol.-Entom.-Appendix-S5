## Little test
## Let's investigate the difference between BC and KLD

## Indeed, it seems that adjustment performs differently if we consider BC or KLD as a performance measure
## Notably for some test datasets

## See the end of Appendix S9 for more details and complementary information

## Let's try this with a variance for which it is clear -> the 'most variable' sugar profile situation
## Var 3 --> Low between, High within

## To run this script, you first need to run
## "Sugar model" (to simulate datasets)
## "Estimating A. ervi parameters" and "Defining variance parameters"

## Parameters we'll use are stored in the "Saved objects" folder
ratio3 <- readRDS(file='Saved objects/ratios.9')
stock3 <- readRDS(file='Saved objects/stock.9')
minGF3 <- readRDS(file='Saved objects/minGF.9')

### Packages

library(RVAideMemoire)
library(randomForest)
library(mclust)
library(reshape2)
library(yarrr)

### Formulas for Random Forest

formula_RF_HPLC  <- as.formula(paste("Tmt.f ~ ", 
                                     paste(c('Glucose','Fructose','Sucrose','Erlose',
                                             'Melezitose','Stacchyose','Maltose','GF_Ratio','H_Ratio'), 
                                           collapse= "+")))


formula_RF_Anth <- as.formula(paste("Tmt.f ~ ", paste(c('Fructose','resSug','RS_Ratio'), collapse= "+")))

#### Let's simulate 100 train datasets with 270 size

set.seed(42)

train.sets <- sugar.simul(categ,n.categ=c(90,90,90),Unf=1,n.Unfed=c(30,30,30),lim=12,stock3,ratio3,
                  Times,nb.simul=100,minGF=minGF3,maxTS=600)

#### Okay, so now let's simulate two test datasets that seem 'problematic'

# 80% Unfed - 10% EFN - 10% ApH
test1 <- sugar.simul(categ,n.categ=c(360,45,45),Unf=1,n.Unfed=c(288,36,36),lim=12,stock3,ratio3,Times,nb.simul=1,minGF3,maxTS=600)[[1]]

# 10% Unfed - 10% EFN - 80% ApH
test2 <- sugar.simul(categ,n.categ=c(45,45,360),Unf=1,n.Unfed=c(5,4,36),lim=12,stock3,ratio3,Times,nb.simul=1,minGF3,maxTS=600)[[1]]

####### Now, let's fit a random forest on each train dataset
####### And then test it on both test datasets
####### We keep confusion matrices and prevalence estimations

##### For now let's do it on 1 dataset

test.set <- test1

confMat.list <- list()
prevEst.unadj.list <- list()
prevEst.adj.list <- list()

acc <- data.frame(set=1:length(train.sets),RF=NA)
kld.unadj <- data.frame(set=1:length(train.sets),RF=NA)
kld.adj <- data.frame(set=1:length(train.sets),RF=NA)
bc.unadj <- kld.unadj 
bc.adj <- kld.adj

levels(test.set$Tmt.f) <- levels(train.sets[[1]]$Tmt.f)
prev.Real <- summary(test.set$Tmt.f)


for (k in 1:length(train.sets))
{
  train.set <- train.sets[[k]]
  print(k)
  
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

  confMat.list[[k]] <- table.RF
  prevEst.unadj.list[[k]] <- unadj.prev.RF
  prevEst.adj.list[[k]] <- adj.prev.RF
  
  
}

##### Some plots

par(mfrow=c(1,1))

plot(kld.unadj$RF,kld.adj$RF,xlim=c(0,0.5),ylim=c(0,0.5),main='Kullback-Leibler Divergence (KLD)',
     ylab='KLD using Adjusted Counting',xlab='KLD using Classify and Count')
abline(0,1)

plot(bc.unadj$RF,bc.adj$RF,xlim=c(0,0.3),ylim=c(0,0.3),main='Bray-Curtis Dissimilarity (BC)',
     ylab='BC using Adjusted Counting',xlab='BC using Classify and Count')
abline(0,1)

##### ---> We can clearly see that ADj is better in terms of BC but not in terms of KLD

KLDRatio <- kld.adj$RF/kld.unadj$RF
BCRatio <- bc.adj$RF/bc.unadj$RF

plot(KLDRatio,BCRatio)

##### ---> The relationship is clearly not linear


## Let's look at cases when AC gives high KLD

plot(kld.unadj$RF[kld.adj$RF>0.3],kld.adj$RF[kld.adj$RF>0.3])
plot(bc.unadj$RF[kld.adj$RF>0.3],bc.adj$RF[kld.adj$RF>0.3])

## Let's plot some of theses cases to see when we have contrasting responses in terms of KLD and BC
## (See Appendix S9 for more details, at the end of the document)

table <- data.frame(bc.unadj=bc.unadj$RF,bc.adj=bc.adj$RF,kld.unadj=kld.unadj$RF,kld.adj=kld.adj$RF)

table[table$kld.adj>0.3,]

barplot(cbind(prev.Real,prevEst.unadj.list[[2]],prevEst.adj.list[[2]]),
        names.arg=c('Real','CC','AC'),legend=T)

text(1.9,300,'BC = 0.15')
text(1.9,250,'KLD = 0.06',col='red')
text(3.1,300,'BC = 0.10')
text(3.1,250,'KLD = 0.36',col='red')


barplot(cbind(prev.Real,prevEst.unadj.list[[6]],prevEst.adj.list[[6]]),
        names.arg=c('Real','CC','AC'))

text(1.9,300,'BC = 0.16')
text(1.9,250,'KLD = 0.10',col='red')
text(3.1,300,'BC = 0.07')
text(3.1,250,'KLD = 0.37',col='red')

barplot(cbind(prev.Real,prevEst.unadj.list[[9]],prevEst.adj.list[[9]]),
        names.arg=c('Real','CC','AC'))

text(1.9,300,'BC = 0.21')
text(1.9,250,'KLD = 0.12',col='red')
text(3.1,300,'BC = 0.13')
text(3.1,250,'KLD = 0.40',col='red')

barplot(cbind(prev.Real,prevEst.unadj.list[[74]],prevEst.adj.list[[74]]),
        names.arg=c('Real','CC','AC'))

text(1.9,300,'BC = 0.17')
text(1.9,250,'KLD = 0.09',col='red')
text(3.1,300,'BC = 0.14')
text(3.1,250,'KLD = 0.41',col='red')
