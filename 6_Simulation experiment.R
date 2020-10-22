library(RVAideMemoire)
library(randomForest)
library(mclust)
library(reshape2)
library(yarrr)

set.seed(42)

#### This script was used to generate the datasets in the simulation experiment
#### to train and test the classifiers and prevalence estimation methods
#### and to compute and save all the metrics on each dataset

#### All metrics were saved under the names "2019-05-20_metrics.3var.HPLC_high" etc.in the file "Saved objects"
#### "2019-05-20_metrics.3var.HPLC_high" refers to "high variance treatments" and "HPLC datasets", etc.

#### To run this script, you need to run first:
## "Sugar model" -> The model
## "Estimating A. ervi parameters" and "Defining variance treatments" -> To get parameters
## "simul.functions" -> The functions used to:
  # train the classifiers on each train dataset
  # test them on each test dataset
  # apply the prevalence estimation methods
  # compute the metrics


### Besides, We'll need these formulas for Random Forest

# HPLC
formula_RF_HPLC  <- as.formula(paste("Tmt.f ~ ", 
                                     paste(c('Glucose','Fructose','Sucrose','Erlose',
                                             'Melezitose','Stacchyose','Maltose','GF_Ratio','H_Ratio'), 
                                           collapse= "+")))


# Anthrone
formula_RF_Anth <- as.formula(paste("Tmt.f ~ ", paste(c('Fructose','resSug','RS_Ratio'), collapse= "+")))



#### LOW BETWEEN-GROUP VARIANCE TREATMENTS ####


### 1 : take the corresponding variance treatments

## All stored in the 'Saved Objects' folder

ratio1 <- readRDS(file='Saved objects/ratios.7')
stock1 <- readRDS(file='Saved objects/stock.7')
minGF1 <- readRDS(file='Saved objects/minGF.7')

ratio2 <- readRDS(file='Saved objects/ratios.8')
stock2 <- readRDS(file='Saved objects/stock.8')
minGF2 <- readRDS(file='Saved objects/minGF.8')

ratio3 <- readRDS(file='Saved objects/ratios.9')
stock3 <- readRDS(file='Saved objects/stock.9')
minGF3 <- readRDS(file='Saved objects/minGF.9')



### 2 : for each variance, simulate 6 dataset sizes (30,60,90,120,180,270) -> 50 JDD/size

#Parameters to be altered: stock, ratio and minGF
param.set <-list( list(stock1,ratio1,minGF1),list(stock2,ratio2,minGF2),list(stock3,ratio3,minGF3) )

simul.list.train <- list()

print(' -- Creating training datasets -- ')

for (i in 1:length(param.set))
{

print(i)
  
t1 <- sugar.simul(categ,n.categ=c(9,9,9),Unf=1,n.Unfed=c(5,2,2),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],
                           Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
t2 <- sugar.simul(categ,n.categ=c(18,20,20),Unf=1,n.Unfed=c(6,6,6),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                       Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
t3 <- sugar.simul(categ,n.categ=c(30,30,30),Unf=1,n.Unfed=c(10,10,10),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                       Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
t4 <- sugar.simul(categ,n.categ=c(39,40,40),Unf=1,n.Unfed=c(13,13,13),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                       Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
t5 <- sugar.simul(categ,n.categ=c(51,50,50),Unf=1,n.Unfed=c(17,17,17),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                  Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
t6 <- sugar.simul(categ,n.categ=c(60,60,60),Unf=1,n.Unfed=c(20,20,20),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                       Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
t7 <- sugar.simul(categ,n.categ=c(69,70,70),Unf=1,n.Unfed=c(23,23,23),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                  Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
t8 <- sugar.simul(categ,n.categ=c(81,80,80),Unf=1,n.Unfed=c(27,27,27),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                  Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
t9 <- sugar.simul(categ,n.categ=c(90,90,90),Unf=1,n.Unfed=c(30,30,30),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                       Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)

simul.list.train[[i]] <- list(t1,t2,t3,t4,t5,t6,t7,t8,t9)

}


### 3 : simulate 12 test datasets (same size) with different ### 3 : simulate 12 test datasets (same size) with different class distributionistribution

simul.list.test <- list()

print(' -- Creating test datasets -- ')

for (i in 1:length(param.set))
{
  
  print(i)

# Balanced: 1/3 - 1/3 - 1/ 3
test1 <- sugar.simul(categ,n.categ=c(150,150,150),Unf=1,n.Unfed=c(50,50,50),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]

#50% Unfed - 25% ApH - 25% EFN
test2 <- sugar.simul(categ,n.categ=c(225,113,112),Unf=1,n.Unfed=c(125,50,50),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]

#25% Unfed - 25% ApH - 50% EFN
test3 <- sugar.simul(categ,n.categ=c(113,225,112),Unf=1,n.Unfed=c(28,57,28),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]

# 1/3 EFN - 2/3 ApH
test4 <-  sugar.simul(categ,n.categ=c(9,141,300),Unf=1,n.Unfed=c(5,2,2),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]

# 2/3 Unfed - 1/3 ApH
test5 <- sugar.simul(categ,n.categ=c(300,5,145),Unf=1,n.Unfed=c(200,2,98),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]

# 2/3 EFN - 1/3 ApH
test6 <- sugar.simul(categ,n.categ=c(9,300,141),Unf=1,n.Unfed=c(5,2,2),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]

# 1/3 Unfed - 2/3 ApH
test7 <- sugar.simul(categ,n.categ=c(150,5,295),Unf=1,n.Unfed=c(48,2,100),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]

# 80% Unfed - 10% EFN - 10% ApH
test8 <- sugar.simul(categ,n.categ=c(360,45,45),Unf=1,n.Unfed=c(288,36,36),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]

# 10% Unfed - 10% EFN - 80% ApH
test9 <- sugar.simul(categ,n.categ=c(45,45,360),Unf=1,n.Unfed=c(5,4,36),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]

# 100% Unfed
test10 <- sugar.simul(categ,n.categ=c(440,5,5),Unf=1,n.Unfed=c(436,2,2),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]

# 100% ApH
test11 <-  sugar.simul(categ,n.categ=c(9,5,436),Unf=1,n.Unfed=c(5,2,2),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]

# 100% EFN
test12 <-  sugar.simul(categ,n.categ=c(9,436,5),Unf=1,n.Unfed=c(5,2,2),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]

simul.list.test[[i]] <- list(test1,test2,test3,test4,test5,test6,
                        test7,test8,test9,test10,test11,test12)

}

### 4 : train models on each train dataset and predict each test dataset, HPLC

print(' -- Training models -- ')

mod.metrics.HPLC <- list()

for (i in 1:length(param.set))
{
  
  cat('Variance Treatment',i,'HPLC \n\n')
  
  #Okay so we have one variance
  #In that variance, we have several sizes
  #We want to test each model of each size on each test set
  
  train.sets <- simul.list.train[[i]]
  test.sets <- simul.list.test[[i]]
  
  #list with results for each test set
  comp <- list()
  
  #There must be a way to do this using mapply
  
  #The first element of comp will be the metrics on all test sets, and so on
  for (j in 1: length(train.sets))
  {
    name <- paste0('size',j)
    comp[[name]] <- lapply(test.sets,function(Y){model.comp.HPLC(train.sets[[j]],Y)})
  }
  
name <- paste0('Var',i)
mod.metrics.HPLC[[name]] <- comp

}

#so in mod.metrics.HPLC, we'll have a list
#1st element is for 1st variance treatment
#In this, 1st element is for 1st size
#In this, 1st element is for 1st test set

saveRDS(mod.metrics.HPLC,file=paste0("Saved objects/",Sys.Date(),'_metrics.3var.HPLC_low'))

### 5 : train models on each train dataset and predict each test dataset, Anthrone

print(' -- Training models -- ')

mod.metrics.Anth <- list()

for (i in 1:length(param.set))
{
  
  cat('Variance Treatment',i,'Anthrone \n\n')
  
  #Okay so we have one variance
  #In that variance, we have several sizes
  #We want to test each model of each size on each test set
  
  train.sets <- simul.list.train[[i]]
  test.sets <- simul.list.test[[i]]
  
  #list with results for each test set
  comp <- list()
  
  #There must be a way to do this using mapply
  
  #The first element of comp will be the metrics on all test sets, and so on
  for (j in 1: length(train.sets))
  {
    name <- paste0('size',j)
    comp[[name]] <- lapply(test.sets,function(Y){model.comp.Anth(train.sets[[j]],Y)})
  }
  
  name <- paste0('Var',i)
  mod.metrics.Anth[[name]] <- comp
  
}

saveRDS(mod.metrics.Anth,file=paste0("Saved objects/",Sys.Date(),'_metrics.3var.Anth_low'))

#### MID BETWEEN-GROUP VARIANCE TREATMENTS ####

### 1 : take the corresponding variance treatments

ratio1 <- readRDS(file='Saved objects/ratios.4')
stock1 <- readRDS(file='Saved objects/stock.4')
minGF1 <- readRDS(file='Saved objects/minGF.4')

ratio2 <- readRDS(file='Saved objects/ratios.5')
stock2 <- readRDS(file='Saved objects/stock.5')
minGF2 <- readRDS(file='Saved objects/minGF.5')

ratio3 <- readRDS(file='Saved objects/ratios.6')
stock3 <- readRDS(file='Saved objects/stock.6')
minGF3 <- readRDS(file='Saved objects/minGF.6')



### 2 : for each variance, simulate 6 dataset sizes (30,60,90,120,180,270) -> 50 JDD/size

#Parameters to be altered: stock, ratio and minGF
param.set <-list( list(stock1,ratio1,minGF1),list(stock2,ratio2,minGF2),list(stock3,ratio3,minGF3) )

simul.list.train <- list()

print(' -- Creating training datasets -- ')

for (i in 1:length(param.set))
{
  
  print(i)
  
  t1 <- sugar.simul(categ,n.categ=c(9,9,9),Unf=1,n.Unfed=c(5,2,2),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t2 <- sugar.simul(categ,n.categ=c(18,20,20),Unf=1,n.Unfed=c(6,6,6),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t3 <- sugar.simul(categ,n.categ=c(30,30,30),Unf=1,n.Unfed=c(10,10,10),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t4 <- sugar.simul(categ,n.categ=c(39,40,40),Unf=1,n.Unfed=c(13,13,13),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t5 <- sugar.simul(categ,n.categ=c(51,50,50),Unf=1,n.Unfed=c(17,17,17),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t6 <- sugar.simul(categ,n.categ=c(60,60,60),Unf=1,n.Unfed=c(20,20,20),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t7 <- sugar.simul(categ,n.categ=c(69,70,70),Unf=1,n.Unfed=c(23,23,23),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t8 <- sugar.simul(categ,n.categ=c(81,80,80),Unf=1,n.Unfed=c(27,27,27),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t9 <- sugar.simul(categ,n.categ=c(90,90,90),Unf=1,n.Unfed=c(30,30,30),lim=12,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  
  simul.list.train[[i]] <- list(t1,t2,t3,t4,t5,t6,t7,t8,t9)
  
}


### 3 : simulate 12 test datasets (same size) with different class distribution

simul.list.test <- list()

print(' -- Creating test datasets -- ')

for (i in 1:length(param.set))
{
  
  print(i)
  
  # Balanced: 1/3 - 1/3 - 1/ 3
  test1 <- sugar.simul(categ,n.categ=c(150,150,150),Unf=1,n.Unfed=c(50,50,50),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  #50% Unfed - 25% ApH - 25% EFN
  test2 <- sugar.simul(categ,n.categ=c(225,113,112),Unf=1,n.Unfed=c(125,50,50),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  #25% Unfed - 25% ApH - 50% EFN
  test3 <- sugar.simul(categ,n.categ=c(113,225,112),Unf=1,n.Unfed=c(28,57,28),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 1/3 EFN - 2/3 ApH
  test4 <-  sugar.simul(categ,n.categ=c(9,141,300),Unf=1,n.Unfed=c(5,2,2),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 2/3 Unfed - 1/3 ApH
  test5 <- sugar.simul(categ,n.categ=c(300,5,145),Unf=1,n.Unfed=c(200,2,98),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 2/3 EFN - 1/3 ApH
  test6 <- sugar.simul(categ,n.categ=c(9,300,141),Unf=1,n.Unfed=c(5,2,2),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 1/3 Unfed - 2/3 ApH
  test7 <- sugar.simul(categ,n.categ=c(150,5,295),Unf=1,n.Unfed=c(48,2,100),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 80% Unfed - 10% EFN - 10% ApH
  test8 <- sugar.simul(categ,n.categ=c(360,45,45),Unf=1,n.Unfed=c(288,36,36),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 10% Unfed - 10% EFN - 80% ApH
  test9 <- sugar.simul(categ,n.categ=c(45,45,360),Unf=1,n.Unfed=c(5,4,36),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 100% Unfed
  test10 <- sugar.simul(categ,n.categ=c(440,5,5),Unf=1,n.Unfed=c(436,2,2),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 100% ApH
  test11 <-  sugar.simul(categ,n.categ=c(9,5,436),Unf=1,n.Unfed=c(5,2,2),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 100% EFN
  test12 <-  sugar.simul(categ,n.categ=c(9,436,5),Unf=1,n.Unfed=c(5,2,2),lim=12,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  simul.list.test[[i]] <- list(test1,test2,test3,test4,test5,test6,
                               test7,test8,test9,test10,test11,test12)
  
}

### 4 : train models on each train dataset and predict each test dataset, HPLC

print(' -- Training models -- ')

mod.metrics.HPLC <- list()

for (i in 1:length(param.set))
{
  
  cat('Variance Treatment',i,'HPLC \n\n')
  
  #Okay so we have one variance
  #In that variance, we have several sizes
  #We want to test each model of each size on each test set
  
  train.sets <- simul.list.train[[i]]
  test.sets <- simul.list.test[[i]]
  
  #list with results for each test set
  comp <- list()
  
  #There must be a way to do this using mapply
  
  #The first element of comp will be the metrics on all test sets, and so on
  for (j in 1: length(train.sets))
  {
    name <- paste0('size',j)
    comp[[name]] <- lapply(test.sets,function(Y){model.comp.HPLC(train.sets[[j]],Y)})
  }
  
  name <- paste0('Var',i)
  mod.metrics.HPLC[[name]] <- comp
  
}

#so in mod.metrics.HPLC, we'll have a list
#1st element is for 1st variance treatment
#In this, 1st element is for 1st size
#In this, 1st element is for 1st test set

saveRDS(mod.metrics.HPLC,file=paste0("Saved objects/",Sys.Date(),'_metrics.3var.HPLC_mid'))

### 5 : train models on each train dataset and predict each test dataset, Anthrone

print(' -- Training models -- ')

mod.metrics.Anth <- list()

for (i in 1:length(param.set))
{
  
  cat('Variance Treatment',i,'Anthrone \n\n')
  
  #Okay so we have one variance
  #In that variance, we have several sizes
  #We want to test each model of each size on each test set
  
  train.sets <- simul.list.train[[i]]
  test.sets <- simul.list.test[[i]]
  
  #list with results for each test set
  comp <- list()
  
  #There must be a way to do this using mapply
  
  #The first element of comp will be the metrics on all test sets, and so on
  for (j in 1: length(train.sets))
  {
    name <- paste0('size',j)
    comp[[name]] <- lapply(test.sets,function(Y){model.comp.Anth(train.sets[[j]],Y)})
  }
  
  name <- paste0('Var',i)
  mod.metrics.Anth[[name]] <- comp
  
}

saveRDS(mod.metrics.Anth,file=paste0("Saved objects/",Sys.Date(),'_metrics.3var.Anth_mid'))

####  HIGH BETWEEN-GROUP VARIANCE TREATMENTS ####

### Formulas for Random Forest

formula_RF_HPLC  <- as.formula(paste("Tmt.f ~ ", 
                                     paste(c('Glucose','Fructose','Sucrose','Erlose',
                                             'Melezitose','Stacchyose','Maltose','GF_Ratio','H_Ratio'), 
                                           collapse= "+")))


formula_RF_Anth <- as.formula(paste("Tmt.f ~ ", paste(c('Fructose','resSug','RS_Ratio'), collapse= "+")))

### 1 : take the corresponding variance treatments

ratio1 <- readRDS(file='Saved objects/ratios.1')
stock1 <- readRDS(file='Saved objects/stock.1')
minGF1 <- readRDS(file='Saved objects/minGF.1')

ratio2 <- readRDS(file='Saved objects/ratios.2')
stock2 <- readRDS(file='Saved objects/stock.2')
minGF2 <- readRDS(file='Saved objects/minGF.2')

ratio3 <- readRDS(file='Saved objects/ratios.3')
stock3 <- readRDS(file='Saved objects/stock.3')
minGF3 <- readRDS(file='Saved objects/minGF.3')



### 2 : for each variance, simulate 6 dataset sizes (30,60,90,120,180,270) -> 50 JDD/size

#Parameters to be altered: stock, ratio and minGF
param.set <-list( list(stock1,ratio1,minGF1),list(stock2,ratio2,minGF2),list(stock3,ratio3,minGF3) )

simul.list.train <- list()

print(' -- Creating training datasets -- ')

for (i in 1:length(param.set))
{
  
  print(i)
  
  t1 <- sugar.simul(categ,n.categ=c(9,9,9),Unf=1,n.Unfed=c(5,2,2),lim=1,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t2 <- sugar.simul(categ,n.categ=c(18,20,20),Unf=1,n.Unfed=c(6,6,6),lim=1,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t3 <- sugar.simul(categ,n.categ=c(30,30,30),Unf=1,n.Unfed=c(10,10,10),lim=1,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t4 <- sugar.simul(categ,n.categ=c(39,40,40),Unf=1,n.Unfed=c(13,13,13),lim=1,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t5 <- sugar.simul(categ,n.categ=c(51,50,50),Unf=1,n.Unfed=c(17,17,17),lim=1,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t6 <- sugar.simul(categ,n.categ=c(60,60,60),Unf=1,n.Unfed=c(20,20,20),lim=1,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t7 <- sugar.simul(categ,n.categ=c(69,70,70),Unf=1,n.Unfed=c(23,23,23),lim=1,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t8 <- sugar.simul(categ,n.categ=c(81,80,80),Unf=1,n.Unfed=c(27,27,27),lim=1,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  t9 <- sugar.simul(categ,n.categ=c(90,90,90),Unf=1,n.Unfed=c(30,30,30),lim=1,param.set[[i]][[1]],param.set[[i]][[2]],
                    Times,nb.simul=50,minGF=param.set[[i]][[3]],maxTS=600)
  
  simul.list.train[[i]] <- list(t1,t2,t3,t4,t5,t6,t7,t8,t9)
  
}


### 3 : simuler 12 JDD de test (même taille) de distrib. différentes

simul.list.test <- list()

print(' -- Creating test datasets -- ')

for (i in 1:length(param.set))
{
  
  print(i)
  
  # Balanced: 1/3 - 1/3 - 1/ 3
  test1 <- sugar.simul(categ,n.categ=c(150,150,150),Unf=1,n.Unfed=c(50,50,50),lim=1,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  #50% Unfed - 25% ApH - 25% EFN
  test2 <- sugar.simul(categ,n.categ=c(225,113,112),Unf=1,n.Unfed=c(125,50,50),lim=1,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  #25% Unfed - 25% ApH - 50% EFN
  test3 <- sugar.simul(categ,n.categ=c(113,225,112),Unf=1,n.Unfed=c(28,57,28),lim=1,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 1/3 EFN - 2/3 ApH
  test4 <-  sugar.simul(categ,n.categ=c(9,141,300),Unf=1,n.Unfed=c(5,2,2),lim=1,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 2/3 Unfed - 1/3 ApH
  test5 <- sugar.simul(categ,n.categ=c(300,5,145),Unf=1,n.Unfed=c(200,2,98),lim=1,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 2/3 EFN - 1/3 ApH
  test6 <- sugar.simul(categ,n.categ=c(9,300,141),Unf=1,n.Unfed=c(5,2,2),lim=1,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 1/3 Unfed - 2/3 ApH
  test7 <- sugar.simul(categ,n.categ=c(150,5,295),Unf=1,n.Unfed=c(48,2,100),lim=1,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 80% Unfed - 10% EFN - 10% ApH
  test8 <- sugar.simul(categ,n.categ=c(360,45,45),Unf=1,n.Unfed=c(288,36,36),lim=1,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 10% Unfed - 10% EFN - 80% ApH
  test9 <- sugar.simul(categ,n.categ=c(45,45,360),Unf=1,n.Unfed=c(5,4,36),lim=1,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 100% Unfed
  test10 <- sugar.simul(categ,n.categ=c(440,5,5),Unf=1,n.Unfed=c(436,2,2),lim=1,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 100% ApH
  test11 <-  sugar.simul(categ,n.categ=c(9,5,436),Unf=1,n.Unfed=c(5,2,2),lim=1,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  # 100% EFN
  test12 <-  sugar.simul(categ,n.categ=c(9,436,5),Unf=1,n.Unfed=c(5,2,2),lim=1,stock=param.set[[i]][[1]],ratios=param.set[[i]][[2]],Times,nb.simul=1,minGF=param.set[[i]][[3]],maxTS=600)[[1]]
  
  simul.list.test[[i]] <- list(test1,test2,test3,test4,test5,test6,
                               test7,test8,test9,test10,test11,test12)
  
}

### 4 : train models on each train dataset and predict each test dataset, HPLC

print(' -- Training models -- ')

mod.metrics.HPLC <- list()

for (i in 1:length(param.set))
{
  
  cat('Variance Treatment',i,'HPLC \n\n')
  
  #Okay so we have one variance
  #In that variance, we have several sizes
  #We want to test each model of each size on each test set
  
  train.sets <- simul.list.train[[i]]
  test.sets <- simul.list.test[[i]]
  
  #list with results for each test set
  comp <- list()
  
  #There must be a way to do this using mapply
  
  #The first element of comp will be the metrics on all test sets, and so on
  for (j in 1: length(train.sets))
  {
    name <- paste0('size',j)
    comp[[name]] <- lapply(test.sets,function(Y){model.comp.HPLC(train.sets[[j]],Y)})
  }
  
  name <- paste0('Var',i)
  mod.metrics.HPLC[[name]] <- comp
  
}

#so in mod.metrics.HPLC, we'll have a list
#1st element is for 1st variance treatment
#In this, 1st element is for 1st size
#In this, 1st element is for 1st test set

saveRDS(mod.metrics.HPLC,file=paste0("Saved objects/",Sys.Date(),'_metrics.3var.HPLC_high'))

### 5 : train models on each train dataset and predict each test dataset, Anthrone

print(' -- Training models -- ')

mod.metrics.Anth <- list()

for (i in 1:length(param.set))
{
  
  cat('Variance Treatment',i,'Anthrone \n\n')
  
  #Okay so we have one variance
  #In that variance, we have several sizes
  #We want to test each model of each size on each test set
  
  train.sets <- simul.list.train[[i]]
  test.sets <- simul.list.test[[i]]
  
  #list with results for each test set
  comp <- list()
  
  #There must be a way to do this using mapply
  
  #The first element of comp will be the metrics on all test sets, and so on
  for (j in 1: length(train.sets))
  {
    name <- paste0('size',j)
    comp[[name]] <- lapply(test.sets,function(Y){model.comp.Anth(train.sets[[j]],Y)})
  }
  
  name <- paste0('Var',i)
  mod.metrics.Anth[[name]] <- comp
  
}

saveRDS(mod.metrics.Anth,file=paste0("Saved objects/",Sys.Date(),'_metrics.3var.Anth_high'))
