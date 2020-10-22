library(randomForest)  
library(ggplot2)
library(vegan)

#### Here are stored and saved all the parameter values that were used for the
#### Simulation experiment
#### There are 9 treatments, corresponding to 
#### 3 (intra-group variance) * 3 (between-group variance) treatments

#### This script must be run after 'Sugar model' and Estimating A. ervi parameters'

#### Objects containing parameters are saved in the 'Saved Objects' folder

#### All parameter values are detailed in Appendix S4

####

#### PARAMETERS ####

categ <- c("Unfed","EFN","ApH")
n.categ <- c(30,30,30)
Unf <- 1
n.Unfed <- c(10,10,10)
lim <- 12
Times <- c(0,1,12,24,48)
nb.simul <- 1
maxTS <- 600

par(mfrow=c(2,1))

formula_RF_HPLC  <- as.formula(paste("Tmt.f ~ ", 
                                     paste(c('Glucose','Fructose','Sucrose','Erlose',
                                             'Melezitose','Stacchyose','Maltose','GF_Ratio','H_Ratio'), 
                                           collapse= "+")))


##### ---- HIGH BETWEEN-GROUP VARIANCE ---- ####

# High between, Low within

## "stock" and "ratios" are estimated from A. ervi dataset ('Estimating A. ervi parameters' script)
## We just round the values a little bit

stock1 <- list(Unfed=c(alpha=3,beta=-0.037,rse=0.1),
               EFN=c(alpha=5,beta=-0.12,rse=0.1),
               ApH=c(alpha=4.2, beta=-0.08, rse=0.1))


ratios1 <- ratios
ratios1$Unfed$`GF_ratio(t)`[,1] <- c(14,13,12,12,12)
ratios1$EFN$`GF_ratio(t)`[,1] <- c(1.2,1.5,6,6,6)
ratios1$ApH$`GF_ratio(t)`[,1] <- c(2.7,2.9,6,6,8)


ratios1$Unfed$`GF_ratio(t)`[,2] <- 1.5
ratios1$EFN$`GF_ratio(t)`[,2] <- 0.1
ratios1$ApH$`GF_ratio(t)`[,2] <- 0.25

minGF1 <- 1.5

table1 <- sugar.simul(categ,n.categ,n.Unfed,Unf=1,lim=1,stock1,ratios1,Times,nb.simul,minGF1,maxTS)[[1]]

plot(table1$GF_Ratio~table1$TotSug,col=table1$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='High between, low within',cex=1.5,pch=16,xlim=c(0,300))

RF1 <- randomForest(formula_RF_HPLC,table1)
RF1$confusion
varImpPlot(RF1)

##saveRDS(stock1,file='Saved objects/stock.1')
##saveRDS(ratios1,file='Saved objects/ratios.1')
##saveRDS(minGF1,file='Saved objects/minGF.1')

# High between, Mid within

  ## (Note: Episyrphus balteatus - like, cf Hogervorst et al. 2007)

stock2 <- list(Unfed=c(alpha=3,beta=-0.037,rse=0.3),
               EFN=c(alpha=5,beta=-0.12,rse=0.3),
               ApH=c(alpha=4.2, beta=-0.08, rse=0.3))

ratios2 <- ratios1
ratios2$Unfed$`GF_ratio(t)`[,1] <- c(14,12,8,8,6)
ratios2$Unfed$`GF_ratio(t)`[,2] <- 5
ratios2$EFN$`GF_ratio(t)`[,2] <- 0.25
ratios2$ApH$`GF_ratio(t)`[,2] <- 1

minGF2 <- 1.25

table2 <- sugar.simul(categ,n.categ,n.Unfed,Unf=1,lim=1,stock2,ratios2,Times,nb.simul,minGF2,maxTS)[[1]]

plot(table2$GF_Ratio~table2$TotSug,col=table2$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='High between, mid within',cex=1.5,pch=16,xlim=c(0,300))

RF2 <- randomForest(formula_RF_HPLC,table2)
RF2$confusion
varImpPlot(RF2)

#saveRDS(stock2,file='Saved objects/stock.2')
#saveRDS(ratios2,file='Saved objects/ratios.2')
#saveRDS(minGF2,file='Saved objects/minGF.2')

# High between, High within

  ## (Note: a bit Venturia canescens like, cf Desouhant et al. 2010)
  ## (But with less y-variance for the 'Unfed' group)

stock3 <- list(Unfed=c(alpha=3,beta=-0.037,rse=0.6),
                EFN=c(alpha=5,beta=-0.12,rse=0.6),
                ApH=c(alpha=4.2, beta=-0.08, rse=0.6))

ratios3 <- ratios2
ratios3$Unfed$`GF_ratio(t)`[,1] <- c(14,11,8,6,5)
ratios3$Unfed$`GF_ratio(t)`[,2] <- 10
ratios3$EFN$`GF_ratio(t)`[,2] <- c(0.5,0.5,1,1,1)
ratios3$ApH$`GF_ratio(t)`[,2] <- 1.5

minGF3 <- 1

table3 <- sugar.simul(categ,n.categ,n.Unfed,Unf=1,lim=1,stock3,ratios3,Times,nb.simul,minGF3,maxTS)[[1]]

plot(table3$GF_Ratio~table3$TotSug,col=table3$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='Low between, low within',cex=1.5,pch=16,xlim=c(0,300))

RF3 <- randomForest(formula_RF_HPLC,table3)
RF3$confusion
varImpPlot(RF3)

#saveRDS(stock3,file='Saved objects/stock.3')
#saveRDS(ratios3,file='Saved objects/ratios.3')
#saveRDS(minGF3,file='Saved objects/minGF.3')

##### ---- MID BETWEEN-GROUP VARIANCE ---- ####

# Mid between, Low within

stock4 <- list(Unfed=c(alpha=2.8,beta=-0.037,rse=0.15),
               EFN=c(alpha=4.8,beta=-0.044,rse=0.15),
               ApH=c(alpha=3.7, beta=-0.04, rse=0.15))


ratios4 <- ratios
ratios4$Unfed$`GF_ratio(t)`[,1] <-  c(12,11,8,7,6)
ratios4$EFN$`GF_ratio(t)`[,1] <- c(1.7,2.1,2.5,4.5,5)
ratios4$ApH$`GF_ratio(t)`[,1] <- c(3,3,3,5,6)


ratios4$Unfed$`GF_ratio(t)`[,2] <- 1.5
ratios4$EFN$`GF_ratio(t)`[,2] <- 0.05
ratios4$ApH$`GF_ratio(t)`[,2] <- 0.5

minGF4 <- 1.5

table4 <- sugar.simul(categ,n.categ,n.Unfed,Unf=1,lim=12,stock4,ratios4,Times,nb.simul,minGF4,maxTS)[[1]]

plot(table4$GF_Ratio~table4$TotSug,col=table4$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='Mid between, low within',cex=1.5,pch=16,xlim=c(0,300))

RF4 <- randomForest(formula_RF_HPLC,table4)
RF4$confusion
varImpPlot(RF4)

#saveRDS(stock4,file='Saved objects/stock.4')
#saveRDS(ratios4,file='Saved objects/ratios.4')
#saveRDS(minGF4,file='Saved objects/minGF.4')

# Mid between, Mid within

  ## (A bit Chrysoperla carnea like, cf Hogervorst et al. 2007)

stock5 <- list(Unfed=c(alpha=2.8,beta=-0.037,rse=0.4),
               EFN=c(alpha=4.8,beta=-0.044,rse=0.4),
               ApH=c(alpha=3.7, beta=-0.04, rse=0.4))

ratios5 <- ratios4
ratios5$Unfed$`GF_ratio(t)`[,2] <- 5
ratios5$EFN$`GF_ratio(t)`[,2] <- 0.25
ratios5$ApH$`GF_ratio(t)`[,2] <- 1

minGF5 <- 1.25

table5 <- sugar.simul(categ,n.categ,n.Unfed,Unf=1,lim=12,stock5,ratios5,Times,nb.simul,minGF5,maxTS)[[1]]

plot(table5$GF_Ratio~table5$TotSug,col=table5$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='Mid between, mid within',cex=1.5,pch=16,xlim=c(0,300))

RF5 <- randomForest(formula_RF_HPLC,table5)
RF5$confusion
varImpPlot(RF5)

#saveRDS(stock5,file='Saved objects/stock.5')
#saveRDS(ratios5,file='Saved objects/ratios.5')
#saveRDS(minGF5,file='Saved objects/minGF.5')

# Mid between, High within

  ## (A bit Aphelinus mali like, cf Ainara's dataset)

stock6 <- list(Unfed=c(alpha=2.8,beta=-0.037,rse=0.7),
               EFN=c(alpha=4.8,beta=-0.044,rse=1),
               ApH=c(alpha=3.7, beta=-0.04, rse=0.8))

ratios6 <- ratios5

ratios6$EFN$`GF_ratio(t)`[,1] <- c(1.7,2.1,2.5,3.5,4.5)
ratios6$ApH$`GF_ratio(t)`[,1] <- c(3,3,3,4,5)

ratios6$Unfed$`GF_ratio(t)`[,2] <- 10
ratios6$EFN$`GF_ratio(t)`[,2] <- c(0.5,0.5,0.5,1,1)
ratios6$ApH$`GF_ratio(t)`[,2] <- 1.5

minGF6 <- 1

table6 <- sugar.simul(categ,n.categ,n.Unfed,Unf=1,lim=12,stock6,ratios6,Times,nb.simul,minGF6,maxTS)[[1]]

plot(table6$GF_Ratio~table6$TotSug,col=table6$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='Mid between, high within',cex=1.5,pch=16,xlim=c(0,300))

RF6 <- randomForest(formula_RF_HPLC,table6)
RF6$confusion
varImpPlot(RF6)

#saveRDS(stock6,file='Saved objects/stock.6')
#saveRDS(ratios6,file='Saved objects/ratios.6')
#saveRDS(minGF6,file='Saved objects/minGF.6')

##### ---- LOW BETWEEN-GROUP VARIANCE ---- ####

# Low between, low within

  #Parameters directly inspired from Aphidius ervi's dataset

stock7 <- list(Unfed=c(alpha=3,beta=-0.037,rse=0.3),
                EFN=c(alpha=4.5,beta=-0.044,rse=0.3),
                ApH=c(alpha=3.7, beta=-0.041, rse=0.3))


ratios7 <- ratios
ratios7$Unfed$`GF_ratio(t)`[,1] <-  c(3.9,12.2,8.5,12.3,6.9)
ratios7$EFN$`GF_ratio(t)`[,1] <- c(1.7,1.9,2.2,2.5,2.8)
ratios7$ApH$`GF_ratio(t)`[,1] <- c(2,3.1,3.6,5.3,5.9)

ratios7$Unfed$`GF_ratio(t)`[,2] <- 1.5
ratios7$EFN$`GF_ratio(t)`[,2] <- 0.05
ratios7$ApH$`GF_ratio(t)`[,2] <- 0.5


minGF7 <- 1.5

table7 <- sugar.simul(categ,n.categ,n.Unfed,Unf=1,lim=14,stock7,ratios7,Times,nb.simul,minGF7,maxTS)[[1]]

plot(table7$GF_Ratio~table7$TotSug,col=table7$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='Low between, low within',cex=1.5,pch=16,xlim=c(0,300))

RF7 <- randomForest(formula_RF_HPLC,table7)
RF7$confusion
varImpPlot(RF7)

#saveRDS(stock7,file='Saved objects/stock.7')
#saveRDS(ratios7,file='Saved objects/ratios.7')
#saveRDS(minGF7,file='Saved objects/minGF.7')

# Low between, mid within

  # A.ervi - like (cf Martin's dataset & Hogervorst et al. 2007)

stock8 <- list(Unfed=c(alpha=3,beta=-0.037,rse=0.7),
                EFN=c(alpha=4.5,beta=-0.044,rse=1),
                ApH=c(alpha=3.7, beta=-0.041, rse=0.7))

ratios8 <- ratios7

ratios8$Unfed$`GF_ratio(t)`[,2] <- 5
ratios8$EFN$`GF_ratio(t)`[,2] <- 0.25
ratios8$ApH$`GF_ratio(t)`[,2] <- 1

minGF8 <- 1.25

table8 <- sugar.simul(categ,n.categ,n.Unfed,Unf=1,lim=13,stock8,ratios8,Times,nb.simul,minGF8,maxTS)[[1]]

plot(table8$GF_Ratio~table8$TotSug,col=table8$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='Low between, mid within',cex=1.5,pch=16,xlim=c(0,300))

RF8 <- randomForest(formula_RF_HPLC,table8)
RF8$confusion
varImpPlot(RF8)

#saveRDS(stock8,file='Saved objects/stock.8')
#saveRDS(ratios8,file='Saved objects/ratios.8')
#saveRDS(minGF8,file='Saved objects/minGF.8')


# Low between, high within

stock9 <- list(Unfed=c(alpha=3,beta=-0.037,rse=1),
                EFN=c(alpha=4.5,beta=-0.044,rse=1.5),
                ApH=c(alpha=3.7, beta=-0.041, rse=1))

ratios9 <- ratios7

ratios9$Unfed$`GF_ratio(t)`[,2] <- 10
ratios9$EFN$`GF_ratio(t)`[,2] <- 0.5
ratios9$ApH$`GF_ratio(t)`[,2] <- 1.5

minGF9 <- 1

table9 <- sugar.simul(categ,n.categ,n.Unfed,Unf=1,lim=12,stock9,ratios9,Times,nb.simul,minGF9,maxTS)[[1]]

plot(table9$GF_Ratio~table9$TotSug,col=table9$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='Low between, high within',cex=1.5,pch=16,xlim=c(0,300))

RF9 <- randomForest(formula_RF_HPLC,table9)
RF9$confusion
varImpPlot(RF9)

#saveRDS(stock9,file='Saved objects/stock.9')
#saveRDS(ratios9,file='Saved objects/ratios.9')
#saveRDS(minGF9,file='Saved objects/minGF.9')

#### ALL GRAPHS ####

par(mfrow=c(3,3))

plot(table1$GF_Ratio~table1$TotSug,col=table1$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='High between, low within',cex=1.5,pch=16,xlim=c(0,500))

plot(table2$GF_Ratio~table2$TotSug,col=table2$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='High between, mid within',cex=1.5,pch=16,xlim=c(0,500))

plot(table3$GF_Ratio~table3$TotSug,col=table3$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='High between, high within',cex=1.5,pch=16,xlim=c(0,500))

plot(table4$GF_Ratio~table4$TotSug,col=table4$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='Mid between, low within',cex=1.5,pch=16,xlim=c(0,500))

plot(table5$GF_Ratio~table5$TotSug,col=table5$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='Mid between, mid within',cex=1.5,pch=16,xlim=c(0,500))

plot(table6$GF_Ratio~table6$TotSug,col=table6$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='Mid between, high within',cex=1.5,pch=16,xlim=c(0,500))

plot(table7$GF_Ratio~table7$TotSug,col=table7$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='Low between, low within',cex=1.5,pch=16,xlim=c(0,500))

plot(table8$GF_Ratio~table8$TotSug,col=table8$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='Low between, mid within',cex=1.5,pch=16,xlim=c(0,500))

plot(table9$GF_Ratio~table9$TotSug,col=table9$Tmt.f,ylim=c(0.5,1),
     xlab='Total Sugars',ylab='GF Ratio',main='Low between, high within',cex=1.5,pch=16,xlim=c(0,500))
