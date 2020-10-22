## Some model checks (more information in Appendix S3)

## This script must be run after "Sugar model" and "Estimating A. ervi parameters"

library(dplyr)
library(randomForest)
library(ggplot2)

### Preparing data table

  aervi.table <- read.table("aervi_data.txt")
  aervi.table$TimeF <- as.factor(aervi.table$Time)
 
  aervi.table$Tmt.f <- as.factor(aervi.table$Tmt2)
  aervi.table[aervi.table$Time>12,]$Tmt.f <- 'Unfed'

  ## Let's prepare a balanced dataset
    
  set.seed(42)

real.tab <- 
  rbind(
    sample_n(aervi.table[aervi.table$Sugar_treatment=="EFN" & aervi.table$Hour_treatment=="0h",]
             , size = 10)  ,
    sample_n(aervi.table[aervi.table$Sugar_treatment=="EFN" & aervi.table$Hour_treatment=="1h",]
             , size = 10)  ,
    sample_n(aervi.table[aervi.table$Sugar_treatment=="EFN" & aervi.table$Hour_treatment=="12h",]
             , size = 10) ,
    sample_n(aervi.table[aervi.table$Sugar_treatment=="EFN" & aervi.table$Hour_treatment=="24h",]
             , size = 5)  ,
    sample_n(aervi.table[aervi.table$Sugar_treatment=="EFN" & aervi.table$Hour_treatment=="48h",]
             , size = 5) ,
    sample_n(aervi.table[aervi.table$Sugar_treatment=="ApH" & aervi.table$Hour_treatment=="0h",]
             , size = 10)  ,
    sample_n(aervi.table[aervi.table$Sugar_treatment=="ApH" & aervi.table$Hour_treatment=="1h",]
             , size = 10)  ,
    sample_n(aervi.table[aervi.table$Sugar_treatment=="ApH" & aervi.table$Hour_treatment=="12h",]
             , size = 10) ,
    sample_n(aervi.table[aervi.table$Sugar_treatment=="ApH" & aervi.table$Hour_treatment=="24h",]
             , size = 5)  ,
    sample_n(aervi.table[aervi.table$Sugar_treatment=="ApH" & aervi.table$Hour_treatment=="48h",]
             , size = 5) ,
    sample_n(aervi.table[aervi.table$Sugar_treatment=="Unfed" & aervi.table$Hour_treatment=="0h",]
             , size = 2)  ,
    sample_n(aervi.table[aervi.table$Sugar_treatment=="Unfed" & aervi.table$Hour_treatment=="1h",]
             , size = 2)  ,
    sample_n(aervi.table[aervi.table$Sugar_treatment=="Unfed" & aervi.table$Hour_treatment=="12h",]
             , size = 2) ,
    sample_n(aervi.table[aervi.table$Sugar_treatment=="Unfed" & aervi.table$Hour_treatment=="24h",]
             , size = 2)  ,
    sample_n(aervi.table[aervi.table$Sugar_treatment=="Unfed" & aervi.table$Hour_treatment=="48h",]
             , size = 2)
      )


### Formula used for random forest

    formula_RF_HPLC  <- as.formula(paste("Tmt.f ~ ", 
                                        paste(c('Glucose','Fructose','Sucrose','Erlose',
                                                'Melezitose','Stacchyose','Maltose','GF_Ratio','H_Ratio'), 
                                              collapse= "+")))

### Random Forest on real data

  RF.real <- randomForest(formula_RF_HPLC,real.tab)
  RF.real$confusion

### Random Forest on simulated data

  ## Simulations
  
  categ <- c("Unfed","EFN","ApH")
  Unf <- 1
  lim <- 12
  Times <- c(0,1,12,24,48)

  simul2 <- sugar.simul(categ,n.categ=c(30,30,30),Unf,n.Unfed=c(10,10,10),lim,stock,ratios,Times,nb.simul=1,minGF=1.3,maxTS=300)[[1]]
  RF.simul2 <- randomForest(formula_RF_HPLC,simul2)
  RF.simul2

  ## Predicting real data with simulated data
  table(real.tab$Tmt.f,predict(RF.simul2,real.tab))
  RF.real$confusion

  ## Predicting simulated data with real data
  ## -- > Bias towards 'honeydew' samples
  table(simul2$Tmt.f,predict(RF.real,simul2))
  RF.simul2$confusion

  ## Relative Importance Plots
  ## Correspondance is quite good
  par(mfrow=c(1,1))
  varImpPlot(RF.real)
  varImpPlot(RF.simul2)

### Now let's do some checks, using multiple datasets

  ## Generating 100 datasets
simul.tables <- sugar.simul(categ,n.categ=c(30,30,30),Unf,n.Unfed=c(10,10,10),lim,stock,ratios,Times,nb.simul=100,minGF=1.3,maxTS=215)
  
  ## RF -> Confusion matrices for all of them
conf.all <- lapply(simul.tables, function(X) randomForest(formula_RF_HPLC,X)$confusion)

  ## RF -> Variable importance for all of them
imp.all <- lapply(simul.tables, function(X) randomForest(formula_RF_HPLC,X)$importance)

  ### Let's have a look at the error rate distribution

  ## Mean and SD values will be stored in those two tables
  simulTab <- RF.real$confusion[,1:3]
  simulTab.sd <- RF.real$confusion[,1:3]

for (i in 1:3){
  for(j in 1:3){
    
    simulTab[i,j] <- mean(unlist(lapply(conf.all,function(X) X[i,j])))
    simulTab.sd[i,j] <- sd(unlist(lapply(conf.all,function(X) X[i,j])))
    
  }
}

  # Adding mean error rates (calculated on means -> but also see below)
  simulTab <- cbind(simulTab, 1 - (diag(simulTab)/colSums(simulTab)))

simulTab
simulTab.sd

  ## If we compare those error rates to the 'real' error rates
  ## It's not perfectly nice in an absolute way, but it's quite okay in a relative way
RF.real$confusion

  ## Let's also get the mean and sd error rates per class

  err.rates <- lapply(simul.tables, function(X) 1 - (diag(randomForest(formula_RF_HPLC,X)$confusion)/30))
  class.err.rates <- do.call("rbind",err.rates)

  # Means
colMeans(class.err.rates)

  # SD
apply(class.err.rates,2,sd)


### Okay, now let's look at variable relative importance

  # Data will be stored here
impTabSug <- list()

for (i in 1:9)
{
  impTabSug[[i]] <- data.frame(Variable=rep(rownames(RF.real$importance)[i],100),Importance=0)
  impTabSug[[i]]$Importance <- unlist(lapply(imp.all,function(X) X[i,]))
  impTabSug[[i]]$Order <- unlist(lapply(imp.all,function(X) rank(-X)[i]))
}

impTab <- do.call("rbind",impTabSug)

  # Just ordering by decreasing order of importance
  impTab$Variable <- reorder(impTab$Variable, -impTab$Importance, mean)
  
  # Let's plot this:
p <- ggplot(impTab,aes(x=Variable,y=Importance))+
      geom_boxplot()+
      ylab('Importance (meanDecreaseGini)')

p
  # Let's add the real data on this:

pointsReal <- data.frame(Variable=rownames(RF.real$importance),ypos=RF.real$importance[,1])
  
p + geom_point(data=pointsReal,aes(x=Variable,y=ypos),shape=18,col='red',size=4)

## Quite nice ! (but see Fructose, Sucrose and H Ratio)

## To save the graph
#VRI_graph <- p + geom_point(data=pointsReal,aes(x=Variable,y=ypos),shape=18,col='red',size=4)
#ggsave('VRI_graph.png', VRI_graph,dpi=300,height=12,width=16,unit='cm')


### Now the graph about the order of importance:

VRI_order <- ggplot(impTab,aes(fill=Variable,x=as.factor(Order)))+
  geom_bar()+
  theme(axis.title.x=element_blank(),legend.position='bottom')

VRI_order

#ggsave('VRI_order.png', VRI_order,dpi=300,height=14,width=16,unit='cm')

## The real order of importance is:
rownames(RF.real$importance)[order(-RF.real$importance)]

### Finally, let's just draw some graphs to compare real and simulated datasets

  ## First let's simulate two more examples

simul3 <- sugar.simul(categ,n.categ=c(30,30,30),Unf,n.Unfed=c(10,10,10),lim,stock,ratios,Times,nb.simul=1,minGF=1.3,maxTS=300)[[1]]
simul4 <- sugar.simul(categ,n.categ=c(30,30,30),Unf,n.Unfed=c(10,10,10),lim,stock,ratios,Times,nb.simul=1,minGF=1.3,maxTS=300)[[1]]

  ## Then let's re-arrange a bit the data (naming the datasets, renaming the variables and treatments to homogeneize, merging the datasets)

real.tab$table <- 'Real'
simul2$table <- 'Simul 1'
simul3$table <- 'Simul 2'
simul4$table <- 'Simul 3'

levels(simul2$Tmt.f) <- c('Honeydew','Nectar','Unfed')
levels(simul3$Tmt.f) <- c('Honeydew','Nectar','Unfed')
levels(simul4$Tmt.f) <- c('Honeydew','Nectar','Unfed')

levels(simul2$Tmt) <- c('Honeydew','Nectar','Unfed')
levels(simul3$Tmt) <- c('Honeydew','Nectar','Unfed')
levels(simul4$Tmt) <- c('Honeydew','Nectar','Unfed')

real.tab$Tmt <- real.tab$Tmt2

allTables <- rbind(real.tab[,c('Glucose','Fructose','Sucrose','Melezitose','Erlose','Maltose','Stacchyose',
                               'GF_Ratio','H_Ratio','Tmt','Tmt.f','Time','table')],
                   simul2[,c('Glucose','Fructose','Sucrose','Melezitose','Erlose','Maltose','Stacchyose',
                             'GF_Ratio','H_Ratio','Tmt','Tmt.f','Time','table')],
                   simul3[,c('Glucose','Fructose','Sucrose','Melezitose','Erlose','Maltose','Stacchyose',
                             'GF_Ratio','H_Ratio','Tmt','Tmt.f','Time','table')],
                   simul4[,c('Glucose','Fructose','Sucrose','Melezitose','Erlose','Maltose','Stacchyose',
                             'GF_Ratio','H_Ratio','Tmt','Tmt.f','Time','table')])

  ## And now the graphs

## Red = Honeydew
## Blue = Unfed
## Green = Nectar

  ## Glucose
a <- ggplot(allTables,aes(x=as.factor(Time),y=Glucose,color=Tmt))+
  facet_wrap(~table)+
  geom_point(position=position_jitterdodge(jitter.width = 0.4,dodge.width=0.4))+
  ggtitle('Glucose')+
  labs(color='Feeding Treatment')+
  theme(legend.position = 'bottom')+
  ylab('Glucose (µg / mg insect)')+
  xlab('Time (h)')
  
  ## Fructose
b <- ggplot(allTables,aes(x=as.factor(Time),y=Fructose,color=Tmt))+
  facet_wrap(~table)+
  geom_point(position=position_jitterdodge(jitter.width = 0.4,dodge.width=0.4))+
  ggtitle('Fructose')+
  labs(color='Feeding Treatment')+
  theme(legend.position = 'bottom')+
  ylab('Fructose (µg / mg insect)')+
  xlab('Time (h)')

  ## Sucrose
c <- ggplot(allTables,aes(x=as.factor(Time),y=Sucrose,color=Tmt))+
  facet_wrap(~table)+
  geom_point(position=position_jitterdodge(jitter.width = 0.4,dodge.width=0.4))+
  ggtitle('Sucrose')+
  labs(color='Feeding Treatment')+
  theme(legend.position = 'bottom')+
  ylab('Sucrose (µg / mg insect)')+
  xlab('Time (h)')

  ## Melezitose
d <- ggplot(allTables,aes(x=as.factor(Time),y=Melezitose,color=Tmt))+
  facet_wrap(~table)+
  geom_point(position=position_jitterdodge(jitter.width = 0.4,dodge.width=0.4))+
  ggtitle('Melezitose')+
  labs(color='Feeding Treatment')+
  theme(legend.position = 'bottom')+
  ylab('Melezitose (µg / mg insect)')+
  xlab('Time (h)')

  ## Erlose
e <- ggplot(allTables,aes(x=as.factor(Time),y=Erlose,color=Tmt))+
  facet_wrap(~table)+
  geom_point(position=position_jitterdodge(jitter.width = 0.4,dodge.width=0.4))+
  ggtitle('Erlose')+
  labs(color='Feeding Treatment')+
  theme(legend.position = 'bottom')+
  ylab('Erlose (µg / mg insect)')+
  xlab('Time (h)')

  ## Maltose
f <- ggplot(allTables,aes(x=as.factor(Time),y=Maltose,color=Tmt))+
  facet_wrap(~table)+
  geom_point(position=position_jitterdodge(jitter.width = 0.4,dodge.width=0.4))+
  ggtitle('Maltose')+
  labs(color='Feeding Treatment')+
  theme(legend.position = 'bottom')+
  ylab('Maltose (µg / mg insect)')+
  xlab('Time (h)')

  ## Stacchyose
g <- ggplot(allTables,aes(x=as.factor(Time),y=Stacchyose,color=Tmt))+
  facet_wrap(~table)+
  geom_point(position=position_jitterdodge(jitter.width = 0.4,dodge.width=0.4))+
  ggtitle('Stacchyose')+
  labs(color='Feeding Treatment')+
  theme(legend.position = 'bottom')+
  ylab('Stacchyose (µg / mg insect)')+
  xlab('Time (h)')

  ## GF Ratio
h <- ggplot(allTables,aes(x=as.factor(Time),y=GF_Ratio,color=Tmt))+
  facet_wrap(~table)+
  geom_point(position=position_jitterdodge(jitter.width = 0.4,dodge.width=0.4))+
  ggtitle('GF Ratio  - Glucose / (Glucose + Fructose)')+
  labs(color='Feeding Treatment')+
  theme(legend.position = 'bottom')+
  xlab('Time (h)')

  ## Honeydew Ratio
i <- ggplot(allTables,aes(x=as.factor(Time),y=H_Ratio,color=Tmt))+
  facet_wrap(~table)+
  geom_point(position=position_jitterdodge(jitter.width = 0.4,dodge.width=0.4))+
  ggtitle('Honeydew Ratio - (Maltose + Erlose + Melezitose) / Total Sugar Amount')+
  labs(color='Feeding Treatment')+
  theme(legend.position = 'bottom', plot.title = element_text(size=12))+
  xlab('Time (h)')