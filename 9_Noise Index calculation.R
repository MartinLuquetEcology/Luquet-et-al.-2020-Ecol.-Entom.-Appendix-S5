##### CALCULATING THE 'NOISE INDEX' ####

## Here, we calculate the noise index for each variance treatment
## The goal is to see if we can define an index that correlates to the performance we obtain
## with each kind of dataset, notably in terms of prevalence estimation
## (see heatmap in the paper)

## The idea is to compute a ratio (inter-group var)/(inter-group var + intra-group var)
## Using RDA
## See Appendix S6

## Must be run after the "Sugar model" script

library(vegan)

## Here we store all the parameter values
## Run "Defining variance treatments" to get these values
stock <- list(stock1,stock2,stock3,stock4,stock5,stock6,stock7,stock8,stock9)
ratios <- list(ratios1,ratios2,ratios3,ratios4,ratios5,ratios6,ratios7,ratios8,ratios9)
minGF <- c(minGF1,minGF2,minGF3,minGF4,minGF5,minGF6,minGF7,minGF8,minGF9)

# Gradient of 'difficulty' according to the heatmap
diff <- c(1,3,6,2,5,7,4,8,9)

## So the order is
## Hb-LW (1), Mb-Lw(2), Hb-Mw(3), Lb-Lw(4), Mb-Mw(5), Hb-Hw(6), Lb-Mw(7), Mb-Hw(8), Lb-Hw(9)

profile <- c('Hb - Lw', 'Mb - Lw', 'Hb - Mw',
             'Lb - Lw','Mb - Mw','Hb - Hw',
             'Mb - Hw','Lb - Mw','Lb - Hw')

## Those are the parameters to determine dataset size and composition
## (We'll work only with balanced datasets here)
n.categL <- list(c(9,9,9),c(18,20,20),c(30,30,30),c(39,40,40),c(51,50,50),
                 c(60,60,60),c(69,70,70),c(81,80,80),c(90,90,90))
n.UnfedL <- list(c(5,2,2),c(6,6,6),c(10,10,10),c(13,13,13),c(17,17,17),
                 c(20,20,20),c(23,23,23),c(27,27,27),c(30,30,30))


#### Noise index calculation for HPLC ####


inert.unconstr.HPLC <- list()

# j = dataset sizes

for (j in 1:9)
{
  
  
  prop.expl <- list()
  
  # i = Variances
  
  for (i in 1:9)
  {
    
    tables <- sugar.simul(categ,n.categL[[j]],n.UnfedL[[j]],Unf=1,lim=1,stock[[i]],ratios[[i]],Times,nb.simul=100,minGF[i],maxTS)
    
    rda.an <- lapply(tables,function(X) rda(X[,c(4:11,13)],X$Tmt.f))
    prop.expl[[i]] <- sapply(rda.an,function(X) summary(X)$unconst.chi/summary(X)$tot.chi)
    
  }
  
  inert.unconstr.HPLC[[j]] <- data.frame(inertP = unlist(prop.expl),
                                       diff = rep(diff,each=100),
                                       size = sum(n.categL[[j]]),
                                       prof = rep(profile,each=100))
  
}

inert.table.HPLC <- do.call("rbind",inert.unconstr.HPLC)

ggplot(inert.table.HPLC,aes(x=as.factor(diff),y=inertP))+
  facet_wrap(~size)+
  geom_boxplot()+
  ggtitle('HPLC')



#### Noise index calculation for Anthrone ####

# Same thing on anthrone datasets

inert.unconstr.Anth <- list()

for (j in 1:9)
{
  
  
  prop.expl <- list()
  
  # i = Variances
  
  for (i in 1:9)
  {
    
    tables <- sugar.simul(categ,n.categL[[j]],n.UnfedL[[j]],Unf=1,lim=1,stock[[i]],ratios[[i]],Times,nb.simul=100,minGF[i],maxTS)
    
    rda.an <- lapply(tables,function(X) rda(X[,c(4,14,15)],X$Tmt.f))
    prop.expl[[i]] <- sapply(rda.an,function(X) summary(X)$unconst.chi/summary(X)$tot.chi)
    
  }
  
  inert.unconstr.Anth[[j]] <- data.frame(inertP = unlist(prop.expl),
                                       diff = rep(diff,each=100),
                                       size = sum(n.categL[[j]]))
  
}

inert.table.Anth <- do.call("rbind",inert.unconstr.Anth)

ggplot(inert.table.Anth,aes(x=as.factor(diff),y=index))+
  facet_wrap(~size)+
  geom_boxplot()+
  ggtitle('Anthrone')+
  scale_y_continuous(breaks = seq(0,1,0.1))
