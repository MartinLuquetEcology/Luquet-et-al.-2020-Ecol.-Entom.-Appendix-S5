#### Plotting the heatmaps
#### To run this script, you first need to run "Plotting functions for simulated data" and "Simulation results & graphs"

### First this function is created to generate the heatmaps

final.heatmap <- function(Anthrone,HPLC,metrics.unadj,metrics.adj,metrics.name,
                          metrics.name2)
  
  
{

listChem <- list(HPLC=HPLC,Anthrone=Anthrone)

  #### Anthrone

allvalues.unadj <- metrics.table(metrics.unadj,listChem[['Anthrone']])
allvalues.adj <- metrics.table(metrics.adj,listChem[['Anthrone']])

  # Unadjusted values
      ## Let's take only the wanted metrics
Anth.unadj <- allvalues.unadj[
  
  allvalues.unadj$size <90 |
    (allvalues.unadj$size<160 & allvalues.unadj$name %in% c('Mb - Hw','Lb - Mw','Lb - Hw'))|
    (allvalues.unadj$size<180 & allvalues.unadj$name %in% c('Mb - Lw'))|
    allvalues.unadj$name %in% c('Hb - Lw')
  
  ,]

Anth.unadj.RF <- Anth.unadj[Anth.unadj$variable=='RF',]

AuRF.means <- aggregate(Anth.unadj.RF$value,list(Anth.unadj.RF$name,Anth.unadj.RF$size),
               function(X)mean(X,na.rm=T))

AuRF.means$Classif <- 'RF'
AuRF.means$Quantif <- 'Unadjusted'

AuRF.sd <- aggregate(Anth.unadj.RF$value,list(Anth.unadj.RF$name,Anth.unadj.RF$size),
          function(X)sd(X,na.rm=T))

  # Adjusted values
      ## Let's take only the wanted metrics

Anth.adj <- allvalues.adj[
  
  allvalues.adj$size >=160 |
    (allvalues.adj$size>=90 & allvalues.adj$name %in% c('Lb - Lw','Mb - Mw',
                                                        'Hb - Mw','Hb - Hw')) |
    (allvalues.adj$size>=150 & allvalues.adj$name %in% c('Lb - Mw','Mb - Hw','Lb - Hw'))
  
  ,]


Anth.adj.RF <- Anth.adj[Anth.adj$variable=='RF',]
Anth.adj.RF1 <- Anth.adj.RF[Anth.adj.RF$name %in% c('Lb - Lw','Hb - Hw','Mb - Mw',
                                                    'Hb - Mw','Mb - Lw'),]

Anth.adj.GMM <- Anth.adj[Anth.adj$variable=='GMM',]
Anth.adj.GMM1 <- Anth.adj.GMM[Anth.adj.GMM$name %in% c('Mb - Hw','Lb - Mw','Lb - Hw'),]


AaRF.means <- aggregate(Anth.adj.RF1$value,list(Anth.adj.RF1$name,Anth.adj.RF1$size),
               function(X)mean(X,na.rm=T))
AaRF.means$Classif <- 'RF'
AaRF.means$Quantif <- 'Adjusted'

AaRF.sd <- aggregate(Anth.adj.RF1$value,list(Anth.adj.RF1$name,Anth.adj.RF1$size),
          function(X)sd(X,na.rm=T))

  
AaGMM.means <- aggregate(Anth.adj.GMM1$value,list(Anth.adj.GMM1$name,Anth.adj.GMM1$size),
               function(X)mean(X,na.rm=T))
AaGMM.means$Classif <- 'GMM'
AaGMM.means$Quantif <- 'Adjusted'

AaGMM.sd <- aggregate(Anth.adj.GMM1$value,list(Anth.adj.GMM1$name,Anth.adj.GMM1$size),
          function(X)sd(X,na.rm=T))


  # Complete Anthrone data

Anth.data <- rbind(AuRF.means,AaRF.means,AaGMM.means)
colnames(Anth.data) <- c('Variance','Size',metrics.name,'Classif','Quantif')
Anth.data$SD <- c(AuRF.sd$x,AaRF.sd$x,AaGMM.sd$x)
Anth.data$Chem <- 'Anthrone'

## Here we order the variance treatments according to the 'difficulty scale' (RDA)

Anth.data$Variance <- factor(Anth.data$Variance, levels = c('Hb - Lw', 'Mb - Lw', 'Hb - Mw',
                                                            'Lb - Lw','Mb - Mw','Hb - Hw',
                                                            'Mb - Hw','Lb - Mw','Lb - Hw'))

  #### HPLC

allvalues.unadj <- metrics.table(metrics.unadj,listChem[['HPLC']])
allvalues.unadj$Quantif <- 'Unadj'

allvalues.adj <- metrics.table(metrics.adj,listChem[['HPLC']])
allvalues.adj$Quantif <- 'Adj'

  # Unadjusted
      ## Let's take only the wanted metrics

HPLC.unadj <- allvalues.unadj[
  
  allvalues.unadj$size <90 |
    allvalues.unadj$name %in% c('Hb - Lw','Mb - Lw')
  
  ,]

HPLC.unadj.RF <- HPLC.unadj[HPLC.unadj$variable=='RF',]



HuRF.means <- aggregate(HPLC.unadj.RF$value,list(HPLC.unadj.RF$name,HPLC.unadj.RF$size),
               function(X)mean(X,na.rm=T))
HuRF.means$Classif <- 'RF'
HuRF.means$Quantif <- 'Unadjusted'

HuRF.sd <- aggregate(HPLC.unadj.RF$value,list(HPLC.unadj.RF$name,HPLC.unadj.RF$size),
                        function(X)sd(X,na.rm=T))

  # Adjusted
      ## Let's take only the wanted metrics

HPLC.adj <- allvalues.adj[
  
  allvalues.adj$size >= 90 &
    allvalues.adj$name != 'Hb - Lw' &
    allvalues.adj$name != 'Mb - Lw'
  ,]


HPLC.adj.RF <- HPLC.adj[HPLC.adj$variable=='RF',]

HaRF.means <- aggregate(HPLC.adj.RF$value,list(HPLC.adj.RF$name,HPLC.adj.RF$size),
               function(X)mean(X,na.rm=T))
HaRF.means$Classif <- 'RF'
HaRF.means$Quantif <- 'Adjusted'

HaRF.sd <- aggregate(HPLC.adj.RF$value,list(HPLC.adj.RF$name,HPLC.adj.RF$size),
                     function(X)sd(X,na.rm=T))


  # Complete HPLC data

HPLC.data <- rbind(HuRF.means,HaRF.means)
colnames(HPLC.data) <- c('Variance','Size',metrics.name,'Classif','Quantif')
HPLC.data$SD <- c(HuRF.sd$x,HaRF.sd$x)
HPLC.data$Chem <- 'HPLC'

## Here we order the variance treatments according to the 'difficulty scale' (RDA)

HPLC.data$Variance <- factor(HPLC.data$Variance, levels = c('Hb - Lw', 'Mb - Lw', 'Hb - Mw',
                                            'Lb - Lw','Mb - Mw','Hb - Hw',
                                            'Mb - Hw','Lb - Mw','Lb - Hw'))


  #### Binding all data
Chem.data <- rbind(Anth.data,HPLC.data)
pretty <- pretty(range(Chem.data[,metrics.name]))
Chem.data$Size <- as.factor(Chem.data$Size)

ggplot(data = Chem.data, aes(x=Size, y=Variance, fill=get(metrics.name)))+
  facet_wrap(~Chem,ncol=1)+
  geom_tile()+
  scale_fill_gradient2(    low = "blue", high = "red", mid = "lightblue", 
                           name=metrics.name2,
                           lim=c(pretty[1],pretty[length(pretty)]),             
                           breaks = pretty)+
  theme(axis.text.y=element_text(angle=45, hjust=1))+
  guides(fill = guide_colourbar(barheight = 25))+
  xlab('Dataset size')+
  ylab('Variance Treatment')

}


### Now we can plot the heatmaps

H <- final.heatmap(mod.metrics.Anth,mod.metrics.HPLC,'BC.unadj','BC.adj','Bray-Curtis Dissimilarity','Bray-Curtis')
H2 <- H + scale_y_discrete(position = "right")

final.heatmap(mod.metrics.Anth,mod.metrics.HPLC,'Kld.unadj','Kld.adj','Kullback-Leibler Divergence','KLD')+
scale_y_discrete(position = "right")


#ggsave("Heatmap.png", plot = H2, width = 16, height = 22,units = c("cm"), dpi = 300)
