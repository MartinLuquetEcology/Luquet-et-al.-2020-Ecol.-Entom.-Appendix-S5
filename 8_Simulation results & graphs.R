library(reshape2)
library(mgcv)
library(ggplot2)
library(gridExtra)
library(ggpubr)

#### This script was used to analyse the simulated datasets
#### Data were generated from the script "Simulation experiment"
#### All loaded data are found in the "Saved objects" folder
#### All functions used can be found in the "plotting functions for simulated data" script

  #### Importing objects ####

mod.metrics.HPLC.high <- readRDS(file='Saved objects/2019-05-20_metrics.3var.HPLC_high')
names(mod.metrics.HPLC.high) <- c('Hb - Lw','Hb - Mw','Hb - Hw')
mod.metrics.HPLC.mid <- readRDS(file='Saved objects/2019-05-27_metrics.3var.HPLC_mid')
names(mod.metrics.HPLC.mid) <- c('Mb - Lw','Mb - Mw','Mb - Hw')
mod.metrics.HPLC.low <- readRDS(file='Saved objects/2019-05-17_metrics.3var.HPLC_low')
names(mod.metrics.HPLC.low) <- c('Lb - Lw','Lb - Mw','Lb - Hw')

mod.metrics.HPLC <- c(mod.metrics.HPLC.high,mod.metrics.HPLC.mid,mod.metrics.HPLC.low)


mod.metrics.Anth.high <- readRDS(file='Saved objects/2019-05-21_metrics.3var.Anth_high')
names(mod.metrics.Anth.high) <- c('Hb - Lw','Hb - Mw','Hb - Hw')
mod.metrics.Anth.mid <- readRDS(file='Saved objects/2019-05-20_metrics.3var.Anth_mid')
names(mod.metrics.Anth.mid) <- c('Mb - Lw','Mb - Mw','Mb - Hw')
mod.metrics.Anth.low <- readRDS(file='Saved objects/2019-05-22_metrics.3var.Anth_low')
names(mod.metrics.Anth.low) <- c('Lb - Lw','Lb - Mw','Lb - Hw')

mod.metrics.Anth <- c(mod.metrics.Anth.high,mod.metrics.Anth.mid,mod.metrics.Anth.low)

sizes <- c(27,58,90,119,151,180,209,241,270)

  #### Comparing all metrics over dataset sizes ####

a <- smooth.size(mod.metrics.Anth,'Anthrone','Accuracy')+ 
  coord_cartesian(ylim=c(0.4,1))

b <- smooth.size(mod.metrics.HPLC,'HPLC','Accuracy')+ 
  coord_cartesian(ylim=c(0.4,1))+
  theme(legend.text=element_text(size=14),legend.title=element_text(size=14))

c <- smooth.size.quant(mod.metrics.Anth,'Bray-Curtis Dissimilarity','BC.unadj','BC.adj','Anthrone')+
  coord_trans(y = "log2",limy = c(0.025,0.6))+
  scale_y_continuous(minor_breaks=c(0.025,0.05,0.1,0.2 ,0.4),
                      breaks=c(0.05,0.1,0.2 ,0.4),
                      labels=c(0.05,0.1,0.2 ,0.4))

d <- smooth.size.quant(mod.metrics.HPLC,'Bray-Curtis Dissimilarity','BC.unadj','BC.adj','HPLC')+
  coord_trans(y = "log2",limy = c(0.025,0.6))+
  scale_y_continuous(minor_breaks=c(0.025,0.05,0.1,0.2 ,0.4),
                      breaks=c(0.05,0.1,0.2 ,0.4),
                      labels=c(0.05,0.1,0.2 ,0.4))

e <- smooth.size.quant(mod.metrics.Anth,'Kullback-Leibler Divergence','Kld.unadj','Kld.adj','Anthrone')+
  coord_trans(y = "log2",limy = c(0.01,3))+
  scale_y_continuous(breaks=c(0.01,0.1,1),labels=c(0.01,0.1,1),
                     minor_breaks=c(0.01,0.03,0.1,0.3,1,3))
  
f <- smooth.size.quant(mod.metrics.HPLC,'Kullback-Leibler Divergence','Kld.unadj','Kld.adj','HPLC')+
  coord_trans(y = "log2",limy = c(0.01,3))+
  scale_y_continuous(breaks=c(0.01,0.1,1),labels=c(0.01,0.1,1),
                      minor_breaks=c(0.01,0.03,0.1,0.3,1,3))


i <- ggarrange(a,b, nrow=1, ncol=2,common.legend=T,legend="left")
j <- ggarrange(c,d,e,f, nrow=2, ncol=2,common.legend=T,legend="left")

#ggsave("Acc_size.png", plot = i, width = 22, height = 10,units = c("cm"),       dpi = 300)

#ggsave("Quantif_size.png",plot = j, width = 22, height = 20,units = c("cm"),       dpi = 300)


  #### Comparing metrics over dataset size and variance ####


  ## Some boxplot over all sizes for all variances

BC.adj.Anth <- metrics.table('BC.adj',mod.metrics.Anth,rm=T)

ggplot(BC.adj.Anth , aes( y=value, colour=variable)) +
  facet_wrap(~name)+
  geom_boxplot(outlier.size = 0.5)

BC.adj.HPLC <- metrics.table('BC.adj',mod.metrics.HPLC,rm=T)

ggplot(BC.adj.HPLC, aes( y=value, colour=variable)) +
  facet_wrap(~name)+
  geom_boxplot(outlier.size = 0.5)

Kld.adj.Anth <- metrics.table('Kld.adj',mod.metrics.Anth,rm=T)

ggplot(Kld.adj.Anth , aes( y=value, colour=variable)) +
  facet_wrap(~name)+
  geom_boxplot(outlier.size = 0.5)

  ## Let's look at size-performances across all variances
  ## With and without data - adding data allows seeing variability quite nicely

# Accuracy

all.acc.anth <- smooth.size(mod.metrics.Anth,'Anthrone','Accuracy')+
  facet_wrap(~name)+  theme(legend.position='bottom')

all.acc.hplc <- smooth.size(mod.metrics.HPLC,'HPLC','Accuracy',k=8)+
  facet_wrap(~name)+  theme(legend.position='bottom')

smooth.size(mod.metrics.HPLC,'HPLC','Accuracy',k=8)+
  facet_wrap(~name)+
  geom_jitter(shape='.')

smooth.size(mod.metrics.Anth,'Anthrone','Accuracy',k=8)+
  facet_wrap(~name)+
  geom_jitter(shape='.')

#ggsave('all.acc.anth.png',all.acc.anth,width=14,height=16,unit='cm')
#ggsave('all.acc.HPLC.png',all.acc.hplc,width=14,height=16,unit='cm')

# It seems estimations are more test dataset-dependent in HPLC (clumps of data)

# Bray-Curtis


u <- smooth.size.quant(mod.metrics.Anth,'Bray-Curtis Dissimilarity','BC.unadj','BC.adj',
                       'Anthrone')+
  facet_wrap(~name)+
  theme(strip.text=element_text(size=14),legend.position='bottom')

all.bc.hplc <- smooth.size.quant(mod.metrics.HPLC,'Bray-Curtis Dissimilarity','BC.unadj','BC.adj',
                  'HPLC',k=8,rt=3)+
  facet_wrap(~name)+
  theme(legend.text=element_text(size=7),legend.position='bottom',
        legend.title=element_text(size=7),
        legend.key.size = unit(0.4, "cm"))

smooth.size.quant(mod.metrics.HPLC,'Bray-Curtis Dissimilarity','BC.unadj','BC.adj',
                  'HPLC',k=8,rt=1)+
  facet_wrap(~name)+
  coord_cartesian(ylim=c(0,0.3))

smooth.size(mod.metrics.HPLC,'HPLC','BC.adj','Bray-Curtis adjusted',k=8,rm=T,rt=3)+
  facet_wrap(~name)+
  geom_jitter(shape='.')

smooth.size(mod.metrics.HPLC,'HPLC','BC.unadj','Bray-Curtis unadjusted',k=8,rm=T)+
  facet_wrap(~name)+
  geom_jitter(shape='.')

smooth.size(mod.metrics.Anth,'Anthrone','BC.adj','Bray-Curtis adjusted',k=8,rm=T)+
  facet_wrap(~name)+
  geom_jitter(shape='.')

smooth.size(mod.metrics.Anth,'Anthrone','BC.unadj','Bray-Curtis unadjusted',k=8,rm=T)+
  facet_wrap(~name)+
  geom_jitter(shape='.')


#ggsave("Bray-Curtis anthrone according to variance.png",plot = u,        width = 22, height = 24,units = c("cm"), dpi = 300)


# KLD

smooth.size.quant(mod.metrics.Anth,'Kullback-Leibler Divergence','Kld.unadj','Kld.adj',
                  'Anthrone')+
  facet_wrap(~name)+
  coord_cartesian(ylim=c(0,0.5))

smooth.size.quant(mod.metrics.HPLC,'Kullback-Leibler Divergence','Kld.unadj','Kld.adj',
                  'HPLC',k=8,rt=2)+
  facet_wrap(~name)+
  coord_cartesian(ylim=c(0,0.5))

smooth.size(mod.metrics.HPLC,'HPLC','Kld.adj','KLD adjusted',k=8,rm=T)+
  facet_wrap(~name)+
  geom_jitter(shape='.')+
  coord_cartesian(ylim=c(0,1))

smooth.size.quant(mod.metrics.HPLC,'Kullback-Leibler Divergence','Kld.unadj','Kld.adj','HPLC',k=5,rt=1)+
  facet_wrap(~name)+
  coord_cartesian(ylim=c(0,0.4))

smooth.size.quant(mod.metrics.Anth,'Kullback-Leibler Divergence','Kld.unadj','Kld.adj','Anthrone',k=5,rt=1)+
  facet_wrap(~name)+
  coord_cartesian(ylim=c(0,0.4))


all.kld.anth1 <- smooth.size.quant(mod.metrics.Anth,'Kullback-Leibler Divergence','Kld.unadj','Kld.adj',
                  'Anthrone')+
  facet_wrap(~name)+
  theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=7),legend.position='bottom',
        legend.title=element_text(size=7),
        legend.key.size = unit(0.4, "cm"))


all.kld.anth2 <- smooth.size.quant(mod.metrics.Anth,'Kullback-Leibler Divergence','Kld.unadj','Kld.adj',
                                   'Anthrone')+
  facet_wrap(~name)+
  coord_cartesian(ylim=c(0,0.5))+
  labs(subtitle='0 to 0.5 scale')+
  theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=7),legend.position='bottom',
        legend.title=element_text(size=7),
        legend.key.size = unit(0.4, "cm"))


all.kld.hplc1 <- smooth.size.quant(mod.metrics.HPLC,'Kullback-Leibler Divergence','Kld.unadj','Kld.adj',
                  'HPLC',k=8,rt=3)+
  facet_wrap(~name)+
  theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=7),legend.position='bottom',
        legend.title=element_text(size=7),
        legend.key.size = unit(0.4, "cm"))


all.kld.hplc2 <- smooth.size.quant(mod.metrics.HPLC,'Kullback-Leibler Divergence','Kld.unadj','Kld.adj',
                  'HPLC',k=8,rt=3)+
  facet_wrap(~name)+
  coord_cartesian(ylim=c(0,0.5))+
  labs(subtitle='0 to 0.5 scale')+
  theme(legend.position = 'bottom')+
  theme(legend.text=element_text(size=7),legend.position='bottom',
        legend.title=element_text(size=7),
        legend.key.size = unit(0.4, "cm"))


#Saving objects
#ggsave('all.bc.hplc.png',all.bc.hplc,height=16,width=14,unit='cm')
#ggsave('all.kld.hplc1.png',all.kld.hplc1,height=16,width=14,unit='cm')
#ggsave('all.kld.hplc2.png',all.kld.hplc2,height=16,width=14,unit='cm')
#ggsave('all.kld.anth1.png',all.kld.anth1,height=16,width=14,unit='cm')
#ggsave('all.kld.anth2.png',all.kld.anth2,height=16,width=14,unit='cm')

  ##### To see better metrics by metrics, we can draw
  ##### some boxplots

for(n in names(mod.metrics.Anth ))
{
  adj.vs.unadj(n,'BC.unadj','BC.adj',mod.metrics.Anth,'Anthrone')
}

for(n in names(mod.metrics.Anth ))
{
  adj.vs.unadj(n,'Kld.unadj','Kld.adj',mod.metrics.Anth,'Anthrone')
}

for(n in names(mod.metrics.HPLC ))
{
  adj.vs.unadj(n,'BC.unadj','BC.adj',mod.metrics.HPLC,'Anthrone')
}

for(n in names(mod.metrics.HPLC ))
{
  adj.vs.unadj(n,'Kld.unadj','Kld.adj',mod.metrics.HPLC,'HPLC')
}


#### Looking more closely, metrics by metrics ####

  ## Let's take some examples: RF and DA

BC.RF.HPLC <- metrics.table.quant(mod.metrics.HPLC,'BC.unadj','BC.adj',
                                   'HPLC','RF',rm=T)

BC.DA.HPLC <- metrics.table.quant(mod.metrics.HPLC,'BC.unadj','BC.adj',
                                  'HPLC','DA',rm=T)

BC.RF.Anth <- metrics.table.quant(mod.metrics.Anth,'BC.unadj','BC.adj',
                                  'Anthrone','RF',rm=T)

BC.GMM.Anth <- metrics.table.quant(mod.metrics.Anth,'BC.unadj','BC.adj',
                                  'Anthrone','GMM',rm=T)

Kld.RF.HPLC <- metrics.table.quant(mod.metrics.HPLC,'Kld.unadj','Kld.adj',
                                   'HPLC','RF',rm=T)

Kld.DA.HPLC <- metrics.table.quant(mod.metrics.HPLC,'Kld.unadj','Kld.adj',
                                   'HPLC','DA',rm=T)

Kld.RF.Anth <- metrics.table.quant(mod.metrics.Anth,'Kld.unadj','Kld.adj',
                                   'Anthrone','RF',rm=T)

Kld.GMM.Anth <- metrics.table.quant(mod.metrics.Anth,'Kld.unadj','Kld.adj',
                                   'Anthrone','GMM',rm=T)

Kld.DA.Anth <- metrics.table.quant(mod.metrics.Anth,'Kld.unadj','Kld.adj',
                                   'Anthrone','DA',rm=T)

Kld.GMM.HPLC <- metrics.table.quant(mod.metrics.HPLC,'Kld.unadj','Kld.adj',
                                    'HPLC','GMM',rm=T)


# 1: boxplots

ggplot(BC.RF.HPLC, aes(x=as.factor(size), y=value, fill=Quantif)) +
  facet_wrap(~name)+
  geom_boxplot(outlier.size=0.5)+
  coord_cartesian(ylim=c(0,0.6))+
  ggtitle('RF - HPLC')

# Looking at BC we have better peformances with adjusted counting -> reduced mean & variance
# Looking at KLD, it's more or less the same, but there are important outliers
# --> See below for further investigation
ggplot(Kld.RF.HPLC, aes(x=as.factor(size), y=value, fill=Prev.est)) +
  facet_wrap(~name)+
  geom_boxplot(outlier.size=0.5)+
  coord_cartesian(ylim=c(0,1))+
  ggtitle('RF - HPLC')

ggplot(Kld.RF.Anth, aes(x=as.factor(size), y=value, fill=Prev.est)) +
  facet_wrap(~name)+
  geom_boxplot(outlier.size=0.5)+
  coord_cartesian(ylim=c(0,1))+
  ggtitle('RF - Anthrone')

# Here we can see this doesn't work well at Var2 & Var3
# No size-performance relationship for those variances
ggplot(Kld.DA.Anth, aes(x=as.factor(size), y=value, fill=Prev.est)) +
  facet_wrap(~name)+
  geom_boxplot(outlier.size=0.5)+
  coord_cartesian(ylim=c(0,0.5))+
  ggtitle('DA - Anthrone')

## Good here !!
ggplot(Kld.GMM.Anth, aes(x=as.factor(size), y=value, fill=Prev.est)) +
  facet_wrap(~name)+
  geom_boxplot(outlier.size=0.5)+
  coord_cartesian(ylim=c(0,1))+
  ggtitle('GMM - Anthrone')

# 2: smoothers, with and without data

ggplot(Kld.RF.HPLC, aes(x=size, y=value, fill=Prev.est,colour=Prev.est)) +
  facet_wrap(~name)+
  geom_smooth(alpha=.2,method='gam',formula=y ~ s(x, k=9,bs = "cs"))+
  ggtitle('Bray-Curtis Dissimilarity - HPLC - Random Forest')

# We can have a nice look at the variability with this graph:
z <- ggplot(BC.RF.HPLC, aes(x=size, y=value, colour=Prev.est)) +
  facet_wrap(~name)+
  ggtitle('Bray-Curtis Dissimilarity - HPLC - Random Forest')+
  geom_jitter(shape='.')+
  geom_smooth(method='gam',formula=y ~ s(x, k=9,bs = "cs"))+
  coord_cartesian(ylim=c(0,0.5))

#ggsave("Bray-Curtis HPLC adj vs unadj.png",plot = z,        width = 22, height = 20,units = c("cm"), dpi = 300)

#Now looking at KLD:

ggplot(Kld.RF.HPLC, aes(x=size, y=value, colour=Prev.est)) +
  facet_wrap(~name)+
  ggtitle('Kullback-Leibler Divergence - HPLC - Random Forest')+
  geom_jitter(shape='.')+
  geom_smooth(method='gam',formula=y ~ s(x, k=9,bs = "cs"))+
  coord_cartesian(ylim=c(0,0.5))

# 3: we can also compare HPLC and Anthrone

BC.RF <- metrics.table.all(mod.metrics.Anth,mod.metrics.HPLC,
                       'BC.unadj','BC.adj',
                       'Anthrone','HPLC',
                       Classif='RF',rm=T)

ggplot(BC.RF, aes(x=size, y=value, colour=chem, lty=Prev.est)) +
  facet_wrap(~name)+
  ggtitle('Bray-Curtis Dissimilarity - HPLC - Random Forest')+
  geom_smooth(method='gam',formula=y ~ s(x, k=9,bs = "cs"))

BC.Kld <- metrics.table.all(mod.metrics.Anth,mod.metrics.HPLC,
                                'Kld.unadj','Kld.adj',
                                'Anthrone','HPLC',
                                Classif='RF',rm=T)

ggplot(BC.Kld, aes(x=size, y=value, colour=chem, lty=Prev.est)) +
  facet_wrap(~name)+
  ggtitle('Kullback-Leibler Divergence - HPLC - Random Forest')+
  geom_smooth(method='gam',formula=y ~ s(x, k=9,bs = "cs"))+
  coord_cartesian(ylim=c(0,0.5))

## Very interesting: at some variances, adjustement doesn't work for DA
## in Anthrone, while it works for RF
## (See next section)

BC.DA <- metrics.table.all(mod.metrics.Anth,mod.metrics.HPLC,
                           'BC.unadj','BC.adj',
                           'Anthrone','HPLC',
                           Classif='DA',rm=T)

ggplot(BC.DA, aes(x=size, y=value, colour=chem, lty=Prev.est)) +
  facet_wrap(~name)+
  ggtitle('Bray-Curtis Dissimilarity - HPLC - Discriminant Analysis')+
  geom_smooth(method='gam',formula=y ~ s(x, k=9,bs = "cs"))

Kld.DA <- metrics.table.all(mod.metrics.Anth,mod.metrics.HPLC,
                           'Kld.unadj','Kld.adj',
                           'Anthrone','HPLC',
                           Classif='RF',rm=T)

ggplot(Kld.DA, aes(x=size, y=value, colour=chem, lty=Prev.est)) +
  facet_wrap(~name)+
  ggtitle('KLD - HPLC - Discriminant Analysis')+
  geom_smooth(method='gam',formula=y ~ s(x, k=9,bs = "cs"))+
  coord_cartesian(ylim=c(0,0.5))

#### Looking a bit more closely at the adjusted/unadjusted relationship ####


##we can evaluate performance according to test dataset distribution

BC.RF.HPLC <- metrics.table.quant(mod.metrics.HPLC,'BC.unadj','BC.adj','HPLC',
                         Classif='RF',rm=T)

Kld.RF.HPLC <- metrics.table.quant(mod.metrics.HPLC,'Kld.unadj','Kld.adj','HPLC',
                                  Classif='RF',rm=T)

# All variances, over all sizes

ggplot(BC.RF.HPLC, aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~name)+
  ggtitle('Bray-Curtis Dissimilarity - HPLC - Random Forest')+
  geom_jitter(shape='.')+
  geom_smooth(method='gam',formula=y ~ s(x, k=8,bs = "cs"))

ggplot(Kld.RF.HPLC , aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~name)+
  ggtitle('KLD - HPLC - Random Forest')+
  geom_jitter(shape='.')+
  geom_smooth(method='gam',formula=y ~ s(x, k=8,bs = "cs"))+
  coord_cartesian(ylim=c(0,0.4))

# Now if we look at only one variance treatment, size by size
# e.g. Var3 in which there is a lot of variability in KLD


## note: funny to just look at Var3 but for all Classifiers
x2 <-  BC.RF.HPLC [BC.RF.HPLC $name=='Lb - Hw',]
x3 <-  Kld.RF.HPLC [Kld.RF.HPLC $name=='Lb - Hw',]
x4 <- BC.RF.Anth [BC.RF.Anth $name=='Lb - Hw',]
x5 <- Kld.RF.Anth [Kld.RF.Anth $name=='Lb - Hw',]
x6 <- Kld.GMM.Anth [Kld.GMM.Anth $name=='Lb - Hw',]
x7 <- BC.GMM.Anth [BC.GMM.Anth $name=='Lb - Hw',]
x8 <- BC.DA.HPLC [BC.DA.HPLC $name=='Lb - Hw',]
x9 <- Kld.DA.HPLC [Kld.DA.HPLC $name=='Lb - Hw',]
x10 <- Kld.GMM.HPLC [Kld.GMM.HPLC $name=='Lb - Hw',]

# Here we can see than generally -> lower mean and lower
# Inter-population variance

m <- ggplot(x2, aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('Bray-Curtis Dissimilarity - HPLC - Random Forest')+
  geom_jitter(size=0.9,alpha=1)+
  geom_smooth(method='gam',formula=y ~ s(x, k=8,bs = "cs"))+
  scale_x_continuous(labels=seq(1,12),breaks=seq(1,12))+
  theme(legend.position='bottom')+
  xlab('Field dataset')+
  coord_cartesian(ylim=c(0,0.4))


ggplot(x3, aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('Kld - HPLC - Random Forest')+
  geom_jitter(size=0.9,alpha=0.5)+
  geom_smooth(method='gam',formula=y ~ s(x, k=8,bs = "cs"))+
  scale_x_continuous(labels=seq(1,12),breaks=seq(1,12))+
  theme(legend.position='bottom')+
  xlab('Field dataset')+
  coord_cartesian(ylim=c(0,0.5))

#ggsave("Adj vs unadj - size.png",plot = m, width = 22, height = 20,units = c("cm"),       dpi = 300)

ggplot(x4, aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('Bray-Curtis Dissimilarity - Anthrone - Random Forest')+
  geom_jitter(size=0.1,alpha=0.5)+
  geom_smooth(method='gam',formula=y ~ s(x, k=8,bs = "cs"))+
  coord_cartesian(ylim=c(0,0.5))

ggplot(x5, aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('KLD - Anthrone - Random Forest')+
  geom_jitter(size=0.1,alpha=0.5)+
  geom_smooth(method='gam',formula=y ~ s(x, k=8,bs = "cs"))+
  coord_cartesian(ylim=c(0,1))

ggplot(x6, aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('KLD - Anthrone - Gaussian Mixture Models')+
  geom_jitter(size=0.1,alpha=0.5)+
  geom_smooth(method='gam',formula=y ~ s(x, k=8,bs = "cs"))+
  coord_cartesian(ylim=c(0,1))

ggplot(x2, aes(x=as.factor(L1), y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  geom_boxplot(outlier.size = 0.01,size=1)+
  ggtitle('KLD - HPLC - Random Forest')+
  geom_jitter(size=0.1,alpha=0.5)+
  coord_cartesian(ylim=c(0,0.5))

ggplot(x3, aes(x=as.factor(L1), y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  geom_boxplot(outlier.size = 0.01,size=1)+
  ggtitle('KLD - HPLC - Random Forest')+
  geom_jitter(size=0.1,alpha=0.5)+
  coord_cartesian(ylim=c(0,0.5))

ggplot(x4, aes(x=as.factor(L1), y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  geom_boxplot(outlier.size = 0.01,size=1)+
  ggtitle('Bray-Curtis - Anthrone - Random Forest')+
  geom_jitter(size=0.1,alpha=0.5)+
  coord_cartesian(ylim=c(0,0.5))

ggplot(x5, aes(x=as.factor(L1), y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  geom_boxplot(outlier.size = 0.01,size=1)+
  ggtitle('KLD - Anthrone - Random Forest')+
  geom_jitter(size=0.1,alpha=0.5)+
  coord_cartesian(ylim=c(0,0.5))

ggplot(x6, aes(x=as.factor(L1), y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  geom_boxplot(outlier.size = 0.01,size=1)+
  ggtitle('KLD - Anthrone - Gaussian Mixture Models')+
  geom_jitter(size=0.1,alpha=0.5)+
  coord_cartesian(ylim=c(0,0.5))


## Also interesting to see it that way :

ggplot(BC.RF.HPLC, aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('Bray-Curtis Dissimilarity - HPLC - Random Forest')+
  geom_jitter(size=0.1,alpha=0.5)+
  geom_boxplot(outlier.size=0.5)+
  coord_cartesian(ylim=c(0,0.5))

ggplot(Kld.RF.HPLC, aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('Bray-Curtis Dissimilarity - HPLC - Random Forest')+
  geom_jitter(size=0.1,alpha=0.5)+
  geom_boxplot(outlier.size=0.5)+
  coord_cartesian(ylim=c(0,0.5))

ggplot(BC.RF.Anth, aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('Bray-Curtis Dissimilarity - HPLC - Random Forest')+
  geom_jitter(size=0.1,alpha=0.5)+
  geom_boxplot(outlier.size=0.5)+
  coord_cartesian(ylim=c(0,0.5))

ggplot(BC.RF.Anth, aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('Bray-Curtis Dissimilarity - Anthrone - Random Forest')+
  geom_jitter(size=0.1,alpha=0.5)+
  geom_boxplot(outlier.size=0.5)+
  coord_cartesian(ylim=c(0,0.5))

ggplot(Kld.RF.Anth, aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('KLD - Anthrone - Random Forest')+
  geom_jitter(size=0.1,alpha=0.5)+
  geom_boxplot(outlier.shape=NA,size=1)+
  coord_cartesian(ylim=c(0,0.5))

ggplot(x6, aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('KLD - Anthrone - Random Forest')+
  geom_boxplot(outlier.shape=NA,size=1)+
  geom_jitter(size=0.5,alpha=0.5)+
  coord_cartesian(ylim=c(0,0.5))


p1 <- ggplot(x4[x4$size>190,], aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('BC - Anthrone - RF')+
  geom_jitter(size=0.9,alpha=0.5)+
  geom_smooth(method='gam',formula=y ~ s(x, k=8,bs = "cs"))+
  coord_cartesian(ylim=c(0,0.7))

p2 <- ggplot(x7[x7$size>190,], aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('BC - Anthrone - GMM')+
  geom_jitter(size=0.9,alpha=0.5)+
  geom_smooth(method='gam',formula=y ~ s(x, k=8,bs = "cs"))+
  coord_cartesian(ylim=c(0,0.7))

p3 <- ggplot(x5[x5$size>190,], aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('Kld - Anthrone - RF')+
  geom_jitter(size=0.9,alpha=0.5)+
  geom_smooth(method='gam',formula=y ~ s(x, k=8,bs = "cs"))+
  coord_cartesian(ylim=c(0,1))

p4 <- ggplot(x6[x6$size>190,], aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('Kld - Anthrone - GMM')+
  geom_jitter(size=0.9,alpha=0.5)+
  geom_smooth(method='gam',formula=y ~ s(x, k=8,bs = "cs"))+
  coord_cartesian(ylim=c(0,1))


ggarrange(p1,p2,p3,p4,common.legend=T)

ggplot(x4[x4$size>190,], aes(x=L1, y=value,colour=Prev.est)) +
  facet_wrap(~size)+
  ggtitle('BC - Anthrone - RF')+
  geom_boxplot(outlier.shape=NA,size=1)+
  geom_jitter(size=0.9,alpha=0.5)+
  coord_cartesian(ylim=c(0,0.7))



#### Analyses of performance variance ####

# Over all variance treatments:

var.plot(mod.metrics.Anth,"Accuracy",'Accuracy','Anthrone')
var.plot(mod.metrics.Anth,"BC.unadj",'Bray-Curtis unadjusted','Anthrone')
var.plot(mod.metrics.Anth,"BC.adj",'Bray-Curtis unadjusted','Anthrone')
var.plot(mod.metrics.Anth,"Kld.unadj",'Kld unadjusted','Anthrone')+
  coord_trans(y = "log2",limy = c(0.1,1.5))
var.plot(mod.metrics.Anth,"Kld.adj",'Kld unadjusted','Anthrone')+
  coord_trans(y = "log2",limy = c(0.1,1.5))

var.plot(mod.metrics.HPLC,"Accuracy",'Accuracy','HPLC')
var.plot(mod.metrics.HPLC,"BC.unadj",'Bray-Curtis unadjusted','HPLC')
var.plot(mod.metrics.HPLC,"BC.adj",'Bray-Curtis unadjusted','HPLC')
var.plot(mod.metrics.HPLC,"Kld.unadj",'Kld unadjusted','HPLC')+
  coord_trans(y = "log2",limy = c(0.05,3))
var.plot(mod.metrics.HPLC,"Kld.adj",'Kld adjusted','HPLC')+
  coord_trans(y = "log2",limy = c(0.05,3))

# On each variance treatment:

  ##Accuracy
var.plot2(mod.metrics.Anth,'Accuracy','Accuracy','Anthrone')
var.plot2(mod.metrics.HPLC,'Accuracy','Accuracy','HPLC')

  ## Bray-Curtis
var.plot3(mod.metrics.Anth,'BC.unadj','BC.adj','Bray-Curtis Dissimilarity','Anthrone')
var.plot3(mod.metrics.HPLC,'BC.unadj','BC.adj','Bray-Curtis Dissimilarity','HPLC')

  ## KLD
# First over the whole range of values, then on the 0 to 0.5 scale

var.plot3(mod.metrics.HPLC,"Kld.adj","Kld.unadj",'Kullback-Leibler Divergence','HPLC')
var.plot3(mod.metrics.HPLC,"Kld.adj","Kld.unadj",'Kullback-Leibler Divergence','HPLC')+
  coord_cartesian(ylim=c(0,0.5))+
  labs(subtitle='0 to 0.5 scale')

var.plot3(mod.metrics.Anth,"Kld.adj","Kld.unadj",'Kullback-Leibler Divergence','Anthrone')
var.plot3(mod.metrics.Anth,"Kld.adj","Kld.unadj",'Kullback-Leibler Divergence','Anthrone')+
  coord_cartesian(ylim=c(0,0.5))+
  labs(subtitle='0 to 0.5 scale')
## It's interesting to note that Adjusted Counting always allows reducing KLD variance, particularly on small size datasets

#### Number of times when AC prevalence estimation is not computable ####

#### HPLC

BC.Thresh.HPLC <- metrics.table.quant(mod.metrics.HPLC,'BC.unadj','BC.adj',
                                      'HPLC','Thresh',rm=F)

length(BC.Thresh.HPLC$value)

length(BC.Thresh.HPLC[BC.Thresh.HPLC$Quantif=='Unadjusted','value'])
sum(is.na(BC.Thresh.HPLC[BC.Thresh.HPLC$Quantif=='Unadjusted','value']))
length(BC.Thresh.HPLC[BC.Thresh.HPLC$Quantif=='Adjusted','value'])
sum(is.na(BC.Thresh.HPLC[BC.Thresh.HPLC$Quantif=='Adjusted','value']))


BC.RF.HPLC <- metrics.table.quant(mod.metrics.HPLC,'BC.unadj','BC.adj',
                                  'HPLC','RF',rm=F)

length(BC.RF.HPLC$value)

length(BC.RF.HPLC[BC.RF.HPLC$Quantif=='Unadjusted','value'])
sum(is.na(BC.RF.HPLC[BC.RF.HPLC$Quantif=='Unadjusted','value']))
length(BC.RF.HPLC[BC.RF.HPLC$Quantif=='Adjusted','value'])
sum(is.na(BC.RF.HPLC[BC.RF.HPLC$Quantif=='Adjusted','value']))


BC.DA.HPLC <- metrics.table.quant(mod.metrics.HPLC,'BC.unadj','BC.adj',
                                  'HPLC','DA',rm=F)

length(BC.DA.HPLC$value)

length(BC.DA.HPLC[BC.DA.HPLC$Quantif=='Unadjusted','value'])
sum(is.na(BC.DA.HPLC[BC.DA.HPLC$Quantif=='Unadjusted','value']))
length(BC.DA.HPLC[BC.DA.HPLC$Quantif=='Adjusted','value'])
sum(is.na(BC.DA.HPLC[BC.DA.HPLC$Quantif=='Adjusted','value']))

BC.GMM.HPLC <- metrics.table.quant(mod.metrics.HPLC,'BC.unadj','BC.adj',
                                   'HPLC','GMM',rm=F)

length(BC.GMM.HPLC$value)

length(BC.GMM.HPLC[BC.GMM.HPLC$Quantif=='Unadjusted','value'])
sum(is.na(BC.GMM.HPLC[BC.GMM.HPLC$Quantif=='Unadjusted','value']))
length(BC.GMM.HPLC[BC.GMM.HPLC$Quantif=='Adjusted','value'])
sum(is.na(BC.GMM.HPLC[BC.GMM.HPLC$Quantif=='Adjusted','value']))


#### Anthrone

BC.Thresh.Anth <- metrics.table.quant(mod.metrics.Anth,'BC.unadj','BC.adj',
                                      'Anth','Thresh',rm=F)

length(BC.Thresh.Anth$value)

length(BC.Thresh.Anth[BC.Thresh.Anth$Quantif=='Unadjusted','value'])
sum(is.na(BC.Thresh.Anth[BC.Thresh.Anth$Quantif=='Unadjusted','value']))
length(BC.Thresh.Anth[BC.Thresh.Anth$Quantif=='Adjusted','value'])
sum(is.na(BC.Thresh.Anth[BC.Thresh.Anth$Quantif=='Adjusted','value']))


BC.RF.Anth <- metrics.table.quant(mod.metrics.Anth,'BC.unadj','BC.adj',
                                  'Anth','RF',rm=F)

length(BC.RF.Anth$value)

length(BC.RF.Anth[BC.RF.Anth$Quantif=='Unadjusted','value'])
sum(is.na(BC.RF.Anth[BC.RF.Anth$Quantif=='Unadjusted','value']))
length(BC.RF.Anth[BC.RF.Anth$Quantif=='Adjusted','value'])
sum(is.na(BC.RF.Anth[BC.RF.Anth$Quantif=='Adjusted','value']))


BC.DA.Anth <- metrics.table.quant(mod.metrics.Anth,'BC.unadj','BC.adj',
                                  'Anth','DA',rm=F)

length(BC.DA.Anth$value)

length(BC.DA.Anth[BC.DA.Anth$Quantif=='Unadjusted','value'])
sum(is.na(BC.DA.Anth[BC.DA.Anth$Quantif=='Unadjusted','value']))
length(BC.DA.Anth[BC.DA.Anth$Quantif=='Adjusted','value'])
sum(is.na(BC.DA.Anth[BC.DA.Anth$Quantif=='Adjusted','value']))

BC.GMM.Anth <- metrics.table.quant(mod.metrics.Anth,'BC.unadj','BC.adj',
                                   'Anth','GMM',rm=F)

length(BC.GMM.Anth$value)

length(BC.GMM.Anth[BC.GMM.Anth$Quantif=='Unadjusted','value'])
sum(is.na(BC.GMM.Anth[BC.GMM.Anth$Quantif=='Unadjusted','value']))
length(BC.GMM.Anth[BC.GMM.Anth$Quantif=='Adjusted','value'])
sum(is.na(BC.GMM.Anth[BC.GMM.Anth$Quantif=='Adjusted','value']))