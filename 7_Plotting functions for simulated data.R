#### Here are some functions to handle, represent and analyse
#### the results from the simulation analyses


#######################################

#### metrics.table.var ####
### Function to pass from a list to a data.frame
### For one variance treatment

metrics.table.var <- function(var,metrics,chem)
{
  var.metrics.table <- data.frame(n=seq(1:27000),variable=NA,value=NA,L1=NA,size=NA,name=var)
  
  mod.metrics <- chem
  
  for (i in names(mod.metrics[[var]]))
  {
    Var1.acc <- lapply(mod.metrics[[var]][[i]], get, x=metrics)
    table <- melt(Var1.acc,id.vars=NULL)
    
    n <- which(names(mod.metrics[[var]])==i)
    nrow <- seq(from=(n-1)*3000+1,to=n*3000,by=1)
    
    var.metrics.table[nrow,2:4] <- table
    var.metrics.table[nrow,'variable'] <- as.character(levels(table$variable))[table$variable]
    
    var.metrics.table[nrow,'size'] <- sizes[n]
    
  }
  
  var.metrics.table <- var.metrics.table[!var.metrics.table$variable == 'set',]
  return(var.metrics.table)
  
}


#######################################

#### metrics.table ####
### Function to pass from a list to a data.frame
### For one metrics across all variance treatments

metrics.table <- function(metrics,chem,rm=F)
{
  
  metrics <- do.call(rbind,
                           lapply(names(chem), function(X){
                             metrics.table.var(X,metrics,chem) })
       
  )
  
  if(rm ==T) {metrics <- na.omit(metrics)}
  
  return(metrics)
  
}

#######################################

#### metrics.table.quant ####
### Function to pass from a list to a data.frame
### Allowing to compare adjusted vs unadjusted

metrics.table.quant <- function(chem,metrics1,metrics2,chemname,
                                Classifier=NULL,rm=F)
  
{
  metrics.unadj <- metrics.table(metrics1,chem)
  
  metrics.adj <-  metrics.table(metrics2,chem)
  
  metrics.unadj$Prev.est <- 'Unadjusted (CC)'
  metrics.adj$Prev.est  <- 'Adjusted (AC)'
  
  metrics <- rbind(metrics.unadj,metrics.adj)
  metrics$chem <- chemname
  
  if(!is.null(Classifier)){metrics <- metrics[metrics$variable==Classifier,]}
  
  if(rm==T) {metrics <- na.omit(metrics)}
  
  return(metrics)
}

#######################################

metrics.table.all <- function(chem1,chem2,
                              metrics1,metrics2,
                              chemname1,chemname2,
                                Classifier=NULL,rm=F)
  
{
  table.metrics1 <- metrics.table.quant(chem1,metrics1,metrics2,
                                  chemname1,Classifier,rm)
  
  table.metrics2 <- metrics.table.quant(chem2,metrics1,metrics2,
                                  chemname2,Classifier,rm)
  
  metrics <- rbind(table.metrics1,table.metrics2)
  
  metrics$Prev.est <- as.factor(metrics$Prev.est)
  metrics$Prev.est <- relevel(metrics$Prev.est, "Unadjusted (CC)")
  
  return(metrics)
}

#######################################

#### smooth.size ####
#### Plot estimated smoothed conditional means for a given metrics
#### According to dataset size

smooth.size <- function(chem,chemname,metrics1,metrics.name=NULL,k=9,rm=T,rt=0)
  
{
  
  if(is.null(metrics.name)){metrics.name <- metrics1}
  
  metrics <- metrics.table(metrics1,chem,rm)
  
  colnames(metrics)[2] <- 'Classifier'
  
  if(rt==1){
    metrics <- metrics[-which(metrics$Classifier=='Thresh'),]
  }
  
  if(rt==2){
    metrics <- metrics[-which(metrics$Classifier=='Thresh'& metrics$name=='Var3'),]
  }

  
  p <- ggplot(metrics, aes(x=size, y=value, colour=Classifier)) +
    geom_smooth(alpha=.2, size=1,method='gam',formula=y ~ s(x, k=k,bs = "cs"),aes(fill=Classifier)) +
    ggtitle(paste(metrics.name,'-',chemname))+
    scale_x_continuous(breaks=c(0,60,120,180,240))+
    scale_fill_manual(values=c(DA="#D55E00",GMM="#009E73",
                               RF="#56B4E9",Thresh="#CC79A7"))+
    scale_color_manual(values=c(DA="#D55E00",GMM="#009E73",
                                RF="#56B4E9",Thresh="#CC79A7"))+
    theme(legend.text=element_text(size=14),legend.title=element_text(size=14))+
    xlab("Dataset size") + ylab(metrics.name)
  
  plot(p)
}

#######################################

#### smooth.size ####
#### Plot estimated smoothed conditional means for quantification
#### According to dataset size

smooth.size.quant <- function(chem,metrics.name,metrics1,metrics2,chemname,k=9,
                              rt =0,lsize=14)
  
{
  
  metrics.unadj <- do.call(rbind,
                           lapply(names(chem), function(X){
                             metrics.table.var(X,metrics1,chem) })
  )
  
  metrics.adj <- do.call(rbind,
                         lapply(names(chem), function(X){
                           metrics.table.var(X,metrics2,chem) })
  )
  
  metrics.unadj$Prev.est <- 'Unadjusted (CC)'
  metrics.adj$Prev.est  <- 'Adjusted (AC)'
  
  metrics <- rbind(metrics.unadj,metrics.adj)
  metrics$Prev.est <- as.factor(metrics$Prev.est)
  metrics$Prev.est <- relevel(metrics$Prev.est, "Unadjusted (CC)")
  
  metrics <- na.omit(metrics)
  colnames(metrics)[2] <- 'Classifier'
  
  if(rt==1){
    metrics <- metrics[-which(metrics$Classifier=='Thresh'),]
  }
  
  if(rt==2){
    metrics <- metrics[-which(metrics$Prev.est=='Adjusted (AC)' & 
                                metrics$Classifier=='Thresh'& metrics$name=='Lb - HW'),]
  }
  
  if(rt==3){
    metrics <-metrics[-which(metrics$name=='Lb - Hw' & 
                         metrics$Classifier=='Thresh' & metrics$Prev.est=='Adjusted (AC)'),]
  }
  
  p <- ggplot(metrics, aes(x=size, y=value, colour=Classifier, lty=Prev.est)) +
    geom_smooth(alpha=.2, size=1,method='gam',formula=y ~ s(x, k=k,bs = "cs"),aes(fill=Classifier)) +
    ggtitle(paste(metrics.name,'-',chemname))+
    scale_x_continuous(breaks=c(0,60,120,180,240))+
    scale_fill_manual(values=c(DA="#D55E00",GMM="#009E73",
                               RF="#56B4E9",Thresh="#CC79A7"))+
    scale_color_manual(values=c(DA="#D55E00",GMM="#009E73",
                                RF="#56B4E9",Thresh="#CC79A7"))+
    theme(legend.text=element_text(size=lsize),legend.title=element_text(size=lsize,face='bold'))+
    xlab("Dataset size") + ylab(metrics.name)
  
  plot(p)
}



#######################################

#### adj.vs.unadj ####
#### some more plots to compare Adjusted Counting and Classify and Count ####

adj.vs.unadj <- function(var,metrics1,metrics2,chem,metrics.name)
  
{
  
  table1 <- metrics.table.var(var,metrics1,chem)
  table2 <- metrics.table.var(var,metrics2,chem)
  
  #Boxplot over all datasets
  
  par(mfrow=c(1,1))
  
  boxplot(table1$value,table2$value,names=c('Unadjusted (CC)','Adjusted (AC)'),outline=F,
          main=metrics.name)

  #Plot over all datasets
  
  inf.adj <- sum(table1$value>table2$value,na.rm=T)
  sup.adj <- sum(table1$value<table2$value,na.rm=T)
  
  col.points <- rep(NA,length(table2$value))
  col.points[which(table1$value>table2$value)] <- "red"
  col.points[which(table1$value<table2$value)] <- "darkgreen"
  
  
  plot(table1$value,table2$value,xlab='Unadjusted (CC)',ylab='Adjusted (AC)',
       main=paste(metrics.name,"-",var),col=col.points,xlim=c(0,1),ylim=c(0,1))
  abline(0,1,lwd=2)
  
  lims <- par("usr")
  
  text(lims[4]*0.1,lims[2]*0.9,paste('Adjusted (AC) > Unadjusted (CC) \n',sup.adj,'times'),col="darkgreen",font=2,cex=0.8)
  text(lims[4]*0.9,lims[2]*0.1,paste('Adjusted (AC) < Unadjusted (CC) \n',inf.adj,'times'),col="red",font=2,cex=0.8)
  
  
  #Boxplot over sizes
  
  par(mfrow=c(1,2))
  boxplot(table1$value~table1$size,main="Unadjusted (CC)")
  boxplot(table2$value~table2$size,main="Adjusted (AC)")
  
  #Plot over all datasets
  
  par(mfrow=c(3,3))
  
  for (i in unique(table1$size))
  {
    subTable1 <- table1[table1$size==i,]
    subTable2 <- table2[table2$size==i,]
    
    inf.adj <- sum(subTable1$value>subTable2$value,na.rm=T)
    sup.adj <- sum(subTable1$value<subTable2$value,na.rm=T)
    
    col.points <- rep(NA,length(subTable2$value))
    col.points[which(subTable1$value>subTable2$value)] <- "red"
    col.points[which(subTable1$value<subTable2$value)] <- "darkgreen"
    
    plot(subTable1$value,subTable2$value,xlab='Unadjusted (CC)',ylab='Adjusted (AC)',
         main=paste0(metrics.name,", dataset size = ",i),col=col.points,xlim=c(0,1),ylim=c(0,1))
    
    lims <- par("usr")
    legend(lims[4]*0.1,lims[2]*0.85,sup.adj,box.col="darkgreen",cex=1.1)
    legend(lims[4]*0.7,lims[2]*0.2,inf.adj,box.col="red",cex=1.1)
    
    
    abline(0,1)
    
    
  }
  
}
 
#######################################

#### var.plot, var.plot2, var.plot3 ####
#### functions to display performance variance across methods and datasets 
#### for each metrics

### var.plot: all variance treatments are pooled

var.plot <- function(chem,metrics,metrics.name,chem.name){
  
  
  allvalues <- do.call(rbind,
                             lapply(names(chem), function(X){
                               metrics.table.var(X,metrics,chem) }))
  
  
  allmeans <- with(allvalues,tapply(value,list(variable,size),
                                             function(x){mean(x,na.rm=T)}))
  
  allsd <- with(allvalues,tapply(value,list(variable,size),
                                          function(x){sd(x,na.rm=T)}))
  
  Std <- melt(allsd)

  ggplot(Std,aes(x=Var2,y=value,colour=Var1))+
    geom_line(lwd=2)+
    ggtitle(paste(metrics.name,'-',chem.name))
  
}

### var.plot2: variance treatment per variance treatment

var.plot2 <- function(chem,metrics,metrics.name,chem.name){
  
  
  allvalues <- do.call(rbind,
                       lapply(names(chem), function(X){
                         metrics.table.var(X,metrics,chem) }))
  
  
  allmeans <- with(allvalues,tapply(value,list(variable,size),
                                    function(x){mean(x,na.rm=T)}))
  
  allsd <- with(allvalues,tapply(value,list(variable,size,name),
                                 function(x){sd(x,na.rm=T)}))
  
  Std <- na.omit(melt(allsd))
  
  ggplot(Std,aes(x=Var2,y=value,colour=Var1))+
    geom_line(lwd=2)+
    ggtitle(paste(metrics.name,'-',chem.name))+facet_wrap(~Var3)+
    ylab(paste(metrics.name,'Std'))
}

### var.plot2: variance treatment per variance treatment with both AC and CC

var.plot3 <- function(chem,metrics1,metrics2,metrics.name,chem.name){
  
  
  allvalues1 <- do.call(rbind,
                       lapply(names(chem), function(X){
                         metrics.table.var(X,metrics1,chem) }))

  allvalues2 <- do.call(rbind,
                        lapply(names(chem), function(X){
                          metrics.table.var(X,metrics2,chem) }))
  
  
  allsd1 <- with(allvalues1,tapply(value,list(variable,size,name),
                                 function(x){sd(x,na.rm=T)}))
  
  allsd2 <- with(allvalues2,tapply(value,list(variable,size,name),
                                   function(x){sd(x,na.rm=T)}))
  
  Std1 <- na.omit(melt(allsd1))
  Std1$Prev.est <- 'Unadjusted'
  
  Std2 <- na.omit(melt(allsd2))
  Std2$Prev.est <- 'Adjusted'
  
  Std <- rbind(Std1,Std2)
  Std$Prev.est <- relevel(as.factor(Std$Prev.est),'Unadjusted')
  
  Std[Std$value==0,'value'] <- min(Std[Std$value!=0,'value'])

    ggplot(Std,aes(x=Var2,y=value,colour=Var1,linetype=Prev.est))+
  geom_line(lwd=1)+
  ggtitle(paste(metrics.name,'-',chem.name))+facet_wrap(~Var3)+
      ylab(paste(metrics.name,'Std'))
}
