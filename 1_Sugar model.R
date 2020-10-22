#### SIMULATION MODEL ####

### See Appendix II for model description

### PARAMETERS
####
#categ = name of categories
  # e.g. categ <- c("Unfed","EFN","ApH")

#n.categ = n per treatment
  # e.g. n.categ <- c(60,60,60)

# Unf = Which treatment corresponds to the "Unfed" treatment?
 # if "Unfed" is in first position in vector "categ" :  Unf <- 1 

#n.Unfed : How many are "Unfed", how many are "starved"?
#(Order : unfed, starved EFN, starved ApH)
  # e.g. n.Unfed <- c(20,20,20)

#lim = detection limit: after what time individuals will be considered as "unfed ?"
  # e.g. lim <- 12 for a 12 'detection time'

#"stock" parameters (glucose + fructose + sucrose)
#list of 3 lists: each has a treatment name and 3 values: alpha, beta and rse for each treatment

#ratio parameters (glucose/fructose, fructose/sucrose)
#list of 3 lists (treatments), each has a treatment name
#and two lists for ratios, each containing a 2-column object: estimate and sd for ratios for each time point and treatment

#Times = time points
  # e.g. Times <- c(0,1,12,24,48)

#nb.simul = number of simulations
  # e.g. nb.simul <- 100

#minGF = Minimum Glucose/Fructose ratio and maxTS = max Total Sugar Amount
#- in order not to be too much unrealistic
  #e.g. minGF <- 1.5
  #e.g. maxTS <- 400

####


#Note: parameters for the 'weak sugars' (maltose, stacchyose, erlose, melezitose) are not alterable here
#such as zeroErloProba, GammaStac, etc...
#Because we'll always use the ones estimated on the A. ervi dataset
#see "Estimating A. ervi parameters" script

sugar.simul <- function(categ,n.categ,Unf,n.Unfed,lim,stock,ratios,Times,nb.simul,minGF,maxTS)
  
{

  #Storing datasets here
  data.list <- list()

  #Preparing dataset format
  untLim <- Times[Times<=lim]
  nrep.f <- n.categ[-Unf]%/%length(untLim)
  
  Treatments.fed <- rep(categ[-Unf],nrep.f*length(untLim))
  
  time_list.fed <- lapply(nrep.f,function(x){sort(rep(untLim ,x))})
  Times.fed <- unlist(time_list.fed)
  
  #### Unfed
  
  if(all.equal(n.categ[Unf],sum(n.Unfed))==F)
  {
    stop("Nb of individuals in n.Unfed must equal number of Unfed in categ")
  }
  
  beyLim <- Times[Times>lim]
  nrep.unf <- n.Unfed[Unf] %/% length(Times)
  nrep.st <- n.Unfed[-Unf] %/% length(beyLim)
  
  nrep.unfTot <- c(nrep.unf*length(Times) , nrep.st*length(beyLim))
  
  Treatments.unfed.pre <- rep( c(categ[Unf],categ[-Unf]),nrep.unfTot)
  Treatments.unfed.post <- rep(categ[Unf],sum(nrep.unfTot))
  
  
  Times.unf <- rep(Times,nrep.unf)
  time_list.st <- lapply(nrep.st,function(x){sort(rep(beyLim ,x))})
  Times.st <- unlist(time_list.st)
  
  
  Treatments <- as.factor(c(Treatments.fed,Treatments.unfed.post))
  Tmt.real <-  as.factor(c(Treatments.fed,Treatments.unfed.pre))
  
  Time <- c(Times.fed,Times.unf,Times.st)


for (k in 1:nb.simul)
{
  
  simul.data <- data.frame(Tmt.f=Treatments,Tmt=Tmt.real,Time=Time,Fructose=0,Glucose=0,Sucrose=0,Stacchyose=0,Erlose=0,Melezitose=0,Maltose=0)

  ### A. Glucose, Fructose and Sucrose values
  
  for (c in categ)
  {
    
    # We'll store data for this category here
    timee <- simul.data[simul.data$Tmt==c,"Time"]
    tmt.data <- data.frame(Time=timee,Glucose=0,Fructose=0,Sucrose=0)

    # Stock parameters for this category
    alpha_stock <- stock[[c]]["alpha"]
    beta_stock <- stock[[c]]["beta"]
    rse_stock <- stock[[c]]["rse"]
    
    # A linear model is used to estimate the log "stock" amount (Glu+Fru+Sucr)
    ymean <- alpha_stock + beta_stock*timee
    ypred <- ymean + rnorm(length(ymean),mean = 0, sd = rse_stock)
    GFSt <- exp(ypred)
    
    # Okay, let's just setup a bound here so that we don't get unrealistic profiles
    if(any(GFSt > maxTS)){GFSt[GFSt > maxTS] <- maxTS - runif(length(GFSt[GFSt > maxTS]),0,maxTS/2)}

    # Now let's get the ratio parameters for this category
    coGF <- ratios[[c]][["GF_ratio(t)"]]
    coFS <- ratios[[c]][["FS_ratio(t)"]]
    
    # let's simulate data with some variation for the individual ratios
    
    for (j in 1:length(unique(timee)))
    {
      u <- unique(timee)[j]
      
      coFS_withvar <- rep(coFS[1],length(GFSt[timee == u])) + rnorm(length(GFSt[timee == u]),mean=0,sd=coFS[2])
      coGF_withvar <- rep(coGF[j,1],length(GFSt[timee == u])) + rnorm(length(GFSt[timee == u]),mean=0,sd=coGF[j,2]) 
      
      # Attempt to avoid negative ratios : If ratio < 0 -> ratio = min(ratio > 0)
      # If there are only negative ratios in this time point, we seek in other time points
      if(sum(coFS_withvar<0)>0)
      {
        if(sum(coFS_withvar<0)<length(coFS_withvar))
          {
            coFS_withvar[coFS_withvar<0] <- min(coFS_withvar[coFS_withvar>0])
        }else {
        coFS_withvar <- min(tmt.data$Fructose/tmt.data$Sucrose)
      }}
        
      if(any(coGF_withvar<minGF))
      {
        while(all(coGF_withvar<minGF))
        {
          coGF_withvar <- rep(coGF[j,1],length(GFSt[timee == u])) + rnorm(length(GFSt[timee == u]),mean=0,sd=coGF[j,2]) 
        }
        coGF_withvar[coGF_withvar<minGF] <- min(coGF_withvar[coGF_withvar>=minGF])
      }
      
      Tot <- (coGF_withvar * coFS_withvar) + coFS_withvar + 1
      
      Sucrose <- (GFSt[timee == u])/Tot
      Fructose <- coFS_withvar*Sucrose
      Glucose <- coGF_withvar*Fructose
      
      tmt.data$Sucrose[timee==u] <- Sucrose
      tmt.data$Fructose[timee==u] <- Fructose
      tmt.data$Glucose[timee==u] <- Glucose
      
    }
    
    # Just checking if sugar calculations are consistent
    test.equal <- all.equal( (tmt.data[,1]+tmt.data[,2]+tmt.data[,3]), GFSt) 
    if(isTRUE(test.equal==F))
    {
      cat("Sums of sugars are not equal to stock for dataset",k,"and treatment",c,"")
    }
      
    
    # Here we put the data in the real dataset
    simul.data[simul.data$Tmt==c,'Sucrose'] <-  tmt.data$Sucrose
    simul.data[simul.data$Tmt==c,'Fructose'] <- tmt.data$Fructose
    simul.data[simul.data$Tmt==c,'Glucose'] <- tmt.data$Glucose
  }

  
  #### 2. 'Weak' sugars
  
  ### MALTOSE
  # For each treatment and each time --> random sample of 0 and 1, taking proba from binomial models as the proba to be "0"
  for (c in categ)
  {
    for (i in 1:length(Times))
    {
      time <- unique(simul.data$Time)[i] 
      proba <- zeroMaltoProba[[c]][i]
      nb.sampl <- nrow(simul.data[simul.data$Time==time & simul.data$Tmt==c,])
      sampl <- sample(c(0,1),  nb.sampl , replace=T, prob=c(proba, 1-proba))
      simul.data$Maltose[simul.data$Time==time & simul.data$Tmt==c] <- sampl
    }
  }

  
  # Now let's set the positive values
  posMalt <- which(simul.data$Maltose!=0)
  simul.data$Maltose[posMalt] <- rgamma(length(posMalt),shapeMalt,rateMalt)
  
  ### STACCHYOSE
  # For each treatment, one proba to be 0
  #For the positive values, Two distribs for each treatment : T0 and T >0
  
  for (c in categ)
  {
    proba <- zeroStacProba[[c]]
    nb.sampl <- length(simul.data[simul.data$Tmt==c,"Stacchyose"] )
    simul.data$Stacchyose[simul.data$Tmt==c] <- sample(c(0,1), nb.sampl, replace=T, prob=c(proba,1-proba))
 
    SubTest <- simul.data[simul.data$Tmt == c & simul.data$Stacchyose!=0,]
    Sub0 <- SubTest[SubTest$Time==0,"Stacchyose"]
    Sub1 <- SubTest[SubTest$Time!=0,"Stacchyose"]
    
    #T = 0
    shape0 <- 1/summary(GammaStac[[c]][["T0"]])$disp
    rate0 <- shape0 * GammaStac[[c]][["T0"]]$coefficients
    sim0 <- rgamma(length(Sub0),shape0,rate0)
    
    #T > 0
    shape1 <- 1/summary(GammaStac[[c]][["T1"]])$disp
    rate1 <- shape1 * GammaStac[[c]][["T1"]]$coefficients
    sim1 <- rgamma(length(Sub1),shape1,rate1)
    
    simul.data[simul.data$Tmt == c & simul.data$Stacchyose!=0 & simul.data$Time==0,"Stacchyose"] <- sim0
    simul.data[simul.data$Tmt == c & simul.data$Stacchyose!=0 & simul.data$Time!=0,"Stacchyose"] <- sim1
    
  }
  
  #MELEZITOSE
  # Same idea than for Stacchyose
  
  for (c in categ)
  {
    proba <- zeroMeleProba[[c]]
    nb.sampl <- length(simul.data[simul.data$Tmt==c,"Melezitose"] )
    simul.data[simul.data$Tmt==c,"Melezitose"] <- sample(c(0,1), nb.sampl, replace=T, prob=c(proba,1-proba))

    SubTest <- simul.data[simul.data$Tmt == c & simul.data$Melezitose!=0,]
    Sub0 <- SubTest[SubTest$Time==0,"Melezitose"]
    Sub1 <- SubTest[SubTest$Time!=0,"Melezitose"]
    
    #T = 0
    shape0 <- 1/summary(GammaMele[[c]][["T0"]])$disp
    rate0 <- shape0 * GammaMele[[c]][["T0"]]$coefficients
    sim0 <- rgamma(length(Sub0),shape0,rate0)
    
    #T > 0
    shape1 <- 1/summary(GammaMele[[c]][["T1"]])$disp
    rate1 <- shape1 * GammaMele[[c]][["T1"]]$coefficients
    sim1 <- rgamma(length(Sub1),shape1,rate1)
    
    simul.data[simul.data$Tmt == c & simul.data$Melezitose!=0 & simul.data$Time==0,"Melezitose"] <- sim0
    simul.data[simul.data$Tmt == c & simul.data$Melezitose!=0 & simul.data$Time!=0,"Melezitose"] <- sim1
    
  }
  
  #ERLOSE

  #let's do this for now:
  # 0/not 0 --> 1 bino/treatment
  #ApH, T <= 12 --> 1 distrib (uniform)
  #Other treatments --> 1 distrib (gamma)
  
  #First, 0/not 0
  for (c in categ)
  {
    proba <- zeroErloProba[[c]]
    nb.sampl <- length(simul.data[simul.data$Tmt==c,"Erlose"] )
    simul.data[simul.data$Tmt==c,"Erlose"] <- sample(c(0,1), nb.sampl, replace=T, prob=c(proba,1-proba))
  }
  
  #Not 0 for ApH and T <= 12: gamma distrib

  shapeErl.H <-  1/(summary(gammaE.H)$disp)
  rateErl.H <- shape * gammaE.H$coefficient
  
  posErl.H <- simul.data[simul.data$Tmt=="ApH" & simul.data$Time <=12 & simul.data$Erlose!=0,"Erlose"]
  SubTest1 <- rgamma(length(posErl.H),shapeErl.H,rateErl.H)

  #Not 0 for the others: gamma distrib too
  
  
  shapeErl.nH <- 1/(summary(gammaE.nH))$disp
  rateErl.nH <- shapeErl.nH * gammaE.nH$coefficients
  
  posErl.nH <- simul.data[(simul.data$Tmt!="ApH" | simul.data$Time > 12) & simul.data$Erlose!=0,"Erlose"]
  SubTest2 <- rgamma(length(posErl.nH),shapeErl.nH,rateErl.nH)
  
  
  simul.data[simul.data$Tmt=="ApH" & simul.data$Time <=12 & simul.data$Erlose!=0,"Erlose"] <- SubTest1
  simul.data[(simul.data$Tmt!="ApH" | simul.data$Time > 12) & simul.data$Erlose!=0,"Erlose"] <- SubTest2 
  
  
 
  # Now let's add the GF Ratio, The H Ratio and the Residual Sugar amount
  simul.data$GF_Ratio <- simul.data$Glucose/(simul.data$Glucose+simul.data$Fructose)
  simul.data$TotSug <- rowSums(simul.data[,4:10])
  simul.data$H_Ratio <- rowSums(simul.data[,8:10]) / simul.data$TotSug
  simul.data$resSug <- simul.data$TotSug-simul.data$Fructose
  simul.data$RS_Ratio <- simul.data$resSug/simul.data$TotSug
  simul.data$F_Ratio <- simul.data$Fructose/simul.data$TotSug
  
  data.list[[k]] <- na.omit(simul.data)
}

return(data.list)

}