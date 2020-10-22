## Scripts to estimate model parameters on the Aphidius ervi dataset

### See Appendix II for model description

library(data.table)
set.seed(42)

# Set working directory to load the A. ervi dataset (aervi_all.txt)
# Loading the data
aervi.table <- read.table("aervi_data.txt")
aervi.table$TimeF <- as.factor(aervi.table$Time)
categ <- c("Unfed","EFN","ApH")


  ####---- MAIN SUGARS ----####

# same script for estimating sugar amounts and ratios
stock <- list()
ratios <- list()

for (c in categ)
{
  SubTest <- aervi.table[aervi.table$Sugar_treatment == c,]
  
  # let's say that : Glu+Fru+Su ~ G+F+S(t=0) * exp(r*t)
  modL <- lm(log(Glucose + Fructose + Sucrose) ~ Time,data=SubTest)
  
  # i simply print but we could define some condition(s),  
  # on the summary(modL)$adj.r.squared for example
  print(summary(modL))
  
  # keeping the results 
  stock[[c]] <- c(modL$coefficients,summary(modL)$sigma)
  names(stock[[c]]) <- c("alpha","beta","rse")
  
  # same here with this print (is time a relevant factor on ratios of interest)
  print(anova(lm(Glucose ~ Fructose,data=SubTest),lm(Glucose ~ Fructose:TimeF,data=SubTest)))
  
  ratiogf <- lm(Glucose ~ Fructose:TimeF-1, data = SubTest)
  
    ## Means
    ratios[[c]][["GF_ratio(t)"]] <- summary(ratiogf)$coefficients
    
    ## Sd: we keep the possibility to change that later but for now, let's consider
    ## a simple homoschedastic model
    ratios[[c]][["GF_ratio(t)"]][,2] <- summary(ratiogf)$sigma
  
  # now for sucrose 
  print(anova(lm(Fructose ~ Sucrose,data=SubTest),lm(Fructose ~ Sucrose:TimeF,data=SubTest)))
  
  ratiofs <- lm(Fructose ~ Sucrose-1, data = SubTest)
  ratios[[c]][["FS_ratio(t)"]] <- summary(ratiofs)$coefficients
}


  ####---- WEAK SUGARS ----####

 ##### Maltose ####
# Time dependence for nb of 0 --> time-dependent binomial
# When sugar amount > 0: marginally significant differences t0 and t>0 but not considered
# & No difference between classes
# --> 1 distrib for amount >0

## A. Nb of 0

zeroMalto <- list() # to store the models
zeroMaltoProba <- list() # to store the proba for each time and treatment

## Note : here I used time as a factor (TimeF), so we work only with 5 time points
## We could also work with a continuous time (Time)

for (c in categ)
{
  #Maltose
  SubTest <- aervi.table[aervi.table$Sugar_treatment == c,]
  SubTest$BM <- SubTest$Maltose == 0
  mod <- glm(BM ~ TimeF-1, family=binomial,data=SubTest)
  zeroMalto[[c]] <- mod
  
  #Coefs
  # If I'm not mistaken, we have P=odds/(1+odds) with odds=exp(logit)
  # Not done here, but we could also estimate the std error of the probas ?
  odds <- exp(zeroMalto[[c]]$coefficients)
  zeroMaltoProba[[c]] <- odds/(1+odds)
  
  # prediction on training set
  binP=predict(zeroMalto[[c]],type = 'response') > 0.5 
  print(paste("Maltose, nb of 0,",c, round(sum(binP == SubTest$BM)/length(SubTest$BM),2) ))
}

## B. Trying to fit non-0 values
# Only one distrib to rule them all ?

#Some t-tests
for (c in categ)
{
  SubTest <- aervi.table[aervi.table$Sugar_treatment == c & aervi.table$Maltose>0,]
  Mt0 <- SubTest[SubTest$Time==0,"Maltose"]
  Mt1 <- SubTest[SubTest$Time!=0,"Maltose"]
  print(c)
  print(t.test(Mt0,Mt1))
}
# Some marginally significant differences for Unfed & ApH
# But let's try without time-dependence

posMalt <-  aervi.table[aervi.table$Maltose!=0,]
GausMalto <- glm(posMalt$Maltose ~ 1,family="gaussian")
gammaMalto <- glm(posMalt$Maltose ~ 1,family="Gamma")

AIC(GausMalto)
AIC(gammaMalto)
## deltaAIC =  - 78
## and gamma predicts only positive values

plot(gammaMalto) #fit is not perfect but way better than lm
print(gammaMalto)

#gamma parameters : shape (alpha) and rate (beta)
#rate = alpha / µ

shapeMalt <- 1/(summary(gammaMalto)$disp)
meanMalt <- 1/gammaMalto$coefficients
rateMalt <- shapeMalt/meanMalt

#predictions
sim <- rgamma(n=nrow(posMalt),shape = shapeMalt, rate=rateMalt)
plot(sim,posMalt$Maltose) 
abline(0,1) #quite symmetrical


#### Stacchyose ####
# No time dependence for nb of 0 --> Binomial with cst proba
# When sugar amount > 0: significant difference t0 and t>0 (all treatments)
# Difference between classes
# --> 6 distrib for amount >0 (2*each categ)

##A. Nb of 0

zeroStac <- list()
zeroStacProba <- list() 

for (c in categ)
{
  SubTest <- aervi.table[aervi.table$Sugar_treatment == c,]
  SubTest$BS <- SubTest$Stacchyose == 0
  
  #mod <- glm(BS ~ Time, family=binomial,data=SubTest) <- Time not signif.
  mod <- glm(BS ~ 1, family=binomial,data=SubTest)
  zeroStac[[c]] <- mod
  
  #probas
  odds <- exp(mod$coefficients)
  proba <- odds/(1+odds)
  zeroStacProba[[c]] <- proba
}
### here the proba is just the % of 0


##B. Positive values

GammaStac <- list() 

for (c in categ)
{
  SubTest <- aervi.table[aervi.table$Sugar_treatment == c & aervi.table$Stacchyose!=0,]
  Sub0 <- SubTest[SubTest$Time==0,]
  Sub1 <- SubTest[SubTest$Time!=0,]
  
  #t-tests -> difference for all groups
  print(c)
  print(t.test(Sub0$Stacchyose,Sub1$Stacchyose))
  
  #Time 0
  #Results not shown here : always better with Gamma than with Gaussian models (for t0 and t>0)
  modGamma0 <- glm(Stacchyose ~ 1, family=Gamma,data=Sub0)
  GammaStac[[c]][["T0"]] <- modGamma0
  #plot(modGamma0)
  
  #Time > 0
  modGamma1 <- glm(Stacchyose ~ 1, family=Gamma,data=Sub1)
  #plot(modGamma1)
  GammaStac[[c]][["T1"]] <- modGamma1
  
  #Predictions
  # Time 0
  shape0 <- 1/(summary(modGamma0)$disp)
  rate0 <- shape0 * modGamma0$coefficients
  
  plot(rgamma(nrow(Sub0),shape0,rate0),Sub0$Stacchyose, main=paste(c,"Stacchyose, Time = 0"),xlab="Simul",ylab="Real")
  abline(0,1)
  
  # Time 1
  shape1 <- 1/(summary(modGamma1)$disp)
  rate1 <- shape1 * modGamma1$coefficients
  
  plot(rgamma(nrow(Sub1),shape1,rate1),Sub1$Stacchyose, main=paste(c,"Stacchyose, Time > 0"),xlab="Simul",ylab="Real")
  abline(0,1)
  
  #Time0
  
  # Time 1
  
}

# Does not seem to work so bad ? Simulated values have +/- the same distrib than real values
# When values are very small, predictions are not so good, but that's no big deal for this sugar


#### Melezitose ####
# No time dependence for nb of 0 --> 1 bino with cst proba / class
# When sugar amount > 0: difference for all groups (but strong only for ApH)
# Difference between classes
# Here I took 6 distribs (3*2) for sugar amount >0
# We could also simplify with 4 distribs (EFN+Unfed+2*ApH)
# Additional note : it could be interesting to have a time-dep nb of 0 for ApH (Wäckers et al. 2006, Hogervorst et al. 2007)

zeroMele <- list()
zeroMeleProba <- list()

for (c in categ)
{
  SubTest <- aervi.table[aervi.table$Sugar_treatment == c,]
  SubTest$BMe <- SubTest$Melezitose == 0
  
  #mod1 <- glm(BMe ~ Time , family=binomial,data=SubTest) <- Time not significant
  mod <- glm(BMe ~ 1, family=binomial,data=SubTest)
  zeroMele[[c]] <- mod
  
  #probas
  odds <- exp(mod$coefficients)
  proba <- odds/(1+odds)
  zeroMeleProba[[c]] <- proba
}  
### here the proba is just the % of 0

GammaMele <- list()

for (c in categ)
{
  SubTest <- aervi.table[aervi.table$Sugar_treatment == c & aervi.table$Melezitose!=0,]
  Sub0 <- SubTest[SubTest$Time==0,]
  Sub1 <- SubTest[SubTest$Time!=0,]
  
  #t-tests -> difference for all groups
  print(c)
  print(t.test(Sub0$Melezitose,Sub1$Melezitose))
  
  #Time 0
  #Results not shown here : always better with Gamma than with Gaussian models (for t0 and t>0)
  modGamma0 <- glm(Melezitose ~ 1, family=Gamma,data=Sub0)
  GammaMele[[c]][["T0"]] <- modGamma0
  #plot(modGamma)
  
  #Time > 0
  modGamma1 <- glm(Melezitose ~ 1, family=Gamma,data=Sub1)
  #plot(modGamma1)
  GammaMele[[c]][["T1"]] <- modGamma1
  
  #Predictions
  # Time 0
  shape0 <- 1/(summary(modGamma0)$disp)
  rate0 <- shape0 * modGamma0$coefficients
  
  plot(rgamma(nrow(Sub0),shape0,rate0),Sub0$Melezitose, main=paste(c,"Melezitose, Time = 0"),xlab="Simul",ylab="Real")
  abline(0,1)
  
  # Time 1
  shape1 <- 1/(summary(modGamma1)$disp)
  rate1 <- shape1 * modGamma1$coefficients
  
  plot(rgamma(nrow(Sub1),shape1,rate1),Sub1$Melezitose, main=paste(c,"Melezitose, Time > 0"),xlab="Simul",ylab="Real")
  abline(0,1)
  
}


#### Erlose ####
# No time dependence for nb of 0
# When sugar amount > 0: some significant differences t0 and t>0 --> see below
# Difference between categories
# Binomial with cst proba (0/not 0)
# 6 distrib for amount >0 (2*each categ) ? --> see below

## A. Nb of 0

zeroErloProba <- list() 

for (c in categ)
{
  SubTest <- aervi.table[aervi.table$Sugar_treatment == c,]
  SubTest$BE <- SubTest$Erlose == 0
  modG <- glm(BE ~ TimeF - 1, family=binomial,data=SubTest)
  # Note : time signif. for the "ApH" treatment, but not considered here
  mod <- glm(BE ~ 1, family=binomial,data=SubTest)
  
  #probas
  odds <- exp(mod$coefficients)
  proba <- odds/(1+odds)
  zeroErloProba[[c]] <- proba
}

# Note : for Unfed & EFN individuals we always have erlose = 0, so we could consider P(erlose = 0) = 1 if T > 0
# but some species can synthetise erlose from sucrose during their lifetime (cf Wäckers et al. 2006, Hogervorst et al. 2007)
# which can be a factor of confusion for predicitons -> interesting
# so let us keep a (low) probability to have a small value for erlose for these treatments, event at T > 0 (so no time dependence)


## B. Trying to fit non-0 values

aervi.table[aervi.table$Erlose!=F,c("Tmt2","Erlose")]
# Very few values >0 for Unfed and EFN
# Suggestion : instead of 3 classes, let's model two classes "ApH" and "not ApH" ?
# (Makes sense : honeydew consumption should change erlose values but not nectar consumption)
# "not ApH" : no time-dependence (we have values only for T=0 anyway)
# "ApH" : no time-dependence (cf t-tests below)
# However, ApH : no time-dependence at all does not seem super realistic (cf Wäckers et al. 2006) <- maybe add later
# Besides here we have only 0 values for this treatment if T > 12
# here if we consider the same distrib for all times for ApH, we'll have the same kind of values at T=24 and T=48
# solution : this distrib until T=12, and for T>12 : take the same distrib than Unfed and EFN ?


## 'Not ApH'
aervi_noH <- aervi.table[aervi.table$Sugar_treatment!="ApH" & aervi.table$Erlose !=0,]

gammaE.nH <- glm(aervi_noH$Erlose~1,family="Gamma")
shape <-  1/(summary(gammaE.nH)$disp)
rate <- shape * gammaE.nH$coefficient

# Some simulations
rgamma(10,shape,rate) # Values are okay

## 'ApH'
SubTest <- aervi.table[aervi.table$Sugar_treatment == "ApH" & aervi.table$Erlose!=0,]
Sub0 <- SubTest[SubTest$Time==0,]
Sub1 <- SubTest[SubTest$Time!=0,]
t.test(Sub0$Erlose,Sub1$Erlose) #Not signif

gaussE <- glm(Erlose ~ 1, family="gaussian",data=SubTest)
gammaE.H <- glm(Erlose ~ 1, family="Gamma",data=SubTest)
print(AIC(gaussE)-AIC(gammaE.H))
#Not such a gain but let's try gamma

print(gammaE.H)
shape <-  1/(summary(gammaE.H)$disp)
rate <- shape * gammaE.H$coefficient

sim <- rgamma(nrow(SubTest),shape,rate)
plot(sim,SubTest$Erlose)
abline(0,1)  # Works okay !!

