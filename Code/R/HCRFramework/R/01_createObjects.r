#-------------------------------------------------------------------------------
#
# Script creates true and observed stock using the output of the assessment
#
# Created by: Thomas Brunel, Niels Hintzen. niels.hintzen@wur.nl
# Project: South Pacific
#
# Rev:
#  1- stock object has 4 areas now, for the 4 fleets
#  2- create the historical stocks by  resampling initial values using the
#     variance covariance matrix
#
#
#  developed under  R version 3.0.1
#            using  FLCore    2.5.0
#
# DO NOT DISTRIBUTE OR COPY THIS CODE WITHOUT PRIOR AGREEMENT OF THE AUTHORS
#
#-------------------------------------------------------------------------------

rm(list=ls())

library(FLCore)
library(FLFleet)
library(PBSadmb)
library(lattice)
library(MASS)

#- Set paths
codePath        <- "D:/Repository/JackMackerel/Code/R/HCRFramework/R/"
dataPath        <- "D:/Repository/JackMackerel/Code/admb/arc/"
outPath         <- "N:/Projecten/SouthPacific/WP4/R/code/submit2GitHub/Results/"

Fname           <- "mod3.1_R.rep"
res             <- readList(file.path(dataPath,Fname))

source(file.path(codePath,"functions.r"))

#-------------------------------------------------------------------------------
# Set key dates
#-------------------------------------------------------------------------------

histMinYr       <- res$N[1,1]
histMaxYr       <- res$N[dim(res$N)[1],1]
nyrs            <- 30 # !!! max 43 years if selected the SGRec regime....
futureMaxYr     <- histMaxYr + nyrs
histPeriod      <- ac(histMinYr:histMaxYr)
projPeriod      <- ac((histMaxYr+1):futureMaxYr)
RecRegime       <- "SGRec" #  either "LTRec", "STRec"  Long term vs short term time series used for rct model, or SGrec : surrowgates time series
if(RecRegime=="SGRec"){
  recrPeriod  <- ac(2002:histMaxYr)
} else {recrPeriod  <- histPeriod}
if(RecRegime=="STRec"){
recrPeriod  <- ac(2002:histMaxYr)
} else {recrPeriod  <- histPeriod}  # chose either 1970:2012 or 2000:2012
selPeriod       <- ac(2002:histMaxYr) #last sel change Fish1=2003, Fish2=2005, Fish3=2002, Fish4=2005/2011
fecYears        <- ac(histMaxYr) #no change in maturity
nits            <- 100
save(list=c("histMinYr","histMaxYr","futureMaxYr","histPeriod","projPeriod","nyrs","nits"),
     file=paste(outPath,"yearParams.RData",sep=""))

settings        <- list(histMinYr=histMinYr,histMaxYr=histMaxYr,futureMaxYr=futureMaxYr,
                        histPeriod=histPeriod,projPeriod=projPeriod,recrPeriod=recrPeriod,
                        nyrs=nyrs,nits=nits,fecYears=fecYears,RecRegime=RecRegime)
datum           <- date()
years           <- res$Yr
ages            <- 1:12
fisheries       <- 1:length(res$Fshry_names)
save(years,ages,fisheries,file=paste(outPath,"dimensions.RData",sep=""))

#-------------------------------------------------------------------------------
# 1): Create the biological object (the assumed true dynamics)
#-------------------------------------------------------------------------------

JMB                       <- FLStock(stock.n=FLQuant(dimnames=
                                      list(age=ages,year=years,unit="unique",season="all",
                                           area=c(res$Fshry_names),iter=1:nits)))
JMB@stock.n[,histPeriod]  <- t(res$N[,-1])
JMB                       <- window(JMB,end=futureMaxYr)
JMB@m.spwn[]              <- 10.5/12
JMB@harvest.spwn[]        <- 10.5/12
JMB@m[]                   <- res$M[1]
JMB@mat[]                 <- res$mature_a
JMB@stock.wt[]            <- res$wt_a_pop
desc(JMB)                 <- "Jack Mackerel Biological Unit"
units(JMB)[1:17]          <- as.list(c(rep(c("tonnes","thousands","kg"),4),
                                         rep("NA",2),"f",rep("NA",2)))

#-------------------------------------------------------------------------------
# 2): Create fisheries object
#-------------------------------------------------------------------------------

dmns                      <- dimnames(m(JMB))
fishery                   <- FLCatch(price=FLQuant(NA,dimnames=dmns))
name(fishery)             <- "Catches"
desc(fishery)             <- "South Pacific Jack Mackerel"
fishery@range             <- range(JMB)

  #-----------------------------------------------------------------------------
  # Simulate weights at age for the two groups (ie fisheries 124
  # and fishery3). (there are only two wt-at-age groups for 4 fisheries)
  #-----------------------------------------------------------------------------

#- for fhs 1
fsh   <- 1 
source(paste(codePath,"01a_generateWeightAtAge.r",sep=""))
for(i in 1:nits)
  fishery@landings.wt[,ac(histMinYr:futureMaxYr),,,c(1),i][]  <- wt[,ac(histMinYr:futureMaxYr),,,,i]
rm(wt)

#- for fhs 2,4
fsh   <- 2
source(paste(codePath,"01a_generateWeightAtAge.r",sep=""))
for(i in 1:nits)
  fishery@landings.wt[,ac(histMinYr:futureMaxYr),,,c(2,4),i][]  <- wt[,ac(histMinYr:futureMaxYr),,,,i]
rm(wt)


#- for fhs 3
fsh   <- 3
source(paste(codePath,"01a_generateWeightAtAge.r",sep=""))
for(i in 1:nits)
  fishery@landings.wt[,ac(histMinYr:futureMaxYr),,,3,i]           <- wt[,ac(histMinYr:futureMaxYr),,,,i]
rm(wt)

#- Take landings from asssesment output
for(iFishery in fisheries)
  fishery@landings.n[,histPeriod,,,iFishery,]  <- t(res[[which(names(res)==paste("C_fsh_",iFishery,sep=""))]][,-1])
landings(fishery)                     <- quantSums(fishery@landings.n * fishery@landings.wt)
catch.q(     fishery)[]               <- 1
discards.sel(fishery)[]               <- 0
fishery@discards.wt[]                 <- 0
fishery@discards.n[]                  <- 0

#- Take selection from assessment output
for(iFishery in fisheries)
  fishery@landings.sel[,histPeriod,,,iFishery,]<- t(res[[which(names(res)==paste("sel_fsh_",iFishery,sep=""))]][,-c(1:2)])

#-------------------------------------------------------------------------------
# 3): update the biological object with catch information
#-------------------------------------------------------------------------------

#- Stock weights in the assessment are the historic average of the weights in CHile
JMB@landings.wt                       <- landings.wt(fishery)
JMB@catch.wt                          <- landings.wt(fishery)
JMB@discards.wt                       <- discards.wt(fishery)

#- Input the historic fishing mortality
for(iFishery in fisheries)
  JMB@harvest[,histPeriod,,,iFishery,]<-  t(res[[which(names(res)==paste("F_age_",iFishery,sep=""))]][,-1])
  # check : res$TotF[10,] is the same as areaSums(JMB@harvest[,ac(1979)])
 
#- Compute the historical catch
JMB@landings.n                        <- fishery@landings.n
JMB@discards.n                        <- fishery@discards.n
JMB@catch.n                           <- catch.n(fishery)
JMB@catch                             <- quantSums(JMB@catch.n * JMB@catch.wt)
JMB@landings                          <- quantSums(JMB@landings.n * JMB@landings.wt)
JMB@discards                          <- quantSums(JMB@discards.n * JMB@discards.wt)
JMB@stock                             <- quantSums(JMB@stock.n * JMB@stock.wt)

  #-------------------------------------------------------------------------------
  # 4):  stock recruitment relationships 
  #-------------------------------------------------------------------------------

#- Load data or compute new data (requires WinBUGS run, so be careful)
newdata                               <- T

# Bayesian method
# run the bayesian estimation?
if (newdata==T & (RecRegime=="STRec" | RecRegime=="LTRec"))
{
  # compute the SR relationships
  # perform the bayesian estimates for each SR model
  source(paste(codePath,"01b_WinBUGSwithR.r",sep=""))
  # prepare a set of SR models based on the probability of each functional relationship
  source(paste(codePath,"01b_WinBUGSwithR_partII.r",sep=""))
}  
#import the SR relationship
if( RecRegime== "STRec" ) SR  <-  read.csv(paste(dataPath,"SRbayeswith R_1000models_short_N_Norm.csv",sep=""))[,-1]
if( RecRegime== "LTRec" ) SR  <-  read.csv(paste(dataPath,"SRbayeswith R_1000models_all_N_Norm.csv",sep=""))[,-1]

if (newdata==T & RecRegime== "SGRec")   source(paste(codePath,'01c_FourierSurrogates.r',sep=""))
if( RecRegime== "SGRec" ) { load(paste(dataPath,"FourrierSurrowgates200RecTS.Rdata",sep="")) ;SR<-Rsg; names(SR) <- ac((histMaxYr+1):(histMaxYr+dim(SR)[2]))}


            
#- Take one model for each iter
cuales                                <- sample(1:200,nits,replace=F)
SR                                    <- SR[cuales,]

  #-------------------------------------------------------------------------------
  # 5):  create replicates for the stock (numbers at age and Fs) resampling values from the corvarcov matrix
  #-------------------------------------------------------------------------------

#- Select the parameters which are need
which.pars            <- c("recruits","fmort",paste("log_selcoffs_fsh[",fisheries,"]",sep=""))
#- Read in the output of the assessment
vcovm                 <- read.table(paste(dataPath,"../jjm.cor",sep=""),skip=1,fill=TRUE,stringsAsFactors=F)
names(vcovm)          <- vcovm[1,]
vcovm                 <- vcovm[-1,]
#- Std devs of the parameters
mu.sd                 <- vcovm[,1:4]
mu.sd                 <- mu.sd[is.element(mu.sd$name,which.pars),]
mu.sd$index           <- an(mu.sd$index)
mu.sd$value           <- an(mu.sd$value)
mu.sd$std.dev         <- an(mu.sd$std.dev)
#- Parameter correlations
cormat                <- vcovm[,-c(1:4)]
cormat                <- cormat[mu.sd$index,mu.sd$index]
cormat[is.na(cormat)] <- 0
cormat                <-cormat + t(cormat)
diag(cormat)          <-1
dims                  <- dim(cormat)[1]
#- Calculate the varcov matrix from correlation matrix and std devs
sd.mat                <- data.frame(matrix(rep(mu.sd$std.dev,dims),nrow=dims))
varcov                <- sd.mat*cormat*t(sd.mat)



#- Resample parameters and recruitments using the varcov matrix
pars                  <- mvrnorm(nits, mu.sd$value, as.matrix(varcov))
recs                  <- pars[,mu.sd$name == "recruits"]
  idx                 <- which(recs <= 0,arr.ind=T)
if(length(idx)>0){
  for(i in 1:nrow(idx))
    recs[idx[i,1],idx[i,2]] <- sample(recs[,idx[i,2]][recs[,idx[i,2]]>0],1,replace=TRUE)
}

#- Need to remove the last point for "2013"
recs                  <- recs[,-dim(recs)[2]]
  # check : plot(recs[2,]) ;lines(stock.n(JMB)[1,,1])   :  OK

Fbar_fsh              <- pars[,mu.sd$name == "fmort"]
Fbar_fsh[Fbar_fsh<0] <- min(Fbar_fsh[Fbar_fsh>0])
Fbar_fsh              <- FLQuant(aperm(array(Fbar_fsh,dim=c(nits,length(years),length(fisheries),1,1,1)),c(6,2,5,4,3,1)),
                                 dimnames=dimnames(catch(JMB)[,histPeriod]))
  # check : plot(Fbar_fsh[1,1:43]) ;lines(Fbar[,,1])   :  OK

#- Get periods where the selection patterns of each fleet is constant
selPeriods            <- list()
for(iFishery in fisheries){
  elem                <- which(names(res)==paste("sel_fsh_",iFishery,sep=""))
  yrs                 <- res[[elem]][,2]; yrs <- c(yrs,rev(yrs)[1]+1)
  blcks               <- which(duplicated(apply(res[[elem]][,-c(1:2)],1,paste,collapse=""))==FALSE); blcks <- c(blcks,length(yrs))
  selPeriods[[iFishery]] <- mapply(seq,yrs[blcks[1:(length(blcks)-1)]],yrs[blcks[2:length(blcks)]]-1)
}

  #-------------------------------------------------------------------------------
  # 6): Update biological object with newly drawn values
  #-------------------------------------------------------------------------------

JMBOrig                               <- JMB

#- Add newly drawn recruits
dmns <- dimnames(JMB@stock.n)
for(iArea in dmns$area)
  stock.n(JMB)[1,histPeriod,,,iArea]  <- t(recs)

#- Scale harvest to newly drawn FBars
for(iArea in dmns$area){
  idx <- unique(which(quantMeans(JMB@harvest[ac(range(JMB)["minfbar"]:range(JMB)["maxfbar"]),histPeriod,,,iArea]) < 0.001,arr.ind=T)[,2])
  print(idx)
  if(length(idx) >= round(0.5*length(histPeriod))) stop("Newly drawn F's are too small")
  if(length(idx)>0){
    harvest(JMB)[,histPeriod[-idx],,,iArea]   <- sweep(harvest(JMB)[,histPeriod[-idx],,,iArea],2:6,Fbar_fsh[,-idx,,,iArea] / quantMeans(JMB@harvest[ac(range(JMB)["minfbar"]:range(JMB)["maxfbar"]),histPeriod[-idx],,,iArea]),"*")
  } else {
    harvest(JMB)[,histPeriod      ,,,iArea]   <- sweep(harvest(JMB)[,histPeriod      ,,,iArea],2:6,Fbar_fsh[,    ,,,iArea] / quantMeans(JMB@harvest[ac(range(JMB)["minfbar"]:range(JMB)["maxfbar"]),histPeriod      ,,,iArea]),"*")
  }
}
  
#- Update all ages stock.n
Z                                                               <- areaSums(harvest(JMB)) + areaMeans(m(JMB))
for(iAge in ages[-1])
  stock.n(JMB)[ac(iAge),ac((histMinYr+1):histMaxYr)]            <- sweep(stock.n(JMB)[ac(iAge-1),ac(histMinYr:(histMaxYr-1))],c(1:4,6),
                                                                         exp(- Z[ac(iAge-1),ac(histMinYr:(histMaxYr-1))]),"*")
stock.n(JMB)[ac(range(JMB)["max"]),ac((histMinYr+1):histMaxYr)] <- (stock.n(JMB)[ac(range(JMB)["max"]),ac((histMinYr+1):histMaxYr)] +
                                                                    sweep(stock.n(JMB)[ac(range(JMB)["max"]),ac(histMinYr:(histMaxYr-1))],c(1:4,6),
                                                                          exp(- Z[ac(range(JMB)["max"]),ac(histMinYr:(histMaxYr-1))]),"*"))


#- Update landings.n
landings.n(JMB)       <- sweep(harvest(JMB),c(1:4,6),(areaMeans(stock.n(JMB)) *
                               (1-exp(-areaSums(harvest(JMB))-areaMeans(m(JMB))))) /
                               (areaSums(harvest(JMB)) + areaMeans(m(JMB))),"*")
catch.n(JMB)          <- landings.n(JMB)

#- Update other slots
JMB@catch             <- quantSums(JMB@catch.n * JMB@catch.wt)
JMB@landings          <- quantSums(JMB@landings.n * JMB@landings.wt)
JMB@discards          <- quantSums(JMB@discards.n * JMB@discards.wt)
JMB@stock             <- quantSums(JMB@stock.n * JMB@stock.wt)
units(JMB)[1:17]<- as.list(c(rep(c("tonnes","thousands","kg"),4),
                             rep("NA",2),"f",rep("NA",2)))

  #-------------------------------------------------------------------------------
  # 7): Update fishery object with newly drawn values
  #-------------------------------------------------------------------------------

#- Add newly drawn selectivities
for(iFishery in fisheries){
  pSel        <- duplicated(landings.sel(fishery)[,1,,,iFishery,1])
  newSel      <- array(t(exp(pars[,mu.sd$name == paste("log_selcoffs_fsh[",iFishery,"]",sep="")])),
                       dim=c(length(which(pSel==FALSE)),length(selPeriods[[iFishery]]),nits))

  #- Define which selectivities are not unique and create repetition vector
  reps        <- rep(1,length(ages))
  for(iAge in an(dimnames(pSel)$age[-1]))
    reps[iAge]<- ifelse(pSel[iAge]==FALSE,reps[iAge-1]+1,reps[iAge-1])

  landings.sel(fishery)[,histPeriod,,,iFishery][]  <-
    newSel[reps,unlist(mapply(rep,1:length(selPeriods[[iFishery]]),lapply(selPeriods[[iFishery]],length))),]
  landings.sel(fishery)[,histPeriod,,,iFishery]    <- sweep(landings.sel(fishery)[,histPeriod,,,iFishery],2:6,
                                                         quantMeans(landings.sel(fishery)[,histPeriod,,,iFishery]),"/")
}

#- Set future selection patterns
#- Assume a constant selectivity for the future years for now
for(iFishery in fisheries)
  fishery@landings.sel[,projPeriod,,,iFishery] <- fishery@landings.sel[,ac(rep(histMaxYr,length(projPeriod))),,,iFishery]


  #-------------------------------------------------------------------------------
  # 8):  create the percieved stocks
  #-------------------------------------------------------------------------------

# only  the N and F at age, adding a cohort correlated noise, of the same magnitude as the uncertainty form the assessment
# no change in the rest (ie weights, catches etc...)

#- Creat the stock object (perceived dynamics)
JMS                 <-  JMB

#- Compute the stocks as observed in the last assessment year
#  the YC stochastic component of the deviation from the true stock
dmns                <- dimnames(stock.n(JMS))
devN                <- FLQuant(NA,dimnames=dmns)
devF                <- FLQuant(NA,dimnames=dmns)
devN[1,,,,1][]      <- rnorm(prod(dim(devN[,,,,1])[-1]),0,1)
devF[1,,,,] []      <- rnorm(prod(dim(devF)       [-1]),0,1)
for(iArea in dmns$area[-1])
  devN[,,,,iArea][] <- devN[,,,,1]
for(iAge in ages[-1]){
  devN[iAge,ac((histMinYr+1):futureMaxYr)] <- devN[iAge-1,ac((histMinYr):(futureMaxYr-1))]
  devF[iAge,ac((histMinYr+1):futureMaxYr)] <- devF[iAge-1,ac((histMinYr):(futureMaxYr-1))]
}
devN[is.na(devN)][] <- 0
devF[is.na(devF)][] <- 0

#- Compute the uncertainty on the N and F at age
CV_stock.n          <-  areaMeans(iterVars(stock.n(window(JMB,2004,histMaxYr)))^0.5/iterMeans(stock.n(window(JMB,2004,histMaxYr))))
CV_harvest          <-  iterVars(harvest(window(JMB,2004,histMaxYr))^0.5/iterMeans(harvest(window(JMB,2004,histMaxYr))))

#- the year effect on the amplitude of the error : historical CV from the assessment output (comes from the resampling from the varcov done on step 5)
# assume that the assessment errors have the same pattern as the assessment uncertainty (reflected by the varcov matrix)

dimnames(CV_stock.n)$area <- 1
CV_stock.n      <- expand(CV_stock.n,area=fisheries)
for(iFisheries in 2:length(fisheries))
  CV_stock.n[,,,,iFisheries] <- CV_stock.n[,,,,1]
CV_stock.n      <- expand(CV_stock.n,year=ac(histMinYr:futureMaxYr)) #- Warning is OK
CV_stock.n      <- propagate(CV_stock.n,iter=nits)
CV_stock.n[,ac(histMinYr:2003)][] <- 0


CV_harvest      <- expand(CV_harvest,year=ac(histMinYr:futureMaxYr)) #- Warning is OK
CV_harvest      <- propagate(CV_harvest,iter=nits)
CV_harvest[,ac(histMinYr:2003)][] <- 0

#- Stock : biological object + deviations
#  with deviation = cohort coherent deviation (devN) and year dependent amplitude (JMB@stock.n*CV_stock.n)
stock.n(JMS)    <-  JMB@stock.n + devN*JMB@stock.n*CV_stock.n
harvest(JMS)    <-  JMB@harvest + devF*JMB@harvest*CV_harvest

#- Don't allow for negative F and N values
stock.n(JMS)[stock.n(JMS)<=0 & !is.na(stock.n(JMS))]  <- 10
harvest(JMS)[harvest(JMS)<0 & !is.na(harvest(JMS))]   <-  0

  #-------------------------------------------------------------------------------
  # 9):  clean and save
  #-------------------------------------------------------------------------------

save(JMB,         file=file.path(outPath,"biol.RData"),compress=T)
save(JMBOrig,     file=file.path(outPath,"biolOrig.RData"),compress=T)
save(JMS,         file=file.path(outPath,"stock.RData"),compress=T)
save(fishery,     file=file.path(outPath,"fishery.RData"),compress=T)
save(CV_stock.n,  file=file.path(outPath,"CVstockn.RData"),compress=T)
save(CV_harvest,  file=file.path(outPath,"CVharvest.RData"),compress=T)
save(devN,        file=file.path(outPath,"devN.RData"),compress=T)
save(devF,        file=file.path(outPath,"devF.RData"),compress=T)
save(SR,          file=file.path(outPath,"SR.RData"),compress=T)
save(settings,    file=file.path(outPath,"settings.RData"),compress=T)


