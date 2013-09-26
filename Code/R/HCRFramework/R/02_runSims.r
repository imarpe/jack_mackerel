#-------------------------------------------------------------------------------
#
# Script runs the simulation projecting the stocks foreward applying
#  the harvest control rule
#
# Created by: Thomas Brunel, Niels Hintzen
# Project: South Pacific
#
#-------------------------------------------------------------------------------

rm(list=ls())

library(FLCore)
library(PBSadmb)
library(lattice)
library(MASS)

#- Set paths
codePath        <- "N:/Projecten/SouthPacific/WP4/R/code/submit2GitHub/R/"
dataPath        <- "N:/Projecten/SouthPacific/WP4/R/code/submit2GitHub/Data/"
outPath         <- "N:/Projecten/SouthPacific/WP4/R/code/submit2GitHub/Results/"


source(file.path(codePath,"functions.r"))

#- Load objects
load(file=file.path(outPath,"biol.RData"))
load(file=file.path(outPath,"stock.RData"))
load(file=file.path(outPath,"fishery.RData"))
load(file=file.path(outPath,"CVstockn.RData"))
load(file=file.path(outPath,"CVharvest.RData"))
load(file=file.path(outPath,"devN.RData"))
load(file=file.path(outPath,"devF.RData"))
load(file=file.path(outPath,"SR.RData"))
load(file=file.path(outPath,"settings.RData"))
load(file=file.path(outPath,"yearParams.RData"))
load(file=file.path(outPath,"dimensions.RData"))


#-  Define management scenario and the reference points which are going to be used for the HCR
#ReF   <- read.csv (paste(PathData,"Refpoints",RecRegime,".csv",sep=""))
Ref   <- list(Bmsy=10556,Btrigger=10566,Bpa=3900,Blim=2800,Fmsy=0.15,Ftarget=0.15,Fmin=0.05)
HCR   <- list(IAV=F,            # whether or not interannual variability of TAC should be limited
              percIAV=15,       # to this given percentage
              Btarget=NULL,     #  give a value to manage towards a Btaraget value
              Ftarget=0.15,     # to manage at constant F
              slope=NULL,       # set to TRUE to have a HCR with slope
              slopeBmin=00,     # B value under which F is contant at Fmin
              slopeBmax=10556,  # B value over which F is constant at Fmax
              slopeFmin=0.05,
              slopeFmax=0.15,
              Recov=NULL        # if TRUE a recovery plan is implemented by decreasing F by 25% annually until B>Blim, then decreasing by 15% until B>Bpa, then F=Fmsy
              )


  #------------------------------------------------------------------------------#
  # 1) create objects, make an assumtion the true F for the first ImY, chose the HCR
  #------------------------------------------------------------------------------#

#- To keep track of the advised TAC
TAC             <- areaSums(catch(JMB))
TAC[,ac(2013)]  <- 360
#- Tp keep track of the perceived stock in the terminal asse'ssment Year
JMSstore        <- JMS



  #------------------------------------------------------------------------------#
  # 2) Start running
  #------------------------------------------------------------------------------#

start.time <- Sys.time()
for (iYr in an(projPeriod)){
  cat(iYr,"\n")
  cat(paste("\n Time running",round(difftime(Sys.time(),start.time,unit="mins"),0),"minutes \n"))

  #----------------------------------------
  # define year names
  #----------------------------------------
  TaY <- ac(iYr-1)  # terminal year in the assessment
  ImY <- ac(iYr)    # intermediate year in the short term forecast, ie, current year
  FcY <- ac(iYr+1)  # year for which the advice is given, ie, forecast two years ahead the last year in the assessment

  #----------------------------------------
  # update the true stock number at age in ImY
  #----------------------------------------

  #- age 1 : recruits
  SSBTaY                              <- ssb(JMB[,TaY])
  stock.n(JMB)[1,ImY,,,1][]           <- B2R(SSBTaY,SR,ImY,settings$RecRegime)
  for(iArea in 2:length(fisheries))
    stock.n(JMB)[1,ImY,,,iArea]       <- stock.n(JMB)[1,ImY,,,1]

  #- survivors
  stock.n(JMB)[,ImY]                  <- calcSurvivors(JMB[,c(TaY,ImY)])

  #- copy to all areas (note that there is only one stock-number-at-age matrix,
  #   but to make calculations easy, we copy it 4 times to match dimensions of area)
  for(iArea in 2:length(fisheries))
    stock.n(JMB)[,ImY,,,iArea]        <- stock.n(JMB)[,ImY,,,1]

  #- Update biological, stock and fisheries object
  fishery@landings.n[,TaY]            <- sweep(harvest(JMB[,TaY]),c(1:4,6),(areaSums(harvest(JMB[,TaY])) + areaMeans(m(JMB[,TaY]))) *
                                               areaMeans(stock.n(JMB[,TaY])) * (1-exp(-areaSums(harvest(JMB[,TaY]))-areaMeans(m(JMB[,TaY])))),"*")
  landings.n(JMB[,TaY])               <- fishery[,TaY]@landings.n
  fishery@landings[,TaY]              <- quantSums(fishery[,TaY]@landings.n * fishery[,TaY]@landings.wt)
  landings(JMB[,TaY])                 <- landings(fishery[,TaY])
  catch.n(JMB[,TaY])                  <- landings.n(JMB[,TaY]) + discards.n(JMB[,TaY])
  JMB@catch[,TaY]                     <- computeCatch(JMB[,TaY])
  stock(JMB[,TaY])                    <- quantSums(JMB@stock.n[,TaY] * JMB@stock.wt[,TaY])
  discards(JMB[,TaY])                 <- quantSums(JMB@discards.n[,TaY] * JMB@discards.wt[,TaY])

  #----------------------------------------
  #- Perform an assessment
  #----------------------------------------
  
  #- Create a new perceived stock by running a pseudo assessment
  JMS                                 <- assessJM(JMB,devN,devF,CV_stock.n,CV_harvest,iYr) #- calculation of TaY inside the function
  #- Keep track of the assessment output for the terminal assessmnet year
  JMSstore[,TaY]                      <- JMS[,TaY]

  #----------------------------------------
  #- Do a short term forcast (stf)
  #----------------------------------------
  if(iYr !=  futureMaxYr){
    stf                               <- JMS[,c(TaY,ImY,FcY)]

    #- Assume same harvest pattern for ImY
    harvest(stf)[,ImY]                <- harvest(JMS)[,TaY]

    #--------------------------------------
    # Update to Intermediate year
    #--------------------------------------

    #- Rec ImY
    stock.n(stf)[1,ImY]               <- exp(yearMeans(log(stock.n(JMS)[1,ac((an(TaY)-9):an(TaY))])))

    #- survivors ImY
    stock.n(stf)[,ImY]                <- calcSurvivors(stf[,c(TaY,ImY)])

    #- copy to all areas
    for(iArea in 2:length(fisheries))
      stock.n(stf)[,ImY,,,iArea]      <- stock.n(stf)[,ImY,,,1]

    #- Scale harvest pattern to the TAC in the ImY
    scaleHarvest  <- numeric(nits)
    for(iTer in 1:nits)
      scaleHarvest[iTer]  <- optimize(f=function(x,tac,stck){ return(
                                         sqrt((tac -
                                              sum(sweep(sweep(stck@harvest[,drop=T] * x,1,(rowSums(stck@harvest[,drop=T],na.rm=T) * x) + rowMeans(stck@m[,drop=T],na.rm=T),"/") *
                                                                  stck@stock.n[,drop=T],1,
                                                                  (1-exp(-rowSums(stck@harvest[,drop=T] * x,na.rm=T) - rowMeans(stck@m[,drop=T],na.rm=T))),"*") * stck@catch.wt[,drop=T],na.rm=T))^2))},
                                      interval=c(0,10),tac=c(iter(TAC[,ImY],iTer)),stck=iter(stf[,ImY],iTer))$minimum
    harvest(stf)[,ImY]    <- sweep(harvest(stf)[,ImY],6,scaleHarvest,"*")

    #--------------------------------------
    # Update to Forecast year
    #--------------------------------------

    #- Rec FcY
    stock.n(stf)[1,FcY]               <- stock.n(stf)[1,ImY]

    #- survivors ImY
    stock.n(stf)[,FcY]                <- calcSurvivors(stf[,c(ImY,FcY)])

    #- copy to all areas
    for(iArea in 2:length(fisheries))
      stock.n(stf)[,FcY,,,iArea]      <- stock.n(stf)[,FcY,,,1]

    #--------------------------------------
    # Find TAC according to HCR
    #--------------------------------------

    TACtemp                           <- findTAC(stf,Ref,HCR,TaY,ImY,FcY)

    # apply the TAC variation limitation
    if(HCR$IAV){
      idxmin                          <- which(TACtemp < ((1-HCR$percIAV/100) * TAC[,ImY]),arr.ind=T)[,6]
      idxmax                          <- which(TACtemp > ((1+HCR$percIAV/100) * TAC[,ImY]),arr.ind=T)[,6]
      if(length(idxmin)>0)
        TACtemp[,,,,,idxmin]          <- (1-HCR$percIAV/100) * TAC[,ImY,,,,idxmin]
      if(length(idxmax)>0)
        TACtemp[,,,,,idxmax]          <- (1+HCR$percIAV/100) * TAC[,ImY,,,,idxmax]
    }
    TAC[,FcY]                         <- TACtemp

    #- Update harvest pattern in ImY for biological object
    JMB                               <- TAC2F(JMB,fishery,TaY,ImY,FcY,TAC)
  }
#save.image(file=paste(outPath,settings$RecRegime,"/HCR Ftarget/simulations results with",settings$RecRegime,".RData",sep=""))
#save.image(file=paste(outPath,settings$RecRegime,"/HCR Hockey stick for report/simulations results with",settings$RecRegime,".RData",sep=""))

#save.image(file=paste(outPath,settings$RecRegime,"/HCR Recov for report/simulations results with",settings$RecRegime,".RData",sep=""))
}  # end of years loop


