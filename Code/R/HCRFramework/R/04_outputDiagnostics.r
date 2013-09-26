#-------------------------------------------------------------------------------
# WKHELP
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
# Performs an MSE of North Sea Herring under different TAC scenario's
#
# Date: 02-Sep-2012
#
# Build for R2.13.2
#-------------------------------------------------------------------------------

rm(list=ls())

library(FLSAM)
library(MASS)
library(msm)
wine <- F

library(FLCore)
library(PBSadmb)
library(lattice)
library(MASS)

ac<-function(x) {return(as.character(x))}
an<-function(x) {return(as.numeric(x))}

codePath        <- "N:/Projecten/SouthPacific/WP4/R/code/submit2GitHub/R/"
dataPath        <- "N:/Projecten/SouthPacific/WP4/R/code/submit2GitHub/Data/"
outPath         <- "N:/Projecten/SouthPacific/WP4/R/code/submit2GitHub/Results/"


setwd(Pathscript)

source(paste(Pathscript,'functions_NTH.r',sep=""))

# load the true and observed stocks at the start of the simulation  from :
RecRegime <-  "LTRec"

#define year ranges
ShortT    <-  ac(2013:2017)
MidT      <-  ac(2018:2027)
LongT     <-  ac(2027:2040)


PathPlots1   <- paste(PathRes,RecRegime,"/Plots input/",sep="")



##-------------------------------------------------------------------------------
## Setup array to save results
##-------------------------------------------------------------------------------



diags<-data.frame(Hcr=NA,RecRegime=NA,Btrigger=NA,Blim=NA,Ftarget=NA,
      Risk2ShortT=NA,Risk2MidT=NA,Risk2LongT=NA,RecovTimeBlim=NA,RecovTimeBpa=NA,
      SSBend=NA,meanSSBLongT=NA,
      Fend=NA,meanF=NA,meanFLongTerm=NA,
      meanYieldShortT=NA,meanYieldMidT=NA,meanYieldLongT=NA,
      meanrelTACIAV=NA,noIAVrestrictup=NA,noIAVrestrictdown=NA,TACup=NA,TACdown=NA)                                                                                                  


##-------------------------------------------------------------------------------
## Load results
##-------------------------------------------------------------------------------

HCRs<-c( " Hockey stick for report" ,
            " Recov for report",
            " Ftarget")


counter <- 1
for(Hcr in HCRs){

load(file=paste(PathRes,RecRegime,"/Hcr",Hcr,"/simulations results with",RecRegime,".RData",sep=""))
projPeriod           <- 2013:(futureMaxYr-1)  
print(counter)
   

diags[counter,"Hcr"]   <- Hcr
diags[counter,"RecRegime"]   <- RecRegime
diags[counter,"Btrigger"]   <-Ref$Btrigger
diags[counter,"Blim"]   <-Ref$Blim
diags[counter,"Ftarget"]   <-Ref$Ftarget


##-------------------------------------------------------------------------------
## Diagnostics on results
##-------------------------------------------------------------------------------
Btrg  <- Ref$Btrigger
Blim  <- Ref$Blim
Bpa  <- Ref$Bpa


JMssb<-ssb(JMB)
JM2ssb<-ssb(JMSstore)
Fbar<-quantMeans(areaSums(harvest(JMB)))
Fbar2<-quantMeans(areaSums(harvest(JMSstore)))

# risk related to Blim
diags[counter,"Risk2ShortT"]             <- length(unique(which(JMssb[,ShortT]<Blim,arr.ind=T)[,"dim6"]))/dims(JMssb)$iter  *100       # percentage of iteration that reach Blim
diags[counter,"Risk2MidT"]               <- length(unique(which(JMssb[,MidT]  <Blim,arr.ind=T)[,"dim6"]))/dims(JMssb)$iter  *100
diags[counter,"Risk2LongT"]              <- length(unique(which(JMssb[,LongT] <Blim,arr.ind=T)[,"dim6"]))/dims(JMssb)$iter  *100
# time before recovery 
diags[counter,"RecovTimeBlim"]        <- round(iterMeans( yearSums(JMssb[,ac(projPeriod)]<Blim)),2)
diags[counter,"RecovTimeBpa"]         <- round(iterMeans( yearSums(JMssb[,ac(projPeriod)]<Bpa)),2)

# stock and fishing mortality
diags[counter,"SSBend"]               <- round(iterMeans(JMssb[,ac(futureMaxYr-1)]))
diags[counter,"meanSSBLongTerm"]      <- round(median(c(apply(JMssb[,LongT],3:6,mean,na.rm=T))))

diags[counter,"meanF"]                <- round(  median(  yearMeans(quantMeans(areaSums(JMB@harvest[,ac(projPeriod)])))@.Data) ,3)
diags[counter,"meanFLongTerm"]        <- round( median(   yearMeans(quantMeans(areaSums(JMB@harvest[,LongT])))@.Data) ,3)
diags[counter,"Fend"]                 <- round( median   (quantMeans(areaSums(JMB@harvest[,ac(futureMaxYr-1)]))@.Data) ,3)

# difference between percieved and true stocks
diags[counter,"SSBabsBias"]           <-  round(apply( yearMeans(100*abs(JM2ssb[,ac(projPeriod)]-JMssb[,ac(projPeriod)])/JMssb[,ac(projPeriod)]),1:5,median,na.rm=T),3)
diags[counter,"SSBBias"]              <-  round(apply( yearMeans(100*(JM2ssb[,ac(projPeriod)]-JMssb[,ac(projPeriod)])/JMssb[,ac(projPeriod)]),1:5,median,na.rm=T),3)
diags[counter,"FbarabsBias"]          <-  round(apply( yearMeans(100*abs((Fbar2[,ac(projPeriod)])-(Fbar[,ac(projPeriod)]))/(Fbar[,ac(projPeriod)])),1:5,median,na.rm=T),3)
diags[counter,"FbarBias"]             <-  round(apply( yearMeans(100*((Fbar2[,ac(projPeriod)])-(Fbar[,ac(projPeriod)]))/(Fbar[,ac(projPeriod)])),1:5,median,na.rm=T),3)



# catches and quotas
diags[counter,"meanYieldShortTerm"]   <- round(median(c(yearMeans((computeCatch(JMB[,ShortT]))))))
diags[counter,"meanYieldMidTerm"]     <- round(median(c(yearMeans((computeCatch(JMB[,MidT]))))))
diags[counter,"meanYieldLongTerm"]    <- round(median(c(yearMeans((computeCatch(JMB[,LongT]))))))




# quota variability
diags[counter,"meanrelTACIAV"]        <- round(median(c(apply(abs(TAC[,ac(projPeriod[2]:rev(projPeriod)[1])] - TAC[,ac(projPeriod[1]:rev(projPeriod)[2])]) /TAC[,ac(projPeriod[2]:rev(projPeriod)[1])] * 100,3:6,mean,na.rm=T))),3)
IAVUp   <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)])] == 1.15* TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],arr.ind=T)
IAVDown <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)])] == 0.85* TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],arr.ind=T)
#  #- Average number of times the IAV rule is applied upwards or downwards
diags[counter,"noIAVrestrictup"]<-0
if((nrow(IAVUp)) > 0 ){
  a <- IAVUp
  diags[counter,"noIAVrestrictup"]    <- max(0,median(aggregate(a[,"dim2"],by=list(a[,"dim6"]),function(x){length(x)})$x),na.rm=T)
}
 diags[counter,"noIAVrestrictdown"]<-0
if((nrow(IAVDown)) > 0 ){
  a <- IAVDown
  diags[counter,"noIAVrestrictdown"]  <- max(0,median(aggregate(a[,"dim2"],by=list(a[,"dim6"]),function(x){length(x)})$x),na.rm=T)
}
#
#  #- Which TAC of the runs go up and which go down
resUp   <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)])] > TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],arr.ind=T)
resDown <- which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)])] < TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-1])],arr.ind=T)
#
#  #- Mean increase in TAC is TAC goes up, or mean decrease in TAC is TAC goes down
diags[counter,"TACup"]   <- round(mean((TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1])] - TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2])])@.Data[which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1])] > TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2])])]))
diags[counter,"TACdown"] <- round(mean((TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2])] - TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1])])@.Data[which(TAC[,ac(projPeriod[2]:projPeriod[length(projPeriod)-1])] < TAC[,ac(projPeriod[1]:projPeriod[length(projPeriod)-2])])]))
#




counter <- counter + 1
}


write.csv(t(diags),file=paste(PathRes,RecRegime,"/tables_diags.csv",sep=""),row.names=T)
