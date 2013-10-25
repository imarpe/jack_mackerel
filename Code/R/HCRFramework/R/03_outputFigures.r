#-------------------------------------------------------------------------------
# WKHERMPII
#
# Author: Niels Hintzen
#         IMARES, The Netherland
#
# Performs an MSE of North Sea Herring under different TAC scenario's
#
# Date: 03-Oct-2011
#
# Build for R2.8.1
#-------------------------------------------------------------------------------

rm(list=ls())

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

source(paste(Pathscript,"functions_NTH.r",sep=""))

# load the true and observed stocks at the start of the simulation  from :
#source('N:/Projecten/SouthPacific/WP4/R/code/1_create Objects02.r')
RecRegime <-  "LTRec"
Hcr <-  " Hockey stick for report"          #  simulations results withLTRec
#Hcr <-  " Recov for report"
#Hcr <-  " Ftarget"

load(file=paste(PathRes,RecRegime,"/Hcr",Hcr,"/simulations results with",RecRegime,".RData",sep=""))


PathPlots2   <- paste(PathRes,RecRegime,"/Hcr",Hcr,"/Plots/",sep="")

# source(paste(codePath,"functions.r",sep=""))

#-------------------------------------------------------------------------------
# Figures on results
#-------------------------------------------------------------------------------
JMBssb<-ssb(JMB)
JMSssb<-ssb(JMS)

rSSBp <- apply(JMBssb@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rSSBs <- apply(JMSssb@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)

#rFstf <- apply(apply(unitSums(fSTF)[ac(2:6),],2:6,mean,na.rm=T)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)

  #- Plot settings
yrs   <- 2001:(futureMaxYr-1)
cl    <- 1.2
ca    <- 1.1
fonts <- 2


#-------------------------------------------------------------------------------
# 3): Plot results of SSB stock and SSB pop
#-------------------------------------------------------------------------------
par(mar=c(5.1,5.1,4.1,2.1))
xrange  <- range(yrs)
yrange  <- c(0,range(pretty(c(rSSBp[,,ac(yrs),,,],rSSBs[,,ac(yrs),,,])),na.rm=T)[2])
  #---------
  #- Biology
  #---------
plot(rSSBp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="Years",ylab="",
     type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
mtext(text="SSB (million tonnes)",side=2,at=0.5,outer=T,cex=cl,line=-1)
grid(); box()
#-Reference level Blim & Bpa
Btrg  <- Ref$Btrigger
Blim  <- Ref$Blim
Bpa  <- Ref$Bpa


abline(h=Blim,col="red",lwd=2,lty=2);
abline(h=Bpa,col="blue",lwd=2,lty=2);
abline(h=Btrg,col="darkgreen",lwd=2,lty=2);
mtext(text="Blim",side=4,at=Blim,las=1,cex=0.8,col="red",line=0.5,font=fonts)
mtext(text="Bpa",side=4,at=Bpa,las=1,cex=0.8,col="blue",line=0.5,font=fonts)
mtext(text="Btrigger",side=4,at=Btrg,las=1,cex=0.8,col="darkgreen",line=0.5,font=fonts)
lines(rSSBp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rSSBp["5%",,ac(yrs),,,]~yrs,lty=3)
lines(rSSBp["95%",,ac(yrs),,,]~yrs,lty=3)

  #---------
  #- Stock
  #---------
lines(rSSBs["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2,col="red")
lines(rSSBs["5%",,ac(yrs),,,]~yrs,lty=3,col="red")
lines(rSSBs["95%",,ac(yrs),,,]~yrs,lty=3,col="red")

  #---------
  #- Legend
  #---------

legend("bottomright",legend=c("SSB assessed stock","SSB true population"),
       col=c("red","black"),lwd=3,lty=1,box.lty=0)
savePlot(file=paste(PathPlots2,"SSB.png",sep=""),type="png");dev.off()
#-------------------------------------------------------------------------------
# 4): Plot results of F stock and true F
#-------------------------------------------------------------------------------
Fbar<-fbar(JMB)
Fbar2<-fbar(JMSstore)

rFp   <- apply(unitSums(Fbar)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
rFs   <- apply(unitSums(Fbar2)@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)


par(mar=c(5.1,5.1,4.1,2.1))
xrange  <- range(yrs)
yrange  <- c(0,range(pretty(c(rFp[,,ac(yrs),,,],rFs[,,ac(yrs),,,])),na.rm=T)[2])

  #---------
  #- Biology
  #---------
plot(rFp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="Years",ylab="",
     type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex=ca,font=fonts)
mtext(text="Fishing mortality (all ages)",side=2,at=0.5,outer=T,cex=cl,line=-1)
grid(); box()
#-Reference level Fpa
Fmsy<-Ref$Fmsy
Fmsy<-0.145
abline(h=Fmsy,col="darkgreen",lwd=2,lty=2);
mtext(text="Fmsyt",side=4,at=Fmsy,las=1,cex=0.8,col="darkgreen",line=0.5,font=fonts)
lines(rFp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
lines(rFp["5%",,ac(yrs),,,]~yrs,lty=3)
lines(rFp["95%",,ac(yrs),,,]~yrs,lty=3)

  #---------
  #- Stock
  #---------
lines(rFs["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2,col="red")
lines(rFs["5%",,ac(yrs),,,]~yrs,lty=3,col="red")
lines(rFs["95%",,ac(yrs),,,]~yrs,lty=3,col="red")

  #---------
  #- STF
  #---------
#lines(rFstf["50%",,ac(2012:rev(projPeriod)[2]),,,]~ac(2012:rev(projPeriod)[2]),lty=1,lwd=2,col="blue")
#lines(rFstf["5%",,ac(2012:rev(projPeriod)[2]),,,]~ac(2012:rev(projPeriod)[2]),lty=3,col="blue")
#lines(rFstf["95%",,ac(2012:rev(projPeriod)[2]),,,]~ac(2012:rev(projPeriod)[2]),lty=3,col="blue")
#
  #---------
  #- Legend
  #---------

legend("bottomright",legend=c("F assessed stock","F true population"),#,"F short term forecast"),
       col=c("red","black","blue"),lwd=3,lty=1,box.lty=0)
savePlot(file=paste(PathPlots2,"F.png",sep=""),type="png");dev.off()




#-------------------------------------------------------------------------------
# 5): Plot results of landings by fleet & TAC on top
#-------------------------------------------------------------------------------

rLandf<- apply(computeCatchperFleet(JMB),1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)
totY  <- apply(((computeCatch(JMB)))@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T) 
#rLands<- apply((unitSums(computeCatch(JM2save)))@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)

rTAC  <- apply(TAC@.Data,1:5,quantile,probs=c(0.05,0.5,0.95),na.rm=T)


par(mar=c(2.1,2.1,2.1,2.1),mfrow=c(2,2),oma=c(3,3,2,0))

  #---------
  #- Fishery
  #---------
for(iFsh in dimnames(fishery@landings.wt)$area){
  xrange  <- range(yrs)
  yrange  <- c(0,range(pretty(c(rLandf[,,ac(yrs),,,iFsh],rTAC[,,ac(yrs),,,])),na.rm=T)[2])

    #- Landings
  plot(rLandf["50%",,ac(yrs),,,iFsh]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
       type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
  axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
  mtext(text="Landings (thousand tonnes)",side=2,at=0.5,outer=T,cex=cl,line=1.5)
  mtext(text="Years",side=1,at=0.5,outer=T,cex=cl,line=1.5)
  grid(); box()
  lines(rLandf["50%",,ac(yrs),,,iFsh]~yrs,lty=1,lwd=2)
  lines(rLandf["5%",,ac(yrs),,,iFsh]~yrs,lty=3)
  lines(rLandf["95%",,ac(yrs),,,iFsh]~yrs,lty=3)
    #- TAC
  lines(rTAC["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2,col="red")
  lines(rTAC["5%",,ac(yrs),,,]~yrs,lty=3,col="red")
  lines(rTAC["95%",,ac(yrs),,,]~yrs,lty=3,col="red")


  text(x=2020,y=yrange[2],adj=c(1,1.5),labels=paste("(Fleet ",iFsh,")",sep=""),cex=0.8,font=fonts)
}
savePlot(file=paste(PathPlots2,"CatchTAC.png",sep=""),type="png");dev.off()
#
##-------------------------------------------------------------------------------
## 6): Plot results total landings fishery and stock
##-------------------------------------------------------------------------------
#par(mar=c(5.1,5.1,4.1,2.1))
#xrange  <- range(yrs)
#yrange  <- c(0,range(pretty(c(rLandf[,,ac(yrs),,,],rLands[,,ac(yrs),,,])),na.rm=T)[2])
#
#  #---------
#  #- Biology
#  #---------
#plot((rLandf["50%",,ac(yrs),,,])~yrs,xlim=xrange,ylim=yrange,xlab="Years",ylab="",
#     type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
#axis(2,las=1,at=pretty(hrayrange),labels=pretty(yrange),cex=ca,font=fonts)
#mtext(text="Landings (tonnes)",side=2,at=0.5,outer=T,cex=cl,line=-1)
#grid(); box()
#lines((rLandf["50%",,ac(yrs),,,])~yrs,lty=1,lwd=2)
#lines((rLandf["5%",,ac(yrs),,,])~yrs,lty=3)
#lines((rLandf["95%",,ac(yrs),,,])~yrs,lty=3)
#
#  #---------
#  #- Stock
#  #---------
#lines(rLands["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2,col="red")
#lines(rLands["5%",,ac(yrs),,,]~yrs,lty=3,col="red")
#lines(rLands["95%",,ac(yrs),,,]~yrs,lty=3,col="red")
#
#  #---------
#  #- Legend
#  #---------
#
#legend("bottomright",legend=c("Catch assessed stock","Catch true population"),
#       col=c("red","black"),lwd=3,lty=1,box.lty=0)
#savePlot(file=paste(PathPlots2,"figures/",scen,"_Catch.png",sep=""),type="png");dev.off()
##
#
##-------------------------------------------------------------------------------
## 7): Plot results of TAC
##-------------------------------------------------------------------------------
#
#par(mar=c(2.1,2.1,2.1,2.1),mfrow=c(2,2),oma=c(3,3,2,0))
#
#  #---------
#  #- Fishery
#  #---------
#for(iFsh in dimnames(fishery@landings.wt)$unit){
#  xrange  <- range(yrs)
#  yrange  <- c(0,range(pretty(rTAC[,,ac(yrs),iFsh,,]),na.rm=T)[2])
#
#  plot(rTAC["50%",,ac(yrs),iFsh,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
#       type="l",lwd=2,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n")
#  axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
#  mtext(text="TAC (thousand tonnes)",side=2,at=0.5,outer=T,cex=cl,line=1.5)
#  mtext(text="Years",side=1,at=0.5,outer=T,cex=cl,line=1.5)
#  grid(); box()
#  lines(rTAC["50%",,ac(yrs),iFsh,,]~yrs,lty=1,lwd=2)
#  lines(rTAC["5%",,ac(yrs),iFsh,,]~yrs,lty=3)
#  lines(rTAC["95%",,ac(yrs),iFsh,,]~yrs,lty=3)
#  text(x=2020,y=yrange[2],adj=c(1,1.5),labels=paste("(Fleet ",iFsh,")",sep=""),cex=0.8,font=fonts)
#}
#savePlot(file=paste(PathPlots2,"figures/",scen,"_TAC.png",sep=""),type="png");dev.off()
#
#
##-------------------------------------------------------------------------------
# 8): Plot trajectories of ssb(biol), rec(biol), f(biol), TAC(A)
#-------------------------------------------------------------------------------
#
#par(mfrow=c(2,2),mar=c(3,3,3,3),oma=c(3.1,3.1,1,1))
#
##- Plot ssb based on biol
#xrange  <- range(yrs)
#yrange  <- c(0,range(pretty(range(iter(JMssb[,ac(yrs)],1:5)/1000)))[2])
#plot(yrs,iter(JMssb[,ac(yrs)],1)/1000,type="l",xlab="Years",ylab="",
#     xlim=xrange,ylim=yrange,lwd=2,font=fonts,cex.lab=cl,cex.axis=ca,las=1)
#mtext(side=2,at=yrange[2]/2,text="Spawning Stock Biomass (kt)",outer=F,line=4,font=fonts,cex=ca)
#grid(); box()
#for(i in 1:10) lines(c(iter(ssbb(biol[,ac(yrs)],unitSums(f[,ac(yrs)])),i))/1000~yrs,col=i,lwd=2)
#
##- Plot rec based on biol
#xrange  <- range(yrs)
#yrange  <- c(0,pretty(range(c(iter(stock.n(JM)[1,ac(yrs)],1:10))/1e6),1)[3])
#plot(c(iter(stock.n(JM)[1,ac(yrs)],1))/1e6~yrs,type="l",xlab="Years",ylab="",
#     xlim=xrange,ylim=yrange,lwd=2,font=fonts,cex.lab=cl,cex.axis=ca,las=1)
#mtext(side=2,at=yrange[2]/2,text="Recruitment (millions)",outer=F,line=4,font=fonts,cex=ca)
#grid();box()
#for(i in 1:10) lines(c(iter(stock.n(JM)[1,ac(yrs)],i))/1e6~yrs,col=i,lwd=2)
#
##- f-ages 2-6
#xrange  <- range(yrs)
#yrange  <- c(0,pretty(range(c(iter(quantMeans(unitSums(f[ac(2:6),ac(yrs)])),1:10))),1)[3])
#plot(c(iter(quantMeans(unitSums(f[ac(2:6),ac(yrs)])),1))~yrs,type="l",xlab="Years",ylab="",
#     xlim=xrange,ylim=yrange,lwd=2,font=fonts,cex.lab=cl,cex.axis=ca,las=1)
#mtext(side=2,at=yrange[2]/2,text="Fishing mortality",outer=F,line=4,font=fonts,cex=ca)
#grid();box()
#for(i in 1:10) lines(c(iter(quantMeans(unitSums(f[ac(2:6),ac(yrs)])),i))~yrs,col=i,lwd=2)
#
##- TAC of fleet A
#xrange  <- range(yrs)
#yrange  <- c(0,pretty(range(c(iter(TAC[,ac(yrs),"A"]/1000,1:10))),1)[3])
#plot(c(iter(TAC[,ac(yrs),"A"]/1000,1))~yrs,type="l",xlab="Years",ylab="",
#     xlim=xrange,ylim=yrange,lwd=2,font=fonts,cex.lab=cl,cex.axis=ca,las=1)
#mtext(side=2,at=yrange[2]/2,text="TAC (kt)",outer=F,line=4,font=fonts,cex=ca)
#grid();box()
#for(i in 1:10) lines(c(iter(TAC[,ac(yrs),"A"]/1000,i))~yrs,col=i,lwd=2)
#
#savePlot(file=paste(PathPlots2,"figures/",scen,"_iterations.png",sep=""),type="png");dev.off()
#
#}
#
##
##-------------------------------------------------------------------------------
## 9): Report figures
##-------------------------------------------------------------------------------
#
#yrs   <- 2010:2020
#cl    <- 1.1
#ca    <- 1
#fonts <- 1
#yrangeSSB <- c(0,2.7e6)
#yrangeLan <- c(0,7e5)
#yrangeLan2<- c(0,2.5e4)
#yrangeF   <- c(0,0.4)
#
#
#par(mfrow=c(3,1),oma=c(6,6,2,3),mar=c(1,0,0,0))
#
#  #---------
#  #- Landings
#  #---------
#
#xrange  <- range(yrs)
#yrange  <- yrangeLan
#
#plot(rLandf["50%",,ac(yrs),1,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
#     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
#axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
#mtext(text=expression(paste("Landings (",10^3," tonnes)",sep="")),side=2,at=(yrange[2]-yrange[1])/2+yrange[1],outer=F,cex=cl,line=4,font=fonts)
#grid(); box()
#lines(rLandf["50%",,ac(yrs),1,,]~yrs,lty=1,lwd=2)
#lines(rLandf["5%",,ac(yrs),1,,]~yrs,lty=3,lwd=1)
#lines(rLandf["95%",,ac(yrs),1,,]~yrs,lty=3,lwd=1)
#par(new=T)
#yrange  <- yrangeLan2
#plot(rLandf["50%",,ac(yrs),2,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
#     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n",col="red")
#axis(4,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts,col="red",col.ticks="red",col.axis="red")
#lines(rLandf["50%",,ac(yrs),2,,]~yrs,lty=1,lwd=2,col="red")
#lines(rLandf["5%",,ac(yrs),2,,]~yrs,lty=3,lwd=1,col="red")
#lines(rLandf["95%",,ac(yrs),2,,]~yrs,lty=3,lwd=1,col="red")
#legend("bottomleft",legend=c("Fleet A","Fleet B"),col=c("black","red"),lwd=2,lty=1,box.lty=0)
#text(x=xrange[1],y=yrange[2],pos=1,labels="(A)",font=fonts,cex=cl)
#
#  #------------------
#  # True F
#  #------------------
#
#xrange  <- range(yrs)
#yrange  <- yrangeF
#plot(rFp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
#     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
#axis(2,las=1,at=pretty(yrange),labels=pretty(yrange),cex=ca,font=fonts)
#mtext(text=expression(paste(F[2-6]," (",year^-1,")",sep="")),side=2,at=(yrange[2]-yrange[1])/2+yrange[1],outer=F,cex=cl,line=4,font=fonts)
#grid(); box()
#lines(rFp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
#lines(rFp["5%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
#lines(rFp["95%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
#
#  #-Reference level Fpa
#abline(h=0.25,col="darkgreen",lwd=1,lty=2);
#mtext(text="Fpa",side=4,at=0.25,las=1,cex=0.65,col="darkgreen",line=0.5,font=fonts)
#text(x=xrange[1],y=yrange[2],pos=1,labels="(B)",font=fonts,cex=cl)
#
#  #------------------
#  # True SSB
#  #------------------
#
#xrange  <- range(yrs)
#yrange  <- yrangeSSB
#plot(rSSBp["50%",,ac(yrs),,,]~yrs,xlim=xrange,ylim=yrange,xlab="",ylab="",
#     type="l",lwd=1,cex.lab=cl,cex.axis=ca,font=fonts,yaxs="i",yaxt="n",xaxt="n")
#axis(2,las=1,at=pretty(yrange),labels=pretty(yrange)/1000,cex=ca,font=fonts)
#axis(1,las=1,at=pretty(xrange),labels=pretty(xrange),cex=ca,font=fonts)
#mtext(text=expression(paste("SSB (",10^3," tonnes)",sep="")),side=2,at=(yrange[2]-yrange[1])/2+yrange[1],outer=F,cex=cl,line=4,font=fonts)
#grid(); box()
#lines(rSSBp["50%",,ac(yrs),,,]~yrs,lty=1,lwd=2)
#lines(rSSBp["5%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
#lines(rSSBp["95%",,ac(yrs),,,]~yrs,lty=3,lwd=1)
#  #-Reference level Blim & Bpa
#abline(h=0.8e6,col="blue",lwd=1,lty=2);
#abline(h=1.3e6,col="darkgreen",lwd=1,lty=2);
#mtext(text="Blim",side=4,at=0.8e6,las=1,cex=0.65,col="blue",line=0.5,font=fonts)
#mtext(text="Bpa",side=4,at=1.3e6,las=1,cex=0.65,col="darkgreen",line=0.5,font=fonts)
#text(x=xrange[1],y=yrange[2],pos=1,labels="(C)",font=fonts,cex=cl)
#
#
#
#  #- Labels x-axis
#mtext(text=expression(Years),side=1,at=(xrange[2]-xrange[1])/2+xrange[1],outer=F,cex=cl,line=4,font=fonts)
#savePlot(paste(PathPlots,scen,"Truth.png",sep=""),type="png")
#