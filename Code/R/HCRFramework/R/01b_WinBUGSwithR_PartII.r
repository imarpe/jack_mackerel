###################################################################
#DM PART 1: Does plots of Winbugs output, works out model probs, creates .dat files of SR parameters for PART 2, etc...
###################################################################
#===========----------------- load libraries -----------------------===========#
# Tested on R 2.8.1
library(coda) # library just used for reading routines for winbugs
library(MASS)

#===========--------------------- setup ----------------------------===========#
# only changes necessary to run the code can be done here

# Species
#DM - Make PLE/SOL compatible - i.e. choose species at top and it runs for both
spp <- "JM"
srComb <- "Three" # "Three", "HS_RK" or "RK_BH"
trnc <- F
                                                                              
#set of results : 
if (RecRegime=="STRec") resset  <-  "short_N"  
if (RecRegime=="LTRec") resset  <-  "all_N"       

# Directory     

#data.dir   <- paste(my.dir,"DATA/",resset,"/",sep="")
#DM Create plotting directory
pathFig <- paste(outPath,RecRegime,"/Plots input/",sep="")
#shell(paste("mkdir ",sep=""))
#pathFig <- paste(pathFig,sep="")                                                                                                                                                     
# Set working directory to data path
#setwd(data.dir)

# Save files
savePlots <- T           # save files or not
#fileFormat <- "png"      # "eps", "wmf", "jpg"
# If not plotting to file, record graph history

## Plot size
wth <- 7                                                                              
hght <- 7
# 
### truncation values from separate analysis            
# (Half increment further out from max and min observed)
#DM hardwired values - different for sole and plaice, depends on SR data. 
# Check "NS PLE and SOL SRs.xls" (in DATA directory) for values (re-calculate in this sheet if SR data changes).
#dats<-read.csv("data.csv")
#  
if (spp=="JM") {
  #       just tryied values a bit larger than the observed rct range / 100 000
  Sru <- max(win.data$Recobs)*1.1
  Srl <- 2
  }   
  

#===========------------------- load data --------------------------===========#
### Read in data files from WINBUGS ---- (2000 values per chain now, 40000 later for plotting)
#File names
Fname=c("BHa","BHb","BinfHS","BinfR","RecHS","RecR","sigmaBH","sigmaHS","sigmaR")

selection<-sample(c(1:length(BevHolt$sims.list$BHa)),2000,replace=F)
BHa<-BevHolt$sims.list$BHa[selection]
BHb<-BevHolt$sims.list$BHb[selection]
# Save orig.
BHasave=BHa; BHbsave=BHb
#Transform to the model spec used by FLR   
#DM different BH formulaiton was used in the Winbugs, tranform back to FLR version
BHb=BHasave/BHbsave
BHa=1/BHbsave
sigmaBH<- (1/BevHolt$sims.list$tauBH[selection])^0.5


#
HSa <-  HockeyStick$sims.list$RecHS[selection] 
HSb <-  HockeyStick$sims.list$BinfHS[selection] 
sigmaHS<- (1/HockeyStick$sims.list$tauHS[selection])^0.5

#
Ra <-  Ricker$sims.list$RecR[selection] 
Rb <-  Ricker$sims.list$BinfR[selection] 
sigmaR<- (1/Ricker$sims.list$tauR[selection])^0.5




BHmods<-data.frame(BHa,BHb,sigmaBH)
Rikmods<-data.frame(Ra,Rb,sigmaR)
HSmods<-data.frame(HSa,HSb,sigmaHS)

###Stock recruit data       
# Must be same data used for Winbugs!!!
# Rec is the (Rec/ 1000), SSB is SSB/1000  .

Rec <- win.data$Recobs
SSB <- win.data$SSB

SSBO=SSB
RecO=Rec
Sssb=sort(SSB,index.return=TRUE)
SSB=Sssb$x
Rec=Rec[Sssb$ix]

#===========------------------ likelihoods -------------------------===========#
#### sections below calculate log likelihood of observations given set of models
### equations will depend on type of model
## here LL is LL fit to log of REC using Normal distribution formulated as 
#DM key part here - getting the probabilities for each model type (and using for proportions later on)

MCMCN=length(Ra)
SRpair=length(SSB)
Rsym=array(0,c(MCMCN,SRpair))
SSBp=seq(0,max(SSB),max(SSB)/50)

## Hockey Stick
PHS=0
P1=c(1:MCMCN)
for (i1 in 1:MCMCN) {
sll=0
sllq=0
for (i in 1:SRpair) {
  mu=(HSa[i1]*HSb[i1]*(SSB[i]>=HSb[i1])+HSa[i1]*SSB[i]*(SSB[i]<HSb[i1])) # Rec
  Rsym[i1,i]=mu
  ll=log(1/sqrt(2*3.14159*sigmaHS[i1]^2)*exp(-1/(2*(sigmaHS[i1]^2))*(mu-Rec[i])^2)) # LL of normal dist
  sll=sll+ll
}
P1[i1]=exp(sll)
}
PHS=sum(1/P1)/MCMCN
PHS<-median(1/P1) 
PHS=1/PHS

# Plot SR data and MCMC fits    #DM units?
if (savePlots) png(file=paste(pathFig,"MCMC_SR-HockeyStick.png",sep=""),height=hght,width=wth,units="in",res=200)
plot(SSB,(Rec), xlim=c(0,1.1*max(SSB)),xlab="SSB", ylab="Recruitment")
for (i1 in seq(1,MCMCN,40)) lines(SSB,(Rsym[i1,]),col=5)
for (qnt in c(0.05,0.25,0.5,0.75,0.95)) lines(SSB,(apply(Rsym,2,quantile,qnt)),col=2)
lines(SSB,(Rsym[which.max(P1),]),col=1,lwd=2)
if (savePlots) dev.off()

## Ricker
PHR=0
P3=c(1:MCMCN)                                                                                                                                                          
for (i1 in 1:MCMCN) {                                                                                                                                            
sll=0                                                                                                                                                           
for (i in 1:SRpair) {                                                                                                                                               
  mu=(Ra[i1]*SSB[i]*exp(-Rb[i1]*SSB[i]))
  Rsym[i1,i]=mu
  ll=log(1/sqrt(2*3.14159*sigmaR[i1]^2)*exp(-1/(2*(sigmaR[i1]^2))*(mu-Rec[i])^2))
  sll=sll+ll                                                                                                                                                    
}                                                                                                                                                               
P3[i1]=exp(sll)                                                                                                                                                      
#PHR=PHR+1/P                                                                                                                                                     
}                                                                                                                                                               
PHR=sum(1/P3)/MCMCN  
PHR<-median(1/P3)                                                                                                                                             
PHR=1/PHR                                                                                                                                                       

# Plot SR data and MCMC fits    #DM units?
if (savePlots) png(file=paste(pathFig,"MCMC_SR-Ricker.png",sep=""),height=hght,width=wth,units="in",res=200)
plot(SSB,(Rec), xlim=c(0,1.1*max(SSB)),xlab="SSB", ylab="Recruitment")
for (i1 in seq(1,MCMCN,40)) lines(SSB,(Rsym[i1,]),col=5)
for (qnt in c(0.05,0.25,0.5,0.75,0.95)) lines(SSB,(apply(Rsym,2,quantile,qnt)),col=2)
lines(SSB,(Rsym[which.max(P1),]),col=1,lwd=2)
if (savePlots) dev.off()

##Bev Holt
PBH=0
P7=c(1:MCMCN)                                                                                                                                                          
for (i1 in 1:MCMCN) {                                                                                                                                            
sll=0                                                                                                                                                           
for (i in 1:SRpair) {
  mu <-(BHa[i1]*SSB[i]/(BHb[i1]+SSB[i]))
  Rsym[i1,i]=mu
  ll=log(1/sqrt(2*3.14159*sigmaBH[i1]^2)*exp(-1/(2*(sigmaBH[i1]^2))*(mu-Rec[i])^2))
  sll=sll+ll                                                                                                                                                    
#print(c(x,log(Rec[i])))
}                                                                                                                                                               
P7[i1]=exp(sll)                                                                                                                                                      
                                                                                                                          
}                                                                                                                                                               
PBH=sum(1/P7)/MCMCN   
PBH<-median(1/P7)                                                                                                                                               
PBH=1/PBH                                                                                                                                                       

# Plot SR data and MCMC fits    #DM units?
if (savePlots) png(file=paste(pathFig,"MCMC_SR-BevHolt.png",sep=""),height=hght,width=wth,units="in",res=200)
plot(SSB,(Rec), xlim=c(0,1.1*max(SSB)),xlab="SSB", ylab="Recruitment")
for (i1 in seq(1,MCMCN,40)) lines(SSB,(Rsym[i1,]),col=5)
for (qnt in c(0.05,0.25,0.5,0.75,0.95)) lines(SSB,(apply(Rsym,2,quantile,qnt)),col=2)
lines(SSB,(Rsym[which.max(P1),]),col=1,lwd=2)
if (savePlots) dev.off()

##DM prints out probabilities of each of the three models picked up   (1,3,7 - random numbers)
pp=seq(1,MCMCN,1)
P1a=sort(P1)[pp]
P3a=sort(P3)[pp]
P7a=sort(P7)[pp]

Total=1
for (nn in seq(1,500,2)) {
pp=seq(nn,MCMCN,1)
P1a=sort(P1/Total)[pp]
P3a=sort(P3/Total)[pp]
P7a=sort(P7/Total)[pp]
P1t=1/mean(1/P1a)
P3t=1/mean(1/P3a)
P7t=1/mean(1/P7a)
}

pp=seq(1,MCMCN,1)
P1a=sort(P1/Total)[pp]
P3a=sort(P3/Total)[pp]
P7a=sort(P7/Total)[pp]

# Plot model likelihoods
if (savePlots) png(file=paste(pathFig,"_Bayes-Lkhd.png",sep=""),height=hght,width=wth,units="in",res=200)
plot(x=pp,y=log(P1a),type="n",main="Likelihood of Bayes model",xlab="model order",ylab="log likelihood") #SOL ylim=c(-75,-55)
lines(x=pp,y=log(P1a),lty=1,col=1)
lines(x=pp,y=log(P3a),lty=2,col=2)
lines(x=pp,y=log(P7a),lty=3,col=3)
legend(x="bottomright",legend=c("HS","Rk","BH"),lty=c(1,2,3),col=c(1,2,3))
if (savePlots) dev.off()

##DM getting coefficients for max likelihood model fit
pp1=which.max(P1)             
pp3=which.max(P3)
pp7=which.max(P7)
maxLFit <- paste('H-Stick : A=',HSa[pp1]," B=",HSb[pp1]," sigma=",sigmaHS[pp1],
';Ricker  : A=',Ra[pp3]," B=",Rb[pp3]," sigma=",sigmaR[pp3],
';Bev-Holt: A=',BHa[pp7]," B=",BHb[pp7]," sigma=",sigmaBH[pp7], sep=" ") 

write(maxLFit, file=paste(pathFig,spp,"_MaxLikelihoodFits.txt",sep=""),sep=",")


##
#===========--------------- Model proportions ----------------------===========#
##################### VALUES SET FROM RESULTYS ABOVE change fractions to fit probabilities
# change start values to get correct total proportions
srComb<-"Three"
if (srComb=="Three") {
  # three models
  srProbs <- c(PHS,PHR,PBH)/sum(PHS,PHR,PBH)
  probHS <- srProbs[1]
  probHR <- srProbs[2]
  probBH <- srProbs[3]
  # HS RK BH = 0.1421944 0.8578056 0.0000000         #DM these values are hard wired
  mod=array("HSL",1000)
  for (i in round(seq(2,1000,(probHS+probHR)/probHR))){mod[i]="RKL"}
  for (i in round(seq(3.4,1000,1/probBH))){mod[i]="BHL"}
  #DM check to see if ended up with right numbers
  print(c(HockeyStick=sum(mod=="HSL"),Ricker=sum(mod=="RKL"),BevHolt=sum(mod=="BHL")))           
  
  } else if (srComb=="HS_RK") {
  # two models: HS and Ricker
  srProbs <- c(PHS,PHR)/sum(PHS,PHR)
  probHS <- srProbs[1]
  probHR <- srProbs[2]
  # HS RK = 0.5832516	0.4167484
  mod=array("HSL",1000)
  for (i in round(seq(1.5,1000,1/probHR))){mod[i]="RKL"}
  #DM check to see if ended up with right numbers
  print(c(sum(mod=="HSL"),sum(mod=="RKL")))           
  
  } else if (srComb=="RK_BH") {
  # two models: Ricker and BevHolt
  srProbs <- c(PHR,PBH)/sum(PHR,PBH)
  probHR <- srProbs[1]
  probBH <- srProbs[2]
  # HS RK = 0.5832516	0.4167484
  mod=array("RKL",1000)
  for (i in round(seq(1.5,1000,1/probBH))){mod[i]="BHL"}
  #DM check to see if ended up with right numbers
  print(c(sum(mod=="RKL"),sum(mod=="BHL")))           
  }
  
#===========-------------------- Merging ---------------------------===========#
###################################################
# setion to put together models to be used
# Unthinned values used, so contours look good
#
## Loop over SR types
#for (sr in c("HS","RH","BH")) {
#
##First save max LL parameter values from input data
#if (sr=="HS") {
## log lilihood vbalues for HS 
#LLa = HSa[pp1]
#LLb = HSb[pp1]
#LLsigma=sigmaHS[pp1]
#i=3
#HSb=read.coda(paste(data.dir,Fname[i],"all1.txt",sep=''),paste(data.dir,Fname[i],"allIndex.txt",sep=''))
#i=5
#HSa=read.coda(paste(data.dir,Fname[i],"all1.txt",sep=''),paste(data.dir,Fname[i],"allIndex.txt",sep=''))
#i=8
#sigmaHS=read.coda(paste(data.dir,Fname[i],"all1.txt",sep=''),paste(data.dir,Fname[i],"allIndex.txt",sep=''))
#AA=HSa
#BB=HSb
############################ jump down
#} else if (sr=="RH") {
## log lilihood vbalues for Ricker 
#LLa = Ra[pp3]
#LLb = Rb[pp3]
#LLsigma=sigmaR[pp3]
#i=4
#Rb=read.coda(paste(data.dir,Fname[i],"all1.txt",sep=''),paste(data.dir,Fname[i],"allIndex.txt",sep=''))
#i=6
#Ra=read.coda(paste(data.dir,Fname[i],"all1.txt",sep=''),paste(data.dir,Fname[i],"allIndex.txt",sep=''))
#i=9
#sigmaR=read.coda(paste(data.dir,Fname[i],"all1.txt",sep=''),paste(data.dir,Fname[i],"allIndex.txt",sep=''))
#AA=Ra
#BB=Rb
####################### jump down
#} else if (sr=="BH") {
## log lilihood vbalues for BH
#LLa = BHasave[pp7]
#LLb = BHbsave[pp7]
#LLsigma=sigmaBH[pp7]
#i=1
#BHa=read.coda(paste(data.dir,Fname[i],"all1.txt",sep=''),paste(data.dir,Fname[i],"allIndex.txt",sep=''))
#i=2
#BHb=read.coda(paste(data.dir,Fname[i],"all1.txt",sep=''),paste(data.dir,Fname[i],"allIndex.txt",sep=''))
#i=7
#sigmaBH=read.coda(paste(data.dir,Fname[i],"all1.txt",sep=''),paste(data.dir,Fname[i],"allIndex.txt",sep=''))
#AA=BHa
#BB=BHb
#}
#
####l run this plotiing bit for each bayes model ------------ Routine borrowed and amended from Mark Payne
#AA.est=median(AA)
#BB.est=median(BB)
#
## some settings normally set in the subroutine call
#n.grid=50    # this and the n below will control smoothness
#show.points=TRUE
#show.ll=TRUE
#do.contours=TRUE
#filled.contours=FALSE
#f.ages=NULL
#margin.plots=TRUE
#xlim= max(BB)
#ylim= max(AA)
#debug=FALSE
#n=40000 # set to match my full data length
#pch="."
#show.grid=TRUE
#alpha=0.05   # sets intervals on pfs
#show.estimate=TRUE
#thin=20
#TH=seq(thin/2,n,thin)
#
## main contouring section
#
#    kern  <-  kde2d(BB,AA,n=n.grid)
#    #Calculate cumulative distribution function
#    kz    <-  as.vector(kern$z)
#    ord   <-  order(kz)
#    cumfrac <-  cumsum(kz[ord])/sum(kz)
#    cumfrac.matrix  <-  matrix(cumfrac[rank(kz)],nrow=nrow(kern$z),ncol=ncol(kern$z))
#    contour.args <- list()
#    contour.args$levels <-   c(0.1,0.25,0.50,0.75,0.90)
#    contour.args$labels <-   NULL
##      if(is.null(contour.args$lty))   contour.args$lty <-   c(1,1,2,3,4)
##      if(is.null(contour.args$lwd))   contour.args$lwd <-   c(1,3,1,1,1)
#    if(is.null(contour.args$method))contour.args$method <-    "edge"
##    if(is.null(contour.args$labels))contour.args$labels <-    NULL
##      if(filled.contours) {      
##       do.call(filled.contour,c(x=list(kern$x),y=list(kern$y),z=list(cumfrac.matrix),add=TRUE,nlevels=100,color.palette=heat.colors))
##    }
#    otolith.obj  <-  c(x=list(kern$x),y=list(kern$y),z=list(cumfrac.matrix),add=TRUE,contour.args)
## this form works but gives odd contour labeling  - alternative crude seting at end to deal with issue
## I cannot work out how to change contour.args the iff line 6 above commented out was an attemp that did not work 
#
#if (savePlots) png(file=paste(pathFig,sr,"_Otolith.png",sep=""),height=hght,width=wth,units="in",res=200)
#
#    if(!show.points) { pch <- NA }
#    if(margin.plots) {
#      layout(matrix(c(1,4,3,2),2,2,byrow=TRUE), c(3,1), c(1,3), TRUE)
#      par(mar=c(0.5,0.5,0.5,0.5),oma=c(5.1,4.1,4.1,2.1))
#    }
#
#
##    x.lab <-  "Q" # not used I think
#   xlim  <- range(pretty(BB))
#   ylim  <- range(pretty(AA))
#    #First the horizontal plot
#    if(margin.plots) { 
#      densF   <-  density(BB)
#      plot(densF,ann=FALSE,xaxt="n",yaxt="n",type="l",xlim=xlim)   
#      if(show.grid) grid()      
#      title(ylab="Prob. Density",xpd=NA,mgp=c(1,1,0))
#      #Calculate 95% confidence intervals
#      cumsumF.fun <-  approxfun(cumsum(densF$y)/sum(densF$y),densF$x)
#      densF.fun   <-  approxfun(densF$x,densF$y)
#      ul.F    <-  cumsumF.fun(1-alpha/2)      
#      ul.dens <-  densF.fun(ul.F)
#      ll.F    <-  cumsumF.fun(alpha/2)      
#      ll.dens <-  densF.fun(ll.F)
#      points(c(ll.F,ul.F),c(ll.dens,ul.dens),pch="|",cex=1.5)
#      text(c(ll.F,ul.F),c(ll.dens,ul.dens),label=sprintf("%.3f",round(c(ll.F,ul.F),3)),pos=4)
#      if(show.estimate) { 
#        points(BB.est,densF.fun(BB.est),pch=19,cex=1.5)
#        text(BB.est,densF.fun(BB.est),label=sprintf("%.3f",round(BB.est,3)),pos=4)
#        }
#      if(show.ll) { 
#        points(LLb,densF.fun(LLb),pch=19,cex=1.5,col=4)
#        text(LLb,densF.fun(LLb),label=sprintf("%.3f",round(LLb,3)),pos=4)
#      }
#    }
#    #Now the vertical plot
#    if(margin.plots) { 
#      densAA <-  density(AA)
#      plot(densAA$y,densAA$x,xaxt="n",yaxt="n",type="l",ylim=ylim)
#      abline(v=0,col="grey")
#      if(show.grid) grid()
#      title(xlab="Prob. Density",xpd=NA,mgp=c(1,1,0))      
#      #Calculate 95% confidence intervals
#      cumsumAA.fun <-  approxfun(cumsum(densAA$y)/sum(densAA$y),densAA$x)
#      densAA.fun   <-  approxfun(densAA$x,densAA$y)
#      ul.AA    <-  cumsumAA.fun(1-alpha/2)      
#      ul.dens <-  densAA.fun(ul.AA)
#      ll.AA    <-  cumsumAA.fun(alpha/2)      
#      ll.dens <-  densAA.fun(ll.AA)
#      points(c(ll.dens,ul.dens),c(ll.AA,ul.AA),pch="-",cex=2)
#      text(c(ll.dens,ul.dens),c(ll.AA,ul.AA),label=round(c(ll.AA,ul.AA),3),pos=4)
#      if(show.estimate) {
#        points(densAA.fun(AA.est),AA.est,pch=19,cex=1.5)      
#        text(densAA.fun(AA.est),AA.est,,label=round(AA.est,3),pos=2)
#        }
#      if(show.ll) { 
#        points(densAA.fun(LLa),LLa,pch=19,cex=1.5,col=4)      
#        text(densAA.fun(LLa),LLa,,label=round(LLa,3),pos=2)
#      }
#    }
#    #Now the main plot
#    plot(0,0,xlim=xlim,ylim=ylim,type="n",xlab="",ylab="")
##   Three sets of points for 3 Bayes MCMC chains - normally one would do
#    if(show.points) points(BB[TH],AA[TH],pch=19,col=2,cex=0.05)
##    if(show.points) points(QC2[thin],QM2[thin]*0.15,pch=19,col=3,cex=0.05)
##    if(show.points) points(QC3[thin],QM3[thin]*0.15,pch=19,col=4,cex=0.05)
#    title(xlab="B",ylab="A",xpd=NA)
##   Crude way to get second graph
#    #if(show.points) points(QC1[thin],MZ1[thin],pch=19,col=2,cex=0.05)
#    #if(show.points) points(QC2[thin],MZ2[thin],pch=19,col=3,cex=0.05)
#    #if(show.points) points(QC3[thin],MZ3[thin],pch=19,col=4,cex=0.05)
#    #title(xlab="Q",ylab="Z",xpd=NA)
#    if(show.estimate) points(BB.est,AA.est,pch=19,cex=1.5)
#    if(show.ll) points(LLb,LLa,pch=19,col=4,cex=1.5)
#    if(show.grid) grid()
#    contour(otolith.obj,levels=contour.args$levels,lables=NULL,add=TRUE)
##    lines below were origonal ---- I uses line above to get correct labled contours
## dont run lines below
#    if(do.contours) {
#      do.call(contour,otolith.obj)
#    }
#    
#if (savePlots) dev.off()
#    
###########################################################################################
#  } #end of sr loop
########### end of plotting ---- use above one for each model type
#

if(1==2)
{

#DM slotting routines into data frame with the right variables for all, then moving across
P1=c(1:2000)
Sr=array(0,c(SRpair,3000))
dev=array(0,c(SRpair,3000))
HSNmu=array(0,c(SRpair,2000))
lowlim=array(0,2000)
uplim=array(0,2000)
RecHSNub=array(0,2000)
mincnt=array(0,2000)
maxcnt=array(0,2000)
RecHSNmn=array(0,2000)
for (i1 in 1:2000) {
sll=0
for (i in 1:SRpair) {
  mu=(HSa[i1]*HSb[i1]*(SSB[i]>=HSb[i1])+HSa[i1]*SSB[i]*(SSB[i]<HSb[i1]))
  dev[i,]=rnorm(3000,0,sigmaHS[i1])
  Sr[i,]=(mu+dev[i,])
  HSNmu[i,i1]=mean(Sr[i,])
}
mincnt[i1]=sum(Sr<Srl)
maxcnt[i1]=sum(Sr>Sru)
lowlim[i1] = min(dev[(Sr>Srl)])
uplim[i1] = min(dev[(Sr>Sru)])    #### why is it min and not max?

RecHSNub[i1]=HSa[i1]+log(mean(HSNmu[,i1])/mean(Sr*(Sr>Srl)*(Sr<Sru)))
RecHSNmn[i1]=mean(HSNmu[,i1])
#print(c(mincnt[i1],maxcnt[i1],lowlim[i1],uplim[i1], mean(HSNmu[,i1]), HSa[i1],RecHSNub[i1]))
P1[i1]=exp(sll)
#PHS=PHS+1/P
}
PHS=median(1/P1)
PHS=1/PHS
mincnt=mincnt/SRpair/3000
maxcnt=maxcnt/SRpair/3000
HSN=as.data.frame(cbind(HSa,HSb,sigmaHS,mincnt,maxcnt,lowlim,uplim,RecHSNmn,RecHSNub))



PHR=0
P2=c(1:2000)                                                                                                                                                          
RNmu=array(0,c(SRpair,2000))
RecRNub=array(0,2000)
RecRNmn=array(0,2000)
for (i1 in 1:2000) {
sll=0                                                                                                                                                           
for (i in 1:SRpair) {
  mu=log(Ra[i1]*SSB[i]*exp(-Rb[i1]*SSB[i]))
  dev[i,]=rnorm(3000,0,sigmaR[i1])
  Sr[i,]=(mu+dev[i,])
  RNmu[i,i1]=mean(Sr[i,])
}
mincnt[i1]=sum(Sr<Srl)
maxcnt[i1]=sum(Sr>Sru)
lowlim[i1] = min(dev[(Sr>Srl)])
uplim[i1] = min(dev[(Sr>Sru)])

RecRNub[i1]=Ra[i1]+log(mean(RNmu[,i1])/mean(Sr*(Sr>Srl)*(Sr<Sru)))
RecRNmn[i1]=mean(RNmu[,i1])
#print(c(mincnt[i1],maxcnt[i1],lowlim[i1],uplim[i1], mean(RNmu[,i1]), Ra[i1],RecRNub[i1]))
P2[i1]=exp(sll)
#PHR=PHR+1/P                                                                                                                                                     
}                                                                                                                                                               
PHR=median(1/P2)                                                                                                                                          
PHR=1/PHR                                                                                                                                                       
mincnt=mincnt/3000/SRpair
maxcnt=maxcnt/3000/SRpair
RkN=as.data.frame(cbind(Ra,Rb,sigmaR,mincnt,maxcnt,lowlim,uplim,RecRNmn,RecRNub))


BHN=0
P4=c(1:2000)                                                                                                                                                          
BHNmu=array(0,c(SRpair,2000))
RecBHNub=array(0,2000)
RecBHNmn=array(0,2000)
for (i1 in 1:2000) {
for (i in 1:SRpair) {
  mu <-(BHa[i1]*SSB[i]/(BHb[i1]+SSB[i]))
  dev[i,]=rnorm(3000,0,sigmaBH[i1])
  Sr[i,]=(mu+dev[i,])
  BHNmu[i,i1]=mean(Sr[i,])
#print(c(x,log(Rec[i])))
}                                                                                                                                                               
mincnt[i1]=sum(Sr<Srl)
maxcnt[i1]=sum(Sr>Sru)
lowlim[i1] = min(dev[(Sr>Srl)])
uplim[i1] = min(dev[(Sr>Sru)])

RecBHNub[i1]=BHa[i1]+log(mean(BHNmu[,i1])/mean(Sr*(Sr>Srl)*(Sr<Sru)))
RecBHNmn[i1]=mean(BHNmu[,i1])
#print(c(mincnt[i1],maxcnt[i1],lowlim[i1],uplim[i1], mean(BHNmu[,i1]),BHa[i1],RecBHNub[i1]))

}                                                                                                                                                               
mincnt=mincnt/3000/SRpair
maxcnt=maxcnt/3000/SRpair
BHN=as.data.frame(cbind(BHa,BHb,sigmaBH,mincnt,maxcnt,lowlim,uplim,RecBHNmn,RecBHNub))

#===========--------------- Create Model Sets ----------------------===========#
####################################################
#This is set up based on results above
#So with different factors it would be different

print(c(Ricker=sum(mod=="RKL"),HockeyStick=sum(mod=="HSL"),BevHolt=sum(mod=="BHL")))
  }
#########################################################
# withOUT tuncation

pp=c(1:1000)
HSLp=2*pp[mod=="HSL"]
RKLp=2*pp[mod=="RKL"]
BHLp=2*pp[mod=="BHL"]

untrncModset=as.data.frame(cbind(mod,BHmods[1:1000,]))

untrncModset[,-1]<-NA

i=RKLp
untrncModset[i/2,-1]<-Rikmods[i,]


i=BHLp
untrncModset[i/2,-1]=BHmods[i,]

i=HSLp
untrncModset[i/2,-1]=HSmods[i,]

names(untrncModset)[2:4]=c("A","B","sigma")


 #end of not truncated loop
#
##########################################################
## with tuncation
#
#pp=c(1:1000)
#HSLp=pp[mod=="HSL"]
#RKLp=pp[mod=="RKL"]
#BHLp=pp[mod=="BHL"]
#
#trncModset=as.data.frame(cbind(mod,BHN[1:1000,1],BHN[1:1000,2:3],BHN[1:1000,6:7],BHN[1:1000,4:5],BHN[1:1000,9]/BHN[1:1000,1]))
#
#i=RKLp
#trncModset[i,2]=RkN[i,9]
#trncModset[i,3:4]=RkN[i,2:3]
#trncModset[i,5:6]=RkN[i,6:7]
#trncModset[i,7:8]=RkN[i,4:5]
#trncModset[i,9]=RkN[i,9]/RkN[i,1]
#
#i=BHLp
#trncModset[i,2]=BHN[i,9]
#trncModset[i,3:4]=BHN[i,2:3]
#trncModset[i,5:6]=BHN[i,6:7]
#trncModset[i,7:8]=BHN[i,4:5]
#trncModset[i,9]=BHN[i,9]/BHN[i,1]
##
#i=HSLp
#trncModset[i,2]=HSN[i,9]
#trncModset[i,3:4]=HSN[i,2:3]
#trncModset[i,5:6]=HSN[i,6:7]
#trncModset[i,7:8]=HSN[i,4:5]
#trncModset[i,9]=HSN[i,9]/HSN[i,1]
#
#
#names(trncModset)[2:4]=c("A","B","sigma")
#names(trncModset)[9]="Cfac"
#
#if (savePlots) png(file=paste(pathFig,"_ModsUpperlim_trunc.png",sep=""),units="in",res=200,height=hght,width=wth)
#hist(trncModset[,8]*100,30,main="No of Models exceeding upper limit", xlab="Percentage of simulated values exceeding upper limit" )
#if (savePlots) dev.off()
#
#if (savePlots) png(file=paste(pathFig,"_ModsLowerlim_trunc.png",sep=""),units="in",res=200,height=hght,width=wth)
#hist(trncModset[,7]*100,30,main="No of Models exceeding lower limit", xlab="Percentage of simulated values exceeding lower limit" )
#if (savePlots) dev.off()
#
 #end of not truncated loop
  #    
##########################################################
##DM saving the model parameters
#write.table(trncModset[,1:6],file=paste("_modelparams_trnc.dat",sep=""),row.names=FALSE,col.names=TRUE,sep=" ")
#write.table(untrncModset,file=paste("_modelparams.dat",sep=""),row.names=FALSE,col.names=TRUE,sep=" ")
#write.table(trncModset,file=paste("_modelparamsfull_trnc.dat",sep="_"),row.names=FALSE,col.names=TRUE,sep=" ")
#write.csv(untrncModset,file=paste("_modelparamsfull.csv",sep="_"))
#
#
#===========--------------- Plot Merged SR fit ---------------------===========#
#DM proceed with plotting according to truncation setting


Modset <- untrncModset


#DM plotting merged SR fit and point scatters. (with and without trunation)
# plot out recruits
mn=40
SSBs=seq(1:mn*200)
Recs=seq(1:mn*200)

### Untruncated
if (!trnc) {
if (savePlots) png(file=paste(pathFig,"_MergedSRfit_unTrunc.png",sep=""),units="in",res=200,height=hght,width=wth)
Xmax <- ceiling(max(SSBO, na.rm=T))
Ymax <- ceiling(max((RecO), na.rm=T))
plot(SSBO,(RecO),xlim=c(0,Xmax),ylim=c(-10,Ymax),type="p",pch=19,col=10,xlab="SSB (1000 of tonnes)",ylab="Recruits (*1E+?)", main=spp)
abline(h=0)
for (i in 1:200) {
SSB = seq(0,7,length.out=mn)
#SSB = seq(0,max(SSBO),length.out=length(SSBO)*40) ## to match SSB values to original
 
 
 if (Modset[i,1]=="HSL") {
  mu=(Modset$A[i]*Modset$B[i]*(SSB>=Modset$B[i])+Modset$A[i]*SSB*(SSB<Modset$B[i]))
  R2=rnorm(mn,0,Modset$sigma[i])
  Rec=mu+R2[1:mn]
 }
 if (Modset[i,1]=="RKL") {
 mu=(Modset$A[i]*SSB*exp(-Modset$B[i]*SSB))
  R2=rnorm(mn,0,Modset$sigma[i])
  Rec=mu+R2[1:mn]
 }
 if (Modset[i,1]=="BHL") {
  mu <-(Modset$A[i]*SSB/(Modset$B[i]+SSB))
  R2=rnorm(mn,0,Modset$sigma[i])
  Rec=mu+R2[1:mn]
 }
points(SSB[1:(mn)],Rec[1:(mn)],type="p",pch=20,col=1,cex=0.0625)
SSBs[((i-1)*mn+1):(i*mn)]=SSB
Recs[((i-1)*mn+1):(i*mn)]=Rec
}
Rsims<-Recs
points(SSBO,(RecO),type="p",pch=19,col=10,cex=1.25)


seqMin <- 0
seqStep <- 0.5

up=seq(seqMin,Xmax,seqStep)
lw=seq(seqMin,Xmax,seqStep)
md=seq(seqMin,Xmax,seqStep)
ssb=seq(seqMin,Xmax,seqStep)
loopNum <- length(ssb)    
for (j in 1:loopNum){     
  up[j]=quantile(Recs[((SSBs>up[j])*(SSBs<up[j+1]))>0],probs=.95,na.rm=TRUE)
  lw[j]=quantile(Recs[((SSBs>lw[j])&(SSBs<lw[j+1]))>0],probs=.05,na.rm=TRUE)
  md[j]=quantile(Recs[((SSBs>md[j])&(SSBs<md[j+1]))>0],probs=.5,na.rm=TRUE)
  }
  
lines(ssb,md[1:loopNum],col="green",lwd=3)
lines(ssb,up[1:loopNum],col=4,lwd=3)
lines(ssb,lw[1:loopNum],col=4,lwd=3)

if (savePlots) dev.off()
#DM ?
#mean(Recs)/mean(exp(RecO))

  } # end of untruncated plot


#===========------------------- Correction -------------------------===========#
## correction factor if requred
#DM work out difference between mean rec from observed SSB levels vs that generated from the models
#mean(Recs)/mean(exp(RecO))          
Correct=mean(Recs)/mean(exp(RecO))
Modset$A=Modset$A/Correct    ####### correction factor just applied on A coefficient (good enough?)
#
##DM writing out corrected values
##write.table(Modset[,1:6],file=paste(spp,"_modelparamsCor",if(trnc) "_trnc",".dat",sep=""),row.names=FALSE,col.names=TRUE,sep=" ")
#write.table(Modset,file=paste(spp,"_modelparamsCorfull",if(trnc) "_trnc",".dat",sep="_"),row.names=FALSE,col.names=TRUE,sep=" ")
#
#===========---------------------- End -----------------------------===========#
#DM end of the bayesian part (mainly plotting)
# Proceed to "proj model combSR verDM PLE.r"
#Modsetsave<-Modset
Modset <- untrncModset
Modset<-Modset[,1:4]
Modset$A[Modset$mod=="HSL"]<-Modset$A[Modset$mod=="HSL"]
Modset$B[Modset$mod=="HSL"]<-Modset$B[Modset$mod=="HSL"]*1000

Modset$A[Modset$mod=="RKL"]<-Modset$A[Modset$mod=="RKL"]/1
Modset$B[Modset$mod=="RKL"]<-Modset$B[Modset$mod=="RKL"]*1e-3

Modset$A[Modset$mod=="BHL"]<-Modset$A[Modset$mod=="BHL"]*1000
Modset$B[Modset$mod=="BHL"]<-Modset$B[Modset$mod=="BHL"]*1000


a<-aggregate(Modset$A,by=list(Modset$mod),FUN=median)
b<-aggregate(Modset$B,by=list(Modset$mod),FUN=median)

ssB<-seq(0,max(SSBO*1000),length.out=50)

plot(SSBO*1000,RecO*1000)
rk<-a[3,2]*ssB*exp(-ssB*b[3,2])
lines(ssB,rk)

hs<-(a[2,2]*ssB*(ssB<b[2,2]) ) + a[2,2]*b[2,2] *(ssB>=b[2,2])
lines(ssB,hs)

bh<- a[1,2]*ssB/(b[1,2]+ssB)
lines(ssB,bh)
write.csv(Modset,file=paste(dataPath,"SRbayeswith R_1000models_",resset,"_Norm.csv",sep=""))




