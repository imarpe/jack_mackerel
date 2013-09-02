#-------------------------------------------------------------------------------
# Script to run the assessment of Jack Mackerel and look at outputs
#
# By: Niels Hintzen
# Upate: 31 Aug 2011
#-------------------------------------------------------------------------------


rm(list=ls())
memory.size(4000)

  # Set libraries & source code
library(lattice)
require(PBSadmb)
library(RColorBrewer)
library(doBy)

  # Set paths
codePath    <- "N:/Projecten/SouthPacific/2012/SWG/R/"
inputPath   <- "D:/Repository/Dropbox/jurel/Codes/Work/"
outputPath  <- "N:/Projecten/SouthPacific/2012/SWG/Results/"
setwd(codePath)

  # Specify control file
controlFile <- "mod1"
  getwd()
  # Run the assessment
source("diagnostics.r")
source("ADMB2R_15102012.r")
source("compareRuns.r")
system(paste('"jjm.exe"','-ind',paste(controlFile,".ctr",sep=""),'-nox'), wait = TRUE)

  # Read in the output of the assessment
jjm.in  <- read.dat(iFilename = "comb.dat",iPath=inputPath)
jjm.in  <- read.dat(iFilename = "north1.dat",iPath=inputPath)
jjm.in  <- read.dat(iFilename = "south.dat",iPath=inputPath)
jjm.in  <- read.dat(iFilename = "mod2.dat",iPath=inputPath)
jjm.out <- readList(file.path(inputPath,"arc/mod7c_r.rep"))
jjm.out <- readList(file.path(inputPath,"arc/comb_r.rep"))
jjm.out <- readList(file.path(inputPath,"arc/north_r.rep"))
jjm.out <- readList(file.path(inputPath,"arc/south_r.rep"))

jjm.ypr <- readYPR(file.path(inputPath,"Fprof.yld"))

  # Create the diagnostics
pdf(paste(outputPath,"summary_test_together.pdf",sep=""),height=16.6,width=12.5,pointsize = 24, bg = "white")
#pdf(paste(outputPath,"mod1 - %02d.pdf"),units = "px", height=1200,width=900,pointsize = 24, bg = "white")
trellis.par.set(fontsize=list(text=24,points=20))
diagnostics(jjm.out,jjm.in,jjm.ypr,what=c("input","fit")) #,"projections","ypr"
dev.off()


  #Compare runs
source("compareRuns.r")
pdf(paste(outputPath,"Compare_1_4.pdf",sep=""),height=16.6,width=12.5,pointsize = 24, bg = "white")
pdf(paste(outputPath,"Compare_0.pdf",sep=""),height=12.5,width=16.6,pointsize = 24, bg = "white")

jjm.mod0.0  <- readList(file.path(inputPath,"arc/mod0.0_r.rep"))
jjm.mod0.1  <- readList(file.path(inputPath,"arc/mod0.1_r.rep"))
jjm.mod0.2  <- readList(file.path(inputPath,"arc/mod0.2_r.rep"))
jjm.mod0.3  <- readList(file.path(inputPath,"arc/mod0.3_r.rep"))
jjm.mod0.4  <- readList(file.path(inputPath,"arc/mod0.4_r.rep"))
jjm.mod0.5  <- readList(file.path(inputPath,"arc/mod0.5_r.rep"))

lstOuts   <- list(Model_0.0=jjm.mod0.0,Model_0.1=jjm.mod0.1,
                  Model_0.2=jjm.mod0.2,Model_0.3=jjm.mod0.3,
                  Model_0.4=jjm.mod0.4,Model_0.5=jjm.mod0.5
)
compareTime(lstOuts,"SSB",SD=T,Sum=NULL)
compareTime(lstOuts,"R",SD=T)
compareTime(lstOuts,"TotBiom",SD=T)
compareMatrix(lstOuts,"TotF",SD=F,Sum=NULL,YrInd=jjm.mod1$Yr,Apply=mean)
compareMatrix(lstOuts,"N",SD=F,YrInd=jjm.mod1$Yr,Apply=sum)
dev.off()

pdf(paste(outputPath,"Compare_1_3.pdf",sep=""),height=12.5,width=16.6,pointsize = 24, bg = "white")

jjm.mod1  <- readList(file.path(inputPath,"arc/mod1_r.rep"))
jjm.mod2  <- readList(file.path(inputPath,"arc/mod2_r.rep"))
jjm.mod3  <- readList(file.path(inputPath,"arc/mod3_r.rep"))
jjm.mod4  <- readList(file.path(inputPath,"arc/mod4_r.rep"))
lstOuts   <- list(Model_1=jjm.mod1,Model_2=jjm.mod2,Model_3=jjm.mod3,Model_4=jjm.mod4)
lstOuts   <- list(Model_1=jjm.mod1,Model_2=jjm.mod2,Model_3=jjm.mod3)

compareTime(lstOuts,"SSB",SD=T,Sum=NULL)
compareTime(lstOuts,"R",SD=T)
compareTime(lstOuts,"TotBiom",SD=T)
compareMatrix(lstOuts,"TotF",SD=F,Sum=NULL,YrInd=jjm.mod1$Yr,Apply=mean)
compareMatrix(lstOuts,"N",SD=F,YrInd=jjm.mod1$Yr,Apply=sum)
dev.off()

compareTimes(lstOuts,"SSB",SD=T,Sum=c("mod2","mod3"))

jjm.comb <- readList(file.path(inputPath,"arc/mod1_r.rep"))
jjm.north <- readList(file.path(inputPath,"arc/north_r.rep"))
jjm.south <- readList(file.path(inputPath,"arc/south_r.rep"))
lstOuts   <- list(comb=jjm.comb,north=jjm.north,south=jjm.south)
windows()
compareTime(lstOuts,"SSB",SD=T,Sum=c("north","south"))
compareTime(lstOuts,"R",SD=T,Sum=c("north","south"))
compareTime(lstOuts,"TotBiom",SD=T,Sum=c("north","south"))

compareMatrix(lstOuts,"TotF",SD=F,Sum=NULL,YrInd=jjm.comb$Yr,Apply=mean)
compareMatrix(lstOuts,"N",SD=F,Sum=c("north","south"),YrInd=jjm.comb$Yr,Apply=sum)

compareTimes(lstOuts,"Obs_catch",SD=F,Sum=c("north","south"),YrInd=jjm.comb$Yr,Apply=sum)
compareTimes(lstOuts,"Obs_catch",SD=F,Sum=NULL,YrInd=jjm.comb$Yr,Apply=sum)

  #Likelihood tables
  # Likelihood table
tab <- cbind(lstOuts[[1]]$Like_Comp_names,do.call(cbind,lapply(lstOuts,function(x){round(x[["Like_Comp"]],2)})))

write.csv(tab,file=file.path(inputPath,"LikelihoodTable2012.csv"),row.names=F)


  #Projection tables
Nfutscens <- length(grep("SSB_fut_",names(lstOuts[[1]])))
for(iScen in 1:Nfutscens){
  res <- as.data.frame(do.call(cbind,lapply(lstOuts,function(x){x[[paste("SSB_fut_",iScen,sep="")]][,2]})))
  res <- cbind(rep(paste("Scenario",iScen),nrow(res)),res)
  colnames(res)[1] <- "Scenario"
  if(iScen == 1) totres <- res
  if(iScen != 1) totres <- rbind(totres,res)
}
write.csv(totres,file=file.path(inputPath,"SSBProjection2012.csv"),row.names=F)

  #Risk tables
jjm.mod7  <- readList(file.path(inputPath,"arc/mod7_r.rep"))
jjm.mod6  <- readList(file.path(inputPath,"arc/mod6_r.rep"))
lstOuts   <- list(Model_6=jjm.mod6,Model_7=jjm.mod7)

  #Get the future SSBs and SDs together in one file
fut       <- do.call(rbind,lapply(lstOuts,function(x){
                     do.call(rbind,lapply(x[grep("SSB_fut_",names(x))],
                                          function(y){return(y[,1:3])}))}))
fut       <- as.data.frame(fut,stringsAsFactors=F)
colnames(fut) <- c("year","SSB","SD")
fut$modelscenario <- paste(rep(names(lstOuts),each=nrow(lstOuts[[1]]$SSB_fut_1) *
                                                   length(grep("SSB_fut_",names(lstOuts[[1]])))),
                           paste("Scen",
                                 rep(1:length(grep("SSB_fut_",names(lstOuts[[1]]))),each=nrow(lstOuts[[1]]$SSB_fut_1)),
                                 sep="_"),
                           sep="_")
  #Get the 2012 SSB and SDs together
ass       <- do.call(rbind,lapply(lstOuts,function(x){
                                  x$SSB[which(x$SSB[,1] == (x$SSB_fut_1[1,1]-1)),1:3]}))
ass       <- as.data.frame(ass,stringsAsFactors=F)
colnames(ass) <- c("year","SSB","SD")
ass$modelscenario <- names(lstOuts)

  #Do the risk calculation:
  # Each final year, what is the risk of SSB(final year) < SSB(ass year)
rsktable    <- matrix(NA,nrow=length(lstOuts),ncol=length(grep("SSB_fut_",names(lstOuts[[1]]))),
                      dimnames=list(names(lstOuts),1:length(grep("SSB_fut_",names(lstOuts[[1]])))))
ratiotable  <- matrix(NA,nrow=length(lstOuts),ncol=length(grep("SSB_fut_",names(lstOuts[[1]]))),
                      dimnames=list(names(lstOuts),1:length(grep("SSB_fut_",names(lstOuts[[1]])))))
for(i in names(lstOuts)){
  futdat <-subset(fut,year==max(fut$year) &
                  paste("Model_",unlist(strsplit(fut$modelscenario,"_"))[seq(2,nrow(fut)*4,4)],sep="") == i)
  assdat <- subset(ass,modelscenario == i)
  rsktable[i,] <- dnorm(assdat$SSB,futdat$SSB,futdat$SD)
  ratiotable[i,] <- round((futdat$SSB - assdat$SSB) / futdat$SSB *100,1)
}

  #Average catch over entire timeseries
meancatch <- do.call(rbind,lapply(lstOuts,function(x){unlist(lapply(x[grep("Catch_fut_",names(x))],function(y){mean(y[,2])}))}))
  