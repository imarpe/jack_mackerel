
#USER INPUT REQUIRED!
#1. Install and load the package R2WinBUGS
library(R2WinBUGS)

MyWinBugsDir <- "C:/WinBUGS14/"


# collect the SR pairs from the assessment output
Rtemp               <-  data.frame(res$R[,1:2])
SSBtemp             <-  data.frame(res$SSB[,1:2])
names(Rtemp)        <- c("year","R")
names(SSBtemp)      <- c("year","SSBtp")
SRpairs             <- merge(Rtemp,SSBtemp,all.y=T)
Rage                <- min(ages)
SRpairs$SSB         <- NA
SRpairs$SSB[-c(1:Rage)] <-  SRpairs$SSBtp[-c(dim(SRpairs)[1]:(dim(SRpairs)[1]+1-Rage))]
SRpairs             <-  SRpairs[!is.na(SRpairs$R),c("year","R","SSB")]

SRpairs             <- SRpairs[is.element(SRpairs$year,as.numeric(recrPeriod)),]
rm(Rtemp,SSBtemp)

set.seed(1)

win.data            <- list(n      = dim(SRpairs)[1],
                            Recobs = SRpairs$R/1000,
                            SSB    = SRpairs$SSB/1000)



#-------------------------------------------------------------------------------
# Beverton and Holt model
#-------------------------------------------------------------------------------

sink("WinBUGStest.txt")
cat("
  model{
    #Diffuse Prior
  BHa~dgamma(1.01,0.125)
  BHb~dgamma(1.2,0.01)
  tauBH~dgamma(0.001,0.001)
  sigmaBH<-sqrt(1/tauBH)
   for(i in 1:n) {
         Recmod[i] <-(SSB[i]/(BHa+SSB[i]*BHb))
         Recobs[i] ~ dnorm(Recmod[i],tauBH)
  }
}
",fill = TRUE)
sink()

#Set the initial values for the betas and sigma
inits <- function () {
   list(BHa=1,
BHb=5,
tauBH=.5)  }

#Parameters to estimate
params <- c("BHa","BHb","tauBH")


#MCMC settings
nc <- 3       #Number of chains
ni <- 1000     #Number of draws from posterior (for each chain)
nb <-  50     #Number of draws to discard as burn-in
nt <-   5     #Thinning rate

BevHolt <- bugs(data = win.data,
            inits = inits,
            parameters = params,
            model = "WinBUGStest.txt",
            n.thin = 5,
            n.chains = 3,
            n.burnin = 50,
            n.iter = 40000,
            n.sims=100,
            debug = FALSE,
            bugs.directory = MyWinBugsDir)


#-------------------------------------------------------------------------------
# Ricker model
#-------------------------------------------------------------------------------
sink("WinBUGStest.txt")
cat("
  model{
    #Diffuse Prior
  BinfR~dgamma(1.2,0.5)
  RecR~dgamma(1.5,.05)
  tauR~dgamma(0.001,0.001)
  sigmaR<-sqrt(1/tauR)

   for(i in 1:n) {
         Recmod[i] <-(RecR *SSB[i]*exp(-BinfR*SSB[i]))
         Recobs[i] ~ dnorm(Recmod[i],tauR)
  }
}
",fill = TRUE)
sink()

#Set the initial values for the betas and sigma
inits <- function () {
   list(BinfR=1,
RecR=5,
tauR=.5)  }

#Parameters to estimate
params <- c("BinfR","RecR","tauR")


#MCMC settings

Ricker <- bugs(data = win.data,
            inits = inits,
            parameters = params,
            model = "WinBUGStest.txt",
            n.thin = 5,
            n.chains = 3,
            n.burnin = 50,
            n.iter = 40000,
            n.sims=100,
            debug = FALSE,
            bugs.directory = MyWinBugsDir)


#-------------------------------------------------------------------------------
# Hockey stick model
#-------------------------------------------------------------------------------

sink("WinBUGStest.txt")
cat("
model{
  #Diffuse Prior
  BinfHS~dgamma(1.1,0.05)I(1.80465,4.16656)
  RecHS~dgamma(1.3,0.05)
  tauHS~dgamma(0.001,0.001)
  sigmaHS<-sqrt(1/tauHS)
   for(i in 1:n) {
            Recmod[i] <-( RecHS *step(SSB[i]-BinfHS)*BinfHS+RecHS*SSB[i]*step(BinfHS-SSB[i]))
          Recobs[i] ~ dnorm(Recmod[i],tauHS)
  }
}
",fill = TRUE)
sink()

#Set the initial values for the betas and sigma
inits <- function () {
   list(BinfHS=2,
RecHS=2,
tauHS=.5)  }

#Parameters to estimate
params <- c("BinfHS","RecHS","tauHS")


#MCMC settings
nc <- 3       #Number of chains
ni <- 1000     #Number of draws from posterior (for each chain)
nb <-  50     #Number of draws to discard as burn-in
nt <-   5     #Thinning rate

HockeyStick <- bugs(data = win.data,
            inits = inits,
            parameters = params,
            model = "WinBUGStest.txt",
            n.thin = 5,
            n.chains = 3,
            n.burnin = 50,
            n.iter = 40000,
            n.sims=100,
            debug = FALSE,
            bugs.directory = MyWinBugsDir)
            
            
            
