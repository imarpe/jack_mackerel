#-------------------------------------------------------------------------------
#
# Script simulates weight-at-age
#
# Created by: Thomas Brunel thomas.brunel@wur.nl
# Project: South Pacific
#
# Rev:
#
#
#  developed under  R version 3.0.1
#            using  FLCore    2.5.0
#
# DO NOT DISTRIBUTE OR COPY THIS CODE WITHOUT PRIOR AGREEMENT OF THE AUTHORS
#
#-------------------------------------------------------------------------------


#- requires R3.0.1 to run the
library(fArma)

#- Read the historical weights
w         <- data.frame(Year=rep(res[[paste("wt_fsh_",fsh,sep="")]]  [,1],length(ages)),Age=sort(rep(ages,length((res[[paste("wt_fsh_",fsh,sep="")]][,1])))),wt=c(res[[paste("wt_fsh_",fsh,sep="")]]  [,-1]))
w$wtp1    <- NA
w$YC      <- w$Year - w$Age
for(i in 1:dim(w)[1])
  try(w$wtp1[i] <- w$wt[w$YC==w$YC[i] & w$Age==w$Age[i]+1],silent=T)
w$deltaW  <- w$wtp1 - w$wt

# Prepare the output object
wt        <- array(data=NA,dim=c(length(ages),length(histMinYr:futureMaxYr),1,1,1,nits),dimnames=list(age=ages,year=histMinYr:futureMaxYr,"all","all","all",iters=1:nits))
wt[,histPeriod,,,,] <-  t(res[[paste("wt_fsh_",fsh,sep="")]]  [,-1])

#- Model weight by an Arma process
#first, age 1 
dats      <- w$wt[w$Age==1]
try       <- armaSearch(dats)
coefs     <- try$model
p         <- coefs[1]
q         <- coefs[2]
ff        <- paste("x~arma(",p,",",q,")",sep="")
ff        <- as.formula(ff)
try       <- armaFit( ff, data=dats)
if(p==0){ aR<-c(0) } else {aR<-coef(try)[1:p]}
if(q==0){ mA<-c(0) } else {mA<-coef(try)[(p+1):(p+q)]}

for (i in 1:nits){
  #- Rescale them
  x       <- armaSim(model = list(ar = aR, ma = mA), n = 300)
  x       <- x@.Data
  #- Standardise the series : same mean and var as the historical series
  x       <- (x-mean(x))/sd(x)
  x       <- sd(dats)*x+mean(dats)

  #- Remove the first simulated years until one is cloase enough to the last historical year to avoid a jump between historical and projection period
  last    <- mean(dats[(length(dats)-2):length(dats)])
  xin     <- (x>0.7*(last)) & (x<1.3*(last))
  xin     <- cumsum(xin)
  x       <- x[xin>0]
  x       <- x[1:length(projPeriod)]
  x[x<min(dats)]  <- min(dats)
  x       <- c(dats,x)
  wt[1,,,,,i]<-x
}

#- For ages >= min age :  modelling the weight increment and add it to the previous year's weight
for(age in ages[-length(ages)]){
  dats    <- w$deltaW[w$Age==age]
  dats    <- dats[!is.na(dats)]
  try     <- armaSearch(dats)
  coefs   <- try$model
  p       <- coefs[1]
  q       <- coefs[2]
  ff      <- paste("x~arma(",p,",",q,")",sep="")
  ff      <- as.formula(ff)
  try     <- armaFit( ff, data=dats)
  if(p==0){ aR<-c(0) } else {aR<-coef(try)[1:p]}
  if(q==0){ mA<-c(0) } else {mA<-coef(try)[(p+1):(p+q)]}

  for(i in 1:nits){
    #- Try to find a way to rescale them
    x     <- armaSim(model = list(ar = aR, ma = mA), n = 300)
    x     <- x@.Data
    #- standardise the series
    x     <- (x-mean(x))/sd(x)
    x     <- sd(dats)*x+mean(dats)

    last  <- mean(dats[(length(dats)-2):length(dats)])
    #- Which simulated values are within +-30% of the last observed?
    xin   <- (x>0.7*(last)) & (x<1.3*(last))
    xin   <- cumsum(xin)
    x     <- x[xin>0]
    x     <- x[1:(length(projPeriod))]
    x[x<min(dats)]  <- min(dats)
    x[x<0]<-0
    wt[age+1,projPeriod,,,,i] <- wt[age,ac(an(projPeriod)-1),,,,i]+x
  }
}


