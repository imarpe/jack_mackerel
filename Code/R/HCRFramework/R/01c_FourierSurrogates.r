################# FOURIER SURROGATES ############################
# This function was originally written by Tristan Rouyer at Ifremer-Sète (tristan.rouyer@ifremer.fr) and is inspired from the following works:
# Theiler J, Eubank S, Longtin A, Galdrikian B, Farmer JD (1992) Testing for nonlinearity in time-series : the method of surrogate data. Physica D 58: 77-94
# Schreiber T, Schmitz, A (2000) Surrogate time series. Physica D 142: 346-382
# From a time-series vector y, the function returns a matrix of n time series with mean, variance and power spectrum identical to the original time series

FourierSurrogates <- function(y,n) {

   L<-length(y)
   z<-fft(y)
   S=matrix(NA,nrow=L,ncol=n)
   for (i in 1:n)
   {
	   if ((length(z)-trunc(length(z)/2)*2)==0) {
	      ph<-2*pi*runif((L-1)/2,0,1)
	      ph<-c(0,ph,0,-ph[length(ph):1])
	   } else {
	      ph<-2*pi*round(runif(floor(length(z)/2),0,1),4);
	      ph<-c(0,ph,-ph[length(ph):1])
	   }
   	   z<-z*exp(1i*ph)
   	   s<-Re(fft(z,inverse=T)/L)
   	   S[,i]=s
   }
return(S)
}

R           <- res$R[,1:2]
R           <- R[is.element(R[,1],as.numeric(recrPeriod)),2]
                                                                                                                                                                                                           
#- Max allowed difference between end of historical time series and start of the simulated one
maxvar      <- quantile( abs(diff(R,1)),0.75)  # the 75% quantile of observed historical variation
nextRlims   <- rev(R)[1] + c(-maxvar,maxvar)
 
#- Create the surrowgates
Rsg         <- FourierSurrogates(R,1000)

#- Select the ones which start close to the last historical R
Rsg         <- Rsg[,Rsg[1,]<max(nextRlims) & Rsg[1,]>min(nextRlims)]
idx         <- sample(1:dim(Rsg)[2],200)
Rsg         <- Rsg[,idx]

#- Shrink the low values to avoid negative rct
Rsg[Rsg<min(R)] <-  min(R)



#png(file=paste(outPath,"SGRec/Plots input/RecSims.png",sep=""),height=6,width=10,res=200,unit="in")
#
#par(mfrow=c(1,2))
#
#plot(c(R,Rsg[,1]),col=c(rep("black",length(R)),rep("red",length(R))),cex=0,main="R sim using Fourrier Surrowgates",ylab="R")
#lines(c(1:length(R)),R,col="red")
#for ( i in 1:200) points(length(R)+c(0:length(R)),c(rev(R)[1],Rsg[,i]),col="black",cex=0.1)
#abline(v=length(R))
#abline(h=0)
#q0.05 <- apply (Rsg,1,FUN=function(x) {quantile(x,0.05)})
#q0.5 <- apply (Rsg,1,FUN=function(x) {quantile(x,0.5)})
#q0.95 <- apply (Rsg,1,FUN=function(x) {quantile(x,0.95)})
#
#lines(length(R)+c(0:length(R)),c(rev(R)[1],q0.05),col="red")
#lines(length(R)+c(0:length(R)),c(rev(R)[1],q0.5),col="red")
#lines(length(R)+c(0:length(R)),c(rev(R)[1],q0.95),col="red")
#
#
#plot(c(R,Rsg[,1]),col=c(rep("black",length(R)),rep("red",length(R))),cex=0,main="R sim using Fourrier Surrowgates",ylab="R")
#lines(c(1:length(R)),R,col="red")
#abline(v=length(R))
#abline(h=0)
#for (i in  1:3)  lines(length(R)+c(0:length(R)),c(rev(R)[1],Rsg[,i]),col=i)
#
#dev.off()
#
#
#
# save the time series
Rsg         <- data.frame(t(Rsg))

save(Rsg,file=paste(dataPath,"FourrierSurrowgates200RecTS.Rdata",sep=""))


