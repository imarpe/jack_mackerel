#-----------------------------------------------------------------------------------------------------------
# Recruitment function => gives the R accoring to the SSB and SR parameters
#-----------------------------------------------------------------------------------------------------------

B2R<-function (Biom, SRp,yr,recreg){

        if( recreg  !=  "SGRec"){
          SRmod <-  SRp[,1]
          A     <-  SRp[,2]
          B     <-  SRp[,3]
          sig   <-  SRp[,4]* 1000  # this rescaling should have been done in the Bayesian protocol


          idxHSL<- which(SRmod=="HSL")
          idxRKL<- which(SRmod=="RKL")
          idxBHL<- which(SRmod=="BHL")
          mu    <- numeric(length(A))
          mu[idxHSL] <- A[idxHSL]*B[idxHSL]*(c(Biom)[idxHSL]>=B[idxHSL])+A[idxHSL]*c(Biom)[idxHSL]*(c(Biom)[idxHSL]<B[idxHSL])
          mu[idxRKL] <- A[idxRKL]*c(Biom)[idxRKL]*exp(-B[idxRKL]*c(Biom)[idxRKL])
          mu[idxBHL] <- A[idxBHL]*c(Biom)[idxBHL]/(B[idxBHL]+c(Biom)[idxBHL])
          REC        <- pmax(rnorm(length(mu),mu,sig),10)
        }
        if (recreg == "SGRec")
          REC       <- SRp[,yr]

    return(REC)}

#-----------------------------------------------------------------------------------------------------------
#  Runs the pseudo assessment
#-----------------------------------------------------------------------------------------------------------

assessJM    <-  function (iStck,iDevN,iDevF,iCVn,iCVf,yr){
                  # years
                  TaY       <- (yr-1)
                  ImY       <- (yr)
                  FcY       <- (yr+1)
                  year1     <- (dimnames(stock.n(iStck))$year)[1]

                  # creates an empty FLquant of the right dimensions
                  quant     <- stock.n(iStck)
                  quant[]   <- NA

                  # update date in the FLquant for error amplitude  on N
                  idx       <- unique(which(is.na(iCVn)==F & iCVn > 0,arr.ind=T)[,2])
                  iCVn      <- iCVn[,idx]
                  # move the error period to the recent years
                  dimnames(iCVn)$year                       <- ac((TaY-(length(idx)-1)):TaY)
                  quant[,dimnames(iCVn)$year]               <- iCVn
                  iCVn      <- quant ;  quant[]             <- NA
                  iCVn[,ac(dimnames(iCVn)$year[1]:(TaY-length(idx)))] <- 0


                  # update date in the FLquant for error amplitude  on F
                  idx       <- unique(which(is.na(iCVf)==F & iCVf > 0,arr.ind=T)[,2])
                  iCVf      <- iCVf[,idx]
                  # move the error period to the recent years
                  dimnames(iCVf)$year                       <- ac((TaY-(length(idx)-1)):TaY)
                  quant[,dimnames(iCVf)$year]               <- iCVf
                  iCVf      <- quant ;  quant[]             <- NA
                  iCVf[,ac(dimnames(iCVf)$year[1]:(TaY-length(idx)))] <- 0

                  # observed stock : true stock + deviations
                  # with deviation = cohort coherent deviation (devN) and year dependent ampliture (JM@stock.n*CV_stock.n)
                  obsstck                                   <- iStck
                  print("Warning: stock assessment object is direct copy from biological object, except for harvest and numbers-at-age")
                  stock.n(obsstck)                          <- iStck@stock.n + iDevN * iStck@stock.n * iCVn
                  harvest(obsstck)                          <- iStck@harvest + iDevF * iStck@harvest * iCVf

                  #  stock.n(obsstck)[1,ac((an(TaY)-9):an(TaY)),1,,,its]
                  stock.n(obsstck)[stock.n(obsstck)<=0 & !is.na(stock.n(obsstck))]  <- 10
                  harvest(obsstck)[harvest(obsstck)<0 & !is.na(harvest(obsstck))]   <-  0
                return(obsstck)}

#-----------------------------------------------------------------------------------------------------------
#  Finds the TAC correspoding to the application of the HCR based on the observed stock
#-----------------------------------------------------------------------------------------------------------

# main function
findTAC <-
function(iStf,iRef,iHCR,iTaY,iImY,iFcY){

# iStf<-stf ; iRef<-Ref  ;  iHCR<-HCR  ; iTaY<-TaY  ;  iImY<-ImY   ; iFcY<-FcY
  
  
      if(is.null(iHCR$Ftarget) & is.null(iHCR$Btarget) & is.null(iHCR$slope) & is.null(iHCR$Recov)) stop("Specify either Ftarget, Btarget or slope")
      if(is.null(iHCR$slope)==F & (is.null(iHCR$slopeBmin) | is.null(iHCR$slopeBmax) | is.null(iHCR$slopeFmin) | is.null(iHCR$slopeFmax))) stop("When specifying slope, also specify slopeBmin, slopeBmax,slopeFmin and slopeFmax")
      if(is.null(iHCR$Ftarget)==F & is.null(iHCR$Btarget)==F) stop("Specify either Ftarget or Btarget, not both")

      # find the F multiplicator for each iteration
      iNits                 <- dims(iStf)$iter

      # Assume similar harvest pattern for forecast year
      iStf@harvest[,FcY]    <- iStf@harvest[,ImY]

      if(is.null(iHCR$Ftarget) == F){
        multip              <- rep(NA,iNits)
        for(iTer in 1:iNits)
          multip[iTer]      <- iHCR$Ftarget/quantMeans(areaSums(harvest(iStf[,ImY,,,,iTer])))
      }
      if(is.null(iHCR$Btarget) == F){
        multip              <- rep(NA,iNits)
        for(iTer in 1:iNits)
          multip[iTer]      <- optimize(f=function(mult,x,target,FcY){x@harvest[,FcY] <- x@harvest[,FcY] * mult;
                                          return(sqrt((ssb(x[,FcY]) - target)^2))},
                                        lower=0,upper=10,x=iter(iStf,iTer),target=iHCR$Btarget,FcY=FcY)$minimum
      }
      if(is.null(iHCR$slope) == F){
        multip              <- rep(NA,iNits)
        for(iTer in 1:iNits)
          multip[iTer]      <- optimize(f=function(mult,x,ImY,FcY){
                                            x@harvest[,FcY] <- x@harvest[,ImY] * mult;
                                            iSSB    <- ssb(x[,FcY])@.Data;
                                            target  <- ifelse(iSSB >= iHCR$slopeBmax,iHCR$slopeFmax,
                                                              ifelse(iSSB<=iHCR$slopeBmin,iHCR$slopeFmin,
                                                                     iHCR$slopeFmin + (iSSB - iHCR$slopeBmin) * (iHCR$slopeFmax - iHCR$slopeFmin) / (iHCR$slopeBmax - iHCR$slopeBmin)))
                                          return(sqrt((fbar(x[,FcY]) - target)^2))},
                                        lower=0.1,upper=3,x=iter(iStf,iTer),ImY=ImY,FcY=FcY)$minimum
      }
      
      
    
    
    
    
     if(is.null(iHCR$Recov) == F){
        multip              <- rep(NA,iNits)
        for(iTer in 1:iNits)
        ifelse( ssb(stf[,ImY,,,,iTer])@.Data<iRef$Blim   ,  multip[iTer] <- 0.75   , 
                                                      ifelse(ssb(iStf[,ImY,,,,iTer])@.Data>iRef$Bpa,multip[iTer] <- iRef$Fmsy/fbar(iStf[,ImY,,,,iTer]),multip[iTer] <- 0.85))
      }  
      
      
      
      
      
      
      

      iStf@harvest[,FcY]  <- sweep(iStf@harvest[,FcY],6,multip,"*")
      iStf@catch.n[,iFcY] <- sweep(sweep(iStf@harvest[,iFcY],c(1:4,6),(areaSums(iStf@harvest[,iFcY]) + areaMeans(iStf@m[,iFcY])),"/") *
                                         iStf@stock.n[,iFcY],c(1:4,6),
                                        (1-exp(-(areaSums(iStf@harvest[,iFcY]))-areaMeans(iStf@m[,iFcY]))),"*")
      iTac                <- quantSums(areaSums(iStf@catch.n[,iFcY]*iStf@catch.wt[,iFcY]))
  return(iTac)}




        
#-----------------------------------------------------------------------------------------------------------
#  Finds the Fs for the true stock for the forecast year, corresponding to the TAC
#-----------------------------------------------------------------------------------------------------------
TAC2F     <- function(iStk,iFleet,iTaY,iImY,iFcY,iTac){
              # stk <-  JMB  ; fleet <-  fishery ; yr  <-  I ; tac <-  TAC
              # find the real Fbar multiplier
              iNits                       <- dims(iStk)$iter
              
              #- Scale the landings.sel to the proportion F by fleet (landings.sel = standardized and therefore each fleet contributes for 1/fleets to the total F)
              iFleet@landings.sel[,iImY]  <- sweep(landings.sel(iFleet[,iImY]),c(5,6),sweep(quantMeans(iStk@harvest[,iTaY]),6,fbar(iStk[,iTaY]),"/"),"*")

              scaleHarvest  <- numeric(nits)
              for(iTer in 1:iNits)
                scaleHarvest[iTer]  <- optimize(f=function(x,iiTac,iiStk,iiFleet,iiTaY,iiImY,iiFcY){
                                                   res <- sqrt((iiTac -
                                                               sum(sweep(sweep(landings.sel(iiFleet)[,drop=T] * x,1,(rowSums(landings.sel(iiFleet)[,drop=T],na.rm=T) * x) + rowMeans(iiStk@m[,drop=T],na.rm=T),"/") *
                                                                         iiStk@stock.n[,drop=T],1,
                                                                         (1-exp(-rowSums(landings.sel(iiFleet)[,drop=T] * x,na.rm=T) - rowMeans(iiStk@m[,drop=T],na.rm=T))),"*") * iiStk@catch.wt[,drop=T],na.rm=T))^2)
                                                  return(res)},
                                                interval=c(0,10),iiTac=c(iter(iTac[,iImY],iTer)),iiStk=iter(iStk[,iImY],iTer),iiFleet=iter(iFleet[,iImY],iTer),iiTaY=iTaY,iiImY=iImY,iiFcY=iFcY)$minimum
              harvest(iStk)[,ImY] <- sweep(landings.sel(iFleet)[,ImY],6,scaleHarvest,"*")
            return(iStk)}

#-------------------------------------------------------------------------------
# areaSums newly defined
#-------------------------------------------------------------------------------

areaSums <- function(x){return(apply(x,c(1:4,6),sum,na.rm=T))}

#-------------------------------------------------------------------------------
# as.numeric
#-------------------------------------------------------------------------------

an <- function(x){return(as.numeric(x))}

#-------------------------------------------------------------------------------
# SSB for the entire stock
#-------------------------------------------------------------------------------

ssb      <- function(x){return(quantSums(areaMeans(stock.n(x))*areaMeans(mat(x))*areaMeans(stock.wt(x)) *
                               exp(-areaSums(x@harvest)*areaMeans(harvest.spwn(x))-areaMeans(x@m)*areaMeans(m.spwn(x)))))}
                               
#-------------------------------------------------------------------------------
# Fbar for the entire stock
#-------------------------------------------------------------------------------

fbar      <- function(x){return(quantMeans(areaSums(harvest(x))[ac(range(JMB)["minfbar"]:range(JMB)["maxfbar"]),]))}

#-------------------------------------------------------------------------------
# Rec for the entire stock
#-------------------------------------------------------------------------------


rec       <- function(x){return(x@stock.n[ac(range(x)["min"]),,,,1])}

#-------------------------------------------------------------------------------
# Calculate survivors
#-------------------------------------------------------------------------------

calcSurvivors <- function(iStk){
                    yrs       <- dimnames(iStk@stock.n)$year[-length(dimnames(iStk@stock.n)$year)]
                    survivors <- areaMeans(stock.n(iStk[,yrs])) * exp(-areaSums(iStk[,yrs]@harvest)-areaMeans(iStk[,yrs]@m))
                    stock.n(iStk)[ac((range(iStk,"min")+1):range(iStk,"max")),ac(an(yrs)+1),,,1] <- survivors[-dim(survivors)[1],,,,,]@.Data
                    for(iArea in 2:length(dimnames(iStk@stock.n)$area))
                      stock.n(iStk)[ac((range(iStk,"min")+1):range(iStk,"max")),ac(an(yrs)+1),,,iArea] <-
                        stock.n(iStk)[ac((range(iStk,"min")+1):range(iStk,"max")),ac(an(yrs)+1),,,1]
                    #- Plusgroup
                    if (!is.na(range(iStk,"plusgroup"))){
                      stock.n(iStk)[ac(range(iStk,"max")),ac(an(yrs)+1),,,1] <- areaMeans(stock.n(iStk)[ac(range(iStk,"max")),ac(an(yrs)+1)]) + survivors[ac(range(iStk,"max"))]
                    }
                    for(iArea in 2:length(dimnames(iStk@stock.n)$area))
                      stock.n(iStk)[ac(range(iStk,"max")),ac(an(yrs)+1),,,iArea] <-
                        stock.n(iStk)[ac(range(iStk,"max")),ac(an(yrs)+1),,,1]
                 return(iStk@stock.n[,ac(an(yrs)+1)])}

#-------------------------------------------------------------------------------
# Calculate catch
#-------------------------------------------------------------------------------

computeCatch <- function(stck){
                    yr <- dimnames(stck@stock.n)$year
                    res <- quantSums(areaSums(sweep(sweep(stck@harvest[,yr],c(1:4,6),(areaSums(stck@harvest[,yr]) + areaMeans(stck@m[,yr])),"/"),c(1:4,6),
                                              areaMeans(stck@stock.n[,yr]) * (1-exp(-areaSums(stck@harvest[,yr]) -areaMeans(stck@m[,yr]))),"*") * stck@catch.wt[,yr]))
                return(res)}
                
                
computeCatchperFleet <- function(stck){
                    yr <- dimnames(stck@stock.n)$year
                    res <- quantSums((sweep(sweep(stck@harvest[,yr],c(1:4,6),(areaSums(stck@harvest[,yr]) + areaMeans(stck@m[,yr])),"/"),c(1:4,6),
                                              areaMeans(stck@stock.n[,yr]) * (1-exp(-areaSums(stck@harvest[,yr]) -areaMeans(stck@m[,yr]))),"*") * stck@catch.wt[,yr]))
                return(res)}    
                

       
                
computeCatchatAgeperFleet <- function(stck){
                    yr <- dimnames(stck@stock.n)$year
                    res <- ((sweep(sweep(stck@harvest[,yr],c(1:4,6),(areaSums(stck@harvest[,yr]) + areaMeans(stck@m[,yr])),"/"),c(1:4,6),
                                              areaMeans(stck@stock.n[,yr]) * (1-exp(-areaSums(stck@harvest[,yr]) -areaMeans(stck@m[,yr]))),"*") * stck@catch.wt[,yr]))
                return(res)}                                

#-------------------------------------------------------------------------------
# Calculate landings
#-------------------------------------------------------------------------------

computeLandings <- function(stck){
                    yr <- dimnames(stck@stock.n)$year
                    res <- quantSums(areaSums(sweep(sweep(stck@harvest[,yr],c(1:4,6),(areaSums(stck@harvest[,yr]) + areaMeans(stck@m[,yr])),"/"),c(1:4,6),
                                              areaMeans(stck@stock.n[,yr]) * (1-exp(-areaSums(stck@harvest[,yr]) -areaMeans(stck@m[,yr]))),"*") * stck@catch.wt[,yr]))
                return(res)}
                
#-------------------------------------------------------------------------------
# Arma search function
#-------------------------------------------------------------------------------


armaSearch = function(
   xx,
   minOrder=c(0,0),
   maxOrder=c(5,5),
   trace=FALSE )
{
   bestAic = 1e9
   len = NROW( xx )
   for( p in minOrder[1]:maxOrder[1] ) for( q in minOrder[2]:maxOrder[2] )
   {
      if( p == 0 && q == 0 )
      {
         next
      }

      formula = as.formula( paste( sep="", "xx ~ arma(", p, ",", q, ")" ) )

      fit = tryCatch( armaFit( formula, data=xx ),
                      error=function( err ) FALSE,
                      warning=function( warn ) FALSE )
      if( !is.logical( fit ) )
      {
         fitAic = fit@fit$aic
         if( fitAic < bestAic )
         {
            bestAic = fitAic
            bestFit = fit
            bestModel = c( p, q )
         }

         if( trace )
         {
            ss = paste( sep="", "(", p, ",", q, "): AIC = ", fitAic )
            print( ss )
         }
      }
      else
      {
         if( trace )
         {
            ss = paste( sep="", "(", p, ",", q, "): None" )
            print( ss )
         }
      }
   }

   if( bestAic < 1e9 )
   {
      return( list( aic=bestAic, fit=bestFit, model=bestModel ) )
   }

   return( FALSE )
}
