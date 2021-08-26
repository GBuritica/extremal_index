#######################################################################
#######################################################################
#######################################################################
#######################################################################
#### Gloria Buriticá
#### Extremal index estimators : 
##                   -Intervals
##                   -Sliding
##                   -Northrop
##                   -Intervals
######### Main function
## plot_estimators plots implemented estimators
##          run plot_estimators(th0=exIndex)     If the return of exIndex has been computed previously and you want to plot results.
##          run plot_estimators(sample0,alpha0)  If estimators + plot needs to be computed.
##          set ei0 = true value of the extremal index (if known).
#########
## exIndex returns a data.frame with all estimators as a function of k : block length.
##                   Intervals estimator.
##                   Sliding estimator.
##                   Northrop estimator.
##                   SCP: spectral cluster estimator for two tunning parameters.

source("/Users/Buritica/Dropbox/Thèse/git/Auxiliar_functions/rank_transform.R")
source("/Users/Buritica/Dropbox/Thèse/git/index_regular_variation/IndexofRV.R")

#######################################################################
#######################################################################
#######################################################################
## Examples of implementation:
## Burr noise:
#set.seed(2806)
#n      <- 10000    
#par    <- 0.6
#alpha  <- 2
#ei     <- 1-par^alpha
#sample <- abs( arima.sim(n = n, list(ar=par, ma=0), rand.gen=function(n) rBurr(n,params = list(b=1,g=2,s=1)) ) )
#alpha  <-  1/alphaestimator(sample,k1=floor(n^0.8))$xi; print(paste0("alpha estimator ", alpha))
#th     <-  plot_estimators(sample,alpha0=alpha, ei0=ei)
#exIndex(sample,alpha0=alpha,k0=32)

## Student noise:
#set.seed(2806)
#n      <- 10000  
#par   <- 0.6
#alpha <- 1
#ei    <- 1-par^alpha
#sample <- abs( arima.sim(n = 3500, list(ar=par, ma=0), rand.gen=function(n) rt(n,df=alpha) ) )
#alpha  <-  1/alphaestimator(sample,k1=floor(n^0.8))$xi; print(alpha)
#th    <- plot_estimators(sample,alpha0=alpha, ei0=ei)
#exIndex(sample,alpha0=alpha,k0=32)
##################################################################
## k0   <- number of order statistics to use. 
##         By default it plots all block 2:sqrt(n) 
##         and then from k=sqrt(n) to the k=n/5 such that 1/k es an integer. 
## alpha0 <- Index of regular variation for the SCP estimator.
## if runs and blocks estimator shall be consider then set runs0 and block0 to T.

#########
## Main functions
plot_estimators<- function(sample0=NA,alpha0=NA,ei0=0,k0=5:floor(n^0.5),runs0=FALSE,blocks0=FALSE,th0=NA){
  if(is.na(th0)) th    <-  exIndex(path0=sample0,k0,alpha0=alpha0)
  else           th    <-  th0
  ## Plot
  plot(1:length(th$k),  rev(th[,1]), type = "l", ylim=c(0,1), col=4, 
       pch=4, cex=0.5, xlab ="Block length", ylab = " ",xaxt = "n", yaxt = "n")                   ## Intervals
  for(i in 2)     lines(1:length(th$k), rev(th[,i]), col = i+3 , lty=1)   
  for(i in 3)     lines(1:length(th$k), rev(th[,i]), col = "lightblue" , lty=1)   ## Sliding and Northrop
  if(runs0)   for(i in 4:5) lines(1:length(th$k), rev(th[,i]), col = "darkgrey" , lty=i+3)   ## runs
  if(blocks0) for(i in 6:7) lines(1:length(th$k), rev(th[,i]), col = "grey" , lty=i+1)       ## blocks
  for(i in 8:9)     lines( 1:length(th$k), rev(th[,i]), col = "darkblue", lty=i-7)
  
  abline(h = ei0, col ="darkgrey", lty=3)  ## true value
  legend( "topright", legend = c("Intervals","Sliding", "Northrop", "SCP: n/bl^2","SCP: n/bl^3"), 
          col =c(4:5,"lightblue" ,"darkblue", "darkblue"), lty=c(rep(1,4),2) ,cex=0.5  )
  ### Setting axes and marks
  ## Draw y-axis.
  axis(side = 2,
       las = 2,           ## Rotate labels perpendicular to y-axis.
       mgp = c(0, 0.5, 0),## Adjust y-axis label positions.
       at = c(0,ei0,1),
       labels = c(0,ei0,1))
  ## Draw the x-axis labels.
  axis(side = 1,
       las  = 1,           ## Rotate labels perpendicular to y-axis.
       mgp  = c(0, 0.5, 0),## Adjust y-axis label positions.
       at   =  seq(1,length(th$k), 10),
       labels = rev(th$k)[seq(1,length(th$k), 10)] )
  return(th)
}
exIndex        <- function(path0,k0=5:floor(n^0.5),alpha0=1,runs0=F,blocks0=F){
  n     <- length(path0)
  theta <- data.frame("Intervals" = NULL, "Sliding" = NULL, "Northrop" = NULL, "runs6"=NULL, "runs7" =NULL, "blocks6"=NULL, "blocks7"=NULL,
                      "Stable2" = NULL,"Stable3" = NULL, "k"=NULL )
  anth <- sort(path0, decreasing=TRUE )
  ## runs and blocks estimator are computed for two different threshold levels: 
  pathrank       <- n/r(path0)
  th6  <- anth[floor(n^0.6)]  
  th7  <- anth[floor(n^0.7)]
  if(min(k0) > 1){
    if(length(k0) > 1){
      re   <- floor(n/max(k0))
      kind <- c(k0, re:2)
    }
    else kind <- k0
    for(k in 1:length(kind)){
      thet <- data.frame("Intervals"  = 0, "Sliding"  = 0, "Northrop" = 0,  "runs6"= 0, "runs7" = 0, "blocks6"= 0, "blocks7"=0,
                         "Stable2" = 0, "Stable3" = 0,"k"=0)
      if( k < length(k0)){
        ## Parameters
        u  <- anth[kind[k]]              ## Defines the quantile
        Ti <- sampleInterexceeences(path0,u)  ## Interexceedence sample
        Nn <- (kind[k])-1                           ## number of exceedences
        bl <- floor(n/kind[k])                             ## block length
      }
      else{
        ## Parameters
        u  <- anth[floor(n/kind[k])]              ## Defines the quantile
        Ti <- sampleInterexceeences(path0,u)  ## Interexceedence sample
        Nn <- ((n/kind[k]))-1                     ## number of exceedences
        bl <- kind[k]                             ## block length
      }
      
      maxsliding      <- Maxsl(path0,bl,n)
      maxdisjoint     <- Maxdj(path0,bl,n)
      if(runs0==TRUE){
        ## Runs Estimator
        thet$runs6   <- thet$runs6 + min(1,(sum( sapply(1:(n-bl), function(l) (path0[l]>th6)*(maxsliding[(l+1)]<=th6) ))/floor(n^0.6)))
        thet$runs7   <- thet$runs7 + min(1,(sum( sapply(1:(n-bl), function(l) (path0[l]>th7)*(maxsliding[(l+1)]<=th7) ))/floor(n^0.7)))
      }
      if(blocks0==TRUE){
        ## Blocks estimator
        thet$blocks6  <- thet$blocks6 + min(1,  sum(maxdisjoint>th6)/floor(n^0.6))
        thet$blocks7  <- thet$blocks7 +  min(1, sum(maxdisjoint>th7)/floor(n^0.7))
      }
      ## Intervals Estimator.
      if(max(Ti) > 2){
        thet$Intervals <- thet$Intervals + Interval2(Ti,Nn)   ## Intervals if Ti > 2
      }
      else thet$Intervals <- thet$Intervals + min(1,Interval1(Ti,Nn))   ## Intervals if Ti <= 2
      
      ## Sliding Estimator
      thet$Sliding  <-  thet$Sliding  + min(1,Sliding(maxsliding,u,n,bl,Nn))      ## Sliding Robert et al.
      
      ## Northrop Estimator
      #thet$Northrop  <-  thet$Northrop + North(maxsliding,path0,n,bl)           ## Northrop
      maxslidingrank <-  Maxsl(pathrank,bl,n) #n/minR_j 
      thet$Northrop  <-  thet$Northrop + min(1,North2(maxslidingrank,n,bl))
      
      ## Stable 
      
      sumdisjoint     <- Sumdj(path0^alpha0,bl,n)
      maxdisjoint     <- maxdisjoint^alpha0
      sortsum        <- sort(sumdisjoint,decreasing=TRUE); n2 <- length(sortsum)
      #thet$Stable1    <- thet$Stable1 + stable2(sumdisjoint,maxdisjoint,bl,sortsum, max(floor( n/bl^(1.5)),2) )  ## proportion n/bl^(1.8) of sample suma
      thet$Stable2    <- thet$Stable2 + stable2(sumdisjoint,maxdisjoint,bl,sortsum, max(floor( n/bl^2),2) )  ## proportion n/bl^2 of sample suma
      thet$Stable3    <- thet$Stable3 + stable2(sumdisjoint,maxdisjoint,bl,sortsum, max(floor( n/bl^3),2) )  ## proportion n/bl^3 of sample suma
      
      #sortmax       <- sort(maxdisjoint,decreasing=TRUE)
      #thet$add  <- thet$add +  min(1, sum(maxdisjoint>sortmax[floor(n2^0.7)])/floor(n^0.7))
      
      thet$k     <- thet$k +bl
      theta      <- rbind(theta, thet)
    }
    return(theta)
  }
}
#########
## Auxiliar functions
sampleInterexceeences <- function(path0,u){
  Extimes    <-  which(path0 > u)     ## Finds index with exceedences
  Interxtimes <- Extimes[2:length(Extimes)] -  Extimes[1:(length(Extimes)-1)] ## Computes the index difference
  return(Interxtimes)
}
Sumdj     <- function(path0,b0,n0) sapply(1:max(1,floor(n0/b0)) , function(i) sum(path0[ ((i-1)*b0 + 1) : (i*b0) ]))
Maxdj     <- function(path0,b0,n0) sapply(1:max(1,floor(n0/b0)) , function(i) max(path0[ ((i-1)*b0 + 1) : (i*b0) ]))
Maxsl     <- function(path0,b0,n0) sapply(1:(n0-b0+1)    , function(i) max(path0[i:(i+b0-1)]))
Sumsl     <- function(path0,b0,n0) sapply(1:(n0-b0+1)    , function(i) sum(path0[i:(i+b0-1)]))
Interval1 <- function(Ti0,Nn0)  min( 1 , (2*sum(Ti0)^2)/( (Nn0-1)*sum(Ti0^2) )  )
Interval2 <- function(Ti0,Nn0)  min( 1 , (2*sum(Ti0-1)^2)/( (Nn0-1)*sum( (Ti0-1)*(Ti0-2) ) ) )
Sliding   <- function(Maxsl0,u,n0,b0,Nn0) -1*n0*log( mean( Maxsl0 <= u ) )/( b0*Nn0 ) 
North     <- function(Maxsl0,path0,n0 ,b0)  1/mean( -1*b0*log( sapply(Maxsl0 , function(m)  mean(path0 <= m ) ))  )
North2    <- function(Maxslrank0,n0 ,b0) 1/mean( -1*b0*log( sapply(Maxslrank0 , function(m)  (n0-(n0/m)+1)/n0  ) )  )
stable2   <- function(sumdj0, maxdj0,rl0,sortsum0,k0) mean( maxdj0[sumdj0 > sortsum0[k0]]/sumdj0[sumdj0 > sortsum0[k0]] )
