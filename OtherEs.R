#######################################################################
#######################################################################
#######################################################################
#######################################################################
#### Gloria Buriticá
#### Estimatiors for SS
#### 1. Robertetal ( Depends on a proportion k of the path )
#### 2. Northrop
#######################################################################
#######################################################################
#######################################################################
#######################################################################
# rn  <- parameter block size.
# num <- parameter  num = floor(length(path)*0.0x)
library(gghighlight)
## Robertetal
Index_Robertelat <- function(rn,num,path){
  path       <- r(path)       ## Getting the rank's X1 > X2 > X3 ..
  n          <- length(path)  ## Pathlength
  pathmax    <- path[1:(n-rn+1)] ## Creating pathmax
  pathmax    <- sapply(1:length(pathmax), function(k) min(path[k:(k+rn-1)]) )
  est        <- -n*log( (1/(n-rn+1))*sum(pathmax > num) )/((num-1)*rn)
  return(min(est,1))
}
## Northrop
Index_Northrop1 <- function(path,rn){
  kn         <- floor(length(path)/rn)             ## Dismissing data at the end
  n          <- kn*rn
  path_rank  <- r(path[1:n],decreasing = FALSE)                       ## Returns the path of Rank's
  ########## Vector of max for disjoint blocks
  vect_max   <- sapply(1:kn , function(k)  max(path_rank[ (( (k-1)*rn)+1) : (k*rn) ]  ))
  ########## Mean-cluster-blocksize
  rec <- (1/kn)*sum( -rn*log( (1/n)*vect_max) )
  return(1/rec)
}
## Intervals Ferro et al. Estimator (u)
# T_i := S_{i+1} - S_{i} exceedence time. There is an exceedence at time i. 



re <- floor(n/max(k0))
for(k in re:2){
  thet <- data.frame("Intervals" = 0, "Sliding" = 0, "Northrop" = 0,  
                     "Stable23" = 0, "Stable24" = 0,"Stable25" =0,"Stable" = 0, "Stablep" = 0, "Stablep2" = 0,"Stablep22"=0,"k"=0)
  ## Parameters
  u  <- anth[floor(n/k)]              ## Defines the quantile
  Ti <- sampleInterexceeences(path0,u)  ## Interexceedence sample
  Nn <- ((n/k)-2)                           ## number of exceedences
  bl <- k                             ## block length
  maxsliding    <- Maxsl(path0,bl,n)
  ## Intervals Estimator.
  if(max(Ti) > 2) thet$Intervals <- thet$Intervals + Interval2(Ti,Nn)   ## Intervals if Ti > 2
  else            thet$Intervals <- thet$Intervals + Interval1(Ti,Nn)   ## Intervals if Ti <= 2
  ## Sliding Estimator
  thet$Sliding  <-  thet$Sliding  + Sliding(maxsliding,u,n,bl,Nn)       ## Sliding Robert et al.
  ## Northrop Estimator
  #thet$Northrop  <-  thet$Northrop + North(maxsliding,path0,n,bl)           ## Northrop
  maxslidingrank <-  Maxsl(pathrank,bl,n) #n/minR_j 
  thet$Northrop  <-  thet$Northrop + North2(maxslidingrank,n,bl)
  ## Stable2
  
  zn            <-  u^alpha0 + bu[1]*mean((path0^alpha0)*(path0 <= u)) 
  zn2           <-  u^alpha0 + bu[2]*mean((path0^alpha0)*(path0 <= u)) 
  zn3           <-  u^alpha0 + bu[3]*mean((path0^alpha0)*(path0 <= u)) 
  zn4           <-  u^alpha0 + bu[4]*mean((path0^alpha0)*(path0 <= u)) 
  
  
  thet$Stablep    <- thet$Stablep +  min((mean( maxsliding2 > u))/(mean( sumsliding2 > zn2)),1)
  thet$Stable     <- thet$Stable +  min((mean( maxsliding1 > u))/(mean( sumsliding1 > zn)),1)
  thet$Stablep2   <- thet$Stablep2 + min((mean( maxsliding3 > u))/(mean( sumsliding3 > zn3)),1)
  thet$Stablep22  <- thet$Stablep22 + min((mean( maxsliding4 > u))/(mean( sumsliding4 > zn4)),1)
  
  k <- bl
  sumdisjoint    <- Sumdj(path0^alpha0,k,n)
  maxdisjoint    <- Maxdj(path0^alpha0,k,n)
  sortsum        <- sort(sumdisjoint); n2 <- length(sortsum)
  thet$Stable23  <- thet$Stable23 + stable2(sumdisjoint,maxdisjoint,k,sortsum, n2*(1-1/n^(0.3) ) ) ## proportion 1/n^(0.3) of sample suma
  thet$Stable24  <- thet$Stable24 + stable2(sumdisjoint,maxdisjoint,k,sortsum, n2*(1-1/n^(0.4) ) ) ## proportion 1/n^(0.4) of sample suma
  thet$Stable25  <- thet$Stable25 + stable2(sumdisjoint,maxdisjoint,k,sortsum, n2*(1-1/n^(0.5) ) ) ## proportion 1/n^(0.5) of sample suma
  ## Stable1log
  #bl <- 2
  #sumsliding     <- Sumsl(path0^alpha0,(bl),n)
  #maxsliding     <- Maxsl(path0,(bl),n)
  #znk            <- u^alpha0 + bl*mean((path0^alpha0)*(path0 <= u)) 
  #thet$Stable    <- thet$Stable +  (mean( maxsliding[(sumsliding > znk)] > u) - mean(maxsliding > u))/(mean( sumsliding[(sumsliding > znk)] > u) - mean((sumsliding > u)))
  
  #thet$Stablep22  <- thet$Stablep22 +
  #    (pevd( u, loc=fit$results$par[1], scale=fit$results$par[2], shape=fit$results$par[3] , type = "GEV" )/pstable(zn4, alpha=fit1[1], beta=1,gamma=fit1[3],delta=fit1[4] )) 
  # thet$Stablep2 <- thet$Stablep2 + (thet$Stable*thet$Stable25^2 +  thet$Stable25*(1-thet$Stable25^2) )
  thet$k <- thet$k +k
  
  theta        <- rbind(theta, thet)
}
## k0 <- number of order statistics to use
## b0 <- block length to use 
exIndex <- function(path0,k0,alpha0=1){
  n     <- length(path0)
  theta <- data.frame("Intervals" = NULL, "Sliding" = NULL, "Northrop" = NULL, "runs6"=NULL, "runs7" =NULL, "blocks6"=NULL, "blocks7"=NULL,
                      "Stable26" = NULL,"Stable27" = NULL, "Stable28" = NULL,"Stable29" = NULL, "k"=NULL )
  anth <- sort(path0, decreasing=TRUE )
  
  pathrank       <- n/r(path0)
  th6  <- anth[floor(n^0.6)]
  th7  <- anth[floor(n^0.7)]
  #print("ok")
  if(min(k0) > 1){
    re   <- floor(n/max(k0))
    kind <- c(k0, re:2)
    for(k in 1:length(kind)){
      thet <- data.frame("Intervals" = 0, "Sliding" = 0, "Northrop" = 0,  "runs6"=0, "runs7" =0, "blocks6"=0, "blocks7"=0,
                           "Stable26" = 0, "Stable27" = 0,"Stable28" = 0,"Stable29" = 0,"k"=0)
      if( k <= length(k0)){
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
      ## Runs Estimator
      thet$runs6   <- thet$runs6 + min(1,(sum( sapply(1:(n-bl), function(l) (path0[l]>th6)*(maxsliding[(l+1)]<=th6) ))/floor(n^0.6)))
      thet$runs7   <- thet$runs7 + min(1,(sum( sapply(1:(n-bl), function(l) (path0[l]>th7)*(maxsliding[(l+1)]<=th7) ))/floor(n^0.7)))
      ## Blocks estimator
      thet$blocks6  <- thet$blocks6 + min(1, sum(maxdisjoint>th6)/floor(n^0.6))
      thet$blocks7  <- thet$blocks7 +  min(1, sum(maxdisjoint>th7)/floor(n^0.7))
     
     
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
     thet$Stable26  <- thet$Stable26 + stable2(sumdisjoint,maxdisjoint,bl,sortsum, max(floor( (n^0.6)/bl),2) )  ## proportion 1/n^(0.5) of sample suma
     thet$Stable27  <- thet$Stable27 + stable2(sumdisjoint,maxdisjoint,bl,sortsum, max(floor( (n^0.7)/bl),2) )  ## proportion 1/n^(0.5) of sample suma
     thet$Stable28  <- thet$Stable28 + stable2(sumdisjoint,maxdisjoint,bl,sortsum, max(floor( (n^0.8)/bl),2) )  ## proportion 1/n^(0.5) of sample suma
     thet$Stable29  <- thet$Stable29 + stable2(sumdisjoint,maxdisjoint,bl,sortsum, max(floor( (n^0.9)/bl),2) )  ## proportion 1/n^(0.5) of sample suma
     
     #sortmax       <- sort(maxdisjoint,decreasing=TRUE)
     #thet$add  <- thet$add +  min(1, sum(maxdisjoint>sortmax[floor(n2^0.7)])/floor(n^0.7))
    
     thet$k     <- thet$k +bl
     theta      <- rbind(theta, thet)
    }
    return(theta)
  }
}
##
par(mfrow=c(2,2))
par(mfrow=c(1,1))
n      <-  10000; ei <- 0.5; kmax <- floor(n^0.8); kmax2 <- floor(n^0.5)
#sample <- ARMAX1(1-ei,n,al=1) #abs(ARCHmodel(lambdaV[2,4],n)); ei <- thetaV[2,4]## ARMAX1(1-ei,n) #n/r(ARMAX1(1-ei,n))#  ARMAX1(1-ei,n)#
sample <- abs( arima.sim(n = 5000, list(ar=0.5, ma=0), rand.gen=function(n)rt(n,df=1)) )
par(mfrow=c(1,1))
{

alpha <- 1/alphaestimator(sample,k0=kmax, R0 = 1)$xi; print(alpha)
th    <- exIndex(path0=sample,k0=5:kmax2,alpha0=alpha)
## Plot

plot(log(th$k),  (th[,1]), type = "p", ylim=c(0,1), col=4, 
     pch=4, cex=0.5, xlab ="log( block length )", ylab = "Extremal Index")   ## Intervals
lines(log(th$k), (th[,1]), col = 4, lty=4)
for(i in 2:3){ 
#  points( log(th$k), th[,i], col = i+3, pch=i+3, cex =0.5)   ## Sliding-Northrop
  lines(log(th$k), th[,i], col = i+3 , lty=i+3)
}
for(i in 4:5){ 
  #points( log(th$k), th[,i], col = 7, pch=i+3, cex =0.5)   ## Sliding-Northrop
  lines(log(th$k), th[,i], col = 7 , lty=i+3)
}
for(i in 6:7){ 
  #points( log(th$k), th[,i], col = 9, pch=i+3, cex =0.5)   ## Sliding-Northrop
  lines(log(th$k), th[,i], col = 9 , lty=i+3)
}
for(i in 8:11){
  lines( log(th$k), (th[,i]), col = 1, lty=i-3);#points( log(th$k), (th[,i]), col = 1, lty=i-3,cex=0.5,pch=i-3)
}
abline(h=0.2)
##
#for(i in 4:6) points( log(th$k), (th[,i]), col = 1, lty=i, pch=16, cex=0.3) ## Stable ^0.3,0.4,0.5

#for(i in 7)   lines( 1:length(th[,1]), th[,i], col = i, lty=i)   ## Stable Biais
#for(i in 8:10) lines( log(th$k), (th[,i]), col = i, lty=i)   ## Rapport + mean
## Legend.
abline(h = ei, col ="red")
legend( "topright", legend = c("Intervals","Sliding", "Northrop", "LD n^0.3","LD n^0.4","LD n^0.5"), 
        col =c(4:6, "black", "black", "black"), lty=c(4:6,1:3),pch=c(4:6,1:3),cex=1  )
}
############################## Additional functions.
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
r <- function(x){
  v <- vector(length = length(x), mode = "integer")
  y <- order(x, decreasing=TRUE)
  for(i in 1:length(x))  v[y[i]] = i 
  return(v)
}
##################################################################
##################################################################
##################################################################
##################################################################
##################################################################
# Simulation study
type          <- c( "Intervals","Sliding","Northrop","Stable23","Stable24","Stable25","Stable","Stablep" ,"Stablep2","Stablep22")
infothetSS    <-  data.frame(  "Intervals"=NULL,"Sliding"=NULL,"Northrop"=NULL,"Stable23"=NULL,"Stable24"=NULL,
                               "Stable25"=NULL,"Stable"=NULL,"Stablep"=NULL ,"Stablep2"=NULL,"Stablep22"=NULL,  ## Estimators
                               "k"=NULL, "thet" = NULL, "alpha" = NULL,"N" = NULL, "id" = NULL)
##################################################################
## ARMAX 
{
n      <-  3500; ei <- 0.2; kmax <- n^0.8; kmax2 <- n^0.5; 
id <- 1; theta <- ei;  km    <-   114# length(th$k)
for(N  in 1:1000){
  sample  <-  (ARMAX1(1-ei,n))
  alpha  <-   1/alphaestimator(sample,k0=kmax1)$xi;# print(alpha)
  th     <-   exIndex(sample,k0=4:kmax2,p0=p,alpha0=alpha)
  th$k  <-    1:km
  th$thet <- rep(theta, floor(km) )
  th$alpha <- rep(alpha, floor(km) )
  th$N     <-  rep(N, floor(km) )
  th$id    <-  rep(id, floor(km) )
  infothetSS <- rbind( infothetSS, th )
  print(N)
}
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "/Users/Buritica/Dropbox/Thèse/0/Projet_EI/SS/SSData/sim0904ARMAX0.2")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.6")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.4")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.2")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.6")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.4")
}
##################################################################
## sqARCH2 
{
n      <-  5000; ei <- thetaV[1,3]; kmax <- n^0.8; 
kmax2 <- n^0.5; id <- 1; theta <- ei;km     <- 137#114
for(N in 1:500){
  sample <-  (squaredARCH(lambdaV[1,3],n))
  alpha  <-  1/alphaestimator(sample,k0=kmax)$xi; #print(alpha)
  th     <-  exIndex(sample,k0=4:kmax2,p0=p,alpha0=alpha)
  th$k   <-  1:km; 
  th$thet <- rep(theta, floor(km) )
  th$alpha <- rep(alpha, floor(km) )
  th$N     <-  rep(N, floor(km) )
  th$id    <-  rep(id, floor(km) )
  infothetSS <- rbind( infothetSS, th )
  print(N)
}
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904sqARCH2")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904sqARCH3")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904sqARCH4")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904sqARCHRank2")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904sqARCHRank3")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904sqARCHRank4")
##
#save(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904sqARCH310000")
#save(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904sqARCH410000")
}
###################################################################
## ARCH
{
n      <-  3500; ei <- thetaV[2,4]; kmax <- n^0.8; 
kmax2 <- n^0.5; id <- 1; theta <- ei;km     <- 114
for(N in 1:500){
  sample <-  (abs(ARCHmodel(lambdaV[2,4],n)))
  alpha  <-  1/alphaestimator(sample,k0=kmax)$xi;# print(alpha)
  th     <-  exIndex(sample,k0=4:kmax2,p0=p,alpha0=alpha)
  th$k  <-    1:km
  th$thet <- rep(theta, floor(km) )
  th$alpha <- rep(alpha, floor(km) )
  th$N     <-  rep(N, floor(km) )
  th$id    <-  rep(id, floor(km) )
  infothetSS <- rbind( infothetSS, th )
  print(N)
}
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARCH2")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARCH3")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARCH4")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARCHRank2")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARCHRank3")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARCHRank4")
##
#save(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARCH210000")
#save(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARCH310000")
#load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARCH410000")
}
###################################################################
##################################################################
##### Plot mean
mean  <- sapply(1:km, function(l) (mean( (infothetSS[infothetSS$k==l ,][ ,1] ) ) ))
var   <- sapply(1:km, function(l) (mean( (infothetSS[infothetSS$k==l ,][ ,1] - mean[l] )^2 ) )) 
biais <- sapply(1:km, function(l) (mean[l] - theta )^2 ) 
MSE   <-  sapply(1:km, function(l) (mean( (infothetSS[infothetSS$k==l ,][ ,1] -theta )^2 ) ))
for(i in 2:10){
  me   <- sapply(1:km, function(l) (mean( (infothetSS[infothetSS$k==l ,][ ,i] ) ) ))
  va   <- sapply(1:km, function(l) (mean( (infothetSS[infothetSS$k==l ,][ ,i] - me[l] )^2 ) )) 
  bi   <- sapply(1:km, function(l) (me[l] - theta )^2 )  
  MSer <-  sapply(1:km, function(l) (mean( (infothetSS[infothetSS$k==l ,][ ,i] -theta )^2 ) ))
  mean <- cbind(mean,me )
  var   <- cbind(var,va)
  biais <- cbind(biais, bi)
  MSE   <- cbind(MSE, MSer)
}  
{
  par(mfrow = c(1,3))
  #plot( log(th$k),  (mean[,1]), type = "l" , ylim = c(0,1),col=4, xlab="block length = log(k)", ylab ="Mean" )
  #points( log(th$k),  (mean[,1]), col=4, pch=16, cex = 0.3)
  #for(i in 2:3) lines( log(th$k), (mean[,i])    , col =i+3)
  #for(i in 4:6) lines( log(th$k), (mean[,i])    , col =1, lty = i)
  #for(i in 7:10) lines((4:kmax2), mean[,i]    ,col ="grey", lty = i)
  #abline(h = theta, col= "red")
  #legend( "bottomright", legend = c("Intervals","Sliding", "Northrop", "LD Threshold", "Mod. Threshold"), 
  #        col =c(4:6, "black", "grey"), lty=1,cex=0.6 , horiz=TRUE )
  #abline(h = theta, col= "red")
  
  plot( log(th$k),  (var[,1]), type = "l" , ylim = c(0,0.04) , col=4, xlab="log( block length )", ylab ="Variance", lty = 4)
  for(i in 2:3){lines(log(th$k), (var[,i])    , col =i+3, lty = i+3)}
  for(i in 4:6) lines( log(th$k), (var[,i])    , col =1, lty = i-3)
  #for(i in 7:10) lines( (4:kmax2), var[,i]   ,col ="grey", lty = i)
  legend( "topright", legend = c("Intervals","Sliding", "Northrop", "LD n^0.3","LD n^0.4","LD n^0.5"), 
          col =c(4:6, "black", "black", "black"), lty=c(4:6,1:3),cex=1  )
  
  plot( log(th$k),  (biais[,1]), type = "l" , ylim = c(0,0.04), col=4, xlab="log( block length )", ylab ="Biais^2",lty=4)
  for(i in 2:3) lines(log(th$k),(biais[,i])    , col =i+3, lty=i+3)
  for(i in 4:6) lines( log(th$k),(biais[,i])    , col =1, lty = i-3)
  #for(i in 7:10) lines((4:kmax2),biais[,i]   ,col ="grey", lty = i)
  legend( "topright", legend = c("Intervals","Sliding", "Northrop", "LD n^0.3","LD n^0.4","LD n^0.5"), 
          col =c(4:6, "black", "black", "black"), lty=c(4:6,1:3),cex=1  )
  
  plot( log(th$k),(MSE[,1]), type = "l" , ylim = c(0,0.04), 
        col=4, xlab="log( block length )", ylab ="MSE",lty=4)
  for(i in 2:3) lines( log(th$k), (MSE[,i]) , col =i+3, lty = i+3)
  for(i in 4:6) lines(log(th$k), (MSE[,i])    , col =1, lty = i-3)
  #for(i in 7:10) lines( (4:kmax2), rev(MSE[,i])   ,col ="grey", lty = i)
  legend( "topright", legend = c("Intervals","Sliding", "Northrop", "LD n^0.3","LD n^0.4","LD n^0.5"), 
          col =c(4:6, "black", "black", "black"), lty=c(4:6,1:3),cex=1  )
  
} ## Plot. 


###################################################################
###################################################################
## SS2
## ARCH(1) - rank
{
  n      <-  5000; ei <- 0.2; kmax <- floor(n^0.8); kmax2 <- floor(n^0.5); 
  id <- 1; theta <- ei;  km    <-   137# length(th$k)
  for(N  in 1:10){
    sample  <-  abs(ARCHm(n))
    sample  <- n/r(sample)
    alpha  <-   1#1/alphaestimator(sample,k1=kmax)$xi ;# print(alpha)
    th     <-   exIndex(path0=sample,k0=4:kmax2,alpha0=alpha)
    th$k  <-    1:km
    th$thet <- rep(theta, floor(km) )
    th$alpha <- rep(alpha, floor(km) )
    th$N     <-  rep(N, floor(km) )
    th$id    <-  rep(id, floor(km) )
    infothetSS <- rbind( infothetSS, th )
    print(N)
  }
  head(infothetSS)
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.2")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.6")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.4")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.2")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.6")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.4")
} #957
#load(n,N,infothetSS, file = "sim1904ARCHRANK.Rdata")
##################################################################
## AR(1) - rank
{
  n      <-  5000; ei <- 0.2; kmax <- floor(n^0.8); kmax2 <- floor(n^0.5); 
  id <- 1; theta <- ei;  km    <-   137# length(th$k)
  for(N  in 10:1000){
    sample  <-  abs( arima.sim(n = 5000, list(ar=0.8, ma=0), rand.gen=function(n)rt(n,df=1)) )
    sample  <-  n/r(sample)
    alpha  <-   1#1/alphaestimator(sample,k1=kmax)$xi ;# print(alpha)
    th     <-   exIndex(path0=sample,k0=4:kmax2,alpha0=alpha)
    th$k  <-    1:km
    th$thet <- rep(theta, floor(km) )
    th$alpha <- rep(alpha, floor(km) )
    th$N     <-  rep(N, floor(km) )
    th$id    <-  rep(id, floor(km) )
    infothetSS <- rbind( infothetSS, th )
    print(N)
  }
  head(infothetSS)
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.2")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.6")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.4")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.2")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.6")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.4")
} 
#load(n,N,infothetSS, file = "sim1904AR0.8RANK.Rdata")
###################################################################
###################################################################
type          <- c( "Intervals","Sliding","Northrop","runs6","runs7","blocks6","blocks7","Stable26","Stable27" )
infothetSS    <-  data.frame(  "Intervals"=NULL,"Sliding"=NULL,"Northrop"=NULL,"runs6"=NULL,"runs7"=NULL,
                               "blocks6"=NULL,"blocks7"=NULL,"Stable26"=NULL,"Stable27"=NULL,"Stable28"=NULL,"Stable29"=NULL, ## Estimators
                               "k"=NULL, "thet" = NULL, "alpha" = NULL,"N" = NULL, "id" = NULL)
##################################################################
## AR(1)
{
  n      <-  5000; ei <- 0.2; kmax <- floor(n^0.8); kmax2 <- floor(n^0.5); 
  id <- 1; theta <- ei;  km    <-   137# length(th$k)
  for(N  in 1:1000){
    sample  <-  abs( arima.sim(n = 5000, list(ar=0.8, ma=0), rand.gen=function(n) rt(n,df=1)) )
    alpha  <-   1/alphaestimator(sample,k1=kmax)$xi ;# print(alpha)
    th     <-   exIndex(path0=sample,k0=4:kmax2,alpha0=alpha)
    #th$k  <-    1:km
    th$thet <- rep(theta, floor(km) )
    th$alpha <- rep(alpha, floor(km) )
    th$N     <-  rep(N, floor(km) )
    th$id    <-  rep(id, floor(km) )
    infothetSS <- rbind( infothetSS, th )
    print(N)
    
  }
  head(infothetSS)
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.2")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.6")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.4")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.2")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.6")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.4")
} #612
#load( n,N,infothetSS, file = "sim2804AR0.82.Rdata")
#load( n,N,infothetSS, file = "sim2804AR0.8.Rdata")
#load(n,N,infothetSS, file = "sim1904AR0.8.Rdata")
## ARCH(1)
{
  n      <-  5000; ei <- 0.2; kmax <- floor(n^0.8); kmax2 <- floor(n^0.5); 
  id <- 1; theta <- ei;  km    <-   137# length(th$k)
  for(N  in 11:1000){
    sample  <-  abs(ARCHm(n))
    alpha  <-   1/alphaestimator(sample,k1=kmax)$xi ;# print(alpha)
    th     <-   exIndex(path0=sample,k0=4:kmax2,alpha0=alpha)
    #th$k  <-    1:km
    th$thet <- rep(theta, floor(km) )
    th$alpha <- rep(alpha, floor(km) )
    th$N     <-  rep(N, floor(km) )
    th$id    <-  rep(id, floor(km) )
    infothetSS <- rbind( infothetSS, th )
    print(N)
  }
  head(infothetSS)
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.2")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.6")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.4")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.2")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.6")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.4")
} #957
#load( n,N,infothetSS, file = "sim2904ARCH.Rdata")
#load( n,N,infothetSS, file = "sim2804ARCH.Rdata")
#load(n,N,infothetSS, file = "sim1904ARCH.Rdata")
## AR(1) - 0.2
{
  n      <-  5000; ei <- 0.8; kmax <- floor(n^0.8); kmax2 <- floor(n^0.5); 
  id <- 1; theta <- ei;  km    <-   137# length(th$k)
  for(N  in 300:1000){
    sample  <-  abs( arima.sim(n = 5000, list(ar=0.2, ma=0), rand.gen=function(n) rt(n,df=1)) )
    alpha  <-   1/alphaestimator(sample,k1=kmax)$xi ;# print(alpha)
    th     <-   exIndex(path0=sample,k0=4:kmax2,alpha0=alpha)
    #th$k  <-    1:km
    th$thet <-  rep(theta, floor(km) )
    th$alpha <- rep(alpha, floor(km) )
    th$N     <-  rep(N, floor(km) )
    th$id    <-  rep(id, floor(km) )
    infothetSS <- rbind( infothetSS, th )
    print(N)
    
  }
  tail(infothetSS)
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.2")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.6")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAX0.4")
  #load("/Users/Buritica/Dropbox/Thèse/0/Projet_EI/SS/SSData/sim0904ARMAXRank0.2")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.6")
  #load(infothetSS,n,ei,kmax,kmax2,id,theta,file = "sim0904ARMAXRank0.4")
} #
load("/Users/Buritica/Dropbox/Thèse/0/Projet_EI/SS/SSData/sim0904ARMAX0.4")
thet <- 0.4
g1 <- ggplot(data=infothetSS, aes(y=(Stable24),x=as.factor(k)))  + geom_boxplot(fill = "#4271AE", colour = "skyblue4", alpha = 0.4,outlier.colour = "#1F3552",outlier.shape = 16, outlier.size = 0.5) + geom_hline(yintercept=thet)
g2 <- ggplot(data=infothetSS, aes(y=(Stable23),x=as.factor(k)))  + geom_boxplot(fill = "#4271AE", colour = "skyblue4", alpha = 0.4,outlier.colour = "#1F3552",outlier.shape = 16, outlier.size = 0.5) + geom_hline(yintercept=thet)
grid.arrange(g1,g2)

#load( n,N,infothetSS, file = "sim2904AR0.22NEW.Rdata")
#load( n,N,infothetSS, file = "sim2904AR0.22.Rdata")
#load(n,N,infothetSS, file = "sim2804AR0.22.Rdata")
#load(n,N,infothetSS, file = "sim2804AR0.2.Rdata")
library(grid)
library(latex2exp)
thet <- 0.2792
thet <- 0.8
{
g1 <- ggplot(data=infothetSS[ (infothetSS$k%in%c(2,4,8,16,32,64,128,250,500)) , ], aes(y=(Intervals), x=as.factor(k))) + geom_boxplot(fill = "#4271AE", colour = "skyblue4", alpha = 0.4,outlier.colour = "#1F3552", outlier.shape = 16, outlier.size = 0.5)
g2 <- ggplot(data=infothetSS[ (infothetSS$k%in%c(2,4,8,16,32,64,128,250,500)) , ], aes(y=(Sliding),   x=as.factor(k)))  + geom_boxplot(fill = "#4271AE", colour = "skyblue4", alpha = 0.4,outlier.colour = "#1F3552",outlier.shape = 16, outlier.size = 0.5)
g3 <- ggplot(data=infothetSS[ (infothetSS$k%in%c(2,4,8,16,32,64,128,250,500)) , ], aes(y=(Northrop), x=as.factor(k))) + geom_boxplot(fill = "#4271AE", colour = "skyblue4", alpha = 0.4,outlier.colour = "#1F3552",outlier.shape = 16, outlier.size = 0.5)
g4 <- ggplot(data=infothetSS[ (infothetSS$k%in%c(2,4,8,16,32,64,128,250,500)) , ], aes(y=(runs6),   x=as.factor(k))) + geom_boxplot(fill = "#4271AE", colour = "skyblue4", alpha = 0.4,outlier.colour = "#1F3552",outlier.shape = 16, outlier.size = 0.5)
g5 <- ggplot(data=infothetSS[ (infothetSS$k%in%c(2,4,8,16,32,64,128,250,500)) , ], aes(y=(runs7),   x=as.factor(k))) + geom_boxplot(fill = "#4271AE", colour = "skyblue4", alpha = 0.4,outlier.colour = "#1F3552",outlier.shape = 16, outlier.size = 0.5)
g6 <- ggplot(data=infothetSS[ (infothetSS$k%in%c(2,4,8,16,32,64,128,250,500)) , ], aes(y=(blocks6), x=as.factor(k))) + geom_boxplot(fill = "#4271AE", colour = "skyblue4", alpha = 0.4,outlier.colour = "#1F3552",outlier.shape = 16, outlier.size = 0.5)
g7 <- ggplot(data=infothetSS[ (infothetSS$k%in%c(2,4,8,16,32,64,128,250,500)) , ], aes(y=(blocks7), x=as.factor(k))) + geom_boxplot(fill = "#4271AE", colour = "skyblue4", alpha = 0.4,outlier.colour = "#1F3552",outlier.shape = 16, outlier.size = 0.5)
g8 <- ggplot(data=infothetSS[ (infothetSS$k%in%c(2,4,8,16,32,64,128,250,500)) , ], aes(y=(Stable26),x=as.factor(k))) + geom_boxplot(fill = "#4271AE", colour = "skyblue4", alpha = 0.4,outlier.colour = "#1F3552",outlier.shape = 16, outlier.size = 0.5)
g9 <- ggplot(data=infothetSS[ (infothetSS$k%in%c(2,4,8,16,32,64,128,250,500)) , ], aes(y=(Stable27),x=as.factor(k)))  + geom_boxplot(fill = "#4271AE", colour = "skyblue4", alpha = 0.4,outlier.colour = "#1F3552",outlier.shape = 16, outlier.size = 0.5)

g1 <- addingtitles(g1, TeX('$x$') , title2=(TeX('$\\widehat{\\theta}^{int}$') ) ) 
g2 <- addingtitles(g2, TeX('$r$') , title2=(TeX('$\\widehat{\\theta}^{slbl}$') ))     
g3 <- addingtitles(g3, TeX('$r$') , title2=(TeX('$\\widehat{\\theta}^{Nsl}$') ))     
g4 <- addingtitles(g4, TeX('$l$') , title2=(TeX('$\\widehat{\\theta}^{runs}$') ))
g5 <- addingtitles(g5, TeX('$r$') , title2=(TeX('$\\widehat{\\theta}^{runs,2}$') ))               
g6 <- addingtitles(g6, TeX('$r$') , title2=(TeX('$\\widehat{\\theta}^{bl}$') ))     
g7 <- addingtitles(g7, " ", title2=(TeX('$\\widehat{\\theta}^{bl,2}$') ))     
g8 <- addingtitles(g8, TeX('$r$') , title2=(TeX('$\\widehat{\\theta}^{scp}$') ))
g9 <- addingtitles(g9, " ", title2=(TeX('$\\widehat{\\theta}^{scp,2}$') ))
}
{
  g11 <- ggplot(data=infothetSS, aes(y=rev(Intervals), x=as.factor(k), group=as.factor(N))) 
  g22 <- ggplot(data=infothetSS, aes(y=rev(Sliding),   x=as.factor(k), group=as.factor(N))) 
  g33 <- ggplot(data=infothetSS, aes(y=rev(Northrop), x=as.factor(k), group=as.factor(N))) 
  g44 <- ggplot(data=infothetSS, aes(y=rev(Stable26), x=as.factor(k), group=as.factor(N))) 
  
  g11 <- g11 + geom_line(colour = alpha("skyblue4", 0.011))+ 
    geom_line(data=infothetSS[((km*10+1):(km*14)), ],colour = alpha("skyblue4", 1), aes(y=rev(Intervals), x=as.factor(k), group=as.factor(N)))
  g22 <- g22 + geom_line(colour = alpha("skyblue4", 0.011))+ 
    geom_line(data=infothetSS[((km*10+1):(km*14)), ],colour = alpha("skyblue4", 1), aes(y=rev(Sliding), x=as.factor(k), group=as.factor(N)))
  g33 <- g33 + geom_line(colour = alpha("skyblue4", 0.011))+ 
    geom_line(data=infothetSS[((km*10+1):(km*14)), ],colour = alpha("skyblue4", 1), aes(y=rev(Northrop), x=as.factor(k), group=as.factor(N)))
  g44 <- g44 + geom_line(colour = alpha("skyblue4", 0.011))+ 
    geom_line(data=infothetSS[((km*10+1):(km*14)), ],colour = alpha("skyblue4",1), aes(y=rev(Stable27), x=as.factor(k), group=as.factor(N)))
  
  g11 <- addingtitles(g11, "","Intervals")               
  g22 <- addingtitle(g22, "Sliding")     
  g33 <- addingtitle(g33, "Northrop")     
  g44 <- addingtitles(g44, "", "LD") 
}
grid.arrange(g1,g2,g3,g4,g6,g8,g5,g7,g9, ncol=3,nrow=3)

grid.arrange(g6,g2,g4,g1,g3,g8, ncol=2,nrow=3)

grid.arrange(g6,g8, ncol=2)

grid.arrange(g1,g11,ncol = 2, nrow=1)
grid.arrange(g2,g22,ncol = 2, nrow=1)
grid.arrange(g3,g33,ncol = 2, nrow=1)
grid.arrange(g4,g44,ncol = 2, nrow=1)

##
grid.arrange(g1,g2,g3,g4, ncol = 2, nrow=2)

#1F3552"
addingtitle <- function(gg,title,title2=" "){
  gg <- gg +
         geom_hline(yintercept=thet,lty=2, col="blue") +
         scale_y_continuous(name =title2 , labels = c(0,thet,1), breaks = c(0,thet,1), limits=c(0,1))+ 
         scale_x_discrete(name = " ",labels=th$k[seq(137, 1, -20)], breaks = seq(1, 137, 20) )+
         theme_bw()+              
         ggtitle(title)+
         theme(#panel.grid.major = element_line(colour = "#d3d3d3"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          plot.title = element_text(size = 14 ,face = "bold", hjust=0.5),
          text=element_text(),
          axis.title = element_text(),
          axis.text.x = element_text(colour="black", size = 9, angle = 45 ),
          axis.text.y = element_text(colour="black", size = 9),
          axis.line =   element_line(size=0.1, colour = "black"))
    
  return(gg)
} 
addingtitles <- function(gg,title,title2=" "){
  gg <- gg +
    geom_hline(yintercept=thet,lty=2, col="blue") +
    scale_y_continuous(name =" " , labels = c(0,thet,1), breaks = c(0,thet,1), limits=c(0,1))+ 
    scale_x_discrete(name = title,labels=c(2,4,8,16,32,64,128,250,500))+
    theme_bw()+              
    ggtitle(title2)+
    theme(#panel.grid.major = element_line(colour = "#d3d3d3"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 24 ,face = "bold", hjust=0),
      text=element_text(),
      axis.title = element_text(size=50),
      axis.title.y= element_text(size=24,angle=0),
      axis.title.x= element_text(size=15,angle=0),
      axis.text.x = element_text(colour="black", size = 15, angle = 15 ),
      axis.text.y = element_text(colour="black", size = 15),
      axis.line =   element_line(size=0.1, colour = "black"))
  
  return(gg)
}  
addingtitlebl <- function(gg,title){
  gg <- gg +
    geom_hline(yintercept=thet,lty=2, col="blue") +
    scale_y_continuous(name =" ", labels = c(0,thet,1), breaks = c(0,thet,1), limits=c(0,1))+ 
    scale_x_discrete(name = "k = n / block length",labels=(floor(n/th$k[seq(137, 1, -15)])), breaks = seq(1,137,15) )+
    theme_bw()+              
    ggtitle(title)+
    theme(#panel.grid.major = element_line(colour = "#d3d3d3"),
      panel.grid.major = element_blank(),
      panel.grid.minor = element_blank(),
      panel.border = element_blank(),
      panel.background = element_blank(),
      plot.title = element_text(size = 14 ,face = "bold", hjust=0.5),
      text=element_text(),
      axis.title = element_text(),
      axis.text.x = element_text(colour="black", size = 9, angle = 45 ),
      axis.text.y = element_text(colour="black", size = 9),
      axis.line =   element_line(size=0.1, colour = "black"))
  
  return(gg)
}  

## B_t = 1

ARCHm <- function(n){
  x0  <- 1
  for(i in 1:(2*n) ) x0  <- c(x0, ( exp(rnorm(1) - 0.5)*(x0[length(x0)]) + 1 ) )
  return(x0[(n+1):(2*n)])
  
}


spec1 <- ugarchspec(mean.model = list(armaOrder = c(0,0), include.mean = FALSE),
                    variance.model = list(model="GARCH",garchOrder = c(1, 0)),
                    fixed.pars = list("omega" = 0, "alpha1" = 1,
                                      "beta1" = 0),
                    distribution.model ='rnorm' )
# simulate the path
path.sgarch <- ugarchpath(spec1, n.sim=3000, n.start=1)
garch.sim(alpha=c(0.5,0.5), beta=1, n = 500, rnd = rnorm , ntrans=100)

fgspec <- garchSpec(model=list(alpha=0.999, 
                               omega=.01, beta=0.0001))
fgsim <-  garchSim(spec=fgspec, extended=TRUE, n=200)
plot.ts(fgsim@.Data[,1])


