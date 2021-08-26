#######################################################################
#######################################################################
#######################################################################
#######################################################################
#### Gloria Buriticá
#### Simulation Study Extremal index
#### 1. Robertetal ( Depends on a proportion k of the path )
#### 2. Northrop
source("/Users/Buritica/Dropbox/Thèse/git/Auxiliar_functions/rank_transform.R")
source("/Users/Buritica/Dropbox/Thèse/git/index_regular_variation/IndexofRV.R")
source("/Users/Buritica/Dropbox/Thèse/git/extremal_index/Estimators.R")

require(ggplot2)
require(gridExtra)
require(latex2exp)
#######################################################################
#######################################################################
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
  
} 
##### Plot
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
#################################################################################
## Additional functions : random paths
ARCHm <- function(n){
  x0  <- 1
  for(i in 1:(2*n) ) x0  <- c(x0, ( exp(rnorm(1) - 0.5)*(x0[length(x0)]) + 1 ) )
  return(x0[(n+1):(2*n)])
}
#################################################################################
## Additional functions : plots
addingtitle   <- function(gg,title,title2=" "){
  gg <- gg +
    geom_hline(yintercept=thet,lty=2, col="blue") +
    scale_y_continuous(name =title2 , labels = c(0,thet,1), breaks = c(0,thet,1), limits=c(0,1))+ 
    scale_x_discrete(name = " ",labels=infothetSS$k[seq(137,1,-20)], breaks = seq(1, 137, 20) )+
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
addingtitles  <- function(gg,title,title2=" "){
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



