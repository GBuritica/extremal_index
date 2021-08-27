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
load("sim2804AR0.82.Rdata") ## n,N,infothetSS
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
load("sim2904ARCH.Rdata")
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
load("sim2904AR0.22NEW.Rdata")

##############################################################
## Create plot
thet <- 0.2
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
  km   <- 137
  g11 <- ggplot(data=infothetSS, aes(y=(Intervals), x=as.factor(k), group=as.factor(N))) 
  g22 <- ggplot(data=infothetSS, aes(y=(Sliding),   x=as.factor(k), group=as.factor(N))) 
  g33 <- ggplot(data=infothetSS, aes(y=(Northrop), x=as.factor(k), group=as.factor(N))) 
  g44 <- ggplot(data=infothetSS, aes(y=(Stable26), x=as.factor(k), group=as.factor(N))) 
  
  g11 <- g11 + geom_line(colour = alpha("skyblue4", 0.011))+ 
    geom_line(data=infothetSS[((km*10+1):(km*14)), ],colour = alpha("skyblue4", 1), aes(y=(Intervals), x=as.factor(k), group=as.factor(N)))
  g22 <- g22 + geom_line(colour = alpha("skyblue4", 0.011))+ 
    geom_line(data=infothetSS[((km*10+1):(km*14)), ],colour = alpha("skyblue4", 1), aes(y=(Sliding), x=as.factor(k), group=as.factor(N)))
  g33 <- g33 + geom_line(colour = alpha("skyblue4", 0.011))+ 
    geom_line(data=infothetSS[((km*10+1):(km*14)), ],colour = alpha("skyblue4", 1), aes(y=(Northrop), x=as.factor(k), group=as.factor(N)))
  g44 <- g44 + geom_line(colour = alpha("skyblue4", 0.011))+ 
    geom_line(data=infothetSS[((km*10+1):(km*14)), ],colour = alpha("skyblue4",1), aes(y=(Stable27), x=as.factor(k), group=as.factor(N)))
  
  g11 <- addingtitle(g11, "","Intervals")               
  g22 <- addingtitle(g22, "Sliding")     
  g33 <- addingtitle(g33, "Northrop")     
  g44 <- addingtitle(g44, "", "LD") 
}
############################### Prints
## Plots results for article
grid.arrange(g6,g2,g4,g1,g3,g8, ncol=2,nrow=3)
## Plots comparison bloccs scp.
grid.arrange(g6,g8, ncol=2)
## Plots all results
grid.arrange(g1,g2,g3,g4,g6,g8,g5,g7,g9, ncol=3,nrow=3)
############################### Trajectory comparison
grid.arrange(g1,g11,ncol = 2, nrow=1)
grid.arrange(g2,g22,ncol = 2, nrow=1)
grid.arrange(g3,g33,ncol = 2, nrow=1)
grid.arrange(g8,g44,ncol = 2, nrow=1)
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
