# Extremal index: 
## Estimators impemented: Intervals
##                   Intervals estimator.
##                   Sliding estimator.
##                   Northrop estimator.
##                   SCP: spectral cluster estimator for two tunning parameters.                      


### Examples of implementation:
#### Burr noise:
set.seed(2806)
par   <- 0.6
alpha <- 2
ei    <- 1-par^alpha
sample <- abs( arima.sim(n = 3500, list(ar=par, ma=0), rand.gen=function(n) rBurr(n,params = list(b=1,g=2,s=1)) ) )
alpha  <-  1/alphaestimator(sample,k1=floor(n^0.8))$xi; print(alpha)
th    <- plot_estimators(sample,alpha0=alpha, ei0=ei)

#### Student noise:
set.seed(2806)
par   <- 0.6
alpha <- 1
ei    <- 1-par^alpha 
sample <- abs( arima.sim(n = 3500, list(ar=par, ma=0), rand.gen=function(n) rt(n,df=alpha) ) )
alpha  <-  1/alphaestimator(sample,k1=floor(n^0.8))$xi; print(alpha)
th       <- plot_estimators(sample,alpha0=alpha, ei0=ei)
