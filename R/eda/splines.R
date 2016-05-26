library(pspline)
library(energy)
library(matie)

## model association using a smoothing spline
g1 <- as.numeric(genes[90,-c(1,96)])
g2 <- as.numeric(genes[94,-c(1,96)])
sp <- sm.spline(x=g1,y=g2)

#plot the results
plot(g1,g2)
lines(sp)
lines(supsmu(g1,g2)) #Friedman's supersmoother
sp$norder

#' Calculate the residuals of a penalised smoothing spline
splineResid <- function(sp, x, y) {
  1-sum((y-predict(sp,x)[,1])^2)/sum((y-mean(y))^2)
}

sp2 <- sm.spline(g2, g1)

splineResid(sp, g1, g2)
splineResid(sp2, g2, g1)
#assymetric

#compare to other estimates of association
cor(g1,g2)^2
dcor(g1,g2)
ma(data.frame(g2,g1))$A


