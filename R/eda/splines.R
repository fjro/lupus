install.packages("pspline")
require(pspline)
example(smooth.Pspline)
## smoother line is given by
xx <- seq(4, 25, length=100)
lines(xx, predict(sm.spline(speed, dist, df=5), xx), col = "red")
## add plots of derivatives
lines(xx, 10*predict(sm.spline(speed, dist), xx, 1), col = "blue")
lines(xx, 100*predict(sm.spline(speed, dist), xx, 2), col = "green")
?sm.spline
g1 = as.numeric(genes[90,-c(1,96)])
g2 = as.numeric(genes[94,-c(1,96)])
plot(g1,g2)
sp = sm.spline(x=g1,y=g2)
lines(sp)
sp$norder
length(sp$ysmth[,1])
length(g2)
splineResid = function(sp,x,y)
{
  1-sum((y-predict(sp,x)[,1])^2)/sum((y-mean(y))^2)
}
splineResid(sp,g1,g2)
sp2 = sm.spline(g2,g1,norder=2)
splineResid(sp2,g2,g1)
#1-(SUM[ (yi-pi)^2 ]/ SUM[ (yi - mean(y))^2])
1-(sum(g2-predict(sp,g1)[,1])^2/sum(g2-mean(g2))^2)
cor(g1,g2)^2
require(energy)
dcor(g1,g2)
require(matie)
ma(data.frame(g2,g1))$A
lines(smooth.spline(g1,g2))
lines(supsmu(g1,g2))
?supersmoo
