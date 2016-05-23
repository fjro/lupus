require(ggplot2)
require(reshape)
library("energy")
library('minerva')
require('matie')
require(pspline)
#install.packages('HHG')

## We write the following function because MIC glitches a small percentage of the time, and we do not wish to average over those trials
hist(rbeta(n,2,8))
notNA.greater <- function(a,b){
  ind <- which(!is.na(a))
  pow <- sum(a[ind] > b)/length(ind)
  return(pow)
}

#pensalised spline
splineR2 = function(x,y)
{
  1-sum((y-predict(sm.spline(x=g1,y=g2),x)[,1])^2)/sum((y-mean(y))^2)
}

x = runif(320)
plot(x, exp(x^2))
plot(x, log(x/(1-x)))
(1 + eβ(0.5−x))−1
(1 + eβ(0.5−x))−1

plot(x, (1 + exp(10*(0.5 - x)))^-1 )
plot(x, sigmoid(x, 0.1, 1, 10))
plot(x, exponential2(x, 0, 1, 1)-1)
plot(x, cubic(x, 1, 1, 1))
plot(x, spike(x, 1, 3, 30))

#define the functions
#define the functions
linear = function(x, noise, noiseLevel, num.noise){x+ noise *(noiseLevel/num.noise)* rnorm(n)}
quadratic = function(x, noise, noiseLevel, num.noise){y=4*(x-.5)^2+  noise * (noiseLevel/num.noise) * rnorm(n)}
cubic = function(x, noise, noiseLevel, num.noise){y=128*(x-1/3)^3-48*(x-1/3)^3-12*(x-1/3)+10* noise  * (noiseLevel/num.noise) *rnorm(n)}

qroot=function(x, noise, noiseLevel, num.noise){y=x^(1/4) + noise * (noiseLevel/num.noise) *rnorm(n)} #x^(1/4) + noise
#exponential2 = function(x, noise, noiseLevel, num.noise){y=exp(5*x) + 10 * noise * (noiseLevel/num.noise) * rnorm(n)}
#exponential2 = function(x, noise, noiseLevel, num.noise){y=exp(x) + noise * (noiseLevel/num.noise) * rnorm(n)}
exponential2 = function(x, noise, noiseLevel, num.noise){y=exp(x^2) + (1.5 *noise * (noiseLevel/num.noise) * rnorm(n))}

logE = function(x, noise, noiseLevel, num.noise){y=log(x) + 2 * (noise * (noiseLevel/num.noise) * rnorm(n))}

sigmoid= function(x, noise, noiseLevel, num.noise){((1 + exp(10*(0.5 - x)))^-1) +( noise * (noiseLevel/num.noise) * rnorm(n))}#their sine + noise
#sigmoid= function(x, noise, noiseLevel, num.noise){y= 0.2*sin(4*(2*x -1)) + 1.1*(2*x-1) + 2* noise * (noiseLevel/num.noise) *rnorm(n)}#their sine + noise

step= function(x, noise, noiseLevel, num.noise){y = (x > 0.5) + noise*5*noiseLevel/num.noise *rnorm(n)}
#spike = function(x, noise, noiseLevel, num.noise) { ifelse (x < 0.05, 20,ifelse (x < 0.1, -18*x + 1.9 ,-x/9 + 1/9 ))+ noise*10*noiseLevel/num.noise *rnorm(n)}
spike = function(x, noise, noiseLevel, num.noise) { ifelse (x < 0.05, 4,ifelse (x < 0.1, -18*x + 1.9 ,-x/9 + 1/9 ))+ noise*5*noiseLevel/num.noise *rnorm(n)}

sin1 = function(x, noise, noiseLevel, num.noise){y=sin(4*pi*x) + 2*noise * (noiseLevel/num.noise) *rnorm(n)}
sin2= function(x, noise, noiseLevel, num.noise){y=sin(16*pi*x) + 2*noise * (noiseLevel/num.noise) *rnorm(n)}#their sine + noise
linearPeriodic= function(x, noise, noiseLevel, num.noise){y= sin(10*pi*x) + x + 4*noise * (noiseLevel/num.noise) *rnorm(n)}
varyingFreq= function(x, noise, noiseLevel, num.noise){y= sin(5*pi*x*(1+x)) + x + 3*noise * (noiseLevel/num.noise) *rnorm(n)}

circle=function(x, noise, noiseLevel, num.noise){y=(2*rbinom(n,1,0.5)-1) * (sqrt(1 - (2*x - 1)^2)) + noise/4*noiseLevel/num.noise *rnorm(n)}
xShaped = function(x, noise, noiseLevel, num.noise){y=((4*(x-.5)^2 + (noiseLevel/num.noise) * rnorm(n)) * sample( c(-1,1), size=n, replace=T ) )}
#xShaped = function(x, noise, noiseLevel, num.noise){y=((4*(x-.5)^2+  noise * (noiseLevel/num.noise) * rnorm(n)) * sample( c(-1,1), size=n, replace=T ) )}

functions = list("linear"=linear,"Quadratic"=parabolic,"Cubic"=cubic, "Fourth Root" = qroot, "Exponential" = exponential2, 
                 "Natural Log" = logE, "Sigmoid" = sigmoid, "Step"=step, "Spike" = spike,
                 "Sine: Low"= sin1, "Sine: High" = sin2, "Linear+Periodic" = linearPeriodic, "Varying Frequency" = varyingFreq,
                 "Circle" = circle, "X" = xShaped)
#functions = list("linear"=linear,"Quadratic"=parabolic)

set.seed(1)

# Here we define parameters which we use to simulate the data

#500
# The number of null datasets we use to estimate our rejection reject regions for an alternative with level 0.05
nsim=50
                         # Number of alternative datasets we use to estimate our power

#30
num.noise <- 30                    # The number of different noise levels used
noise <- 3                          # A constant to determine the amount of noise

#320
# Number of data points per simulation
n=320

# Vectors holding the null "correlations" (for pearson, dcor and mic respectively) for each of the nsim null datasets at a given noise level
val.cor=val.dcor=val.mine=val.A=val.spear=rep(NA,nsim)        
# Vectors holding the alternative "correlations" (for pearson, dcor and mic respectively) for each of the nsim2 alternative datasets at a given noise level
val.cor2=val.dcor2=val.mine2=val.A2=val.spear2= rep(NA,nsim)        

power.cor=power.dcor=power.mine=power.A=power.spear= array(NA, c(length(functions),num.noise))                # Arrays holding the estimated power for each of the "correlation" types, for each data type (linear, parabolic, etc...) with each noise level

## We loop through the noise level and functional form; each time we estimate a null distribution based on the marginals of the data, and then use that null distribution to estimate power

## We use a uniformly distributed x, because in the original paper the authors used the same

for(l in 1:num.noise)
{
  for(typ in 1:length(functions))
  {
    ## This next loop simulates data under the null with the correct marginals (x is uniform, and y is a function of a uniform with gaussian noise)
    
    cat(sprintf("Noise Null %s; Type %s\n",  l, typ))
    for(ii in 1:nsim)
    {
      x=runif(n)
      #x = rbeta(n,2,5)
      y=functions[[typ]](x, noise, l, num.noise)
      
      # We resimulate x so that we have the null scenario
      x <- runif(n) 
      #x = rbeta(n,2,5)
      
      # calculate the association
      val.cor[ii]=(cor(x,y))^2            # Calculate the correlation
      val.dcor[ii]=dcor(x,y)              # Calculate dcor
      val.mine[ii]=mine(x,y)$MIC            # Calculate mic
      val.A[ii]=ma(data.frame(x,y))$A             # Calculate A
      val.spear[ii]=(cor(x,y,method='spearman'))^2             # Calculate spearman
    }
    
    val.mine <- val.mine[which(!is.na(val.mine))]                 # we remove the mic trials which glitch
    
    ## Next we calculate our rejection cutoffs
    cut.cor=quantile(val.cor,.95)
    cut.dcor=quantile(val.dcor,.95)
    cut.mine=quantile(val.mine,.95)
    cut.A=quantile(val.A,.95)
    cut.spear=quantile(val.spear,.95)
    #cut.spline=quantile(val.spline,.95)
    
    ## Next we simulate the data again, this time under the alternative
    cat(sprintf("Noise Alt %s; Type %s;\n",  l, typ))
    for(ii in 1:nsim)
    {
      x=runif(n)
      #x = rbeta(n,2,5)
      
      y=functions[[typ]](x, noise, l, num.noise)
      
      ## We again calculate our "correlations"
      val.cor2[ii]=(cor(x,y))^2
      val.dcor2[ii]=dcor(x,y)
      val.mine2[ii]=mine(x,y)$MIC
      val.A2[ii]=ma(data.frame(x,y))$A
      val.spear2[ii]=(cor(x,y,method='spearman'))^2
      #val.spline2[ii]=splineR2(x,y)
    }
    
    ## Now we estimate the power as the number of alternative statistics exceeding our estimated cutoffs
    power.cor[typ,l] <- sum(val.cor2 > cut.cor)/nsim
    power.dcor[typ,l] <- sum(val.dcor2 > cut.dcor)/nsim
    power.mine[typ,l] <- notNA.greater(val.mine2, cut.mine)
    power.A[typ,l] <- sum(val.A2 > cut.A)/nsim
    power.spear[typ,l] <- sum(val.spear2 > cut.spear)/nsim
  }
}

#Plot the results
powerdf= data.frame(power.cor)
powerdf = t(powerdf)
colnames(powerdf) = names(functions)

ff = melt(powerdf)
ff = cbind(ff,rep('Pearson',30))
ff = cbind(ff,(1:30)/10)
colnames(ff) = c('var','Form','Power','Statistic','Noise')


allCors = function(df, x, stat)
{
  tf= data.frame(x)
  tf = t(tf)
  colnames(tf) = names(functions)
  tf = melt(tf)
  tf = cbind(tf,rep(stat,30))
  tf = cbind(tf,(1:30)/10)
  colnames(tf) = c('var','Form','Power','Statistic','Noise')
  df=rbind(df,tf)
  df
}

ff = allCors(ff,power.dcor,'dcor')
ff = allCors(ff,power.mine,'MIC')
ff = allCors(ff,power.A,'A')
ff = allCors(ff,power.spear,'Spearman')
#ff = allCors(ff,power.spline,'Spline')

#reorder based on form
ff$Form <- reorder(ff$Form, new.order=names(functions))
off = ff[ order(ff$Form), ]
ggplot(off, aes(x=Noise, y=Power,group=Statistic,colour=Statistic)) +geom_line(size=1.1) + facet_wrap(~ Form, ncol=3) + theme(legend.position="bottom")  



write.csv(ff,'powerNoise320Beta25.csv')


#plot the function forms
n=320
noise=1
num.noise=1
l=0.1

#plot(x,10*(128*(x-1/3)^3-48*(x-1/3)^3-12*(x-1/3))+10* noise *10 * (noiseLevel/num.noise) *rnorm(n))
#plot(mm(x), 41*(4*x^3 + x^2 -4*x))
#cos1 = function(x, noise, noiseLevel, num.noise){y=runif(n) + 2*noise * (noiseLevel/num.noise) *rnorm(n)}
#cos2= function(x, noise, noiseLevel, num.noise){y= 0.2*sin(4*(2*x -1)) + 1.1*(2*x-1) + noise * (noiseLevel/num.noise) *rnorm(n)}#their sine + noise
sample(c(-1,1),1) *
  plot(x,x/sqrt(1+x^4))
plot(x,tanh(x))
plot(x, x^0.25)
plot(x,2/pi * atan(pi/2 *x))
#builds up the frame
appendForm = function(df, fun, type, noiseLevel)
{
  y=fun(x, 3, noiseLevel, 30)
  ndf=data.frame(x,y, rep(type,n),rep(noiseLevel,n))
  colnames(ndf)[3:4]=c('Form','Noise')
  df = rbind(df, ndf)
  df
}

#hist(rbeta(10000,15,2))
hist(rbeta(1000,2,5))
#hist(rbeta(10000,5,5))

#define the frame
x=runif(n)
require(ggplot2)
#hist(x)
#x=rnorm
#x = rbeta(n,2,5)
#x=rexp(n)
#appendForm = function(df, fun, type, noiseLevel)
n=320
df=data.frame(x,linear(x,3,0.1,30),rep('Linear',n), rep(0.1,n))
colnames(df) = c('x', 'y', 'Form','Noise')

#build it up
df = appendForm(df,linear,'Linear',1)
df = appendForm(df,linear,'Linear',10)
df = appendForm(df,linear,'Linear',20)
df = appendForm(df,linear,'Linear',30)

df = appendForm(df,parabolic,'Quadratic',0.1)
df = appendForm(df,parabolic,'Quadratic',1)
df = appendForm(df,parabolic,'Quadratic',10)
df = appendForm(df,parabolic,'Quadratic',20)
df = appendForm(df,parabolic,'Quadratic',30)

df = appendForm(df,cubic,'Cubic',0.1)
df = appendForm(df,cubic,'Cubic',1)
df = appendForm(df,cubic,'Cubic',10)
df = appendForm(df,cubic,'Cubic',20)
df = appendForm(df,cubic,'Cubic',30)

df = appendForm(df,qroot,'X^(1/4)',0.1)
df = appendForm(df,qroot,'X^(1/4)',1)
df = appendForm(df,qroot,'X^(1/4)',10)
df = appendForm(df,qroot,'X^(1/4)',20)
df = appendForm(df,qroot,'X^(1/4)',30)

df = appendForm(df,exponential2,'Exponential',0.1)
df = appendForm(df,exponential2,'Exponential',1)
df = appendForm(df,exponential2,'Exponential',10)
df = appendForm(df,exponential2,'Exponential',20)
df = appendForm(df,exponential2,'Exponential',30)

df = appendForm(df,logE,'Log',0.1)
df = appendForm(df,logE,'Log',1)
df = appendForm(df,logE,'Log',10)
df = appendForm(df,logE,'Log',20)
df = appendForm(df,logE,'Log',30)

df = appendForm(df,sigmoid,'Sigmoid',0.1)
df = appendForm(df,sigmoid,'Sigmoid',1)
df = appendForm(df,sigmoid,'Sigmoid',10)
df = appendForm(df,sigmoid,'Sigmoid',20)
df = appendForm(df,sigmoid,'Sigmoid',30)

df = appendForm(df,step,'Step Function',0.1)
df = appendForm(df,step,'Step Function',1)
df = appendForm(df,step,'Step Function',10)
df = appendForm(df,step,'Step Function',20)
df = appendForm(df,step,'Step Function',30)

df = appendForm(df,spike,'Spike',0.1)
df = appendForm(df,spike,'Spike',1)
df = appendForm(df,spike,'Spike',10)
df = appendForm(df,spike,'Spike',20)
df = appendForm(df,spike,'Spike',30)

df = appendForm(df,sin1,'Sine: period 1/8',0.1)
df = appendForm(df,sin1,'Sine: period 1/8',1)
df = appendForm(df,sin1,'Sine: period 1/8',10)
df = appendForm(df,sin1,'Sine: period 1/8',20)
df = appendForm(df,sin1,'Sine: period 1/8',30)

df = appendForm(df,sin2,'Sine: period 1/4',0.1)
df = appendForm(df,sin2,'Sine: period 1/4',1)
df = appendForm(df,sin2,'Sine: period 1/4',10)
df = appendForm(df,sin2,'Sine: period 1/4',20)
df = appendForm(df,sin2,'Sine: period 1/4',30)

df = appendForm(df,linearPeriodic,'Linear+Periodic',0.1)
df = appendForm(df,linearPeriodic,'Linear+Periodic',1)
df = appendForm(df,linearPeriodic,'Linear+Periodic',10)
df = appendForm(df,linearPeriodic,'Linear+Periodic',20)
df = appendForm(df,linearPeriodic,'Linear+Periodic',30)

df = appendForm(df,varyingFreq,'Varying Frequency',0.1)
df = appendForm(df,varyingFreq,'Varying Frequency',1)
df = appendForm(df,varyingFreq,'Varying Frequency',10)
df = appendForm(df,varyingFreq,'Varying Frequency',20)
df = appendForm(df,varyingFreq,'Varying Frequency',30)

df = appendForm(df,circle,'Circle',0.1)
df = appendForm(df,circle,'Circle',1)
df = appendForm(df,circle,'Circle',10)
df = appendForm(df,circle,'Circle',20)
df = appendForm(df,circle,'Circle',30)

df = appendForm(df,xShaped,'X',0.1)
df = appendForm(df,xShaped,'X',1)
df = appendForm(df,xShaped,'X',10)
df = appendForm(df,xShaped,'X',20)
df = appendForm(df,xShaped,'X',30)
df$Noise = df$Noise/10
ggplot(df, aes(x=x, y=y, colour=Noise)) +geom_point(alpha=0.2, size=1) + facet_wrap(~ Form,scales="free_y", ncol=3) + theme_bw() + theme(legend.position = "bottom")

ggplot(df[df$Noise <=0.01,], aes(x=x, y=y)) +geom_point() + facet_wrap(~ Form,scales="free_y", ncol=3 )+ theme(legend.position = "bottom")

