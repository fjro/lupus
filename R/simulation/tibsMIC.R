require(ggplot2)
require(reshape)
## We write the following function because MIC glitches a small percentage of the time, and we do not wish to average over those trials

notNA.greater <- function(a,b){
  ind <- which(!is.na(a))
  pow <- sum(a[ind] > b)/length(ind)
  return(pow)
}

# download from  mine.jar exploredata.net
#install.packages("proxy")
library("energy")
library('minerva')
require('matie')
## This is a short R wrapper which uses system calls to call MIC
# (we wrote our own wrapper, as the one provided by the authors was difficult to use)

mymine=function(x,y){
  xx=cbind(x,y)
  write("x,y",file="test.csv")
  write(t(xx),sep=",",file="test.csv",ncol=2,append=T)
  command <- 'java -jar MINE.jar "test.csv" -allPairs'
  system(command)
  res=scan("test.csv,B=n^0.6,k=15.0x,Results.csv",what="",sep=",")
  val=as.numeric(res[11])
  return(val)
}

ma(data.frame(x,y)) 
mymine(x,y)
cor(x,y,method='spearman')
?cor
mymine=function(x,y){ mine(x,y)$MIC }

?mine

set.seed(1)

# Here we define parameters which we use to simulate the data

#500
nsim=50                         # The number of null datasets we use to estimate our rejection reject regions for an alternative with level 0.05
nsim2=50                          # Number of alternative datasets we use to estimate our power

#30
num.noise <- 30                    # The number of different noise levels used
noise <- 3                          # A constant to determine the amount of noise

#320
n=100                              # Number of data points per simulation

val.cor=val.dcor=val.mine=val.A=val.spear=rep(NA,nsim)              # Vectors holding the null "correlations" (for pearson, dcor and mic respectively) for each of the nsim null datasets at a given noise level

val.cor2=val.dcor2=val.mine2=val.A2=val.spear2= rep(NA,nsim2)              # Vectors holding the alternative "correlations" (for pearson, dcor and mic respectively) for each of the nsim2 alternative datasets at a given noise level

power.cor=power.dcor=power.mine=power.A=power.spear= array(NA, c(8,num.noise))                # Arrays holding the estimated power for each of the "correlation" types, for each data type (linear, parabolic, etc...) with each noise level

## We loop through the noise level and functional form; each time we estimate a null distribution based on the marginals of the data, and then use that null distribution to estimate power

## We use a uniformly distributed x, because in the original paper the authors used the same

functions = list("linear"=3,"quadratic"=4)
functions$linear
functions[2]
length(functions)

for(l in 1:num.noise)
{
  for(typ in 1:8)
  {
    ## This next loop simulates data under the null with the correct marginals (x is uniform, and y is a function of a uniform with gaussian noise)
    
    for(ii in 1:nsim){
      cat(sprintf("Noise %s; Type %s; sim = %s\"\n",  l, typ,ii))
      x=runif(n)
      
      if(typ==1){
        y=x+ noise *(l/num.noise)* rnorm(n)
      }
      #parabolic+noise
      if(typ==2){
        y=4*(x-.5)^2+  noise * (l/num.noise) * rnorm(n)
      }
      #cubic+noise
      if(typ==3){
        y=128*(x-1/3)^3-48*(x-1/3)^3-12*(x-1/3)+ 10* noise  * (l/num.noise) *rnorm(n)
      }
      #sin+noise
      if(typ==4){
        y=sin(4*pi*x) + 2*noise * (l/num.noise) *rnorm(n)
      }
      #their sine + noise
      if(typ==5){
        y=sin(16*pi*x) + noise * (l/num.noise) *rnorm(n)
      }
      #x^(1/4) + noise
      if(typ==6){
        y=x^(1/4) + noise * (l/num.noise) *rnorm(n)
      }
      #circle
      if(typ==7){
        y=(2*rbinom(n,1,0.5)-1) * (sqrt(1 - (2*x - 1)^2)) + noise/4*l/num.noise *rnorm(n)
      }
      #step function
      if(typ==8){
        y = (x > 0.5) + noise*5*l/num.noise *rnorm(n)
      }
      
      x <- runif(n)                       # We resimulate x so that we have the null scenario
      
      val.cor[ii]=(cor(x,y))^2            # Calculate the correlation
      val.dcor[ii]=dcor(x,y)              # Calculate dcor
      val.mine[ii]=mymine(x,y)            # Calculate mic
      val.A[ii]=ma(data.frame(x,y))$A             # Calculate A
      val.spear[ii]=(cor(x,y,method='spearman'))^2             # Calculate mic
    }
    
    val.mine <- val.mine[which(!is.na(val.mine))]                 # we remove the mic trials which glitch
    
    ## Next we calculate our 3 rejection cutoffs
    
    cut.cor=quantile(val.cor,.95)
    cut.dcor=quantile(val.dcor,.95)
    cut.mine=quantile(val.mine,.95)
    cut.A=quantile(val.A,.95)
    cut.spear=quantile(val.spear,.95)
    ## Next we simulate the data again, this time under the alternative
    
    for(ii in 1:nsim2){
      cat(sprintf("Noise 2 %s; Type %s; sim = %s\"\n",  l, typ,ii))
      x=runif(n)
      
      #lin+noise
      if(typ==1){
        y=x+ noise * (l/num.noise) *rnorm(n)
      }
      #parabolic+noise
      if(typ==2){
        y=4*(x-.5)^2+  noise * (l/num.noise)*rnorm(n)
      }
      #cubic+noise
      if(typ==3){
        y=128*(x-1/3)^3-48*(x-1/3)^3-12*(x-1/3)+10* noise * (l/num.noise) *rnorm(n)
      }
      #sin+noise
      if(typ==4){
        y=sin(4*pi*x) + 2*noise * (l/num.noise) *rnorm(n)
      }
      #their sine + noise
      if(typ==5){
        y=sin(16*pi*x) + noise * (l/num.noise) *rnorm(n)
      }
      #x^(1/4) + noise
      if(typ==6){
        y=x^(1/4) + noise * (l/num.noise) *rnorm(n)
      }
      #circle
      if(typ==7){
        y=(2*rbinom(n,1,0.5)-1) * (sqrt(1 - (2*x - 1)^2)) + noise/4*l/num.noise *rnorm(n)
      }
      #step function
      if(typ==8){
        y = (x > 0.5) + noise*5*l/num.noise *rnorm(n)
      }
      
      ## We again calculate our "correlations"
      
      val.cor2[ii]=(cor(x,y))^2
      val.dcor2[ii]=dcor(x,y)
      val.mine2[ii]=mymine(x,y)
      val.A2[ii]=ma(data.frame(x,y))$A
      val.spear2[ii]=(cor(x,y,method='spearman'))^2
    }
    
    ## Now we estimate the power as the number of alternative statistics exceeding our estimated cutoffs
    
    power.cor[typ,l] <- sum(val.cor2 > cut.cor)/nsim2
    power.dcor[typ,l] <- sum(val.dcor2 > cut.dcor)/nsim2
    power.mine[typ,l] <- notNA.greater(val.mine2, cut.mine)
    power.A[typ,l] <- sum(val.A2 > cut.A)/nsim2
    power.spear[typ,l] <- sum(val.spear2 > cut.spear)/nsim2
  }
}

#save.image()

## The rest of the code is for plotting the image
#?pdf
#pdf("power.pdf")
#dev.off()
par(mfrow = c(4,2), cex = 0.45)
plot((1:30)/10, power.cor[1,], ylim = c(0,1), main = "Linear", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[1,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[1,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","A"), pch = c(1,2,3), col = c("black","green","red"))

plot((1:30)/10, power.cor[2,], ylim = c(0,1), main = "Quadratic", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[2,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[2,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","A"), pch = c(1,2,3), col = c("black","green","red"))

plot((1:30)/10, power.cor[3,], ylim = c(0,1), main = "Cubic", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[3,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[3,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","A"), pch = c(1,2,3), col = c("black","green","red"))

plot((1:30)/10, power.cor[5,], ylim = c(0,1), main = "Sine: period 1/8", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[5,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[5,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","A"), pch = c(1,2,3), col = c("black","green","red"))

plot((1:30)/10, power.cor[4,], ylim = c(0,1), main = "Sine: period 1/2", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[4,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[4,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","A"), pch = c(1,2,3), col = c("black","green","red"))

plot((1:30)/10, power.cor[6,], ylim = c(0,1), main = "X^(1/4)", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[6,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[6,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","A"), pch = c(1,2,3), col = c("black","green","red"))

plot((1:30)/10, power.cor[7,], ylim = c(0,1), main = "Circle", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[7,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[7,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","A"), pch = c(1,2,3), col = c("black","green","red"))

plot((1:30)/10, power.cor[8,], ylim = c(0,1), main = "Step function", xlab = "Noise Level", ylab = "Power", pch = 1, col = "black", type = 'b')
points((1:30)/10, power.dcor[8,], pch = 2, col = "green", type = 'b')
points((1:30)/10, power.mine[8,], pch = 3, col = "red", type = 'b')
legend("topright",c("cor","dcor","A"), pch = c(1,2,3), col = c("black","green","red"))

#ggplot
powerdf= data.frame(power.cor)
powerdf = t(powerdf)
#powerdf = cbind(powerdf,rep('Pearson',30))
colnames(powerdf) = c('Linear','Quadratic','Cubic','Sine: period 1/8',
                          'Sine: period 1/4', 'X^(1/4)','Circle','Step Function')

require(ggplot2)
ff = melt(powerdf)
ff = cbind(ff,rep('Pearson',30))
ff = cbind(ff,(1:30)/10)
colnames(ff) = c('var','Form','Power','Statistic','Noise')

allCors = function(df, x, stat)
{
  tf= data.frame(x)
  tf = t(tf)
  colnames(tf) = c('Linear','Quadratic','Cubic','Sine: period 1/8',
                        'Sine: period 1/4', 'X^(1/4)','Circle','Step Function')
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
write.csv(ff,'power1000.csv')
ggplot(ff, aes(x=Noise, y=Power,group=Statistic,colour=Statistic)) +geom_line() + facet_wrap(~ Form,scales="free_y")



#plot the function forms
n=1000 
noise =1
num.noise=1
l=0.1
#define the functions
linear = function(x, noise, noiseLevel, num.noise){x+ noise *(noiseLevel/num.noise)* rnorm(n)}
parabolic = function(x, noise, noiseLevel, num.noise){y=4*(x-.5)^2+  noise * (noiseLevel/num.noise) * rnorm(n)}
cubic = function(x, noise, noiseLevel, num.noise){y=128*(x-1/3)^3-48*(x-1/3)^3-12*(x-1/3)+10* noise  * (noiseLevel/num.noise) *rnorm(n)}


sin1 = function(x, noise, noiseLevel, num.noise){y=sin(4*pi*x) + 2*noise * (noiseLevel/num.noise) *rnorm(n)}
sin2= function(x, noise, noiseLevel, num.noise){y=sin(16*pi*x) + noise * (noiseLevel/num.noise) *rnorm(n)}#their sine + noise
qroot=function(x, noise, noiseLevel, num.noise){y=x^(1/4) + noise * (noiseLevel/num.noise) *rnorm(n)} #x^(1/4) + noise
circle=function(x, noise, noiseLevel, num.noise){y=(2*rbinom(n,1,0.5)-1) * (sqrt(1 - (2*x - 1)^2)) + noise/4*noiseLevel/num.noise *rnorm(n)}
step= function(x, noise, noiseLevel, num.noise){y = (x > 0.5) + noise*5*noiseLevel/num.noise *rnorm(n)}

spike = function(x, noise, noiseLevel, num.noise) { ifelse (x < 0.05, 20,ifelse (x < 0.1, -18*x + 1.9 ,-x/9 + 1/9 ))+ noise*10*noiseLevel/num.noise *rnorm(n)}
#sigmoid = function(x, noise, noiseLevel, num.noise) { ifelse (x < 0.49, 0,ifelse (x < 0.51, 50*(x-0.5) + 0.5 ,1))+ noise*noiseLevel/num.noise *rnorm(n)}
#sigmoid = function(x, noise, noiseLevel, num.noise) { x/sqrt(1+x^2)+ noise*noiseLevel/num.noise *rnorm(n)}
sigmoid= function(x, noise, noiseLevel, num.noise){y= 0.2*sin(4*(2*x -1)) + 1.1*(2*x-1) + noise * (noiseLevel/num.noise) *rnorm(n)}#their sine + noise

xShaped = function(x, noise, noiseLevel, num.noise){y=((4*(x-.5)^2+  noise * (noiseLevel/num.noise) * rnorm(n)) * sample( c(-1,1), size=n, replace=T ) )}

linearPeriodic= function(x, noise, noiseLevel, num.noise){y= sin(10*pi*x) + x + noise * (noiseLevel/num.noise) *rnorm(n)}
varyingFreq= function(x, noise, noiseLevel, num.noise){y= sin(5*pi*x*(1+x)) + x + noise * (noiseLevel/num.noise) *rnorm(n)}
logE = function(x, noise, noiseLevel, num.noise){y=log(x) + (noise * (noiseLevel/num.noise) * rnorm(n))}
plot(x, log(x) + (noise * (noiseLevel/num.noise) * rnorm(n)))

exponential2 = function(x, noise, noiseLevel, num.noise){y=exp(5*x) + 10 * noise * (noiseLevel/num.noise) * rnorm(n)}
exponential10 = function(x, noise, noiseLevel, num.noise){y=10^(5*x) + 10* noise * (noiseLevel/num.noise) * rnorm(n)}
cubicY = function(x, noise, noiseLevel, num.noise){y=48*(128*(x-1/3)^3-48*(x-1/3)^3-12*(x-1/3))+10* noise *10 * (noiseLevel/num.noise) *rnorm(n)}
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
appendForm = function(df, fun, type, noise, noiseLevel, num.noise)
{
  y=fun(x, noise, noiseLevel, num.noise)
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
#hist(x)
#x=rnorm
#x = rbeta(n,2,5)
#x=rexp(n)
df=data.frame(x,linear(x,0.1,0.1,1),rep('Linear',n), rep(0.1,n))
colnames(df) = c('x', 'y', 'Form','Noise')
 
#build it up
df = appendForm(df,linear,'Linear',0.1,1,1)
df = appendForm(df,linear,'Linear',0.1,3,1)

df = appendForm(df,parabolic,'Quadratic',0.1,0.1,1)
df = appendForm(df,parabolic,'Quadratic',0.1,1,1)
df = appendForm(df,parabolic,'Quadratic',0.1,3,1)

df = appendForm(df,cubic,'Cubic',0.1,0.1,1)
df = appendForm(df,cubic,'Cubic',0.1,1,1)
df = appendForm(df,cubic,'Cubic',0.1,3,1)

df = appendForm(df,cubicY,'Cubic Y Stretched',0.1,0.1,1)
df = appendForm(df,cubicY,'Cubic Y Stretched',0.1,1,1)
df = appendForm(df,cubicY,'Cubic Y Stretched',0.1,3,1)

df = appendForm(df,qroot,'X^(1/4)',0.1,0.1,1)
df = appendForm(df,qroot,'X^(1/4)',0.1,1,1)
df = appendForm(df,qroot,'X^(1/4)',0.1,3,1)

df = appendForm(df,exponential2,'Exponential',0.1,0.1,1)
df = appendForm(df,exponential2,'Exponential',0.1,1,1)
df = appendForm(df,exponential2,'Exponential',0.1,3,1)

df = appendForm(df,logE,'Log',0.1,0.1,1)
df = appendForm(df,logE,'Log',0.1,1,1)
df = appendForm(df,logE,'Log',0.1,3,1)

df = appendForm(df,sigmoid,'Sigmoid',0.1,0.1,1)
df = appendForm(df,sigmoid,'Sigmoid',0.1,1,1)
df = appendForm(df,sigmoid,'Sigmoid',0.1,3,1)

df = appendForm(df,step,'Step Function',0.1,0.1,1)
df = appendForm(df,step,'Step Function',0.1,1,1)
df = appendForm(df,step,'Step Function',0.1,3,1)

df = appendForm(df,spike,'Spike',0.1,0.1,1)
df = appendForm(df,spike,'Spike',0.1,1,1)
df = appendForm(df,spike,'Spike',0.1,3,1)

df = appendForm(df,sin1,'Sine: period 1/8',0.1,0.1,1)
df = appendForm(df,sin1,'Sine: period 1/8',0.1,1,1)
df = appendForm(df,sin1,'Sine: period 1/8',0.1,3,1)

df = appendForm(df,sin2,'Sine: period 1/4',0.1,0.1,1)
df = appendForm(df,sin2,'Sine: period 1/4',0.1,1,1)
df = appendForm(df,sin2,'Sine: period 1/4',0.1,3,1)

df = appendForm(df,medFreq,'Linear+Periodic',0.1,0.1,1)
df = appendForm(df,medFreq,'Linear+Periodic',0.1,1,1)
df = appendForm(df,medFreq,'Linear+Periodic',0.1,3,1)

df = appendForm(df,varyingFreq,'Varying Frequency',0.1,0.1,1)
df = appendForm(df,varyingFreq,'Varying Frequency',0.1,1,1)
df = appendForm(df,varyingFreq,'Varying Frequency',0.1,3,1)


df = appendForm(df,circle,'Circle',0.1,0.1,1)
df = appendForm(df,circle,'Circle',0.1,1,1)
df = appendForm(df,circle,'Circle',0.1,3,1)

df = appendForm(df,xShaped,'X',0.1,0.1,1)
df = appendForm(df,xShaped,'X',0.1,1,1)
df = appendForm(df,xShaped,'X',0.1,3,1)

plot(x, log(x))

ggplot(df, aes(x=x, y=y, colour=Noise)) +geom_point(alpha=0.2) + facet_wrap(~ Form,scales="free_y", ncol=4) + theme(legend.position="bottom")  

##understand alg
for(l in 1:num.noise){
  for(typ in 1:8){
    ## This next loop simulates data under the null with the correct marginals (x is uniform, and y is a function of a uniform with gaussian noise)
    
    for(ii in 1:nsim){
    }
x=runif(n)
y=x+ noise *(l/num.noise)* rnorm(n)
x <- runif(n)                       # We resimulate x so that we have the null scenario
val.cor[ii]=(cor(x,y))^2   # Calculate the correlation
cut.cor=quantile(val.cor,.95)

## Next we simulate the data again, this time under the alternative
xa=runif(n)
ya=xa+ noise * (l/num.noise) *rnorm(n)
val.cor2[ii]=(cor(xa,ya))^2
  

## Now we estimate the power as the number of alternative statistics exceeding our estimated cutoffs

power.cor[typ,l] <- sum(val.cor2 > cut.cor)/nsim2
power.dcor[typ,l] <- sum(val.dcor2 > cut.dcor)/nsim2
power.mine[typ,l] <- notNA.greater(val.mine2, cut.mine)
power.A[typ,l] <- sum(val.A2 > cut.A)/nsim2
power.spear[typ,l] <- sum(val.spear2 > cut.spear)/nsim2
