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



ma(data.frame(x,y)) 
mymine(x,y)
rA = function(y,x) { ma(cbind(x,lm(y~x)$residuals))$A }
rDcor = function(y,x) { dcor(x,lm(y~x)$residuals) }
rMIC = function(y,x){ mine(x,lm(y~x)$residuals)$MIC}
mymine=function(x,y){ mine(x,y)$MIC }

set.seed(1)

# Here we define parameters which we use to simulate the data

#500
nsim=500                          # The number of null datasets we use to estimate our rejection reject regions for an alternative with level 0.05
nsim2=500                           # Number of alternative datasets we use to estimate our power

#30
num.noise <- 30                    # The number of different noise levels used
noise <- 3                          # A constant to determine the amount of noise

n=320                               # Number of data points per simulation

val.rdcor=val.rmine=val.rA=rep(NA,nsim)              # Vectors holding the null "correlations" (for pearson, dcor and mic respectively) for each of the nsim null datasets at a given noise level

val.rdcor2=val.rmine2=val.rA2= rep(NA,nsim2)              # Vectors holding the alternative "correlations" (for pearson, dcor and mic respectively) for each of the nsim2 alternative datasets at a given noise level

power.rdcor=power.rmine=power.rA= array(NA, c(8,num.noise))                # Arrays holding the estimated power for each of the "correlation" types, for each data type (linear, parabolic, etc...) with each noise level

## We loop through the noise level and functional form; each time we estimate a null distribution based on the marginals of the data, and then use that null distribution to estimate power

## We use a uniformly distributed x, because in the original paper the authors used the same

for(l in 1:num.noise){
  for(typ in 1:8){
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
        y=128*(x-1/3)^3-48*(x-1/3)^3-12*(x-1/3)+10* noise  * (l/num.noise) *rnorm(n)
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
      
      val.rdcor[ii]=rDcor(x,y)              # Calculate dcor
      val.rmine[ii]=rMIC(x,y)            # Calculate mic
      val.rA[ii]=rA(x,y)           # Calculate A
    }
    
    val.rmine <- val.mine[which(!is.na(val.mine))]                 # we remove the mic trials which glitch
    
    ## Next we calculate our 3 rejection cutoffs
    
    cut.rdcor=quantile(val.rdcor,.95)
    cut.rmine=quantile(val.rmine,.95)
    cut.rA=quantile(val.rA,.95)
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
      
      val.rdcor2[ii]=rDcor(x,y)
      val.rmine2[ii]=rMIC(x,y)
      val.rA2[ii]=rA(x,y)
    }
    
    ## Now we estimate the power as the number of alternative statistics exceeding our estimated cutoffs
    power.rdcor[typ,l] <- sum(val.rdcor2 > cut.rdcor)/nsim2
    power.rmine[typ,l] <- notNA.greater(val.rmine2, cut.rmine)
    power.rA[typ,l] <- sum(val.rA2 > cut.rA)/nsim2
  }
}


#ggplot
powerdf= data.frame(power.rdcor)
powerdf = t(powerdf)
#powerdf = cbind(powerdf,rep('Pearson',30))
colnames(powerdf) = c('Linear','Quadratic','Cubic','Sine: period 1/8',
                      'Sine: period 1/4', 'X^(1/4)','Circle','Step Function')

require(ggplot2)
ff = melt(powerdf)
ff = cbind(ff,rep('rDor',30))
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

#ff = allCors(ff,power.rdcor,'dcor')
ff = allCors(ff,power.rmine,'rMIC')
ff = allCors(ff,power.rA,'rA')

#write.csv(ff,'power1000.csv')
ggplot(ff, aes(x=Noise, y=Power,group=Statistic,colour=Statistic)) +geom_line() + facet_wrap(~ Form,scales="free_y")

