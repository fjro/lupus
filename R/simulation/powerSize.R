require(ggplot2)
require(reshape)
library("energy")
library('minerva')
require('matie')

## We write the following function because MIC glitches a small percentage of the time, and we do not wish to average over those trials

notNA.greater <- function(a,b){
  ind <- which(!is.na(a))
  pow <- sum(a[ind] > b)/length(ind)
  return(pow)
}

functions = list("linear"=linear,"Quadratic"=parabolic,"Cubic"=cubic, "Fourth Root" = qroot, "Exponential" = exponential2, 
                 "Natural Log" = logE, "Sigmoid" = sigmoid, "Step"=step, "Spike" = spike,
                 "Sine: Low"= sin1, "Sine: High" = sin2, "Linear+Periodic" = linearPeriodic, "Varying Frequency" = varyingFreq,
                 "Circle" = circle, "X" = xShaped)
#functions = list("linear"=linear,"Quadratic"=parabolic)

set.seed(1)

# Here we define parameters which we use to simulate the data

#500
# The number of null datasets we use to estimate our rejection reject regions for an alternative with level 0.05
nsim=500
# Number of alternative datasets we use to estimate our power

#30
#num.noise <- 30                    # The number of different noise levels used
noise <- 3                          # A constant to determine the amount of noise

#320
# Number of data points per simulation
sizes=c(10,20,40,80,160,320,640,1280)

# Vectors holding the null "correlations" (for pearson, dcor and mic respectively) for each of the nsim null datasets at a given noise level
val.cor=val.dcor=val.mine=val.A=val.spear=rep(NA,nsim)        
# Vectors holding the alternative "correlations" (for pearson, dcor and mic respectively) for each of the nsim2 alternative datasets at a given noise level
val.cor2=val.dcor2=val.mine2=val.A2=val.spear2= rep(NA,nsim)        

power.cor=power.dcor=power.mine=power.A=power.spear= array(NA, c(length(functions),num.noise))                # Arrays holding the estimated power for each of the "correlation" types, for each data type (linear, parabolic, etc...) with each noise level

## We loop through the noise level and functional form; each time we estimate a null distribution based on the marginals of the data, and then use that null distribution to estimate power

## We use a uniformly distributed x, because in the original paper the authors used the same

for(k in 1:length(sizes))
{
  n = sizes[k]
  for(typ in 1:length(functions))
  {
    ## This next loop simulates data under the null with the correct marginals (x is uniform, and y is a function of a uniform with gaussian noise)
    
    for(ii in 1:nsim)
      {
      cat(sprintf("Noise Null %s; Type %s; sim = %s\"\n",  l, typ,ii))
      
      x=runif(n)
      y=functions[[typ]](x, 3, 10, 30)
      
      # We resimulate x so that we have the null scenario
      x <- runif(n)                       
      
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

    
    ## Next we simulate the data again, this time under the alternative
    
    for(ii in 1:nsim)
    {
      #cat(sprintf("Noise Alt %s; Type %s; sim = %s\"\n",  l, typ,ii))
      x=runif(n)
      
      y=functions[[typ]](x, 3, 10, 30)
      
      ## We again calculate our "correlations"
      val.cor2[ii]=(cor(x,y))^2
      val.dcor2[ii]=dcor(x,y)
      val.mine2[ii]=mine(x,y)$MIC
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

#Plot the results
powersizedf= data.frame(power.cor)
powersizedf = t(powersizedf)
colnames(powersizedf) = names(functions)

fff = melt(powersizedf)
fff = cbind(fff,rep('Pearson',30))
fff = cbind(fff,(1:30)/10)
colnames(fff) = c('var','Form','Power','Statistic','Noise')

allCors = function(df, x, stat)
{
  tf= data.frame(x)
  tf = t(tf)
  colnames(tf) = names(functions)
  tf = melt(tf)
  tf = cbind(tf,rep(stat,30))
  tf = cbind(tf,(1:30)/10)
  colnames(tf) = c('var','Form','Power','Statistic','Size')
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
puwerU = ggplot(off, aes(x=Noise, y=Power,group=Statistic,colour=Statistic)) +geom_line(size=1.1) + facet_wrap(~ Form,scales="free_y") + theme(legend.position="bottom")  



write.csv(ff,'power1000ULast.csv')


    