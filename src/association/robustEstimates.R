#Bootstrapped estimates of mean and sd of 2 and 3 way associations. 

library(plyr)
library(energy)
library(minerva)
library(matie)

tmpfun <- function() 
{
  fff = function(df) 
  {
    data.frame( Petal.Width=df$Petal.Width, Sepal.Width=sample(df$Sepal.Width) ) 
  } 
  tmp <- ddply( iris, .(Species),fff)
  dcor(tmp$Petal.Width, tmp$Sepal.Width)
}
?sample
out <- c( dcor(iris$Petal.Width, iris$Sepal.Width), replicate(999, tmpfun()) )

bdcor = function(x, y, n = 1000)
{
  n = 1000
  jj = cbind(x,y)
  
  res = rep(0,n)
  
  for (i in 1:n)
  {
    bs = jj[sample(nrow(jj), replace = T), ]
    res[i]= dcor(bs[,1], bs[,2])
    
  }
  
  br= list('Mean' = mean(res), 'SD' = sd(res))
  return(br)
}

bcor = function(x, y, n = 1000)
{
  n = 1000
  jj = cbind(x,y)
  
  res = rep(0,n)
  
  for (i in 1:n)
  {
    bs = jj[sample(nrow(jj), replace = T), ]
    res[i]= cor(bs[,1], bs[,2])^2
    
  }
  
  br= list('Mean' = mean(res), 'SD' = sd(res))
  return(br)
}

bmic = function(x, y, n = 1000)
{
  n = 1000
  jj = cbind(x,y)
  
  res = rep(0,n)
  
  for (i in 1:n)
  {
    bs = jj[sample(nrow(jj), replace = T), ]
    res[i]= mine(bs[,1], bs[,2])$MIC
  }
  
  br= list('Mean' = mean(res), 'SD' = sd(res))
  return(br)
}

bA = function(x, y, n = 1000)
{
  n = 1000
  jj = cbind(x,y)
  
  res = rep(0,n)
  
  for (i in 1:n)
  {
    bs = jj[sample(nrow(jj), replace = T), ]
    res[i]= ma(data.frame(bs[,1], bs[,2]))$A
  }
  
  br= list('Mean' = mean(res), 'SD' = sd(res), 'res' = res)
  return(br)
}

bdcor3 = function(x, y, z, n = 1000)
{
  n = 1000
  jj = cbind(x,y, z)
  
  res = rep(0,n)
  
  for (i in 1:n)
  {
    bs = jj[sample(nrow(jj), replace = T), ]
    res[i]= dcor(bs[,1], bs[,2:3])
    
  }
  
  br= list('Mean' = mean(res), 'SD' = sd(res))
  return(br)
}

bR23 = function(x, y, z, n = 1000)
{
  n = 1000
  jj = cbind(x,y, z)
  
  res = rep(0,n)
  
  for (i in 1:n)
  {
    bs = jj[sample(nrow(jj), replace = T), ]
    lm.res = lm(bs[,1] ~ bs[,2] * bs[,3])
    res[i]= summary(lm.res)$r.squared
  }
  
  br= list('Mean' = mean(res), 'SD' = sd(res))
  return(br)
}

bA3 = function(x, y, z, n = 1000)
{
  n = 1000
  jj = cbind(x,y, z)
  
  res = rep(0,n)
  
  for (i in 1:n)
  {
    bs = jj[sample(nrow(jj), replace = T), ]
    res[i]= ma(data.frame(bs[,2],bs[,3],bs[,1]))$A
  }
  
  br= list('Mean' = mean(res), 'SD' = sd(res))
  return(br)
}


dcor(lupus[,4],lupus[,6:7])
ffg = bdcor3(lupus[,4],lupus[,6], lupus[,7])
ffg$Mean
ffg$SD

ffg = bcor(lupus[,4],lupus[,6])
ffg$Mean
ffg$SD

ffg = bmic(lupus[,4],lupus[,6])
ffg$Mean
ffg$SD

ffg = bA(lupus[,4],lupus[,6])
ffg$Mean
ffg$SD

ffg = bA3(lupus[,4],lupus[,6], lupus[,7])
ffg$Mean
ffg$SD

ffg = bR23(lupus[,4],lupus[,6], lupus[,7])
ffg$Mean
ffg$SD

dim(lupus)
lll = data.frame(lll)
colnames(lll)[1:2] = c('X','Y')
bresults = data.frame(rep(0,15),rep(0,15),rep(0,15),rep(0,15),rep(0,15),rep(0,15),rep(0,15),rep(0,15))
colnames(bresults) = c('dcorMean', 'dcorSD', 'R2Mean', 'R2SD', 'AMean', 'ASD','micMean','micSD')
for(i in 1:nrow(lll))
{
  print(i)
  res1 = bdcor(lupus[,lll$X[i]+1], lupus[,lll$Y[i]+1])
  bresults[i,1]= res1$Mean
  bresults[i,2]= res1$SD
   
  print('dcor')
    res1 = bcor(lupus[,lll$X[i]+1], lupus[,lll$Y[i]+1])
    bresults[i,3]= res1$Mean
    bresults[i,4]= res1$SD
   
    print('cor') 
    res1 = bA(lupus[,lll$X[i]+1], lupus[,lll$Y[i]+1])
    bresults[i,5]= res1$Mean
    bresults[i,6]= res1$SD
    print('A') 
    res1 = bmic(lupus[,lll$X[i]+1], lupus[,lll$Y[i]+1])
    bresults[i,7]= res1$Mean
    bresults[i,8]= res1$SD
    print('mic')
}


lll = cbind(lll, bresults)
lll= cbind(lll, lll$dcorMean - lll$R2Mean)
colnames(lll)[15] = 'rDiff'
require(xtable)

x.small <- xtable(lll[,-c(3:6,15)], label = 'tabsmall', caption = 'Strong Associations', digits = 3)

print(x.small,latex.environments = "",table.placement = 'h')
###############################################################################test

bresults = data.frame(rep(0,15),rep(0,15),rep(0,15),rep(0,15),rep(0,15),rep(0,15),rep(0,15),rep(0,15))
colnames(bresults) = c('dcorMean', 'dcorSD', 'R2Mean', 'R2SD', 'AMean', 'ASD','micMean','micSD')
for(i in 4:nrow(lllBig))
{
  print(i)
  res1 = bdcor(subLupusBD[,lllBig$X[i]], subLupusBD[,lllBig$Y[i]])
  bresults[i,1]= res1$Mean
  bresults[i,2]= res1$SD
  
  print('dcor')
  res1 = bcor(subLupusBD[,lllBig$X[i]], subLupusBD[,lllBig$Y[i]])
  bresults[i,3]= res1$Mean
  bresults[i,4]= res1$SD
  
  print('cor') 
  res1 = bA(subLupusBD[,lllBig$X[i]], subLupusBD[,lllBig$Y[i]])
  bresults[i,5]= res1$Mean
  bresults[i,6]= res1$SD
  print('A') 
  res1 = bmic(subLupusBD[,lllBig$X[i]], subLupusBD[,lllBig$Y[i]])
  bresults[i,7]= res1$Mean
  bresults[i,8]= res1$SD
  print('mic')
}


lllBig = cbind(lllBig, bresults)

require(xtable)

lastresults = arrange(lllBig[,-c(1:7)], desc(dcorMean))
x.small <- xtable(lastresults, label = 'tabsmall', caption = 'Strong Associations', digits = 3)

print(x.small,latex.environments = "",table.placement = 'h',include.rownames=FALSE)

#######################################################################
#3 way

bresults3 = data.frame(rep(0,15),rep(0,15),rep(0,15),rep(0,15),rep(0,15),rep(0,15))
colnames(bresults3) = c('dcorMean', 'dcorSD', 'R2Mean', 'R2SD', 'AMean', 'ASD')
for(i in 1:nrow(lll5))
{
  print(i)
  res1 = bdcor3(lupus[,lll5$X[i]+1], lupus[,lll5$Y[i]+1], lupus[,lll5$Z[i]+1])
  bresults3[i,1]= res1$Mean
  bresults3[i,2]= res1$SD
  
  print('dcor')
  res1 = bR23(lupus[,lll5$X[i]+1], lupus[,lll5$Y[i]+1], lupus[,lll5$Z[i]+1])
  bresults3[i,3]= res1$Mean
  bresults3[i,4]= res1$SD
  
  print('cor') 
  res1 = bA3(lupus[,lll5$X[i]+1], lupus[,lll5$Y[i]+1], lupus[,lll5$Z[i]+1])
  bresults3[i,5]= res1$Mean
  bresults3[i,6]= res1$SD
  print('A') 

}
lll5 = cbind(lll5, bresults3[1:12,])
x.small <- xtable(lll5[,c(6,7,8,16:21)], label = 'tabsmall', caption = 'Strong Associations', digits = 3)
print(x.small,latex.environments = "",table.placement = 'h')
ffg = bA3(lupus[,4],lupus[,6], lupus[,7])

lll5 = cbind(lll5, bresults3)