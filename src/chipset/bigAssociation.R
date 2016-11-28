require(energy)
require(matie)
require(minerva)
myMA = function(x, y) { ma(data.frame(x, y))$A}

bdAssociation = data.frame(rep())
addAssociation = function(x, y, cx, cy, lx, ly, bigResults, alpha, beta, gamma)
{
  cvalues = c(cor(cx, cy), myMA(cx, cy), dcor(cx, cy), mine(cx, cy)$MIC)
  lvalues = c(cor(lx, ly), myMA(lx, ly), dcor(lx, ly), mine(lx, ly)$MIC)

  # 
#   cR2 = cor(cx, cy)
#   cA = myMA(cx, cy)
#   cDcor = dcor(cx, cy)
#   cMIC = mine(cx, cy)$MIC
#   
#   lR2 = cor(lx, ly)
#   lA = myMA(lx, ly)
#   lDcor = dcor(lx, ly)
#   lMIC = mine(lx, ly)$MIC
  
  if((any(cvalues >= alpha) | any(lvalues >= alpha)) & 
      any(lvalues[2:4] - lvalues[1] >= beta) &
      any(lvalues - cvalues >= gamma))
    {
    print('Adding...')
    bigResults = rbind(bigResults, c(x, y, cvalues,lvalues))
  }
}

addAssociation(4, 5, controlBD[,4], controlBD[,5], lupusBD[,8], lupusBD[,5], bigResults, 0.8, 0.1,0.4)





maBD = function(sampleC, sampleL, alpha = 0.8, beta=0.1, gamma=0.4)
{
  names <- names(sampleC)
  lr <- ncol(sampleC)
  total = lr * (lr - 1)/2
  print(paste("Processing", toString(lr), "variables.", toString(total), "pairs"))

  bigResults = data.frame()
  for (i in 1:(lr - 1)) 
  {
    for (j in (i + 1):lr) 
    {
      #addAssociation(i, j, sampleC[,i],sampleC[,j], sampleL[,i],sampleL[,j], bigResults, alpha, beta, gamma)
      cvalues = c(cor(sampleC[,i], sampleC[,j]), myMA(sampleC[,i], sampleC[,j]), dcor(sampleC[,i], sampleC[,j]), mine(sampleC[,i], sampleC[,j])$MIC)
      lvalues = c(cor(sampleL[,i], sampleL[,j]), myMA(sampleL[,i], sampleL[,j]), dcor(sampleL[,i], sampleL[,j]), mine(sampleL[,i], sampleL[,j])$MIC)
      
      if((any(cvalues >= alpha) | any(lvalues >= alpha)) & 
         any(lvalues[2:4] - lvalues[1] >= beta) &
         any(lvalues - cvalues >= gamma))
      {
        print('Adding...')
        bigResults = rbind(bigResults, c(i, j, cvalues,lvalues))
      }
      if(j %% 100 == 0)
      {
        print(j)
      }
    }
    
    #if(i %% 50 == 0)
    #{
      print(paste(toString(i), " variables of ", toString(total), sep = " "))
    #}
  }
  colnames(bigResults) = c('X', 'Y', 'cR2', 'cA', 'cDcor', 'cMIC', 'lR2', 'lA', 'lDcor', 'lMIC')
  
  bigResults
}


maBD2 = function(from, to, sampleC, sampleL, alpha = 0.8, beta=0.1, gamma=0.4)
{
  names <- names(sampleC)
  lr <- ncol(sampleC)
  total = lr * (lr - 1)/2
  #print(paste("Processing", toString(lr), "variables.", toString(total), "pairs"))
  
  bigResults = data.frame()
  for (i in from:to) 
  {
    for (j in (i + 1):lr) 
    {
      #addAssociation(i, j, sampleC[,i],sampleC[,j], sampleL[,i],sampleL[,j], bigResults, alpha, beta, gamma)
      cvalues = c(cor(sampleC[,i], sampleC[,j]), dcor(sampleC[,i], sampleC[,j]))
      lvalues = c(cor(sampleL[,i], sampleL[,j]), dcor(sampleL[,i], sampleL[,j]))
      
      if((any(cvalues >= alpha) | any(lvalues >= alpha)) & 
         any(lvalues[2] - lvalues[1] >= beta) &
         any(lvalues - cvalues >= gamma))
      {
        #print('Adding...')
        bigResults = rbind(bigResults, c(i, j, cvalues,lvalues))
      }
     
    }
    
    print(paste(toString(i), " variables of ", toString(total), sep = " "))
    
  }
  colnames(bigResults) = c('X', 'Y', 'cR2', 'cDcor', 'lR2', 'lDcor')
  
  bigResults
}

maBD3 = function(sampleC, sampleL, alpha = 0.8, beta=0.1, gamma=0.4)
{
  names <- names(sampleC)
  lr <- ncol(sampleC)
  total = lr * (lr - 1)/2
  print(paste("Processing", toString(lr), "variables.", toString(total), "pairs"))
  
  bigResults = data.frame()
  for (i in 1:(lr - 1)) 
  {
    for (j in (i + 1):lr) 
    {
      #addAssociation(i, j, sampleC[,i],sampleC[,j], sampleL[,i],sampleL[,j], bigResults, alpha, beta, gamma)
      cvalues = c(cor(sampleC[,i], sampleC[,j])^2, dcor(sampleC[,i], sampleC[,j]))
      lvalues = c(cor(sampleL[,i], sampleL[,j])^2, dcor(sampleL[,i], sampleL[,j]))
      
      if(any(cvalues >= alpha) | any(lvalues >= alpha)) 
      {
        print(paste("Adding ", toString(i) ,":", toString(j), sep = ""))
        bigResults = rbind(bigResults, c(i, j, cvalues,lvalues))
      }
      if(j %% 500 == 0)
      {
        print(paste("j = ", toString(j), sep = " "))
      }
    }
    
    #if(i %% 50 == 0)
    #{
    print(paste(toString(i), " variables of ", toString(lr-1), sep = " "))
    #}
  }
  colnames(bigResults) = c('X', 'Y', 'cR2', 'cDcor', 'lR2', 'lDcor')
  
  bigResults
}

maBD4 = function(sampleL, alpha = 0.8)
{
  names <- names(sampleL)
  lr <- ncol(sampleL)
  total = lr * (lr - 1)/2
  print(paste("Processing", toString(lr), "variables.", toString(total), "pairs"))
  
  bigResults = data.frame()
  for (i in 1:(lr - 1)) 
  {
    for (j in (i + 1):lr) 
    {
      value =  .dcov(sampleL[,i], sampleL[,j])
      
      if(value >= alpha)
      {
        print(paste("Adding ", toString(i) ,":", toString(j), sep = ""))
        bigResults = rbind(bigResults, c(i, j, value))
      }
      if(j %% 500 == 0)
      {
        print(paste("j = ", toString(j), sep = " "))
      }
    }
    
    #if(i %% 50 == 0)
    #{
    print(paste(toString(i), " variables of ", toString(lr-1), sep = " "))
    #}
  }
  colnames(bigResults) = c('X', 'Y', 'lDcor')
  
  bigResults
}

x <- iris[1:50, 1:4]
y <- iris[51:100, 1:4]
set.seed(111)
## R version
system.time(replicate(1000, DCOR(x, y)))
set.seed(111)
## C version
system.time(replicate(1000, .dcov(x, y)))

dim(controlBD)
dim(lupusBD)

dcor(controlBD[,3:4])
myMA(controlBD[,3], controlBD[,4])
dcor(controlBD[,3], controlBD[,4])
mine(controlBD[,3], controlBD[,4])$MIC

cor(lupusBD[,3], lupusBD[,4])
myMA(lupusBD[,3], lupusBD[,4])
dcor(lupusBD[,3], lupusBD[,4])
mine(lupusBD[,3], lupusBD[,4])$MIC

bigResults = maBD3(controlBD, lupusBD)
bigResults = maBD3(controls[, -c(1,96)], lupus[,-c(1,96)])
bigResults = maBD4(lupusBD)

bigResults = mtl2(ichunks[1], controls[1:20, 3:30], lupus[1:20,3:30])

dim(allBD)

























library(snow)

mtl2 = function(ichunk, ml, mc)
{
  maBD2(unlist(ichunk)[1], unlist(ichunk)[2], ml, mc)
}

mutlinks2 = function(cls, ml, mc)
{
  options(warn=-1)
  ichunks = list('1' = c(1, 20), '2'=c(21, 40))
  options(warn=0)
  bdR = clusterApply(cls, ichunks, mtl2, ml, mc)
  bdR
}


?clusterApply
unlist(ichunks[2])[2] 
mtl = function(ichunk, m)
{
  n = ncol(m)
  matches = 0
  for(i in ichunk)
  {
    if(i < n)
    {
      rowi = m[i,]
      matches = matches + sum(m[(i+1):n,] %*% rowi)
    }

  }
  matches
}

mutlinks = function(cls, m)
{
  m = testm
  n = nrow(m)
  nc = length(cls)
  nc = 2
  options(warn=-1)
  ichunks =  split(1:n, 1:nc)
  options(warn=0)
  counts = clusterApply(cls, ichunks, mtl, m)
  do.call(sum, counts)/ (n*(n-1)/2)
}

?makeCluster
cls = makeCluster(type="SOCK", rep('localhost',4))
mm = 4
testm = matrix(sample(0:1, mm^2, replace = T), nrow=mm)
mutlinks(cl,testm)
