#install.packages('ff')
require(ff)

gc()


fff = cor(controls[,3:10])
#the max matrix size
sqrt(.Machine$integer.max)

#filter by sd
dim(controlBD)
bdAllSD =  apply(allBD, 2, sd)
bdControlSD =  apply(controlBD, 2, sd)
bdLupusSD =  apply(lupusBD, 2, sd)

bdAllIQR=  apply(allBD, 2, IQR)
bdControlIQR =  apply(controlBD, 2, IQR)
bdLupusIQR =  apply(lupusBD, 2, IQR)

par(mfrow=c(2,3))
hist(bdAllSD)
hist(bdControlSD)
hist(bdLupusSD)

hist(bdAllIQR)
hist(bdControlIQR)
hist(bdLupusIQR)

probs = (1:54675)/54675



quantile(bdAllIQR, probs=(53675/54675))
#probs[9657]
?quantile
#subBD = allBD[,bdAllIQR >= quantile(bdAllIQR, probs[9675])]
length(which(bdAllSD >= quantile(bdAllSD, probs)))
subBD = allBD[,bdAllIQR >= quantile(bdAllIQR, probs=(53675/54675))]
subBD2 = allBD[,bdAllSD >= quantile(bdAllSD, probs=(53675/54675))]
dim(subBD)
dim(subBD2)


subBD3 = allBD[,bdAllSD >= quantile(bdAllSD, probs=(9675/54675))]
subBD3 = allBD[,bdAllSD >= quantile(bdAllSD, probs=(44675/54675))]

subLupusBD = lupusBD[,bdLupusSD >= quantile(bdLupusSD, probs=4675/54675)]
dim(subLupusBD)




subControlBD = controlBD[,bdControlSD >= quantile(bdControlSD, probs=4675/54675)]
dim(subControlBD)

######################################################################
subBDAll = safeBDAll[,resultsBig$ProbeID[1:50000]]
dim(subBDAll)

dim(safeControlBD)
subControlBD = safeControlBD[,resultsBig$ProbeID[1:50000]]
dim(subControlBD)
subLupusBD = safeLupusBD[,resultsBig$ProbeID[1:50000]]
dim(subLupusBD)
all(colnames(subControlBD) == colnames(subLupusBD))
colnames(subControlBD)[1:20]
colnames(subLupusBD)[1:20]
dim(safeBDAll)
safeBD = safeBDAll[,resultsBig$ProbeID[1:1000]]
dim(safeBD)
bb = bigcor(subBD3)
length(bb)
any(COR > 0.940)
COR <- as.ffdf(bb)
length(bb)
bb[100]
write.table(COR, file = "bigcor.csv")
help("memory.size")
memory.size()
?ff
dimnames(bb)

x <- ff(1:12, dim=c(3,4), dimnames=list(letters[1:3], NULL))
x[ 1, 1]
typeof(COR)


myCor = cor(lupusBD[,1:20000])
myCor[lower.tri(myCor, diag = T)] <- 0
length(which(myCor > 0.7, arr.ind = T))

x = subBD3[,1:5000]


findAssociations <- function(x, alpha =0.7, nblocks = 10, ...)
{
  n <- ncol(x)

  #for computational conveniance we need to break the matrix up into eqaul sized blocks  
  if (n %% nblocks != 0) stop("Choose s different 'nblocks' so that ncol(x) %% nblocks = 0!")
  
  # split column numbers into 'nblocks' groups
  splits <- split(1:n, rep(1:nblocks, each = n/nblocks))
  
  # create all unique combinations of blocks
  blocks <- expand.grid(1:length(splits), 1:length(splits))
  blocks <- t(apply(blocks, 1, sort))
  blocks <- unique(blocks)
  
  # iterate through each block combination, calculate association matrix
  # and identify the global matrix indices of any associations >= alpha
 
  indices = data.frame()

  for (i in 1:nrow(blocks)) 
  {
    block <- blocks[i, ]
    g1 <- blocks[[block[1]]]
    g2 <- blocks[[block[2]]]
    cat("Searching block", block[1], "with Block", block[2], "\n")
    flush.console()
    
    
    amat <- cor(x[, g1], x[, g2],...)
    amat[lower.tri(amat, diag = T)] <- 0
    amat = amat^2
    
    #identify the associations of interest
    tempRes = which(amat >= alpha, arr.ind = T)
    tempRes[,1] = tempRes[,1] + g1[1] -1
    tempRes[,2] = tempRes[,2] + g2[1] -1
    indices = rbind(indices, tempRes)
  }

  gc()
  return(indices)
}

dim(subLupusBD)
myResL1 = searchCor(subLupusBD, 0.7, nblocks=10)
myResC2 = searchCor(subControlBD, 0.7, nblocks=10)

myRes = rbind(myResL1, myResC2)
myRes = cbind(myRes, paste(myRes[,1], myRes[,2], sep=":"))
colnames(myRes) = c('X','Y','Pair')


dim(myRes)
myRes = myRes[!duplicated(myRes$Pair), ]

require(energy)
require(minerva)
require(matie)


twoWayBig = function(model, data)
{
  n = nrow(model)
#   resultsBig = data.frame(rep(0, n), rep(0, n))#, rep(0, n), rep(0, n))
#   colnames(resultsBig) = c('dcor', 'A', 'R2', 'MIC')
  resultsBig = data.frame(rep(0, n), rep(0, n))
  colnames(resultsBig) = c('dcor', 'R2')
  for(i in 1:n)
  {
    X = model[i,1]
    Y = model[i,2]
   
    #linear and nonlinear
    resultsBig[i, 1] = dcor(data[,X],data[,Y])
    #resultsBig[i, 2] = ma(data.frame(data[,X],data[,Y]))$A
    resultsBig[i, 2] = cor(data[,X],data[,Y])^2
    #resultsBig[i, 4] = mine(data[,X],data[,Y])$MIC

    if(i %% 100 == 0)
    {
      print(paste(toString(i), " tests complete of ", n,  sep = ""))
    }
  }
  resultsBig
}

myResL1 = searchCor(subLupusBD, 0.7, nblocks=10)
myResC2 = searchCor(subControlBD, 0.7, nblocks=10)
dim(myResL1)
dim(myResC2)
nrow(myResC2) + nrow(myResL1)

myRes = rbind(myResL1, myResC2)
myRes = cbind(myRes, paste(myRes[,1], myRes[,2], sep=":"))
colnames(myRes) = c('X','Y','Pair')
myRes = myRes[!duplicated(myRes$Pair), ]

controlBD2way = twoWayBig(myRes, subControlBD) 
lupusBD2way = twoWayBig(myRes, subLupusBD) 


length(which((lupusBD2way$dcor > 0.90 & lupusBD2way$R2 < 0.75)))

write.csv(lupusBD2way, 'lupusBD2way.csv')
write.csv(controlBD2way, 'controlBD2way.csv')


bd2wayResults = read.csv('bd2wayResults.csv', header = T)
lupusBD2way = read.csv('lupusBD2way.csv', header = T)
controlBD2way = read.csv('controlBD2way.csv', header = T) 
lupusBD2way = lupusBD2way[,-1]
controlBD2way = controlBD2way[,-1]

controlBD2way = cbind(bd2wayResults[,2:4], controlBD2way)
lupusBD2way = cbind(bd2wayResults[,2:4], lupusBD2way)

dim(lupusBD2way)

head(myRes)

myRes = cbind(myRes, controlBD2way)
colnames(myRes)[4:5] = c('Ldcor', 'LR2')

#colnames(myRes)[4:7] = c('Ldcor', 'LA', 'LR2', 'LMIC')

myRes = cbind(myRes, lupusBD2way)
colnames(myRes)[6:7] = c('Cdcor', 'CR2')
write.csv(myRes, 'bd2wayResults.csv')

plot(table(myRes$X[which(lupusBD2way$dcor > 0.90 & lupusBD2way$R2 < 0.8)]))


controlBD2way = controlBD2way[,c(1:3,15:18)]
dim(controlBD2way)
lupusBD2way = cbind(myRes[,1:3], lupusBD2way)
controlBD2way = cbind(myRes[,1:3], controlBD2way)
############################################# bug checking

for (i in 1:50)
{
  xc = subControlBD[,myResC2[i,1]]
  yc = subControlBD[,myResC2[i,2]]
  #if(cor(xc,yc)^2 < 0.7)
  {
    print(cor(xc,yc)^2)
  }
}

for (i in 1:200)
{
  xc = subLupusBD[,myResL1[i,1]]
  yc = subLupusBD[,myResL1[i,2]]
  
}


dim(controlBD2way)
dim(myRes)

controlBD2way = cbind(myRes[,1:3], controlBD2way)
lupusBD2way = cbind(myRes[,1:3], lupusBD2way[,1:4])

#swap!
temp1 = controlBD2way
temp2 = lupusBD2way

controlBD2way = lupusBD2way
lupusBD2way = temp1


require(energy)
xc = subControlBD[,49]
yc = subControlBD[,50]
dcor(xc,yc)

xl = subLupusBD[,49]
yl = subLupusBD[,50]
dcor(xl,yl)

controlBD2way[1:3,]
dim(myRes)
