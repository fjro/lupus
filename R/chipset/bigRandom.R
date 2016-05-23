library(plyr)

#############################################################################
# sign based randomisaations
binom.test(x=8, n=10, p=0.5, alternative="two.sided")$p.value

testDiffSign = function(x1, y1, x2, y2, n = 1000)
{
  
  jj = cbind(x1,y1)
  kk = cbind(x2,y2)
  res1 = rep(0,n)
  res2 = rep(0,n)
  tres = rep(0,n)
  
  tObs = dcor(x1,y1) - dcor(x2, y2)
  
  for (i in 1:n)
  {
    
    bs1 = jj[sample(nrow(jj), replace = T), ]
    res1[i]= dcor(bs1[,1], bs1[,2])
    
    bs2 = kk[sample(nrow(kk), replace = T), ]
    res2[i]= dcor(bs2[,1], bs2[,2])
  }
  
  
  br= list('Mean1' = mean(res1), 'SD1' = sd(res1), 'Mean2' = mean(res2), 'SD2' = sd(res2), 'sign' = binom.test(x=sum(res1 >= res2), n=n, p=0.5, alternative="two.sided")$p.value, 'pvalueA' = 1-mean(abs(res1 -res2) >= 0), 'Res1' = res1, 'Res2' = res2, 'Wilcox'=wilcox.test(res1, res2, paired =T)$p.value)
  return(br)
}

testDiffSign(safeLupusBD[, lllBig$X[1]], safeLupusBD[, lllBig$Y[1]], safeControlBD[, lllBig$X[1]], safeControlBD[, lllBig$Y[1]])


?binom.test
dcor(lupus[,31], lupus[,32])
dcor(controls[,31], controls[,32])
testDiff(lupus[,42], lupus[,43], controls[,42], controls[,43])

testDiff(lupus[,31], lupus[,32], controls[,31], controls[,32])
testDiffSign(lupus[,42], lupus[,43], controls[,42], controls[,43]) #97 and 98

tdRes = data.frame(rep(2, nrow(lllBig)), rep(2, nrow(lllBig)), rep(2, nrow(lllBig)))
colnames(tdRes) = c('b','w', 's')
for(i in 1:nrow(lllBig))
{ 
  print(i)
  temp =  testDiffSign(safeLupusBD[, lllBig$X[i]], safeLupusBD[, lllBig$Y[i]], safeControlBD[, lllBig$X[i]], safeControlBD[, lllBig$Y[i]])
  tdRes[i,1] = temp$pvalueA
  tdRes[i,2] = temp$Wilcox
  tdRes[i,3] = temp$sign
}

cdRes = data.frame(rep(2, nrow(cccBig)), rep(2, nrow(cccBig)), rep(2, nrow(cccBig)))
colnames(cdRes) = c('b','w', 's')
for(i in 1:nrow(cccBig))
{ 
  print(i)
  temp =  testDiffSign(safeLupusBD[, cccBig$X[i]], safeLupusBD[, cccBig$Y[i]], safeControlBD[, cccBig$X[i]], safeControlBD[, cccBig$Y[i]])
  cdRes[i,1] = temp$pvalueA
  cdRes[i,2] = temp$Wilcox
  cdRes[i,3] = temp$sign
}


######################################