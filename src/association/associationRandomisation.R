#' sign based randomisation tests for differences in 2 way association
testDiffSign <- function(x1, y1, x2, y2, n = 1000)
{

  jj = cbind(x1,y1)
  kk = cbind(x2,y2)
  res1 = rep(0,n)
  res2 = rep(0,n)
  tres = rep(0,n)
  
  tObs = dcor(x1,y1) - dcor(x2, y2)
  
  for (i in 1:n) {
    bs1 = jj[sample(nrow(jj), replace = T), ]
    res1[i]= dcor(bs1[,1], bs1[,2])
    
    bs2 = kk[sample(nrow(kk), replace = T), ]
    res2[i]= dcor(bs2[,1], bs2[,2])
  }
  
 
  br= list('Mean1' = mean(res1), 'SD1' = sd(res1), 
           'Mean2' = mean(res2), 'SD2' = sd(res2), 
           'sign' = binom.test(x=sum(res1 >= res2), n=n, p=0.5, 
                               alternative="two.sided")$p.value, 'pvalueA' = 1-mean(abs(res1 -res2) >= 0), 
           'Res1' = res1, 
           'Res2' = res2, 
           'Wilcox'=wilcox.test(res1, res2, paired =T)$p.value)
  return(br)
}

#' sign based randomisation tests for differences in 3 way association
testDiffSign3 <- function(x1, y1, z1, x2, y2, z2, n = 1000)
{
  n = 1000
  jj = cbind(x1,y1, z1)
  kk = cbind(x2,y2, z2)
  res1 = rep(0,n)
  res2 = rep(0,n)
  
  for (i in 1:n)
  {
    
    bs1 = jj[sample(nrow(jj), replace = T), ]
    res1[i]= dcor(bs1[,1], bs1[,2:3])
    
    bs2 = kk[sample(nrow(kk), replace = T), ]
    res2[i]= dcor(bs2[,1], bs2[,2:3])
  }
  
  br= list('Mean1' = mean(res1), 'SD1' = sd(res1), 'Mean2' = mean(res2), 'SD2' = sd(res2), 'sign' = binom.test(x=sum(res1 >= res2), n=n, p=0.5, alternative="two.sided")$p.value, 'pvalueA' = 1-mean(abs(res1 -res2) >= 0), 'Res1' = res1, 'Res2' = res2, 'Wilcox'=wilcox.test(res1, res2, paired =T)$p.value)
  return(br)
}


dcor(lupus[,31], lupus[,32])
dcor(controls[,31], controls[,32])
testDiff(lupus[,42], lupus[,43], controls[,42], controls[,43])

testDiff(lupus[,31], lupus[,32], controls[,31], controls[,32])
testDiffSign(lupus[,42], lupus[,43], controls[,42], controls[,43]) #97 and 98

x = 31
y = 32
x = 83
y = 90
testDiff(lupus[,x], lupus[,y], controls[,x], controls[,y])
bres = testDiffSign(lupus[,x], lupus[,y], controls[,x], controls[,y]) #97 and 98
bres$pvalueA
bres$Wilcox
bres$sign

tdRes = data.frame(rep(2, nrow(lll)), rep(2, nrow(lll)), rep(2, nrow(lll)))
colnames(tdRes) = c('b','w', 's')
for(i in 1:nrow(lll))
{ 
  print(i)
  temp =  testDiffSign(lupus[,lll$X[i]+1], lupus[,lll$Y[i]+1], controls[,lll$X[i]+1], controls[,lll$Y[i]+1])
  tdRes[i,1] = temp$pvalueA
  tdRes[i,2] = temp$Wilcox
  tdRes[i,3] = temp$sign
}

i=1
dcor(lupus[,42], lupus[,43])
dcor(controls[,42], controls[,43])

testDiff(lupus[,31], lupus[,32], controls[,31], controls[,32])
testDiff(lupus[,45], lupus[,70], controls[,45], controls[,70])
testDiff(lupus[,42], lupus[,43], controls[,42], controls[,43])

testDiff3(lupus[,lll5$X[i]+1], lupus[,lll5$Y[i]+1], lupus[,lll5$Z[i]+1], controls[,lll5$X[i]+1], controls[,lll5$Y[i]+1], controls[,lll5$Z[i]+1])
which(controlDcor>0.90 & lupusDcor <0.98, arr.ind = T)
which(controlDcor == lupusDcor & lupusDcor < 1, arr.ind = T)

tdRes5 = data.frame(rep(2, nrow(lll5)),rep(2, nrow(lll5)),rep(2, nrow(lll5)),rep(2, nrow(lll5)),rep(2, nrow(lll5)),rep(2, nrow(lll5)),rep(2, nrow(lll5)))
for(i in 1:nrow(lll5))
{ 
  print(i)
  td = testDiffSign3(lupus[,lll5$X[i]+1], lupus[,lll5$Y[i]+1], lupus[,lll5$Z[i]+1], controls[,lll5$X[i]+1], controls[,lll5$Y[i]+1], controls[,lll5$Z[i]+1])
  tdRes5[i, 1] = td$Mean1
  tdRes5[i, 2] = td$SD1
  tdRes5[i, 3] = td$Mean2
  tdRes5[i, 4] = td$SD2
  tdRes5[i, 5] = td$pvalueA
  tdRes5[i,6] = td$Wilcox
  tdRes5[i,7] = td$sign
}
colnames(tdRes5) = c('meanDcorL', 'sdDcorL', 'meanDcorC', 'sdDcorC', 'pvalue','Wilcox','sign')

tdRes6 = data.frame(rep(2, nrow(lll6)),rep(2, nrow(lll6)),rep(2, nrow(lll6)),rep(2, nrow(lll6)),rep(2, nrow(lll6)))
for(i in 1:nrow(lll6))
{ 
  print(i)
  td = testDiff3(lupus[,lll6$X[i]+1], lupus[,lll6$Y[i]+1], lupus[,lll6$Z[i]+1], controls[,lll6$X[i]+1], controls[,lll6$Y[i]+1], controls[,lll6$Z[i]+1])
  tdRes6[i, 1] = td$Mean1
  tdRes6[i, 2] = td$SD1
  tdRes6[i, 3] = td$Mean2
  tdRes6[i, 4] = td$SD2
  tdRes6[i, 5] = td$pvalueA
}

tdRes3 = data.frame(rep(2, nrow(lll3)),rep(2, nrow(lll3)),rep(2, nrow(lll3)),rep(2, nrow(lll3)),rep(2, nrow(lll3)))
for(i in 1:nrow(lll3))
{ 
  print(i)
  td = testDiff3(lupus[,lll3$X[i]+1], lupus[,lll3$Y[i]+1], lupus[,lll3$Z[i]+1], controls[,lll3$X[i]+1], controls[,lll3$Y[i]+1], controls[,lll3$Z[i]+1])
  tdRes3[i, 1] = td$Mean1
  tdRes3[i, 2] = td$SD1
  tdRes3[i, 3] = td$Mean2
  tdRes3[i, 4] = td$SD2
  tdRes3[i, 5] = td$pvalueA
}

######################################