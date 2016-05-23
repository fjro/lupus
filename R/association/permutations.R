library(gtools)
library(stringi)
library(energy)
library(matie)
pr = permutations(4, 3, letters[1:4], repeats.allowed = F)

pr = data.frame(permutations(94, 3, 1:94, repeats.allowed = F))
pr = cbind(pr, rep(0,nrow(pr)))
#pr = cbind(pr, paste(pr[,2], pr[,3], sep=":"))
#pr[,4] = order(pr[,4])

colnames(pr) = c('X','Y','Z','Pair')
  

for(i in 1:nrow(pr))
{
  ss = sort(c(pr[i,2],pr[i,3]))
  pr[i,4] = paste(pr[i,1],ss[1], ss[2], sep=":")
}

dim(pr)
dd = pr[!duplicated(pr$Pair), ]
dim(dd)
head(dd)
write.csv(dd, "3way.csv")
dcor(controls[,20], controls[,3:4])
dcov.test(controls[,20], controls[,3:4])$estimates[2]

threeWay = function(model, data)
{
  n = nrow(model)
  results = data.frame(rep(0, n), rep(0, n), rep(0, n), rep(0, n), rep(0, n))
  colnames(results) = c('dcor', 'A', 'R2', 'nldcor', 'nlA')
  for(i in 1:n)
  {
    X = model[i,1]
    Y = model[i,2]
    Z = model[i,3]
    #linear and nonlinear
    results[i, 1] = dcov.test(data[,X],data[,c(Y,Z)])$estimates[2]
    results[i, 2] = ma(data.frame(data[,c(Y,Z)],data[,X]))$A
    res = lm(data[,X] ~ data[,Y] * data[,Z])
    results[i, 3] = summary(res)$r.squared
    
    #residual association
    results[i, 4] = dcov.test(res$residuals,data[,c(Y,Z)])$estimates[2]
    results[i, 5] = ma(data.frame(data[,c(Y,Z)],res$residuals))$A
    
    if(i %% 100 == 0)
    {
      print(paste(toString(i), " tests complete", sep = ""))
    }
  }
  results
}

control3way = threeWay(dd, controls[,-c(1,96)]) 
lupus3way = threeWay(dd, lupus[,-c(1,96)])


summary(lm(genes[,2] ~ genes[,3] * genes[,4]))

dim(dd)
colnames(dd)
colnames(control3way)
control3way = cbind(control3way, dd)
control3way = cbind(control3way, rep('Control', nrow(dd)))
colnames(control3way)[10] = 'Disease'

lupus3way = cbind(lupus3way, dd)
lupus3way = cbind(lupus3way, rep('Lupus', nrow(dd)))
colnames(lupus3way)[10] = 'Disease'

both3way = rbind(control3way, lupus3way)

write.csv(control3way, "control3way.csv")
write.csv(lupus3way, "lupus3way.csv")
write.csv(both3way, "both3way.csv")

lupus3way = read.csv('lupus3way.csv')
lupus3way = lupus3way[,-1]

control3way = read.csv('control3way.csv')
control3way = control3way[,-1]
