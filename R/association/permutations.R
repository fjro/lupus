#Workflow for measuring all possible 3 way associations (X ~ Y, Z etc.) 

library(gtools)
library(stringi)
library(energy)
library(matie)

#first evaluate the indices of all possible 3 way associations
pr <- data.frame(permutations(94, 3, 1:94, repeats.allowed = F))
colnames(pr) <- c('X','Y','Z')
pr$Triplet <- rep(0,nrow(pr))

#' Forms symmetric triplets of indices
#' 
#' @param row The dataframe row
#' @return The triplet indices
triplet <- function(row) {
  ss <- sort(row[2:3])
  paste(row[1],ss[1], ss[2], sep=":")
}

pr$Triplet <- apply(pr, 1, triplet)

#now deduplicate them, recall that for some estimate of assoction, A, A(X ~ Y, Z) == A(X ~ Z, Y)
dd <- pr[!duplicated(pr$Triplet), ]
dim(dd)

write.csv(dd, "output/3way.csv")
dcor(controls[,20], controls[,3:4])
dcov.test(controls[,20], controls[,3:4])$estimates[2]

#' Computes all 3 way association using dcor, A and R2 and also computes the residual
#' non-linear association.
#' 
#' @param model A dataframe containing the indices of the 3 way associations.
#' @param data The data.
#' @return A dataframe with the results.
threeWay <- function(model, data)
{
  n <- nrow(model)
  results <- data.frame(rep(0, n), rep(0, n), rep(0, n), rep(0, n), rep(0, n))
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

#evaluate the 3 way associations separately for each group and also combined.
control3way <- threeWay(dd, controls[,-c(1,96)]) 
lupus3way <- threeWay(dd, lupus[,-c(1,96)])


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

write.csv(control3way, "output/control3way.csv")
write.csv(lupus3way, "output/lupus3way.csv")
write.csv(both3way, "output/both3way.csv")

lupus3way = read.csv('output/lupus3way.csv')
lupus3way = lupus3way[,-1]

control3way = read.csv('output/control3way.csv')
control3way = control3way[,-1]
