#Workflow for measuring all possible 3 way associations (X ~ Y, Z etc.) 

library(gtools)
library(stringi)
library(energy)
library(matie)
library(dplyr)
library(microbenchmark)

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

#' Computes all 3 way associations using dcor, A and R2 and also computes the residual
#' non-linear association.
#' 
#' @param model A dataframe containing the indices of the 3 way associations.
#' @param data The data.
#' @return A dataframe with the results.
threeWayAll <- function(models, data, ncores = 'all') {
  require(data.table, quietly = T)
  require(snowfall, quietly = T)
  require(parallel, quietly = T)
  cat(format(Sys.time(),usetz = TRUE), ": Calculating ", nrow(models), " associations \n")
  
  if (ncores == 'all') {
    ncores <- detectCores() -1
  }
  
  sfInit( parallel=TRUE, cpus=ncores, slaveOutfile="log.txt")
  sfExport('threeWayAssociation')
  
  result <- rbindlist(sfLapply(models,
                               threeWayAssociation,
                               data = data))
  sfStop()
  cat("Finished at ", format(Sys.time(),usetz = TRUE), '\n')
  as.data.frame(result)
}

#' Computes a single set of 3 way association estimates for a given model.
#' 
#' @param model The triplet of column indices to use.
#' @param data The dataframe.
#' @return A dataframe.
threeWayAssociation <- function(model, data){
  results <- data.frame(t(rep(0, 5)))
  colnames(results) = c('dcor', 'A', 'R2', 'nldcor', 'nlA')
  X <- unlist(model[1])
  Y <- unlist(model[2])
  Z <- unlist(model[3])
    
  #linear and nonlinear association
  results$dcor <- dcov.test(data[,X],data[,c(Y,Z)])$estimates[2]
  results$A <- ma(data.frame(data[,c(Y,Z)],data[,X]))$A
  res <- lm(data[,X] ~ data[,Y] * data[,Z])
  results$R2 <- summary(res)$r.squared
    
  #residual nonlinear association
  results$nldcor <- dcov.test(res$residuals,data[,c(Y,Z)])$estimates[2]
  results$nlA <- ma(data.frame(data[,c(Y,Z)],res$residuals))$A
    
  results
}


threeWayAssociation(pr[1,], data = data.frame(genes[, -c(1,96,97)]))
threeWayAll(pr[1:10,], data = data.frame(genes[, -c(1,96,97)]))

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
