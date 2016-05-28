#Workflow for measuring all possible 3 way associations (X ~ Y, Z etc.) between variables. A number of
#association estimates are used and the non-linear residual association is also computed.

library(gtools)
library(stringi)
library(energy)
library(matie)
library(dplyr)
library(microbenchmark)
library(snowfall)

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
  result <- rbindlist(sfApply(models, 1,
                               threeWayAssociation,
                               data = data))
  sfStop()
  cat("Finished at ", format(Sys.time(),usetz = TRUE), '\n')
  result <- as.data.frame(result)
  result <- cbind(models, result)
  result
}

#' Computes a single set of 3 way association estimates for a given model.
#' 
#' @param model The triplet of column indices to use.
#' @param data The dataframe.
#' @return A dataframe.
threeWayAssociation <- function(model, data) {
  require(energy)
  require(matie)
  results <- data.frame(t(rep(0, 5)))
  colnames(results) = c('dcor', 'A', 'R2', 'nldcor', 'nlA')
  X <- as.numeric(model[1])
  Y <- as.numeric(model[2])
  Z <- as.numeric(model[3])
 
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

#run a benchmark to get a feel for the execution time
microbenchmark(times = 1, thousand = threeWayAll(pr[1:1000,], data = data.frame(genes[, -c(1,96,97)])))
#1000 associations takes about 136 seconds using 7 cores on an Intel I7 using Microsoft R Open

#evaluate the 3 way associations separately for each group and also combined.
#there are about 400000 3 way associations so this needs to run overnight
control3way <- threeWayAll(dd, data.frame(controls[,-c(1,96, 97)]))
lupus3way <- threeWayAll(dd, data.frame(lupus[,-c(1,96, 97)]))

#build up the results and save them to csv for later use.
control3way$Disease <- rep('Control', nrow(dd))
lupus3way$Disease <- rep('Lupus', nrow(dd))
both3way <- rbind(control3way, lupus3way)

write.csv(control3way, "output/control3way.csv", row.names = F)
write.csv(lupus3way, "output/lupus3way.csv", row.names = F)
write.csv(both3way, "output/both3way.csv", row.names = F)
