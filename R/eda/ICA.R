install.packages('fastICA')
require(fastICA)
#---------------------------------------------------
#Example 1: un-mixing two mixed independent uniforms
#---------------------------------------------------
S <- matrix(runif(10000), 5000, 2)
A <- matrix(c(1, 1, -1, 3), 2, 2, byrow = TRUE)
X <- S %*% A
a <- fastICA(X, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
             method = "C", row.norm = FALSE, maxit = 200,
             tol = 0.0001, verbose = TRUE)
par(mfrow = c(1, 3))
plot(a$X, main = "Pre-processed data")
plot(a$X %*% a$K, main = "PCA components")
plot(a$S, main = "ICA components")
#--------------------------------------------
#Example 2: un-mixing two independent signals
#--------------------------------------------
S <- cbind(sin((1:1000)/20), rep((((1:200)-100)/100), 5))
A <- matrix(c(0.291, 0.6557, -0.5439, 0.5572), 2, 2)
X <- S %*% A
a <- fastICA(X, 2, alg.typ = "parallel", fun = "logcosh", alpha = 1,
             method = "R", row.norm = FALSE, maxit = 200,
             tol = 0.0001, verbose = TRUE)
par(mfcol = c(2, 3))
plot(1:1000, S[,1 ], type = "l", main = "Original Signals",
     xlab = "", ylab = "")
plot(1:1000, S[,2 ], type = "l", xlab = "", ylab = "")
plot(1:1000, X[,1 ], type = "l", main = "Mixed Signals",
     xlab = "", ylab = "")
plot(1:1000, X[,2 ], type = "l", xlab = "", ylab = "")
plot(1:1000, a$S[,1 ], type = "l", main = "ICA source estimates",
     xlab = "", ylab = "")
plot(1:1000, a$S[, 2], type = "l", xlab = "", ylab = "")
#-----------------------------------------------------------
#Example 3: using FastICA to perform projection pursuit on a
# mixture of bivariate normal distributions
#-----------------------------------------------------------
if(require(MASS)){
  x <- mvrnorm(n = 1000, mu = c(0, 0), Sigma = matrix(c(10, 3, 3, 1), 2, 2))
  x1 <- mvrnorm(n = 1000, mu = c(-1, 2), Sigma = matrix(c(10, 3, 3, 1), 2, 2))
  X <- rbind(x, x1)
  a <- fastICA(X, 2, alg.typ = "deflation", fun = "logcosh", alpha = 1,
               method = "R", row.norm = FALSE, maxit = 200,
               tol = 0.0001, verbose = TRUE)
  par(mfrow = c(1, 3))
  plot(a$X, main = "Pre-processed data")
  plot(a$X %*% a$K, main = "PCA components")
  plot(a$S, main = "ICA components")
}

X <- genes[,c(-1,-ncol(genes))]
a <- fastICA(X, 6, alg.typ = "deflation", fun = "logcosh", alpha = 1,
             method = "R", row.norm = FALSE, maxit = 200,
             tol = 0.0001, verbose = TRUE)
par(mfrow = c(1, 3))
plot(a)
plot(a$X, main = "Pre-processed data")
plot(a$X %*% a$K, main = "PCA components")
plot(a$S, main = "ICA components")
res = as.data.frame(a$S)
colnames(res) = c('V1','V2','V3','V4','V5','V6')
ggplot(res, aes(x=V1, y=V2,color=genes$diseased, shape=genes$diseased,size=5)) +geom_point()
ggplot(res, aes(x=V3, y=V6,color=genes$diseased, shape=genes$diseased,size=5)) +geom_point()
a$S
