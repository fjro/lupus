## Generate some data from the Circle example
n = 50
X = hhg.example.datagen(n, 'Circle')
plot(X[1,], X[2,])

## Compute distance matrices, on which the HHG test will be based
Dx = as.matrix(dist((X[1,]), diag = TRUE, upper = TRUE))
Dy = as.matrix(dist((X[2,]), diag = TRUE, upper = TRUE))

## Compute HHG statistics, and p-values using 1000 random permutations
set.seed(1) #set the seed for the random permutations
hhg = hhg.test(Dx, Dy, nr.perm = 1000)

## Print the  statistics and their permutation p-value

hhg

## 1.2. A null univariate example

n = 50
X = hhg.example.datagen(n, '4indclouds') 

Dx = as.matrix(dist((X[1,]), diag = TRUE, upper = TRUE))
Dy = as.matrix(dist((X[2,]), diag = TRUE, upper = TRUE))

set.seed(1) #set the seed for the random permutations
hhg = hhg.test(Dx, Dy, nr.perm = 1000)

hhg

## 1.3. A multivariate example
library(MASS)

n = 50
p = 5
x = t(mvrnorm(n, rep(0, p), diag(1, p)))
y = log(x ^ 2)
Dx = as.matrix(dist((t(x)), diag = TRUE, upper = TRUE))
Dy = as.matrix(dist((t(y)), diag = TRUE, upper = TRUE))

set.seed(1) #set the seed for the random permutations
hhg = hhg.test(Dx, Dy, nr.perm = 1000)

hhg

## 2. The k-sample test

n = 50
D = hhg.example.datagen(n, 'FourClassUniv')
Dx = as.matrix(dist(D$x, diag = TRUE, upper = TRUE))

set.seed(1) #set the seed for the random permutations
hhg = hhg.test.k.sample(Dx, D$y, nr.perm = 1000)

hhg

## End(Not run)
