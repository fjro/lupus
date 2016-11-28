require(energy)
x <- matrix(rnorm(100), 10, 10)
y <- matrix(runif(100), 10, 10)
dx <- dist(x)
dy <- dist(y)

x <- genes[genes$diseased=='Control',-c(1,96,97)]
y <- genes[genes$diseased=='Lupus',-c(1,96,97)]
tadcor(x)
nl = nltadcor(y)
diag(nl) = 0
any(nl > 0.45)
which(nl == max(nl), arr.ind = T)

#nonlinear dcor
nldcor = function(x, y) { dcor(lm(y~x)$residuals, x) }
nlMIC = function(x, y) { mine(lm(y~x)$residuals, x)$MIC }
nlA = function(x, y) { ma.nl(x,y)$rA }

tadcor = function(dataset) 
{
  names <- names(dataset)
  lr <- length(names)
  print(paste("Processing", toString(lr), "variables.", toString(lr * (lr - 1)/2), "pairs"))
  m = matrix(nrow=lr,ncol=lr)
  diag(m) = 1
  print(paste("m = ", toString(dim(m))))
  for (i in 1:(lr - 1)) 
  {
    for (j in (i + 1):lr) 
    {
      dc = dcor(dataset[,i],dataset[,j])
      m[i,j] = dc
      m[j,i] = dc
    }
    print(paste(toString(i), " variable(s) complete", sep = ""))
  }
  colnames(m) <- names
  rownames(m) <= names
  return(m)
}

nltadcor = function(dataset) 
{
  names <- names(dataset)
  lr <- length(names)
  print(paste("Processing", toString(lr), "variables.", toString(lr * (lr - 1)/2), "pairs"))
  m = matrix(nrow=lr,ncol=lr)
  diag(m) = 1
  print(paste("m = ", toString(dim(m))))
  for (i in 1:(lr - 1)) 
  {
    for (j in (i + 1):lr) 
    {
      dc = nldcor(dataset[,i],dataset[,j])
      m[i,j] = dc
      m[j,i] = dc
    }
    print(paste(toString(i), " variable(s) complete", sep = ""))
  }
  colnames(m) <- names
  rownames(m) <= names
  return(m)
}

nltaMIC = function(dataset) 
{
  names <- names(dataset)
  lr <- length(names)
  print(paste("Processing", toString(lr), "variables.", toString(lr * (lr - 1)/2), "pairs"))
  m = matrix(nrow=lr,ncol=lr)
  diag(m) = 1
  print(paste("m = ", toString(dim(m))))
  for (i in 1:(lr - 1)) 
  {
    for (j in (i + 1):lr) 
    {
      dc = nlMIC(dataset[,i],dataset[,j])
      m[i,j] = dc
      m[j,i] = dc
    }
    print(paste(toString(i), " variable(s) complete", sep = ""))
  }
  colnames(m) <- names
  rownames(m) <= names
  return(m)
}

nltap = function(dataset) 
{
  names <- names(dataset)
  lr <- length(names)
  print(paste("Processing", toString(lr), "variables.", toString(lr * (lr - 1)/2), "pairs"))
  m = matrix(nrow=lr,ncol=lr)
  diag(m) = 1
  print(paste("m = ", toString(dim(m))))
  for (i in 1:(lr - 1)) 
  {
    for (j in (i + 1):lr) 
    {
      dc = nlA(dataset[,i],dataset[,j])
      m[i,j] = dc
      m[j,i] = dc
    }
    print(paste(toString(i), " variable(s) complete", sep = ""))
  }
  colnames(m) <- names
  rownames(m) <= names
  return(m)
}
 


require(matie)
require(minerva)
?tap
x = rnorm(100)
y = x^2
nlMIC(x,y)
plot(x,y)
plot(lm(y~x)$residuals,x)
ma(data.frame(y, x))



fit$residuals
#multivariate independence
dcor.t(x, y)
dx <- dist(t(x))
dy <- dist(t(y))
dcor.t(t(x), t(y))
?dcor.t
bcdcor(dx, dy, distance=TRUE)
dcor.ttest(t(x), t(y))


