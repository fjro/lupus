#A Bayesian assesment of non-linear associations in the lupus data
library(spikeslab)

colnames(genes)
y = genes$diseased == "Control"
x = genes[,-c(1,96)]
obj <- spikeslab(y ~ ., data=genes[,-c(1,96)])
summary(obj)

#crossvalidation
cv.obj <- cv.spikeslab(x = x, y = y, K = 20)
#about 15 folds should do

cv.stb <- as.data.frame(cv.obj$stability)
gnet <- cv.stb$gnet
stability <- cv.stb$stability
plot(cv.obj, plot.type = "cv")
plot(gnet, stability)
which(stability >80)
cv.stb$bma
#multiclass
sparsePC.out <- sparsePC(x = genes[, -c(1,96)], y = genes$diseased, n.rep = 3)
rf.obj <- sparsePC.out$rf.obj
varImpPlot(rf.obj)

##############################################
#spike slab gam
library(spikeSlabGAM)
set.seed(100)
n <- 400
d <- data.frame(f1 = sample(gl(3, n/3)), f2 = sample(gl(4,
                                                        n/4)), x1 = runif(n), x2 = runif(n), x3 = runif(n))
# true model:
#   - interaction between f1 and x1
#   - smooth interaction between x1 and x2,
#   - x3 and f2 are noise variables without influence on y
nf1 <- as.numeric(d$f1)
d$f <- with(d, 5 * (nf1 + 2 * sin(4 * pi * (x1 - 0.2) *
                                    (x2 - 0.7)) - nf1 * x1))
d$y <- with(d, scale(f + rnorm(n)))
d$yp <- with(d, rpois(n, exp(f/10)))

# fit & display the model:
m1 <- spikeSlabGAM(y ~ x1 * f1 + f1 * f2 + x3 * f1 +
                     x1 * x2, data = d)
summary(m1)

# visualize estimates:
plot(m1)
plot(m1, cumulative = FALSE)
(plotTerm("fct(f1):fct(f2)", m1, aggregate = median))
print(p <- plotTerm("sm(x1):sm(x2)", m1, quantiles = c(0.25,
                                                       0.75), cumulative = FALSE))

# change MCMC settings and priors:
mcmc <- list(nChains = 3, burnin = 100, chainLength = 1000,
             thin = 3, reduceRet = TRUE)
hyper <- list(gamma = c(v0 = 5e-04), tau2 = c(10,
                                              r30), w = c(2, 1))

# complicated formula example, poisson response:
m2 <- spikeSlabGAM(yp ~ x1 * (x2 + f1) + (x2 + x3 + f2)^2 -
                     sm(x2):sm(x3), data = d,
                   family = "poisson", mcmc = mcmc,
                   hyperparameters = hyper)
summary(m2)
plot(m2)

# quick & dirty convergence diagnostics:
print(b <- ssGAM2Bugs(m1))
plot(b)


# fit & display the model:
xnam = c('x1', 'x2', 'f1', 'f2')
(fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"))))
m1 <- spikeSlabGAM(y ~ x1 * f1 + f1 * f2 + x3 * f1 +   x1 * x2, data = d)
m1 <- spikeSlabGAM(fmla, data = d)
summary(m1)


        
g2 = genes[,-c(1,96)]        
colnames(g2) = paste('x', 1:94)
colnames(g2) = gsub(" ", "", colnames(g2))
colnames(g2) = gsub("t", "", colnames(g2))
xnam = colnames(g2)
(fmla <- as.formula(paste("y ~ ", paste(xnam, collapse= "+"))))
m1 <- spikeSlabGAM(y ~ x1 * f1 + f1 * f2 + x3 * f1 +   x1 * x2, data = d)
m1 <- spikeSlabGAM(fmla, data = g2)
summary(m1)
y
m2 <- spikeSlabGAM(y ~ x5 + x33 + x55 + x84 + x5*x33 + x5 * x55 + x5*x84 * x33 * x55 + x33*x84 + x55*x84, data = g2)
summary(m2)
plot(m2)                                                       
plot(g2$x33, g2$x84)

yt = y[testg]
m3 <- spikeSlabGAM(yt ~ x82+x89+ x82*x83, data = g2[testg,])
summary(m3)
plot(m3)                                                       

predict(m3, newdata= g2[-testg,],type='terms')

testg = c(sample(1:53, 37), sample(54:420, 257))
