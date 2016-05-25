library(fastICA)

a <- fastICA(genes[,-c(1,96,97)], 6, alg.typ = "deflation", fun = "logcosh", alpha = 1,
             method = "R", row.norm = FALSE, maxit = 200,
             tol = 0.0001, verbose = TRUE)
par(mfrow = c(1, 3))
plot(a)
plot(a$X, main = "Pre-processed data")
plot(a$X %*% a$K, main = "PCA components")
plot(a$S, main = "ICA components")

#build up the frame for plotting
res <- as.data.frame(a$S)
colnames(res) = c('V1','V2','V3','V4','V5','V6')
ggplot(res, aes(x=V1, y=V2,color=genes$diseased, shape=genes$diseased,size=5)) +geom_point()
ggplot(res, aes(x=V1, y=V3,color=genes$diseased, shape=genes$diseased,size=5)) +geom_point()

#no effective class separation