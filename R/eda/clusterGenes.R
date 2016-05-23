require(cluster)
pc = cor(as.matrix(genes[,c(-1,-ncol(genes))]))
R2 = pc^2
aa = aScore[-c(1,96),-c(1,2,96)]
diag(aa)= 1
rownames(aa)=colnames(aa)
nlr = aa - R2

plot(pam(as.dist(nlr),k=2)
 
par(mfrow=c(2,1))  

r2 = cor(as.matrix(genes[genes$diseased=="RED",c(-1,-ncol(genes))]))
adlupus = abs(lupusA[-c(1,96),-c(1,2,97)] - r2)
diag(adlupus) = 1
plot(agnes(as.dist(t(1-adlupus)),stand=F,method='complete'), which.plots=c(2),cex=0.7,main='',xlab='Lupus')

adcontrol = abs(controlA[-c(1,96),-c(1,2,97)] -cor(as.matrix(genes[genes$diseased=="GREEN",c(-1,-ncol(genes))])))
diag(adcontrol) = 1
plot(diana(as.dist(t(1-adcontrol))), which.plots=c(2),cex=0.7,main='',xlab='Control')
 ?agnes
require(vegan)
plot(cmdscale(as.dist(t(1-adcontrol))))
?spantree
plot(spantree(as.dist(t(1-adcontrol))))
find_k <- function(data,m=20)
{
  asw <- numeric(m)
  for (k in 2:m)
    asw[[k]] <- pam(data, k) $ silinfo $ avg.width
  k.best <- which.max(asw)
  cat("silhouette-optimal number of clusters:", k.best, "\n")
  #optimal silhoute widt at k=11 but fit not great
  plot(asw)
}
find_k(as.dist(t(1/adlupus)))
plot(pam(as.dist(t(1/adcontrol)),k=2), which.plots=c(2))

max(r2)
which(adlupus == max(adlupus), arr.ind = TRUE)
