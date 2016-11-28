library(cluster)
library(vegan)

#get the Pearson correlation and R^2
pc <- cor(as.matrix(genes[,-c(1,96,97)]))
R2 <- pc^2
aa <- aScore[-c(1,96),-c(1,2,96)]
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


plot(cmdscale(as.dist(t(1-adcontrol))))

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


require(devtools)
library(rafalib)
mypar2(1,1)
hc <- hclust(d)

plot(hc,labels=genes[,96])
as.character()
myplclust(hc, labels=genes[,96], lab.col=as.fumeric(as.character(genes[,96])))

abline(h=12)
hclusters <- cutree(hc, h=10)
table(true=genes[,96], cluster=hclusters)

#k-means
plot(e[1,], e[2,])
set.seed(1)
km <- kmeans(t(e[1:2,]), centers=7)
names(km)
plot(e[1,], e[2,], col=km$cluster, pch=16)
plot(e[1,], e[2,], col=as.fumeric(tissue), pch=16)
table(true=tissue,cluster=km$cluster)

#kmeans with cmd
mds <- cmdscale(d)
plot(mds[,1], mds[,2]) 
km <- kmeans(genes[,-c(1)], centers=7)
plot(mds[,1], mds[,2], col=km$cluster, pch=16)
table(true=tissue,cluster=km$cluster)
