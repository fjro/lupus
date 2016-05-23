require(MASS)
require(vegan)
require(energy)
require(grid)
require(ggplot2)
#get the data
pc = cor(as.matrix(genes[,c(-1,-ncol(genes))]))
R2 = pc^2
aa = aScore[-c(1,96),-c(1,2,96)]
diag(aa)= 1
rownames(aa)=colnames(aa)
nlr = aa - R2
diag(nlr)=1

?dcor


#nonlineat distance scaled
#not much to see
mds1 = cmdscale(as.dist(1-nlr),eig=T)
plot(mds1$points[,1],mds1$points[,2])

#plot mds of R2 and A side by side
samScree(R2,20)
samScree(aa,20)

#mds all samples
#mdsR2 = cmdscale(as.dist(1-R2),eig=T)
mdsR2 = sammon(as.dist(1-R2))
mdsdfR2 = as.data.frame(mdsR2$points)
mdsdfR2 = cbind(mdsdfR2, colnames(genes[-c(1,96,97)]))
mdsdfR2 = cbind(mdsdfR2, rep('R2',94))
colnames(mdsdfR2)[3:4] = c("Gene","Association")

#mdsA = cmdscale(as.dist(1-aa),eig=T)
mdsA = sammon(as.dist(1-aa))
mdsdfA = as.data.frame(mdsA$points)
mdsdfA = cbind(mdsdfA, colnames(genes[-c(1,96,97)]))
mdsdfA = cbind(mdsdfA, rep('A',94))
colnames(mdsdfA)[3:4] = c("Gene","Association")
mdsdf = rbind(mdsdfR2, mdsdfA)

ll= rep(1:94, 2)
ggplot(mdsdf, aes(x=V1, y=V2, label=ll)) +geom_text() + facet_grid(Association~.)  

#by association type and disease status 
#mdsR2 = cmdscale(as.dist(1-R2),eig=T)
R2C = controlCM^2
mdsCR2 = sammon(as.dist(1-R2C))
mdsdfCR2 = as.data.frame(mdsCR2$points)
mdsdfCR2 = cbind(mdsdfCR2, colnames(genes[-c(1,96,97)]))
mdsdfCR2 = cbind(mdsdfCR2, rep('R2',94))
mdsdfCR2 = cbind(mdsdfCR2, rep('Control',94))
colnames(mdsdfCR2)[3:5] = c("Gene","Association","Disease")

R2L = lupusCM^2
mdsLR2 = sammon(as.dist(1-R2L))
mdsdfLR2 = as.data.frame(mdsLR2$points)
mdsdfLR2 = cbind(mdsdfLR2, colnames(genes[-c(1,96,97)]))
mdsdfLR2 = cbind(mdsdfLR2, rep('R2',94))
mdsdfLR2 = cbind(mdsdfLR2, rep('Lupus',94))
colnames(mdsdfLR2)[3:5] = c("Gene","Association","Disease")
mdsdf = rbind(mdsdfCR2, mdsdfLR2)

#A
mdsAC = sammon(as.dist(1-controlA))
mdsdfAC = as.data.frame(mdsAC$points)
mdsdfAC = cbind(mdsdfAC, colnames(genes[-c(1,96,97)]))
mdsdfAC = cbind(mdsdfAC, rep('A',94))
mdsdfAC = cbind(mdsdfAC, rep('Control',94))
colnames(mdsdfAC)[3:5] = c("Gene","Association","Disease")
mdsdf = rbind(mdsdf, mdsdfAC)

mdsAL = sammon(as.dist(1-lupusA))
mdsdfAL = as.data.frame(mdsAL$points)
mdsdfAL = cbind(mdsdfAL, colnames(genes[-c(1,96,97)]))
mdsdfAL = cbind(mdsdfAL, rep('A',94))
mdsdfAL = cbind(mdsdfAL, rep('Lupus',94))
colnames(mdsdfAL)[3:5] = c("Gene","Association","Disease")
mdsdf = rbind(mdsdf, mdsdfAL)

#dcor
mdsDC = sammon(as.dist(1-controlDcor))
mdsdfDC = as.data.frame(mdsDC$points)
mdsdfDC = cbind(mdsdfDC, colnames(genes[-c(1,96,97)]))
mdsdfDC = cbind(mdsdfDC, rep('dcor',94))
mdsdfDC = cbind(mdsdfDC, rep('Control',94))
colnames(mdsdfDC)[3:5] = c("Gene","Association","Disease")
mdsdf = rbind(mdsdf, mdsdfDC)

mdsDL = sammon(as.dist(1-lupusDcor))
mdsdfDL = as.data.frame(mdsDL$points)
mdsdfDL = cbind(mdsdfDL, colnames(genes[-c(1,96,97)]))
mdsdfDL = cbind(mdsdfDL, rep('dcor',94))
mdsdfDL = cbind(mdsdfDL, rep('Lupus',94))
colnames(mdsdfDL)[3:5] = c("Gene","Association","Disease")
mdsdf = rbind(mdsdf, mdsdfDL)

ll= rep(1:94, 6)
ggplot(mdsdf, aes(x=V1, y=V2, label=ll)) +geom_text(size=2) + facet_grid(Association~Disease)  

## look at lower components
par(mfrow=c(3,1))
mdsDL = sammon(as.dist(1-lupusCM), k=4)
plot(mdsDL$points[,2], mdsDL$points[,3])

mdsDL = sammon(as.dist(1-lupusA), k=4)
plot(mdsDL$points[,2], mdsDL$points[,3])

mdsDL = sammon(as.dist(1-lupusDcor), k=4)
plot(mdsDL$points[,2], mdsDL$points[,3])


samScree(lupusDcor,20)
samScree(aa,20)













#clustered genes based on R2
colnames(genes[which((mdsdfR2$V1 > 0.4) & (mdsdfR2$V2 < 0.25) & (mdsdfR2$V2 > -0.25))])

#clustered genes based on A
colnames(genes[which((mdsdfA$V1 > 0.42) & (mdsdfA$V2 < 0.25) & (mdsdfA$V2 > -0.25))])
colnames(genes[c(20,21,23,25)])
colnames(genes[c(58,78,53,9,10,11,12)])
colnames(genes[c(7,54)])


sp1 = spantree(as.dist(1-R2))
plot(sp1, sammon(as.dist(1-R2),k=2),pch=16,col=2,cex=1.5)
#text(cmdscale(dm,k=k),names(data),cex=0.8,adj=c(0.5,-0.6))

sp2 = spantree(as.dist(1-aa))
plot(sp2, sammon(as.dist(1-aa),k=2),pch=16,col=2,cex=1.5)

plot(mdsA$points[,1],mdsA$points[,2])
plot(mdsA$points[,3],mdsA$points[,4])
plot(mdsA$points[,5],mdsA$points[,6])
plot(mdsA$points[,7],mdsA$points[,8])

sammon(as.dist(1-aa))

plot(mds0$points[,1],mds0$points[,2])

ggplot(as.data.frame(nfl.pca$x), aes(x=PC1, y=PC2, label=rownames(nfl2000))) +
  geom_point(alpha=.3) +geom_text()  + coord_fixed(ratio = 1,xlim = c(-6,6), ylim = c(-6, 6))

CMDscreeplot(as.dist(1-aa),maxcomp=10)

plot(mds2$points[,1],mds2$points[,2])
plot(1:94,mds2$eig)
#maybe something in this

k=8
mds2 = cmdscale(as.dist(1-aa),eig=T,k=k)
mdsdf = as.data.frame(mds2$points)
mdsdf = cbind(mdsdf, colnames(genes[-c(1,96,97)]))
colnames(mdsdf)[k+1] = "Gene"
ggplot(mdsdf, aes(x=V1, y=V2, label=Gene)) +geom_text() 
#other projects dont discrinate
#ggplot(mdsdf, aes(x=V3, y=V4)) +geom_point(size =4) + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))
#ggplot(mdsdf, aes(x=V5, y=V6)) +geom_point(size =4) + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))
#ggplot(mdsdf, aes(x=V7, y=V8)) +geom_point(size =4) + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))

#try lupus and contorls separately
mdsL = sammon(as.dist(1-adlupus),k=2)
mdsdfL = as.data.frame(mdsL$points)
mdsdfL = cbind(mdsdfL, colnames(genes[-c(1,96,97)]))
colnames(mdsdfL)[k+1] = "Gene"
ggplot(mdsdfL, aes(x=V1, y=V2, label=Gene)) +geom_text() 
ggplot(mdsdfL, aes(x=V3, y=V4, label=Gene)) +geom_text() 
ggplot(mdsdfL, aes(x=V5, y=V6, label=Gene)) +geom_text() 

mdsC = cmdscale(as.dist(1-adcontrol),eig=T,k=k)
mdsdfC = as.data.frame(mdsC$points)
mdsdfC = cbind(mdsdfC, colnames(genes[-c(1,96,97)]))
colnames(mdsdfC)[k+1] = "Gene"
ggplot(mdsdfC, aes(x=V1, y=V2, label=Gene)) +geom_text() 
ggplot(mdsdfC, aes(x=V3, y=V4, label=Gene)) +geom_text() 
ggplot(mdsdfC, aes(x=V5, y=V6, label=Gene)) +geom_text() 


#variations of standard disatnces and krusakal or sammon scaling failed to discriminate
dg = dist(genes[,-c(1,97)],method="minkowski")
mds1 = isoMDS(dg, k=6)
mdsdf = as.data.frame(mds1$points)
mdsdf = cbind(mdsdf, genes$diseased)
mdsdf = cbind(mdsdf, genes$Type)
colnames(mdsdf)[7:8] = c("Disease","Type")
ggplot(mdsdf, aes(x=V1, y=V2,color=Type, shape=Disease)) +geom_point(size =4) + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))
ggplot(mdsdf, aes(x=V3, y=V4,color=Type, shape=Disease)) +geom_point(size =4) + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))
ggplot(mdsdf, aes(x=V5, y=V6,color=Type, shape=Disease)) +geom_point(size =4) + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))

#try A based sample distance
tt = t(genes[,-c(1,97)])
aSamples = tap(as.data.frame(tt))

#does not work
mds2 = cmdscale(as.dist(aSamples), k=6)
mdsdf = as.data.frame(mds2)
mdsdf = cbind(mdsdf, genes$diseased)
mdsdf = cbind(mdsdf, genes$Type)
colnames(mdsdf)[7:8] = c("Disease","Type")
ggplot(mdsdf, aes(x=V1, y=V2,color=Type, shape=Disease)) +geom_point(size =4) + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))
ggplot(mdsdf, aes(x=V3, y=V4,color=Type, shape=Disease)) +geom_point(size =4) + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))
ggplot(mdsdf, aes(x=V5, y=V6,color=Type, shape=Disease)) +geom_point(size =4) + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))


CMDscreeplot(as.matrix(genes[,-c(1,97)]),maxcomp=20)
image(as.matrix(nlr))

samScree = function(data,k=10)
{ 
  dm = as.dist(1-data)
  stress = rep(0,k)
  #k= ncol(data)-2
  for(i in 1:k)
  {
    stress[i] = sammon(dm,k=i)$stress
  }
  
  plot(1:k,stress)
}



cmdScree = function(mydata,raw=T,abs=T,maxcomp=10) {
  n=min(dim(mydata))-1
  #x<-cmdscale(as.matrix(mydata),k=n,eig=TRUE)
  x<-cmdscale(mydata,k=n,eig=TRUE)
  #x<-cmdscale(as.matrix(mydata))
  
  
  k<-min(n,maxcomp)
  w<-c(0:k)
  
  if (raw==T) y<-x$eig
  
  
  if (raw==F & abs==T) y<-abs(x$eig)
  if(raw==F & abs==F) y<- x$eig*x$eig
  y<-c(0,y)
  z<-100*cumsum(y)/sum(y)
  
  if (raw==T) { 
    plot(w,z[1:(k+1)],type="l",xlab="number of dimensions",
         cex.main=1.5, lwd=3, col="red",
         ylim=c(0,100),
         ylab="cumulative percentage of total eigenvlaues",
         main="Scree plot of eigenvalues",
         xaxt="n", yaxt="n") 
  }
  
  if (raw==F & abs==T) { plot(w,z[1:(k+1)],type="l",xlab="number of dimensions",
                              cex.main=1.5, lwd=3, col="red",
                              ylim=c(0,100),
                              ylab="cumulative percentage of total absolutes eigenvalues",
                              main="Scree plot of absolute values of eigenvalues",
                              xaxt="n", yaxt="n") }
  
  if (raw==F & abs==F) { plot(w,z[1:(k+1)],type="l",xlab="number of dimensions",
                              cex.main=1.5, lwd=3, col="red",
                              ylim=c(0,100),
                              ylab="cumulative percentage of total variance",
                              main="Scree plot of squared eigenvalues",
                              xaxt="n", yaxt="n") }
  
  axis(1,at=w,lwd=2)
  axis(2,at=c(0,20,40,60,80,100),lwd=2)
  abline(a=100,b=0,lwd=2,lty="dashed",col="orange")
  text(w,z[1:(k+1)],labels=w,cex=0.8,adj=c(1.2,-.1),col="blue")
}
