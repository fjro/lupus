require(ggplot2)
require(gridExtra)
getwd()
lupus.pca <- prcomp(lupus[,c(-1,-ncol(genes))], retx=TRUE,
                    center = TRUE,
                    scale. = TRUE) 

plot(lupus.pca)

screeplot(lupus[,c(-1,-ncol(genes))],maxcomp=25)
screeplot(genes[,c(-1,-ncol(genes))],maxcomp=25)




#project samples
#plot PC1 v PC2
pca.1 <- prcomp(genes[,-c(1,96,97)], retx=TRUE,center = TRUE, scale. = TRUE) 

screeplot(as.matrix(pca.1))
pcadf = as.data.frame(pca.1$x)
pcadf = cbind(pcadf, genes$diseased)
pcadf = cbind(pcadf, genes$Source)
colnames(pcadf)[95:96] = c("Patient","Source")
ggplot(pcadf, aes(x=PC1, y=PC2,color=Source, shape=Patient)) +geom_point(size =3) + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))
ggplot(pcadf, aes(x=PC3, y=PC4,color=Source, shape=Patient)) +geom_point(size =3)
ggplot(pcadf, aes(x=PC5, y=PC6,color=Source, shape=Patient)) +geom_point(size =3)
ggplot(pcadf, aes(x=PC7, y=PC8,color=Source, shape=Patient)) +geom_point(size =3)

loadings = as.data.frame(pca.1$rotation)
loadings = cbind(rownames(loadings),loadings)
colnames(loadings)[1] = 'Probe'
loadings = cbind(1:94, loadings)
colnames(loadings)[1] = 'ProbeNo'
ggplot(loadings, aes(y=PC1,x=ProbeNo))+ scale_x_discrete(limits = loadings$ProbeNo)+ 
  geom_bar(stat = "identity",position="dodge")+ theme(axis.text.x=element_text(angle=270, vjust = 0.4,size=6))

#project genes
pca.2 <- prcomp(t(genes[,-c(1,96,97)]), retx=TRUE,center = TRUE, scale. = TRUE) 

screeplot(as.matrix(t(genes[,-c(1,96,97)])))
pcadf2 = as.data.frame(pca.2$x)
pcadf2 = cbind(pcadf2, colnames(genes[-c(1,96,97)]))
#pcadf2 = cbind(pcadf2, genes$Type)
colnames(pcadf2)[95] = "Gene"
ggplot(pcadf, aes(x=PC1, y=PC2,color=Type, shape=Disease)) +geom_point(size =4) + theme(legend.position='bottom',legend.key.size = unit(.4, "cm"))
ggplot(pcadf2, aes(x=PC1, y=PC2, label=Gene)) + geom_point(alpha=.3) +geom_text() 
ggplot(pcadf2, aes(x=PC3, y=PC4, label=Gene)) + geom_point(alpha=.3) +geom_text() 
ggplot(pcadf2, aes(x=PC5, y=PC6, label=Gene)) + geom_point(alpha=.3) +geom_text() 

loadings2 = as.data.frame(pca.2$rotation)
loadings2 = cbind(genes$Var.1,loadings2)
colnames(loadings2)[1] = 'Sample'
ggplot(loadings2, aes(y=PC1,x=Sample))+ scale_x_discrete(limits = loadings2$Sample)+ geom_bar(stat = "identity",position="dodge")+ theme(axis.text.x=element_text(angle=270, vjust = 0.5))

