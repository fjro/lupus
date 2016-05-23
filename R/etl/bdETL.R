require(gdata)

#read the data
genes = read.csv('genes.csv', header =T)
options( stringsAsFactors=F ) 
bd <- read.csv("Big Data.csv", header = T)
rownames(bd) = bd[,1]
bd = bd[,-1]
?read.csv2
typeof(bd[3,3])
table(genes$diseased)
samples <- read.table("Sample Map.csv", as.is = T)
dim(bd)
bd = t(bd)
bd = data.frame(t(bd), stringsAsFactors = F)
colnames(bd) = bd[1,]
bd=bd[-1,]
dim(bd)
length(bd[,1])
b

#pull out the sample names
controlNames = table(substr(samples[samples$Color=='GREEN',]$Sample,0,6))
lupusNames = table(substr(samples[samples$Color=='RED',]$Sample,0,6))
affyColNames = table(substr(colnames(bd),0,6))

genes$X.1
intersect(rownames(allBD), genes$X.1)
length(intersect(rownames(allBD), genes$Var.1[1:53]))
length(intersect(rownames(allBD), genes$Var.1[54:420]))

length(intersect(colnames(bd), genes$Var.1[1:53]))
length(intersect(colnames(bd), genes$Var.1[54:420]))

safe = which(rownames(bd) %in% intersect(rownames(bd), genes$X.1))
safeBDAll = bd[safe,]
dim(safeBDAll)
colnames(safeBDAll) = colnames(bd)
typeof(safeBDAll[1,1])

safeControlBD = bd[which(rownames(bd) %in% intersect(rownames(bd), genes$X.1[1:53])),]
dim(safeControlBD)
colnames(safeControlBD) = colnames(bd)

safeLupusBD = bd[which(rownames(bd) %in% intersect(rownames(bd), genes$X.1[54:420])),]
dim(safeLupusBD)
colnames(safeLupusBD) = colnames(bd)

rownames(safeControlBD) == rownames(safeBDAll)[1:29]
rownames(safeLupusBD) == rownames(safeBDAll)[30:397]

#separate controls and lupus
controlBD = sort(intersect(names(controlNames), names(affyColNames)))
lupusBD = sort(intersect(names(lupusNames), names(affyColNames)))


lupusSamples = c()
for(i in 1:length(lupusBD))
{
  lupusSamples = union(lupusSamples, colnames(bd)[which(startsWith(colnames(bd), lupusBD[i], trim=FALSE, ignore.case=FALSE))])
}
length(lupusSamples)

controlSamples = c()
for(i in 1:length(controlBD))
{
  controlSamples = union(controlSamples, colnames(bd)[which(startsWith(colnames(bd), controlBD[i], trim=FALSE, ignore.case=FALSE))])
}
length(controlSamples)
length(lupusSamples)

allBD = bd[,c(controlSamples, lupusSamples)]
allBD = t(allBD)
colnames(allBD) = bd$X
dim(allBD)
controlBD = bd[,controlSamples]
controlBD = t(controlBD)
colnames(controlBD) = bd$X
lupusBD = bd[,lupusSamples]
lupusBD = t(lupusBD)
colnames(lupusBD) = bd$X

#look at distributions
require(plyr)
require(reshape)
require(ggplot2)
rrr = cbind(safeBDAll, c(rep('Control', 29), rep('Lupus', nrow(safeLupusBD))))
colnames(rrr)[ncol(rrr)] = 'Type'
gM <- melt(rrr[,1:1000])
colnames(gM)[c(2,3)] <- c("Gene","Expression")
dp = ggplot(gM, aes(x=Expression)) + geom_line(aes(group=Gene),stat="density", alpha=0.15)
dp + facet_grid(Type ~ .)

gM <- melt(safeControlBD)
colnames(gM)[c(2,3)] <- c("Gene","Expression")
ggplot(gM, aes(x=Expression)) + geom_line(aes(group=Gene),stat="density", alpha=0.15)

